use clap::Parser;
use flate2::read::GzDecoder;
use std::collections::BTreeMap;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Write};
use thiserror::Error;

mod write;
use write::parallel_read_write;

type Result<T> = std::result::Result<T, CombineError>;

/// Entries below this value are never written to the output, and are not counted in row sums.
pub(crate) const STATIC_THRESHOLD: f64 = 0.000005;

/// Default for `--restrict` with no value: the largest number of non-zero entries a sparse matrix
/// can hold in R's Matrix package (.Machine$integer.max, 2^31-1).
const R_MAX_ENTRIES: &str = "2147483647";

/// The NOTICE file: this program's licence, and the notices of its third-party dependencies
const NOTICE: &str = include_str!("../NOTICE");

/// The opening statement of NOTICE, i.e. everything before the first third-party notice
fn get_notice_summary() -> &'static str {
    NOTICE
        .split_once("\n---")
        .map_or(NOTICE, |(summary, _)| summary)
}

#[derive(Debug, Error)]
enum CombineError {
    #[error(transparent)]
    IOError(#[from] io::Error),
    #[error(transparent)]
    ParseIntError(#[from] std::num::ParseIntError),
    #[error(transparent)]
    ParseFloatError(#[from] std::num::ParseFloatError),
    #[error(transparent)]
    GZipError(#[from] gzp::GzpError),
    #[error(transparent)]
    ThreadPoolBuildError(#[from] rayon::ThreadPoolBuildError),
    #[error("{0}")]
    Message(String),
}

// ---------------------------------------------------------------------------
// Sparse matrix
//
// HMat is a sparse matrix: a (BTree)Map of HVec rows.
// HVec is a sparse vector: only non-default entries are stored in a (BTree)Map.
// ---------------------------------------------------------------------------

struct HMat {
    m: BTreeMap<usize, HVec>,
}

impl HMat {
    fn new() -> Self {
        HMat { m: BTreeMap::new() }
    }

    fn get(&mut self, row: usize) -> &mut HVec {
        self.m.entry(row).or_insert_with(HVec::new)
    }
}

struct HVec {
    v: BTreeMap<usize, f64>,
}

impl HVec {
    fn new() -> Self {
        HVec { v: BTreeMap::new() }
    }

    fn add(&mut self, p: usize, delta: f64) {
        let current = self.v.entry(p).or_insert(0.0);
        *current += delta;
    }
}

/// How `--restrict` chooses which entries to drop.
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
enum RestrictMode {
    /// Drop the smallest values in the whole matrix (one cut-off for every row)
    #[default]
    Global,
    /// Keep the same fraction of every row's entries (each row's largest values)
    RowFraction,
}

// ---------------------------------------------------------------------------
// CLI entry point
// ---------------------------------------------------------------------------

#[derive(Parser, Debug)]
#[command(
    version,
    disable_version_flag = true,
    override_usage = "combine --chunkpathsfile <chunkpathsfile> --snpcountsfile <snpcountsfile> --maplengthsfile <maplengthsfile> [OPTIONS]"
)]
struct Args {
    /// A file containing the paths to the chunk files.
    #[arg(short, long, required_unless_present_any = ["licences", "licences_summary"], default_value = "")]
    chunkpathsfile: String,

    /// A file containing the SNP counts for each chunk, in the same order as the chunk paths file.
    #[arg(short, long, required_unless_present_any = ["licences", "licences_summary"], default_value = "")]
    snpcountsfile: String,

    /// A file containing the map lengths for each chunk, in the same order as the chunk paths file.
    #[arg(short, long, required_unless_present_any = ["licences", "licences_summary"], default_value = "")]
    maplengthsfile: String,

    /// The maximum number of entries (lines) to write to the output file. Which entries are dropped is set by --restrict-mode [Default with flag but no value: 2^31-1, the most R's Matrix package can hold]
    #[arg(short, long, default_missing_value = R_MAX_ENTRIES, num_args = 0..=1, require_equals = true)]
    restrict: Option<usize>,

    /// How --restrict drops entries when there are too many: "global" drops the smallest values in the whole matrix; "row-fraction" keeps the same fraction of each row's entries (its largest values), so that rows with small values are not dropped disproportionately
    #[arg(long, value_enum, default_value_t = RestrictMode::Global, requires = "restrict")]
    restrict_mode: RestrictMode,

    /// Write the sum of each row (index1) to <OUT>.sums, as "index1 sum", counting all values >= 0.000005 including any dropped by --restrict. Always written if --restrict is used.
    #[arg(short, long, required = false, default_value_t = false)]
    writesums: bool,

    /// The number of threads to use for writing.
    #[arg(short, long, required = false, default_value_t = 8)]
    threads: usize,

    /// The prefix for the output file(s).
    #[arg(short, long, required = false, default_value = "combined")]
    out: String,

    /// Print a summary of licence information
    #[arg(short = 'l', default_value_t = false)]
    licences_summary: bool,

    /// Print full licence information, including all third-party notices
    #[arg(long, aliases = ["licenses", "license"], default_value_t = false)]
    licences: bool,

    /// Print version
    #[arg(short = 'v', long, action = clap::ArgAction::Version)]
    version: Option<bool>,
}

fn main() {
    let args = Args::parse();

    if args.licences {
        print!("{NOTICE}");
        std::process::exit(0);
    }

    if args.licences_summary {
        println!("{}", get_notice_summary());
        println!(
            "Use \"combine --licences\" to print the full text, including all third-party notices."
        );
        std::process::exit(0);
    }

    run(args).unwrap_or_else(|e| {
        eprintln!("{}", e);
        std::process::exit(1);
    });
}

fn run(args: Args) -> Result<()> {
    let out = args.out.clone();
    combine(args, |writesums| {
        let hap_writer = BufWriter::new(File::create(format!("{out}.txt.gz"))?);
        let sums_writer = if writesums {
            Some(BufWriter::new(File::create(format!("{out}.sums"))?))
        } else {
            None
        };
        Ok((hap_writer, sums_writer))
    })
}

/// Does all the work of `run`, except creating the output writers. Once the input has been loaded,
/// `open_outputs` is called with whether the sums file is wanted, and returns the writers for the
/// matrix and (if wanted) the sums. Tests use this to write to memory instead of files.
fn combine<W, F>(mut args: Args, open_outputs: F) -> Result<()>
where
    W: Write + Send + 'static,
    F: FnOnce(bool) -> Result<(W, Option<W>)>,
{
    if args.threads < 1 {
        return Err(CombineError::Message(
            "invalid value for --threads <THREADS>: must be at least 1".to_string(),
        ));
    }

    let nrows = args.restrict;

    if nrows.is_some() {
        args.writesums = true;
    }

    let chr_filenames = read_values_per_line::<String>(&args.chunkpathsfile)?;
    let total_snp = read_values_per_line::<f64>(&args.snpcountsfile)?;
    let total_gd = read_values_per_line::<f64>(&args.maplengthsfile)?;

    if chr_filenames.len() != total_snp.len() || chr_filenames.len() != total_gd.len() {
        return Err(CombineError::Message(format!(
            "{}, {}, and {} must all have the same number of lines (one per input file)",
            &args.chunkpathsfile, &args.snpcountsfile, &args.maplengthsfile
        )));
    }

    // Per-chromosome weight: genetic map length / SNP count, rounded to 6 d.p.
    let weights: Vec<f64> = total_gd
        .iter()
        .zip(total_snp.iter())
        .map(|(gd, snp)| (gd / snp * 1e6).round() / 1e6)
        .collect();

    // Accumulate weighted chunk lengths across all chromosomes.
    let mut cl = HMat::new();
    for (i, path) in chr_filenames.iter().enumerate() {
        eprintln!("Loading data from {}", filename_from_path(path));
        read_pbwt_out(path, &mut cl, weights[i])?;
    }
    eprintln!("Loading data complete");

    // -----------------------------------------------------------------------
    // Determine dynamic thresholds (one per row, in row order).
    //
    // If --restrict=N is set and more than N entries survive the static
    // threshold, find cut-offs so that only entries strictly above them are
    // written, leaving at most N entries: one cut-off for the whole matrix
    // (--restrict-mode=global), or one per row (--restrict-mode=row-fraction).
    // If there are ties at a cut-off, fewer than N entries may be written.
    // -----------------------------------------------------------------------
    let dynamic_thresholds: Vec<f64> = match (nrows, args.restrict_mode) {
        (None, _) => {
            eprintln!("No dynamic threshold needed");
            vec![f64::NEG_INFINITY; cl.m.len()]
        }
        (Some(n), RestrictMode::Global) => vec![get_dynamic_threshold(&cl, n); cl.m.len()],
        (Some(n), RestrictMode::RowFraction) => get_row_thresholds(&cl, n),
    };

    // -----------------------------------------------------------------------
    // Write the output matrix, applying both thresholds.
    //
    // An entry is written iff:
    //   val >= STATIC_THRESHOLD  (always applied)
    //   val > its row's dynamic threshold  (-inf unless --restrict needs to drop entries)
    // -----------------------------------------------------------------------

    let (hap_writer, sums_writer) = open_outputs(args.writesums)?;

    eprintln!("Writing output matrix");
    parallel_read_write(
        args.threads,
        cl,
        hap_writer,
        sums_writer,
        dynamic_thresholds,
    )?;
    Ok(())
}

/// Returns the value at or below which entries are excluded from the output, so that at most
/// `max_entries` entries are written. Returns -infinity if no entries need to be excluded.
fn get_dynamic_threshold(cl: &HMat, max_entries: usize) -> f64 {
    eprintln!("Determining dynamic threshold");
    // Only count entries that will survive the static threshold. This also excludes NaN.
    let mut all_values: Vec<f64> =
        cl.m.values()
            .flat_map(|row| row.v.values())
            .copied()
            .filter(|&v| v >= STATIC_THRESHOLD)
            .collect();

    let total_length = all_values.len();

    let dynamic_threshold = if total_length > max_entries {
        let n_to_drop = total_length - max_entries;
        // After select_nth_unstable_by(k), the element at k is the one that would be at sorted
        // position k. Keeping only entries > it drops at least the n_to_drop smallest entries.
        let threshold = *all_values
            .select_nth_unstable_by(n_to_drop - 1, |a, b| a.total_cmp(b))
            .1;
        eprintln!(
            "Note: {} entries survive the static threshold, which is more than the {} allowed by \
             --restrict; dropping entries <= {:.6} (at least {} entries).",
            total_length, max_entries, threshold, n_to_drop
        );
        threshold
    } else {
        eprintln!(
            "{} entries survive the static threshold, no more than the {} allowed by --restrict; \
             no dynamic threshold needed.",
            total_length, max_entries
        );
        f64::NEG_INFINITY
    };

    eprintln!("Determining dynamic threshold complete");
    dynamic_threshold
}

/// Returns a cut-off for each row (in row order) for --restrict-mode=row-fraction. If more than
/// `max_entries` entries survive the static threshold, a row with n such entries keeps its largest
/// floor(n * max_entries / total) entries (those strictly above its cut-off, so ties at the cut-off
/// can make it fewer), so that at most `max_entries` are written in total. Otherwise every cut-off
/// is -infinity.
fn get_row_thresholds(cl: &HMat, max_entries: usize) -> Vec<f64> {
    eprintln!("Determining per-row thresholds");
    let passing = |row: &HVec| -> Vec<f64> {
        row.v
            .values()
            .copied()
            .filter(|&v| v >= STATIC_THRESHOLD)
            .collect()
    };
    let total_length: usize =
        cl.m.values()
            .map(|row| row.v.values().filter(|&&v| v >= STATIC_THRESHOLD).count())
            .sum();

    if total_length <= max_entries {
        eprintln!(
            "{} entries survive the static threshold, no more than the {} allowed by --restrict; \
             no dynamic thresholds needed.",
            total_length, max_entries
        );
        return vec![f64::NEG_INFINITY; cl.m.len()];
    }

    let mut n_kept = 0;
    let thresholds =
        cl.m.values()
            .map(|row| {
                let mut values = passing(row);
                let n = values.len();
                let keep = (n as u128 * max_entries as u128 / total_length as u128) as usize;
                n_kept += keep;
                if keep >= n {
                    f64::NEG_INFINITY
                } else {
                    // Keeping only entries > the (n - keep)th smallest drops at least n - keep entries.
                    *values
                        .select_nth_unstable_by(n - keep - 1, |a, b| a.total_cmp(b))
                        .1
                }
            })
            .collect();

    eprintln!(
        "Note: {} entries survive the static threshold, which is more than the {} allowed by \
         --restrict; keeping the largest {:.4}% of each row's entries (at most {} entries).",
        total_length,
        max_entries,
        100.0 * max_entries as f64 / total_length as f64,
        n_kept
    );
    eprintln!("Determining per-row thresholds complete");
    thresholds
}

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

/// Read one whitespace-trimmed token per non-blank line into a Vec<T>.
fn read_values_per_line<T>(path: &str) -> Result<Vec<T>>
where
    T: std::str::FromStr,
    T::Err: std::fmt::Debug,
{
    let file = open_file_or_err(path)?;
    let reader = BufReader::new(file);
    let mut values = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let value: T = trimmed.parse().map_err(|e| {
            CombineError::Message(format!(
                "could not parse '{}' from {}: {:?}",
                trimmed, path, e
            ))
        })?;
        values.push(value);
    }
    Ok(values)
}

/// Read a per-chromosome painting file, accumulating value * weight
/// into existing entries or creating new ones.
fn read_pbwt_out(path: &str, cl: &mut HMat, weight: f64) -> Result<()> {
    let file = open_file_or_err(path)?;
    let reader = BufReader::new(GzDecoder::new(file));

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<_> = line.split_ascii_whitespace().collect();
        if parts.len() != 3 {
            return Err(CombineError::Message(format!(
                "invalid line: {}, in {}",
                line,
                filename_from_path(path)
            )));
        }
        let ind1: usize = parts[0].parse()?;
        let ind2: usize = parts[1].parse()?;
        let value: f64 = parts[2].parse()?;

        cl.get(ind1 - 1).add(ind2 - 1, value * weight);
    }

    Ok(())
}

fn open_file_or_err(path: &str) -> Result<File> {
    File::open(path).map_err(|e| {
        CombineError::Message(format!(
            "could not open {}: {}",
            filename_from_path(path),
            e
        ))
    })
}

fn filename_from_path(path: &str) -> &str {
    path.split('/').next_back().unwrap_or(path)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use std::cell::Cell;
    use std::io::BufReader;
    use std::path::PathBuf;
    use std::sync::{Arc, Mutex};

    #[test]
    fn verify_cli() {
        use clap::CommandFactory;
        Args::command().debug_assert();
    }

    #[test]
    fn test_notice_summary_stops_before_first_third_party_notice() {
        let summary = get_notice_summary();
        assert!(summary.starts_with("combine is licensed under the GNU General Public License"));
        assert!(!summary.contains("---"));
        assert!(summary.len() < NOTICE.len());
    }

    fn testdata(filename: &str) -> String {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("testdata")
            .join(filename)
            .to_string_lossy()
            .into_owned()
    }

    /// Read bytes (optionally gzipped) and return the lines.
    fn bytes_2_lines(bytes: &[u8], gz: bool) -> Vec<String> {
        if gz {
            BufReader::new(GzDecoder::new(bytes))
                .lines()
                .map(|l| l.unwrap())
                .collect()
        } else {
            BufReader::new(bytes).lines().map(|l| l.unwrap()).collect()
        }
    }

    fn testdata_lines(filename: &str) -> Vec<String> {
        bytes_2_lines(
            &std::fs::read(testdata(filename)).unwrap(),
            filename.ends_with(".gz"),
        )
    }

    /// An in-memory writer that can be handed to `combine` (which needs `Send + 'static`) and read
    /// back afterwards through a clone.
    #[derive(Clone, Default)]
    struct SharedBuf(Arc<Mutex<Vec<u8>>>);

    impl Write for SharedBuf {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buf);
            Ok(buf.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// The lines of the matrix and, if it was written, the sums file.
    struct Output {
        entries: Vec<String>,
        sums: Option<Vec<String>>,
    }

    /// Run combine on testdata with the given extra arguments, writing to memory.
    fn run_in_memory(extra_args: &[&str]) -> Output {
        let mut argv = vec![
            "combine".to_string(),
            "--chunkpathsfile".to_string(),
            testdata("chunklength.files.txt"),
            "--snpcountsfile".to_string(),
            testdata("nsnps.txt"),
            "--maplengthsfile".to_string(),
            testdata("map_lengths.txt"),
        ];
        argv.extend(extra_args.iter().map(|a| a.to_string()));
        let args = Args::parse_from(argv);

        let hap = SharedBuf::default();
        let sums = SharedBuf::default();
        let sums_wanted = Cell::new(false);
        combine(args, |writesums| {
            sums_wanted.set(writesums);
            Ok((hap.clone(), writesums.then(|| sums.clone())))
        })
        .unwrap();

        Output {
            entries: bytes_2_lines(&hap.0.lock().unwrap(), true),
            sums: sums_wanted
                .get()
                .then(|| bytes_2_lines(&sums.0.lock().unwrap(), false)),
        }
    }

    /// Return the number of unique row indices in the output data.
    fn get_nrows(data: &[String]) -> usize {
        data.iter()
            .map(|s| {
                s.split_ascii_whitespace()
                    .next()
                    .unwrap()
                    .parse::<usize>()
                    .unwrap()
            })
            .unique()
            .count()
    }

    /// Number of individuals (rows) in testdata
    const N_INDS: usize = 3202;

    /// Number of entries in testdata/expected.rust.txt.gz
    const N_EXPECTED_ENTRIES: usize = 1929576;

    /// Parse "index sum" lines from a .sums file.
    fn parse_sums(lines: &[String]) -> Vec<(usize, f64)> {
        lines
            .iter()
            .map(|l| {
                let (i, s) = l.split_once(' ').unwrap();
                (i.parse().unwrap(), s.parse().unwrap())
            })
            .collect()
    }

    /// Check that a .sums file has one line for every row in testdata, in order, and that the sums
    /// match those calculated from the C++ implementation's output (up to floating point error).
    fn check_sums(lines: &[String]) {
        let sums = parse_sums(lines);
        let expected = testdata_lines("expected.cpp.sorted.sums");
        assert_eq!(sums.len(), N_INDS);
        assert_eq!(expected.len(), N_INDS);
        for (k, ((i, s), e)) in sums.iter().zip(&expected).enumerate() {
            assert_eq!(*i, k + 1);
            let e: f64 = e.parse().unwrap();
            assert!((s - e).abs() < 1e-3, "row {i}: sum {s}, expected {e}");
        }
    }

    #[test]
    fn test_combine_matches_expected() {
        let out = run_in_memory(&[]);
        assert_eq!(out.entries, testdata_lines("expected.rust.txt.gz"));
        assert!(out.sums.is_none());
    }

    #[test]
    fn test_combine_matches_expected_with_sums() {
        let out = run_in_memory(&["--writesums"]);
        assert_eq!(out.entries, testdata_lines("expected.rust.txt.gz"));
        let sums = out.sums.unwrap();
        assert_eq!(get_nrows(&out.entries), sums.len());
        check_sums(&sums);
    }

    #[test]
    fn test_restrict_10000() {
        let out = run_in_memory(&["--restrict=10000"]);
        assert!(out.entries.len() <= 10000);
        // rows with no entries left after --restrict are still in the sums file, with their full sums
        assert!(get_nrows(&out.entries) < N_INDS);
        check_sums(&out.sums.unwrap());
    }

    #[test]
    fn test_restrict_equal_to_entry_count_keeps_everything() {
        let out = run_in_memory(&[&format!("--restrict={N_EXPECTED_ENTRIES}")]);
        assert_eq!(out.entries, testdata_lines("expected.rust.txt.gz"));
    }

    #[test]
    fn test_restrict_one_below_entry_count_drops_smallest() {
        let out = run_in_memory(&[&format!("--restrict={}", N_EXPECTED_ENTRIES - 1)]);
        let expected = testdata_lines("expected.rust.txt.gz");
        assert!(out.entries.len() < N_EXPECTED_ENTRIES);
        check_sums(&out.sums.unwrap());

        let value = |l: &String| l.split_whitespace().nth(2).unwrap().parse::<f64>().unwrap();
        let min_kept = out.entries.iter().map(value).fold(f64::INFINITY, f64::min);
        let kept: std::collections::HashSet<&String> = out.entries.iter().collect();
        let max_expected_dropped = expected
            .iter()
            .filter(|l| !kept.contains(l))
            .map(value)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(max_expected_dropped <= min_kept);
    }

    #[test]
    fn test_restrict_zero_writes_no_entries_but_all_sums() {
        let out = run_in_memory(&["--restrict=0"]);
        assert!(out.entries.is_empty());
        check_sums(&out.sums.unwrap());
    }

    #[test]
    fn test_combine_matches_expected_cpp() {
        let out = run_in_memory(&["--writesums"]);
        let expected = testdata_lines("expected.cpp.sorted.txt.gz");
        assert_eq!(out.entries.len(), expected.len());

        // there are small differences due to floating point precision between the rust
        // implementation and the cpp implementation
        let split = |l: &String| {
            let f: Vec<&str> = l.split_whitespace().collect();
            let ij: (usize, usize) = (f[0].parse().unwrap(), f[1].parse().unwrap());
            (ij, f[2].parse::<f64>().unwrap())
        };
        for (a, e) in out.entries.iter().zip(&expected) {
            let (a_ij, a_val) = split(a);
            let (e_ij, e_val) = split(e);
            assert_eq!(a_ij, e_ij);
            assert!((a_val - e_val).abs() < 1e-3, "{a} vs {e}");
        }

        check_sums(&out.sums.unwrap());
    }

    #[test]
    fn test_restrict_mode_requires_restrict() {
        let base = ["combine", "-c", "a", "-s", "b", "-m", "c"];
        let without = Args::try_parse_from(base.iter().chain(&["--restrict-mode=row-fraction"]));
        assert!(without.is_err());
        let with = Args::try_parse_from(
            base.iter()
                .chain(&["--restrict=10", "--restrict-mode=row-fraction"]),
        )
        .unwrap();
        assert_eq!(with.restrict_mode, RestrictMode::RowFraction);
        assert_eq!(
            Args::try_parse_from(base).unwrap().restrict_mode,
            RestrictMode::Global
        );
    }

    #[test]
    fn test_restrict_mode_global_is_the_default() {
        let explicit = run_in_memory(&["--restrict=10000", "--restrict-mode=global"]);
        let default = run_in_memory(&["--restrict=10000"]);
        assert_eq!(explicit.entries, default.entries);
        assert_eq!(explicit.sums, default.sums);
    }

    #[test]
    fn test_row_fraction_keeps_everything_when_under_the_limit() {
        let out = run_in_memory(&[
            &format!("--restrict={N_EXPECTED_ENTRIES}"),
            "--restrict-mode=row-fraction",
        ]);
        assert_eq!(out.entries, testdata_lines("expected.rust.txt.gz"));
        check_sums(&out.sums.unwrap());
    }

    #[test]
    fn test_row_fraction_zero_writes_no_entries_but_all_sums() {
        let out = run_in_memory(&["--restrict=0", "--restrict-mode=row-fraction"]);
        assert!(out.entries.is_empty());
        check_sums(&out.sums.unwrap());
    }

    /// Group "i j value" lines by row: row -> values.
    fn values_by_row(lines: &[String]) -> BTreeMap<usize, Vec<f64>> {
        let mut rows: BTreeMap<usize, Vec<f64>> = BTreeMap::new();
        for l in lines {
            let f: Vec<&str> = l.split_whitespace().collect();
            rows.entry(f[0].parse().unwrap())
                .or_default()
                .push(f[2].parse().unwrap());
        }
        rows
    }

    /// With --restrict-mode=row-fraction, each row keeps at most floor(n * N / total) of its n
    /// entries, the ones it keeps are its largest, and the total is at most N.
    fn check_row_fraction(max_entries: usize) -> Output {
        let out = run_in_memory(&[
            &format!("--restrict={max_entries}"),
            "--restrict-mode=row-fraction",
        ]);
        let expected = testdata_lines("expected.rust.txt.gz");
        assert!(out.entries.len() <= max_entries);
        let kept_set: std::collections::HashSet<&String> = out.entries.iter().collect();
        assert!(out.entries.iter().all(|l| kept_set.contains(l)));

        let all = values_by_row(&expected);
        let kept = values_by_row(&out.entries);
        let dropped = values_by_row(
            &expected
                .iter()
                .filter(|l| !kept_set.contains(l))
                .cloned()
                .collect::<Vec<_>>(),
        );
        let mut total_budget = 0;
        for (i, values) in &all {
            let budget = values.len() * max_entries / N_EXPECTED_ENTRIES;
            total_budget += budget;
            let k = kept.get(i).map_or(&[][..], |v| v.as_slice());
            assert!(k.len() <= budget, "row {i}: kept {} > {budget}", k.len());
            let min_kept = k.iter().copied().fold(f64::INFINITY, f64::min);
            let max_dropped = dropped.get(i).map_or(f64::NEG_INFINITY, |v| {
                v.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            });
            assert!(
                max_dropped <= min_kept,
                "row {i}: dropped {max_dropped} > kept {min_kept}"
            );
        }
        // rows keep fewer than their budget only when values tie at their cut-off: rare
        let shortfall = total_budget - out.entries.len();
        assert!(
            shortfall * 1000 < total_budget,
            "{shortfall} of {total_budget} entries lost to ties"
        );
        check_sums(out.sums.as_ref().unwrap());
        out
    }

    #[test]
    fn test_row_fraction_half() {
        check_row_fraction(N_EXPECTED_ENTRIES / 2);
    }

    #[test]
    fn test_row_fraction_10000_keeps_more_rows_than_global() {
        let row_fraction = check_row_fraction(10000);
        let global = run_in_memory(&["--restrict=10000"]);
        // every row with at least N_EXPECTED_ENTRIES / 10000 entries keeps at least one
        let all = values_by_row(&testdata_lines("expected.rust.txt.gz"));
        let kept = values_by_row(&row_fraction.entries);
        for (i, values) in &all {
            if values.len() * 10000 >= N_EXPECTED_ENTRIES {
                assert!(kept.contains_key(i), "row {i} lost all its entries");
            }
        }
        assert!(get_nrows(&row_fraction.entries) > get_nrows(&global.entries));
    }

    /// End-to-end check that `run` writes `<out>.txt.gz` and `<out>.sums`, in a temporary directory.
    #[test]
    fn test_run_writes_output_files() {
        let dir = std::env::temp_dir().join(format!("combine-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("combined").to_string_lossy().into_owned();
        let args = Args::parse_from([
            "combine",
            "--chunkpathsfile",
            &testdata("chunklength.files.txt"),
            "--snpcountsfile",
            &testdata("nsnps.txt"),
            "--maplengthsfile",
            &testdata("map_lengths.txt"),
            "--restrict=10000",
            "--out",
            &prefix,
        ]);
        let result = run(args);
        let entries = std::fs::read(format!("{prefix}.txt.gz"));
        let sums = std::fs::read(format!("{prefix}.sums"));
        std::fs::remove_dir_all(&dir).unwrap();

        result.unwrap();
        let in_memory = run_in_memory(&["--restrict=10000"]);
        assert_eq!(bytes_2_lines(&entries.unwrap(), true), in_memory.entries);
        assert_eq!(
            bytes_2_lines(&sums.unwrap(), false),
            in_memory.sums.unwrap()
        );
    }
}
