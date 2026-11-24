use clap::Parser;
use flate2::read::GzDecoder;
use std::collections::BTreeMap;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter};
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

    /// The maximum number of entries (lines) to write to the output file. The smallest values are dropped first [Default with flag but no value: 2^31-1, the most R's Matrix package can hold]
    #[arg(short, long, default_missing_value = R_MAX_ENTRIES, num_args = 0..=1, require_equals = true)]
    restrict: Option<usize>,

    /// Write the sums across the matrix before any filtering. If --restrict is in effect, this file will be written anyway.
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

fn run(mut args: Args) -> Result<()> {
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
    // Determine dynamic threshold.
    //
    // If --restrict=N is set and more than N entries survive the static
    // threshold, find a cut-off so that only entries strictly above it are
    // written, leaving at most N entries. If there are ties at the cut-off,
    // fewer than N entries may be written.
    // -----------------------------------------------------------------------
    let dynamic_threshold: f64 = if let Some(n) = nrows {
        get_dynamic_threshold(&cl, n)
    } else {
        eprintln!("No dynamic threshold needed");
        f64::NEG_INFINITY
    };

    // -----------------------------------------------------------------------
    // Write the output matrix, applying both thresholds.
    //
    // An entry is written iff:
    //   val >= STATIC_THRESHOLD  (always applied)
    //   val > dynamic_threshold  (-inf unless --restrict needs to drop entries)
    // -----------------------------------------------------------------------

    let output_path = format!("{}.txt.gz", &args.out);
    let output_file = File::create(&output_path)?;
    let hap_writer = BufWriter::new(Box::new(output_file));

    let sums_writer = if args.writesums {
        Some(BufWriter::new(Box::new(File::create(format!(
            "{}.sums",
            &args.out
        ))?)))
    } else {
        None
    };

    eprintln!("Writing output matrix");
    parallel_read_write(args.threads, cl, hap_writer, sums_writer, dynamic_threshold)?;
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
    use std::io::BufReader;
    use std::path::PathBuf;

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

    /// Read a file (optionally gzipped) and return the lines as a Vec<String>.
    fn file_2_lines(bytes: &[u8], gz: bool) -> Vec<String> {
        if gz {
            BufReader::new(GzDecoder::new(bytes))
                .lines()
                .map(|l| l.unwrap())
                .collect()
        } else {
            BufReader::new(bytes).lines().map(|l| l.unwrap()).collect()
        }
    }

    /// Return the number of unique row indices in the output data.
    fn get_nrows(data: Vec<String>) -> usize {
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

    #[test]
    #[rustfmt::skip]
    fn test_combine_matches_expected() {
        let args = Args::parse_from([
            "combine",
            "--chunkpathsfile", &testdata("chunklength.files.txt"),
            "--snpcountsfile", &testdata("nsnps.txt"),
            "--maplengthsfile", &testdata("map_lengths.txt"),
            "--out", &testdata("combined"),
        ]);

        run(args).unwrap();

        let actual = file_2_lines(&std::fs::read(testdata("combined.txt.gz")).unwrap(), true);
        let expected = file_2_lines(&std::fs::read(testdata("expected.rust.txt.gz")).unwrap(), true);

        assert_eq!(actual, expected);

        std::fs::remove_file(testdata("combined.txt.gz")).unwrap();
    }

    #[test]
    #[rustfmt::skip]
    fn test_combine_matches_expected_r() {
        let args = Args::parse_from([
            "combine",
            "--chunkpathsfile", &testdata("chunklength.files.txt"),
            "--snpcountsfile", &testdata("nsnps.txt"),
            "--maplengthsfile", &testdata("map_lengths.txt"),
            "--writesums",
            "--out", &testdata("combined.r"),
        ]);

        run(args).unwrap();

        let actual = file_2_lines(&std::fs::read(testdata("combined.r.txt.gz")).unwrap(), true);
        let expected = file_2_lines(&std::fs::read(testdata("expected.rust.txt.gz")).unwrap(), true);

        assert_eq!(actual, expected);

        let actual_sums = file_2_lines(&std::fs::read(testdata("combined.r.sums")).unwrap(), false);

        assert_eq!(get_nrows(actual), actual_sums.len());

        std::fs::remove_file(testdata("combined.r.txt.gz")).unwrap();
        std::fs::remove_file(testdata("combined.r.sums")).unwrap();
    }

    #[test]
    #[rustfmt::skip]
    fn test_combine_matches_expected_r_10000() {
        let args = Args::parse_from([
            "combine",
            "--chunkpathsfile", &testdata("chunklength.files.txt"),
            "--snpcountsfile", &testdata("nsnps.txt"),
            "--maplengthsfile", &testdata("map_lengths.txt"),
            "--restrict=10000",
            "--out", &testdata("combined.r.10000"),
        ]);

        run(args).unwrap();

        let actual = file_2_lines(&std::fs::read(testdata("combined.r.10000.txt.gz")).unwrap(), true);
        assert!(actual.len() <= 10000);

        let actual_sums = file_2_lines(&std::fs::read(testdata("combined.r.10000.sums")).unwrap(), false);

        assert_eq!(get_nrows(actual), actual_sums.len());

        std::fs::remove_file(testdata("combined.r.10000.txt.gz")).unwrap();
        std::fs::remove_file(testdata("combined.r.10000.sums")).unwrap();
    }

    /// Number of entries in testdata/expected.rust.txt.gz
    const N_EXPECTED_ENTRIES: usize = 1929576;

    fn run_restricted(restrict: usize, prefix: &str) -> (Vec<String>, Vec<String>) {
        let restrict_arg = format!("--restrict={restrict}");
        let args = Args::parse_from([
            "combine",
            "--chunkpathsfile",
            &testdata("chunklength.files.txt"),
            "--snpcountsfile",
            &testdata("nsnps.txt"),
            "--maplengthsfile",
            &testdata("map_lengths.txt"),
            &restrict_arg,
            "--out",
            &testdata(prefix),
        ]);
        run(args).unwrap();
        let out = testdata(&format!("{prefix}.txt.gz"));
        let sums = testdata(&format!("{prefix}.sums"));
        let actual = file_2_lines(&std::fs::read(&out).unwrap(), true);
        let actual_sums = file_2_lines(&std::fs::read(&sums).unwrap(), false);
        std::fs::remove_file(out).unwrap();
        std::fs::remove_file(sums).unwrap();
        (actual, actual_sums)
    }

    #[test]
    fn test_restrict_equal_to_entry_count_keeps_everything() {
        let (actual, _) = run_restricted(N_EXPECTED_ENTRIES, "combined.r.exact");
        let expected = file_2_lines(
            &std::fs::read(testdata("expected.rust.txt.gz")).unwrap(),
            true,
        );
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_restrict_one_below_entry_count_drops_smallest() {
        let (actual, actual_sums) = run_restricted(N_EXPECTED_ENTRIES - 1, "combined.r.minus1");
        let expected = file_2_lines(
            &std::fs::read(testdata("expected.rust.txt.gz")).unwrap(),
            true,
        );
        assert!(actual.len() < N_EXPECTED_ENTRIES);
        assert_eq!(get_nrows(actual.clone()), actual_sums.len());

        let value = |l: &String| l.split_whitespace().nth(2).unwrap().parse::<f64>().unwrap();
        let min_kept = actual.iter().map(value).fold(f64::INFINITY, f64::min);
        let kept: std::collections::HashSet<&String> = actual.iter().collect();
        let max_expected_dropped = expected
            .iter()
            .filter(|l| !kept.contains(l))
            .map(value)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(max_expected_dropped <= min_kept);
    }

    #[test]
    fn test_restrict_zero_writes_nothing() {
        let (actual, actual_sums) = run_restricted(0, "combined.r.zero");
        assert!(actual.is_empty());
        assert!(actual_sums.is_empty());
    }

    #[test]
    #[rustfmt::skip]
    fn test_combine_matches_expected_cpp() {
        let args = Args::parse_from([
            "combine",
            "--chunkpathsfile", &testdata("chunklength.files.txt"),
            "--snpcountsfile", &testdata("nsnps.txt"),
            "--maplengthsfile", &testdata("map_lengths.txt"),
            "--writesums",
            "--out", &testdata("combined.cpp"),
        ]);

        run(args).unwrap();

        let actual = file_2_lines(&std::fs::read(testdata("combined.cpp.txt.gz")).unwrap(), true);
        let expected = file_2_lines(&std::fs::read(testdata("expected.cpp.sorted.txt.gz")).unwrap(), true);

        assert_eq!(actual.len(), expected.len());

        // there are rounding differences due to floating point precision between the rust implementation and the cpp implementation
        let actual_values_rounded = actual.iter().map(|l| l.split_whitespace().next().unwrap().parse::<f64>().unwrap().round()).collect::<Vec<_>>();
        let expected_values_rounded = expected.iter().map(|l| l.split_whitespace().next().unwrap().parse::<f64>().unwrap().round()).collect::<Vec<_>>();

        assert_eq!(actual_values_rounded, expected_values_rounded);

        let actual_sums = file_2_lines(&std::fs::read(testdata("combined.cpp.sums")).unwrap(), false);
        let expected_sums = file_2_lines(&std::fs::read(testdata("expected.cpp.sorted.sums")).unwrap(), false);

        assert_eq!(actual_sums.len(), expected_sums.len());

        // there are rounding differences due to floating point precision between the rust implementation and the cpp implementation
        let actual_sums_rounded = actual_sums.iter().map(|l| l.parse::<f64>().unwrap().round()).collect::<Vec<_>>();
        let expected_sums_rounded = expected_sums.iter().map(|l| l.parse::<f64>().unwrap().round()).collect::<Vec<_>>();

        assert_eq!(actual_sums_rounded, expected_sums_rounded);

        std::fs::remove_file(testdata("combined.cpp.txt.gz")).unwrap();
        std::fs::remove_file(testdata("combined.cpp.sums")).unwrap();
    }
}
