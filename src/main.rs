use clap::Parser;
use flate2::read::GzDecoder;
use std::collections::BTreeMap;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter};
use thiserror::Error;

mod licences;
use licences::print_licences;

mod write;
use write::parallel_read_write;

type Result<T> = std::result::Result<T, CombineError>;

#[derive(Debug, Error)]
enum CombineError {
    #[error(transparent)]
    IOError(#[from] io::Error),
    #[error(transparent)]
    ParseIntError(#[from] std::num::ParseIntError),
    #[error(transparent)]
    ParseFloatError(#[from] std::num::ParseFloatError),
    #[error(transparent)]
    ChanSendRowLocationErr(#[from] crossbeam_channel::SendError<write::RowLocation>),
    #[error(transparent)]
    ChanSendWriteDataErr(#[from] crossbeam_channel::SendError<write::WriteData>),
    #[error("lock poisoned")]
    LockPoisoned,
    #[error(transparent)]
    GZipError(#[from] gzp::GzpError),
    #[error("thread panicked: {0}")]
    ThreadPanic(String),
    #[error("{0}")]
    Message(String),
}

impl<T> From<std::sync::PoisonError<T>> for CombineError {
    fn from(_: std::sync::PoisonError<T>) -> Self {
        CombineError::LockPoisoned
    }
}

impl From<Box<dyn std::any::Any + Send + 'static>> for CombineError {
    fn from(e: Box<dyn std::any::Any + Send + 'static>) -> Self {
        let msg = e
            .downcast_ref::<&str>()
            .map(|s| s.to_string())
            .or_else(|| e.downcast_ref::<String>().cloned())
            .unwrap_or_else(|| "unknown panic".to_string());
        CombineError::ThreadPanic(msg)
    }
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
    override_usage = "combine --chunkpathsfile <chunkpathsfile> --snpcountsfile <snpcountsfile> --maplengthsfile <maplengthsfile> [OPTIONS]"
)]
struct Args {
    /// A file containing the paths to the chunk files.
    #[arg(short, long, required_unless_present = "licences", default_value = "")]
    chunkpathsfile: String,

    /// A file containing the SNP counts for each chunk, in the same order as the chunk paths file.
    #[arg(short, long, required_unless_present = "licences", default_value = "")]
    snpcountsfile: String,

    /// A file containing the map lengths for each chunk, in the same order as the chunk paths file.
    #[arg(short, long, required_unless_present = "licences", default_value = "")]
    maplengthsfile: String,

    /// The maximum number of rows to write to the output file [Default with flag but no value: 2^31-1]
    #[arg(short, long, default_missing_value = "2147483648", num_args = 0..=1, require_equals = true)]
    restrictrows: Option<usize>,

    /// Write the row sums. If --restrictrows is in effect, this file will be written anyway.
    #[arg(short, long, required = false, default_value_t = false)]
    writerowsums: bool,

    /// The number of threads to use for writing.
    #[arg(short, long, required = false, default_value_t = 8)]
    threads: usize,

    /// The prefix for the output file(s).
    #[arg(short, long, required = false, default_value = "combined")]
    out: String,

    /// Print the licences and exit.
    #[arg(short, long, alias = "license", default_value_t = false)]
    licences: bool,
}

fn main() {
    let args = Args::parse();

    if args.licences {
        print_licences();
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

    let nrows = args.restrictrows;

    if nrows.is_some() {
        args.writerowsums = true;
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
    // If --restrictrows is set, and the number of entries surviving the static
    // threshold is >= args.restrictrows, we must drop additional entries (the
    // smallest ones) until fewer than args.restrictrows remain. We find the exact
    // cut-off value using a linear-time selection algorithm, then keep only
    // entries strictly above it.
    //
    // Note: if there are ties at the threshold value, more than the minimum
    // number of entries may be dropped, but the count is guaranteed < args.restrictrows.
    // -----------------------------------------------------------------------
    let dynamic_threshold: f64 = if let Some(n) = nrows {
        get_dynamic_threshold(&cl, n)
    } else {
        eprintln!("No dynamic threshold needed");
        f64::NEG_INFINITY
    };

    // -----------------------------------------------------------------------
    // Pass 2: write the output matrix, applying both thresholds.
    //
    // An entry is written iff:
    //   val >= 0.000005          (static threshold, always applied)
    //   val > dynamic_threshold  (only active when total entries >= args.restrictrows)
    // -----------------------------------------------------------------------

    let output_path = format!("{}.txt.gz", &args.out);
    let output_file = File::create(&output_path)?;
    let hap_writer = BufWriter::new(Box::new(output_file));

    let rowsums_writer = if args.writerowsums {
        Some(BufWriter::new(Box::new(File::create(format!(
            "{}.rowsums",
            &args.out
        ))?)))
    } else {
        None
    };

    eprintln!("Writing output matrix");
    parallel_read_write(
        args.threads,
        cl,
        hap_writer,
        rowsums_writer,
        dynamic_threshold,
        args.writerowsums,
    )?;

    Ok(())
}

/// Returns threshold below which to exclude values from the output if there are more
/// rows in the input than are permitted
fn get_dynamic_threshold(cl: &HMat, nrows: usize) -> f64 {
    eprintln!("Determining dynamic threshold");
    let mut all_values: Vec<f64> =
        cl.m.values()
            .flat_map(|row| row.v.values())
            .copied()
            .collect();

    let total_length = all_values.len();

    let dynamic_threshold: f64 = if all_values.len() >= nrows {
        let n_to_drop = total_length - nrows;
        // select_nth_unstable_by(k) rearranges all_values so that the element
        // at sorted position k is in place. Elements before it are <= it and
        // elements after are >= it. We then keep entries strictly > threshold,
        // guaranteeing at most MAX_ENTRIES entries in the output.
        // f64 doesn't implement Ord (due to NaN), so we use partial_cmp.
        // The unwrap is safe: all values here are finite (>= 0.000005).
        let threshold = all_values
            .select_nth_unstable_by(n_to_drop - 1, |a, b| a.partial_cmp(b).unwrap())
            .1;
        eprintln!(
            "Note: {} entries survive the static threshold (>= 2^31); \
             applying dynamic threshold {:.6} to drop at least {} entries.",
            total_length, threshold, n_to_drop
        );
        *threshold
    } else {
        // No dynamic exclusion needed; use -infinity so the condition
        // `val > dynamic_threshold` is always satisfied for finite values.
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
                    .into_iter()
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
            "--writerowsums",
            "--out", &testdata("combined.r"),
        ]);

        run(args).unwrap();

        let actual = file_2_lines(&std::fs::read(testdata("combined.r.txt.gz")).unwrap(), true);
        let expected = file_2_lines(&std::fs::read(testdata("expected.rust.txt.gz")).unwrap(), true);

        assert_eq!(actual, expected);

        let actual_rowsums = file_2_lines(&std::fs::read(testdata("combined.r.rowsums")).unwrap(), false);

        assert_eq!(get_nrows(actual), actual_rowsums.len());

        std::fs::remove_file(testdata("combined.r.txt.gz")).unwrap();
        std::fs::remove_file(testdata("combined.r.rowsums")).unwrap();
    }

    #[test]
    #[rustfmt::skip]
    fn test_combine_matches_expected_r_10000() {
        let args = Args::parse_from([
            "combine",
            "--chunkpathsfile", &testdata("chunklength.files.txt"),
            "--snpcountsfile", &testdata("nsnps.txt"),
            "--maplengthsfile", &testdata("map_lengths.txt"),
            "--restrictrows=10000",
            "--out", &testdata("combined.r.10000"),
        ]);

        run(args).unwrap();

        let actual = file_2_lines(&std::fs::read(testdata("combined.r.10000.txt.gz")).unwrap(), true);
        assert!(actual.len() <= 10000);

        let actual_rowsums = file_2_lines(&std::fs::read(testdata("combined.r.10000.rowsums")).unwrap(), false);

        assert_eq!(get_nrows(actual), actual_rowsums.len());

        std::fs::remove_file(testdata("combined.r.10000.txt.gz")).unwrap();
        std::fs::remove_file(testdata("combined.r.10000.rowsums")).unwrap();
    }

    #[test]
    #[rustfmt::skip]
    fn test_combine_matches_expected_cpp() {
        let args = Args::parse_from([
            "combine",
            "--chunkpathsfile", &testdata("chunklength.files.txt"),
            "--snpcountsfile", &testdata("nsnps.txt"),
            "--maplengthsfile", &testdata("map_lengths.txt"),
            "--writerowsums",
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

        let actual_rowsums = file_2_lines(&std::fs::read(testdata("combined.cpp.rowsums")).unwrap(), false);
        let expected_rowsums = file_2_lines(&std::fs::read(testdata("expected.cpp.sorted.rowsums")).unwrap(), false);

        assert_eq!(actual_rowsums.len(), expected_rowsums.len());

        // there are rounding differences due to floating point precision between the rust implementation and the cpp implementation
        let actual_rowsums_rounded = actual_rowsums.iter().map(|l| l.parse::<f64>().unwrap().round()).collect::<Vec<_>>();
        let expected_rowsums_rounded = expected_rowsums.iter().map(|l| l.parse::<f64>().unwrap().round()).collect::<Vec<_>>();

        assert_eq!(actual_rowsums_rounded, expected_rowsums_rounded);

        std::fs::remove_file(testdata("combined.cpp.txt.gz")).unwrap();
        std::fs::remove_file(testdata("combined.cpp.rowsums")).unwrap();
    }
}
