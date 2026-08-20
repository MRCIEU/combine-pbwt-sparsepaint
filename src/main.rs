use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use std::collections::BTreeMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

// ---------------------------------------------------------------------------
// Sparse matrix
//
// HVec is a sparse vector: only non-default entries are stored in a HashMap.
// HMat is a row-major sparse matrix: a Vec of HVec rows.
// ---------------------------------------------------------------------------

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

struct HMat {
    m: BTreeMap<usize, HVec>,
}

impl HMat {
    fn new() -> Self {
        HMat { m: BTreeMap::new() }
    }

    fn get(&mut self, row: usize) -> &mut HVec {
        self.m.entry(row).or_insert_with(|| HVec::new())
    }
}

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

/// Read one whitespace-trimmed token per non-blank line into a Vec<T>.
fn read_values_per_line<T>(filename: &str) -> Vec<T>
where
    T: std::str::FromStr,
    T::Err: std::fmt::Debug,
{
    let file =
        File::open(filename).unwrap_or_else(|e| panic!("could not open file {}: {}", filename, e));
    let reader = BufReader::new(file);
    let mut values = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap();
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let value: T = trimmed
            .parse()
            .unwrap_or_else(|e| panic!("could not parse '{}' from {}: {:?}", trimmed, filename, e));
        values.push(value);
    }
    values
}

/// Read a per-chromosome painting file, accumulating value * weight
/// into existing entries or creating new ones.
/// (Equivalent to C++ readdata: reads existing value and adds.)
fn read_data(filename: &str, cl: &mut HMat, weight: f64) {
    let file =
        File::open(filename).unwrap_or_else(|e| panic!("could not open {}: {}", filename, e));
    let reader = BufReader::new(GzDecoder::new(file));

    // let mut q = 1usize;
    for line in reader.lines() {
        let line = line.unwrap();
        let mut parts = line.split_ascii_whitespace();
        let ind1: usize = parts.next().unwrap().parse().unwrap();
        let ind2: usize = parts.next().unwrap().parse().unwrap();
        let value: f64 = parts.next().unwrap().parse().unwrap();

        // if ind1 == q {
        //     if q % 1000 == 0 {
        //         println!("ind{}", q);
        //     }
        //     q += 1;
        // }

        cl.get(ind1 - 1).add(ind2 - 1, value * weight);
    }
}

// ---------------------------------------------------------------------------
// Core logic
// ---------------------------------------------------------------------------

/// Maximum number of entries permitted in the output matrix (2^31 - 1).
/// This is a downstream constraint: indices must fit in a signed 32-bit integer.
const MAX_ENTRIES: usize = (1usize << 31) - 1;

/// Accumulate per-chromosome chunk-length painting files into a single
/// genome-wide weighted matrix and write the result to `out` (gzipped).
///
/// Per-row sums (computed before any dynamic exclusion) are written to
/// `rowsums_out` so downstream consumers can use the correct normalisation
/// denominator even when the dynamic threshold is active.
fn run<W: Write>(
    chr_filenames: &[String],
    total_snp: &[f64],
    total_gd: &[f64],
    nind: usize,
    out: W,
    rowsums_out: W,
) {
    // Per-chromosome weight: genetic map length / SNP count, rounded to 6 d.p.
    let weights: Vec<f64> = total_gd
        .iter()
        .zip(total_snp.iter())
        .map(|(gd, snp)| (gd / snp * 1e6).round() / 1e6)
        .collect();

    // Accumulate weighted chunk lengths across all chromosomes.
    let mut cl = HMat::new();
    for i in 0..chr_filenames.len() {
        println!("Processing chromosome {}", i + 1);
        read_data(&chr_filenames[i], &mut cl, weights[i]);
    }

    // -----------------------------------------------------------------------
    // Pass 1: collect all values that survive the static threshold (>= 0.000005)
    // and compute per-row sums before any dynamic exclusion.
    //
    // The row sums represent the true total weighted copying for each individual
    // (subject only to the static floor). They are written to a separate file so
    // that downstream consumers can use them as the correct normalisation
    // denominator even when the dynamic threshold further reduces the matrix.
    // -----------------------------------------------------------------------
    let mut all_values: Vec<f64> = Vec::new();
    let mut row_sums: Vec<f64> = vec![0.0; nind];

    for i in 0..nind {
        for (&_j, &val) in &cl.get(i).v {
            if val >= 0.000005 {
                row_sums[i] += val;
                all_values.push(val);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Determine dynamic threshold.
    //
    // If the number of entries surviving the static threshold is >= 2^31, we
    // must drop additional entries (the smallest ones) until fewer than 2^31
    // remain. We find the exact cut-off value using a linear-time selection
    // algorithm, then keep only entries strictly above it.
    //
    // Note: if there are ties at the threshold value, more than the minimum
    // number of entries may be dropped, but the count is guaranteed < 2^31.
    // -----------------------------------------------------------------------
    let dynamic_threshold: f64 = if all_values.len() >= (1usize << 31) {
        let n_to_drop = all_values.len() - MAX_ENTRIES;
        // select_nth_unstable_by(k) rearranges all_values so that the element
        // at sorted position k is in place. Elements before it are <= it and
        // elements after are >= it. We then keep entries strictly > threshold,
        // guaranteeing at most MAX_ENTRIES entries in the output.
        // f64 doesn't implement Ord (due to NaN), so we use partial_cmp.
        // The unwrap is safe: all values here are finite (>= 0.000005).
        let threshold = *all_values
            .select_nth_unstable_by(n_to_drop - 1, |a, b| a.partial_cmp(b).unwrap())
            .1;
        eprintln!(
            "Note: {} entries survive the static threshold (>= 2^31); \
             applying dynamic threshold {:.6} to drop at least {} entries.",
            all_values.len(),
            threshold,
            n_to_drop
        );
        threshold
    } else {
        // No dynamic exclusion needed; use -infinity so the condition
        // `val > dynamic_threshold` is always satisfied for finite values.
        f64::NEG_INFINITY
    };

    // -----------------------------------------------------------------------
    // Pass 2: write the output matrix, applying both thresholds.
    //
    // An entry is written iff:
    //   val >= 0.000005          (static threshold, always applied)
    //   val > dynamic_threshold  (only active when total entries >= 2^31)
    // -----------------------------------------------------------------------
    let mut writer = BufWriter::new(GzEncoder::new(out, Compression::default()));

    // for i in 0..nind {
    //     // TODO: sort by indices to print in order?
    //     for (&j, &val) in &cl.get(i).v {
    //         if val >= 0.000005 && val > dynamic_threshold {
    //             writeln!(writer, "{} {} {:.5}", i + 1, j + 1, val).unwrap();
    //         }
    //     }
    // }

    for (i, row) in cl.m {
        for (j, val) in row.v {
            if val >= 0.000005 && val > dynamic_threshold {
                writeln!(writer, "{} {} {:.5}", i + 1, j + 1, val).unwrap();
            }
        }
    }

    // -----------------------------------------------------------------------
    // Write row sums.
    //
    // Contains one value per line (individual 1 on line 1, etc.) representing
    // the sum of all entries in that row that survived the static threshold,
    // before any dynamic exclusion. Downstream code should use these as the
    // normalisation denominator rather than recomputing from the (possibly
    // truncated) output matrix.
    // -----------------------------------------------------------------------
    let mut rowsums_writer = BufWriter::new(rowsums_out);
    for &rs in &row_sums {
        writeln!(rowsums_writer, "{:.5}", rs).unwrap();
    }
}

// ---------------------------------------------------------------------------
// CLI entry point
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.iter().any(|a| a == "-h" || a == "--help") {
        println!(
            "Usage: {} <chr_chunks_filenames_list> <chr_snp_counts_list> <chr_map_lengths_list> <nind> <outfile>",
            args[0]
        );
        println!("  Row sums (pre-dynamic-exclusion) are written to <outfile>.rowsums");
        return;
    }

    if args.len() != 6 {
        eprintln!(
            "Usage: {} <chr_chunks_filenames_list> <chr_snp_counts_list> <chr_map_lengths_list> <nind> <outfile>",
            args[0]
        );
        std::process::exit(1);
    }

    let filenames_list = &args[1];
    let snp_counts_list = &args[2];
    let map_lengths_list = &args[3];
    let nind: usize = args[4].parse().unwrap_or_else(|_| {
        eprintln!("Error: <nind> must be a positive integer");
        std::process::exit(1);
    });
    let outfile = &args[5];

    let chr_filenames = read_values_per_line::<String>(filenames_list);
    let total_snp = read_values_per_line::<f64>(snp_counts_list);
    let total_gd = read_values_per_line::<f64>(map_lengths_list);

    if chr_filenames.len() != total_snp.len() || chr_filenames.len() != total_gd.len() {
        eprintln!(
            "Error: {}, {} and {} must all have the same number of lines (one per chromosome)",
            filenames_list, snp_counts_list, map_lengths_list
        );
        std::process::exit(1);
    }

    let out = File::create(outfile).unwrap_or_else(|e| {
        eprintln!("Error: could not create {}: {}", outfile, e);
        std::process::exit(1)
    });

    let rowsums_path = format!("{}.rowsums", outfile);
    let rowsums_out = File::create(&rowsums_path).unwrap_or_else(|e| {
        eprintln!("Error: could not create {}: {}", rowsums_path, e);
        std::process::exit(1)
    });

    run(
        &chr_filenames,
        &total_snp,
        &total_gd,
        nind,
        out,
        rowsums_out,
    );
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;
    use std::path::PathBuf;

    fn testdata(filename: &str) -> String {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("testdata")
            .join("n10000")
            .join(filename)
            .to_string_lossy()
            .into_owned()
    }

    /// Decompress a gzipped byte slice and return its lines sorted.
    /// Sorting is necessary because HashMap iteration order is non-deterministic.
    fn decode_gz(bytes: &[u8]) -> Vec<String> {
        let reader = BufReader::new(GzDecoder::new(bytes));
        let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();
        lines
    }

    fn load_f64_vec(filename: &str) -> Vec<f64> {
        let mut vec: Vec<f64> = Vec::new();
        let mut file = std::fs::File::open(filename).unwrap();
        let reader = std::io::BufReader::new(&mut file);
        for line in reader.lines() {
            let line = line.unwrap();
            vec.push(line.parse().unwrap());
        }
        vec
    }

    #[test]
    fn test_combine_matches_expected() {
        let chr_filenames = vec![
            testdata("chr01.chunklengths.s.out.gz"),
            testdata("chr02.chunklengths.s.out.gz"),
            testdata("chr03.chunklengths.s.out.gz"),
            testdata("chr04.chunklengths.s.out.gz"),
            testdata("chr05.chunklengths.s.out.gz"),
        ];
        let total_snp = load_f64_vec(&testdata("nsnps.txt"));
        let total_gd = load_f64_vec(&testdata("chr1-5.maplengths.txt"));
        let nind = 10000;

        let mut out_buf: Vec<u8> = Vec::new();
        let mut rowsums_buf: Vec<u8> = Vec::new();
        run(
            &chr_filenames,
            &total_snp,
            &total_gd,
            nind,
            &mut out_buf,
            &mut rowsums_buf,
        );

        let actual = decode_gz(&out_buf);
        let expected = decode_gz(&std::fs::read(testdata("combined.chunklengths.txt.gz")).unwrap());

        assert_eq!(actual, expected);
    }
}
