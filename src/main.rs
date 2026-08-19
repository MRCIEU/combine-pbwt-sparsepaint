use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use std::collections::HashMap;
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
    v: HashMap<usize, f64>,
    x0: f64, // default value for missing entries
}

impl HVec {
    fn new(x0: f64) -> Self {
        HVec {
            v: HashMap::new(),
            x0,
        }
    }

    fn get(&self, p: usize) -> f64 {
        *self.v.get(&p).unwrap_or(&self.x0)
    }

    fn set(&mut self, p: usize, val: f64) {
        self.v.insert(p, val);
    }

    fn add(&mut self, p: usize, delta: f64) {
        let current = self.get(p);
        self.v.insert(p, current + delta);
    }
}

struct HMat {
    m: Vec<HVec>,
}

impl HMat {
    fn new(nrows: usize, x0: f64) -> Self {
        let m = (0..nrows).map(|_| HVec::new(x0)).collect();
        HMat { m }
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
    let file = File::open(filename).unwrap_or_else(|e| {
        eprintln!("Error: could not open file {}: {}", filename, e);
        std::process::exit(1);
    });
    let reader = BufReader::new(file);
    let mut values = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap();
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let value: T = trimmed.parse().unwrap_or_else(|e| {
            eprintln!(
                "Error: could not parse '{}' from {}: {:?}",
                trimmed, filename, e
            );
            std::process::exit(1);
        });
        values.push(value);
    }
    values
}

/// Read the first per-chromosome painting file, setting entries to value * weight.
/// (Equivalent to C++ readdatafirst: uses set, not accumulate.)
fn read_data_first(filename: &str, cl: &mut HMat, weight: f64) {
    let file = File::open(filename).unwrap_or_else(|e| {
        eprintln!("Error: could not open {}: {}", filename, e);
        std::process::exit(1);
    });
    let reader = BufReader::new(GzDecoder::new(file));

    let mut q = 1usize;
    for line in reader.lines() {
        let line = line.unwrap();
        let mut parts = line.split_ascii_whitespace();
        let ind1: usize = parts.next().unwrap().parse().unwrap();
        let ind2: usize = parts.next().unwrap().parse().unwrap();
        let value: f64 = parts.next().unwrap().parse().unwrap();

        if ind1 == q {
            if q % 1000 == 0 {
                println!("ind{}", q);
            }
            q += 1;
        }

        cl.m[ind1 - 1].set(ind2 - 1, value * weight);
    }
}

/// Read a subsequent per-chromosome painting file, accumulating value * weight
/// into existing entries.
/// (Equivalent to C++ readdata: reads existing value and adds.)
fn read_data(filename: &str, cl: &mut HMat, weight: f64) {
    let file = File::open(filename).unwrap_or_else(|e| {
        eprintln!("Error: could not open {}: {}", filename, e);
        std::process::exit(1);
    });
    let reader = BufReader::new(GzDecoder::new(file));

    let mut q = 1usize;
    for line in reader.lines() {
        let line = line.unwrap();
        let mut parts = line.split_ascii_whitespace();
        let ind1: usize = parts.next().unwrap().parse().unwrap();
        let ind2: usize = parts.next().unwrap().parse().unwrap();
        let value: f64 = parts.next().unwrap().parse().unwrap();

        if ind1 == q {
            if q % 1000 == 0 {
                println!("ind{}", q);
            }
            q += 1;
        }

        cl.m[ind1 - 1].add(ind2 - 1, value * weight);
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

/// Maximum number of entries permitted in the output matrix (2^31 - 1).
/// This is a downstream constraint: indices must fit in a signed 32-bit integer.
const MAX_ENTRIES: usize = (1usize << 31) - 1;

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

    // Per-chromosome weight: genetic map length / SNP count, rounded to 6 d.p.
    let weights: Vec<f64> = total_gd
        .iter()
        .zip(total_snp.iter())
        .map(|(gd, snp)| (gd / snp * 1e6).round() / 1e6)
        .collect();

    // Accumulate weighted chunk lengths across all chromosomes.
    let mut cl = HMat::new(nind, 0.0);
    read_data_first(&chr_filenames[0], &mut cl, weights[0]);
    for i in 1..chr_filenames.len() {
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
        for (&_j, &val) in &cl.m[i].v {
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
        // select_nth_unstable(k) rearranges all_values so that the element at
        // sorted position k is in place. Elements before it are <= it and
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
    //   val >= 0.000005   (static threshold, always applied)
    //   val > dynamic_threshold  (dynamic threshold, only active when >= 2^31 entries)
    // -----------------------------------------------------------------------
    let out = File::create(outfile).unwrap_or_else(|e| {
        eprintln!("Error: could not create output file {}: {}", outfile, e);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(GzEncoder::new(out, Compression::default()));

    for i in 0..nind {
        for (&j, &val) in &cl.m[i].v {
            if val >= 0.000005 && val > dynamic_threshold {
                writeln!(writer, "{} {} {:.5}", i + 1, j + 1, val).unwrap();
            }
        }
    }

    // -----------------------------------------------------------------------
    // Write row sums file.
    //
    // Contains one value per line (individual 1 on line 1, etc.) representing
    // the sum of all entries in that row that survived the static threshold,
    // before any dynamic exclusion. Downstream code should use these as the
    // normalisation denominator rather than recomputing from the (possibly
    // truncated) output matrix.
    // -----------------------------------------------------------------------
    let rowsums_path = format!("{}.rowsums", outfile);
    let rowsums_file = File::create(&rowsums_path).unwrap_or_else(|e| {
        eprintln!(
            "Error: could not create row sums file {}: {}",
            rowsums_path, e
        );
        std::process::exit(1);
    });
    let mut rowsums_writer = BufWriter::new(rowsums_file);
    for &rs in &row_sums {
        writeln!(rowsums_writer, "{:.5}", rs).unwrap();
    }
}
