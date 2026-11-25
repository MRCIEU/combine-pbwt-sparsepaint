use gzp::{ZBuilder, deflate::Gzip};
use rayon::prelude::*;
use std::io::Write;

use crate::{HMat, HVec, Result, STATIC_THRESHOLD};

/// The first 8 bytes of a binary row stream (`--binary`), followed by the number of rows (u64).
pub(crate) const BINARY_MAGIC: &[u8; 8] = b"CMBROWS1";

/// One formatted row: its index, its text output, its binary output (if wanted) and its sum.
type FormattedRow = (usize, Vec<u8>, Option<Vec<u8>>, f64);

// The number of rows formatted in parallel per formatting thread before being written out.
// This bounds the amount of formatted (uncompressed) output held in memory at once.
const ROWS_PER_THREAD: usize = 16;

// Parallel sparse matrix reading + parallel gzip writing.
// Rows are formatted in parallel in chunks on a rayon thread pool. While one chunk is being
// formatted, the previous (already formatted, in-order) chunk is written by the calling thread
// to the multi-threaded gzip writer from the gzp crate.
pub(crate) fn parallel_read_write<BW: Write + Send + 'static>(
    nthreads: usize,
    cl: HMat,
    hap_writer: BW,
    mut rowsum_writer: Option<BW>,
    binary: Option<(BW, usize)>,
    dynamic_thresholds: Vec<f64>,
) -> Result<()> {
    // binary row stream: header now; rows (with empty rows for samples that have none) as they come
    let mut binary = match binary {
        Some((mut w, n_rows)) => {
            w.write_all(BINARY_MAGIC)?;
            w.write_all(&(n_rows as u64).to_le_bytes())?;
            Some(BinaryRows {
                w,
                next_row: 0,
                n_rows,
            })
        }
        None => None,
    };
    let want_binary = binary.is_some();
    assert_eq!(
        dynamic_thresholds.len(),
        cl.m.len(),
        "one threshold per row"
    );

    // Threads are split between gzip and formatting, with the gzip count rounded down.
    // This still works with nthreads = 1: (gzip_threads = 0, read_mat_threads = 1)
    let gzip_threads = nthreads / 2;
    let read_mat_threads = nthreads - gzip_threads;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(read_mat_threads)
        .build()?;

    let mut parz = ZBuilder::<Gzip, _>::new()
        .num_threads(gzip_threads)
        .from_writer(hap_writer);

    let rows: Vec<(&usize, &HVec, f64)> =
        cl.m.iter()
            .zip(dynamic_thresholds)
            .map(|((key, row), threshold)| (key, row, threshold))
            .collect();

    let mut pending = Vec::new();
    for chunk in rows.chunks(read_mat_threads * ROWS_PER_THREAD) {
        let mut formatted = None;
        pool.in_place_scope(|s| {
            s.spawn(|_| {
                formatted = Some(
                    chunk
                        .par_iter()
                        .map(|(key, row, threshold)| {
                            format_row(**key, row, *threshold, want_binary)
                        })
                        .collect::<Result<Vec<_>>>(),
                )
            });
            write_rows(&pending, &mut parz, &mut rowsum_writer, &mut binary)
        })?;
        pending = formatted.unwrap()?;
    }
    write_rows(&pending, &mut parz, &mut rowsum_writer, &mut binary)?;

    parz.finish()?;
    if let Some(rowsum_writer) = &mut rowsum_writer {
        rowsum_writer.flush()?;
    }
    if let Some(b) = &mut binary {
        b.finish()?;
    }

    Ok(())
}

/// The binary row stream's writer: rows must be written for every index 0..n_rows, in order, so
/// rows absent from the matrix are written as empty rows.
struct BinaryRows<W: Write> {
    w: W,
    next_row: usize,
    n_rows: usize,
}

impl<W: Write> BinaryRows<W> {
    fn empty_rows_until(&mut self, row: usize) -> Result<()> {
        while self.next_row < row {
            self.w.write_all(&0u32.to_le_bytes())?;
            self.next_row += 1;
        }
        Ok(())
    }

    fn write_row(&mut self, row: usize, data: &[u8]) -> Result<()> {
        self.empty_rows_until(row)?;
        self.w.write_all(data)?;
        self.next_row += 1;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.empty_rows_until(self.n_rows)?;
        self.w.flush()?;
        Ok(())
    }
}

// Writes formatted rows (in order) to the gzip writer, "index sum" for every row to the rowsum
// writer (including rows with no entries written because of --restrict), and the binary rows.
fn write_rows<W: Write>(
    rows: &[FormattedRow],
    parz: &mut impl Write,
    rowsum_writer: &mut Option<W>,
    binary: &mut Option<BinaryRows<W>>,
) -> Result<()> {
    for (key, data, bin, rowsum) in rows {
        parz.write_all(data)?;
        if let Some(rowsum_writer) = rowsum_writer {
            writeln!(rowsum_writer, "{} {:.5}", key + 1, rowsum)?;
        }
        if let (Some(b), Some(bin)) = (binary.as_mut(), bin) {
            b.write_row(*key, bin)?;
        }
    }
    Ok(())
}

// Formats one row of the sparse matrix for writing (text, and binary if wanted), and computes its
// sum over all values that pass the static threshold (whether or not they pass the dynamic
// threshold). A binary row is: the number of entries n (u32), then n column indices (u32, 0-based,
// ascending), then n values (f32), all little-endian.
fn format_row(
    key: usize,
    row: &HVec,
    dynamic_threshold: f64,
    binary: bool,
) -> Result<FormattedRow> {
    let mut buffer = Vec::new();
    let mut rowsum = 0.0;
    let mut cols: Vec<u32> = Vec::new();
    let mut vals: Vec<f32> = Vec::new();

    for (j, &val) in &row.v {
        if val >= STATIC_THRESHOLD {
            rowsum += val;
            if val > dynamic_threshold {
                writeln!(buffer, "{} {} {:.5}", key + 1, j + 1, val)?;
                if binary {
                    cols.push(*j as u32);
                    vals.push(val as f32);
                }
            }
        }
    }

    let bin = binary.then(|| {
        let mut b = Vec::with_capacity(4 + 8 * cols.len());
        b.extend_from_slice(&(cols.len() as u32).to_le_bytes());
        cols.iter()
            .for_each(|c| b.extend_from_slice(&c.to_le_bytes()));
        vals.iter()
            .for_each(|v| b.extend_from_slice(&v.to_le_bytes()));
        b
    });

    Ok((key, buffer, bin, rowsum))
}
