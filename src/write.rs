use gzp::{ZBuilder, deflate::Gzip};
use rayon::prelude::*;
use std::io::Write;

use crate::{HMat, HVec, Result};

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
    dynamic_threshold: f64,
) -> Result<()> {
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

    let rows: Vec<(&usize, &HVec)> = cl.m.iter().collect();

    let mut pending = Vec::new();
    for chunk in rows.chunks(read_mat_threads * ROWS_PER_THREAD) {
        let mut formatted = None;
        pool.in_place_scope(|s| {
            s.spawn(|_| {
                formatted = Some(
                    chunk
                        .par_iter()
                        .map(|(key, row)| format_row(**key, row, dynamic_threshold))
                        .collect::<Result<Vec<_>>>(),
                )
            });
            write_rows(&pending, &mut parz, &mut rowsum_writer)
        })?;
        pending = formatted.unwrap()?;
    }
    write_rows(&pending, &mut parz, &mut rowsum_writer)?;

    parz.finish()?;
    if let Some(rowsum_writer) = &mut rowsum_writer {
        rowsum_writer.flush()?;
    }

    Ok(())
}

// Writes formatted rows (in order) to the gzip writer, and the sums of non-empty rows to the
// rowsum writer.
fn write_rows<W: Write>(
    rows: &[(Vec<u8>, f64)],
    parz: &mut impl Write,
    rowsum_writer: &mut Option<W>,
) -> Result<()> {
    for (data, rowsum) in rows {
        parz.write_all(data)?;
        if let Some(rowsum_writer) = rowsum_writer
            && !data.is_empty()
        {
            writeln!(rowsum_writer, "{:.5}", rowsum)?;
        }
    }
    Ok(())
}

// Formats one row of the sparse matrix for writing, and computes its sum.
fn format_row(key: usize, row: &HVec, dynamic_threshold: f64) -> Result<(Vec<u8>, f64)> {
    let mut buffer = Vec::new();
    let mut rowsum = 0.0;

    for (j, &val) in &row.v {
        if val >= 0.000005 {
            rowsum += val;
            if val > dynamic_threshold {
                writeln!(buffer, "{} {} {:.5}", key + 1, j + 1, val)?;
            }
        }
    }

    Ok((buffer, rowsum))
}
