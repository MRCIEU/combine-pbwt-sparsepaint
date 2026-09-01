use crossbeam_channel::{Receiver, Sender, bounded};
use crossbeam_utils::sync::WaitGroup;
use gzp::{ZBuilder, deflate::Gzip};
use std::collections::HashMap;
use std::io::Write;
use std::sync::{Arc, RwLock};
use std::thread;

use crate::{CombineError, HMat, Result};

// The location of one row in the sparse matrix.
// Used to communicate the location of a row to a formatting thread.
pub(crate) struct RowLocation {
    idx: usize, // index for ordering
    key: usize, // the actual entry in the map (e.g. the sample's ID)
}

// One row of the sparse matrix, formatted for writing.
// Used to communicate a row of formatted data to a gzip writer.
pub(crate) struct WriteData {
    idx: usize,
    data: Vec<u8>,
    rowsum: Option<f64>,
}

// Splits the number of threads into gzip and read matrix threads.
// For odd numbers of threads, the gzip thread count is rounded down.
// This still works with nthreads = 1: (gzip_threads = 0, read_mat_threads = 1)
fn split_n_threads(nthreads: usize) -> (usize, usize) {
    let gzip_threads = nthreads / 2;
    (gzip_threads, nthreads - gzip_threads)
}

// Top level logic for parallel sparse matrix reading + parallel gzip writing (using multiple threads).
pub(crate) fn parallel_read_write<BW: Write + Send + 'static>(
    nthreads: usize,
    cl: HMat,
    hap_writer: BW,
    rowsum_writer: Option<BW>,
    dynamic_threshold: f64,
    write_rowsums: bool,
) -> Result<()> {
    let (txrl, rxrl) = bounded::<RowLocation>(nthreads);
    let (txwd, rxwd) = bounded::<WriteData>(nthreads);

    let (gzip_threads, read_mat_threads) = split_n_threads(nthreads);

    let writer_thread = thread::spawn(move || {
        gather_parallel_write(&rxwd, hap_writer, rowsum_writer, gzip_threads)
    });

    let keys = cl.m.keys().cloned().collect::<Vec<usize>>();

    let keys_thread = thread::spawn({
        move || -> Result<()> {
            for (counter, k) in keys.into_iter().enumerate() {
                let rl = RowLocation {
                    idx: counter,
                    key: k,
                };
                txrl.send(rl)?;
            }
            Ok(())
        }
    });

    let mut workers = Vec::new();
    for _ in 0..read_mat_threads {
        workers.push((rxrl.clone(), txwd.clone()))
    }

    let waitgroup = WaitGroup::new();
    let arc = Arc::new(RwLock::new(cl));
    let handles: Vec<_> = workers
        .into_iter()
        .map(|(rxrl, txwd)| {
            let wg = waitgroup.clone();
            let arc = arc.clone();
            thread::spawn(move || {
                let result = get_row(&rxrl, &txwd, arc, dynamic_threshold, write_rowsums);
                drop(wg);
                result
            })
        })
        .collect();

    keys_thread.join()??;

    waitgroup.wait();
    drop(txwd);

    writer_thread.join()??;

    for handle in handles {
        handle.join()??;
    }

    Ok(())
}

// Worker thread that reads rows from the sparse matrix and writes them (formatted) to a channel,
// to be passed to the (parallel) gzip writer
fn get_row(
    rxrl: &Receiver<RowLocation>,
    txwd: &Sender<WriteData>,
    cl: Arc<RwLock<HMat>>,
    dynamic_threshold: f64,
    write_rowsums: bool,
) -> Result<()> {
    for rl in rxrl.iter() {
        let row = cl.read()?.m.get(&rl.key).unwrap().v.clone();
        let mut buffer = Vec::new();
        let mut write_this_rowsum = write_rowsums;
        let mut rowsum = 0.0;

        for (j, val) in row {
            if val >= 0.000005 {
                if write_rowsums {
                    rowsum += val;
                }
                if val > dynamic_threshold {
                    writeln!(buffer, "{} {} {:.5}", rl.key + 1, j + 1, val)?;
                }
            }
        }
        if buffer.is_empty() {
            write_this_rowsum = false;
        }
        let wd = WriteData {
            idx: rl.idx,
            data: buffer,
            rowsum: if write_this_rowsum {
                Some(rowsum)
            } else {
                None
            },
        };
        txwd.send(wd)?;
    }
    Ok(())
}

// Worker thread that gathers rows from a channel and passes them in order to the gzip writer.
// The gzip writer built here is a multi-threaded gzip writer from the gzp crate
fn gather_parallel_write<W: Write + Send + 'static>(
    rxwd: &Receiver<WriteData>,
    hap_writer: W,
    mut rowsum_writer: Option<W>,
    nthreads: usize,
) -> Result<()> {
    let mut m: HashMap<usize, WriteData> = HashMap::new();

    let mut counter: usize = 0;

    let mut parz = ZBuilder::<Gzip, _>::new()
        .num_threads(nthreads)
        .from_writer(hap_writer);

    for r in rxwd.iter() {
        m.insert(r.idx, r);
        while m.contains_key(&counter) {
            if let Some(rv) = m.remove(&counter) {
                parz.write_all(&rv.data)?;
                if let Some(rowsum_writer) = &mut rowsum_writer
                    && let Some(rowsum) = rv.rowsum
                {
                    writeln!(rowsum_writer, "{:.5}", rowsum)?;
                }
                counter += 1;
            } else {
                return Err(CombineError::Message(
                    "couldn't find index in writing map".to_string(),
                ));
            }
        }
    }

    parz.finish()?;
    if let Some(rowsum_writer) = &mut rowsum_writer {
        rowsum_writer.flush()?;
    }

    Ok(())
}
