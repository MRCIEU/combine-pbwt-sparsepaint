use gzp::{ZBuilder, ZWriter, deflate::Gzip};
use std::collections::HashMap;
use std::io::Write;
use std::sync::{Arc, RwLock};
use std::thread;

use crossbeam_channel::{Receiver, Sender, bounded};
use crossbeam_utils::sync::WaitGroup;

use crate::HMat;

#[derive(Clone)]
struct RowLocation {
    idx: usize, // index for ordering
    key: usize, // the actual entry in the map (e.g. the sample's ID)
}

struct WriteData {
    idx: usize,
    data: Vec<u8>,
}

pub fn parallel_read_write<BW: Write + Send + 'static>(
    nthreads: usize,
    cl: HMat,
    writer: BW,
    dynamic_threshold: f64,
) {
    let (txrl, rxrl) = bounded::<RowLocation>(nthreads);
    let (txwd, rxwd) = bounded::<WriteData>(nthreads);

    let writer_thread = thread::spawn({
        move || {
            // gather_write(&rxwd, writer);
            gather_parallel_write(&rxwd, writer);
        }
    });

    let keys = cl.m.keys().cloned().collect::<Vec<usize>>();

    let keys_thread = thread::spawn({
        move || {
            let mut counter = 0;
            for k in keys {
                let rl = RowLocation {
                    idx: counter,
                    key: k,
                };
                txrl.send(rl).unwrap();
                counter += 1;
            }
        }
    });

    let mut workers = Vec::new();
    for _ in 0..nthreads {
        workers.push((rxrl.clone(), txwd.clone()))
    }

    let waitgroup = WaitGroup::new();
    let arc = Arc::new(RwLock::new(cl));
    for (rxrl, txwd) in workers {
        let wg = waitgroup.clone();
        let arc = arc.clone();
        thread::spawn(move || {
            get_row(&rxrl, &txwd, arc, dynamic_threshold);
            drop(wg);
        });
    }

    keys_thread.join().unwrap();

    waitgroup.wait();
    drop(txwd);

    writer_thread.join().unwrap();
}

fn get_row(
    rxrl: &Receiver<RowLocation>,
    txwd: &Sender<WriteData>,
    cl: Arc<RwLock<HMat>>,
    dynamic_threshold: f64,
) {
    for rl in rxrl.iter() {
        let row = cl.read().unwrap().m.get(&rl.key).unwrap().v.clone();
        let mut buffer = Vec::new();
        for (j, val) in row {
            if val >= 0.000005 && val > dynamic_threshold {
                writeln!(buffer, "{} {} {:.5}", rl.key + 1, j + 1, val).unwrap();
            }
        }
        let wd = WriteData {
            idx: rl.idx,
            data: buffer,
        };
        txwd.send(wd).unwrap();
    }
}

fn gather_write<W: Write>(rxwd: &Receiver<WriteData>, mut writer: W) {
    let mut m: HashMap<usize, WriteData> = HashMap::new();

    let mut counter: usize = 0;

    for r in rxwd.iter() {
        m.insert(r.idx, r);
        while m.contains_key(&counter) {
            let rv = m.remove(&counter).unwrap();
            writer.write_all(&rv.data).unwrap();
            counter += 1;
        }
    }

    _ = writer.flush();
}

fn gather_parallel_write<W: Write + Send + 'static>(rxwd: &Receiver<WriteData>, writer: W) {
    let mut m: HashMap<usize, WriteData> = HashMap::new();

    let mut counter: usize = 0;

    let mut parz = ZBuilder::<Gzip, _>::new()
        .num_threads(8)
        .from_writer(writer);

    for r in rxwd.iter() {
        m.insert(r.idx, r);
        while m.contains_key(&counter) {
            let rv = m.remove(&counter).unwrap();
            parz.write_all(&rv.data).unwrap();
            counter += 1;
        }
    }

    _ = parz.finish().unwrap();
}
