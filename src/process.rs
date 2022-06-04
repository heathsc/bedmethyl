use std::{
	fs,
	collections::VecDeque,
};

use crossbeam_channel::bounded;
use crossbeam_utils::thread;
use r_htslib::*;

mod writer;
mod reader;
mod print_value;
mod distance;
mod hist;
mod heatmaps;
mod smooth;

use crate::config::*;
use crate::process::reader::MultiRecord;

const BLOCK_SIZE:usize = 1024;

pub type MsgBlock = VecDeque<MultiRecord>;

pub fn process(cfg: &Config, mut hts_vec: Vec<Hts>) -> anyhow::Result<()> {

	// Create output directory if necessary
	if let Some(dir) = cfg.dir() {
		fs::create_dir_all(dir)?
	}

	// Create thread scope so that we can share references across threads
	thread::scope(|scope| {

		// Spawn writer / processing thread if required
		let (out_tx, rx) = bounded(64);
		let writer_task = scope.spawn(move |_| { writer::writer_thread(cfg, rx) });

		// Spawn reader threads. Partition input files across reader threads
		// Each thread gets a vector of Hts structures and a channel
		// to send the input records
		let nthr = cfg.threads().min(cfg.n_samples());
		let mut channels = Vec::new();
		let mut read_tasks = Vec::new();
		let mut sample_idx = 0;

		for ix in 0..nthr {
			let ns = (hts_vec.len() / (nthr - ix)).max(1);
			let hvec: Vec<_> = hts_vec.drain(..ns).collect();
			let tx = if nthr > 1 {
				let (tx, rx) = bounded(64);
				channels.push((rx, ns));
				tx
			} else {
				// If we only have one thread we skip the merge step and send the output from the read thread directly to the next stage
				out_tx.clone()
			};
			read_tasks.push(scope.spawn(move |_| { reader::reader_thread(cfg, hvec, sample_idx, tx) }));
			sample_idx += ns;
		}

		// If multiple read threads are used, set up merge thread
		let mut merge_task = if nthr > 1 {
			Some(
				scope.spawn(move |_| { reader::merge_thread(cfg, channels, out_tx)})
			)
		} else {
			drop(out_tx);
			None
		};

		// Wait for read threads
		for jh in read_tasks.drain(..) { let _ = jh.join(); }

		// Wait for merge thread if used
		if let Some(jh) = merge_task.take() { let _ = jh.join(); }

		// Wait for writer thread
		let _ = writer_task.join();

	}).expect("Error creating thread scope");

	Ok(())

}