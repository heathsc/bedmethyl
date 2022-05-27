use std::{
	sync::Arc,
	io::{self, Write},
	ops::Deref,
	collections::HashMap,
	thread::{self, JoinHandle},
};

use crossbeam_channel::{bounded, Sender, Receiver};

use crate::config::Config;
use crate::sample::Sample;
use super::hist::LnP;

enum Dist<'a> {
	DistPtr(&'a [f64]),
	Dist(Vec<f64>),
}

impl<'a> Deref for Dist<'a> {
	type Target = [f64];

	fn deref(&self) -> &Self::Target {
		match self {
			Self::DistPtr(p) => p,
			Self::Dist(d) => d.as_slice(),
		}
	}
}

struct BetaCache {
	cache: Vec<Box<[f64]>>,
	min_n: u16,
	max_n: u16,
}

impl BetaCache {
	fn new(min_n: u16, max_n: u16, lnp: &LnP) -> Self {
		 
		let sumx = |x| x * (x + 1) / 2;
		
		let cache = if max_n >= min_n {
			let n_bins = lnp.n_bins();
			let size = sumx((max_n + 1) as usize) - sumx(min_n as usize);
			debug!("Setting up cache for n from {} to {}: {} elements", min_n, max_n, size);
			let mut v = Vec::with_capacity(size);
			
			for n in min_n..=max_n {
				for a in 0..=n {
					let b = n - a;
					let mut wk = Vec::with_capacity(n_bins);
					let sum = lnp.get_dist(a as u32, b as u32, &mut wk);
					for elem in wk.iter_mut() { *elem /= sum }
					v.push(wk.into_boxed_slice())
				}
			}
			debug!("Set up cache with {} elements", v.len());
			v
		} else {
			Vec::new()
		};
		
		Self{ cache, min_n, max_n }
	}
	
	fn get_dist(&self, lnp: &LnP, a: u32, b: u32) -> Dist {
		let n = a + b;
		if n > self.max_n as u32 {
			let mut dst = Vec::with_capacity(lnp.n_bins());
			let sum = lnp.get_dist(a, b, &mut dst);
			for elem in dst.iter_mut() { *elem /= sum }
			Dist::Dist(dst)
		} else {
			let sumx = |x| x * (x + 1) / 2;
			assert!(n >= self.min_n as u32);
			let ix = sumx(n as usize) - sumx(self.min_n as usize) + a as usize;
			Dist::DistPtr(&self.cache[ix])
		}
	}	
}

pub struct Obs {
	pub(super) idx: usize,
	pub(super) a: u32,
	pub(super) b: u32,
}

pub struct PairObs {
	pub(super) idx: usize,
	pub(super) n: usize,
	pub(super) cts: [u32; 4],
}

type VObs = Vec<Obs>;

struct HeatDist {
	dist: Box<[f64]>,
	n: usize,
}

impl HeatDist {
	fn new(n_bins: usize) -> Self {
		Self{dist: vec!(0.0; n_bins * n_bins).into_boxed_slice(), n: 0}
	}
}

const BLOCKSIZE: usize = 256;

type CacheHashKey = [u8; 4];
type CacheHash = HashMap<CacheHashKey, u32>;

struct HeatMapPair {
	cache: CacheHash,
	dist: HeatDist,
}

#[derive(Default)]
struct HeatMapTask {
	hash: HashMap<usize, HeatMapPair>,	
}

impl HeatMapTask {
	fn add_pair(&mut self, ix: usize, pair: HeatMapPair) {
		if self.hash.insert(ix, pair).is_some() {
			panic!("Pair with id {} inserted twice", ix);
		}
	}

	fn pair(&mut self, ix: usize) -> Option<&mut HeatMapPair> {
		self.hash.get_mut(&ix)
	}
}

pub(crate) struct HeatMaps {
	n_bins: usize,
	tx: Option<Sender<Vec<VObs>>>,
	task: Option<JoinHandle<Vec<HeatDist>>>,
	tvec: Option<Vec<VObs>>,
}

struct HeatMapCfg {
	n_bins: usize,
	max_n: u16,
	min_n: u16,
	full_limit: u32,
	threads: usize,
}

impl HeatMapCfg {
	fn new(cfg: &Config) -> Self {
		Self {
			n_bins: cfg.n_bins(),
			max_n: cfg.kde_cache_limit(),
			min_n: cfg.min_samples().min(u16::MAX as usize) as u16,
			full_limit: cfg.heatmaps_full_limit() as u32,
			threads: cfg.threads(),
		}
	}
}

impl HeatMaps {
	pub(crate) fn new(n_samples: usize, cfg: &Config) -> Self {
		assert!(n_samples > 0);
		let hm_cfg = HeatMapCfg::new(cfg);
		let nt = hm_cfg.threads;
		let n_bins = hm_cfg.n_bins;
		let (tx, rx) = bounded(nt * 4);
		let task = thread::spawn(move || { heatmap_thread(n_samples, hm_cfg, rx) });
		let tvec = Some(Vec::with_capacity(BLOCKSIZE));
		Self{n_bins, tx: Some(tx), task: Some(task), tvec}

	}
	
	pub(crate) fn add_vec_obs(&mut self, v: Vec<Obs>) {
		
		let mut tv = self.tvec.take().expect("add_vec_obs() called when tvec is None");
		tv.push(v);
		if tv.len() >= BLOCKSIZE {
			self.tx.as_mut().unwrap().send(tv).expect("Error sending observations for heatmap calculation");
			trace!("Sent message to heatmaps thread");
			tv = Vec::with_capacity(BLOCKSIZE);
		}	
		self.tvec = Some(tv)
	}
	
	// We take() tx so that afterwards it will be dropped and the channel will be closed
	fn flush_buffer(&mut self) {
		trace!("heathmaps: flushing buffer");
		let tx = self.tx.take().unwrap();
		if let Some(tv) = self.tvec.take() {
			tx.send(tv).expect("Error sending last observations for heatmap calculation");
		}		
	}
	
	pub(crate) fn print_table<W: Write>(&mut self, wrt: &mut W, samples: &[Sample], blank_lines: bool) -> io::Result<()> {
		
		trace!("heatmaps: printing table");
		// Flush buffers and close channel so signalling the child thread to exit.
		self.flush_buffer();
		
		// Wait for child to exit
		let task = self.task.take().unwrap();
		let mut hdist = task.join().expect("Error joining heatmaps thread");

		write!(wrt, "x\ty")?;
		for (i, s1) in samples[1..].iter().enumerate() {
			for s2 in &samples[0..=i] {
				write!(wrt, "\t{}_{}", s1.name(), s2.name())?;
			}
		}
		writeln!(wrt)?;

		let nb = self.n_bins as f64;

		for hd in hdist.iter_mut() {
			if hd.n > 0 {
				let konst = nb * nb / (hd.n as f64);
				for elem in hd.dist.iter_mut() { *elem *= konst }
			}
		}

		for ix in 0..self.n_bins {
			let x = ((ix as f64) + 0.5) / nb;
			for iy in 0..self.n_bins {
				let y = ((iy as f64) + 0.5) / nb;
				write!(wrt, "{}\t{}", x, y)?;
				for hd in hdist.iter() {
					write!(wrt, "\t{:.4e}", hd.dist[ix * self.n_bins + iy])?;
				}
				writeln!(wrt)?;
			}
			if blank_lines { writeln!(wrt)? }
		}
		Ok(())
	}
}

fn heatmap_thread(n_samples: usize, cfg: HeatMapCfg, recv: Receiver<Vec<VObs>>) -> Vec<HeatDist> {

	assert!(n_samples > 1);
	let n_bins = cfg.n_bins;
	let lnp = Arc::new(LnP::new(n_bins));
		
	let max_n = cfg.max_n;
	let min_n = cfg.min_n;
	let full_limit = cfg.full_limit as u32;
	let cache = Arc::new(BetaCache::new(min_n, max_n, &lnp));

	let ssize = n_samples * (n_samples - 1) / 2;
	let nt = cfg.threads;

	let ext_thr = if ssize < nt {
		let ext_th = nt;
		let (tx, rx) = bounded(ext_th * 64);
		let jobs: Vec<_> = (0..ext_th).map(|_| {
			let rx = rx.clone();
			let cache = cache.clone();
			let lnp = lnp.clone();
			thread::spawn(move || { heatmap_calc(ssize, cache, lnp, rx)})
		}).collect();
		Some((tx, jobs))
	} else {
		None
	};

	let nt = nt.min(ssize);
	let mut hm_tasks: Vec<_> = (0..nt).map(|_| HeatMapTask::default()).collect();
	let mut task_idx = Vec::with_capacity(ssize);
	for (i, j) in (0..ssize).zip((0..nt).cycle()) {
		let pair = HeatMapPair{
			cache: HashMap::new(),
			dist: HeatDist::new(n_bins),
		};
		task_idx.push(j);
		hm_tasks[j].add_pair(i, pair);
	}
	
	let mut senders = Vec::with_capacity(nt);
	let mut tasks: Vec<_> = hm_tasks.drain(..).enumerate().map(|(ix, hm)| {
		let (tx, rx) = bounded(256);
		let cache = cache.clone();
		let lnp = lnp.clone();
		senders.push(tx);
		let ext_tx = ext_thr.as_ref().map(|(t, _)| t.clone());
		thread::spawn(move || { heatmap_subtask(ix, full_limit, cache, lnp, hm, rx, ext_tx) })
	}).collect();

	let mut ext_thr = ext_thr.map(|(_, jobs)| jobs);

	let mut pobs_vec: Vec<_> = (0..nt).map(|_| Some(Vec::with_capacity(BLOCKSIZE))).collect();
	for rec in recv.iter() {
		for vrec in rec.iter() {
			for (i, obs1) in vrec[1..].iter().enumerate() {
				let ix = obs1.idx * (obs1.idx - 1) / 2;	
				for obs2 in vrec[0..=i].iter() {
					let idx = ix + obs2.idx;
					let task_id = task_idx[idx];
					pobs_vec[task_id].as_mut().unwrap().push(PairObs{idx, n: 1, cts: [obs1.a, obs1.b, obs2.a, obs2.b] });
					if pobs_vec[task_id].as_ref().unwrap().len() >= BLOCKSIZE {
						let pv = pobs_vec[task_id].take().unwrap();
						pobs_vec[task_id] = Some(Vec::with_capacity(BLOCKSIZE));
						senders[task_id].send(pv).expect("Error sending message to sub task");
					}
				}
			}			
		}
	}
	for (pv, tx) in pobs_vec.drain(..).zip(senders.drain(..)) {
		let pv = pv.unwrap();
		if !pv.is_empty() {
			tx.send(pv).expect("Error sending message to sub task");
		}
	}

	info!("End of input.  Processing cache");
	debug!("heatmap_thread(): waiting for child tasks to end");

	let mut hdist: Vec<Option<HeatDist>> = (0..ssize).map(|_| None).collect();
		
	for task in tasks.drain(..) {
		let mut hm = task.join().expect("Error waiting for heatmap subtask");
		for (ix, hmp) in hm.hash.drain() {
			let HeatMapPair{cache: _, dist} = hmp;
			hdist[ix] = Some(dist);
		}
	}
	let mut hdist: Vec<HeatDist> = hdist.drain(..).map(|hm| hm.expect("Missing element!")).collect();

	if let Some(mut tvec) = ext_thr.take() {
		debug!("Adding contributions from extra heatmap compute threads");
		for task in tvec.drain(..) {
			let mut vdist = task.join().expect("Error waiting for extra heatmap compute threads");
			for (dist, dist1) in vdist.drain(..).zip(hdist.iter_mut()) {
				dist1.n += dist.n;
				for (elem, elem1) in dist.dist.iter().zip(dist1.dist.iter_mut()) {
					*elem1 += *elem
				}
			}
		}
	}
	debug!("heatmap_thread(): closing down");
	hdist
}

fn heatmap_subtask(ix: usize, full_limit: u32, cache: Arc<BetaCache>, lnp: Arc<LnP>, mut hm: HeatMapTask, rx: Receiver<Vec<PairObs>>,mut tx: Option<Sender<Vec<PairObs>>>) -> HeatMapTask {

	let mut pvec = if tx.is_some() { Some(Vec::with_capacity(BLOCKSIZE)) } else { None };

	debug!("heatmap_subtask({}): waiting for input", ix);
	for mut rec in rx.iter() {
		for pobs in rec.drain(..) {
			let pair = hm.pair(pobs.idx).expect("Pair idx not found");
			if pobs.cts.iter().any(|x| *x > full_limit) {
				if let Some(v) = pvec.as_mut() {
					v.push(pobs);
					if v.len() >= BLOCKSIZE {
						let pv = pvec.take().unwrap();
						tx.as_mut().unwrap().send(pv).expect("Error sending message from subtask");
						pvec = Some(Vec::with_capacity(BLOCKSIZE));
					}
				} else {
					add_contrib(pobs.n, &pobs.cts, &cache, &lnp, &mut pair.dist)
				}
			} else {
				let cts = &pobs.cts;
				let key = [cts[0] as u8, cts[1] as u8, cts[2] as u8, cts[3] as u8];
				let elem = pair.cache.entry(key).or_insert(0);
				*elem += 1;
			}
		}
	}

	debug!("heatmap_subtask({}): end of input; processing caches", ix);

	for (i, hmp) in hm.hash.iter_mut() {
		let hdist = &mut hmp.dist;
		for(cts, n) in hmp.cache.iter() {
			let cts = [cts[0] as u32, cts[1] as u32, cts[2] as u32, cts[3] as u32];
			if let Some(v) = pvec.as_mut() {
				let pobs = PairObs{ idx: *i, n: *n as usize, cts };
				v.push(pobs);
				if v.len() >= BLOCKSIZE {
					let pv = pvec.take().unwrap();
					tx.as_mut().unwrap().send(pv).expect("Error sending message from subtask");
					pvec = Some(Vec::with_capacity(BLOCKSIZE));
				}
			} else {
				add_contrib(*n as usize, &cts, &cache, &lnp, hdist)
			}
/*
			hdist.n += *n as usize;
			let td1 = cache.get_dist(&lnp, cts[0] as u32, cts[1] as u32);
			let td2 = cache.get_dist(&lnp, cts[2] as u32, cts[3] as u32);
			let mut ix = 0;
			let n = *n as f64;
			for x in td1.iter()  {
				for (y, elem) in td2.iter().zip(hdist.dist[ix..ix+n_bins].iter_mut()) {
					*elem += x * y * n;
				}
				ix += n_bins;
			} */
		}
	}

	if let Some(pv) = pvec.take() {
		if !pv.is_empty() {
			tx.as_mut().unwrap().send(pv).expect("Error sending message from subtask")
		}
	}

	debug!("heatmap_subtask({}): closing down", ix);
	hm
}

fn heatmap_calc(ssize: usize, cache: Arc<BetaCache>, lnp: Arc<LnP>, rx: Receiver<Vec<PairObs>>) -> Vec<HeatDist> {

	let n_bins = lnp.n_bins();
	let mut hdists: Vec<_> = (0..ssize).map(|_| HeatDist::new(n_bins)).collect();

	for mut rec in rx.iter() {
		for pobs in rec.drain(..) {
			add_contrib(pobs.n, &pobs.cts, &cache, &lnp, &mut hdists[pobs.idx]);
		}
	}
	hdists
}

fn add_contrib(n: usize, cts: &[u32; 4], cache: &BetaCache, lnp: &LnP, hd: &mut HeatDist) {
	let n_bins = lnp.n_bins();
	let td1 = cache.get_dist(&lnp, cts[0], cts[1]);
	let td2 = cache.get_dist(&lnp, cts[2], cts[3]);
	hd.n += n;
	let nn = n as f64;
	let mut ix = 0;
	for x in td1.iter() {
		for (y, elem) in td2.iter().zip(hd.dist[ix..ix + n_bins].iter_mut()) {
			*elem += x * y * nn
		}
		ix += n_bins;
	}
}