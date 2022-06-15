use std::{
	collections::HashMap,
	fmt,
	io::{self, Write},
};


/// Struct for accumulating variance in one pass
/// avoiding round off errors
#[derive(Default, Debug, Copy, Clone)]
pub struct VarStore {
	n: usize,
	mean: f64,
	sum_sq: f64,
}

impl VarStore {
	pub fn new() -> Self { Self::default() }
	
	pub fn add(&mut self, x: f64) {
		let delta = x - self.mean;
		self.n += 1;
		self.mean += delta / (self.n as f64);
		self.sum_sq += (x - self.mean) * delta;
	}
	
	pub fn n(&self) -> usize { self.n }
	
	pub fn mean(&self) -> f64 { self.mean }
	
//	pub fn sum_sq(&self) -> f64 { self.sum_sq }
	
	pub fn var(&self) -> f64 {
		if self.n < 2 {
			0.0
		} else {
			self.sum_sq / ((self.n - 1) as f64)
		}
	}
}

pub struct SampleStats {
	sample: String,
	barcode: String,
	hist: HashMap<u32, u32>,
	var: VarStore,
}

impl SampleStats {
	pub fn new<P: AsRef<str>, Q: AsRef<str>>(barcode: P, sample: Q) -> Self { 
		Self {
			sample: sample.as_ref().to_owned(),
			barcode: barcode.as_ref().to_owned(),
			hist: HashMap::new(),
			var: VarStore::new(),
		}
	}
	
	pub fn add_site(&mut self, a: u32, b: u32) {
		let ct = self.hist.entry(a + b).or_insert(0);
		*ct += 1;
		let x = ((a + 1) as f64) / ((a + b + 2) as f64);
		self.var.add(x);
	}
}

impl fmt::Display for SampleStats {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		
		write!(f, "{}\t{}\t{}\t{:.4}\t{:.4}", self.barcode, self.sample, self.var.n(), self.var.mean(), self.var.var())?;
		if self.hist.is_empty() {
			write!(f, "\tNA\tNA\tNA\tNA\tNA")
		} else {
			let mut v: Vec<(u32, u32)> = self.hist.iter().map(|(a,b)| (*a, *b)).collect();
			v.sort_unstable_by_key(|(ct, _)| *ct);
			let tot = v.iter().fold(0, |s, (ct, n)| s + ct * n);
			let min = v[0].0;
			let max = v[v.len() - 1].0;
			let c1 = tot >> 2;
			let c2 = tot >> 1;
			let c3 = (tot * 3) >> 2;
			let mut s = 0;
			let mut q1 = None;
			let mut q2 = None;
			let mut q3 = None;
			for (ct, n) in v.iter() {
				s += ct * n;
				if q1.is_none() && s >= c1 {
					q1 = Some(*ct)
				}
				if q2.is_none() && s >= c2 {
					q2 = Some(*ct)
				}
				if q3.is_none() && s >= c3 {
					q3 = Some(*ct)
				}
			}
			write!(f, "\t{}\t{}\t{}\t{}\t{}", min, q1.unwrap(), q2.unwrap(), q3.unwrap(), max)
		}	
	}
}

const NBINS:usize = 5000;
const BFACTOR:usize = NBINS << 1;
const NBINSF:f64 = NBINS as f64;

pub struct SiteStats {
	sdev_dist: Vec<usize>,
}

impl Default for SiteStats {
	fn default() -> Self {
		Self{ sdev_dist: vec!(0; NBINS + 1)}	
	}
}

impl SiteStats {
	pub fn add_sdev(&mut self, sd: f64) {
		self.sdev_dist[(sd * NBINSF).round() as usize] += 1
	}
	
	pub fn write_report<W: Write>(&self, w: &mut W) -> io::Result<()> {
		let (n_sites, sum) = self.sdev_dist.iter().enumerate().fold((0, 0), |(s1, s2), (i, ct)| (s1 + ct, s2 + i * ct));
		if sum == 0 { return Ok(()) }
		writeln!(w, "Total sites: {}", n_sites)?;
		let mut s = 0;
		let mut q1 = None;
		let mut q2 = None;
		let mut q3 = None;
		let mut max = 0;
		let c1 = sum >> 2;
		let c2 = sum >> 1;
		let c3 = (sum * 3) >> 2;
		for (i, ct) in self.sdev_dist.iter().enumerate() {
			s += i * ct;
			if q1.is_none() && s >= c1 {
				q1 = Some(i)
			}
			if q2.is_none() && s >= c2 {
				q2 = Some(i)
			}
			if q3.is_none() && s >= c3 {
				q3 = Some(i)
			}
			if *ct > 0 {
				max = i 
			}
		}
		let factor = BFACTOR as f64;
		let den = factor * (n_sites as f64);
		writeln!(w, "Mean sdev: {:.4}", (sum as f64)/ den)?;
		writeln!(w, "Max sdev: {:.4}", (max as f64) / factor)?;
		writeln!(w, "Median sdev: {:.4}, Q1: {:.4}, Q3: {:.4}", (q2.unwrap() as f64) / factor, (q1.unwrap() as f64) / factor, (q3.unwrap() as f64) / factor)?;
		writeln!(w, "\nPercentiles of sdev distribution\npct\tsdev")?;
		let mut pc = 0;
		s = 0;
		for (i, ct) in self.sdev_dist.iter().enumerate() {
			s += i * ct;
			while s >= (sum * pc) / 100 {
				writeln!(w, "{}\t{:.4}", pc, (i as f64) / factor)?;
				pc += 1;
			}
		}		
		writeln!(w, "\nDistribution of sdev\nsdev\tfreq\tcum. freq")?;
		s = 0;
		for (i, ct) in self.sdev_dist.iter().enumerate() {
			s += ct;
			writeln!(w, "{:.4}\t{:.6}\t{:.6}", (i as f64) / factor, (*ct as f64) / (n_sites as f64), (s as f64) / (n_sites as f64))?;
			if i == max {
				break
			}
		}			
		Ok(())
	}
}
