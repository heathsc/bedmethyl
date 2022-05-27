use std::{
	io::BufRead,
	cmp::{Ord, Ordering},
	fmt,
	ops::Deref,
};

use regex::{Regex, RegexSet};
use lazy_static::lazy_static;

use compress_io::compress;

#[derive(Debug)]
pub enum Reg<'a> {
	Chrom(&'a str, bool),
	Open(&'a str, usize, bool),
	Closed(&'a str, usize, usize, bool),
}

fn write_ctg(f: &mut fmt::Formatter, s: &str, quoted: bool) -> fmt::Result {
	if quoted && f.alternate() {
		write!(f, "{{{}}}", s)
	} else {
		write!(f, "{}", s)
	}
}

impl <'a>fmt::Display for Reg<'a> {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		match self {
			Reg::Chrom(a, quoted) => write_ctg(f, *a, *quoted),
			Reg::Open(a, x, quoted) => {
				write_ctg(f, *a, *quoted)?;
				write!(f, ":{}", *x)
			},
			Reg::Closed(a, x, y, quoted) if *x == 0 => {
				write_ctg(f, *a, *quoted)?;
				write!(f, ":-{}", *y)
			},
			Reg::Closed(a, x, y, quoted) => {
				write_ctg(f, *a, *quoted)?;
				write!(f, ":{}-{}", *x, *y)
			},
		}
	}
}

impl <'a>Ord for Reg<'a> {
	fn cmp(&self, other: &Self) -> Ordering {
		match (self, other) {
			(Reg::Chrom(a, _), Reg::Chrom(b, _)) => a.cmp(b),
			(Reg::Chrom(a, _), Reg::Open(b, ..)) | (Reg::Chrom(a, _), Reg::Closed(b, ..)) => if a == b {Ordering::Less } else { a.cmp(b) },
			(Reg::Open(a, ..), Reg::Chrom(b, _)) | (Reg::Closed(a, ..), Reg::Chrom(b, _)) => if a == b { Ordering::Greater } else { a.cmp(b) },
			(Reg::Open(a, x1, ..), Reg::Open(b, x2, ..)) => if a == b { x1.cmp(x2) } else { a.cmp(b) },
			(Reg::Open(a, x1, ..), Reg::Closed(b, x2, ..)) => if a == b { if x1 == x2 { Ordering::Greater } else {	x1.cmp(x2) }} else { a.cmp(b) },
			(Reg::Closed(a, x1, ..), Reg::Open(b, x2, ..)) => if a == b { if x1 == x2 { Ordering::Less } else {	x1.cmp(x2) }} else { a.cmp(b) },
			(Reg::Closed(a, x1, y1, ..), Reg::Closed(b, x2, y2, ..)) => if a == b { if x1 == x2 { y1.cmp(y2) } else {	x1.cmp(x2) }} else { a.cmp(b) },
		}
	}	
}

impl <'a>PartialOrd for Reg<'a> {
	fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
		Some(self.cmp(other))
	}	
}

impl <'a>PartialEq for Reg<'a> {
	fn eq(&self, other: &Self) -> bool {
		match (self, other) {
			(Reg::Chrom(a, ..), Reg::Chrom(b, ..)) => a == b,
			(Reg::Open(a, x1, ..), Reg::Open(b, x2, ..)) => a == b && x1 == x2,
			(Reg::Closed(a, x1, y1, ..), Reg::Closed(b, x2, y2, ..)) => a == b && x1 == x2 && y1 == y2,
			_ => false,
		}
	}
}

lazy_static!{
	// Matches when the contig is disambiguated using brackets i.e.., {chr2}:20000-50000
	// The Regex for the contig name comes from the VCF4.3 spec
	static ref RE_REGION1: Regex = Regex::new(r#"^[{]([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)[}]:?([0-9,]+)?-?([0-9,]+)?"#).unwrap();

	// Matches when the contig is present without brackets i.e., chr2:20000-30000
	static ref RE_REGION2: Regex = Regex::new(r#"^([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./;=?@^_|~-]*):?([0-9,]+)?-?([0-9,]+)?"#).unwrap();

	// RegexSet to check for valid contig
	static ref RE_SET_CONTIG: RegexSet = RegexSet::new(&[
		r"^[{]([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)[}]$",
		r"^([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./;=?@^_|~-]*)$"
	]).unwrap();
}

#[derive(Debug)]
enum ChkRegOverlap<'a> {
	Yes(Reg<'a>),
	No(Reg<'a>, Reg<'a>),	
}

fn parse_usize_with_commas(s: &str) -> Option<usize> {
	s.replace(',', "").parse::<usize>().ok()
}

fn parse_range(s1: &str, s2: &str) -> Option<(usize, usize)> {
	match (parse_usize_with_commas(s1), parse_usize_with_commas(s2)) {
		(Some(x), Some(y)) => Some((x, y)),
		_ => None,
	}	
}

impl <'a> Reg<'a> {
	// Assume that self <= other
	fn check_overlap(mut self, other: Self) -> ChkRegOverlap<'a> {
		match (&self, &other) {
			(Reg::Chrom(a, ..), Reg::Chrom(b, ..)) | (Reg::Chrom(a, ..), Reg::Open(b, ..)) | (Reg::Chrom(a, ..), Reg::Closed(b, ..))
				| (Reg::Open(a, ..), Reg::Open(b, ..)) | (Reg::Open(a, ..), Reg::Closed(b, ..)) => if a == b {
				ChkRegOverlap::Yes(self)
			} else {
				ChkRegOverlap::No(self, other)
			},
			(Reg::Open(a, ..), Reg::Chrom(b, ..)) | (Reg::Closed(a, ..), Reg::Chrom(b, ..)) => if a == b {
				ChkRegOverlap::Yes(other)
			} else {
				ChkRegOverlap::No(self, other)
			},
			(Reg::Closed(a, x1, y1, quoted), Reg::Open(b, x2, ..)) => if a == b && y1 + 1 >= *x2 {
				self = Reg::Open(a, *x1, *quoted);
				ChkRegOverlap::Yes(self)
			} else {
				ChkRegOverlap::No(self, other)
			},	
			(Reg::Closed(a, x1, y1, quoted), Reg::Closed(b, x2, y2, ..)) => if a == b {
				if x1 == x2 {
					ChkRegOverlap::Yes(other)
				} else {
					self = Reg::Closed(a, *x1, *y1.max(y2), *quoted);
					ChkRegOverlap::Yes(self)
				}
			} else {
				ChkRegOverlap::No(self, other)
			}
		}
	}
	
	fn from_region_str(s: &'a str) -> Option<Self> {
		if let Some((cap, quoted)) = RE_REGION1.captures(s)
			.map(|c| (c, true)).or_else(|| RE_REGION2.captures(s).map(|c| (c, false))) {
			debug!("{:?} {:?} {:?} {}", cap.get(1), cap.get(2), cap.get(2), quoted);
			match (cap.get(1), cap.get(2), cap.get(3)) {
				(Some(c), None, None) => Some(Reg::Chrom(c.as_str(), quoted)),
				(Some(c), Some(p), None) => parse_usize_with_commas(p.as_str()).map(|x| Reg::Open(c.as_str(), x, quoted)),
				(Some(c), None, Some(q)) => parse_usize_with_commas(q.as_str()).map(|x| Reg::Closed(c.as_str(), 0, x, quoted)),
				(Some(c), Some(p), Some(q)) => parse_range(p.as_str(), q.as_str()).map(|(x, y)| Reg::Closed(c.as_str(), x, y, quoted)),
				_ => None,
			}
		} else {
			None
		}
	} 
}

impl <'a>Eq for Reg<'a> { }

pub struct Region {
	contig: Box<str>,
	start: usize,
	end: Option<usize>,
	quoted: bool,
}

impl fmt::Display for Region {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		match (self.start, self.end) {
			(0, None) => write_ctg(f, self.contig.as_ref(), self.quoted),
			(0, Some(y)) => {
			write_ctg(f, self.contig.as_ref(), self.quoted)?;
				write!(f, ":-{}", y)
			},
			(x, None) => {
				write_ctg(f, self.contig.as_ref(), self.quoted)?;
				write!(f, ":{}-", x)
			},
			(x, Some(y)) => {
				write_ctg(f, self.contig.as_ref(), self.quoted)?;
				write!(f, ":{}-{}", x, y)
			},
		}
	}
}

impl Region {
	pub fn from_reg(reg: &Reg) -> Self {
		debug!("Converting reg {} to region", reg);
		match reg {
			Reg::Chrom(s, quoted) => {
				debug!("Chrom: {}", s);
				Self {
					contig: s.to_string().into_boxed_str(),
					start: 0,
					end: None,
					quoted: *quoted,
				}
			},
			Reg::Open(s, x, quoted) => {
				debug!("Open: {} {}", s, x);
				Self {
					contig: s.to_string().into_boxed_str(),
					start: *x,
					end: None,
					quoted: *quoted,
				}
			},
			Reg::Closed(s, x, y, quoted) => {
				debug!("Closed: {} {} {}", s, x, y);
				Self {
					contig: s.to_string().into_boxed_str(),
					start: *x,
					end: Some(*y),
					quoted: *quoted,
				}
			},
		}
	}
	pub fn from_ctg(ctg: &str) -> Self {
		Self {
			contig: ctg.to_string().into_boxed_str(),
			start: 0,
			end: None,
			quoted: ctg.contains(':'),
		}
	}
	pub fn contig(&self) -> &str { self.contig.as_ref() }
}

pub struct Regions {
	regs: Vec<Region>,
}

impl Deref for Regions {
	type Target = [Region];

	fn deref(&self) -> &Self::Target {
		&self.regs
	}
}

impl Regions {
	pub fn from_regions<'a, I: IntoIterator<Item = &'a str>>(regions: I) -> Self {	
		let v: Vec<_> = regions.into_iter().filter_map(Reg::from_region_str).collect();
		let regs: Vec<_> = merge_regions(v).iter().map(|r| Region::from_reg(r)).collect();
		Self{regs}
	}	
	
	pub fn from_file<P: AsRef<str>>(file: P) -> anyhow::Result<Self> {
		
		let fname = file.as_ref();
		let mut buf = String::new();
		let mut f = compress::CompressIo::new().path(fname).bufreader()?;
		
		let mut reg_str = Vec::new();
		loop {
			buf.clear();
		
			if f.read_line(&mut buf)? == 0 {
				break;
			}
			let fd: Vec<_> = buf.trim().split('\t').collect();
			
			if fd.len() > 1 {
				let m = RE_SET_CONTIG.matches(fd[0]);
				let quoted = match (m.matched(0), m.matched(1)) {
					(true, _) => true,
					(false, true) => false,
					_ => return Err(anyhow!("Illegal contig name {}", fd[0]))
				};
				let start = fd.get(1).and_then(|s| s.parse::<usize>().ok());
				let stop = fd.get(2).and_then(|s| s.parse::<usize>().ok());
				reg_str.push((fd[0].to_owned(), start, stop, quoted));
			}
		}
		let mut regs = Vec::new();
		for (chr, start, stop, quoted) in reg_str.iter() {
			regs.push(match (start, stop) {
				(Some(x), Some(y)) => Reg::Closed(chr, *x, *y, *quoted),
				(Some(x), None) => Reg::Closed(chr, *x, *x, *quoted),
				_ => Reg::Chrom(chr, *quoted),
			});
		}	
		let regs: Vec<_> = merge_regions(regs).iter().map(|r| Region::from_reg(r)).collect();
		Ok(Self{regs})
	}

	pub fn from_str_vec(ctgs: &[&str]) -> Self {
		let regs: Vec<_> = ctgs.iter().map(|s| Region::from_ctg(s)).collect();
		Self { regs }
	}
}

fn merge_regions(mut regs: Vec<Reg>) -> Vec<Reg> {
		regs.sort_unstable();
		let mut merged = Vec::with_capacity(regs.len());
		let mut prev: Option<Reg> = None;
		for r in regs.into_iter() {
			prev = Some(if let Some(p) = prev.take() {
				match p.check_overlap(r) {
					ChkRegOverlap::Yes(a) => a,
					ChkRegOverlap::No(a, b) => {
						merged.push(a);
						b
					},
				}
			} else {
				r
			});
		}
		if let Some(p) = prev.take() {
			merged.push(p)
		}
		merged	
}
