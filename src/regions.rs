use std::num::NonZeroUsize;
use std::{
    cmp::{Ord, Ordering},
    collections::{btree_map, BTreeMap, HashMap},
    fmt,
    io::BufRead,
    slice,
    sync::Arc,
};

use anyhow::Context;

use lazy_static::lazy_static;
use regex::{Regex, RegexSet};

use compress_io::compress;

#[derive(Debug)]
pub enum Reg<'a> {
    Chrom(&'a str),
    Open(&'a str, usize),
    Closed(&'a str, usize, usize),
}

fn write_ctg(f: &mut fmt::Formatter, s: &str) -> fmt::Result {
    if f.alternate() && s.contains(':') {
        write!(f, "{{{}}}", s)
    } else {
        write!(f, "{}", s)
    }
}

impl<'a> fmt::Display for Reg<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Reg::Chrom(a) => write_ctg(f, *a),
            Reg::Open(a, x) => {
                write_ctg(f, *a)?;
                write!(f, ":{}", *x)
            }
            Reg::Closed(a, x, y) if *x == 0 => {
                write_ctg(f, *a)?;
                write!(f, ":-{}", *y)
            }
            Reg::Closed(a, x, y) => {
                write_ctg(f, *a)?;
                write!(f, ":{}-{}", *x, *y)
            }
        }
    }
}

lazy_static! {
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
enum ChkRegOverlap {
    Yes(Region),
    No(Region, Region),
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

impl<'a> Reg<'a> {
    fn from_region_str(s: &'a str) -> Option<Self> {
        if let Some(cap) = RE_REGION1.captures(s).or_else(|| RE_REGION2.captures(s)) {
            match (cap.get(1), cap.get(2), cap.get(3)) {
                (Some(c), None, None) => Some(Reg::Chrom(c.as_str())),
                (Some(c), Some(p), None) => {
                    parse_usize_with_commas(p.as_str()).map(|x| Reg::Open(c.as_str(), x))
                }
                (Some(c), None, Some(q)) => {
                    parse_usize_with_commas(q.as_str()).map(|x| Reg::Closed(c.as_str(), 0, x))
                }
                (Some(c), Some(p), Some(q)) => {
                    parse_range(p.as_str(), q.as_str()).map(|(x, y)| Reg::Closed(c.as_str(), x, y))
                }
                _ => None,
            }
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Region {
    contig: Arc<str>,
    start: usize,
    end: Option<usize>,
    sites: Option<Vec<usize>>,
}

impl Ord for Region {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.contig.cmp(&other.contig) {
            Ordering::Equal => match self.start.cmp(&other.start) {
                Ordering::Equal => match (self.end.as_ref(), other.end.as_ref()) {
                    (Some(x), Some(y)) => x.cmp(y),
                    (Some(_), None) => Ordering::Less,
                    (None, Some(_)) => Ordering::Greater,
                    _ => Ordering::Equal,
                },
                x => x,
            },
            x => x,
        }
    }
}

impl PartialOrd for Region {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.start, self.end) {
            (0, None) => write_ctg(f, self.contig.as_ref()),
            (0, Some(y)) => {
                write_ctg(f, self.contig.as_ref())?;
                write!(f, ":-{}", y)
            }
            (x, None) => {
                write_ctg(f, self.contig.as_ref())?;
                write!(f, ":{}-", x)
            }
            (x, Some(y)) => {
                write_ctg(f, self.contig.as_ref())?;
                write!(f, ":{}-{}", x, y)
            }
        }
    }
}

impl Region {
    pub fn contig(&self) -> &str {
        self.contig.as_ref()
    }
    pub fn sites(&self) -> Option<&Vec<usize>> {
        self.sites.as_ref()
    }

    // Assume that self <= other
    fn check_overlap(mut self, other: Self) -> ChkRegOverlap {
        if self.contig != other.contig {
            ChkRegOverlap::No(self, other)
        } else {
            match (self.end.as_ref(), other.end.as_ref()) {
                (Some(x), Some(y)) => {
                    if other.start <= x + 1 {
                        self.end = Some(*x.max(y));
                        ChkRegOverlap::Yes(self)
                    } else {
                        ChkRegOverlap::No(self, other)
                    }
                }
                (Some(x), None) => {
                    if other.start <= x + 1 {
                        self.end = None;
                        ChkRegOverlap::Yes(self)
                    } else {
                        ChkRegOverlap::No(self, other)
                    }
                }
                (None, Some(_)) => ChkRegOverlap::Yes(self),
                (None, None) => ChkRegOverlap::Yes(self),
            }
        }
    }

    fn overlap_interval(&self, x: usize, y: usize) -> Option<(usize, usize)> {
        if self.start > y {
            None
        } else {
            match self.end.as_ref() {
                Some(z) => {
                    if *z < x {
                        None
                    } else {
                        Some((self.start.max(x), y.min(*z)))
                    }
                }
                None => Some((self.start.max(x), y)),
            }
        }
    }
}

#[derive(Default)]
pub struct Regions {
    ctg_reg: BTreeMap<Arc<str>, Vec<Region>>,
    all_regs: Option<Vec<Region>>,
}

impl Regions {
    pub fn region_from_reg(&mut self, reg: &Reg) -> Region {
        debug!("Converting reg {} to region", reg);
        match reg {
            Reg::Chrom(s) => {
                debug!("Chrom: {}", s);
                Region {
                    contig: self.add_contig(s),
                    start: 0,
                    end: None,
                    sites: None,
                }
            }
            Reg::Open(s, x) => {
                debug!("Open: {} {}", s, x);
                Region {
                    contig: self.add_contig(s),
                    start: *x,
                    end: None,
                    sites: None,
                }
            }
            Reg::Closed(s, x, y) => {
                debug!("Closed: {} {} {}", s, x, y);
                Region {
                    contig: self.add_contig(s),
                    start: *x,
                    end: Some(*y),
                    sites: None,
                }
            }
        }
    }

    pub fn iter(&self) -> RegionIter {
        let mut reg_vec = self.ctg_reg.values();
        let regs = reg_vec.next().map(|v| v.iter());

        RegionIter { reg_vec, regs }
    }

    fn mk_all_regs(&mut self) {
        if self.all_regs.is_none() {
            let n = self.ctg_reg.values().map(|v| v.len()).sum();
            let mut reg_vec = Vec::with_capacity(n);
            for v in self.ctg_reg.values() {
                reg_vec.extend(v.iter().cloned())
            }
            self.all_regs = Some(reg_vec)
        }
    }

    pub fn get(&self, idx: usize) -> Option<&Region> {
        self.all_regs.as_ref().and_then(|v| v.get(idx))
    }

    fn add_contig(&mut self, ctg: &str) -> Arc<str> {
        let ctg = Arc::from(ctg);
        if !self.ctg_reg.contains_key(&*ctg) {
            self.ctg_reg.insert(Arc::clone(&ctg), Vec::new());
            self.all_regs = None;
        }
        ctg
    }

    fn add_region(&mut self, region: Region) {
        self.all_regs = None;
        self.ctg_reg
            .get_mut(region.contig.as_ref())
            .unwrap()
            .push(region)
    }

    fn add_reg(&mut self, reg: Reg) {
        let region = self.region_from_reg(&reg);
        self.add_region(region)
    }

    pub fn from_regions<'b, I: IntoIterator<Item = &'b str>>(reg_str: I) -> Self {
        let mut regions = Self::default();

        for reg in reg_str.into_iter().filter_map(Reg::from_region_str) {
            regions.add_reg(reg)
        }
        regions.merge()
    }

    pub fn from_file<P: AsRef<str>>(file: P) -> anyhow::Result<Self> {
        let fname = file.as_ref();
        let mut buf = String::new();
        let mut f = compress::CompressIo::new().path(fname).bufreader()?;
        let mut regions = Self::default();

        loop {
            buf.clear();

            if f.read_line(&mut buf)? == 0 {
                break;
            }
            let fd: Vec<_> = buf.trim().split('\t').collect();

            if !fd.is_empty() {
                if !RE_SET_CONTIG.is_match(fd[0]) {
                    return Err(anyhow!("Illegal contig name {}", fd[0]));
                }
                let start = fd.get(1).and_then(|s| s.parse::<usize>().ok());
                let stop = fd.get(2).and_then(|s| s.parse::<usize>().ok());
                regions.add_reg(match (start, stop) {
                    (Some(x), Some(y)) => Reg::Closed(fd[0], x, y),
                    (Some(x), None) => Reg::Closed(fd[0], x, x),
                    _ => Reg::Chrom(fd[0]),
                })
            }
        }
        Ok(regions.merge())
    }

    pub fn from_str_vec(ctgs: &[&str]) -> Self {
        let mut regions = Regions::default();
        for s in ctgs.iter() {
            let region = Region {
                contig: regions.add_contig(*s),
                start: 0,
                end: None,
                sites: None,
            };
            regions.add_region(region)
        }
        regions.merge()
    }

    fn merge(mut self) -> Self {
        for v in self.ctg_reg.values_mut() {
            v.sort_unstable();
            let mut nv = Vec::new();
            let mut prev: Option<Region> = None;
            for r in v.drain(..) {
                prev = Some(if let Some(p) = prev.take() {
                    match p.check_overlap(r) {
                        ChkRegOverlap::Yes(a) => a,
                        ChkRegOverlap::No(a, b) => {
                            nv.push(a);
                            b
                        }
                    }
                } else {
                    r
                });
            }
            if let Some(p) = prev.take() {
                nv.push(p)
            }
            *v = nv;
        }
        self.mk_all_regs();
        self
    }

    pub fn n_regions(&self) -> usize {
        self.ctg_reg.values().map(|v| v.len()).sum()
    }

    pub fn trim(&mut self, sites: &Sites, win_size: usize) {
        debug!("Trimming regions to match sites");
        let mut nreg = Regions::default();
        for (ctg, v) in self.ctg_reg.iter() {
            if let Some(inter) = sites.intervals(win_size, ctg) {
                let mut site_it = sites.sites(ctg).unwrap().iter().peekable();
                let mut it = v.iter();
                let mut curr_reg = it.next();
                for (x, y) in inter.iter() {
                    if let Some(region) = curr_reg {
                        if let Some((a, b)) = region.overlap_interval(*x, *y) {
                            let contig = nreg.add_contig(ctg);
                            let mut vs = Vec::new();
                            loop {
                                if site_it.peek().map(|&&x| x <= b).unwrap_or(false) {
                                    let x = site_it.next().unwrap();
                                    if *x >= a {
                                        vs.push(*x)
                                    }
                                } else {
                                    break;
                                }
                            }
                            assert!(!vs.is_empty());
                            let r = Region {
                                contig,
                                start: a,
                                end: Some(b),
                                sites: Some(vs),
                            };
                            nreg.add_region(r)
                        }
                    } else {
                        curr_reg = it.next();
                        if curr_reg.is_none() {
                            break;
                        }
                    }
                }
            }
        }
        nreg.mk_all_regs();
        self.ctg_reg = nreg.ctg_reg;
        self.all_regs = nreg.all_regs;
    }
}

pub struct RegionIter<'a> {
    reg_vec: btree_map::Values<'a, Arc<str>, Vec<Region>>,
    regs: Option<slice::Iter<'a, Region>>,
}

impl<'a> Iterator for RegionIter<'a> {
    type Item = &'a Region;

    fn next(&mut self) -> Option<Self::Item> {
        if self.regs.is_none() {
            None
        } else {
            match self.regs.as_mut().and_then(|v| v.next()) {
                Some(x) => Some(x),
                None => {
                    self.regs = self.reg_vec.next().map(|v| v.iter());
                    self.next()
                }
            }
        }
    }
}

pub struct Sites {
    contig_site: HashMap<Box<str>, Vec<usize>>,
}

impl Sites {
    pub fn from_file<P: AsRef<str>>(file: P) -> anyhow::Result<Self> {
        let fname = file.as_ref();
        debug!("Reading in site list from {fname}");
        let mut buf = String::new();
        let mut f = compress::CompressIo::new().path(fname).bufreader()?;
        let mut contig_site = HashMap::new();
        let mut line = 0;
        let mut n_sites = 0;
        loop {
            buf.clear();
            if f.read_line(&mut buf)? == 0 {
                break;
            }
            line += 1;
            let fd: Vec<_> = buf.trim().split('\t').collect();
            // Expecting at laest 2 columns with contig and position (1 based)
            if fd.len() > 1 {
                // Try to parse position
                let x = fd[1].parse::<NonZeroUsize>().with_context(|| {
                    format!(
                        "{}: {} - Could not parse position for site: {}",
                        fname, line, fd[1]
                    )
                })?;
                if !contig_site.contains_key(fd[0]) {
                    let ctg = fd[0].to_owned().into_boxed_str();
                    contig_site.insert(ctg, Vec::new());
                }
                contig_site.get_mut(fd[0]).unwrap().push(usize::from(x) - 1);
                n_sites += 1;
            }
        }
        debug!("Read in {n_sites} sites on {} contigs", contig_site.len());
        // Sort sites within contigs
        for v in contig_site.values_mut() {
            v.sort_unstable();
            // Remove duplicates
            v.dedup()
        }

        Ok(Sites { contig_site })
    }

    pub fn sites(&self, ctg: &str) -> Option<&Vec<usize>> {
        self.contig_site.get(ctg)
    }

    pub fn intervals(&self, win_size: usize, ctg: &str) -> Option<Vec<(usize, usize)>> {
        self.contig_site.get(ctg).map(|v| {
            let mut nv = Vec::new();
            let mut min = v[0].saturating_sub(win_size);
            let mut max = v[0].saturating_add(win_size);
            for x in v[1..].iter() {
                let x1 = x.saturating_sub(win_size);
                if x1 > max + 1 {
                    nv.push((min, max));
                    min = x1;
                }
                max = x.saturating_add(win_size);
            }
            nv.push((min, max));
            nv
        })
    }
}
