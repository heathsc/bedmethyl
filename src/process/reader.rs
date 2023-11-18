use anyhow::Context;
use crossbeam_channel::{Receiver, Sender};
use lazy_static::lazy_static;
use regex::Regex;

use std::{collections::VecDeque, ops::AddAssign};

use r_htslib::*;

use crate::config::*;

use super::{smooth::*, MsgBlock, BLOCK_SIZE};

use crate::regions::Regions;
use crate::sample::Sample;

/// Observed counts
#[derive(Debug, Copy, Clone, Default)]
pub struct Counts {
    pub non_converted: u32,
    pub converted: u32,
}

impl AddAssign for Counts {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            non_converted: self.non_converted + other.non_converted,
            converted: self.converted + other.converted,
        };
    }
}

impl Counts {
    pub fn is_zero(&self) -> bool {
        self.non_converted + self.converted == 0
    }
    pub fn n(&self) -> u32 {
        self.non_converted + self.converted
    }
}

/// Imputed counts (if smoothing is used)
#[derive(Debug, Copy, Clone, Default)]
pub struct ImpCounts {
    pub non_converted: f64,
    pub converted: f64,
}

impl ImpCounts {
    pub fn imputed(&self) -> bool {
        self.converted + self.non_converted > 0.000001
    }
}

lazy_static! {
    static ref RE_DESC: Regex = Regex::new(r#"description="([^"]*)""#).unwrap();
    static ref RE_RGB: Regex = Regex::new(r#"^[0-9]+,[0-9]+,[0-9]+$"#).unwrap();
    static ref RE_BC: Regex = Regex::new(r#"^(.*)_c[ph][gh].bed"#).unwrap();
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
    Any,
}

fn send_block(blk: MsgBlock, tx: &mut Sender<MsgBlock>) -> anyhow::Result<()> {
    tx.send(blk).with_context(|| "Error sending block")
}

pub(crate) struct InputFile<'a, 'b, 'c> {
    trec: TbxRec,
    iter: HtsItrReader<'a, 'b, TbxRec>,
    pending: Option<Record>,
    current: Option<Record>,
    sample: &'c Sample,
    regions: &'c Regions,
    smooth: Option<SmoothFile<'c>>,
    last_smooth_fit: Option<(u32, SmoothFit)>,
    line: usize,
    merge_strands: bool,
    eof: bool,
    specific_sites: bool,
    format_bug: Option<bool>,
}

impl<'a, 'b, 'c> InputFile<'a, 'b, 'c> {
    fn update(&mut self) -> anyhow::Result<Option<&Record>> {
        if self.current.is_none() {
            if let Some(mut sm) = self.smooth.take() {
                self.current = sm.get_rec(self)?;
                self.smooth = Some(sm);
            } else {
                loop {
                    self.current = self.read_rec()?;
                    if self.specific_sites {
                        if let Some(r) = self.current.as_ref() {
                            let v = self
                                .regions
                                .get(r.reg_idx as usize)
                                .expect("Missing region")
                                .sites()
                                .expect("No sites for region");
                            if v.binary_search(&r.pos).is_ok() {
                                break;
                            }
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                }
            }
        }
        Ok(self.current.as_ref())
    }

    fn set_current(&mut self, rec: Record) {
        self.current = Some(rec);
    }

    // Skip over records that have zero counts (not observed)
    pub(super) fn read_rec(&mut self) -> anyhow::Result<Option<Record>> {
        loop {
            if let Some(rec) = self._read_rec()? {
                if rec.is_observed() {
                    return Ok(Some(rec));
                }
            } else {
                return Ok(None);
            }
        }
    }

    fn _read_rec(&mut self) -> anyhow::Result<Option<Record>> {
        while !self.eof {
            // First see if there is a pending record that can be issued immediately
            if let Some(mut prec) = self.pending.take() {
                match (prec.strand, prec.cg) {
                    (Strand::Plus, Some(true)) | (Strand::Plus, None) => self.pending = Some(prec),
                    (Strand::Minus, Some(true)) => {
                        prec.pos -= 1;
                        prec.strand = Strand::Plus;
                        prec.two_base = true;
                        return Ok(Some(prec));
                    }
                    _ => return Ok(Some(prec)),
                }
            }

            let res = self
                .iter
                .read(&mut self.trec)
                .with_context(|| format!("Error reading from {}", self.sample.path().display()))?;
            self.line += 1;
            if !res {
                // EOF
                self.eof = true;
                break;
            }
            let idx = self
                .iter
                .current_region()
                .and_then(|r| r.idx())
                .expect("Missing region index");

            let rec =
                Record::from_tbx_record(&self.trec, idx, &mut self.format_bug, self.specific_sites)
                    .with_context(|| {
                        format!(
                            "Error reading from {}:{}",
                            self.sample.path().display(),
                            self.line
                        )
                    })?;

            if !self.merge_strands {
                return Ok(Some(rec));
            }

            // If we are merging strands for CpGs we have to see if we have two records from the same CpG

            if let Some(mut prec) = self.pending.take() {
                if rec.strand == Strand::Minus
                    && rec.reg_idx == prec.reg_idx
                    && rec.pos == prec.pos + 1
                {
                    // Pair found so we merge records
                    prec.counts += rec.counts;
                    prec.two_base = true;
                } else {
                    self.pending = Some(rec);
                    if prec.strand == Strand::Plus {
                        // && prec.cg == Some(true) {
                        prec.two_base = true;
                    }
                }
                return Ok(Some(prec));
            } else {
                self.pending = Some(rec)
            }
        }
        // We only get here if we are the end of the input
        if let Some(mut rec) = self.pending.take() {
            if rec.cg == Some(true) {
                if rec.strand == Strand::Minus {
                    rec.strand = Strand::Plus;
                    rec.pos -= 1
                }
                rec.two_base = true;
            }
            return Ok(Some(rec));
        } else {
            self.current = None;
        }
        Ok(None)
    }
}

#[derive(Debug, Copy, Clone)]
pub(super) struct Record {
    pub(super) reg_idx: u32, // Index into cfg.regions
    pub(super) counts: Counts,
    pub(super) smooth_fit: Option<SmoothFit>,
    pub(super) pos: usize,
    strand: Strand,
    cg: Option<bool>,
    two_base: bool,
    pub(super) keep: bool,
}

impl Record {
    pub(super) fn is_observed(&self) -> bool {
        !self.counts.is_zero()
    }

    fn from_tbx_record(
        trec: &TbxRec,
        reg_idx: u32,
        format_bug: &mut Option<bool>,
        specific_sites: bool,
    ) -> anyhow::Result<Self> {
        let mut it = trec
            .to_str()
            .expect("Null tbx record")
            .trim_end()
            .split('\t');

        let (s1, s2, s3, s4) = if format_bug.is_none() {
            let fd: Vec<_> = it.collect();
            let s = fd.get(9).ok_or_else(|| anyhow!("Short input line"))?;
            let bug = RE_RGB.is_match(s);
            let kludge = if bug { 1 } else { 0 };
            *format_bug = Some(bug);
            (
                fd.get(5 + kludge).copied(),
                fd.get(9 + kludge).copied(),
                fd.get(10 + kludge).copied(),
                fd.get(11 + kludge).copied(),
            )
        } else {
            let kludge = if format_bug.unwrap() { 1 } else { 0 };
            let s1 = it.nth(5 + kludge);
            let s2 = it.nth(3);
            let s3 = it.next();
            let s4 = it.next();
            (s1, s2, s3, s4)
        };

        // If the BedMethyl file is from gemBS we can see if this is part of a CpG or not
        let cg = s4.and_then(|s| match s {
            "CG" => Some(true),
            "CHG" | "CHH" => Some(false),
            _ => None,
        });

        // Get strand, non-converted and converted counts
        if let (Some(s1), Some(s2), Some(s3)) = (s1, s2, s3) {
            let cov = s2.parse::<u32>().with_context(|| "parsing coverage")?;
            let meth = s3
                .parse::<f64>()
                .with_context(|| "Error parsing methylation")?;
            if (0.0..=100.0).contains(&meth) {
                let non_converted = ((cov as f64) * 0.01 * meth).round() as u32;
                let converted = cov - non_converted;
                let strand = match s1 {
                    "+" => Strand::Plus,
                    "-" => Strand::Minus,
                    _ => Strand::Any,
                };
                Ok(Record {
                    reg_idx,
                    counts: Counts {
                        non_converted,
                        converted,
                    },
                    smooth_fit: None,
                    pos: (trec.begin() + 1) as usize,
                    strand,
                    cg,
                    two_base: false,
                    keep: !specific_sites,
                })
            } else {
                Err(anyhow!("Illegal methylation value {} (out of range)", meth))
            }
        } else {
            Err(anyhow!("Short line"))
        }
    }
}

#[derive(Debug)]
pub struct MultiRecord {
    pub reg_idx: u32,
    pub pos: usize,
    pub two_base: bool,
    pub strand: Strand,
    pub counts: Vec<Counts>,
    pub imp_counts: Option<Vec<ImpCounts>>,
    pub smooth_fit: Option<Vec<Option<SmoothFit>>>,
}

pub(super) fn reader_thread(
    cfg: &Config,
    mut hts_vec: Vec<Hts>,
    sample_idx: usize,
    mut sender: Sender<MsgBlock>,
    final_level: bool,
) -> anyhow::Result<()> {
    let nf = hts_vec.len();

    let tid = unsafe { libc::syscall(libc::SYS_gettid) };
    //	let tid = unsafe { libc::gettid() };
    debug!("Reader thread starting {}", tid);

    debug!("Creating sample specific region lists");
    let reg_vec: Vec<_> = hts_vec
        .iter_mut()
        .map(|h| {
            let v: Vec<_> = cfg.regions().iter().map(|r| format!("{:#}", r)).collect();
            h.make_region_list(&v)
        })
        .collect();

    // Should we merge the two strands?
    let merge_strands = cfg.combine_cpgs();

    let samples = cfg.sample_info().samples();
    // Make InputFile structs
    let mut ifiles: Vec<_> = hts_vec
        .iter_mut()
        .zip(reg_vec.iter())
        .enumerate()
        .map(|(k, (h, rlist))| {
            let smooth = if cfg.smooth() {
                Some(SmoothFile::new(
                    cfg.min_sites(),
                    cfg.max_distance(),
                    cfg.window_size(),
                    cfg.regions(),
                    cfg.specific_sites(),
                ))
            } else {
                None
            };
            InputFile {
                trec: TbxRec::new(),
                iter: h.itr_reader(rlist),
                pending: None,
                current: None,
                sample: &samples[k + sample_idx],
                regions: cfg.regions(),
                line: 0,
                smooth,
                merge_strands,
                eof: false,
                specific_sites: cfg.specific_sites(),
                format_bug: None,
                last_smooth_fit: None,
            }
        })
        .collect();

    let mut rec_vec = VecDeque::with_capacity(BLOCK_SIZE);
    loop {
        // Find left most entry amongst all files
        let mut curr = None;
        let mut two_base = false;
        for ifile in ifiles.iter_mut() {
            if let Some(rec) = ifile.update()? {
                if let Some((idx, x, _)) = curr {
                    if idx > rec.reg_idx || (idx == rec.reg_idx && x > rec.pos) {
                        curr = Some((rec.reg_idx, rec.pos, rec.strand));
                        two_base = rec.two_base
                    } else if idx == rec.reg_idx && x == rec.pos {
                        two_base = two_base || rec.two_base
                    }
                } else {
                    curr = Some((rec.reg_idx, rec.pos, rec.strand));
                    two_base = rec.two_base
                }
            }
        }
        // Create MultiRecord from matching individual records

        if let Some((idx, x, strand)) = curr {
            trace!("Reader thread {tid} has a multirecord");
            let mut counts = Vec::with_capacity(nf);
            let mut imp_counts = if cfg.smooth() {
                Some(Vec::with_capacity(nf))
            } else {
                None
            };
            let mut sm_vec = if cfg.smooth() && !final_level {
                Some(Vec::with_capacity(nf))
            } else {
                None
            };

            for ifile in ifiles.iter_mut() {
                let mut sm: Option<SmoothFit> = None;
                if let Some(rec) = ifile.current.take().and_then(|mut r| {
                    if let Some((r_idx, _)) = ifile.last_smooth_fit {
                        if r_idx != idx {
                            ifile.last_smooth_fit = None
                        }
                    }
                    if r.reg_idx == idx {
                        if r.pos == x {
                            Some(r)
                        } else if two_base
                            && strand == Strand::Plus
                            && r.pos == x + 1
                            && r.strand == Strand::Minus
                        {
                            // Recover partially observed CpGs where we only have the minus strand
                            r.pos -= 1;
                            r.two_base = true;
                            r.strand = Strand::Plus;
                            ifile.last_smooth_fit =
                                r.smooth_fit.and_then(|sf| Some((r.reg_idx, sf)));
                            Some(r)
                        } else {
                            assert!(x < r.pos);
                            // See if we have a smooth fit close enough to use
                            let dist = ifile.last_smooth_fit.and_then(|(_, sf)| {
                                if x.abs_diff(sf.pos) as f64 <= sf.bandwidth {
                                    Some((x.abs_diff(sf.pos), sf))
                                } else {
                                    None
                                }
                            });
                            let dist1 = r.smooth_fit.and_then(|sf| {
                                if x.abs_diff(sf.pos) as f64 <= sf.bandwidth {
                                    Some(x.abs_diff(sf.pos))
                                } else {
                                    None
                                }
                            });
                            sm = match (dist, dist1) {
                                (Some((x, sf)), Some(y)) => {
                                    if x <= y {
                                        Some(sf)
                                    } else {
                                        r.smooth_fit
                                    }
                                }
                                (Some((_, sf)), None) => Some(sf),
                                (None, Some(_)) => r.smooth_fit,
                                _ => None,
                            };
                            ifile.set_current(r);
                            None
                        }
                    } else {
                        sm = ifile.last_smooth_fit.and_then(|(_, sf)| {
                            if x.abs_diff(sf.pos) as f64 <= sf.bandwidth {
                                Some(sf)
                            } else {
                                None
                            }
                        });
                        ifile.set_current(r);
                        None
                    }
                }) {
                    counts.push(rec.counts);

                    if cfg.smooth() {
                        let ic: ImpCounts = rec
                            .smooth_fit
                            .as_ref()
                            .map(|sf| sf.impute(x))
                            .unwrap_or_default();
                        imp_counts.as_mut().unwrap().push(ic);
                    }
                    if let Some(v) = sm_vec.as_mut() {
                        v.push(rec.smooth_fit);
                    }
                } else {
                    if cfg.smooth() {
                        let ic: ImpCounts = sm.as_ref().map(|sf| sf.impute(x)).unwrap_or_default();
                        imp_counts.as_mut().unwrap().push(ic);
                    }
                    counts.push(Counts::default());
                    if let Some(v) = sm_vec.as_mut() {
                        v.push(sm);
                    }
                }
            }
            let mrec = MultiRecord {
                reg_idx: idx,
                pos: x,
                counts,
                imp_counts,
                two_base,
                strand,
                smooth_fit: sm_vec,
            };
            rec_vec.push_back(mrec);
            if rec_vec.len() >= BLOCK_SIZE {
                send_block(rec_vec, &mut sender)?;
                rec_vec = VecDeque::with_capacity(BLOCK_SIZE);
            }
        } else {
            if !rec_vec.is_empty() {
                send_block(rec_vec, &mut sender)?
            }
            break;
        }
    }

    debug!("Reader thread returning");
    Ok(())
}

struct Channel {
    r: Receiver<MsgBlock>,
    curr_rec: Option<MultiRecord>,
    curr_mb: Option<MsgBlock>,
    last_smooth_fit: Vec<Option<(u32, SmoothFit)>>,
    ns: usize, // number of samples being handled by this channel
    eof: bool,
}

impl Channel {
    fn update(&mut self) -> Option<&MultiRecord> {
        if self.curr_rec.is_none() {
            self.next_rec();
        }
        self.curr_rec.as_ref()
    }

    fn next_rec(&mut self) {
        if let Some(mb) = self.curr_mb.as_mut().and_then(|mb| mb.pop_front()) {
            self.curr_rec = Some(mb);
        } else if !self.eof {
            match self.r.recv() {
                Ok(mb) => {
                    self.curr_mb = Some(mb);
                    self.next_rec()
                }
                Err(_) => {
                    self.eof = true;
                }
            }
        }
    }
}

pub(super) fn merge_thread(
    cfg: &Config,
    mut recv_vec: Vec<(Receiver<MsgBlock>, usize)>,
    mut sender: Sender<MsgBlock>,
    final_level: bool,
) -> anyhow::Result<()> {
    let tid = unsafe { libc::syscall(libc::SYS_gettid) };
    //	let tid = unsafe { libc::gettid() };
    debug!("Starting merge thread {}", tid);

    let mut channels: Vec<_> = recv_vec
        .drain(..)
        .map(|(r, ns)| Channel {
            r,
            curr_rec: None,
            curr_mb: None,
            last_smooth_fit: vec![None; ns],
            ns,
            eof: false,
        })
        .collect();

    let mut rec_vec = VecDeque::with_capacity(BLOCK_SIZE);
    loop {
        let mut curr = None;
        let mut two_base = cfg.combine_cpgs() && cfg.assume_cpg();
        for chan in channels.iter_mut() {
            if let Some(rec) = chan.update() {
                if let Some((idx, x, _)) = curr {
                    if idx > rec.reg_idx || (idx == rec.reg_idx && x > rec.pos) {
                        curr = Some((rec.reg_idx, rec.pos, rec.strand));
                        two_base = rec.two_base
                    } else if idx == rec.reg_idx && x == rec.pos {
                        two_base = two_base || rec.two_base
                    }
                } else {
                    curr = Some((rec.reg_idx, rec.pos, rec.strand));
                    two_base = rec.two_base
                }
            }
        }
        if let Some((idx, x, strand)) = curr {
            let mut counts: Vec<Counts> = Vec::new();
            let mut imp_counts = if cfg.smooth() { Some(Vec::new()) } else { None };
            let mut sm_vec = if cfg.smooth() && !final_level {
                Some(Vec::new())
            } else {
                None
            };

            for chan in channels.iter_mut() {
                let mut sm: Vec<Option<SmoothFit>> = vec![None; chan.ns];
                if let Some(rec) = chan.curr_rec.take().and_then(|mut r| {
                    // Remove last_smooth_fit for any channels where we have changed contig
                    for sf in chan.last_smooth_fit.iter_mut() {
                        if let Some((r_idx, _)) = sf {
                            if *r_idx != idx {
                                *sf = None
                            }
                        }
                    }
                    if r.reg_idx == idx {
                        if r.pos == x {
                            Some(r)
                        } else if two_base
                            && strand == Strand::Plus
                            && r.pos == x + 1
                            && r.strand == Strand::Minus
                        {
                            r.pos -= 1;
                            r.two_base = true;
                            r.strand = Strand::Plus;

                            if cfg.smooth() {
                                // Copy Smooth fit if existing to channel store
                                if let Some(rsf) = r.smooth_fit.as_ref() {
                                    for (sf1, sf2) in
                                        chan.last_smooth_fit.iter_mut().zip(rsf.iter())
                                    {
                                        if let Some(s) = sf2 {
                                            *sf1 = Some((idx, *s))
                                        }
                                    }
                                }
                            }
                            Some(r)
                        } else {
                            assert!(x < r.pos);
                            // see if we have a smooth fit close enough to use for each of the samples in the channel
                            if cfg.smooth() {
                                for ix in 0..chan.ns {
                                    let dist = chan.last_smooth_fit[ix].and_then(|(_, sf)| {
                                        if x.abs_diff(sf.pos) as f64 <= sf.bandwidth {
                                            Some((x.abs_diff(sf.pos), sf))
                                        } else {
                                            None
                                        }
                                    });
                                    let dist1 =
                                        r.smooth_fit.as_ref().and_then(|v| v[ix]).and_then(|sf| {
                                            if x.abs_diff(sf.pos) as f64 <= sf.bandwidth {
                                                Some(x.abs_diff(sf.pos))
                                            } else {
                                                None
                                            }
                                        });
                                    sm[ix] = match (dist, dist1) {
                                        (Some((x, sf)), Some(y)) => {
                                            if x <= y {
                                                Some(sf)
                                            } else {
                                                r.smooth_fit.as_ref().unwrap()[ix]
                                            }
                                        }
                                        (Some((_, sf)), None) => Some(sf),
                                        (None, Some(_)) => r.smooth_fit.as_ref().unwrap()[ix],
                                        _ => None,
                                    };
                                }
                            }
                            chan.curr_rec = Some(r);
                            None
                        }
                    } else {
                        if cfg.smooth() {
                            for (sf1, o) in sm.iter_mut().zip(chan.last_smooth_fit.iter()) {
                                *sf1 = o.and_then(|(_, sf)| {
                                    if x.abs_diff(sf.pos) as f64 <= sf.bandwidth {
                                        Some(sf)
                                    } else {
                                        None
                                    }
                                })
                            }
                        }
                        chan.curr_rec = Some(r);
                        None
                    }
                }) {
                    let MultiRecord {
                        counts: mut cn,
                        imp_counts: mut icn,
                        smooth_fit: svec,
                        ..
                    } = rec;
                    counts.append(&mut cn);
                    if let Some(v) = sm_vec.as_mut() {
                        v.append(&mut svec.unwrap());
                    }
                    if let Some(ic) = imp_counts.as_mut() {
                        ic.append(icn.as_mut().unwrap())
                    }
                } else {
                    for _ in 0..chan.ns {
                        counts.push(Counts::default())
                    }
                    if let Some(ic) = imp_counts.as_mut() {
                        for ix in 0..chan.ns {
                            let icn = sm[ix].as_ref().map(|sf| sf.impute(x)).unwrap_or_default();
                            ic.push(icn);
                        }
                    }
                    if let Some(v) = sm_vec.as_mut() {
                        v.append(&mut sm);
                    }
                }
            }
            let mrec = MultiRecord {
                reg_idx: idx,
                pos: x,
                two_base,
                counts,
                imp_counts,
                strand,
                smooth_fit: sm_vec, // We don't need these anymore
            };
            rec_vec.push_back(mrec);
            if rec_vec.len() >= BLOCK_SIZE {
                send_block(rec_vec, &mut sender).with_context(|| "Error sending merged block")?;
                rec_vec = VecDeque::with_capacity(BLOCK_SIZE);
            }
        } else {
            if !rec_vec.is_empty() {
                send_block(rec_vec, &mut sender)?
            }
            break;
        }
    }
    debug!("Merge thread returning");

    Ok(())
}
