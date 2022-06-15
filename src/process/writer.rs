use std::{
	io::{Write, BufWriter},
	fs::File,
};
use anyhow::Context;

use crossbeam_channel::{Sender, Receiver};

use crate::{
	config::*,
	stats::{SampleStats, SiteStats, VarStore},
	sample::Sample,
};
use crate::process::print_value::PrintValue;

use super::{
	reader::{ImpCounts, Counts},
	print_value::{get_print_value_u32, get_print_value_f64},
	distance::Distance,
	hist::Hist,
	heatmaps::{HeatMaps, Obs},
	MsgBlock,
};

fn print_header<T>(pv: Option<&mut Box<dyn PrintValue<T>>>, samples: &[Sample]) -> anyhow::Result<()> {
	if let Some(p) = pv {
		p.print_str("chrom\tstart\tstop\tn_pass\tsdev")?;
		for s in samples.iter().map(|sample| sample.name()) { p.print_header(s)? }
		p.print_str("\n")?;
	}
	Ok(())
}

pub(super) fn smooth_writer_thread(cfg: &Config, recv: Receiver<MsgBlock>) -> anyhow::Result<()> {

	let mut pv = match get_print_value_f64(&cfg, true)? {
		Some(p) => p,
		None => return Ok(())
	};

	let tid = unsafe { libc::gettid() };
	debug!("Smooth writer thread: {:?}", tid);

	info!("Starting writing of smoothed values:");
	let min_sdev = cfg.min_sdev();

	let samples = cfg.sample_info().samples();

	// Print header lines if required
	print_header(Some(&mut pv), samples)?;

	while let Ok(mblk) = recv.recv() {
		for rec in mblk.into_iter() {
			let mut var = VarStore::new();
			for (a, b) in rec.imp_counts.as_ref().unwrap().iter().map(|ic| {
				let ImpCounts { non_converted: ncnv, converted: cnv } = ic;
				(*ncnv, *cnv)
			}) {
				if a + b > 0.0001 {
					var.add((a + 1.0) as f64 / ((a + b + 2.0) as f64))
				}
			}
			let sd = var.var().sqrt();
			let chrom = cfg.regions().get(rec.reg_idx as usize).expect("Missing contig!").contig();
			if sd >= min_sdev {
				pv.print_str(format!("{}\t{}\t{}\t{}\t{:.4}", chrom, rec.pos, rec.pos + if rec.two_base { 2 } else { 1 }, var.n(), var.var().sqrt()).as_str())?;

				if let Some(ict) = rec.imp_counts.as_ref() {
					if cfg.round_counts() {
						for ct in ict.iter() {
							if ct.imputed() {
								let a = ct.non_converted.round();
								let b = ct.converted.round();
								pv.print_value(a, b)?
							} else { pv.print_missing()? }
						}
					} else {
						for ct in ict.iter() {
							if ct.imputed() {
								pv.print_value(ct.non_converted, ct.converted)?
							} else { pv.print_missing()? }
						}
					}
				}
				pv.print_str("\n")?
			}
		}
	}
	info!("Finished smooth output");
	Ok(())
}

pub(super) fn writer_thread(cfg: &Config, recv: Receiver<MsgBlock>, mut smooth_tr: Option<Sender<MsgBlock>>) -> anyhow::Result<()> {

	let tid = unsafe { libc::gettid() };
	debug!("Writer thread: {:?}", tid);

	info!("Starting processing:");

	// Keep track of the current contig so that we can update the
	let mut curr_contig: Option<&str> = None;

	let mut pv = get_print_value_u32(&cfg, false)?;

	let min_counts = cfg.min_counts();
	let min_samples = cfg.min_samples();
	let min_sdev = cfg.min_sdev();

	let mut site_stats = SiteStats::default();
	let samples = cfg.sample_info().samples();

	// Print header lines if required
	print_header(pv.as_mut(), samples)?;

	let mut stats: Vec<_> = samples.iter().map(|s| SampleStats::new(s.name(), s.name())).collect();
	let mut hist = if cfg.distribution() {
		Some(Hist::new(samples.len(), cfg))
	} else { None };

	let mut heatmaps = if cfg.heatmaps() {
		Some(HeatMaps::new(samples.len(), cfg))
	} else { None };

	let mut distance = cfg.similarity_type().map(|stype| Distance::new(samples.len(), stype));

	while let Ok(mblk) = recv.recv() {
		for rec in mblk.iter() {
			let mut var = VarStore::new();
			for (a, b) in rec.counts.iter().map(|c| {
				let Counts{non_converted: ncnv, converted: cnv} = c;
				(*ncnv, *cnv)
			}) {
				if a + b >= min_counts {
					var.add((a + 1) as f64 / ((a + b + 2) as f64));
				}
			}
			let chrom = cfg.regions().get(rec.reg_idx as usize).expect("Missing contig").contig();

			if curr_contig.as_ref().map(|&s| s != chrom).unwrap_or(true) {
				curr_contig = Some(chrom);
				info!("At contig {}", chrom)
			}

			if var.n() >= min_samples {
				let sd = var.var().sqrt();
				if sd >= min_sdev {
					if let Some(p) = pv.as_mut() {
						p.print_str(format!("{}\t{}\t{}\t{}\t{:.4}", chrom, rec.pos, rec.pos + if rec.two_base { 2 } else { 1 }, var.n(), sd).as_str())?
					}
					if let Some(d) = distance.as_mut() { d.clear_obs() };
					let mut ovec = if heatmaps.is_some() && var.n() > 1 { Some(Vec::with_capacity(var.n())) } else { None };
					for (i, ct) in rec.counts.iter().enumerate() {
						let Counts{non_converted: a, converted: b} = *ct;
						if a + b >= min_counts {
							stats[i].add_site(a, b);
							if let Some(p) = pv.as_mut() { p.print_value(a, b)? }
							if let Some(d) = distance.as_mut() { d.add_obs(i, ((a + 1) as f64) / ((a + b + 2) as f64)) };
							if let Some(hst) = hist.as_mut() { hst.add_obs(i, a, b) }
							if let Some(v) = ovec.as_mut()  { v.push(Obs{idx: i, a, b})}
						} else if let Some(p) = pv.as_mut() {
							p.print_missing()?
						}
					}
					if let Some(p) = pv.as_mut() { p.print_str("\n")? }
					site_stats.add_sdev(sd);
					if let Some(v) = ovec.take() { heatmaps.as_mut().unwrap().add_vec_obs(v) }
				}
			}
		}
		if let Some(tr) = smooth_tr.as_mut() {
			tr.send(mblk).with_context(|| "Error sending message block to smooth writer thread")?;
		}
	}
	info!("Finishing up");

	if cfg.summary() {
	// Write out sample summary statistics
		let fname = cfg.mk_path("sample_stats.txt", false, false);
		debug!("Writing out sample statistics to {:?}", fname);
		let mut summ_wrt = BufWriter::new(File::create(fname)?);

		writeln!(summ_wrt, "barcode\tsample\tn_sites_passed\tmeth_mean\tmeth_var\tcounts:min\tcounts:Q1\tcounts:Q2\tcounts:Q3\tcounts:max")?;
		for st in stats.iter() {
			writeln!(summ_wrt, "{}", st)?;
		}

		if samples.len() > 1 {
			// Write out site summary statistics
			let fname = cfg.mk_path("site_stats.txt", false, false);
			debug!("Writing out site statistics to {:?}", fname);
			let mut summ_wrt = BufWriter::new(File::create(fname)?);
			site_stats.write_report(&mut summ_wrt)?;
		}
	}

	// Write out similarity matrix if required
	if let Some(dist) = distance {
		let fname = cfg.mk_path(match cfg.similarity_type().unwrap() {
			SimilarityType::Distance => "distance.txt",
			SimilarityType::Correlation => "corr.txt",
			SimilarityType::Covariance => "covar.txt",
		}, false, false);
		let mut summ_wrt = BufWriter::new(File::create(fname)?);
		dist.print_table(&mut summ_wrt, samples)?;
	}
	
	// Write out histograms if required

	if let Some(mut hst) = hist {
		hst.handle_cache();
		let fname = cfg.mk_path("dist.txt", false, false);
		let mut summ_wrt = BufWriter::new(File::create(fname)?);
		hst.print_table(&mut summ_wrt, samples)?;
	}
	
	// Write out heatmap data if required

	if let Some(mut hm) = heatmaps {
		if samples.len() > 1 {
			let fname = cfg.mk_path("heatmaps.txt", false, false);
			let mut summ_wrt = BufWriter::new(File::create(fname)?);
			hm.print_table(&mut summ_wrt, samples, cfg.blank_lines())?;
		}
	}
	Ok(())
}