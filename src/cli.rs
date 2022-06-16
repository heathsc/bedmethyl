mod cli_model;

use std::{
	num::NonZeroUsize,
	str::FromStr,
	fmt::Display,
};

use anyhow::{Context, Error};
use clap::ArgMatches;
use log::Level::Debug;

use r_htslib::*;

use crate::{
	config::*,
	regions::{Regions, Sites},
	sample
};

fn parse<T>(s: &str) -> anyhow::Result<T>
where
	T: FromStr,
	T::Err: Display + Send + Sync + 'static,
{
	s.parse::<T>().map_err(|e| Error::msg(format!("Parse error {}", e)))
}

fn parse_percent<S: AsRef<str>>(s: S) -> anyhow::Result<f64> {
	let x = s.as_ref().parse::<f64>().with_context(|| "Invalid argument to parse_percent")?;
	if (0.0..=100.0).contains(&x) {
		Ok(x)
	} else {
		Err(anyhow!("Invalid argument to parse_percent:  '{}' out of range", s.as_ref()))
	}
}

// Create Config structure and handle options in common between the sub-commands
fn handle_common(m: &ArgMatches, msub: &ArgMatches) -> anyhow::Result<(Config, Vec<Hts>)> {

	// Parse sample file
	let sample_info = sample::read_sample_file(msub.value_of("sample_file"))?;

	let nt = if let Some(s) = m.value_of("threads") {
		Some(parse::<NonZeroUsize>(s)?)
	} else { None };

	let hts_vec = sample::open_input_files(sample_info.samples(), nt.map(usize::from).unwrap_or_else(num_cpus::get))?;

	// Parse region options
	let region = if let Some(f)= m.value_of("regions") {
		Some(Regions::from_file(f)?)
	} else {
		m.values_of("region").map(Regions::from_regions)
	}.unwrap_or_else(|| {
		let mut ctgs = sample::get_contig_list(&hts_vec);
		ctgs.sort();
		Regions::from_str_vec(&ctgs)
	});
	
	let mut cfg = Config::new( sample_info, region)?;

	if let Some(n) = nt { cfg.set_threads(n) }

	if let Some(s) = m.value_of("min_counts") {
		cfg.set_min_counts(parse::<u32>(s)?)
	}

	if let Some(s) = m.value_of("min_sdev") {
		cfg.set_min_sdev(parse::<f64>(s)?)
	}
	
	if let Some(s) = m.value_of("min_samples") {
		cfg.set_min_samples(parse::<NonZeroUsize>(s)?)
	} else if let Some(s) = m.value_of("min_pcnt") {
		let p = parse_percent(s)?;
		cfg.set_min_samples(NonZeroUsize::new(((p * 0.01 * (cfg.n_samples() as f64)).round() as usize).max(1)).unwrap())
	} 

	if let Some(s) = m.value_of("dir") {
		cfg.set_dir(s)
	}		

	cfg.set_prefix(m.value_of("prefix").unwrap());

	cfg.set_combine_cpgs(!m.is_present("no_combine_cpgs"));

	cfg.set_assume_cpg(!m.is_present("assume_cpg"));

	cfg.set_kde(!m.is_present("no_kde"));
	if let Some(s) = m.value_of("bins") {
		let nb = parse::<NonZeroUsize>(s)?;
		cfg.set_n_bins(nb)
	}
	cfg.set_kde_cache_limit(parse::<u16>(m.value_of("kde_cache_limit").unwrap())?);

	Ok((cfg, hts_vec))
}

// Handle region/site options
fn handle_sites(cfg: &mut Config, m: &ArgMatches) -> anyhow::Result<()> {
	// Parse site options
	if let Some(f) = m.value_of("sites") {
		let sites = Sites::from_file(f)?;
		let win_size = if cfg.smooth() {
			cfg.min_sites() * cfg.max_distance()
		} else { 10000 };
		cfg.regions_mut().trim(&sites, win_size);
		cfg.set_specific_sites(true);
	}
	Ok(())
}

// Handle options for merge sub-command
fn merge_com(m: &ArgMatches, msub: &ArgMatches) -> anyhow::Result<(Config, Vec<Hts>)> {

	let (mut cfg, hts_vec) = handle_common(m, msub)?;

	// Handle output options
	let value_delim = if msub.is_present("split_output") {
		Some(ValueDelim::Split)
	} else {
		msub.value_of("value_delim").map(ValueDelim::from_str).flatten()
	};

	let output_type = msub.value_of("output_type").map(OutputType::from_str).flatten();

	cfg.set_merge_outputs(output_type, value_delim, msub.values_of("output_names"));

	cfg.set_compress(msub.is_present("compress"));

	cfg.set_summary(!msub.is_present("no_summary"));

	if msub.is_present("similarity") {
		handle_similarity_opts(&mut cfg, msub)
	}

	cfg.set_distribution(msub.is_present("distribution"));

	if msub.is_present("heatmaps") {
		handle_heatmaps_opts(&mut cfg, msub)?
	}

	if msub.is_present("smooth") {
		handle_smooth_opts(&mut cfg, msub, true)?
	} else {

	}

	Ok((cfg, hts_vec))
}

fn summary_com(m: &ArgMatches, msub: &ArgMatches) -> anyhow::Result<(Config, Vec<Hts>)> {

	let (mut cfg, hts_vec)  = handle_common(m, msub)?;
	
	cfg.set_summary(true);
	
	Ok((cfg, hts_vec))
}

fn handle_similarity_opts(cfg: &mut Config, msub: &ArgMatches) {

	// Handle similarity type
	let similarity_type = if let Some(s) = msub.value_of("type") {
		SimilarityType::from_str(s).expect("Unknown Similarity Type")
	} else {
		SimilarityType::Correlation
	};
	cfg.set_similarity_type(similarity_type);	
}


fn similarity_com(m: &ArgMatches, msub: &ArgMatches) -> anyhow::Result<(Config, Vec<Hts>)> {

	let (mut cfg, hts_vec) = handle_common(m, msub)?;

	handle_similarity_opts(&mut cfg, msub);	

	Ok((cfg, hts_vec))
}

fn distribution_com(m: &ArgMatches, msub: &ArgMatches) -> anyhow::Result<(Config, Vec<Hts>)> {

	let (mut cfg, hts_vec) = handle_common(m, msub)?;
	
	cfg.set_distribution(true);
	
	Ok((cfg, hts_vec))
}

fn handle_heatmaps_opts(cfg: &mut Config, msub: &ArgMatches) -> anyhow::Result<()> {

	// Handle heatmaps options
	cfg.set_heatmaps_full_limit(parse::<u8>(msub.value_of("full_limit").unwrap())?);
	cfg.set_blank_lines(!msub.is_present("no_blank_lines"));

	cfg.set_heatmaps(true);
	
	Ok(())
}

fn handle_smooth_opts(cfg: &mut Config, msub: &ArgMatches, output: bool) -> anyhow::Result<()> {

	// Handle heatmaps options
	cfg.set_window_size(parse::<NonZeroUsize>(msub.value_of("window_size").unwrap())?);
	cfg.set_min_sites(parse::<NonZeroUsize>(msub.value_of("min_sites").unwrap())?);
	cfg.set_max_distance(parse::<usize>(msub.value_of("max_distance").unwrap())?);
	cfg.set_round_counts(msub.is_present("round_counts"));
	cfg.set_compress(msub.is_present("compress"));

	if output {
		let value_delim = if msub.is_present("split_output") {
			Some(ValueDelim::Split)
		} else {
			msub.value_of("value_delim").map(ValueDelim::from_str).flatten()
		};

		let output_type = msub.value_of("smooth_output_type")
			.or_else(|| msub.value_of("output_type"))
			.or(Some("non_conv-meth"))
			.and_then(OutputType::from_str);

		cfg.set_smooth_outputs(output_type, value_delim);
	}

	cfg.set_smooth(true);

	Ok(())
}

fn heatmaps_com(m: &ArgMatches, msub: &ArgMatches) -> anyhow::Result<(Config, Vec<Hts>)> {

	let (mut cfg, hts_vec) = handle_common(m, msub)?;
	
	handle_heatmaps_opts(&mut cfg, msub)?;	
	
	Ok((cfg, hts_vec))
}

fn smooth_com(m: &ArgMatches, msub: &ArgMatches) -> anyhow::Result<(Config, Vec<Hts>)> {

	let (mut cfg, hts_vec) = handle_common(m, msub)?;

	handle_smooth_opts(&mut cfg, msub, true)?;

	Ok((cfg, hts_vec))
}

fn model_com(m: &ArgMatches, msub: &ArgMatches) -> anyhow::Result<(Config, Vec<Hts>)> {

	let (mut cfg, hts_vec) = handle_common(m, msub)?;

	if !msub.is_present("no-smooth") {
		handle_smooth_opts(&mut cfg, msub, false)?;
	}
	cfg.set_model(true);

	Ok((cfg, hts_vec))
}

pub fn handle_cli() -> anyhow::Result<(Config, Vec<Hts>)> {

	let m = cli_model::cli_model();
	super::log_utils::init_log(&m);
	
	debug!("Handling command line inputs");

	let (mut cfg, hts_vec) = match m.subcommand() {
		Some(("merge", m_sub)) => {
			merge_com(&m, m_sub)
		},
		Some(("summary", m_sub)) => {
			summary_com(&m, m_sub)
		},
		Some(("similarity", m_sub)) => {
			similarity_com(&m, m_sub)
		},
		Some(("distribution", m_sub)) => {
			distribution_com(&m, m_sub)
		},
		Some(("heatmaps", m_sub)) => {
			heatmaps_com(&m, m_sub)
		},
		Some(("smooth", m_sub)) => {
			smooth_com(&m, m_sub)
		},
		Some(("model", m_sub)) => {
			model_com(&m, m_sub)
		},
		_ => {
			Err(anyhow!("Unknown subcommand"))
		},
	}?;

	handle_sites(&mut cfg, &m)?;

	info!("Configuration");
	info!("  Number of input files: {}", cfg.n_samples());
	let si = cfg.sample_info();
	if !si.factors().is_empty() {
		info!("  Number of factors: {}", si.factors().len());
		for f in si.factors().iter() {
			info!("    Name: {}, No. levels: {}", f.name(), f.levels().len());
		}
	}
	if !si.covariates().is_empty() {
		info!("  Number of covariates: {}", si.covariates().len());
		for c in si.covariates().iter() {
			info!("    Name: {}", c.name());
		}
	}
	info!("  No. Regions: {}", cfg.regions().n_regions());
	if log_enabled!(Debug) {
		for reg in cfg.regions().iter() {
			let x =reg.sites()
				.map(|s| format!("Sites: {}", s.len()))
				.unwrap_or_default();
			debug!("{}\t{}", reg, x)
		}
	}
	info!("  Output compression: {}", cfg.compress());
	info!("  Prefix set to {}", cfg.prefix());
	if let Some(d) = cfg.dir() { info!("  Dir set to {:?}", d) };
	info!("  Minimum counts set to {}", cfg.min_counts());	
	info!("  Minimum samples set to {}", cfg.min_samples());	
	info!("  Minimum methylation stdev set to {}", cfg.min_sdev());	
	info!("  Combined CpGs: {}", cfg.combine_cpgs());
	if let Some(mo) = cfg.merge_output() {
		info!("  Creating merged output");
		info!("    Output Type: {:?}", mo.output_type());
		info!("    Value Delim: {:?}", mo.value_delim());
		for p in cfg.outputs(false).iter() {
			info!("    Output file: {:?}", p)
		}
	}
	if let Some(mo) = cfg.smooth_output() {
		info!("  Creating smoothed output");
		info!("    Output Type: {:?}", mo.output_type());
		info!("    Value Delim: {:?}", mo.value_delim());
		for p in cfg.outputs(true).iter() {
			info!("    Output file: {:?}", p)
		}
	}
	if cfg.summary() { info!("  Creating summary files") }
	if cfg.distribution() {
		info!("  Creating methylation distribution file");
		info!("    Number of bins: {}", cfg.n_bins());
		info!("    Kernel density estimation: {:?}", cfg.kde());
		if cfg.kde() { debug!("   kde_cache_limit: {}", cfg.kde_cache_limit()) }
	}
	if let Some(stype) = cfg.similarity_type() {
		info!("  Creating similarity matrix");
		info!("  Similarity type: {:?}", stype);
	}
	
	Ok((cfg, hts_vec))
}
