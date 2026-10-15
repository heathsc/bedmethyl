mod cli_model;
mod log_level;
mod log_utils;
pub mod output_type;
pub mod similarity_type;
pub mod value_delim;

use std::{num::NonZeroUsize, path::PathBuf};

use clap::ArgMatches;
use log::Level::Debug;

use r_htslib::*;

use crate::{
    config::*,
    regions::{Regions, Sites},
    sample::self,
};

use output_type::OutputType;
use similarity_type::SimilarityType;
use value_delim::ValueDelim;

fn check_percent(x: f64) -> anyhow::Result<f64> {
    if (0.0..=100.0).contains(&x) {
        Ok(x)
    } else {
        Err(anyhow!(
            "Invalid argument to parse_percent:  '{x}' out of range"
        ))
    }
}

// Create Config structure and handle options in common between the sub-commands
fn handle_common(m: &ArgMatches, msub: &ArgMatches) -> anyhow::Result<(Config, Vec<Hts>)> {
    // Parse sample file
    let sample_info = sample::read_sample_file(msub.get_one::<PathBuf>("sample_file"))?;

    let nt = m
        .get_one::<u32>("threads")
        .and_then(|s| NonZeroUsize::new(*s as usize));

    let hts_vec = sample::open_input_files(
        sample_info.samples(),
        nt.map(usize::from).unwrap_or_else(num_cpus::get),
    )?;

    // Parse region options
    let region = if let Some(f) = m.get_one::<PathBuf>("regions") {
        Some(Regions::from_file(f)?)
    } else {
        m.get_many::<String>("region").map(Regions::from_regions)
    }
    .unwrap_or_else(|| {
        let mut ctgs = sample::get_contig_list(&hts_vec);
        ctgs.sort();
        Regions::from_str_vec(&ctgs)
    });

    let mut cfg = Config::new(sample_info, region)?;

    if let Some(n) = nt {
        cfg.set_threads(n)
    }

    cfg.set_min_counts(m.get_one::<u32>("min_counts").map(|x| *x).unwrap());
    cfg.set_min_sdev(m.get_one::<f64>("min_sdev").copied().unwrap());

    if let Some(s) = m.get_one::<u32>("min_samples") {
        cfg.set_min_samples(NonZeroUsize::new(*s as usize).unwrap())
    } else if let Some(s) = m.get_one::<f64>("min_pcnt") {
        let p = check_percent(*s)?;
        cfg.set_min_samples(
            NonZeroUsize::new(((p * 0.01 * (cfg.n_samples() as f64)).round() as usize).max(1))
                .unwrap(),
        )
    }

    if let Some(s) = m.get_one::<PathBuf>("dir") {
        cfg.set_dir(s)
    }

    cfg.set_prefix(m.get_one::<String>("prefix").unwrap());

    cfg.set_combine_cpgs(!m.get_flag("no_combine_cpgs"));

    cfg.set_assume_cpg(!m.get_flag("assume_cpg"));

    cfg.set_kde(!m.get_flag("no_kde"));

    cfg.set_n_bins(
        m.get_one::<u32>("bins")
            .map(|x| NonZeroUsize::new(*x as usize).unwrap())
            .unwrap(),
    );

    cfg.set_kde_cache_limit(m.get_one::<u16>("kde_cache_limit").copied().unwrap());

    Ok((cfg, hts_vec))
}

// Handle region/site options
fn handle_sites(cfg: &mut Config, m: &ArgMatches) -> anyhow::Result<()> {
    // Parse site options
    if let Some(f) = m.get_one::<PathBuf>("sites") {
        let sites = Sites::from_file(f)?;
        let win_size = if cfg.smooth() {
            cfg.min_sites() * cfg.max_distance()
        } else {
            10000
        };
        cfg.regions_mut().trim(&sites, win_size);
        cfg.set_specific_sites(true);
    }
    Ok(())
}

// Handle options for merge sub-command
fn merge_com(m: &ArgMatches, msub: &ArgMatches) -> anyhow::Result<(Config, Vec<Hts>)> {
    let (mut cfg, hts_vec) = handle_common(m, msub)?;

    // Handle output options
    let value_delim = if msub.get_flag("split_output") {
        Some(ValueDelim::Split)
    } else {
        msub.get_one::<ValueDelim>("value_delim").copied()
    };

    let output_type = msub.get_one::<OutputType>("output_type").copied();

    cfg.set_merge_outputs(
        output_type,
        value_delim,
        msub.get_many::<PathBuf>("output_names"),
    );

    cfg.set_compress(msub.get_flag("compress"));

    cfg.set_summary(!msub.get_flag("no_summary"));

    if msub.get_flag("similarity") {
        handle_similarity_opts(&mut cfg, msub)
    }

    cfg.set_distribution(msub.get_flag("distribution"));

    if msub.get_flag("heatmaps") {
        handle_heatmaps_opts(&mut cfg, msub)?
    }

    if msub.get_flag("smooth") {
        handle_smooth_opts(&mut cfg, msub, true)?
    }

    Ok((cfg, hts_vec))
}

fn summary_com(m: &ArgMatches, msub: &ArgMatches) -> anyhow::Result<(Config, Vec<Hts>)> {
    let (mut cfg, hts_vec) = handle_common(m, msub)?;

    cfg.set_summary(true);

    Ok((cfg, hts_vec))
}

fn handle_similarity_opts(cfg: &mut Config, msub: &ArgMatches) {
    // Handle similarity type
    let similarity_type = if let Some(s) = msub.get_one::<SimilarityType>("stype") {
        *s
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
    cfg.set_heatmaps_full_limit(*msub.get_one("full_limit").unwrap());
    cfg.set_blank_lines(!msub.get_flag("no_blank_lines"));

    cfg.set_heatmaps(true);

    Ok(())
}

fn handle_smooth_opts(cfg: &mut Config, msub: &ArgMatches, output: bool) -> anyhow::Result<()> {
    // Handle heatmaps options
    cfg.set_window_size(
        msub.get_one::<u32>("window_size")
            .and_then(|x| NonZeroUsize::new(*x as usize))
            .unwrap(),
    );
    cfg.set_min_sites(
        msub.get_one::<u32>("min_sites")
            .and_then(|x| NonZeroUsize::new(*x as usize))
            .unwrap(),
    );
    cfg.set_max_distance(*msub.get_one::<usize>("max_distance").unwrap());
    cfg.set_round_counts(msub.get_flag("round_counts"));
    cfg.set_compress(msub.get_flag("compress"));

    if output {
        let value_delim = if msub.get_flag("split_output") {
            Some(ValueDelim::Split)
        } else {
            msub.get_one::<ValueDelim>("value_delim").copied()
        };

        let output_type = msub
            .get_one::<OutputType>("smooth_output_type")
            .copied()
            .or_else(|| msub.get_one::<OutputType>("output_type").copied())
            .or(Some(OutputType::NonconvConv));

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

    if !msub.get_flag("no-smooth") {
        handle_smooth_opts(&mut cfg, msub, false)?;
    }
    cfg.set_model(true);

    Ok((cfg, hts_vec))
}

pub fn handle_cli() -> anyhow::Result<(Config, Vec<Hts>)> {
    let m = cli_model::cli_model().get_matches();
    log_utils::init_log(&m);

    debug!("Handling command line inputs");

    let (mut cfg, hts_vec) = match m.subcommand() {
        Some(("merge", m_sub)) => merge_com(&m, m_sub),
        Some(("summary", m_sub)) => summary_com(&m, m_sub),
        Some(("similarity", m_sub)) => similarity_com(&m, m_sub),
        Some(("distribution", m_sub)) => distribution_com(&m, m_sub),
        Some(("heatmaps", m_sub)) => heatmaps_com(&m, m_sub),
        Some(("smooth", m_sub)) => smooth_com(&m, m_sub),
        Some(("model", m_sub)) => model_com(&m, m_sub),
        _ => Err(anyhow!("Unknown subcommand")),
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
            let x = reg
                .sites()
                .map(|s| format!("Sites: {}", s.len()))
                .unwrap_or_default();
            debug!("{}\t{}", reg, x)
        }
    }
    info!("  Output compression: {}", cfg.compress());
    info!("  Prefix set to {}", cfg.prefix());
    if let Some(d) = cfg.dir() {
        info!("  Dir set to {:?}", d)
    };
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
    if cfg.summary() {
        info!("  Creating summary files")
    }
    if cfg.distribution() {
        info!("  Creating methylation distribution file");
        info!("    Number of bins: {}", cfg.n_bins());
        info!("    Kernel density estimation: {:?}", cfg.kde());
        if cfg.kde() {
            debug!("   kde_cache_limit: {}", cfg.kde_cache_limit())
        }
    }
    if let Some(stype) = cfg.similarity_type() {
        info!("  Creating similarity matrix");
        info!("  Similarity type: {:?}", stype);
    }

    Ok((cfg, hts_vec))
}
