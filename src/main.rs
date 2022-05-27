#[macro_use]
extern crate log;
#[macro_use]
extern crate anyhow;

mod cli;
mod config;
mod regions;
mod log_utils;
mod stats;
mod process;
mod sample;

fn main() -> anyhow::Result<()> {
	let (cfg, hts_vec) = cli::handle_cli()?;
	process::process(&cfg, hts_vec)
}
