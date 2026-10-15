use clap::ArgMatches;

use super::log_level::LogLevel;

/// Set up stderrlog using options from clap::ArgMatches
/// 
/// Relevant options:
/// 
///   loglevel: can be error, warn, info, trace or none.  Default is info
///   quiet: turns off all loggin irrespective of loglevel
///   timestamp: sets timestamp options for logger
/// 
///   loglevel none is a synonym for the quiet option
/// 
pub fn init_log(m: &ArgMatches) {
    let verbose = m
        .get_one::<LogLevel>("loglevel")
        .copied()
        .expect("Missing default log level");
	
    let quiet = verbose.is_none() || m.get_flag("quiet");
    let ts = m
        .get_one::<stderrlog::Timestamp>("timestamp")
        .copied()
        .unwrap_or(stderrlog::Timestamp::Off);

    stderrlog::new()
        .quiet(quiet)
        .verbosity(verbose.get_level())
        .timestamp(ts)
        .init()
        .unwrap();
}
