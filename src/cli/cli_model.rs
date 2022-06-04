use clap::{
   Command, Arg, ArgMatches, crate_version,
};

pub(super) fn cli_model() -> ArgMatches {
   Command::new("bed_methyl").version(crate_version!()).author("Simon Heath")
      .about("bedmethyl is a toolkit for handling bedMethyl files (ENCODE WGBS output files)")
      .arg(
         Arg::new("quiet")
            .short('q').long("quiet")
            .help("Silence all output")
            .conflicts_with("loglevel")
      )
      .arg(
         Arg::new("timestamp")
            .long("timestamp")
            .takes_value(true).value_name("GRANULARITY")
            .possible_values(&["none", "error", "warn", "info", "debug", "trace"])
            .ignore_case(true).default_value("none")
            .help("Prepend log entries with a timestamp")
      )
      .arg(
         Arg::new("loglevel")
            .short('l').long("loglevel")
            .takes_value(true).value_name("LOGLEVEL")
            .possible_values(&["none", "error", "warn", "info", "debug", "trace"])
            .ignore_case(true).default_value("info")
            .help("Set log level")
      )
      .arg(
         Arg::new("region")
            .short('r').long("region").global(true)
            .takes_value(true).value_name("REGION").multiple_occurrences(true)
            .help("Genomic region")
      )
      .arg(
         Arg::new("regions")
            .short('R').long("regions").global(true)
            .takes_value(true).value_name("FILE")
            .conflicts_with("region")
            .help("BED file with genomic regions")
      )
      .arg(
         Arg::new("min_counts")
            .short('n').long("min-counts").global(true)
            .takes_value(true).value_name("INT").default_value("1")
            .help("Minimum counts for a sample to pass filter")
      )
      .arg(
         Arg::new("min_samples")
            .short('m').long("min-samples").global(true)
            .takes_value(true).value_name("INT").default_value("1")
            .help("Minimum samples passing filter for entry to be output")
      )
      .arg(
         Arg::new("min_pcnt")
            .short('p').long("min-pcnt").global(true)
            .takes_value(true).value_name("NUMBER").default_value("0.0")
            .help("Minimum % samples passing filter for entry to be output")
      )
      .arg(
         Arg::new("min_sdev")
            .short('v').long("min-sdev").global(true)
            .takes_value(true).value_name("NUMBER").default_value("0.0")
            .help("Minimum methylation std. dev. in samples for entry to be output")
      )
      .arg(
         Arg::new("combine_cpgs")
            .short('c').long("combine-cpgs").global(true)
            .help("Combine +/- strand entries for CpGs [default]")
      )
      .arg(
         Arg::new("no_combine_cpgs")
            .short('C').long("no-combine-cpgs").global(true)
            .conflicts_with("combine_cpgs")
            .help("Do not combine +/- strand entries for CpGs")
      )
      .arg(
         Arg::new("assume_cpg")
            .long("assume-cpg").global(true)
            .help("Assume all records are in CpG context")
      )
      .arg(
         Arg::new("threads")
            .short('t').long("threads").global(true)
            .takes_value(true).value_name("INT")
            .help("Maximum number of read threads [default: number of available threads]")
      )
      .arg(
         Arg::new("prefix")
            .short('P').long("prefix").global(true)
            .takes_value(true).value_name("STRING")
            .default_value("bedmethyl")
            .help("Set prefix for output files")
      )
      .arg(
         Arg::new("dir")
            .short('P').long("dir").global(true)
            .takes_value(true).value_name("DIR")
            .help("Set output directory [default: current directory]")
      )
      .arg(
         Arg::new("bins")
            .short('b').long("bins").global(true)
            .takes_value(true).value_name("INT").default_value("100")
            .help("Number of bins for methylation distributions")
      )
      .arg(
         Arg::new("kde")
            .short('k').long("kde").global(true)
            .help("Use kernel density estimation for distributions [default]")
      )
      .arg(
         Arg::new("no_kde")
            .short('K').long("no-kde").global(true)
            .conflicts_with("kde")
            .help("Do not use kernel density estimation for distributions")
      )
      .arg(
         Arg::new("kde_cache_limit")
            .long("kde-cache-limit").global(true).hide(true)
            .takes_value(true).value_name("INT").default_value("512")
            .help("kde cache usage where number of counts <= kde_cache_limit")
      )
      .subcommand(
         Command::new("merge")
            .about("Generate methylation count data from multiple bedMethyl input files")
            .arg(
               Arg::new("output_names")
                  .short('o').long("output-names")
                  .takes_value(true).value_name("FILE")
                  .multiple_values(true).require_value_delimiter(true)
                  .max_values(2)
                  .help("Set merged output files names")
            )
            .arg(
               Arg::new("output_type")
                  .short('T').long("output-type")
                  .takes_value(true).value_name("OUTPUT TYPE")
                  .possible_values(&["non_conv-conv", "non_conv-cov", "meth-cov", "meth"])
                  .ignore_case(true).default_value("non_conv-conv")
                  .help("Set output type")
            )
            // tab, space, comma, semicolon, slash
            .arg(
               Arg::new("value_delim")
                  .short('D').long("value-delim")
                  .takes_value(true).value_name("DELIM")
                  .possible_values(&["tab", "space", "comma", "semicolon", "slash"])
                  .ignore_case(true).default_value("tab")
                  .help("Set value delimiter")
            )
            .arg(
               Arg::new("split_output")
                  .short('L').long("split-output")
                  .conflicts_with("value_delim")
                  .help("Generate two output files for the two values per sample")
            )
            .arg(
               Arg::new("compress")
                  .short('z').long("compress")
                  .help("Compress output files (with bgzip)")
            )
            .arg(
               Arg::new("summary")
                  .short('s').long("summary")
                  .help("Generate summary files [default]")
            )
            .arg(
               Arg::new("no_summary")
                  .short('S').long("no_summary")
                  .conflicts_with("summary")
                  .help("Do not generate summary files")
            )
            .arg(
               Arg::new("stype")
                  .short('X').long("similarity-type")
                  .takes_value(true).value_name("TYPE")
                  .possible_values(&["correlation", "covariance", "distance"])
                  .ignore_case(true).default_value("correlation")
                  .help("Set similarity type")
            )
            .arg(
               Arg::new("distribution")
                  .long("distribution")
                  .help("Generate per sample methylation distributions")
            )
            .arg(
               Arg::new("heatmaps")
                  .long("heatmaps")
                  .help("Generate joint methylation distribution heatmaps per sample pair")
            )
            .arg(
               Arg::new("full_limit")
                  .long("full-limit").hide(true)
                  .takes_value(true).value_name("INT").default_value("127")
                  .help("Cache usage pairs where individual counts <= full_limit")
            )
            .arg(
               Arg::new("blank_lines")
                  .long("blank_lines")
                  .help("Generate blank lines between iso-lines [default]")
            )
            .arg(
               Arg::new("no_blank_lines")
                  .short('B').long("no_blank_lines")
                  .help("Do not generate blank lines between iso-lines")
            )
            .arg(
               Arg::new("smooth")
                  .long("smooth")
                  .help("Generate smoothed counts using local regression")
            )
            .arg(
               Arg::new("window_size")
                  .short('w').long("window-size")
                  .takes_value(true).value_name("BASE PAIRS").default_value("1000")
                  .help("Minimum window size for smoothing")
            )
            .arg(
               Arg::new("min_sites")
                  .long("min-sites")
                  .takes_value(true).value_name("INT").default_value("70")
                  .help("Minimum number of observed sites in window for smoothing")
            )
            .arg(
               Arg::new("max_distance")
                  .long("max-distance")
                  .takes_value(true).value_name("BASE PAIRS").default_value("50000")
                  .help("Maximum distance between adjacent observed sites for smoothing (0 = no limit)")
            )
            .arg(
               Arg::new("smooth_output_type")
                  .long("smooth-output-type")
                  .takes_value(true).value_name("OUTPUT TYPE")
                  .possible_values(&["non_conv-conv", "non_conv-cov", "meth-cov", "meth"])
                  .ignore_case(true)
                  .help("Set output type")
            )
            .arg(
               Arg::new("sample_file")
                  .takes_value(true).value_name("FILE")
                  .help("Sample input file [default: <stdin>]")
            )
      )
      .subcommand(
         Command::new("summary")
            .about("Generate summary files from multiple BedMethyl input files")
            .arg(
               Arg::new("sample_file")
                  .takes_value(true).value_name("FILE")
                  .help("Sample input file [default: <stdin>]")
            )
      )
      .subcommand(
         Command::new("similarity")
            .about("Generate similarity matrix from multiple BedMethyl input files")
            .arg(
               Arg::new("stype")
                  .short('X').long("similarity-type")
                  .takes_value(true).value_name("TYPE")
                  .possible_values(&["correlation", "covariance", "distance"])
                  .ignore_case(true).default_value("correlation")
                  .help("Set similarity type")
            )
            .arg(
               Arg::new("sample_file")
                  .takes_value(true).value_name("FILE")
                  .help("Sample input file [default: <stdin>]")
            )
      )
      .subcommand(
         Command::new("distribution")
            .about("Generate per sample methylation distributions from multiple BedMethyl input files")
            .arg(
               Arg::new("sample_file")
                  .takes_value(true).value_name("FILE")
                  .help("Sample input file [default: <stdin>]")
            )
      )
      .subcommand(
         Command::new("heatmaps")
            .about("Generate heatmaps showing the joint methylation distributions for all sample pairs from multiple BedMethyl input files")
            .arg(
               Arg::new("full_limit")
                  .long("full-limit").hide(true)
                  .takes_value(true).value_name("INT").default_value("127")
                  .help("Cache usage pairs where individual counts <= full_limit")
            )
            .arg(
               Arg::new("blank_lines")
                  .long("blank_lines")
                  .help("Generate blank lines between iso-lines [default]")
            )
            .arg(
               Arg::new("no_blank_lines")
                  .short('B').long("no_blank_lines")
                  .help("Do not generate blank lines between iso-lines")
            )
            .arg(
               Arg::new("sample_file")
                  .takes_value(true).value_name("FILE")
                  .help("Sample input file [default: <stdin>]")
            )
      )
      .subcommand(
         Command::new("smooth")
            .about("Generate smoothed counts using local regression")
            .arg(
               Arg::new("window_size")
                  .long("window-size")
                  .takes_value(true).value_name("BASE PAIRS").default_value("1000")
                  .help("Minimum window size for smoothing")
            )
            .arg(
               Arg::new("min_sites")
                  .long("min-sites")
                  .takes_value(true).value_name("INT").default_value("70")
                  .help("Minimum number of observed sites in window for smoothing")
            )
            .arg(
               Arg::new("max_distance")
                  .long("max-distance")
                  .takes_value(true).value_name("BASE PAIRS").default_value("50000")
                  .help("Maximum distance between adjacent observed sites for smoothing (0 = no_limit)")
            )
            .arg(
               Arg::new("smooth_output_type")
                  .short('T').long("smooth-output-type")
                  .takes_value(true).value_name("OUTPUT TYPE")
                  .possible_values(&["non_conv-conv", "non_conv-cov", "meth-cov", "meth"])
                  .ignore_case(true).default_value("non_conv-conv")
                  .help("Set output type")
            )
            // tab, space, comma, semicolon, slash
            .arg(
               Arg::new("value_delim")
                  .short('D').long("value-delim")
                  .takes_value(true).value_name("DELIM")
                  .possible_values(&["tab", "space", "comma", "semicolon", "slash"])
                  .ignore_case(true).default_value("tab")
                  .help("Set value delimiter")
            )
            .arg(
               Arg::new("round_counts")
                  .long("round-counts")
                  .help("Round smoothed counts")
            )
            .arg(
               Arg::new("compress")
                  .short('z').long("compress")
                  .help("Compress output files (with bgzip)")
            )
            .arg(
               Arg::new("sample_file")
                  .takes_value(true).value_name("FILE")
                  .help("Sample input file [default: <stdin>]")
            )
      )
      .get_matches()
}