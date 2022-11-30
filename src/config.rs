use std::{
	io,
	path::{Path, PathBuf},
	num::NonZeroUsize,
};

use crate::{
	regions::Regions,
	sample::SampleInfo,
};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum OutputType {
	NonconvConv,
	NonconvCov,
	MethCov,
	Meth
}

impl OutputType {
	pub fn from_str(s: &str) -> Option<Self> {
		match s.to_ascii_lowercase().as_str() {
			"non_conv-conv" => Some(Self::NonconvConv),
			"non_conv-cov" => Some(Self::NonconvCov),
			"meth-cov" => Some(Self::MethCov),
			"meth" => Some(Self::Meth),
			_ => None,
		}
	}	
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum ValueDelim {
	Tab,
	Space,
	Semicolon,
	Comma,
	Slash,
	Split,
	None,
}

impl ValueDelim {
	pub fn from_str(s: &str) -> Option<Self> {
		match s.to_ascii_lowercase().as_str() {
			"tab" => Some(Self::Tab),
			"space" => Some(Self::Space),
			"comma" => Some(Self::Comma),
			"semicolon" => Some(Self::Semicolon),
			"slash" => Some(Self::Slash),
			_ => None,
		}
	}

	pub fn get_delim_char(&self) -> char {
		match self {
			ValueDelim::Space => ' ',
			ValueDelim::Comma => ',',
			ValueDelim::Semicolon => ';',
			_ => '\t',
		}
	}
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum SimilarityType {
	Correlation,
	Covariance,
	Distance
}

impl SimilarityType {
	pub fn from_str(s: &str) -> Option<Self> {
		match s.to_ascii_lowercase().as_str() {
			"correlation" => Some(Self::Correlation),
			"covariance" => Some(Self::Covariance),
			"distance" => Some(Self::Distance),
			_ => None,
		}
	}	
}

pub struct MergeOutput {
	outputs: [Option<PathBuf>; 2],
	output_type: OutputType,
	value_delim: ValueDelim,	
}

impl MergeOutput {
	fn new(otype: Option<OutputType>, delim: Option<ValueDelim>) -> Self {
		let output_type = otype.unwrap_or(OutputType::NonconvConv); 
		let mut value_delim = delim.unwrap_or_else(|| {
			if output_type == OutputType::Meth { ValueDelim::None } else { ValueDelim::Tab }	
		}); 
		if output_type == OutputType::Meth && value_delim != ValueDelim::None {
			warn!("value_delim set to None for OutPutType::Meth");
			value_delim = ValueDelim::None
		}
		debug!("Output Type: {:?}, Value Delim: {:?}", output_type, value_delim);
		Self {
			outputs: [None, None],
			output_type,
			value_delim,
		}
	}

	fn set_outputs<P: AsRef<Path>>(&mut self, output1: Option<P>, output2: Option<P>) {
		let output1 = output1.map(|x| x.as_ref().to_owned());
		let mut output2 = output2.map(|x| x.as_ref().to_owned());
		match self.output_type {
			OutputType::Meth => {
				if output2.is_some() { 
					error!("Second output file should not be set for Meth output option");
					output2 = None;
				 }
			},
			_ => {
				if output1.is_none() && output2.is_some() {
					error!("Second output file should not be set if first output file is missing");
					output2 = None
				}
				if output1.is_some() {
					if output2.is_some() { self.value_delim = ValueDelim::Split }
					else if self.value_delim == ValueDelim::Split {
						panic!("To split output files, zero or two output files must be given");						
					}
				} 
			}
		}
		self.outputs = [output1, output2]
	}
	
	pub fn value_delim(&self) -> ValueDelim { self.value_delim }

	pub fn output_type(&self) -> OutputType { self.output_type }
	
	fn output_paths(&self, cfg: &Config, smooth: bool) -> Vec<PathBuf> {
		let mut v = Vec::new();

		// Add directory if present
		let add_dir = |cfg: &Config, s: &PathBuf| {
			if let Some(d) = cfg.dir() {
				d.join(s)
			} else {
				s.clone()
			}
		};
		
		match &self.outputs {
			[Some(s1), Some(s2)] => {
				assert_eq!(self.value_delim, ValueDelim::Split);
				v.push(add_dir(cfg, s1));
				v.push(add_dir(cfg, s2));
			},
			[Some(s1), None] => {	
				assert_ne!(self.value_delim, ValueDelim::Split);
				v.push(add_dir(cfg, s1));
			},	
			// No outputs provided so we build them
			_ => {
				let split = self.value_delim == ValueDelim::Split;
				 
				match self.output_type {
					OutputType::NonconvConv => {
						if split {
							v.push(cfg.mk_path("non_conv.txt", true, smooth));
							v.push(cfg.mk_path("conv.txt", true, smooth));
						} else {
							v.push(cfg.mk_path("nconv_conv.txt", true, smooth));
						}
					},
					OutputType::NonconvCov => {
						if split {
							v.push(cfg.mk_path("non_conv.txt", true, smooth));
							v.push(cfg.mk_path("cov.txt", true, smooth));
						} else {
							v.push(cfg.mk_path("nconv_cov.txt", true, smooth));
						}
					},					
					OutputType::MethCov => {
						if split {
							v.push(cfg.mk_path("meth.txt", true, smooth));
							v.push(cfg.mk_path("cov.txt", true, smooth));
						} else {
							v.push(cfg.mk_path("meth_cov.txt", true, smooth));
						}
					},
					OutputType::Meth => {
						assert!(!split);
						v.push(cfg.mk_path("meth.txt", true, smooth));
					},					
				};
			},
		}
		v
	}
}

#[derive(Default)]
pub struct ConfigCore {
	// Common output options
	prefix: Option<String>,
	dir: Option<PathBuf>,
	
	// Common process options
	threads: Option<NonZeroUsize>,

	// Common selection options
	min_counts: u32,
	min_samples: Option<NonZeroUsize>,
	min_sdev: f64,
	specific_sites: bool,

	// Options for merge subcommand
	merge_output: Option<MergeOutput>,
	combine_cpgs: bool,
	assume_cpg: bool,
	compress: bool,

	// Options for similarity subcommand
	similarity_type: Option<SimilarityType>,

	// Options for summary subcommand
	summary: bool,

	// Options for distribution subcommand
	kde: bool,
	kde_cache_limit: u16, // Limit on a + b
	n_bins: Option<NonZeroUsize>,
	distribution: bool,
	
	// Options for heatmaps subcommand
	heatmaps: bool,
	blank_lines: bool,
	heatmaps_full_limit: u8, // Limit on a or b

	// Options for smoothing
	smooth: bool,
	round_counts: bool,
	window_size: usize,
	min_sites: usize,
	max_distance: usize,
	smooth_output: Option<MergeOutput>,

	// Options for model fitting
	model: bool,
}

pub struct Config {
	core: ConfigCore,
	// Sample information
	sample_info: SampleInfo,
	// Region information
	regions: Regions,
}

impl Config {
	pub fn new(sample_info: SampleInfo, regions: Regions) -> io::Result<Self> {
		let core = ConfigCore::default();
		Ok(Self{core, sample_info, regions})
	}
	
	fn mk_merge_output(&mut self, output_type: Option<OutputType>, value_delim: Option<ValueDelim>) {
		if self.core.merge_output.is_none() {
			self.core.merge_output = Some(MergeOutput::new(output_type, value_delim))
		}	
	}

	pub fn set_smooth_outputs(&mut self, output_type: Option<OutputType>, value_delim: Option<ValueDelim>) {
		if self.core.smooth_output.is_none() {
			self.core.smooth_output = Some(MergeOutput::new(output_type, value_delim))
		}
	}

	pub fn set_merge_outputs<P: AsRef<Path>, I: IntoIterator<Item = P>>(&mut self, output_type: Option<OutputType>, value_delim: Option<ValueDelim>, outputs: Option<I>) {
		self.mk_merge_output(output_type, value_delim);
		let (output1, output2) = if let Some(x) = outputs {
			let mut it = x.into_iter();
			(it.next(), it.next())
		} else {
			(None, None)
		};
		self.core.merge_output.as_mut().unwrap().set_outputs(output1, output2)
	}
	
	pub fn set_min_counts(&mut self, counts: u32) { self.core.min_counts = counts }
	
	pub fn set_min_samples(&mut self, samples: NonZeroUsize) { self.core.min_samples = Some(samples) }

	pub fn set_min_sdev(&mut self, sdev: f64) { self.core.min_sdev = sdev }
	
	pub fn set_threads(&mut self, threads: NonZeroUsize) { self.core.threads = Some(threads) }
	
	pub fn set_similarity_type(&mut self, st: SimilarityType) { self.core.similarity_type = Some(st) }
	
	pub fn set_prefix<S: AsRef<str>>(&mut self, prefix: S) {
		let prefix = prefix.as_ref().to_owned();
		self.core.prefix = Some(prefix)
	}

	pub fn set_dir<P: AsRef<Path>>(&mut self, dir: P) {
		let dir = dir.as_ref().to_owned();
		self.core.dir = Some(dir)
	}

	pub fn set_compress(&mut self, compress: bool) { self.core.compress = compress }

	pub fn set_model(&mut self, x: bool) { self.core.model = x }

	pub fn set_combine_cpgs(&mut self, combine: bool) {	self.core.combine_cpgs = combine }

	pub fn set_assume_cpg(&mut self, assume: bool) {	self.core.assume_cpg = assume }

	pub fn set_summary(&mut self, summary: bool) { self.core.summary = summary }

	pub fn set_kde(&mut self, kde: bool) { self.core.kde = kde }

	pub fn set_distribution(&mut self, x: bool) { self.core.distribution = x }

	pub fn set_heatmaps(&mut self, x: bool) { self.core.heatmaps = x }

	pub fn set_blank_lines(&mut self, x: bool) { self.core.blank_lines = x }

	pub fn set_n_bins(&mut self, n: NonZeroUsize) { self.core.n_bins = Some(n) }

	pub fn set_smooth(&mut self, x: bool) { self.core.smooth = x }

	pub fn set_specific_sites(&mut self, x: bool) { self.core.specific_sites = x }

	pub fn set_round_counts(&mut self, x: bool) { self.core.round_counts = x }

	pub fn set_window_size(&mut self, n: NonZeroUsize) { self.core.window_size = usize::from(n) }

	pub fn set_min_sites(&mut self, n: NonZeroUsize) { self.core.min_sites = usize::from(n) }

	pub fn set_max_distance(&mut self, n: usize) { self.core.max_distance = n }

	pub fn set_kde_cache_limit(&mut self, n: u16) { self.core.kde_cache_limit = n }

	pub fn set_heatmaps_full_limit(&mut self, n: u8) { self.core.heatmaps_full_limit = n }

	pub fn merge_output(&self) -> Option<&MergeOutput> { self.core.merge_output.as_ref() }

	pub fn outputs(&self, smooth: bool) -> Option<Vec<PathBuf>> {
		if smooth {
			self.core.smooth_output.as_ref()
		} else {
			self.core.merge_output.as_ref()
		}.map(|m| m.output_paths(self, smooth))
	}

	pub fn smooth_output(&self) -> Option<&MergeOutput> { self.core.smooth_output.as_ref() }

	pub fn summary(&self) -> bool { self.core.summary }

	pub fn regions(&self) -> &Regions { &self.regions }

	pub fn regions_mut(&mut self) -> &mut Regions { &mut self.regions }

	pub fn min_counts(&self) -> u32 { self.core.min_counts }

	pub fn min_samples(&self) -> usize { self.core.min_samples.map(usize::from).unwrap_or(1) }

	pub fn min_sdev(&self) -> f64 { self.core.min_sdev }
	
	pub fn threads(&self) -> usize {
		self.core.threads.map(usize::from).unwrap_or_else(num_cpus::get)
	}

	pub fn n_bins(&self) -> usize { self.core.n_bins.map(usize::from).unwrap_or(100) }

	pub fn similarity_type(&self) -> Option<SimilarityType> { self.core.similarity_type }
	
	pub fn prefix(&self) -> &str { self.core.prefix.as_deref().unwrap() }
	
	pub fn dir(&self) -> Option<&Path> { self.core.dir.as_deref() }

	pub fn compress(&self) -> bool { self.core.compress }

	pub fn specific_sites(&self) -> bool { self.core.specific_sites }

//	pub fn model(&self) -> bool { self.core.model }

	pub fn distribution(&self) -> bool { self.core.distribution }

	pub fn heatmaps(&self) -> bool { self.core.heatmaps }

	pub fn blank_lines(&self) -> bool { self.core.blank_lines }

	pub fn smooth(&self) -> bool { self.core.smooth }

	pub fn round_counts(&self) -> bool { self.core.round_counts }

	pub fn window_size(&self) -> usize { self.core.window_size }

	pub fn min_sites(&self) -> usize{ self.core.min_sites }

	pub fn max_distance(&self) -> usize{ self.core.max_distance }

	pub fn kde(&self) -> bool { self.core.kde }

	pub fn kde_cache_limit(&self) -> u16 { self.core.kde_cache_limit }

	pub fn heatmaps_full_limit(&self) -> u8 { self.core.heatmaps_full_limit }

	pub fn combine_cpgs(&self) -> bool { self.core.combine_cpgs }

	pub fn assume_cpg(&self) -> bool { self.core.assume_cpg }

	pub fn n_samples(&self) -> usize { self.sample_info.n_samples() }

	pub fn sample_info(&self) -> &SampleInfo { &self.sample_info }

	pub fn mk_path(&self, name: &str, add_suffix: bool, smooth: bool) -> PathBuf {
		let mut s = if smooth {
			format!("{}_smooth_{}", self.prefix(), name)
		} else {
			format!("{}_{}", self.prefix(), name)
		};
		if add_suffix && self.compress() && !s.ends_with(".gz") {
			s.push_str(".gz")
		}
		let mut fname = PathBuf::from(s);
		if let Some(d) = self.dir() {
			fname = d.join(&fname);
		}
		fname
	}
}

