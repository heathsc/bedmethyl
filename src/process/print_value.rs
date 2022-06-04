use std::io::{self, Write, BufWriter};

use r_htslib::*;

use crate::config::*;

type Wrt = BufWriter<OwnedWriter>;

pub(super) trait PrintValue {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()>;
	fn print_missing(&mut self) -> io::Result<()>;
	fn print_header(&mut self, id: &str) -> io::Result<()>;
	fn print_str(&mut self, s: &str) -> io::Result<()>;
}

struct NoOutput {}

impl PrintValue for NoOutput {
	fn print_value(&mut self, _a: u32, _b: u32) -> io::Result<()> { Ok(()) }
	fn print_missing(&mut self) -> io::Result<()> { Ok(()) }
	fn print_header(&mut self, _id: &str) -> io::Result<()> { Ok(()) }
	fn print_str(&mut self, _s: &str) -> io::Result<()> { Ok(()) }
}

struct NonconvConv {
	wrt: Wrt,
	delim: char,
}

impl PrintValue for NonconvConv {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()> { write!(self.wrt, "\t{}{}{}", a, self.delim, b) }
	fn print_missing(&mut self) -> io::Result<()> { write!(self.wrt, "\tNA{}NA", self.delim) }
	fn print_header(&mut self, id: &str) -> io::Result<()> { write!(self.wrt, "\t{}:nc{}{}:c", id, self.delim, id) }
	fn print_str(&mut self, s: &str) -> io::Result<()> { write!(self.wrt, "{}", s) }
}

struct SplitNonconvConv {
	pub(super) wrt: Vec<Wrt>,
}

impl PrintValue for SplitNonconvConv {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()> {
		write!(self.wrt[0], "\t{}", a)?;
		write!(self.wrt[1], "\t{}", b)
	}
	fn print_missing(&mut self) -> io::Result<()> {
		write!(self.wrt[0], "\tNA")?;
		write!(self.wrt[1], "\tNA")
	}
	fn print_header(&mut self, id: &str) -> io::Result<()> {
		write!(self.wrt[0], "\t{}:nc", id)?;
		write!(self.wrt[1], "\t{}:c", id)
	}
	fn print_str(&mut self, s: &str) -> io::Result<()> {
		write!(self.wrt[0], "{}", s)?;
		write!(self.wrt[1], "{}", s)
	}
}

struct NonconvCov {
	wrt: Wrt,
	delim: char,
}

impl PrintValue for NonconvCov {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()> { write!(self.wrt, "\t{}{}{}", a, self.delim, a + b) }
	fn print_missing(&mut self) -> io::Result<()> { write!(self.wrt, "\tNA{}NA", self.delim) }
	fn print_header(&mut self, id: &str) -> io::Result<()> { write!(self.wrt, "\t{}:nc{}{}:cov", id, self.delim, id) }
	fn print_str(&mut self, s: &str) -> io::Result<()> { write!(self.wrt, "{}", s) }
}

struct SplitNonconvCov {
	wrt: Vec<Wrt>,
}

impl PrintValue for SplitNonconvCov {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()> {
		write!(self.wrt[0], "\t{}", a)?;
		write!(self.wrt[1], "\t{}", a + b)
	}
	fn print_missing(&mut self) -> io::Result<()> {
		write!(self.wrt[0], "\tNA")?;
		write!(self.wrt[1], "\tNA")
	}
	fn print_header(&mut self, id: &str) -> io::Result<()> {
		write!(self.wrt[0], "\t{}:nc", id)?;
		write!(self.wrt[1], "\t{}:cov", id)
	}
	fn print_str(&mut self, s: &str) -> io::Result<()> {
		write!(self.wrt[0], "{}", s)?;
		write!(self.wrt[1], "{}", s)
	}
}

struct MethCov {
	wrt: Wrt,
	delim: char,
}

impl PrintValue for MethCov {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()> { write!(self.wrt, "\t{:.4}{}{}", (a as f64) / ((a + b) as f64), self.delim, a + b) }
	fn print_missing(&mut self) -> io::Result<()> { write!(self.wrt, "\tNA{}NA", self.delim) }
	fn print_header(&mut self, id: &str) -> io::Result<()> { write!(self.wrt, "\t{}:meth{}{}:cov", id, self.delim, id) }
	fn print_str(&mut self, s: &str) -> io::Result<()> { write!(self.wrt, "{}", s) }
}

struct SplitMethCov {
	wrt: Vec<Wrt>,
}

impl PrintValue for SplitMethCov {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()> {
		write!(self.wrt[0], "\t{:.4}", (a as f64) / ((a + b) as f64))?;
		write!(self.wrt[1],"\t{}", a + b)
	}
	fn print_missing(&mut self) -> io::Result<()> {
		write!(self.wrt[0], "\tNA")?;
		write!(self.wrt[1], "\tNA")
	}
	fn print_header(&mut self, id: &str) -> io::Result<()> {
		write!(self.wrt[0], "\t{}:meth", id)?;
		write!(self.wrt[1], "\t{}:cov", id)
	}
	fn print_str(&mut self, s: &str) -> io::Result<()> {
		write!(self.wrt[0], "{}", s)?;
		write!(self.wrt[1], "{}", s)
	}
}

struct Meth {
	wrt: Wrt,
}

impl PrintValue for Meth {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()> { write!(self.wrt,"\t{:.4}", (a as f64) / ((a + b) as f64)) }
	fn print_missing(&mut self) -> io::Result<()> { write!(self.wrt, "\tNA") }
	fn print_header(&mut self, id: &str) -> io::Result<()> { write!(self.wrt, "\t{}:meth", id) }
	fn print_str(&mut self, s: &str) -> io::Result<()> { write!(self.wrt, "{}", s) }
}

pub(super) fn get_print_value(cfg: &Config, smooth: bool) -> io::Result<Option<Box<dyn PrintValue>>> {
	let x = if smooth { cfg.smooth_output() } else { cfg.merge_output() };
	Ok(if let Some(mo) = x {
		// Compression mode
		let mode = if cfg.compress() { "wz" } else { "w" };
		let mut wrt = Vec::new();

		for s in cfg.outputs(smooth).expect("No outputs").iter() {
			debug!("Opening output file {:?} ", s);
			let h = Hts::open(s, mode)?;
			let mut w = OwnedWriter::new(h)?;
			if cfg.compress() {
				w.set_mt(num_cpus::get_physical(), 256);
			}
			wrt.push(BufWriter::new(w))
		}
		
		let output_type = mo.output_type();
		Some(if wrt.len() > 1 {
			match output_type {
				OutputType::NonconvConv => Box::new(SplitNonconvConv{wrt}),
				OutputType::NonconvCov => Box::new(SplitNonconvCov{wrt}),
				OutputType::MethCov => Box::new(SplitMethCov{wrt}),
				_ => panic!("Output type can not be split"),	
			}
		} else {
			let w = wrt.into_iter().next().unwrap();
			let delim = mo.value_delim().get_delim_char();
			match output_type {
				OutputType::NonconvConv => Box::new(NonconvConv{wrt: w, delim}),
				OutputType::NonconvCov => Box::new(NonconvCov{wrt: w, delim}),
				OutputType::MethCov => Box::new(MethCov{wrt: w, delim}),
				OutputType::Meth => Box::new(Meth{wrt: w}),
			}
		})
	} else {
		None
	})
}
