use std::io::{self, Write, BufWriter};

use compress_io::{
	compress_type::{CompressThreads, CompressType},
};

use crate::config::*;

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

struct NonconvConv<W: Write> {
	wrt: BufWriter<W>,
	delim: char,
}

impl <W: Write>PrintValue for NonconvConv<W> {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()> { write!(self.wrt, "\t{}{}{}", a, self.delim, b) }
	fn print_missing(&mut self) -> io::Result<()> { write!(self.wrt, "\tNA{}NA", self.delim) }
	fn print_header(&mut self, id: &str) -> io::Result<()> { write!(self.wrt, "\t{}:nc{}{}:c", id, self.delim, id) }
	fn print_str(&mut self, s: &str) -> io::Result<()> { write!(self.wrt, "{}", s) }
}

struct SplitNonconvConv<W: Write> {
	pub(super) wrt: Vec<BufWriter<W>>,
}

impl <W: Write>PrintValue for SplitNonconvConv<W> {
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

struct NonconvCov<W: Write> {
	wrt: BufWriter<W>,
	delim: char,
}

impl <W: Write>PrintValue for NonconvCov<W> {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()> { write!(self.wrt, "\t{}{}{}", a, self.delim, a + b) }
	fn print_missing(&mut self) -> io::Result<()> { write!(self.wrt, "\tNA{}NA", self.delim) }
	fn print_header(&mut self, id: &str) -> io::Result<()> { write!(self.wrt, "\t{}:nc{}{}:cov", id, self.delim, id) }
	fn print_str(&mut self, s: &str) -> io::Result<()> { write!(self.wrt, "{}", s) }
}

struct SplitNonconvCov<W: Write> {
	wrt: Vec<BufWriter<W>>,
}

impl <W: Write>PrintValue for SplitNonconvCov<W> {
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

struct MethCov<W: Write> {
	wrt: BufWriter<W>,
	delim: char,
}

impl <W: Write>PrintValue for MethCov<W> {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()> { write!(self.wrt, "\t{:.4}{}{}", (a as f64) / ((a + b) as f64), self.delim, a + b) }
	fn print_missing(&mut self) -> io::Result<()> { write!(self.wrt, "\tNA{}NA", self.delim) }
	fn print_header(&mut self, id: &str) -> io::Result<()> { write!(self.wrt, "\t{}:meth{}{}:cov", id, self.delim, id) }
	fn print_str(&mut self, s: &str) -> io::Result<()> { write!(self.wrt, "{}", s) }
}

struct SplitMethCov<W: Write> {
	wrt: Vec<BufWriter<W>>,
}

impl <W: Write>PrintValue for SplitMethCov<W> {
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

struct Meth<W: Write> {
	wrt: BufWriter<W>,
}

impl <W: Write>PrintValue for Meth<W> {
	fn print_value(&mut self, a: u32, b: u32) -> io::Result<()> { write!(self.wrt,"\t{:.4}", (a as f64) / ((a + b) as f64)) }
	fn print_missing(&mut self) -> io::Result<()> { write!(self.wrt, "\tNA") }
	fn print_header(&mut self, id: &str) -> io::Result<()> { write!(self.wrt, "\t{}:meth", id) }
	fn print_str(&mut self, s: &str) -> io::Result<()> { write!(self.wrt, "{}", s) }
}

pub(super) fn get_print_value(cfg: &Config) -> io::Result<Option<Box<dyn PrintValue>>> {
	Ok(if let Some(mo) = cfg.merge_output() {
		let ctype = if cfg.compress() {
			Some(CompressType::Bgzip)
		} else {
			None
		};
		let mut wrt = Vec::new();
		for s in cfg.merge_outputs().expect("No outputs").iter() {
			debug!("Opening output file {:?} with compress type {:?}", s, ctype);
			let mut c = compress_io::compress::CompressIo::new();
			if let Some(ct) = ctype {
				c.ctype(ct).cthreads(CompressThreads::NPhysCores);
			}
			wrt.push(c.path(s).bufwriter()?)
		}
		let delim = match mo.value_delim() {
			ValueDelim::Space => ' ',
			ValueDelim::Comma => ',',
			ValueDelim::Semicolon => ';',
			_ => '\t',
		};
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
