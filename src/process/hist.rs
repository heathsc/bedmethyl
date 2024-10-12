use std::{
    collections::HashMap,
    io::{self, Write},
};

use crate::config::Config;
use crate::sample::Sample;

// Stores ln(x) and ln(1x) for
// the midpoints of all bins
pub(super) struct LnP(Vec<(f64, f64)>);

impl LnP {
    pub(super) fn new(n: usize) -> Self {
        let z = 1.0 / (n as f64);
        let lnp: Vec<_> = (0..n)
            .map(|i| {
                let x = ((i as f64) + 0.5) * z;
                (x.ln(), (1.0 - x).ln())
            })
            .collect();

        Self(lnp)
    }

    pub(super) fn get_dist(&self, a: u32, b: u32, wk: &mut Vec<f64>) -> f64 {
        let a = a as f64;
        let b = b as f64;
        let mut max = -f64::INFINITY;
        wk.clear();
        for (lp, lp1) in self.0.iter() {
            let z = a * lp + b * lp1;
            wk.push(z);
            max = max.max(z);
        }
        let mut sum = 0.0;
        for z in wk.iter_mut() {
            *z = (*z - max).exp();
            sum += *z;
        }
        sum
    }

    pub(super) fn n_bins(&self) -> usize {
        self.0.len()
    }
}

struct Beta {
    lnp: LnP,
    tmp: Vec<f64>,
    cache: HashMap<(u32, u32), Box<[u32]>>,
}

impl Beta {
    fn new(n: usize) -> Self {
        let lnp = LnP::new(n);
        let tmp = Vec::with_capacity(n);

        Self {
            lnp,
            tmp,
            cache: HashMap::new(),
        }
    }

    fn get_dist(&mut self, a: u32, b: u32) -> (f64, &[f64]) {
        let sum = self.lnp.get_dist(a, b, &mut self.tmp);
        (sum, &self.tmp)
    }

    fn add_to_cache(&mut self, ns: usize, ix: usize, a: u32, b: u32) {
        let ct = self
            .cache
            .entry((a, b))
            .or_insert_with(|| vec![0; ns].into_boxed_slice());
        ct[ix] += 1;
    }

    fn handle_cache(&mut self, bins: &mut Vec<Vec<f64>>) {
        for ((a, b), v) in self.cache.drain() {
            let sum = self.lnp.get_dist(a, b, &mut self.tmp);
            for (ix, ct) in v.iter().enumerate() {
                if *ct > 0 {
                    let z = (*ct as f64) / sum;
                    for (src, bin) in self.tmp.iter().zip(bins[ix].iter_mut()) {
                        *bin += src * z
                    }
                }
            }
        }
    }
}

pub(super) struct Hist {
    bins: Vec<Vec<f64>>,
    beta: Option<Beta>,
    kde_cache_limit: u32,
    n_obs: Vec<usize>,
    n_bins: usize,
    n_samples: usize,
}

impl Hist {
    pub(super) fn new(n_samples: usize, cfg: &Config) -> Self {
        let beta = if cfg.kde() {
            Some(Beta::new(cfg.n_bins()))
        } else {
            None
        };

        let bins: Vec<_> = (0..n_samples).map(|_| vec![0.0; cfg.n_bins()]).collect();
        debug!(
            "Setting up Hist for {} samples and {} bins",
            n_samples,
            cfg.n_bins()
        );
        Self {
            bins,
            beta,
            kde_cache_limit: cfg.kde_cache_limit() as u32,
            n_obs: vec![0; n_samples],
            n_bins: cfg.n_bins(),
            n_samples,
        }
    }

    pub(super) fn add_obs<F: Fn(u32, u32) -> f64>(&mut self, ix: usize, a: u32, b: u32, f: F) {
        if let Some(beta) = self.beta.as_mut() {
            if a + b < self.kde_cache_limit {
                beta.add_to_cache(self.n_samples, ix, a, b)
            } else {
                let (sum, v) = beta.get_dist(a, b);
                for (src, bin) in v.iter().zip(self.bins[ix].iter_mut()) {
                    *bin += src / sum
                }
            }
        } else {
            let m = f(a, b);
            let i = ((m * (self.n_bins as f64)).trunc() as usize).min(self.n_bins - 1);
            self.bins[ix][i] += 1.0;
        }
        self.n_obs[ix] += 1;
    }

    pub(super) fn print_table<W: Write>(&self, wrt: &mut W, samples: &[Sample]) -> io::Result<()> {
        
        let mut tot = vec!(0.0; samples.len());
        let nb = self.n_bins as f64;
        
        for bin in 0..self.n_bins {
            let x = ((bin as f64) + 0.5) / nb;
            for (((v, sm), n_obs),t) in self.bins.iter().zip(samples.iter()).zip(self.n_obs.iter()).zip(tot.iter_mut()) {
                let pos = match sm.conversion.as_ref() {
                    Some(c) => (x - 1.0 + c.conversion) / (c.conversion - c.over_conversion),
                    None => x,
                };
                if (0.0..=1.0).contains(&pos) {
                    let z = if *n_obs > 0 {
                        v[bin] / (*n_obs as f64)
                    } else {
                        0.0
                    };
                    *t += z
                } 
            }
        }
        
        write!(wrt, "index")?;
        for (ix, sm) in samples.iter().enumerate() {
            let id = sm.name();
            write!(wrt, "\t{id}_pos\t{id}_freq")?;
            debug!("Sample {id}, n_obs: {}", self.n_obs[ix]);
        }
        writeln!(wrt)?;

        
        for bin in 0..self.n_bins {
            let x = ((bin as f64) + 0.5) / nb;
            write!(wrt, "{}", bin)?;
            for (((v, sm), n_obs), t) in self.bins.iter().zip(samples.iter()).zip(self.n_obs.iter()).zip(tot.iter()) {
                let pos = match sm.conversion.as_ref() {
                    Some(c) => (x - 1.0 + c.conversion) / (c.conversion - c.over_conversion),
                    None => x,
                };
                if (0.0..=1.0).contains(&pos) {
                    let z = if *n_obs > 0 {
                        v[bin] * nb / (*t * (*n_obs as f64))
                    } else {
                        0.0
                    };
                    write!(wrt, "\t{pos}\t{z}")?;
                } else {
                    write!(wrt, "\t{}\t-", pos.clamp(0.0,1.0))?;
                }
            }
            writeln!(wrt)?;
        }
        Ok(())
    }

    pub(super) fn handle_cache(&mut self) {
        if let Some(beta) = self.beta.as_mut() {
            debug!("Adding contributions from cached entries");
            beta.handle_cache(&mut self.bins);
            debug!("Finished adding contributions from cached entries");
        }
    }
}
