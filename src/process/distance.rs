use std::io::{self, Write};

use crate::cli::similarity_type::SimilarityType;
use crate::sample::Sample;

type PItem = (f64, [f64; 5]);

pub(super) struct Distance {
    n_samples: usize,
    accum: Vec<PItem>,
    obs: Vec<(usize, f64, f64)>,
    stype: SimilarityType,
}

impl Distance {
    pub(super) fn new(n_samples: usize, stype: SimilarityType) -> Self {
        let sz = n_samples * (n_samples + 1) / 2;
        Self {
            n_samples,
            accum: vec![(0.0, [0.0; 5]); sz],
            obs: Vec::new(),
            stype,
        }
    }

    fn update<F>(&mut self, mut f: F)
    where
        F: FnMut(&mut PItem, &f64, &f64, f64),
    {
        let (ix, x, ni) = self.obs.last().unwrap();
        let k = ix * (ix + 1) / 2;
        for (iy, y, mi) in self.obs.iter() {
            // wt is the harmonic mean of n and m (note ni = 1/n; mi = 1/m)
            let wt = 2.0 / (ni + mi);
            f(&mut self.accum[k + iy], x, y, wt)
        }
    }

    pub(super) fn add_obs(&mut self, ix: usize, x: f64, n: u32) {
        assert!(ix < self.n_samples);
        assert!(n > 0);

        if let Some(prev) = self.obs.last() {
            assert!(prev.0 < ix);
        }
        self.obs.push((ix, x, 1.0 / (n.min(100) as f64)));
        match self.stype {
            SimilarityType::Distance => self.update(|item, x, y, wt| {
                item.0 += wt;
                item.1[0] += wt * (x - y) * (x - y);
            }),
            SimilarityType::Covariance => self.update(|item, x, y, wt| {
                item.0 += wt;
                item.1[0] += wt * x;
                item.1[1] += wt * y;
                item.1[2] += wt * x * y;
            }),
            _ => self.update(|item, x, y, wt| {
                item.0 += wt;
                item.1[0] += wt * x;
                item.1[1] += wt * y;
                item.1[2] += wt * x * y;
                item.1[3] += wt * x * x;
                item.1[4] += wt * y * y;
            }),
        }
    }

    pub(super) fn clear_obs(&mut self) {
        self.obs.clear()
    }

    pub(super) fn print_table<W: Write>(&self, wrt: &mut W, samples: &[Sample]) -> io::Result<()> {
        write!(wrt, "Barcode")?;
        for id in samples.iter() {
            write!(wrt, "\t{}", id.name())?;
        }
        writeln!(wrt)?;

        let prt: Box<dyn Fn(&PItem) -> f64> = match self.stype {
            SimilarityType::Distance => Box::new(|item| {
                if item.0 > 0.0 {
                    item.1[0].sqrt() / item.0
                } else {
                    0.0
                }
            }),
            SimilarityType::Covariance => Box::new(|item| {
                if item.0 > 1.0 {
                    (item.1[2] - item.1[0] * item.1[1] / item.0) / (item.0 - 1.0)
                } else {
                    0.0
                }
            }),
            SimilarityType::Correlation => Box::new(|item| {
                if item.0 > 1.0 {
                    let sxy = item.1[2] - item.1[0] * item.1[1] / item.0;
                    if sxy.abs() > 1.0e-12 {
                        let sx = (item.1[3] - item.1[0] * item.1[0] / item.0)
                            .max(0.0)
                            .sqrt();
                        let sy = (item.1[4] - item.1[1] * item.1[1] / item.0)
                            .max(0.0)
                            .sqrt();
                        if sx * sy > 0.0 {
                            sxy / (sx * sy)
                        } else {
                            0.0
                        }
                    } else {
                        0.0
                    }
                } else {
                    0.0
                }
            }),
        };

        for (ix, id) in samples.iter().enumerate() {
            write!(wrt, "{}", id.name())?;
            let k = ix * (ix + 1) / 2;
            for iy in 0..=ix {
                let z: f64 = prt(&self.accum[k + iy]);
                write!(wrt, "\t{}", z)?;
            }
            for iy in ix + 1..self.n_samples {
                let z: f64 = prt(&self.accum[iy * (iy + 1) / 2 + ix]);
                write!(wrt, "\t{}", z)?;
            }
            writeln!(wrt)?;
        }
        Ok(())
    }
}
