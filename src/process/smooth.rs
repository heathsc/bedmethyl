use std::{
   collections::VecDeque,
};

use super::reader::{Record, ImpCounts, InputFile};

const TRICUBE_CONST: f64 = 70.0 / 81.0;
// Used to bound the methylation away from 0 or 1
const FIT_CONST: f64 = 0.1;

/// Results of fitting smoothing curve
#[derive(Debug, Copy, Clone, Default)]
pub struct SmoothFit {
   effects: [f64; 3], // [a0, a1, a2]: we fit a local quadratic, so y = a0 + a1 * x + a2 * x^2
   cov: [f64; 6], // Var-cov matrix (lower triangle) for effect estimates
   pub(super) bandwidth: f64,
   pub(super) pos: usize,
}

impl SmoothFit {
   pub(super) fn impute(&self, pos: usize) -> ImpCounts {
      let eff = &self.effects;
      let cov = &self.cov;
      let (z, var) = if pos == self.pos {
         (eff[0], cov[0])
      } else {
         let x = ((pos as f64) - (self.pos as f64))  / self.bandwidth;
         let x2 = x * x;
         (
            eff[0] + x * eff[1] + x2 * eff[2],
            cov[0] + x * cov[1] + x2 * cov[3]
            + x * (cov[1] + x * cov[2] + x2 * cov[4])
            + x2 * (cov[3] + x * cov[4] + x2 * cov[5])
         )
      };
      let m = (z.sin() + 1.0) * 0.5;
      let n = 1.0 / var;
      ImpCounts {
         non_converted: m * n,
         converted: (1.0 - m) * n,
      }
   }
}


#[derive(Debug, Copy, Clone, Eq, PartialEq)]
enum SmState {
   Filling,
   InBlock,
   Draining,
   Finished,
}

pub(super) struct SmoothFile {
   n_sites: usize, // Minimum number of observed sites
   max_distance: usize, // Maximum distance between adjacent sites (0 = no limit, except same chromosome)
   bandwidth: usize, // Minimum bandwidth (bp)
   curr_idx: usize,
   cache: Option<Record>,
   data: VecDeque<Record>,
   work: Vec<FitVal>,
   state: SmState,
   eof: bool,
}

impl SmoothFile {
   pub(super) fn new(n_sites: usize, dist: usize, bandwidth: usize) -> Self {
      trace!("Setting up new SmoothFile");
      let max_distance = if dist > 0 { dist } else { usize::MAX };
      Self { n_sites, max_distance, bandwidth, curr_idx: 0, cache: None,
         data: VecDeque::with_capacity(n_sites),
         work: Vec::with_capacity(n_sites),
         state: SmState::Filling, eof: false }
   }

   fn read_rec(&mut self, f: &mut InputFile) -> anyhow::Result<Option<Record>> {
      match self.cache.take() {
         None => f.read_rec(),
         x => Ok(x),
      }
   }

   pub(super) fn get_rec(&mut self, f: &mut InputFile) -> anyhow::Result<Option<Record>> {

      loop {
         match self.state {

            SmState::Finished => {
               assert!(self.data.is_empty());
               return Ok(None)
            },

            SmState::Draining => {
               if let Some(rec) = self.data.pop_front() {
                  return Ok(Some(rec))
               } else {
                  self.curr_idx = 0;
                  self.state = if self.eof { SmState::Finished } else { SmState::Filling };
               }
            },

            SmState::Filling => {
               if self.try_add_record(f)? && self.range_ok() {
                  self.fit();
                  self.state = SmState::InBlock
               }
            },

            SmState::InBlock => {
               // Trim left end of block if possible
               if self.data.len() > self.n_sites && self.curr_idx > 0 && self.data[self.curr_idx].pos - self.data[1].pos >= self.bandwidth {
                  self.curr_idx -= 1;
                  return Ok(self.data.pop_front())
               }
               // Add to right end if possible
               if !self.check_right_span() {
                  self.try_add_record(f)?;
               } else {
                  // See if the next record (if it exists) is closer to the index position than the left most record
                  // and if the following record on the left still maintains the desired bandwidth, we can remove
                  // the leftmost record and add a new record on the right
                  if self.curr_idx > 0 {
                     if let Some(rec) = self.check_next_rec(f)? {
                        let x = self.data[self.curr_idx].pos;
                        if x - self.data[1].pos >= self.bandwidth && rec.pos - x <= x - self.data[0].pos {
                           self.data.push_back(rec);
                           self.curr_idx -= 1;
                           return Ok(self.data.pop_front())
                        } else {
                           self.cache = Some(rec)
                        }
                        self.fit();
                     } else {
                        self.drain()
                     }
                  } else {
                     self.fit();
                  }
               }
            },
         }
      }
   }

   // Fit model at curr_idx position
   fn fit(&mut self) -> bool {
      assert!(self.curr_idx < self.data.len() && self.data.len() >= self.n_sites, "{} {}", self.data.len(), self.curr_idx);
      let pos = self.data[self.curr_idx].pos;
      let h = (self.data.back().unwrap().pos - pos).max(pos - self.data[0].pos) as f64;
      let zpos = pos as f64;
      self.work.clear();
      let mut xwx = [0.0; 6]; // Lower triangle of X'WX matrix
      let mut xwz = [0.0; 3]; // X'WZ vector
      let mut denom = 0.0;

      // Add data to work, build up X'WX and X'YZ
      for r in self.data.iter() {
         let x = (((r.pos as f64) - zpos) as f64) / h;
         let n = (r.counts.n() as f64) + 2.0 * FIT_CONST;
         let m = ((r.counts.non_converted as f64) + FIT_CONST) / n;
         denom += n - 1.0;
         let val = FitVal::new(n, x, m);
         val.accum(&mut xwx, &mut xwz);
         self.work.push(val)
      }
      xwx[3] = xwx[2]; // These are both equal to Sigma w * x^2

      // Cholesky decomposition
      let c = chol(&mut xwx);
      // Get effects
      let b = chol_solve(&c, &mut xwz);
      // Get residual SS
      let rss: f64 = self.work.iter().map(|v| v.rss(b)).sum();
      // Estimate phi
      let s2 = rss / ((self.work.len() - 3) as f64);
      let phi = ((self.work.len() as f64) * (s2 - 1.0) / denom).max(0.001).min(0.999);

      // Update GLS matrices
      xwx = [0.0; 6];
      xwz = [0.0; 3];
      for v in self.work.iter_mut() {
         v.update_wt(v[0] / (1.0 + (v[0] - 1.0) * phi));
         v.accum(&mut xwx, &mut xwz);
      }
      xwx[3] = xwx[2]; // These are both equal to Sigma w * x^2

      // Cholesky decomposition
      let c = chol(&mut xwx);
      // Get updated effects
      let b = chol_solve(&c, &mut xwz);
      // Get Var-cov matrix of effecys
      let inv = chol_inv(c);

      // And store results
      self.data[self.curr_idx].smooth_fit = Some(SmoothFit {
         effects: *b,
         cov: inv,
         bandwidth: h,
         pos,
      });

      self.curr_idx += 1;
      self.curr_idx < self.data.len()
   }

   // Returns true if we have enough sites and the bandwidth is sufficient for the current position
   fn range_ok(&self) -> bool {
      if self.data.len() < self.n_sites || self.curr_idx >= self.data.len() {
         false
      } else {
         self.check_right_span()
      }
   }

   fn check_right_span(&self) -> bool {
      self.data.back().unwrap().pos - self.data[self.curr_idx].pos >= self.bandwidth
   }

   fn check_next_rec(&mut self, f: &mut InputFile) -> anyhow::Result<Option<Record>> {
      if let Some(rec) = self.read_rec(f)? {
         if self.data.back().map(|r| r.reg_idx == rec.reg_idx && rec.pos - r.pos <= self.max_distance).unwrap_or(true) {
            Ok(Some(rec))
         } else {
            self.cache = Some(rec);
            Ok(None)
         }
      } else {
         self.eof = true;
         Ok(None)
      }
   }

   fn try_add_record(&mut self, f: &mut InputFile) -> anyhow::Result<bool> {

      Ok(if let Some(rec) = self.check_next_rec(f)? {
         self.data.push_back(rec);
         true
      } else {
         self.drain();
         false
      })
   }

   fn drain(&mut self) {
      while self.data.len() >= self.n_sites && self.curr_idx < self.data.len() {
          self.fit();
      }
      self.state = SmState::Draining;
   }
}

// Cholesky decomposition for a 3 x 3 symmetric PD matrix
// supplied with the lower triangle as a 6 element vector
fn chol(x: &mut [f64; 6]) -> &[f64; 6] {
   assert!(x[0] > 0.0, "Matrix not PD");
   x[0] = x[0].sqrt(); // L[0,0]
   x[1] /= x[0]; // L[1,0]
   let z = x[2] - x[1] * x[1];
   assert!(z > 0.0, "Matrix not PD");
   x[2] = z.sqrt(); // L[1,1]
   x[3] /= x[0]; // L[2,0]
   x[4] = (x[4] - x[1] * x[3]) / x[2]; // L[2,1]
   let z = x[5] - x[3] * x[3] - x[4] * x[4];
   assert!(z > 0.0, "Matrix not PD");
   x[5] = z.sqrt(); // L[2;2]
   x
}

fn chol_solve<'a>(c: &[f64; 6], y: &'a mut [f64; 3]) -> &'a [f64; 3] {
   let a0 = y[0] / c[0];
   let a1 = (y[1] - c[1] * a0) / c[2];
   let a2 = (y[2] - c[3] * a0 - c[4] * a1) / c[5];
   y[2] = a2 / c[5];
   y[1] = (a1 - c[4] * y[2]) / c[2];
   y[0] = (a0 - c[3] * y[2] - c[1] * y[1]) / c[0];
   y
}

// Calculate inverse of 3 x 3 symmetric PD matrix given
// the cholesky decomposition (lower triangular)
fn chol_inv(c: &[f64; 6]) -> [f64; 6] {
   // First col
   let a0 = 1.0 / c[0];
   let a1 = -c[1] * a0 / c[2];
   let a2 = (-c[3] * a0 - c[4] * a1) / c[5];
   let i3 = a2 / c[5];
   let i1 = (a1 - c[4] * i3) / c[2];
   let i0 = (a0 - c[3] * i3 - c[1] * i1) / c[0];
   // Second col
   let a1 = 1.0 / c[2];
   let a2 = (-c[4] * a1) / c[5];
   let i4 = a2 / c[5];
   let i2 = (a1 - c[4] * i4) / c[2];
   // last col
   let i5 = 1.0 / (c[5] * c[5]);
   [i0, i1, i2, i3, i4, i5]
}

#[derive(Copy, Clone)]
struct FitVal([f64;7]); // wt, x, x^2, x^3, x^4, z, tc

impl FitVal {
   fn new(wt: f64, x: f64, m: f64) -> Self {
      let x1 = x.abs();
      assert!(x1 <= 1.0, "x1 = {}", x1);
      let z = x * x;
      let z1 = 1.0 - z * x1;

      // Tricube kernel
      let tc = TRICUBE_CONST * z1 * z1 * z1;

      // Transform meth value
      let y = (2.0 * m - 1.0).asin();
      Self([wt, x, z, x * z, z * z, y, tc])
   }

   fn accum(&self, xwx: &mut [f64; 6], xwz: &mut [f64; 3]) {
      let w = self[0] * self[6];
      // Accumulate lower triangle of XWX
      xwx[0] += w; // Sigma w
      xwx[1] += w * self[1]; // w * x
      xwx[2] += w * self[2]; // w * x^2 - note the xwx[6] == xwx[4]
      xwx[4] += w * self[3]; // w * x^3
      xwx[5] += w * self[4]; // w * x^4
      xwz[0] += w * self[5]; // w * z
      xwz[1] += w * self[1] * self[5]; // w * x * z
      xwz[2] += w * self[2] * self[5]; // w * x^2 * z
   }

   fn update_wt(&mut self, wt: f64) { self.0[0] = wt }

   // Calculate z - X * Beta
   fn resid(&self, beta: &[f64; 3]) -> f64 {
      self[5] - (beta[0] + beta[1] * self[1] + beta[2] * self[2])
   }

   // Calculate n * (z - X * Beta)^2
   fn rss(&self, beta: &[f64; 3]) -> f64 {
      let e = self.resid(beta);
      self[0] * e * e
   }
}

impl std::ops::Deref for FitVal {
   type Target = [f64];

   fn deref(&self) -> &Self::Target {
      &self.0
   }
}
