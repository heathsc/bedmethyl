use std::{
   collections::VecDeque,
};

use anyhow::Context;

use super::reader::{Record, Counts, InputFile};

/// Results of fitting smoothing curve
#[derive(Debug, Copy, Clone, Default)]
pub struct SmoothFit {
   effects: [f64; 3], // [a0, a1, a2]: we fit a local quadratic, so y = a0 + a1 * x + a2 * x^2
   cov: [f64; 6], // Var-cov matrix (lower triangle) for effect estimates
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
   state: SmState,
   eof: bool,
}

impl SmoothFile {
   pub(super) fn new(n_sites: usize, max_distance: usize, bandwidth: usize) -> Self {
      trace!("Setting up new SmoothFile");
      Self { n_sites, max_distance, bandwidth, curr_idx: 0, cache: None, data: VecDeque::with_capacity(n_sites), state: SmState::Filling, eof: false }
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
               trace!("State: Finished");
               return Ok(None)
            },

            SmState::Draining => {
               trace!("State: Draining");
               if let Some(rec) = self.data.pop_front() {
                  return Ok(Some(rec))
               } else {
                  self.state = if self.eof { SmState::Finished } else { SmState::Filling };
               }
            },

            SmState::Filling => {
               trace!("State: Filling");
               if self.try_add_record(f)? && self.range_ok() {
                  self.impute();
                  self.state = SmState::InBlock
               }
            },

            SmState::InBlock => {
               trace!("State: InBlock");
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
                        self.impute();
                     } else {
                        self.drain()
                     }
                  } else {
                     self.impute();
                  }
               }
            },
         }
      }
   }

   // Impute at current position
   fn impute(&mut self) -> bool {
      // Do imputation magic here
      let r0 = &self.data[0];
      let r1 = self.data.back().unwrap();
      let rix = &self.data[self.curr_idx];

      debug!("Impute: {}:{} {}:{}-{}:{} left: {} {}, right: {} {}, ns: {}",
         rix.reg_idx, rix.pos, r0.reg_idx, r0.pos, r1.reg_idx, r1.pos,
         self.curr_idx, rix.pos - r0.pos,
         self.data.len() - self.curr_idx, r1.pos - rix.pos,
         self.data.len()
      );

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
      trace!("try_add_record");

      Ok(if let Some(rec) = self.check_next_rec(f)? {
         self.data.push_back(rec);
         trace!("try_add_record - got new record.  State = {:?}", self.state);
         true
      } else {
         self.drain();
         false
      })
   }

   fn drain(&mut self) {
      while self.impute() {};
      self.state = SmState::Draining;
   }
}