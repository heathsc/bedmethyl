use std::{
    collections::{HashMap, HashSet},
    io::{BufRead, BufReader},
    ops::Deref,
    path::{Path, PathBuf},
    str::FromStr,
};

use anyhow::Context;
use crossbeam_channel::{unbounded, Receiver, Sender};
use crossbeam_utils::thread;
use lazy_static::lazy_static;
use regex::Regex;

use compress_io::compress::{CompressIo, Reader};
use r_htslib::*;

lazy_static! {
    static ref RE: Regex = Regex::new(r#"^@(\S+)\s*<([^>]*)>"#).unwrap();
    static ref RE1: Regex = Regex::new(r#"^\s*(\S+)\s*=\s*"([^"]*)"\s*$"#).unwrap();
    static ref RE2: Regex = Regex::new(r#"^\s*(\S+)\s*=\s*(\S+)\s*$"#).unwrap();
}

pub type SCov = SVar<f64>;

pub struct SampleInfo {
    samples: Vec<Sample>,
    factors: Vec<SFactor>,
    covariates: Vec<SCov>,
}

impl SampleInfo {
    pub fn samples(&self) -> &[Sample] {
        &self.samples
    }
    pub fn factors(&self) -> &[SFactor] {
        &self.factors
    }
    pub fn covariates(&self) -> &[SCov] {
        &self.covariates
    }

    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }
}

pub struct SFactor {
    var: SVar<u32>,
    levels: Vec<String>,
}

impl SFactor {
    //   pub fn var(&self) -> &SVar<u32> { &self.var }
    pub fn levels(&self) -> &[String] {
        &self.levels
    }
    pub fn name(&self) -> &str {
        self.var.name()
    }
}

pub struct SVar<T> {
    name: String,
    values: Vec<Option<T>>,
}

impl<T> Deref for SVar<T> {
    type Target = [Option<T>];

    fn deref(&self) -> &Self::Target {
        &self.values
    }
}

impl<T> SVar<T> {
    pub fn new(name: String) -> Self {
        Self {
            name,
            values: Vec::new(),
        }
    }

    pub fn push(&mut self, x: Option<T>) {
        self.values.push(x)
    }

    pub fn name(&self) -> &str {
        &self.name
    }
}

#[derive(Copy, Clone)]
pub struct Conversion {
    pub conversion: f64,
    pub over_conversion: f64,
}

impl Default for Conversion {
    fn default() -> Self {
        Self {
            conversion: 1.0,
            over_conversion: 0.0,
        }
    }
}

impl Conversion {
    fn new(x1: f64, x2: f64) -> anyhow::Result<Self> {
        if (0.0..=1.0).contains(&x1) && (0.0..=1.0).contains(&x2) && x2 < x1 {
            Ok(Self {
                conversion: x1,
                over_conversion: x2,
            })
        } else {
            Err(anyhow!("Invalid conversion rates {x1} {x2}"))
        }
    }
}

pub struct Sample {
    pub name: String,
    pub path: PathBuf, // Path to bedmethyl file
    pub conversion: Option<Conversion>,
}

impl Sample {
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn path(&self) -> &Path {
        self.path.as_path()
    }
}

// Get union of contig lists from all input files
pub fn get_contig_list(hts_vec: &[Hts]) -> Vec<&str> {
    let mut ctg_hash = HashSet::new();
    for hts in hts_vec {
        for ctg in hts.seq_names() {
            ctg_hash.insert(ctg);
        }
    }
    ctg_hash.iter().copied().collect()
}

fn open_file_thread(
    recv: Receiver<(usize, PathBuf)>,
    send: Sender<(usize, Hts)>,
) -> anyhow::Result<()> {
    for (ix, input_file) in recv {
        trace!("Checking file {}", input_file.display());
        if !input_file.exists() {
            return Err(anyhow!("File {} is not accessible", input_file.display()));
        }
        let hts = Hts::open(Some(&input_file), "r")
            .with_context(|| format!("Error opening file {}", input_file.display()))?;
        trace!("Sending htsfile {} {}", ix, hts.name());
        send.send((ix, hts)).expect("Error sending Sample");
    }
    Ok(())
}

struct FactorInfo {
    vinfo: VInfo,
    values: SVar<u32>,
    hash: HashMap<String, u32>,
}

impl FactorInfo {
    fn process_obs(&mut self, fd: &[&str]) {
        self.values.push(self.vinfo.process_obs(fd).and_then(|s| {
            self.hash.get(s).copied().or_else(|| {
                let i = self.hash.len() as u32;
                self.hash.insert(s.to_owned(), i);
                Some(i)
            })
        }));
    }
}

struct ConversionInfo {
    vinfo: VInfo,
}

impl ConversionInfo {
    fn process_obs(&self, fd: &[&str]) -> anyhow::Result<Option<Conversion>> {
        let x = if let Some(s) = self.vinfo.process_obs(fd) {
            let v: Vec<_> = s.split(',').map(|t| t.trim()).collect();
            if v.len() != 2 {
                return Err(anyhow!("Error getting conversion values from {s} - expected two comma separated values"));
            }
            match (f64::from_str(v[0]), f64::from_str(v[1])) {
                (Ok(x1), Ok(x2)) => Some(Conversion::new(x1, x2)?),
                _ => return Err(anyhow!("Error parsing conversion values from {s}")),
            }
        } else {
            None
        };
        Ok(x)
    }
}

struct CovInfo {
    vinfo: VInfo,
    values: SVar<f64>,
}

impl CovInfo {
    fn process_obs(&mut self, fd: &[&str]) -> anyhow::Result<()> {
        let x = if let Some(s) = self.vinfo.process_obs(fd) {
            Some(f64::from_str(s).with_context(|| {
                format!(
                    "Error getting float for variate {}: '{}'",
                    self.values.name(),
                    s
                )
            })?)
        } else {
            None
        };
        self.values.push(x);
        Ok(())
    }
}

#[derive(Debug)]
struct VInfo {
    col: usize,
    missing: Option<String>,
}

impl VInfo {
    fn process_obs<'a>(&self, fd: &[&'a str]) -> Option<&'a str> {
        let s = fd.get(self.col - 1).and_then(|s| {
            let t = s.trim();
            if t.is_empty() {
                None
            } else {
                Some(t)
            }
        });
        match (s, self.missing.as_ref()) {
            (Some(t), Some(m)) => {
                if t == m {
                    None
                } else {
                    Some(t)
                }
            }
            (Some(t), None) => Some(t),
            _ => None,
        }
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
enum VType {
    Factor,
    Covariate,
    Conversion,
}

/// Process header line that has been split into the type and argument fields
/// i.e., for a header line:
///    @Factor<name="group",col=5>
/// s will be "@Factor" and t will be "name=\"group\",col=5"
fn proc_header_line(
    s: &str,
    t: &str,
) -> anyhow::Result<Option<(VType, String, Option<String>, Option<usize>)>> {
    let vt = match s.to_ascii_lowercase().as_str() {
        "factor" => VType::Factor,
        "covariate" => VType::Covariate,
        "conversion" => VType::Conversion,
        _ => return Ok(None),
    };
    let mut name: Option<String> = None;
    let mut missing: Option<String> = None;
    let mut col: Option<usize> = None;
    for field in t.split(',') {
        if let Some(cap) = RE1.captures(field).or_else(|| RE2.captures(field)) {
            match cap[1].to_ascii_lowercase().as_str() {
                "name" => name = Some(cap[2].to_owned()),
                "missing" => missing = Some(cap[2].to_owned()),
                "col" => {
                    col = {
                        let x = usize::from_str(&cap[2])
                            .with_context(|| "Illegal numeric argument to col")?;
                        Some(x)
                    }
                }
                _ => (),
            }
        }
    }
    if matches!(vt, VType::Conversion) {
        Ok(Some((vt, String::from("conversion"), missing, col)))
    } else if let Some(n) = name {
        Ok(Some((vt, n, missing, col)))
    } else {
        Err(anyhow!("Missing name field"))
    }
}

struct SHeader {
    factors: Vec<FactorInfo>,
    covs: Vec<CovInfo>,
    conversion: Option<ConversionInfo>,
}

fn read_header(
    rd: &mut BufReader<Reader>,
    buf: &mut String,
    line: &mut usize,
) -> anyhow::Result<(usize, SHeader)> {
    let mut fact = Vec::new();
    let mut cov = Vec::new();
    let mut conversion = None;

    let mut current_col = 3;
    let mut name_hash: HashSet<String> = HashSet::new();
    let mut col_hash: HashMap<usize, String> = HashMap::new();
    let mut l: usize;
    loop {
        l = rd
            .read_line(buf)
            .unwrap_or_else(|_| panic!("Error reading from sample file"));
        if l == 0 {
            break;
        }
        *line += 1;
        if buf.starts_with('@') {
            if let Some(cap) = RE.captures(buf) {
                if let Some((vt, name, missing, col)) = proc_header_line(&cap[1], &cap[2])
                    .with_context(|| format!("Error in header line {}: {}", line, buf))?
                {
                    if !name_hash.insert(name.clone()) {
                        return Err(anyhow!(
                            "Variable {} listed multiple times at line {}",
                            name,
                            line
                        ));
                    }
                    let col = col.unwrap_or(current_col);
                    if let Some(s) = col_hash.insert(col, name.clone()) {
                        return Err(anyhow!(
                            "Line {}: Column {} for variable {} already in use with variable {}",
                            line,
                            col,
                            name,
                            s
                        ));
                    }
                    let vinfo = VInfo { missing, col };
                    match vt {
                        VType::Factor => fact.push(FactorInfo {
                            vinfo,
                            values: SVar::new(name),
                            hash: Default::default(),
                        }),
                        VType::Covariate => cov.push(CovInfo {
                            vinfo,
                            values: SVar::new(name),
                        }),
                        VType::Conversion => {
                            if conversion.is_none() {
                                conversion = Some(ConversionInfo { vinfo })
                            } else {
                                return Err(anyhow!(
                                    "Multiple conversion columns listed in header"
                                ));
                            }
                        }
                    }
                    current_col += 1;
                }
            } else {
                return Err(anyhow!("Unrecognized header line at {}: {}", line, buf));
            }
        } else {
            break;
        }
        buf.clear();
    }
    Ok((
        l,
        SHeader {
            factors: fact,
            covs: cov,
            conversion,
        },
    ))
}

/// Read in sample descriptions from sample file.
///
/// The file should be a tab separated text file with at least two columns and optional header lines.
///
/// Column 1 should have the sample name and column 2 the path to the bedmethyl file for that sample.
/// Subsequent columns contain optional sample information (i.e., groups, covariates etc.)
/// The header lines if present give information about the extra columns
pub fn read_sample_file(name: Option<&str>) -> anyhow::Result<SampleInfo> {
    let mut rd = CompressIo::new()
        .opt_path(name)
        .bufreader()
        .with_context(|| "Error opening sample file for input")?;

    debug!("Reading sample information");

    // Read header if present
    let mut line = 0;
    let mut buf = String::new();
    let (l, mut sheader) = read_header(&mut rd, &mut buf, &mut line)
        .with_context(|| "Error reading header from sample file")?;

    if l == 0 {
        return Err(anyhow!("Sample file has no content"));
    }

    // Get largest column required for sample file
    let max_col = sheader
        .factors
        .iter()
        .map(|f| f.vinfo.col)
        .chain(sheader.covs.iter().map(|v| v.vinfo.col))
        .chain(sheader.conversion.iter().map(|v| v.vinfo.col))
        .max()
        .unwrap_or(2);

    let mut samples = Vec::new();

    loop {
        // buf already contains the next line at this point so we can process it directly
        let fd: Vec<&str> = buf.trim().split('\t').collect();
        if fd.len() < max_col {
            warn!("Short line in sample file (line {}) - skipping", line)
        } else {
            let name = fd[0].to_owned();
            let path = PathBuf::from_str(fd[1])
                .unwrap_or_else(|_| panic!("Illegal characters in filename {}", fd[1]));
            for f in sheader.factors.iter_mut() {
                f.process_obs(&fd)
            }
            for c in sheader.covs.iter_mut() {
                c.process_obs(&fd)
                    .with_context(|| format!("Error reading sample file at line {line}"))?
            }
            let conversion = match sheader.conversion.as_mut() {
                Some(c) => c.process_obs(&fd).with_context(|| {
                    format!("Error reading conversion values from sample file at line {line}")
                })?,
                None => None,
            };
            samples.push(Sample {
                name,
                path,
                conversion,
            });
        }

        // Read in next line
        buf.clear();
        let l = rd
            .read_line(&mut buf)
            .unwrap_or_else(|_| panic!("Error reading from sample file"));
        if l == 0 {
            break;
        }
        line += 1;
    }

    // Make SFactor vector
    let factors = sheader
        .factors
        .drain(..)
        .map(|f| {
            let FactorInfo {
                values, mut hash, ..
            } = f;
            let mut v: Vec<_> = hash.drain().collect();
            v.sort_unstable_by_key(|(_, x)| *x);
            let levels: Vec<_> = v.drain(..).map(|(s, _)| s).collect();
            SFactor {
                var: values,
                levels,
            }
        })
        .collect();

    // Make SCov vector
    let covariates = sheader
        .covs
        .drain(..)
        .map(|c| {
            let CovInfo { values, .. } = c;
            values
        })
        .collect();

    Ok(SampleInfo {
        samples,
        factors,
        covariates,
    })
}

/// Open all input files along with associated index files
pub fn open_input_files(samples: &[Sample], nthr: usize) -> anyhow::Result<Vec<Hts>> {
    debug!("Opening all input files + index files");
    let (res_s, res_r) = unbounded();
    let mut open_thread_results = Vec::new();
    thread::scope(|s| {
        let (msg_s, msg_r) = unbounded();
        let jh: Vec<_> = (0..nthr)
            .map(|_| {
                let msg = msg_r.clone();
                let res = res_s.clone();
                s.spawn(move |_| open_file_thread(msg, res))
            })
            .collect();

        for (ix, sample) in samples.iter().enumerate() {
            trace!("Sending {} to open index thread", sample.path().display());
            let pb = sample.path().to_owned();
            msg_s
                .send((ix, pb))
                .expect("Error sending message to open index threads");
        }
        drop(msg_s);
        debug!("Waiting for open file threads to finish");
        for j in jh {
            open_thread_results.push(j.join().expect("Error from thread"))
        }
        drop(res_s);
    })
    .expect("Error creating thread scope");

    let mut failed = false;
    for res in open_thread_results.drain(..) {
        if let Err(e) = res {
            error!("Error when opening input file: {}", e);
            failed = true;
        }
    }
    if failed {
        Err(anyhow!(
            "One or more of the input files could not be opened"
        ))
    } else {
        let mut v: Vec<Option<Hts>> = (0..samples.len()).map(|_| None).collect();
        for (ix, h) in res_r.iter() {
            trace!("Received hts file {} {}", ix, h.name());
            v[ix] = Some(h)
        }
        let hts_vec: Vec<_> = v.drain(..).map(|x| x.unwrap()).collect();

        debug!("All {} sample files opened successfully", samples.len());

        Ok(hts_vec)
    }
}
