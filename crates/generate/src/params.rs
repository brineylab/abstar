use std::collections::BTreeMap;
use std::path::PathBuf;
use std::str::FromStr;

use rand::Rng;
use rand_distr::{Beta, Distribution};

use crate::GenerateError;

/// A per-sequence draw: `"0.1"` | `"fixed:0.1"` | `"uniform:0.15-0.35"` | `"beta:2,10"`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Spec {
    Fixed(f64),
    Uniform(f64, f64),
    Beta(f64, f64),
}

impl Spec {
    pub fn sample(&self, rng: &mut impl Rng) -> f64 {
        match *self {
            Spec::Fixed(x) => x,
            Spec::Uniform(a, b) => a + (b - a) * rng.random::<f64>(),
            Spec::Beta(a, b) => Beta::new(a, b).expect("validated at parse").sample(rng),
        }
    }

    /// Smallest and largest value a draw can take.
    pub(crate) fn range(&self) -> (f64, f64) {
        match *self {
            Spec::Fixed(x) => (x, x),
            Spec::Uniform(a, b) => (a, b),
            Spec::Beta(..) => (0.0, 1.0),
        }
    }
}

/// The form `FromStr` reads: `fixed:X`, `uniform:LO-HI`, `beta:A,B`.
impl std::fmt::Display for Spec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            Spec::Fixed(x) => write!(f, "fixed:{x}"),
            Spec::Uniform(a, b) => write!(f, "uniform:{a}-{b}"),
            Spec::Beta(a, b) => write!(f, "beta:{a},{b}"),
        }
    }
}

impl FromStr for Spec {
    type Err = GenerateError;

    fn from_str(s: &str) -> Result<Self, GenerateError> {
        let bad = || {
            GenerateError::InvalidArgument(format!(
                "bad distribution '{s}': expected a number, fixed:X, uniform:LO-HI or beta:A,B"
            ))
        };
        let num = |x: &str| x.trim().parse::<f64>().ok().filter(|v| v.is_finite()).ok_or_else(bad);
        let (kind, rest) = s.split_once(':').unwrap_or(("fixed", s));
        match kind.trim() {
            "fixed" => Ok(Spec::Fixed(num(rest)?)),
            "uniform" => {
                // the `-` between two numbers, so `1e-3-0.1` and `-0.1-0.2` parse
                let (a, b) = rest
                    .match_indices('-')
                    .find_map(|(i, _)| Some((num(&rest[..i]).ok()?, num(&rest[i + 1..]).ok()?)))
                    .ok_or_else(bad)?;
                if a > b {
                    return Err(bad());
                }
                Ok(Spec::Uniform(a, b))
            }
            "beta" => {
                let (a, b) = rest.split_once(',').ok_or_else(bad)?;
                let (a, b) = (num(a)?, num(b)?);
                if a <= 0.0 || b <= 0.0 {
                    return Err(bad());
                }
                Ok(Spec::Beta(a, b))
            }
            _ => Err(bad()),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct GenerateParams {
    pub n: usize,
    pub seed: u64,
    /// Target realized V divergence per sequence.
    pub mutation: Spec,
    /// Per junction, probability of drawing the N length from `insertions_long_tail`.
    pub long_cdr3_weight: f64,
    /// Per-sequence probability of one indel.
    pub indel_rate: f64,
    pub frameshift_indels: bool,
    pub indel_duplication: bool,
    /// nt removed from the 5' / 3' end of the finished read.
    pub truncate_5p: Spec,
    pub truncate_3p: Spec,
    pub constant_len: Option<usize>,
    /// C gene name → weight; empty ⇒ uniform over the locus's C alleles.
    pub isotype_weights: BTreeMap<String, f64>,
    pub productive_only: bool,
    /// Overrides the model's `shm_model` ("motif" | "none").
    pub shm_model: Option<String>,
    /// S5F directory; replaces the motif model wherever the model mutates.
    pub shm_table: Option<PathBuf>,
    pub threads: Option<usize>,
    pub batch_size: usize,
}

pub const PRESETS: [&str; 4] = ["naive", "memory", "bnab", "tcr"];

pub fn preset(name: &str) -> Result<GenerateParams, GenerateError> {
    let base = GenerateParams {
        n: 1000,
        seed: 0,
        mutation: Spec::Fixed(0.0),
        long_cdr3_weight: 0.0,
        indel_rate: 0.0,
        frameshift_indels: false,
        indel_duplication: false,
        truncate_5p: Spec::Fixed(0.0),
        truncate_3p: Spec::Fixed(0.0),
        constant_len: None,
        isotype_weights: BTreeMap::new(),
        productive_only: false,
        shm_model: None,
        shm_table: None,
        threads: None,
        batch_size: 4096,
    };
    Ok(match name {
        "naive" | "tcr" => base,
        "memory" => GenerateParams { mutation: Spec::Uniform(0.02, 0.12), ..base },
        "bnab" => GenerateParams {
            mutation: Spec::Uniform(0.15, 0.35),
            // Tuned on productive output (`productive_only`, the preset's default), over
            // productive human-igh reads (tests/validation.rs `bnab_targets`, 50,000 reads,
            // seed 2026). long_cdr3_weight: smallest round weight with margin for ≥ 25% CDRH3
            // ≥ 24 aa and ≥ 2% ≥ 35 aa; 0.08 gives 28.0% / 2.43% (7 SE above 2%); 0 gives
            // 21.6% / 0.34%, 0.07 gives 27.2% / 2.18% (3 SE). indel_rate: 0.110 centers
            // has_indel at 10.0% of productive reads (0.108 gave 9.87%: indel reads are less
            // often productive). Mean 6.1 attempts per read.
            long_cdr3_weight: 0.08,
            indel_rate: 0.110,
            productive_only: true,
            ..base
        },
        _ => {
            return Err(GenerateError::InvalidArgument(format!(
                "unknown preset '{name}'; valid presets: {}",
                PRESETS.join(", ")
            )))
        }
    })
}
