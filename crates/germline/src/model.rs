use crate::GermlineError;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;
use std::str::FromStr;

macro_rules! str_enum {
    ($t:ident { $($v:ident = $s:literal),+ }) => {
        impl fmt::Display for $t {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(match self { $(Self::$v => $s),+ })
            }
        }
        impl FromStr for $t {
            type Err = GermlineError;
            fn from_str(s: &str) -> Result<Self, GermlineError> {
                match s {
                    $($s => Ok(Self::$v),)+
                    _ => Err(GermlineError::Invalid(format!(
                        "unknown {} {:?}", stringify!($t).to_lowercase(), s))),
                }
            }
        }
    };
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Locus { IGH, IGK, IGL, TRA, TRB, TRG, TRD }
str_enum!(Locus { IGH = "IGH", IGK = "IGK", IGL = "IGL", TRA = "TRA", TRB = "TRB", TRG = "TRG", TRD = "TRD" });

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Segment { V, D, J, C }
str_enum!(Segment { V = "V", D = "D", J = "J", C = "C" });

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Functionality { F, ORF, P, NF, Unknown }
str_enum!(Functionality { F = "F", ORF = "ORF", P = "P", NF = "NF", Unknown = "unknown" });

#[derive(Debug, Clone, PartialEq)]
pub struct Allele {
    pub id: String,
    pub name: String,
    pub gene: String,
    pub locus: Locus,
    pub segment: Segment,
    pub species: String,
    pub functionality: Functionality,
    pub sequence_gapped: String,
    /// 0-based offset into the ungapped sequence: first base of the V Cys104 or J Trp/Phe118 codon.
    pub anchor: Option<usize>,
    pub source: String,
    pub usage: Option<f64>,
}

impl Allele {
    /// Gapped sequence with '.' removed.
    pub fn sequence(&self) -> String {
        self.sequence_gapped.chars().filter(|&c| c != '.').collect()
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Provenance {
    pub species: String,
    pub files: Vec<String>,
    pub sha256: Vec<String>,
    pub description: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BuildWarning {
    pub allele: String,
    pub message: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Manifest {
    pub format_version: u32,
    pub name: String,
    pub description: String,
    pub abstar_version: String,
    /// YYYY-MM-DD, supplied by the caller.
    pub built: String,
    pub alleles_sha256: String,
    pub species: Vec<String>,
    pub loci: Vec<Locus>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub prefixes: BTreeMap<String, String>,
    #[serde(default)]
    pub provenance: Vec<Provenance>,
    #[serde(default)]
    pub warnings: Vec<BuildWarning>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Database {
    pub manifest: Manifest,
    pub alleles: Vec<Allele>,
}

impl Database {
    /// Allele counts keyed "LOCUS/SEGMENT/species", e.g. "IGH/V/human".
    pub fn counts(&self) -> BTreeMap<String, usize> {
        let mut counts = BTreeMap::new();
        for a in &self.alleles {
            *counts.entry(format!("{}/{}/{}", a.locus, a.segment, a.species)).or_insert(0) += 1;
        }
        counts
    }
}
