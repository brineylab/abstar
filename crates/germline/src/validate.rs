use crate::{anchors::find_v_anchor, Allele, BuildWarning, GermlineError, Segment};
use std::collections::HashSet;

pub fn validate(alleles: &[Allele]) -> Result<Vec<BuildWarning>, GermlineError> {
    let bad = |a: &Allele, m: String| GermlineError::Invalid(format!("allele {}: {m}", a.id));
    let mut seen = HashSet::new();
    let mut warnings = Vec::new();
    for a in alleles {
        if !seen.insert(a.id.as_str()) {
            return Err(bad(a, "duplicate id".into()));
        }
        for (field, v) in [
            ("id", &a.id), ("name", &a.name), ("gene", &a.gene), ("species", &a.species),
            ("sequence_gapped", &a.sequence_gapped), ("source", &a.source),
        ] {
            if v.contains(['\t', '\n', '\r']) {
                return Err(bad(a, format!("{field} contains a tab or newline")));
            }
        }
        if let Some(c) = a
            .sequence_gapped
            .chars()
            // ACGTN + IUPAC ambiguity codes (same alphabet as queries) + IMGT gap '.'
            .find(|c| !matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'N' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' | '.'))
        {
            return Err(bad(a, format!("invalid sequence character {c:?}")));
        }
        let len = a.sequence().len();
        if len == 0 {
            return Err(bad(a, "empty sequence".into()));
        }
        match a.anchor {
            Some(p) if p >= len => {
                return Err(bad(a, format!("anchor {p} out of range for length {len}")))
            }
            None if matches!(a.segment, Segment::V | Segment::J) => warnings.push(BuildWarning {
                allele: a.id.clone(),
                message: format!("{} allele has no anchor", a.segment),
            }),
            _ => {}
        }
    }
    Ok(warnings)
}

/// Share of V alleles that must carry Cys104 at IMGT 310-312 for a database to count as IMGT-gapped.
const MIN_ANCHORED_V: f64 = 0.95;

/// Refuses a database whose V alleles are not IMGT-gapped (record 0036). The rule is
/// database-level: a few real alleles lack Cys104 (mutated or partial), most must not.
pub(crate) fn check_imgt_gapped(name: &str, alleles: &[Allele]) -> Result<(), GermlineError> {
    let v: Vec<&Allele> = alleles.iter().filter(|a| a.segment == Segment::V).collect();
    let bad: Vec<&Allele> = v
        .iter()
        .copied()
        .filter(|a| !a.sequence_gapped.contains('.') || find_v_anchor(&a.sequence_gapped).is_none())
        .collect();
    let ungapped = v.iter().any(|a| !a.sequence_gapped.contains('.'));
    if v.is_empty() || (!ungapped && (v.len() - bad.len()) as f64 >= MIN_ANCHORED_V * v.len() as f64) {
        return Ok(());
    }
    let ex: Vec<&str> = bad.iter().take(5).map(|a| a.name.as_str()).collect();
    Err(GermlineError::Invalid(format!(
        "database {name} is not IMGT-gapped: {} of {} V alleles lack IMGT gaps or Cys104 (TGT/TGC) at \
         IMGT 310-312 (at least {:.0}% must have it, and every V allele needs gaps; e.g. {}); \
         IMGT-gapped V sequences are required",
        bad.len(),
        v.len(),
        MIN_ANCHORED_V * 100.0,
        ex.join(", ")
    )))
}
