use crate::{Allele, BuildWarning, GermlineError, Segment};
use std::collections::HashSet;

pub fn validate(alleles: &[Allele]) -> Result<Vec<BuildWarning>, GermlineError> {
    let bad = |a: &Allele, m: String| GermlineError::Invalid(format!("allele {}: {m}", a.id));
    let mut seen = HashSet::new();
    let mut warnings = Vec::new();
    for a in alleles {
        if !seen.insert(a.id.as_str()) {
            return Err(bad(a, "duplicate id".into()));
        }
        if let Some(c) = a
            .sequence_gapped
            .chars()
            .find(|c| !matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'N' | '.'))
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
