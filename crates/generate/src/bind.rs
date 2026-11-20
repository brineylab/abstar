use std::collections::{BTreeMap, BTreeSet};

use abstar_germline::{Allele, Database, Segment};

use crate::{GenerateError, Model};

/// Minimum retained model probability mass for each of V, D and J.
const MIN_COVERAGE: f64 = 0.90;

#[derive(Debug, Clone)]
pub struct Candidates {
    pub v: Vec<Allele>,
    pub d: Vec<Allele>,
    pub j: Vec<Allele>,
    pub c: Vec<Allele>,
}

#[derive(Debug, Clone)]
pub struct Bound {
    pub model: Model,
    pub v: Vec<(Allele, f64)>,
    pub d_given_j: BTreeMap<String, Vec<(Allele, f64)>>,
    pub j: Vec<(Allele, f64)>,
    pub candidates: Candidates,
    pub warnings: Vec<String>,
}

struct Ctx<'a> {
    model: &'a Model,
    alleles: Vec<&'a Allele>,
}

impl Ctx<'_> {
    /// Generatable alleles of `gene` in `seg`, each with its share of the gene.
    fn gene_alleles(&self, gene: &str, seg: Segment) -> Vec<(Allele, f64)> {
        let mut out: Vec<(Allele, f64)> = self
            .alleles
            .iter()
            .filter(|a| a.gene == gene && a.segment == seg && (seg == Segment::D || a.anchor.is_some()))
            .map(|a| ((*a).clone(), self.model.allele_usage.get(&a.name).copied().unwrap_or(0.0)))
            .collect();
        let has_usage = out.iter().any(|(a, _)| self.model.allele_usage.contains_key(&a.name));
        let total: f64 = out.iter().map(|x| x.1).sum();
        if !has_usage || total <= 0.0 {
            out.iter_mut().for_each(|x| x.1 = 1.0);
        }
        let total: f64 = out.iter().map(|x| x.1).sum();
        out.iter_mut().for_each(|x| x.1 /= total);
        out
    }

    /// Bind a gene -> mass table; returns (alleles with renormalized mass, retained mass, dropped genes).
    fn table(&self, usage: &BTreeMap<String, f64>, seg: Segment) -> (Vec<(Allele, f64)>, f64, Vec<String>) {
        let (mut out, mut kept, mut dropped) = (Vec::new(), 0.0, Vec::new());
        for (gene, &p) in usage {
            let alleles = self.gene_alleles(gene, seg);
            if alleles.is_empty() {
                dropped.push(gene.clone());
                continue;
            }
            kept += p;
            out.extend(alleles.into_iter().map(|(a, s)| (a, p * s)));
        }
        if kept > 0.0 {
            out.iter_mut().for_each(|x| x.1 /= kept);
        }
        (out, kept, dropped)
    }
}

fn check_coverage(seg: &str, kept: f64, dropped: &BTreeSet<String>) -> Result<(), GenerateError> {
    if kept < MIN_COVERAGE {
        return Err(GenerateError::ModelInvalid(format!(
            "only {:.1}% of {seg} probability mass maps to alleles in the germline database \
             (need {:.0}%); dropped genes: {}",
            kept * 100.0,
            MIN_COVERAGE * 100.0,
            dropped.iter().cloned().collect::<Vec<_>>().join(", ")
        )));
    }
    Ok(())
}

fn warn(warnings: &mut Vec<String>, seg: &str, dropped: &BTreeSet<String>) {
    for g in dropped {
        warnings.push(format!("{seg} gene {g} has no generatable allele in the germline database; dropped"));
    }
}

pub fn bind(model: Model, db: &Database) -> Result<Bound, GenerateError> {
    let alleles: Vec<&Allele> =
        db.alleles.iter().filter(|a| a.species == model.species && a.locus == model.locus).collect();
    let of = |s: Segment| alleles.iter().filter(|a| a.segment == s).map(|a| (*a).clone()).collect::<Vec<_>>();
    let candidates = Candidates { v: of(Segment::V), d: of(Segment::D), j: of(Segment::J), c: of(Segment::C) };
    if let Some(a) = candidates.v.iter().find(|a| !a.sequence_gapped.contains('.')) {
        return Err(GenerateError::ModelInvalid(format!(
            "V allele {} has no IMGT gaps in its gapped sequence; the germline database is not usable for generation",
            a.id
        )));
    }

    let ctx = Ctx { model: &model, alleles };
    let mut warnings = Vec::new();

    let (v, kept, dropped) = ctx.table(&model.v_usage, Segment::V);
    let dropped: BTreeSet<_> = dropped.into_iter().collect();
    check_coverage("V", kept, &dropped)?;
    warn(&mut warnings, "V", &dropped);

    let (j, kept, dropped) = ctx.table(&model.j_usage, Segment::J);
    let dropped: BTreeSet<_> = dropped.into_iter().collect();
    check_coverage("J", kept, &dropped)?;
    warn(&mut warnings, "J", &dropped);

    let mut d_given_j = BTreeMap::new();
    if !model.d_given_j.is_empty() {
        let kept_j: BTreeSet<&str> = j.iter().map(|(a, _)| a.gene.as_str()).collect();
        let (mut dropped, mut mass, mut weight) = (BTreeSet::new(), 0.0, 0.0);
        for (jg, row) in &model.d_given_j {
            if !kept_j.contains(jg.as_str()) {
                continue;
            }
            let (d, kept, dr) = ctx.table(row, Segment::D);
            if d.is_empty() {
                return Err(GenerateError::ModelInvalid(format!(
                    "d_given_j.{jg}: none of its D genes ({}) have alleles in the germline database",
                    row.keys().cloned().collect::<Vec<_>>().join(", ")
                )));
            }
            let w = model.j_usage.get(jg).copied().unwrap_or(0.0);
            mass += w * kept;
            weight += w;
            dropped.extend(dr);
            d_given_j.insert(jg.clone(), d);
        }
        if let Some(jg) = kept_j.iter().find(|g| !d_given_j.contains_key(**g)) {
            return Err(GenerateError::ModelInvalid(format!("d_given_j.{jg}: missing row for J gene {jg}")));
        }
        check_coverage("D", if weight > 0.0 { mass / weight } else { 1.0 }, &dropped)?;
        warn(&mut warnings, "D", &dropped);
    }

    Ok(Bound { v, d_given_j, j, candidates, warnings, model })
}
