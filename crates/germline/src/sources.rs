use crate::{Allele, Functionality, GermlineError, Locus, Segment};
use serde::Deserialize;

fn allele(name: &str, locus: Locus, segment: Segment, functionality: Functionality, seq: String) -> Allele {
    Allele {
        id: name.into(), // replaced once the species mix is known
        name: name.into(),
        gene: name.split('*').next().unwrap_or(name).into(),
        locus,
        segment,
        species: String::new(),
        functionality,
        sequence_gapped: seq,
        anchor: None,
        source: String::new(),
        usage: None,
    }
}

#[derive(Deserialize)]
struct OgrdbFile {
    #[serde(rename = "GermlineSet")]
    sets: Vec<GermlineSet>,
}

#[derive(Deserialize)]
struct GermlineSet {
    germline_set_name: String,
    release_version: serde_json::Number,
    release_date: String,
    allele_descriptions: Vec<AlleleDescription>,
}

#[derive(Deserialize)]
struct AlleleDescription {
    label: String,
    sequence_type: String,
    locus: String,
    functional: Option<bool>,
    coding_sequence: Option<String>,
    v_gene_delineations: Option<Vec<Delineation>>,
    j_cdr3_end: Option<usize>,
}

#[derive(Deserialize)]
struct Delineation {
    delineation_scheme: Option<String>,
    aligned_sequence: Option<String>,
}

/// Alleles from an OGRDB germline-set JSON, plus one description per germline set.
/// J alleles carry the OGRDB anchor (`j_cdr3_end - 1`) in `anchor`.
pub(crate) fn ogrdb(bytes: &[u8], file: &str) -> Result<(Vec<Allele>, Vec<String>), GermlineError> {
    let bad = |m: String| GermlineError::RecipeInvalid(format!("{file}: {m}"));
    let parsed: OgrdbFile = serde_json::from_slice(bytes).map_err(|e| bad(e.to_string()))?;
    let mut alleles = Vec::new();
    let mut descriptions = Vec::new();
    for set in parsed.sets {
        descriptions.push(format!(
            "OGRDB {} release {} ({})",
            set.germline_set_name, set.release_version, set.release_date
        ));
        for d in set.allele_descriptions {
            let bad = |m: &str| bad(format!("allele {}: {m}", d.label));
            let segment: Segment = d.sequence_type.parse().map_err(|_| bad("unknown sequence_type"))?;
            let locus: Locus = d.locus.parse().map_err(|_| bad("unknown locus"))?;
            let seq = match segment {
                Segment::V => d
                    .v_gene_delineations
                    .as_deref()
                    .and_then(<[_]>::first)
                    .filter(|x| x.delineation_scheme.as_deref() == Some("IMGT"))
                    .and_then(|x| x.aligned_sequence.clone())
                    .ok_or_else(|| bad("no IMGT v_gene_delineations[0].aligned_sequence"))?,
                _ => d.coding_sequence.clone().ok_or_else(|| bad("no coding_sequence"))?,
            };
            let func = if d.functional == Some(true) { Functionality::F } else { Functionality::NF };
            let mut a = allele(&d.label, locus, segment, func, seq);
            if segment == Segment::J {
                a.anchor = match d.j_cdr3_end {
                    Some(0) => return Err(bad("j_cdr3_end is 0")),
                    end => end.map(|e| e - 1),
                };
            }
            alleles.push(a);
        }
    }
    Ok((alleles, descriptions))
}

/// Alleles from a (gapped) FASTA; the segment comes from the recipe key, never the name.
pub(crate) fn fasta(bytes: &[u8], file: &str, segment: Segment) -> Result<Vec<Allele>, GermlineError> {
    let bad = |m: String| GermlineError::RecipeInvalid(format!("{file}: {m}"));
    let mut reader = needletail::parse_fastx_reader(bytes).map_err(|e| bad(e.to_string()))?;
    let mut out = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| bad(e.to_string()))?;
        let header = String::from_utf8_lossy(rec.id()).into_owned();
        let (name, func) = if header.contains('|') {
            let f: Vec<&str> = header.split('|').collect();
            let func = match f.get(3).map(|s| s.trim().replace(['(', ')', '[', ']'], "")).as_deref() {
                Some("F") => Functionality::F,
                Some("ORF") => Functionality::ORF,
                Some("P") => Functionality::P,
                _ => Functionality::Unknown,
            };
            (f.get(1).map_or("", |s| s.trim()), func)
        } else {
            (header.split_whitespace().next().unwrap_or(""), Functionality::Unknown)
        };
        if name.is_empty() {
            return Err(bad(format!("no allele name in header {header:?}")));
        }
        let locus: Locus = name.get(..3).and_then(|l| l.parse().ok()).ok_or_else(|| {
            bad(format!("allele {name}: name does not start with a known locus"))
        })?;
        let seq = String::from_utf8_lossy(&rec.seq()).into_owned();
        out.push(allele(name, locus, segment, func, seq));
    }
    Ok(out)
}
