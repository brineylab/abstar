#![allow(dead_code)]
use abstar_germline::{Allele, Database, Functionality, Locus, Manifest, Segment};

pub fn allele(name: &str, locus: Locus, segment: Segment, species: &str, anchored: bool) -> Allele {
    let gene = name.split('*').next().unwrap().to_string();
    let seq = "ACGTACGTACGTACGTACGTACGTTGTGCGCGT";
    let gapped = if segment == Segment::V { "ACGTACGTAC..GTACGTACGTACGTTGTGCGCGT" } else { seq };
    Allele {
        id: name.to_string(),
        name: name.to_string(),
        gene,
        locus,
        segment,
        species: species.to_string(),
        functionality: Functionality::F,
        sequence_gapped: gapped.to_string(),
        anchor: anchored.then_some(24),
        source: "test".into(),
        usage: None,
    }
}

pub fn db(alleles: Vec<Allele>) -> Database {
    Database {
        manifest: Manifest {
            format_version: 1,
            name: "tiny".into(),
            description: String::new(),
            abstar_version: "0".into(),
            built: "2026-11-20".into(),
            alleles_sha256: String::new(),
            species: vec!["human".into()],
            loci: vec![Locus::IGH],
            prefixes: Default::default(),
            provenance: vec![],
            warnings: vec![],
        },
        alleles,
    }
}

/// Alleles matching tests/data/tiny_igh.model.toml, human IGH.
pub fn tiny_igh_alleles() -> Vec<Allele> {
    use Segment::*;
    let h = |n, s| allele(n, Locus::IGH, s, "human", true);
    vec![
        h("IGHV1-1*01", V), h("IGHV1-1*02", V), h("IGHV1-2*01", V), h("IGHV2-1*01", V),
        h("IGHD1-1*01", D), h("IGHD2-1*01", D),
        h("IGHJ1*01", J), h("IGHJ2*01", J),
        allele("IGHM*01", Locus::IGH, C, "human", false),
    ]
}

pub fn tiny_igh_db() -> Database {
    db(tiny_igh_alleles())
}
