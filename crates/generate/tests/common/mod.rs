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

/// Random stop-free codons (deterministic in `seed`).
fn codons(n: usize, seed: u64) -> Vec<String> {
    use rand::{Rng, SeedableRng};
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);
    let b = b"ACGT";
    (0..n)
        .map(|_| loop {
            let c: String = (0..3).map(|_| b[rng.random_range(0..4)] as char).collect();
            if !["TAA", "TAG", "TGA"].contains(&c.as_str()) {
                break c;
            }
        })
        .collect()
}

/// IMGT-gapped V: 104 codons with gap codons at IMGT 10, 31-34 and 60-61, Cys104 = TGT
/// at gapped 309..312, then `tail`. Seed picks the sequence; `snp` changes one base.
fn sim_v(seed: u64, tail: &str, snp: Option<usize>) -> String {
    let mut c = codons(104, seed);
    for i in [9, 30, 31, 32, 33, 59, 60] {
        c[i] = "...".into();
    }
    c[103] = "TGT".into();
    let mut s = c.concat() + tail;
    if let Some(p) = snp {
        let b = if &s[p..p + 1] == "A" { "G" } else { "A" };
        s.replace_range(p..p + 1, b);
    }
    s
}

fn sim(name: &str, locus: Locus, segment: Segment, gapped: String) -> Allele {
    let ungapped: String = gapped.chars().filter(|&c| c != '.').collect();
    let anchor = match segment {
        Segment::V => abstar_germline::find_v_anchor(&gapped),
        Segment::J => abstar_germline::find_j_anchor(&ungapped),
        _ => None,
    };
    Allele { anchor, sequence_gapped: gapped, ..allele(name, locus, segment, "human", false) }
}

/// Realistic-length human IGH alleles with the gene names of tiny_igh.model.toml:
/// IMGT-gapped stop-free V (IGHV1-1*01/*02 differ by one base), real-sequence D and J,
/// and two C genes (IGHM, IGHG1).
pub fn sim_igh_alleles() -> Vec<Allele> {
    use Segment::*;
    let h = |n, s, g: String| sim(n, Locus::IGH, s, g);
    let tail = "GCGAGAGA";
    vec![
        h("IGHV1-1*01", V, sim_v(11, tail, None)),
        h("IGHV1-1*02", V, sim_v(11, tail, Some(150))),
        h("IGHV1-2*01", V, sim_v(12, tail, None)),
        h("IGHV2-1*01", V, sim_v(13, "GCACACAG", None)),
        h("IGHD1-1*01", D, "GGTACAACTGGAACGAC".into()),
        h("IGHD2-1*01", D, "AGGATATTGTGGTGGTGACTGCTATTCC".into()),
        h("IGHJ1*01", J, "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG".into()),
        h("IGHJ2*01", J, "CTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG".into()),
        h("IGHM*01", C, "GGGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGATACG".into()),
        h("IGHG1*01", C, "CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACA".into()),
    ]
}

/// Realistic-length human IGK alleles with the gene names of tiny_igk.model.toml.
pub fn sim_igk_alleles() -> Vec<Allele> {
    use Segment::*;
    let k = |n, s, g: String| sim(n, Locus::IGK, s, g);
    vec![
        k("IGKV1-1*01", V, sim_v(21, "CAACAGTATAATAGTTACCC", None)),
        k("IGKJ1*01", J, "GTGGACGTTCGGCCAAGGGACCAAGGTGGAAATCAAAC".into()),
        k("IGKC*01", C, "CGAACTGTGGCTGCACCATCTGTCTTCATCTTCCCGCCATCTGATGAGCAGTTGAAATCTGGAACTGCC".into()),
    ]
}
