mod common;
use abstar_annotate::aa::assign::{assign_aa, AaCalls};
use abstar_annotate::aa::germline::AaSegment;
use abstar_annotate::align::ProteinScoring;
use abstar_germline::{Allele, Database, Segment};
use std::sync::OnceLock;

const IGH_AA: &str = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLGWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARDRITIFGVVIIPDYFDYWGQGTLVTVSS";

struct Ref {
    db: Database,
    v: AaSegment,
    j: AaSegment,
}

fn r() -> &'static Ref {
    static R: OnceLock<Ref> = OnceLock::new();
    R.get_or_init(|| {
        let db = common::human_db();
        let v = AaSegment::new(&db.alleles, Segment::V, None, false);
        let j = AaSegment::new(&db.alleles, Segment::J, None, false);
        Ref { db, v, j }
    })
}

fn call(q: &[u8]) -> Option<AaCalls> {
    let r = r();
    assign_aa(q, &r.v, &r.j, &r.db.alleles, &ProteinScoring::default())
}

fn a(i: usize) -> &'static Allele {
    &r().db.alleles[i]
}

fn names(ids: &[usize]) -> Vec<&'static str> {
    ids.iter().map(|&i| a(i).name.as_str()).collect()
}

fn seq_of(seg: &AaSegment, name: &str) -> Vec<u8> {
    seg.groups.iter().find(|g| g.alleles.iter().any(|&i| a(i).name == name)).unwrap().seq.clone()
}

#[test]
fn igh_calls() {
    let c = call(IGH_AA.as_bytes()).unwrap();
    assert_eq!(a(c.v.primary).gene, "IGHV1-2");
    assert!(names(&c.v.ties).contains(&"IGHV1-2*02"), "{:?}", names(&c.v.ties));
    assert_eq!(c.v.ties[0], c.v.primary);
    let j = c.j.as_ref().unwrap();
    assert_eq!(a(j.primary).gene, "IGHJ4");
    assert!(names(&j.ties).contains(&"IGHJ4*02"), "{:?}", names(&j.ties));
    assert_eq!(j.ties[0], j.primary);
    assert_eq!((c.v_end, c.j_start), (98, 111));
}

#[test]
fn synonymous_alleles_tie() {
    let c = call(IGH_AA.as_bytes()).unwrap();
    let g = &r().v.groups[c.v.group];
    assert!(g.alleles.contains(&c.v.primary));
    for i in &g.alleles {
        assert!(c.v.ties.contains(i), "{} not in ties {:?}", a(*i).name, names(&c.v.ties));
    }
}

#[test]
fn no_j() {
    let c = call(&IGH_AA.as_bytes()[..96]).unwrap();
    assert!(c.j.is_none());
    assert_eq!(c.j_start, 96);
}

#[test]
fn unrelated_protein_no_v() {
    assert!(call(b"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG").is_none());
}

#[test]
fn overlap_cut_earliest() {
    // IGHV3-15*01 ends `...YYCTT`; every IGHJ6 allele ends `...GTTVTVSS`. The read is the whole
    // V germline translation followed by `VTVSS`, so its last V residues `TT` are also the
    // `TT` of J6's `TTVTVSS`. J is searched from the Cys (`CTTVTVSS`); J6 aligns `TTVTVSS`
    // (the `C` against `G` is trimmed), so its alignment starts 2 residues before V's ends.
    // Each `T` scores T/T whichever segment takes it, so every cut in n-2..=n scores the same
    // and the earliest, n-2, wins.
    let v = seq_of(&r().v, "IGHV3-15*01");
    let j6 = seq_of(&r().j, "IGHJ6*02");
    assert!(v.ends_with(b"YYCTT"));
    assert!(j6.ends_with(b"GTTVTVSS"));
    let n = v.len();
    let mut q = v.clone();
    q.extend_from_slice(b"VTVSS");
    let c = call(&q).unwrap();
    let j = c.j.as_ref().unwrap();
    assert_eq!(a(c.v.primary).gene, "IGHV3-15");
    assert_eq!(a(j.primary).gene, "IGHJ6");
    assert_eq!((c.v.aln.q_end, j.aln.q_start), (n, n - 2));
    assert_eq!((c.v_end, c.j_start), (n - 2, n - 2));
}

#[test]
fn leader_not_pulled_into_v() {
    let q = format!("MDWTWRILFLVAAATGAHS{IGH_AA}");
    assert_eq!(call(q.as_bytes()).unwrap().v.aln.q_start, 19);
}

#[test]
fn long_extension_without_identity_not_taken() {
    // QVQL -> WWWW: every column scores below zero, so the trim cuts all four, and a
    // 4-residue extension with no residue identical to the allele is not taken.
    let q = format!("WWWW{}", &IGH_AA[4..]);
    let c = call(q.as_bytes()).unwrap();
    assert_eq!(a(c.v.primary).gene, "IGHV1-2");
    assert_eq!(c.v.aln.q_start, 4);
}
