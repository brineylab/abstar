mod common;
use abstar_annotate::aa::assign::{assign_aa, AaCalls};
use abstar_annotate::aa::derive::{aa_alignments, aa_masks, aa_mutations, aa_regions, AaRegions};
use abstar_annotate::aa::germline::{translate_v, AaSegment};
use abstar_annotate::align::ProteinScoring;
use abstar_germline::{Database, Segment};
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

fn call(q: &[u8]) -> AaCalls {
    let r = r();
    assign_aa(q, &r.v, &r.j, &r.db.alleles, &ProteinScoring::default()).unwrap()
}

fn regions(q: &[u8]) -> (AaCalls, AaRegions) {
    let c = call(q);
    let reg = aa_regions(q, &c, &r().v, &r().j);
    (c, reg)
}

fn slice(q: &[u8], x: Option<(usize, usize)>) -> &str {
    let (a, b) = x.unwrap();
    std::str::from_utf8(&q[a..b]).unwrap()
}

/// IGHV1-2*02's ungapped residues and their 1-based IMGT positions, from its gapped translation.
fn ighv1_2() -> (Vec<u8>, Vec<usize>) {
    let a = r().db.alleles.iter().find(|a| a.name == "IGHV1-2*02").unwrap();
    let g = translate_v(&a.sequence_gapped).unwrap();
    g.iter().enumerate().filter(|(_, &x)| x != b'.').map(|(k, &x)| (x, k + 1)).unzip()
}

#[test]
fn igh_regions() {
    let q = IGH_AA.as_bytes();
    let (_, reg) = regions(q);
    let got: Vec<&str> = [reg.fwr1, reg.cdr1, reg.fwr2, reg.cdr2, reg.fwr3, reg.cdr3, reg.fwr4].iter().map(|&x| slice(q, x)).collect();
    assert_eq!(
        got,
        [
            "QVQLVQSGAEVKKPGASVKVSCKAS",
            "GYTFTGYY",
            "MHWVRQAPGQGLGWMGW",
            "INPNSGGT",
            "NYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC",
            "ARDRITIFGVVIIPDYFDY",
            "WGQGTLVTVSS"
        ]
    );
    assert_eq!(slice(q, reg.junction), "CARDRITIFGVVIIPDYFDYW");
    assert_eq!(reg.vdj, (0, q.len()));
}

#[test]
fn igh_mutations() {
    let q = IGH_AA.as_bytes();
    let c = call(q);
    let [v, j] = aa_mutations(q, &c, &r().v, &r().j);
    let (v, j) = (v.unwrap(), j.unwrap());
    assert_eq!(v.nt, ["IMGT51:E>G"]);
    assert_eq!(v.count, 1);
    assert!(v.insertions.is_empty() && v.deletions.is_empty());
    assert!(j.nt.is_empty() && j.insertions.is_empty() && j.deletions.is_empty());
    assert_eq!(j.count, 0);
}

#[test]
fn igh_masks() {
    let q = IGH_AA.as_bytes();
    let (c, reg) = regions(q);
    let m = aa_masks(q, &c, &r().v, &r().j, &reg);
    assert_eq!(m.segment, format!("{}{}{}", "V".repeat(98), "N".repeat(13), "J".repeat(15)));
    let cdr = ["0".repeat(25), "1".repeat(8), "0".repeat(17), "2".repeat(8), "0".repeat(38), "3".repeat(19), "0".repeat(11)].concat();
    assert_eq!(m.cdr.unwrap(), cdr);
    let ones: Vec<usize> = m.nongermline.bytes().enumerate().filter(|x| x.1 == b'1').map(|x| x.0).collect();
    assert_eq!(ones.len(), 14);
    assert_eq!(ones.iter().filter(|&&i| (98..111).contains(&i)).count(), 13);
}

#[test]
fn v_insertion_and_deletion() {
    // IGH_AA covers IGHV1-2*02 from its first residue (query index k = germline index k). `GS`
    // goes after the last germline residue at or before IMGT 61 (IMGT 60-61 may be gaps in an
    // 8-residue CDR2), whose IMGT position is `p`; the first FWR3 `E` (`...YMELS...`, IMGT `d`) is
    // deleted. Neither indel can shift: `GS` is not followed by `GS`, and `E`'s neighbours differ
    // from it.
    let (seq, imgt) = ighv1_2();
    let k = imgt.iter().rposition(|&x| x <= 61).unwrap();
    let e = (0..seq.len()).find(|&u| imgt[u] >= 66 && seq[u] == b'E').unwrap();
    let (p, d) = (imgt[k], imgt[e]);
    assert_eq!(IGH_AA.as_bytes()[e], b'E');
    let mut q = IGH_AA.as_bytes().to_vec();
    q.remove(e);
    q.splice(k + 1..k + 1, *b"GS");
    let c = call(&q);
    assert_eq!((c.v.aln.q_start, c.v.aln.g_start), (0, 0));
    let [v, _] = aa_mutations(&q, &c, &r().v, &r().j);
    let v = v.unwrap();
    assert_eq!(v.insertions, [format!("IMGT{p}:ins:GS")]);
    assert_eq!(v.deletions, [format!("IMGT{d}:del:E")]);
}

#[test]
fn projected_cys() {
    // The 8 residues before Cys104 (`SDDTAVYY` -> `WWWWWWPP`) score -25 against IGHV1-2*02, more
    // than the following `CAR` regains (+18), so V's trimmed alignment ends before them and
    // Cys104 is projected 8 residues along its last diagonal. (The worst BLOSUM62 substitutions
    // of only the last 6, `DTAVYY`, total -18: a tie, which the trim resolves to the longer V.)
    let mut q = IGH_AA.as_bytes().to_vec();
    q[87..95].copy_from_slice(b"WWWWWWPP");
    let (c, reg) = regions(&q);
    let g = &r().v.groups[c.v.group];
    let cys = g.anchor.unwrap();
    let last = c.v.aln.g_end - 1;
    assert!(last < cys, "V's alignment reaches Cys104");
    assert!(c.v.aln.q_end - 1 + (cys - last) == 95 && cys - last <= 16);
    assert_eq!(q[95], b'C');
    assert_eq!(reg.cdr3.unwrap().0, 96);
    assert_eq!(slice(&q, reg.junction), "CARDRITIFGVVIIPDYFDYW");
}

#[test]
fn no_cys_no_cdr3() {
    let mut q = IGH_AA.as_bytes().to_vec();
    q.drain(90..98);
    let (c, reg) = regions(&q);
    assert!(c.v.aln.g_end <= r().v.groups[c.v.group].anchor.unwrap());
    assert_eq!(reg.cdr3, None);
    assert_eq!(reg.junction, None);
    assert_eq!(aa_masks(&q, &c, &r().v, &r().j, &reg).cdr, None);
}

#[test]
fn alignment_strings() {
    let q = IGH_AA.as_bytes();
    let c = call(q);
    let (s, g) = aa_alignments(q, &c, &r().v, &r().j);
    assert_eq!(s.replace('.', ""), IGH_AA);
    assert_eq!(s.len(), g.len());
    let n0 = s.find("DRITIF").unwrap();
    let dashes: Vec<usize> = g.bytes().enumerate().filter(|x| x.1 == b'-').map(|x| x.0).collect();
    assert_eq!(dashes, (n0..n0 + 13).collect::<Vec<_>>());
}

#[test]
fn mutated_first_residue_extends_v() {
    // Q>P at IMGT 1 is cut by the trim; the 5' extension (record 0029) brings it back.
    let q = format!("P{}", &IGH_AA[1..]);
    let q = q.as_bytes();
    let (c, reg) = regions(q);
    assert_eq!(c.v.aln.q_start, 0);
    assert_eq!(reg.fwr1, Some((0, 25)));
    let [v, _] = aa_mutations(q, &c, &r().v, &r().j);
    assert!(v.unwrap().nt.contains(&"IMGT1:Q>P".to_string()));
}
