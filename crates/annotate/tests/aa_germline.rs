mod common;
use abstar_annotate::aa::germline::{translate_j, translate_v, AaSegment};
use abstar_germline::{Locus, Segment};

fn find<'a>(db: &'a abstar_germline::Database, name: &str) -> (usize, &'a abstar_germline::Allele) {
    db.alleles.iter().enumerate().find(|(_, a)| a.name == name).unwrap()
}

#[test]
fn translate_v_keeps_imgt_numbering() {
    let db = common::human_db();
    let (_, a) = find(&db, "IGHV1-2*02");
    let aa = translate_v(&a.sequence_gapped).unwrap();
    assert_eq!(aa[103], b'C');
    let g = a.sequence_gapped.as_bytes();
    for (i, &r) in aa.iter().enumerate() {
        assert_eq!(r == b'.', &g[i * 3..i * 3 + 3] == b"...");
    }
    let ungapped: String = aa.iter().filter(|&&r| r != b'.').map(|&r| r as char).collect();
    assert!(ungapped.starts_with("QVQLVQSGAEVKKPGASVKVSCKAS"));
}

#[test]
fn translate_v_rejects_ungapped() {
    assert!(translate_v("CAGGTGCAG").is_none());
}

#[test]
fn translate_j_from_anchor() {
    let db = common::human_db();
    let (_, a) = find(&db, "IGHJ4*02");
    let (aa, idx) = translate_j(&a.sequence(), a.anchor.unwrap());
    assert_eq!(String::from_utf8(aa.clone()).unwrap(), "YFDYWGQGTLVTVSS");
    assert_eq!(idx, 4);
    assert_eq!(aa[idx], b'W');
}

#[test]
fn synonymous_alleles_share_a_group() {
    let db = common::human_db();
    let seg = AaSegment::new(&db.alleles, Segment::V, Some(&[Locus::IGH]), false);
    let (i, a) = find(&db, "IGHV1-2*02");
    let g = seg.groups.iter().find(|g| g.alleles.contains(&i)).unwrap();
    let want = translate_v(&a.sequence_gapped).unwrap();
    let ungapped: Vec<u8> = want.iter().copied().filter(|&r| r != b'.').collect();
    assert_eq!(g.seq, ungapped);
    for (j, b) in db.alleles.iter().enumerate() {
        if b.name.starts_with("IGHV1-2*") && !b.functionality.is_nonfunctional() {
            assert_eq!(g.alleles.contains(&j), translate_v(&b.sequence_gapped) == Some(want.clone()), "{}", b.name);
        }
    }
    for g in &seg.groups {
        assert!(g.alleles.iter().all(|&j| db.alleles[j].locus == g.locus));
        assert!(g.alleles.windows(2).all(|w| w[0] < w[1]));
    }
    let g2 = seg.groups.iter().find(|g| g.alleles.contains(&i)).unwrap();
    assert_eq!(g2.seq[g2.anchor.unwrap()], b'C');
}

#[test]
fn nonfunctional_excluded_by_default() {
    let db = common::human_db();
    let (i, _) = find(&db, "IGKV2-29*01");
    let has = |nf| AaSegment::new(&db.alleles, Segment::V, Some(&[Locus::IGK]), nf).groups.iter().any(|g| g.alleles.contains(&i));
    assert!(!has(false));
    assert!(has(true));
}
