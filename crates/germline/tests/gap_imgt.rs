use abstar_germline::gap::imgt::*;
use abstar_germline::{load, Segment};
use std::path::PathBuf;

fn db_dir(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../python/abstar/germline_dbs")
        .join(name)
}

#[test]
fn cdr_rule_table() {
    assert_eq!(
        cdr_positions(Cdr::One, 8),
        Some(vec![27, 28, 29, 30, 35, 36, 37, 38])
    );
    assert_eq!(cdr_positions(Cdr::One, 12), Some((27..=38).collect()));
    assert_eq!(
        cdr_positions(Cdr::One, 11),
        Some((27..=38).filter(|&p| p != 33).collect())
    );
    assert_eq!(cdr_positions(Cdr::Two, 3), Some(vec![56, 57, 65]));
    assert_eq!(
        cdr_positions(Cdr::Two, 8),
        Some((56..=65).filter(|&p| p != 60 && p != 61).collect())
    );
    assert_eq!(cdr_positions(Cdr::One, 13), None);
    assert_eq!(cdr_positions(Cdr::Two, 11), None);
    assert_eq!(cdr_positions(Cdr::One, 0), Some(vec![]));
}

#[test]
fn cdr_rule_matches_shipped_databases() {
    let mut checked = 0;
    for name in ["human", "c57bl6", "balbc"] {
        let db = load(&db_dir(name)).unwrap();
        for a in db
            .alleles
            .iter()
            .filter(|a| a.segment == Segment::V && a.anchor.is_some())
        {
            let s = &a.sequence_gapped;
            for (cdr, lo, hi, dots_before) in
                [(Cdr::One, 27u16, 38u16, 114), (Cdr::Two, 56, 65, 195)]
            {
                // complete: starts with a base (not 5'-partial) and no mixed codons through the CDR
                let head = &s.as_bytes()[..dots_before.min(s.len())];
                if s.len() < 3 * hi as usize
                    || s.starts_with('.')
                    || head.chunks(3).any(|c| c.contains(&b'.') && c != b"...")
                {
                    continue;
                }
                let got: Vec<u16> = template_residues(s)
                    .into_iter()
                    .map(|r| r.pos)
                    .filter(|p| (lo..=hi).contains(p))
                    .collect();
                assert_eq!(
                    Some(got.clone()),
                    cdr_positions(cdr, got.len()),
                    "{}",
                    a.name
                );
                checked += 1;
            }
        }
    }
    assert!(checked >= 1000, "only {checked} checked");
}

#[test]
fn template_residues_human_ighv1_2() {
    let db = load(&db_dir("human")).unwrap();
    let a = db.alleles.iter().find(|a| a.name == "IGHV1-2*02").unwrap();
    let rs = template_residues(&a.sequence_gapped);
    let at = |p: u16| rs.iter().find(|r| r.pos == p);
    assert_eq!(at(23).unwrap().aa, b'C');
    assert_eq!(at(41).unwrap().aa, b'W');
    assert_eq!(at(104).unwrap().aa, b'C');
    assert_eq!(at(104).unwrap().nt, 309..312);
    assert!(at(10).is_none() && at(73).is_none());
}

#[test]
fn template_residues_skip_partial_codons() {
    let rs = template_residues(&format!("{}ACGTGTACG", ".".repeat(16)));
    assert_eq!((rs[0].pos, rs[0].nt.clone()), (7, 18..21));
}
