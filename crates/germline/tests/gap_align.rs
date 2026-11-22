use abstar_germline::gap::align::*;

fn pair_sum(q: &[u8], t: &[u8], cols: &[Col]) -> i32 {
    cols.iter()
        .map(|c| match *c {
            Col::Pair(i, j) => blosum62(q[i], t[j]),
            _ => 0,
        })
        .sum()
}

#[test]
fn identical_sequences_align_fully() {
    let s = b"QVQLVQSGAEVKKPGASVKVSC";
    let a = align(s, s, &[], &Scoring::DEFAULT);
    assert_eq!(a.cols, (0..22).map(|i| Col::Pair(i, i)).collect::<Vec<_>>());
    assert_eq!(a.score, s.iter().map(|&c| blosum62(c, c)).sum::<i32>());
}

#[test]
fn end_gaps_are_free() {
    let a = align(
        b"VKKPGASVKV",
        b"QVQLVQSGAEVKKPGASVKVSCKAS",
        &[],
        &Scoring::DEFAULT,
    );
    assert_eq!(
        a.cols,
        (0..10).map(|i| Col::Pair(i, i + 10)).collect::<Vec<_>>()
    );
}

#[test]
fn internal_deletion() {
    let t = b"QVQLVQSGAEVKKPGASVKVSCKASGYTFT";
    let mut q = t[..14].to_vec();
    q.extend_from_slice(&t[16..]);
    let a = align(&q, t, &[], &Scoring::DEFAULT);
    let runs = a
        .cols
        .iter()
        .filter(|c| matches!(c, Col::TemplateOnly(_)))
        .count();
    assert_eq!(runs, 2);
    assert!(a
        .cols
        .windows(2)
        .any(|w| matches!((w[0], w[1]), (Col::TemplateOnly(_), Col::TemplateOnly(_)))));
    assert!(!a.cols.iter().any(|c| matches!(c, Col::QueryOnly(_))));
    let s = Scoring::DEFAULT;
    assert_eq!(
        a.score,
        pair_sum(&q, t, &a.cols) - s.gap_open - s.gap_extend
    );
}

#[test]
fn blosum62_spot_values() {
    assert_eq!(blosum62(b'W', b'W'), 11);
    assert_eq!(blosum62(b'C', b'C'), 9);
    assert_eq!(blosum62(b'A', b'R'), -1);
    assert_eq!(blosum62(b'*', b'A'), -4);
    assert_eq!(blosum62(b'X', b'A'), -1);
}

#[test]
fn gap_in_cdr_columns_opens_cheaper() {
    let t = b"QVQLVQSGAEVKKPGASVKVSCKASGYTFT";
    let mut q = t[..14].to_vec();
    q.extend_from_slice(&t[16..]);
    let cdr: Vec<bool> = (0..t.len()).map(|j| (14..16).contains(&j)).collect();
    let s = Scoring::DEFAULT;
    let a = align(&q, t, &cdr, &s);
    assert!(a.cols.contains(&Col::TemplateOnly(14)) && a.cols.contains(&Col::TemplateOnly(15)));
    assert_eq!(
        a.score,
        pair_sum(&q, t, &a.cols) - s.cdr_gap_open - s.gap_extend
    );
    assert!(s.cdr_gap_open < s.gap_open);
}
