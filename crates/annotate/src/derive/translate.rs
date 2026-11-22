/// Translates `nt[frame..]` with the standard genetic code. A partial terminal codon is
/// dropped; a codon with a base other than ACGT (either case) is `X`.
pub fn translate(nt: &[u8], frame: usize) -> String {
    const AA: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
    let code = |b: u8| match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    };
    nt.get(frame..)
        .unwrap_or_default()
        .as_chunks::<3>()
        .0
        .iter()
        .map(|c| match (code(c[0]), code(c[1]), code(c[2])) {
            (Some(a), Some(b), Some(d)) => AA[a * 16 + b * 4 + d] as char,
            _ => 'X',
        })
        .collect()
}

/// Translation of read region `a..b` in the reading frame whose codons start at `frame`
/// (`frame` and every third position after it): the codons whose third base lies in the
/// region, so a codon crossing a region boundary belongs to the downstream region.
pub fn region_aa(seq: &[u8], frame: usize, (a, b): (usize, usize)) -> String {
    let first = frame + a.saturating_sub(frame + 2).div_ceil(3) * 3;
    translate(&seq[first.min(b)..b], 0)
}
