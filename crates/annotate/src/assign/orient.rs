use crate::align::KmerIndex;

/// The read on the strand chosen by `orient`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Oriented {
    /// The uppercased read, reverse-complemented when `rev_comp`.
    pub seq: Vec<u8>,
    pub rev_comp: bool,
}

/// IUPAC complement; other bytes are kept.
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        o => o,
    }
}

/// Picks the strand with more k-mer evidence, the best V allele's hits plus the best J allele's
/// hits (record 0029); a tie goes forward. `None` (no V hit) when neither strand shares a k-mer
/// with any V allele.
pub fn orient(read: &[u8], v_index: &KmerIndex, j_index: &KmerIndex) -> Option<Oriented> {
    let fwd = read.to_ascii_uppercase();
    let rev: Vec<u8> = fwd.iter().rev().map(|&b| complement(b)).collect();
    let best = |index: &KmerIndex, s: &[u8]| index.hits(s).iter().map(|&(_, n)| n).max().unwrap_or(0);
    let (vf, vr) = (best(v_index, &fwd), best(v_index, &rev));
    if vf == 0 && vr == 0 {
        return None;
    }
    Some(if vr + best(j_index, &rev) > vf + best(j_index, &fwd) {
        Oriented { seq: rev, rev_comp: true }
    } else {
        Oriented { seq: fwd, rev_comp: false }
    })
}
