//! IMGT position helpers for V gapping.
use std::ops::Range;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cdr {
    One,
    Two,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Region {
    Fr1,
    Cdr1,
    Fr2,
    Cdr2,
    Fr3,
    Tail,
}

/// CDR1 positions in the order gaps open (from the centre).
const CDR1_GAP_ORDER: [u16; 12] = [33, 32, 34, 31, 35, 30, 36, 29, 37, 28, 38, 27];
const CDR2_GAP_ORDER: [u16; 10] = [61, 60, 62, 59, 63, 58, 64, 57, 65, 56];

/// IMGT positions a CDR of `len` residues occupies, ascending; `None` if it cannot be represented.
pub fn cdr_positions(cdr: Cdr, len: usize) -> Option<Vec<u16>> {
    let (full, order): (Range<u16>, &[u16]) = match cdr {
        Cdr::One => (27..39, &CDR1_GAP_ORDER),
        Cdr::Two => (56..66, &CDR2_GAP_ORDER),
    };
    let n_gaps = order.len().checked_sub(len)?;
    Some(full.filter(|p| !order[..n_gaps].contains(p)).collect())
}

/// Standard-code translation of one codon; non-ACGT gives `X`, stops give `*`.
pub fn translate(codon: &[u8]) -> u8 {
    const AA: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
    let mut idx = 0;
    for &b in codon {
        let v = match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return b'X',
        };
        idx = idx * 4 + v;
    }
    if codon.len() == 3 {
        AA[idx]
    } else {
        b'X'
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Residue {
    pub pos: u16,
    pub aa: u8,
    /// Range in the gapped string.
    pub nt: Range<usize>,
}

/// One residue per IMGT position whose codon is three bases; gap and partial-gap codons are skipped.
pub fn template_residues(gapped: &str) -> Vec<Residue> {
    let b = gapped.as_bytes();
    (1..=b.len() / 3)
        .filter_map(|p| {
            let nt = 3 * (p - 1)..3 * p;
            let codon = &b[nt.clone()];
            (!codon.contains(&b'.')).then(|| Residue {
                pos: p as u16,
                aa: translate(codon),
                nt,
            })
        })
        .collect()
}

pub fn region(pos: u16) -> Region {
    match pos {
        ..=26 => Region::Fr1,
        27..=38 => Region::Cdr1,
        39..=55 => Region::Fr2,
        56..=65 => Region::Cdr2,
        66..=104 => Region::Fr3,
        _ => Region::Tail,
    }
}
