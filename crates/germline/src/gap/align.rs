//! Semi-global affine-gap protein aligner (BLOSUM62) used to place V germline residues.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Scoring {
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Scoring {
    /// A gap of length k costs `open + (k-1) * extend`.
    pub const DEFAULT: Scoring = Scoring {
        gap_open: 11,
        gap_extend: 1,
    };
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Col {
    Pair(usize, usize),
    QueryOnly(usize),
    TemplateOnly(usize),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinAlignment {
    pub score: i32,
    pub cols: Vec<Col>,
}

const SYMBOLS: &[u8; 24] = b"ARNDCQEGHILKMFPSTWYVBZX*";

#[rustfmt::skip]
const BLOSUM62: [[i8; 24]; 24] = [
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4],
    [-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4],
    [-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
    [ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4],
    [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1],
];

fn idx(c: u8) -> usize {
    SYMBOLS
        .iter()
        .position(|&s| s == c)
        .unwrap_or(SYMBOLS.len() - 2) // X
}

/// Standard BLOSUM62, except that `X` scores -1 against every residue (as in EMBOSS), including
/// `A`, `S`, `T` and `P`, which NCBI scores 0 or -2; `X`/`*` stays -4.
pub fn blosum62(a: u8, b: u8) -> i32 {
    let (i, j) = (idx(a), idx(b));
    if (i == 22 || j == 22) && i != 23 && j != 23 {
        return -1;
    }
    BLOSUM62[i][j] as i32
}

const NEG: i32 = i32::MIN / 2;

#[derive(Clone, Copy, PartialEq)]
enum St {
    M,
    X, // TemplateOnly column
    Y, // QueryOnly column
    Start,
}

/// Gotoh DP, free end gaps on both sequences. Ties prefer diagonal, then a gap in the query
/// (`TemplateOnly`), then `QueryOnly`; tracing back from the end keeps pairs as long as
/// possible, which places gaps leftmost.
pub fn align(query: &[u8], template: &[u8], s: &Scoring) -> ProteinAlignment {
    let (n, m) = (query.len(), template.len());
    let w = m + 1;
    let at = |i: usize, j: usize| i * w + j;
    let size = (n + 1) * w;
    // M: ends in a pair; X: ends in TemplateOnly; Y: ends in QueryOnly.
    let (mut mm, mut xx, mut yy) = (vec![NEG; size], vec![NEG; size], vec![NEG; size]);
    for i in 1..=n {
        for j in 1..=m {
            let prev = at(i - 1, j - 1);
            let mut best = if i == 1 || j == 1 { 0 } else { NEG };
            best = best.max(mm[prev]).max(xx[prev]).max(yy[prev]);
            mm[at(i, j)] = best + blosum62(query[i - 1], template[j - 1]);
            let l = at(i, j - 1);
            xx[at(i, j)] = (mm[l] - s.gap_open).max(xx[l] - s.gap_extend);
            let u = at(i - 1, j);
            yy[at(i, j)] = (mm[u] - s.gap_open).max(yy[u] - s.gap_extend);
        }
    }
    // Best end: last row or last column, ending in a pair.
    let mut end = None;
    let mut score = NEG;
    for i in 1..=n {
        for j in 1..=m {
            if (i == n || j == m) && mm[at(i, j)] > score {
                score = mm[at(i, j)];
                end = Some((i, j));
            }
        }
    }
    let Some((mut i, mut j)) = end else {
        return ProteinAlignment {
            score: 0,
            cols: Vec::new(),
        };
    };
    let mut cols = Vec::new();
    let mut st = St::M;
    loop {
        match st {
            St::M => {
                cols.push(Col::Pair(i - 1, j - 1));
                let target = mm[at(i, j)] - blosum62(query[i - 1], template[j - 1]);
                let prev = at(i - 1, j - 1);
                st = if (i > 1 && j > 1) && mm[prev] == target {
                    St::M
                } else if (i > 1 && j > 1) && xx[prev] == target {
                    St::X
                } else if (i > 1 && j > 1) && yy[prev] == target {
                    St::Y
                } else {
                    St::Start
                };
                if st != St::Start {
                    i -= 1;
                    j -= 1;
                }
            }
            St::X => {
                cols.push(Col::TemplateOnly(j - 1));
                let l = at(i, j - 1);
                st = if xx[at(i, j)] == mm[l] - s.gap_open {
                    St::M
                } else {
                    St::X
                };
                j -= 1;
            }
            St::Y => {
                cols.push(Col::QueryOnly(i - 1));
                let u = at(i - 1, j);
                st = if yy[at(i, j)] == mm[u] - s.gap_open {
                    St::M
                } else {
                    St::Y
                };
                i -= 1;
            }
            St::Start => break,
        }
    }
    cols.reverse();
    ProteinAlignment { score, cols }
}
