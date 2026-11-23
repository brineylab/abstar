//! Scalar Gotoh overlap aligner for protein sequences (BLOSUM62, affine gaps).
use super::aligner::push;
use super::{Alignment, Op};

/// Affine gap scores; a gap of length n costs `gap_open + (n - 1) * gap_extend`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ProteinScoring {
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for ProteinScoring {
    fn default() -> Self {
        ProteinScoring { gap_open: -11, gap_extend: -1 }
    }
}

const ORDER: &[u8; 24] = b"ARNDCQEGHILKMFPSTWYVBZX*";
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

/// Byte -> `ORDER` index, case-insensitive; any other byte -> 22 (`X`).
const IDX: [u8; 256] = {
    let mut t = [22u8; 256];
    let mut k = 0;
    while k < ORDER.len() {
        t[ORDER[k] as usize] = k as u8;
        t[ORDER[k].to_ascii_lowercase() as usize] = k as u8;
        k += 1;
    }
    t
};

fn idx(c: u8) -> usize {
    IDX[c as usize] as usize
}

/// NCBI BLOSUM62, case-insensitive; `J`, `O`, `U` and any other byte score as `X`.
pub fn blosum62(a: u8, b: u8) -> i32 {
    BLOSUM62[idx(a)][idx(b)] as i32
}

fn is_match(q: u8, g: u8) -> bool {
    q.eq_ignore_ascii_case(&g) && !q.eq_ignore_ascii_case(&b'X')
}

const NEG: i32 = i32::MIN / 4;

/// Gotoh overlap alignment: both query and germline ends are free; terminal gaps are scored
/// then dropped from `ops`. Empty when the optimum is <= 0. Traceback prefers the diagonal, so
/// gaps land leftmost.
// ponytail: scalar O(n·m) per allele; AA k-mer prefilter + block-aligner AAMatrix if protein volume grows
pub fn align_protein(query: &[u8], germline: &[u8], s: &ProteinScoring) -> Alignment {
    let (m, n) = (query.len(), germline.len());
    let w = m + 1;
    let at = |i: usize, j: usize| i * w + j;
    // h: best; e: ends in Del (germline-only); f: ends in Ins (query-only). Row/col 0 free (h=0).
    let mut h = vec![0i32; (n + 1) * w];
    let mut e = vec![NEG; (n + 1) * w];
    let mut f = vec![NEG; (n + 1) * w];
    for i in 1..=n {
        for j in 1..=m {
            e[at(i, j)] = (h[at(i - 1, j)] + s.gap_open).max(e[at(i - 1, j)] + s.gap_extend);
            f[at(i, j)] = (h[at(i, j - 1)] + s.gap_open).max(f[at(i, j - 1)] + s.gap_extend);
            let diag = h[at(i - 1, j - 1)] + blosum62(query[j - 1], germline[i - 1]);
            h[at(i, j)] = diag.max(e[at(i, j)]).max(f[at(i, j)]);
        }
    }
    let (mut score, mut ei, mut ej) = (0, 0, 0);
    for j in 0..=m {
        if h[at(n, j)] > score {
            (score, ei, ej) = (h[at(n, j)], n, j);
        }
    }
    for i in 0..n {
        if h[at(i, m)] > score {
            (score, ei, ej) = (h[at(i, m)], i, m);
        }
    }
    if score <= 0 {
        return Alignment { score: 0, q_start: 0, q_end: 0, g_start: 0, g_end: 0, ops: vec![] };
    }
    let (mut i, mut j, mut state) = (ei, ej, 0u8); // 0 h, 1 e (Del), 2 f (Ins)
    let mut rev = Vec::new();
    while i > 0 && j > 0 {
        let k = at(i, j);
        match state {
            0 if h[k] == h[at(i - 1, j - 1)] + blosum62(query[j - 1], germline[i - 1]) => {
                rev.push(if is_match(query[j - 1], germline[i - 1]) { Op::Match } else { Op::Mismatch });
                i -= 1;
                j -= 1;
            }
            0 if h[k] == e[k] => state = 1,
            0 => state = 2,
            1 => {
                rev.push(Op::Del);
                state = if e[k] == h[at(i - 1, j)] + s.gap_open { 0 } else { 1 };
                i -= 1;
            }
            _ => {
                rev.push(Op::Ins);
                state = if f[k] == h[at(i, j - 1)] + s.gap_open { 0 } else { 2 };
                j -= 1;
            }
        }
    }
    // Walk forward from the traceback start, dropping leading and trailing gaps.
    let mut ops = Vec::new();
    let (mut start, mut end, mut keep) = (None, (0, 0), 0);
    for op in rev.into_iter().rev() {
        match op {
            Op::Del => {
                if start.is_some() {
                    push(&mut ops, Op::Del);
                }
                i += 1;
            }
            Op::Ins => {
                if start.is_some() {
                    push(&mut ops, Op::Ins);
                }
                j += 1;
            }
            _ => {
                start.get_or_insert((j, i));
                push(&mut ops, op);
                i += 1;
                j += 1;
                end = (j, i);
                keep = ops.len();
            }
        }
    }
    let Some((q_start, g_start)) = start else {
        return Alignment { score: 0, q_start: 0, q_end: 0, g_start: 0, g_end: 0, ops: vec![] };
    };
    ops.truncate(keep);
    Alignment { score, q_start, q_end: end.0, g_start, g_end: end.1, ops }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn blosum62_is_symmetric() {
        for a in 0..24 {
            for b in 0..24 {
                assert_eq!(BLOSUM62[a][b], BLOSUM62[b][a]);
            }
        }
    }

    #[test]
    fn idx_matches_a_case_insensitive_scan() {
        for c in 0..=255u8 {
            assert_eq!(idx(c), ORDER.iter().position(|&o| o == c.to_ascii_uppercase()).unwrap_or(22), "{c}");
        }
    }
}
