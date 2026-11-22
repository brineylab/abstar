//! Portable backend: Gotoh overlap DP with the SIMD backend's semantics (free read ends, free
//! germline ends, terminal gaps penalized then dropped, empty when the optimum is <= 0) and its
//! tie-breaking (end cell and traceback), so both backends give the same alignment.
use super::aligner::{is_match, push, sim_dims};
use super::{Alignment, Op, Scoring};

const NEG: i32 = i32::MIN / 4;

/// Overlap DP, rows = germline, columns = query. Returns (optimum, end row, end column). With
/// `tr`, records per cell: bit 0 best == Ins state, bit 1 best == Del state, bit 2 Ins opened
/// here, bit 3 Del opened here.
fn dp(query: &[u8], germline: &[u8], s: &Scoring, mut tr: Option<&mut Vec<u8>>) -> (i32, usize, usize) {
    let (n, m) = (germline.len(), query.len());
    let (o, e) = (s.gap_open as i32, s.gap_extend as i32);
    // best and Del (vertical) values of the previous and current rows; row 0 and column 0 are
    // free starts (0)
    let (mut best_p, mut del_p) = (vec![0i32; m + 1], vec![NEG; m + 1]);
    let (mut best_c, mut del_c) = (vec![0i32; m + 1], vec![NEG; m + 1]);
    let mut last_col = vec![0i32; n + 1];
    if let Some(t) = tr.as_deref_mut() {
        t.clear();
        t.resize((n + 1) * (m + 1), 0);
    }
    for i in 1..=n {
        let mut ins = NEG;
        for j in 1..=m {
            let sub = if is_match(query[j - 1], germline[i - 1]) { s.matched } else { s.mismatch } as i32;
            let diag = best_p[j - 1] + sub;
            let (d_open, d_ext) = (best_p[j] + o, del_p[j] + e);
            let (i_open, i_ext) = (best_c[j - 1] + o, ins + e);
            let del = d_open.max(d_ext);
            ins = i_open.max(i_ext);
            let b = diag.max(del).max(ins);
            best_c[j] = b;
            del_c[j] = del;
            if let Some(t) = tr.as_deref_mut() {
                t[i * (m + 1) + j] = (ins == b) as u8 | ((del == b) as u8) << 1 | ((i_open >= i_ext) as u8) << 2 | ((d_open >= d_ext) as u8) << 3;
            }
        }
        last_col[i] = best_c[m];
        std::mem::swap(&mut best_p, &mut best_c);
        std::mem::swap(&mut del_p, &mut del_c);
    }
    // The SIMD backend's choice among tied end cells: it scans the padded matrix in blocks of
    // columns (the first `bs` wide, then 8 at a time), keeps the first block reaching the best
    // score, and within it the last column. An end cell's column there is the read and germline
    // consumed: padding, then the read (last row) or the read plus the germline tail laid on the
    // padding (last column), then the sentinels.
    let (k, bs) = sim_dims(n, s);
    let pos = |(i, j): (usize, usize)| {
        let c = n + j + (n - i) + k;
        (if c < bs { 0 } else { 1 + (c - bs) / 8 }, c)
    };
    let cells = (0..=m).map(|j| (best_p[j], n, j)).chain((0..n).map(|i| (last_col[i], i, m)));
    let mut end = (0, n, 0);
    let mut key = (usize::MAX, 0);
    for (v, i, j) in cells {
        let pk = pos((i, j));
        if v > end.0 || (v == end.0 && (pk.0 < key.0 || (pk.0 == key.0 && pk.1 > key.1))) {
            end = (v, i, j);
            key = pk;
        }
    }
    end
}

pub(super) fn score(query: &[u8], germline: &[u8], s: &Scoring) -> (i32, usize, usize) {
    match dp(query, germline, s, None) {
        (sc, i, j) if sc > 0 => (sc, j, i),
        _ => (0, 0, 0),
    }
}

// ponytail: full (n+1)*(m+1) traceback bytes per call; the scalar path is only the fallback.
pub(super) fn align(query: &[u8], germline: &[u8], s: &Scoring) -> Option<Alignment> {
    let mut tr = Vec::new();
    let (score, ei, ej) = dp(query, germline, s, Some(&mut tr));
    if score <= 0 {
        return None;
    }
    let m1 = query.len() + 1;
    // Traceback with the SIMD backend's preferences: Ins, then Del, then the diagonal; a gap
    // opens as soon as it can.
    let (mut i, mut j, mut state) = (ei, ej, 0u8); // 0 best, 1 Del, 2 Ins
    let mut rev = Vec::new();
    while i > 0 && j > 0 {
        let t = tr[i * m1 + j];
        match state {
            0 if t & 1 != 0 => state = 2,
            0 if t & 2 != 0 => state = 1,
            0 => {
                rev.push(Op::Match);
                i -= 1;
                j -= 1;
            }
            1 => {
                rev.push(Op::Del);
                state = if t & 8 != 0 { 0 } else { 1 };
                i -= 1;
            }
            _ => {
                rev.push(Op::Ins);
                state = if t & 4 != 0 { 0 } else { 2 };
                j -= 1;
            }
        }
    }
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
                push(&mut ops, if is_match(query[j], germline[i]) { Op::Match } else { Op::Mismatch });
                end = (j + 1, i + 1);
                keep = ops.len();
                i += 1;
                j += 1;
            }
        }
    }
    let (q_start, g_start) = start?;
    ops.truncate(keep);
    Some(Alignment { score, q_start, q_end: end.0, g_start, g_end: end.1, ops })
}
