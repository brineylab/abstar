use super::{Alignment, Op, Scoring};
use std::ops::Range;

pub(super) fn norm(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        c @ (b'A' | b'C' | b'G' | b'T') => c,
        _ => b'N',
    }
}

/// A non-ACGT query base never matches.
pub(super) fn is_match(q: u8, g: u8) -> bool {
    let q = q.to_ascii_uppercase();
    matches!(q, b'A' | b'C' | b'G' | b'T') && q == g.to_ascii_uppercase()
}

pub(super) fn push(ops: &mut Vec<(Op, u32)>, op: Op) {
    match ops.last_mut() {
        Some((last, n)) if *last == op => *n += 1,
        _ => ops.push((op, 1)),
    }
}

pub(super) fn score_ops(ops: &[(Op, u32)], s: &Scoring) -> i32 {
    ops.iter()
        .map(|&(op, n)| {
            let n = n as i32;
            match op {
                Op::Match => n * s.matched as i32,
                Op::Mismatch => n * s.mismatch as i32,
                Op::Ins | Op::Del => s.gap_open as i32 + (n - 1) * s.gap_extend as i32,
            }
        })
        .sum()
}

/// Sentinel bonus of the SIMD backend (see `simd.rs`).
pub(super) const BONUS: i8 = 127;

/// The SIMD backend's sentinel count and block size for a germline of length `gl`. K sentinels
/// outscore the costliest way to reach the last row from any cell (deleting the rest of the
/// germline). A block taller than the germline covers every row, so the DP is exact; the size
/// depends only on the germline length, so results never depend on buffer history. The scalar
/// backend uses these to break ties as the SIMD backend does.
pub(super) fn sim_dims(gl: usize, s: &Scoring) -> (usize, usize) {
    let worst = -(s.gap_open as i64) - (s.gap_extend as i64) * gl as i64;
    let k = (worst.max(0) / BONUS as i64 + 1) as usize;
    (k, (gl + k + 1).next_power_of_two().max(32))
}

/// SIMD backend available on this CPU (std caches the CPUID probe).
#[cfg(target_arch = "x86_64")]
fn simd() -> bool {
    std::arch::is_x86_feature_detected!("avx2")
}

#[cfg(target_arch = "aarch64")]
fn simd() -> bool {
    true
}

/// Overlap alignment before gap placement: `score` is the DP optimum, ops exclude terminal gaps.
/// `None` when the optimum is <= 0.
fn raw_align(query: &[u8], germline: &[u8], s: &Scoring) -> Option<Alignment> {
    #[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
    if simd() {
        return super::simd::align(query, germline, s);
    }
    super::scalar::align(query, germline, s)
}

/// Semi-global (overlap) affine-gap alignment: the query's (read's) ends are free and the
/// germline's 5' and 3' ends are free; the alignment starts on the first base of one sequence and
/// ends on the last base of one. Uses AVX2 (x86_64, detected at run time) or NEON (aarch64), and
/// a scalar Gotoh aligner with the same semantics elsewhere.
///
/// Terminal gaps: the alignment is chosen with terminal gaps penalized like any gap, then they are
/// dropped from `ops` and the coordinates, so it starts and ends on aligned bases. `score` is the
/// DP optimum, so it is the ops' score plus the dropped gaps' cost (it equals the ops' score
/// unless a terminal gap was dropped) and equals `score_semiglobal`'s score.
///
/// When the optimum is <= 0 (nothing aligns usefully) the result is empty: score 0, all
/// coordinates 0, no ops. Gaps are placed by `place_gaps` (no protected range), so gap placement
/// never depends on the backend.
pub fn align_semiglobal(query: &[u8], germline: &[u8], s: &Scoring) -> Alignment {
    match raw_align(query, germline, s) {
        Some(aln) => place_gaps(&aln, query, germline, None, s),
        None => Alignment { score: 0, q_start: 0, q_end: 0, g_start: 0, g_end: 0, ops: vec![] },
    }
}

/// Score-only `align_semiglobal`: (score, q_end, g_end), without traceback. The score is
/// `align_semiglobal`'s; the ends are the DP's end cell, equal to `align_semiglobal`'s unless it
/// dropped a trailing gap. (0, 0, 0) when the optimum is <= 0.
pub fn score_semiglobal(query: &[u8], germline: &[u8], s: &Scoring) -> (i32, usize, usize) {
    #[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
    if simd() {
        return super::simd::score(query, germline, s);
    }
    super::scalar::score(query, germline, s)
}

/// Moves each gap (left to right, the others fixed) to the placement of its kind and length with
/// the fewest query-vs-germline mismatches over the aligned span, leftmost on ties. A placement
/// needs an aligned base on each side and may not overlap another gap of its kind; its flanks are
/// the germline positions of those two bases, and it is rejected if `[left flank, right flank]`
/// intersects `protect` (germline coordinates).
/// A gap with no allowed placement stays where it is. The span is unchanged and the ops are
/// rebuilt; the score changes by the change in the ops' score under `s`, so an offset already in
/// `aln.score` (a dropped terminal gap's cost, see `align_semiglobal`) is kept.
///
/// Mismatches count every aligned pair that is not a match, so a read N against a germline N is a
/// mismatch (the generator counts equal bytes as a match; this differs only if a germline has N).
pub fn place_gaps(aln: &Alignment, query: &[u8], germline: &[u8], protect: Option<Range<usize>>, s: &Scoring) -> Alignment {
    // (is insertion, start, length): an insertion starts at a query position, a deletion at a
    // germline position.
    let mut gaps: Vec<(bool, usize, usize)> = Vec::new();
    let (mut qi, mut gi) = (aln.q_start, aln.g_start);
    for &(op, n) in &aln.ops {
        let n = n as usize;
        match op {
            Op::Ins => {
                gaps.push((true, qi, n));
                qi += n;
            }
            Op::Del => {
                gaps.push((false, gi, n));
                gi += n;
            }
            Op::Match | Op::Mismatch => {
                qi += n;
                gi += n;
            }
        }
    }
    let (gs, gn, qs, qn) = (aln.g_start, aln.g_end - aln.g_start, aln.q_start, aln.q_end - aln.q_start);
    // at[u - gs] = query position aligned to germline u, walking from the left (or right) end of
    // the span with every gap but `skip`.
    let walk = |gaps: &[(bool, usize, usize)], skip: usize, from_right: bool, at: &mut Vec<Option<usize>>| {
        let covers = |ins: bool, x: usize| gaps.iter().enumerate().any(|(k, &(i, p, n))| k != skip && i == ins && (p..p + n).contains(&x));
        at.clear();
        at.resize(gn, None);
        let mut q = 0; // query bases consumed from the walk's end
        for k in 0..gn {
            let u = if from_right { gn - 1 - k } else { k };
            if covers(false, gs + u) {
                continue;
            }
            let qpos = |q: usize| if from_right { qs + qn - 1 - q } else { qs + q };
            while q < qn && covers(true, qpos(q)) {
                q += 1;
            }
            if q < qn {
                at[u] = Some(qpos(q));
            }
            q += 1;
        }
    };
    let mm = |at: Option<usize>, u: usize| at.is_some_and(|q| !is_match(query[q], germline[gs + u])) as u32;
    let (mut left, mut right) = (Vec::new(), Vec::new());
    let (mut inv_l, mut inv_r) = (vec![None; qn], vec![None; qn]);
    let (mut pl, mut sr) = (vec![0u32; gn + 1], vec![0u32; gn + 1]);
    for j in 0..gaps.len() {
        let (ins, p0, len) = gaps[j];
        // With gap j at a placement whose aligned flanks are germline l and r, the columns left
        // of it follow the left walk and those right of it the right walk, so its mismatches are
        // pl[l + 1] + sr[r] (prefix and suffix sums).
        walk(&gaps, j, false, &mut left);
        walk(&gaps, j, true, &mut right);
        inv_l.fill(None);
        inv_r.fill(None);
        for u in 0..gn {
            pl[u + 1] = pl[u] + mm(left[u], u);
            if let Some(q) = left[u] {
                inv_l[q - qs] = Some(gs + u);
            }
            if let Some(q) = right[u] {
                inv_r[q - qs] = Some(gs + u);
            }
        }
        for u in (0..gn).rev() {
            sr[u] = sr[u + 1] + mm(right[u], u);
        }
        let (lo, hi) = if ins { (qs, qs + qn) } else { (gs, gs + gn) };
        let mut best: Option<(u32, usize)> = None;
        for p in lo + 1..hi.saturating_sub(len) {
            if gaps.iter().enumerate().any(|(k, &(i, x, n))| k != j && i == ins && x < p + len && p < x + n) {
                continue;
            }
            let flanks = if ins {
                inv_l[p - 1 - qs].zip(inv_r[p + len - qs])
            } else {
                left[p - 1 - gs].and(right[p + len - gs]).map(|_| (p - 1, p + len))
            };
            let Some((l, r)) = flanks else { continue };
            if protect.as_ref().is_some_and(|pr| pr.start < pr.end && l < pr.end && pr.start <= r) {
                continue;
            }
            let n = pl[l - gs + 1] + sr[r - gs];
            if best.is_none_or(|(b, _)| n < b) {
                best = Some((n, p));
            }
        }
        gaps[j] = (ins, best.map_or(p0, |(_, p)| p), len);
    }
    let mut at = Vec::new();
    walk(&gaps, usize::MAX, false, &mut at);
    let mut ops = Vec::new();
    let (mut q, mut u) = (aln.q_start, 0);
    while q < aln.q_end || u < at.len() {
        match at.get(u) {
            Some(None) => {
                push(&mut ops, Op::Del);
                u += 1;
            }
            Some(&Some(x)) if x == q => {
                let m = is_match(query[q], germline[aln.g_start + u]);
                push(&mut ops, if m { Op::Match } else { Op::Mismatch });
                q += 1;
                u += 1;
            }
            _ => {
                push(&mut ops, Op::Ins);
                q += 1;
            }
        }
    }
    let score = aln.score - score_ops(&aln.ops, s) + score_ops(&ops, s);
    Alignment { score, q_start: aln.q_start, q_end: aln.q_end, g_start: aln.g_start, g_end: aln.g_end, ops }
}

#[cfg(all(test, any(target_arch = "x86_64", target_arch = "aarch64")))]
mod tests {
    use super::*;

    fn rand_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    /// The scalar fallback gives the SIMD backend's alignments and scores.
    #[test]
    fn scalar_matches_simd() {
        if !simd() {
            return;
        }
        let s = Scoring::default();
        for t in 0..300u64 {
            let v = rand_seq(60 + (t as usize * 13) % 260, 5000 + t);
            let cut = (t as usize * 7) % (v.len() / 2);
            let mut read = match t % 4 {
                0 => [rand_seq(t as usize % 40, 6000 + t), v[cut..].to_vec()].concat(), // 5' lead / truncation
                1 => v[..v.len() - cut].to_vec(),                                          // 3' truncation
                2 => [&v[..v.len() / 2], &rand_seq(1 + t as usize % 30, 7000 + t)[..], &v[v.len() / 2..]].concat(),
                _ => [&v[..v.len() / 3], &v[v.len() / 3 + 1 + t as usize % 30..]].concat(),
            };
            let mut x = t;
            for b in read.iter_mut() {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                if (x >> 60) < 2 + t % 3 {
                    *b = b"ACGTN"[(x >> 40) as usize % 5];
                }
            }
            read.extend(rand_seq(t as usize % 25, 8000 + t));
            if t % 10 == 9 {
                read = rand_seq(40 + t as usize % 50, 9000 + t); // unrelated
            }
            let place = |a: Option<Alignment>| a.map(|a| place_gaps(&a, &read, &v, None, &s));
            assert_eq!(place(super::super::scalar::align(&read, &v, &s)), place(super::super::simd::align(&read, &v, &s)), "case {t}");
            assert_eq!(super::super::scalar::score(&read, &v, &s), super::super::simd::score(&read, &v, &s), "case {t}");
        }
    }
}
