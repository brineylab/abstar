use super::{Alignment, Op, Scoring};
use block_aligner::cigar::{Cigar, Operation};
use block_aligner::scan_block::{Block, PaddedBytes};
use block_aligner::scores::{Gaps, Matrix, NucMatrix};
use std::cell::RefCell;
use std::ops::Range;

// ponytail: block-aligner has no scalar backend; add one (or simd_sse2/simd_wasm) for other targets.
#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
compile_error!("the aligner needs x86_64 (AVX2) or aarch64 (NEON)");

// Overlap alignment on block-aligner. It frees only the start of its reference (the read):
// FREE_QUERY_START_GAPS with the germline as its query. Its FREE_QUERY_END_GAPS end search is
// unusable in 0.5.1: the last-row maximum is a per-lane maximum over every row congruent to the
// query length mod the SIMD width, so interior rows leak into the score and the end column is
// wrong. Instead:
// - the read is padded on both sides with PAD, which scores 0 against every base, so germline
//   bases overhanging either read end are free;
// - the germline gets K trailing SENT bases scoring +127 against everything, enough that the
//   best cell of the whole matrix is on the last row, and X-drop mode (with no drop) reports
//   that best cell exactly. The bonus is the same for every complete path, so it never changes
//   which path wins.
// PAD (`E`) and SENT (`H`) collide with no ACGTN or block-aligner's `Z` in NucMatrix's index bits.
const PAD: u8 = b'E';
const SENT: u8 = b'H';
const BONUS: i8 = 127;

fn norm(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        c @ (b'A' | b'C' | b'G' | b'T') => c,
        _ => b'N',
    }
}

/// A non-ACGT query base never matches.
fn is_match(q: u8, g: u8) -> bool {
    let q = q.to_ascii_uppercase();
    matches!(q, b'A' | b'C' | b'G' | b'T') && q == g.to_ascii_uppercase()
}

fn push(ops: &mut Vec<(Op, u32)>, op: Op) {
    match ops.last_mut() {
        Some((last, n)) if *last == op => *n += 1,
        _ => ops.push((op, 1)),
    }
}

fn score_ops(ops: &[(Op, u32)], s: &Scoring) -> i32 {
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

/// Per-thread buffers, grown on demand and reused across calls.
struct Scratch {
    block: Block<true, true, false, true, false>,
    cap: (usize, usize, usize), // germline + sentinel length, padded read length, block size
    g: PaddedBytes,
    r: PaddedBytes,
    gbuf: Vec<u8>,
    rbuf: Vec<u8>,
    cigar: Cigar,
}

impl Scratch {
    fn new((g, r, b): (usize, usize, usize)) -> Scratch {
        Scratch {
            block: Block::new(g, r, b),
            cap: (g, r, b),
            g: PaddedBytes::new::<NucMatrix>(g, b),
            r: PaddedBytes::new::<NucMatrix>(r, b),
            gbuf: Vec::with_capacity(g),
            rbuf: Vec::with_capacity(r),
            cigar: Cigar::new(g, r),
        }
    }
}

thread_local! {
    static SCRATCH: RefCell<Option<Scratch>> = const { RefCell::new(None) };
}

/// Semi-global affine-gap alignment: query (read) ends are free, and germline 5' and 3' ends
/// are free (overlap alignment: it starts on the first base of one sequence and ends on the
/// last base of one). Terminal gaps are dropped, so the alignment starts and ends on aligned
/// bases. Gaps are then placed by `place_gaps` (no protected range), so ties never depend on the
/// SIMD backend.
pub fn align_semiglobal(query: &[u8], germline: &[u8], s: &Scoring) -> Alignment {
    let empty = Alignment { score: 0, q_start: 0, q_end: 0, g_start: 0, g_end: 0, ops: vec![] };
    if query.is_empty() || germline.is_empty() {
        return empty;
    }
    let (gl, ql) = (germline.len(), query.len());
    // K sentinels outscore the costliest way to reach the last row from any cell (deleting the
    // rest of the germline), so the best cell is on it.
    let worst = -(s.gap_open as i64) - (s.gap_extend as i64) * gl as i64;
    let k = (worst.max(0) / BONUS as i64 + 1) as usize;
    // ponytail: pads by a whole germline on each side (~2x the DP width, ~140 us per 300-nt V on
    // AVX2); pad only the overhang a caller allows if alignment time matters.
    let (hl, rl) = (gl + k, ql + 2 * gl + k);
    // A block taller than the germline covers every row, so the DP is exact. The size depends
    // only on the germline length, so results never depend on buffer history.
    let bs = (hl + 1).next_power_of_two().max(32);
    let mut matrix = NucMatrix::new_simple(s.matched, s.mismatch);
    for c in [b'A', b'C', b'G', b'T', b'N', PAD] {
        matrix.set(PAD, c, 0);
        matrix.set(SENT, c, BONUS);
    }
    matrix.set(b'N', b'N', s.mismatch);
    let gaps = Gaps { open: s.gap_open, extend: s.gap_extend };

    let mut ops = Vec::new();
    let (mut start, mut end, mut keep) = (None, (0, 0), 0);
    SCRATCH.with_borrow_mut(|slot| {
        let sc = match slot {
            Some(sc) if sc.cap.0 >= hl && sc.cap.1 >= rl && sc.cap.2 >= bs => sc,
            _ => {
                let cap = slot.as_ref().map_or((hl, rl, bs), |o| (o.cap.0.max(hl), o.cap.1.max(rl), o.cap.2.max(bs)));
                slot.insert(Scratch::new(cap))
            }
        };
        sc.gbuf.clear();
        sc.gbuf.extend(germline.iter().map(|&b| norm(b)));
        sc.gbuf.resize(hl, SENT);
        sc.rbuf.clear();
        sc.rbuf.resize(gl, PAD);
        sc.rbuf.extend(query.iter().map(|&b| norm(b)));
        sc.rbuf.resize(rl, PAD);
        sc.g.set_bytes::<NucMatrix>(&sc.gbuf, bs);
        sc.r.set_bytes::<NucMatrix>(&sc.rbuf, bs);
        sc.block.align(&sc.g, &sc.r, &matrix, gaps, bs..=bs, i32::MAX);
        let res = sc.block.res();
        debug_assert_eq!(res.query_idx, hl);
        sc.block.trace().cigar(res.query_idx, res.reference_idx, &mut sc.cigar);

        // block-aligner's I = germline only (our Del), D = read only (our Ins).
        let (mut i, mut j) = (res.query_idx, res.reference_idx);
        for k in 0..sc.cigar.len() {
            let o = sc.cigar.get(k);
            if matches!(o.op, Operation::M | Operation::Eq | Operation::X | Operation::I) {
                i -= o.len;
            }
            if matches!(o.op, Operation::M | Operation::Eq | Operation::X | Operation::D) {
                j -= o.len;
            }
        }
        let read = gl..gl + ql;
        for k in 0..sc.cigar.len() {
            let o = sc.cigar.get(k);
            for _ in 0..o.len {
                match o.op {
                    Operation::I => {
                        if start.is_some() {
                            push(&mut ops, Op::Del);
                        }
                        i += 1;
                    }
                    Operation::D => {
                        if start.is_some() && read.contains(&j) {
                            push(&mut ops, Op::Ins);
                        }
                        j += 1;
                    }
                    _ => {
                        if i < gl && read.contains(&j) {
                            let q = j - gl;
                            start.get_or_insert((q, i));
                            push(&mut ops, if is_match(query[q], germline[i]) { Op::Match } else { Op::Mismatch });
                            end = (q + 1, i + 1);
                            keep = ops.len();
                        }
                        i += 1;
                        j += 1;
                    }
                }
            }
        }
    });
    let Some((q_start, g_start)) = start else { return empty };
    ops.truncate(keep);
    let aln = Alignment { score: 0, q_start, q_end: end.0, g_start, g_end: end.1, ops };
    place_gaps(&aln, query, germline, None, s)
}

/// Moves each gap (left to right, the others fixed) to the placement of its kind and length with
/// the fewest query-vs-germline mismatches over the aligned span, leftmost on ties. A placement
/// needs an aligned base on each side; its flanks are the germline positions of those two bases,
/// and it is rejected if `[left flank, right flank]` intersects `protect` (germline coordinates).
/// A gap with no allowed placement stays where it is. The span is unchanged; ops and score are
/// rebuilt (the score with `s`).
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
    Alignment { score: score_ops(&ops, s), q_start: aln.q_start, q_end: aln.q_end, g_start: aln.g_start, g_end: aln.g_end, ops }
}
