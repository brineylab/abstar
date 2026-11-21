//! block-aligner backend (AVX2 on x86_64, NEON on aarch64).
//!
//! Overlap alignment on block-aligner. It frees only the start of its reference (the read):
//! FREE_QUERY_START_GAPS with the germline as its query. Its FREE_QUERY_END_GAPS end search is
//! unusable in 0.5.1: the last-row maximum is a per-lane maximum over every row congruent to the
//! query length mod the SIMD width, so interior rows leak into the score and the end column is
//! wrong. Instead:
//! - the read is padded on both sides with PAD, which scores 0 against every base, so germline
//!   bases overhanging either read end are free;
//! - the germline gets K trailing SENT bases scoring +127 against everything, enough that the
//!   best cell of the whole matrix is on the last row, and X-drop mode (with no drop) reports
//!   that best cell exactly. The bonus is the same for every complete path, so it never changes
//!   which path wins.
//!
//! PAD (`E`) and SENT (`H`) collide with no ACGTN or block-aligner's `Z` in NucMatrix's index bits.
use super::aligner::{is_match, norm, push, sim_dims, BONUS};
use super::{Alignment, Op, Scoring};
use block_aligner::cigar::{Cigar, Operation};
use block_aligner::scan_block::{AlignResult, Block, PaddedBytes};
use block_aligner::scores::{Gaps, Matrix, NucMatrix};
use std::cell::RefCell;

const PAD: u8 = b'E';
const SENT: u8 = b'H';

/// Per-thread buffers, grown on demand and reused across calls.
struct Scratch {
    trace: Block<true, true, false, true, false>,
    score: Block<false, true, false, true, false>,
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
            trace: Block::new(g, r, b),
            score: Block::new(g, r, b),
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

/// Germline length, read length, sentinel count.
struct Dims {
    gl: usize,
    ql: usize,
    k: usize,
}

/// Runs block-aligner on the padded problem (TRACE picks the block) and hands `f` the scratch and
/// the result. Returns `None` for empty input.
fn run<T>(query: &[u8], germline: &[u8], s: &Scoring, trace: bool, f: impl FnOnce(&mut Scratch, &Dims, AlignResult) -> T) -> Option<T> {
    if query.is_empty() || germline.is_empty() {
        return None;
    }
    let (gl, ql) = (germline.len(), query.len());
    let (k, bs) = sim_dims(gl, s);
    // ponytail: pads by a whole germline on each side (~2x the DP width); the overhang can be up
    // to a germline long, so narrowing needs a caller-supplied bound.
    let (hl, rl) = (gl + k, ql + 2 * gl + k);
    let mut matrix = NucMatrix::new_simple(s.matched, s.mismatch);
    for c in [b'A', b'C', b'G', b'T', b'N', PAD] {
        matrix.set(PAD, c, 0);
        matrix.set(SENT, c, BONUS);
    }
    matrix.set(b'N', b'N', s.mismatch);
    let gaps = Gaps { open: s.gap_open, extend: s.gap_extend };
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
        let res = if trace {
            sc.trace.align(&sc.g, &sc.r, &matrix, gaps, bs..=bs, i32::MAX);
            sc.trace.res()
        } else {
            sc.score.align(&sc.g, &sc.r, &matrix, gaps, bs..=bs, i32::MAX);
            sc.score.res()
        };
        debug_assert_eq!(res.query_idx, hl);
        let res = AlignResult { score: res.score - k as i32 * BONUS as i32, ..res };
        Some(f(sc, &Dims { gl, ql, k }, res))
    })
}

pub(super) fn score(query: &[u8], germline: &[u8], s: &Scoring) -> (i32, usize, usize) {
    run(query, germline, s, false, |_, d, res| {
        if res.score <= 0 {
            return (0, 0, 0);
        }
        // the germline's last row is reached K columns before the sentinel row's end; past the
        // read, the germline's tail lies on padding
        let c = res.reference_idx - d.k;
        if c <= d.gl + d.ql {
            (res.score, c - d.gl, d.gl)
        } else {
            (res.score, d.ql, d.gl - (c - d.gl - d.ql))
        }
    })
    .unwrap_or((0, 0, 0))
}

pub(super) fn align(query: &[u8], germline: &[u8], s: &Scoring) -> Option<Alignment> {
    run(query, germline, s, true, |sc, d, res| {
        if res.score <= 0 {
            return None;
        }
        let gl = d.gl;
        sc.trace.trace().cigar(res.query_idx, res.reference_idx, &mut sc.cigar);
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
        let read = gl..gl + d.ql;
        let mut ops = Vec::new();
        let (mut start, mut end, mut keep) = (None, (0, 0), 0);
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
        let (q_start, g_start) = start?;
        ops.truncate(keep);
        Some(Alignment { score: res.score, q_start, q_end: end.0, g_start, g_end: end.1, ops })
    })
    .flatten()
}
