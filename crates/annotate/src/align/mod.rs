//! Search and alignment primitives: a k-mer prefilter index and a SIMD affine-gap
//! semi-global aligner with deterministic gap placement.
//!
//! Coordinates are 0-based half-open. `query` is the (oriented) read, `germline` the allele.
mod aligner;
mod kmer;

pub use aligner::{align_semiglobal, place_gaps};
pub use kmer::KmerIndex;

/// One alignment column type. `Ins` = base in the query only; `Del` = base in the germline only.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    Match,
    Mismatch,
    Ins,
    Del,
}

/// A query-vs-germline alignment: `ops` (run-length, adjacent runs of one kind merged) cover
/// `query[q_start..q_end]` and `germline[g_start..g_end]`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignment {
    pub score: i32,
    pub q_start: usize,
    pub q_end: usize,
    pub g_start: usize,
    pub g_end: usize,
    pub ops: Vec<(Op, u32)>,
}

/// Linear scores; a gap of length n costs `gap_open + (n - 1) * gap_extend`.
/// A non-ACGT query base never scores as a match (it scores `mismatch`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Scoring {
    pub matched: i8,
    pub mismatch: i8,
    pub gap_open: i8,
    pub gap_extend: i8,
}

impl Default for Scoring {
    /// Record 0029.
    fn default() -> Self {
        Scoring { matched: 2, mismatch: -3, gap_open: -10, gap_extend: -1 }
    }
}

/// AIRR CIGAR: leading `S` (query before the alignment), leading `N` (germline before it),
/// `M` (match or mismatch), `I`, `D`, trailing `S` (query after it).
pub fn cigar(a: &Alignment, query_len: usize) -> String {
    let mut out = String::new();
    let mut push = |n: usize, c: char| {
        if n > 0 {
            out.push_str(&format!("{n}{c}"));
        }
    };
    push(a.q_start, 'S');
    push(a.g_start, 'N');
    let mut i = 0;
    while i < a.ops.len() {
        let c = match a.ops[i].0 {
            Op::Match | Op::Mismatch => 'M',
            Op::Ins => 'I',
            Op::Del => 'D',
        };
        let mut n = 0;
        while i < a.ops.len() && matches!((c, a.ops[i].0), ('M', Op::Match | Op::Mismatch) | ('I', Op::Ins) | ('D', Op::Del)) {
            n += a.ops[i].1 as usize;
            i += 1;
        }
        push(n, c);
    }
    push(query_len - a.q_end, 'S');
    out
}
