//! Throwaway V-assignment throughput spike (Task 1). Not production code.
//!
//! usage: vsearch <v.fasta> <seqs.fasta[.gz]> <k> <n_candidates> <reps> <threads,csv> [preds.tsv]
//!
//! Per record: k-mer prefilter over both orientations -> top-N (allele, strand)
//! candidates -> score-only SIMD semi-global affine alignment of each candidate
//! (germline end-to-end, read ends free) -> best score wins -> one traceback
//! alignment (CIGAR) for the winner, as a real engine would need.

use block_aligner::{cigar::Cigar, scan_block::*, scores::*};
use rayon::prelude::*;
use std::io::{BufRead, BufReader, Read, Write};
use std::time::Instant;

const BLOCK: usize = 512; // must exceed germline length (max 317) for FREE_QUERY_END_GAPS
const GAPS: Gaps = Gaps { open: -12, extend: -2 };

fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let f = std::fs::File::open(path).expect("open fasta");
    let r: Box<dyn Read> = if path.ends_with(".gz") { Box::new(flate2::read::MultiGzDecoder::new(f)) } else { Box::new(f) };
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in BufReader::new(r).lines() {
        let line = line.unwrap();
        if let Some(name) = line.strip_prefix('>') {
            out.push((name.to_string(), Vec::new()));
        } else {
            // Uppercase; anything that is not ACGT becomes N (block-aligner's NucMatrix knows only ATCGN).
            out.last_mut().unwrap().1.extend(line.trim().bytes().map(|b| match b.to_ascii_uppercase() {
                c @ (b'A' | b'C' | b'G' | b'T') => c,
                _ => b'N',
            }));
        }
    }
    out
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| match b { b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', _ => b'N' }).collect()
}

fn code(b: u8) -> Option<u32> {
    match b { b'A' => Some(0), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None }
}

/// Calls `f(kmer)` for every N-free k-mer of `s`.
fn for_kmers(s: &[u8], k: usize, mut f: impl FnMut(usize)) {
    let mask = (1u32 << (2 * k)) - 1;
    let (mut km, mut valid) = (0u32, 0usize);
    for &b in s {
        match code(b) {
            Some(c) => { km = ((km << 2) | c) & mask; valid += 1; }
            None => valid = 0,
        }
        if valid >= k { f(km as usize); }
    }
}

/// CSR k-mer -> allele index (each allele listed once per k-mer).
struct Index { k: usize, offsets: Vec<u32>, postings: Vec<u16> }

impl Index {
    fn new(germ: &[(String, Vec<u8>)], k: usize) -> Self {
        let mut pairs: Vec<(u32, u16)> = Vec::new();
        for (a, (_, s)) in germ.iter().enumerate() {
            for_kmers(s, k, |km| pairs.push((km as u32, a as u16)));
        }
        pairs.sort_unstable();
        pairs.dedup();
        let mut offsets = vec![0u32; (1 << (2 * k)) + 1];
        for &(km, _) in &pairs { offsets[km as usize + 1] += 1; }
        for i in 1..offsets.len() { offsets[i] += offsets[i - 1]; }
        Index { k, offsets, postings: pairs.into_iter().map(|p| p.1).collect() }
    }
}

struct Ctx<'a> {
    germ: &'a [(String, Vec<u8>)],
    germ_pb: Vec<PaddedBytes>,
    index: Index,
    n_cand: usize,
    matrix: NucMatrix,
}

struct Scratch {
    score_block: Block<false, false, false, true, true>,
    trace_block: Block<true, false, false, true, true>,
    read_pb: [PaddedBytes; 2],
    counts: Vec<u16>,
    cand: Vec<(u16, u32)>,
    cigar: Cigar,
}

struct Call { best: Vec<u16>, strand: u8, score: i32, cigar_ops: usize }

impl Ctx<'_> {
    fn scratch(&self) -> Scratch {
        let ml = 2048;
        Scratch {
            score_block: Block::new(BLOCK, ml, BLOCK),
            trace_block: Block::new(BLOCK, ml, BLOCK),
            read_pb: [PaddedBytes::new::<NucMatrix>(ml, BLOCK), PaddedBytes::new::<NucMatrix>(ml, BLOCK)],
            counts: vec![0; 2 * self.germ.len()],
            cand: Vec::new(),
            cigar: Cigar::new(BLOCK, ml),
        }
    }

    fn assign(&self, read: &[u8], sc: &mut Scratch) -> Call {
        let n = self.germ.len();
        let rc = revcomp(read);
        let strands = [read, &rc[..]];
        // 1. k-mer prefilter, both strands.
        sc.counts.fill(0);
        for (st, s) in strands.iter().enumerate() {
            let counts = &mut sc.counts[st * n..(st + 1) * n];
            for_kmers(s, self.index.k, |km| {
                let (a, b) = (self.index.offsets[km] as usize, self.index.offsets[km + 1] as usize);
                for &al in &self.index.postings[a..b] { counts[al as usize] += 1; }
            });
        }
        sc.cand.clear();
        sc.cand.extend(sc.counts.iter().enumerate().filter(|(_, c)| **c > 0).map(|(i, &c)| (c, i as u32)));
        let take = self.n_cand.min(sc.cand.len());
        if take < sc.cand.len() {
            sc.cand.select_nth_unstable_by(take, |a, b| b.0.cmp(&a.0));
            sc.cand.truncate(take);
        }
        // 2. Score-only SIMD alignment of each candidate. Query = germline (end to end),
        //    reference = read (free leading/trailing read sequence).
        let mut loaded = [false; 2];
        let (mut best, mut best_score, mut best_strand) = (Vec::new(), i32::MIN, 0u8);
        for &(_, ci) in &sc.cand {
            let (st, al) = (ci as usize / n, ci as usize % n);
            if !loaded[st] { sc.read_pb[st].set_bytes::<NucMatrix>(strands[st], BLOCK); loaded[st] = true; }
            sc.score_block.align(&self.germ_pb[al], &sc.read_pb[st], &self.matrix, GAPS, BLOCK..=BLOCK, 0);
            let s = sc.score_block.res().score;
            if s > best_score { best_score = s; best.clear(); best_strand = st as u8; }
            if s == best_score && st as u8 == best_strand { best.push(al as u16); }
        }
        // 3. Traceback for the winner (needed downstream for V coordinates).
        let mut cigar_ops = 0;
        if let Some(&al) = best.first() {
            let st = best_strand as usize;
            sc.trace_block.align(&self.germ_pb[al as usize], &sc.read_pb[st], &self.matrix, GAPS, BLOCK..=BLOCK, 0);
            let r = sc.trace_block.res();
            sc.trace_block.trace().cigar(r.query_idx, r.reference_idx, &mut sc.cigar);
            cigar_ops = sc.cigar.len();
        }
        Call { best, strand: best_strand, score: best_score, cigar_ops }
    }
}

fn gene(name: &str) -> &str {
    name.split("__").next().unwrap().split('*').next().unwrap()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 7 {
        eprintln!("usage: vsearch <v.fasta> <seqs.fasta[.gz]> <k> <n_candidates> <reps> <threads,csv> [preds.tsv]");
        std::process::exit(2);
    }
    let germ = read_fasta(&a[1]);
    let reads: Vec<Vec<u8>> = read_fasta(&a[2]).into_iter().map(|r| r.1).collect();
    let (k, n_cand, reps): (usize, usize, usize) = (a[3].parse().unwrap(), a[4].parse().unwrap(), a[5].parse().unwrap());
    let threads: Vec<usize> = a[6].split(',').map(|t| t.parse().unwrap()).collect();
    assert!(k <= 13 && germ.len() < u16::MAX as usize && germ.iter().all(|g| g.1.len() < BLOCK));
    assert!(reads.iter().all(|r| r.len() < 2048));

    let mut matrix = NucMatrix::new_simple(2, -3);
    for c in [b'A', b'C', b'G', b'T', b'N'] { matrix.set(b'N', c, 0); }
    let ctx = Ctx {
        germ: &germ,
        germ_pb: germ.iter().map(|g| PaddedBytes::from_bytes::<NucMatrix>(&g.1, BLOCK)).collect(),
        index: Index::new(&germ, k),
        n_cand,
        matrix,
    };
    let total = reads.len() * reps;
    eprintln!("alleles={} reads={} k={} candidates={} reps={} total={}", germ.len(), reads.len(), k, n_cand, reps, total);

    // Optional prediction dump (one pass over the corpus, not timed).
    if let Some(out) = a.get(7) {
        let calls: Vec<Call> = reads.par_iter().map_init(|| ctx.scratch(), |sc, r| ctx.assign(r, sc)).collect();
        let mut w = std::io::BufWriter::new(std::fs::File::create(out).unwrap());
        writeln!(w, "corpus_ordinal\tstrand\tscore\tbest_alleles\tbest_genes").unwrap();
        for (i, c) in calls.iter().enumerate() {
            let alleles: Vec<&str> = c.best.iter().map(|&b| germ[b as usize].0.split("__").next().unwrap()).collect();
            let mut genes: Vec<&str> = c.best.iter().map(|&b| gene(&germ[b as usize].0)).collect();
            genes.dedup();
            writeln!(w, "{i}\t{}\t{}\t{}\t{}", if c.strand == 0 { "+" } else { "-" }, c.score, alleles.join(","), genes.join(",")).unwrap();
        }
    }

    // Timed runs. Replication = cycling the in-memory corpus `reps` times (index modulo corpus size).
    for &t in &threads {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(t).build().unwrap();
        let start = Instant::now();
        let checksum: i64 = pool.install(|| {
            (0..total).into_par_iter()
                .map_init(|| ctx.scratch(), |sc, i| { let c = ctx.assign(&reads[i % reads.len()], sc); c.score as i64 + c.cigar_ops as i64 })
                .sum()
        });
        let secs = start.elapsed().as_secs_f64();
        println!("threads={t}\trecords={total}\tseconds={secs:.2}\trecords_per_sec={:.0}\tper_thread={:.0}\tchecksum={checksum}",
            total as f64 / secs, total as f64 / secs / t as f64);
    }
}
