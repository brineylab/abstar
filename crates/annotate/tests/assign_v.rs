//! Orientation and V assignment (record 0029) on the built-in human database.
mod common;
use abstar_annotate::align::{Op, Scoring};
use abstar_annotate::assign::{assign_v, orient, primary_by_tiebreak, Oriented, SegmentIndex, VHit};
use abstar_generate::Spec;
use abstar_germline::{Allele, Functionality, Locus, Segment};
use common::{generated, human_db, rand_seq, revcomp, Generated};
use rayon::prelude::*;
use std::collections::BTreeSet;
use std::sync::OnceLock;
use std::time::Instant;

struct Ref {
    alleles: Vec<Allele>,
    v: SegmentIndex,
    j: SegmentIndex,
}

impl Ref {
    fn new(alleles: Vec<Allele>) -> Ref {
        let v = SegmentIndex::new(&alleles, Segment::V);
        let j = SegmentIndex::new(&alleles, Segment::J);
        Ref { alleles, v, j }
    }
    fn call(&self, read: &[u8]) -> Option<(Oriented, VHit)> {
        let o = orient(read, &self.v.index, &self.j.index)?;
        let h = assign_v(&o, &self.alleles, &self.v, &Scoring::default())?;
        Some((o, h))
    }
    fn idx(&self, name: &str) -> usize {
        self.alleles.iter().position(|a| a.name == name).unwrap_or_else(|| panic!("no allele {name}"))
    }
    fn names(&self, ids: &[usize]) -> BTreeSet<String> {
        ids.iter().map(|&i| self.alleles[i].name.clone()).collect()
    }
}

fn human() -> &'static Ref {
    static R: OnceLock<Ref> = OnceLock::new();
    R.get_or_init(|| Ref::new(human_db().alleles))
}

fn ungapped(a: &Allele) -> Vec<u8> {
    a.sequence().into_bytes()
}

/// Same-locus V alleles with no mismatch against `a` over `a`'s bases at gapped positions
/// `..upto` (a gap or a missing position where `a` has a base is a mismatch).
fn zero_mismatch(r: &Ref, a: usize, upto: usize) -> Vec<usize> {
    let ag = r.alleles[a].sequence_gapped.as_bytes();
    (0..r.alleles.len())
        .filter(|&c| {
            let ca = &r.alleles[c];
            ca.segment == Segment::V
                && ca.locus == r.alleles[a].locus
                && (0..upto.min(ag.len())).all(|g| ag[g] == b'.' || ca.sequence_gapped.as_bytes().get(g) == Some(&ag[g]))
        })
        .collect()
}

const NP_J: &[u8] = b"GCGAAAGATCGGGGATATAGTAGTGGCTGGTAC";

#[test]
fn exact_v_read_called() {
    let r = human();
    let a = r.idx("IGHV3-23*01");
    let v = ungapped(&r.alleles[a]);
    let read = [&v[..], NP_J, &ungapped(&r.alleles[r.idx("IGHJ4*02")])].concat();
    let (o, h) = r.call(&read).expect("V hit");
    assert!(!o.rev_comp);
    assert_eq!(o.seq, read);
    let expected = zero_mismatch(r, a, usize::MAX);
    assert_eq!(r.names(&h.ties), r.names(&expected));
    assert_eq!(r.alleles[h.primary].name, "IGHV3-23*01");
    assert_eq!(h.ties[0], h.primary);
    assert_eq!(h.mismatches, 0);
    assert_eq!((h.aln.q_start, h.aln.q_end, h.aln.g_start, h.aln.g_end), (0, v.len(), 0, v.len()));
    assert_eq!(h.aln.ops, vec![(Op::Match, v.len() as u32)]);
    // Read position -> IMGT gapped index of the called allele; None past the V.
    let gapped: Vec<usize> =
        r.alleles[a].sequence_gapped.bytes().enumerate().filter(|&(_, b)| b != b'.').map(|(g, _)| g).collect();
    assert_eq!(h.imgt.len(), read.len());
    for (i, g) in h.imgt.iter().enumerate() {
        assert_eq!(*g, gapped.get(i).copied(), "read position {i}");
    }
}

#[test]
fn sibling_tie_reported() {
    let r = human();
    // First pair (database order) of same-gene IGH V alleles identical over at least 150 gapped
    // positions; the read is the first allele's bases before their first difference.
    let vs: Vec<usize> = (0..r.alleles.len()).filter(|&i| r.alleles[i].segment == Segment::V && r.alleles[i].locus == Locus::IGH).collect();
    let (a, b, p) = vs
        .iter()
        .flat_map(|&a| vs.iter().map(move |&b| (a, b)))
        .filter(|&(a, b)| a != b && r.alleles[a].gene == r.alleles[b].gene)
        .find_map(|(a, b)| {
            let (x, y) = (r.alleles[a].sequence_gapped.as_bytes(), r.alleles[b].sequence_gapped.as_bytes());
            let p = (0..x.len().min(y.len())).find(|&g| x[g] != y[g])?;
            (p >= 150 && x[..p].iter().filter(|&&c| c != b'.').count() >= 150).then_some((a, b, p))
        })
        .expect("a sibling pair");
    let read: Vec<u8> = r.alleles[a].sequence_gapped.as_bytes()[..p].iter().copied().filter(|&c| c != b'.').collect();
    let (_, h) = r.call(&read).expect("V hit");
    let expected = zero_mismatch(r, a, p);
    assert!(expected.contains(&b) && expected.len() >= 2);
    assert_eq!(r.names(&h.ties), r.names(&expected), "read from {} to gapped {p}", r.alleles[a].name);
    assert_eq!(h.primary, primary_by_tiebreak(&expected, &r.alleles));
    assert_eq!(h.mismatches, 0);
}

#[test]
fn primary_by_tiebreak_order() {
    let base = human_db().alleles.into_iter().find(|a| a.segment == Segment::V).unwrap();
    let mk = |name: &str, f: Functionality, usage: Option<f64>| Allele { name: name.into(), functionality: f, usage, ..base.clone() };
    use Functionality::*;
    let pick = |v: &[Allele]| v[primary_by_tiebreak(&(0..v.len()).collect::<Vec<_>>(), v)].name.clone();
    // Functionality: F > ORF > P > NF > unknown, before usage and name.
    let order = [Unknown, NF, P, ORF, F];
    for k in 1..order.len() {
        let v: Vec<Allele> = (0..=k).map(|i| mk(&format!("X*{:02}", i + 1), order[i], Some(1.0 - i as f64 / 10.0))).collect();
        assert_eq!(pick(&v), format!("X*{:02}", k + 1));
    }
    // Usage: higher first; a missing usage ranks after any given usage.
    assert_eq!(pick(&[mk("A*01", F, Some(0.1)), mk("B*01", F, Some(0.3))]), "B*01");
    assert_eq!(pick(&[mk("A*01", F, None), mk("B*01", F, Some(0.0))]), "B*01");
    // Natural sort of names.
    assert_eq!(pick(&[mk("IGHV1-2*10", F, None), mk("IGHV1-2*9", F, None)]), "IGHV1-2*9");
    assert_eq!(pick(&[mk("IGHV1-18*01", F, None), mk("IGHV1-2*01", F, None)]), "IGHV1-2*01");
    assert_eq!(pick(&[mk("IGHV1-2*01", F, None), mk("IGHV1-2*01", F, None)]), "IGHV1-2*01");
}

#[test]
fn reverse_complement_identical() {
    // Review Focus 1: a read and its reverse complement give the same oriented read and V hit.
    let r = human();
    let g = generated(7, "memory", 200, |p| p.indel_rate = 0.3);
    let mut called = 0;
    for i in 0..g.len() {
        let s = g.sequence(i).as_bytes();
        let fwd = r.call(s);
        let rev = r.call(&revcomp(s));
        match (fwd, rev) {
            (Some((of, hf)), Some((or, hr))) => {
                assert!(!of.rev_comp && or.rev_comp, "read {i}");
                assert_eq!(of.seq, or.seq, "read {i}");
                assert_eq!(hf, hr, "read {i}");
                called += 1;
            }
            (None, None) => {}
            _ => panic!("read {i}: one strand called, the other not"),
        }
    }
    assert_eq!(called, g.len());
}

#[test]
fn orientation_tie_goes_forward() {
    let r = human();
    let v = ungapped(&r.alleles[r.idx("IGHV3-23*01")]);
    let read = [&v[..], &revcomp(&v)].concat(); // its own reverse complement: equal evidence
    let o = orient(&read, &r.v.index, &r.j.index).unwrap();
    assert!(!o.rev_comp);
    assert_eq!(o.seq, read);
}

#[test]
fn iupac_bases_no_crash() {
    // Review Focus 2: a non-ACGT base is a mismatch for every candidate.
    let r = human();
    let a = r.idx("IGHV3-23*01");
    let clean = [&ungapped(&r.alleles[a])[..], NP_J, &ungapped(&r.alleles[r.idx("IGHJ4*02")])].concat();
    let (_, hc) = r.call(&clean).unwrap();
    let mut read = clean.to_ascii_lowercase();
    for (pos, b) in [(10, b'N'), (50, b'r'), (100, b'Y'), (150, b'k'), (200, b'w')] {
        read[pos] = b;
    }
    let (o, h) = r.call(&read).expect("V hit");
    assert_eq!(o.seq, read.to_ascii_uppercase());
    assert_eq!(h.mismatches, 5);
    assert_eq!(h.primary, hc.primary);
    assert!(r.names(&h.ties).is_superset(&r.names(&hc.ties)));
    // Reverse complement with IUPAC codes, and reads with no usable bases.
    assert_eq!(r.call(&revcomp(&read)).unwrap().1, h);
    assert!(r.call(&vec![b'N'; 400]).is_none());
    assert!(r.call(b"").is_none());
    assert!(r.call(b"RYKMSWBDHVN").is_none());
}

#[test]
fn random_dna_no_hit() {
    // Review Focus 3: random DNA, adapters and very short reads give no V hit.
    let r = human();
    let hits: Vec<usize> = (0..2000u64)
        .into_par_iter()
        .filter(|&s| r.call(&rand_seq(30 + (s as usize * 7) % 600, s + 1)).is_some())
        .map(|s| s as usize)
        .collect();
    assert!(hits.is_empty(), "random reads called: {hits:?}");
    let adapters: [&[u8]; 3] = [
        b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
        b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        b"CTGTCTCTTATACACATCTCCGAGCCCACGAGACAAGCAGAAGACGGCATACGAGAT",
    ];
    for ad in adapters {
        assert!(r.call(ad).is_none(), "{}", String::from_utf8_lossy(ad));
    }
    let v = ungapped(&r.alleles[r.idx("IGHV3-23*01")]);
    assert!(r.call(&v[100..120]).is_none(), "20 nt of V");
}

#[test]
fn tcr_read_against_ig_only_db_no_hit() {
    // Review Focus 3: every TR V allele, as a VJ read, against an IG-only database.
    let ig = Ref::new(human_db().alleles.into_iter().filter(|a| matches!(a.locus, Locus::IGH | Locus::IGK | Locus::IGL)).collect());
    let r = human();
    let trbj = ungapped(&r.alleles[r.idx("TRBJ2-7*01")]);
    let tr: Vec<usize> = (0..r.alleles.len())
        .filter(|&i| r.alleles[i].segment == Segment::V && matches!(r.alleles[i].locus, Locus::TRA | Locus::TRB | Locus::TRG | Locus::TRD))
        .collect();
    assert!(tr.len() > 100);
    let forced: Vec<String> = tr
        .par_iter()
        .filter_map(|&i| {
            let read = [&ungapped(&r.alleles[i])[..], b"GCCAGCAGTTTAGGGACAGG", &trbj].concat();
            // The full database calls it on its own locus.
            let (_, h) = r.call(&read).unwrap_or_else(|| panic!("{} not called on the full database", r.alleles[i].name));
            assert_eq!(r.alleles[h.primary].locus, r.alleles[i].locus, "{}", r.alleles[i].name);
            ig.call(&read).map(|(_, h)| format!("{} -> {}", r.alleles[i].name, ig.alleles[h.primary].name))
        })
        .collect();
    assert!(forced.is_empty(), "TR reads given an IG V: {forced:?}");
}

#[test]
fn insertion_positions_skipped() {
    // A 9-nt insertion in FWR1: inserted read bases have no IMGT position and are not counted.
    let r = human();
    let a = r.idx("IGHV3-23*01");
    let v = ungapped(&r.alleles[a]);
    let read = [&v[..60], b"TTTAAACCC", &v[60..], NP_J].concat();
    let (_, h) = r.call(&read).unwrap();
    assert_eq!(r.alleles[h.primary].name, "IGHV3-23*01");
    assert_eq!(h.mismatches, 0);
    let ins: Vec<usize> = (0..read.len()).filter(|&i| i < v.len() + 9 && h.imgt[i].is_none()).collect();
    assert_eq!(ins.len(), 9);
    assert!(h.aln.ops.contains(&(Op::Ins, 9)));
}

/// Recoverable reads: the true V allele is in its non-empty valid set.
fn recoverable(g: &Generated, i: usize) -> bool {
    !g.list("v_call_valid", i).is_empty() && !g.bool("unrecoverable_v", i).unwrap()
}

#[test]
fn prefilter_never_drops_a_tie() {
    // The tie set holds every allele of v_call_valid on recoverable reads (seed 1, 2,000 reads,
    // mutation uniform 0-0.3). Record 0029: the read alone cannot always reproduce truth's
    // origin-based sets (N bases matching a sibling's 3' tail, spurious compensating gaps at
    // > 20% mutation), so the bound is 99%, and the exact rate is printed.
    let r = human();
    let g = generated(1, "naive", 2000, |p| p.mutation = Spec::Uniform(0.0, 0.3));
    let rows: Vec<usize> = (0..g.len()).filter(|&i| recoverable(&g, i)).collect();
    let missed: Vec<String> = rows
        .par_iter()
        .filter_map(|&i| {
            let valid: BTreeSet<String> = g.list("v_call_valid", i).into_iter().collect();
            let ties = r.call(g.sequence(i).as_bytes()).map(|(_, h)| r.names(&h.ties)).unwrap_or_default();
            (!ties.is_superset(&valid)).then(|| format!("{}: valid {valid:?} ties {ties:?}", g.str("sequence_id", i).unwrap()))
        })
        .collect();
    let rate = 1.0 - missed.len() as f64 / rows.len() as f64;
    println!("V tie coverage (mutation 0-0.3): {}/{} = {:.4}%", rows.len() - missed.len(), rows.len(), 100.0 * rate);
    assert!(rate >= 0.99, "{} reads dropped a valid allele: {:#?}", missed.len(), &missed[..missed.len().min(10)]);
}

#[test]
fn v_call_accuracy_naive() {
    let r = human();
    let g = generated(1, "naive", 2000, |_| {});
    let rows: Vec<usize> = (0..g.len()).filter(|&i| recoverable(&g, i)).collect();
    let calls: Vec<Option<VHit>> = rows.par_iter().map(|&i| r.call(g.sequence(i).as_bytes()).map(|(_, h)| h)).collect();
    let (mut correct, mut covered) = (0, 0);
    for (&i, h) in rows.iter().zip(&calls) {
        let valid: BTreeSet<String> = g.list("v_call_valid", i).into_iter().collect();
        if let Some(h) = h {
            correct += valid.contains(&r.alleles[h.primary].name) as usize;
            covered += r.names(&h.ties).is_superset(&valid) as usize;
        }
    }
    // Per-read time of the V step (orientation and V assignment), one thread, 200 reads.
    let t = Instant::now();
    for &i in rows.iter().take(200) {
        r.call(g.sequence(i).as_bytes());
    }
    let per_read = t.elapsed().as_secs_f64() / rows.len().min(200) as f64;
    let rate = correct as f64 / rows.len() as f64;
    println!(
        "V call accuracy (naive): {correct}/{} = {:.4}%; tie coverage {covered}/{}; V step {:.1} us/read ({})",
        rows.len(),
        100.0 * rate,
        rows.len(),
        per_read * 1e6,
        if cfg!(debug_assertions) { "debug" } else { "release" }
    );
    assert!(rate >= 0.99, "V call accuracy {rate}");
}
