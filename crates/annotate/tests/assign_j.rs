//! J assignment (record 0029) on the built-in human database.
mod common;
use abstar_annotate::align::{Op, Scoring};
use abstar_annotate::assign::{assign_j, primary_by_tiebreak, JHit};
use abstar_generate::Spec;
use abstar_germline::{Locus, Segment};
use common::{generated, human, rand_seq, Generated, Ref};
use rayon::prelude::*;
use std::collections::BTreeSet;
use std::time::Instant;

fn ungapped(r: &Ref, name: &str) -> Vec<u8> {
    r.alleles[r.idx(name)].sequence().into_bytes()
}

/// `locus` J alleles containing `bases` (no mismatch at some shift, every base covered).
fn containing(r: &Ref, locus: Locus, bases: &[u8]) -> Vec<usize> {
    (0..r.alleles.len())
        .filter(|&i| {
            let a = &r.alleles[i];
            a.segment == Segment::J && a.locus == locus && a.sequence().as_bytes().windows(bases.len()).any(|w| w == bases)
        })
        .collect()
}

const NP: &[u8] = b"GCGAAAGATCGGGGATATAGTAGTGGCTGGTAC";

#[test]
fn exact_j_called() {
    let r = human();
    let j = ungapped(r, "IGHJ4*02");
    let read = [&ungapped(r, "IGHV3-23*01")[..], NP, &j].concat();
    let h = r.call_j(&read).expect("J hit");
    let expected = containing(r, Locus::IGH, &j);
    assert_eq!(r.names(&h.ties), r.names(&expected));
    assert_eq!(r.alleles[h.primary].name, "IGHJ4*02");
    assert_eq!(h.ties[0], h.primary);
    assert_eq!(h.mismatches, 0);
    let start = read.len() - j.len();
    assert_eq!((h.aln.q_start, h.aln.q_end, h.aln.g_start, h.aln.g_end), (start, read.len(), 0, j.len()));
    assert_eq!(h.aln.ops, vec![(Op::Match, j.len() as u32)]);
}

#[test]
fn short_j_sibling_tie() {
    // IGKJ2*01, *03 and *04 share their last 32 nt. A read keeping only those ties them; the
    // full-length IGKJ2*01 does not, since the siblings' different 5' bases are mismatches.
    let r = human();
    let full = ungapped(r, "IGKJ2*01");
    let v = ungapped(r, "IGKV3-20*01");
    let short = &full[full.len() - 32..];
    // N bases unlike every sibling's base before the shared stretch.
    let read = [&v[..], b"GGGCCCTTT", short].concat();
    let h = r.call_j(&read).expect("J hit");
    let expected = containing(r, Locus::IGK, short);
    assert!(r.names(&expected).is_superset(&["IGKJ2*01", "IGKJ2*03", "IGKJ2*04"].map(String::from).into()));
    assert_eq!(r.names(&h.ties), r.names(&expected));
    assert_eq!(h.primary, primary_by_tiebreak(&expected, &r.alleles));
    assert_eq!(h.mismatches, 0);
    let h = r.call_j(&[&v[..], b"GGGCCCTTT", &full].concat()).expect("J hit");
    assert_eq!(r.names(&h.ties), r.names(&containing(r, Locus::IGK, &full)));
    assert!(!r.names(&h.ties).contains("IGKJ2*03"));
    assert_eq!(r.alleles[h.primary].name, "IGKJ2*01");
}

#[test]
fn j_from_v_locus_only() {
    // An IGK V followed by an IGH J: the J call comes from IGK, or there is none.
    let r = human();
    let read = [&ungapped(r, "IGKV3-20*01")[..], NP, &ungapped(r, "IGHJ4*02")].concat();
    let (o, v) = r.call(&read).expect("V hit");
    assert_eq!(r.alleles[v.primary].locus, Locus::IGK);
    if let Some(h) = assign_j(&o, &v, &r.alleles, &r.j, &Scoring::default()) {
        let loci: BTreeSet<Locus> = h.ties.iter().map(|&i| r.alleles[i].locus).collect();
        assert_eq!(loci, [Locus::IGK].into(), "{:?}", r.names(&h.ties));
    }
}

#[test]
fn missing_j_is_none() {
    let r = human();
    let v = ungapped(r, "IGHV3-23*01");
    assert!(r.call_j(&v).is_none(), "V only");
    // V followed by random bases of 1-120 nt.
    let called: Vec<String> = (0..1000u64)
        .into_par_iter()
        .filter_map(|s| {
            let read = [&v[..], &rand_seq(1 + s as usize % 120, s + 1)].concat();
            r.call_j(&read).map(|h| format!("seed {s}: {} score {}", r.alleles[h.primary].name, h.aln.score))
        })
        .collect();
    assert!(called.is_empty(), "random tails given a J: {called:?}");
}

/// Recoverable reads: the true J allele is in its non-empty valid set.
fn recoverable(g: &Generated, i: usize) -> bool {
    !g.list("j_call_valid", i).is_empty() && !g.bool("unrecoverable_j", i).unwrap()
}

/// J primary in `j_call_valid`, `j_call_valid` ⊆ ties and no J hit, over recoverable reads,
/// and the J step's per-read time (one thread, 200 reads); printed. Returns the primary rate
/// and up to 10 misses.
fn grade(label: &str, g: &Generated) -> (f64, Vec<String>) {
    let r = human();
    let rows: Vec<usize> = (0..g.len()).filter(|&i| recoverable(g, i)).collect();
    let calls: Vec<Option<JHit>> = rows.par_iter().map(|&i| r.call_j(g.sequence(i).as_bytes())).collect();
    let (mut correct, mut covered, mut none) = (0, 0, 0);
    let mut missed = Vec::new();
    for (&i, h) in rows.iter().zip(&calls) {
        let valid: BTreeSet<String> = g.list("j_call_valid", i).into_iter().collect();
        match h {
            Some(h) => {
                let ok = valid.contains(&r.alleles[h.primary].name);
                correct += ok as usize;
                covered += r.names(&h.ties).is_superset(&valid) as usize;
                if !ok && missed.len() < 10 {
                    missed.push(format!("{}: valid {valid:?} ties {:?}", g.str("sequence_id", i).unwrap(), r.names(&h.ties)));
                }
            }
            None => none += 1,
        }
    }
    let vs: Vec<_> = rows.iter().take(200).filter_map(|&i| r.call(g.sequence(i).as_bytes())).collect();
    let t = Instant::now();
    for (o, v) in &vs {
        assign_j(o, v, &r.alleles, &r.j, &Scoring::default());
    }
    let per_read = t.elapsed().as_secs_f64() / vs.len().max(1) as f64;
    let n = rows.len();
    let pct = |x: usize| 100.0 * x as f64 / n as f64;
    println!(
        "J {label}: primary in valid {correct}/{n} = {:.2}%; valid in ties {covered}/{n} = {:.2}%; no J hit {none} ({:.2}%); J step {:.1} us/read ({})",
        pct(correct),
        pct(covered),
        pct(none),
        per_read * 1e6,
        if cfg!(debug_assertions) { "debug" } else { "release" }
    );
    (correct as f64 / n as f64, missed)
}

#[test]
fn j_call_accuracy_naive() {
    let (rate, missed) = grade("naive, seed 1", &generated(1, "naive", 2000, |_| {}));
    assert!(rate >= 0.99, "J call accuracy {rate}: {missed:#?}");
}

#[test]
fn j_call_accuracy_report() {
    // Printed only (not asserted): the Task 5 sets.
    grade("mutation 0-30%, seed 1", &generated(1, "naive", 2000, |p| p.mutation = Spec::Uniform(0.0, 0.3)));
    grade("memory + 20% indels, seed 11", &generated(11, "memory", 3000, |p| p.indel_rate = 0.2));
    grade("bnab, seed 1", &generated(1, "bnab", 2000, |_| {}));
}

#[test]
fn short_3p_fragment_called() {
    // Ruling B1: a read ending after the first 14 J nt (exact) gets a J call; 11 nt (below IGH's
    // 12-nt minimum) does not.
    let r = human();
    let (v, j) = (ungapped(r, "IGHV3-23*01"), ungapped(r, "IGHJ4*02"));
    let h = r.call_j(&[&v[..], NP, &j[..14]].concat()).expect("J hit");
    assert!(r.names(&h.ties).contains("IGHJ4*02"), "{:?}", r.names(&h.ties));
    assert_eq!(h.aln.q_end - h.aln.q_start, 14);
    assert!(r.call_j(&[&v[..], NP, &j[..11]].concat()).is_none());
}
