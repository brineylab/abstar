//! Joint junction resolution (record 0030), D in both orientations (record 0031) and C
//! assignment on the built-in human database.
mod common;
use abstar_annotate::align::Scoring;
use abstar_annotate::assign::{assign_c, assign_j, primary_by_tiebreak, resolve_junction};
use abstar_generate::Spec;
use abstar_germline::{Locus, Segment};
use common::{generated, generated_inverting, human, revcomp, Calls, Generated, Ref};
use rayon::prelude::*;
use std::collections::BTreeSet;
use std::time::Instant;

fn ungapped(r: &Ref, name: &str) -> Vec<u8> {
    r.alleles[r.idx(name)].sequence().into_bytes()
}

/// A base other than `b`.
fn other(b: u8) -> u8 {
    if b == b'A' { b'C' } else { b'A' }
}

fn calls(r: &Ref, read: &[u8]) -> Calls {
    let c = r.call_all(read).expect("V and J hits");
    assert!(!c.read.rev_comp);
    c
}

#[test]
fn n_matching_v_extends_v_end() {
    // Record 0025 (G5): an N base matching V's next germline base moves the V end onto it.
    let r = human();
    let v = ungapped(r, "IGKV3-20*01");
    let j = ungapped(r, "IGKJ1*01");
    let cut = v.len() - 3;
    let n = [v[cut], other(v[cut + 1]), other(v[cut + 2])];
    let read = [&v[..cut], &n, &j].concat();
    let c = calls(r, &read);
    assert_eq!(c.junction.v_end, cut + 1);
    assert_eq!(c.junction.j_start, cut + 3);
    assert!(c.junction.d.is_none());
}

#[test]
fn mutated_edge_base_goes_to_n() {
    // Record 0025: a mutated base at V's edge goes to N (0 beats -1).
    let r = human();
    let mut v = ungapped(r, "IGKV3-20*01");
    let j = ungapped(r, "IGKJ1*01");
    let last = v.len() - 1;
    v[last] = other(v[last]);
    let read = [&v[..], b"TTCCA", &j].concat();
    let c = calls(r, &read);
    assert_eq!(c.junction.v_end, last);
    assert_eq!(c.junction.j_start, v.len() + 5);
}

#[test]
fn short_d_is_null() {
    // 3 nt of D between full V and J: the best D run is 3 nt, so no D; the no-D cuts stand.
    let r = human();
    let v = ungapped(r, "IGHV3-23*01");
    let d = ungapped(r, "IGHD3-3*01");
    let read = [&v[..], &d[5..8], &ungapped(r, "IGHJ4*02")].concat();
    let c = calls(r, &read);
    assert!(c.junction.d.is_none(), "{:?}", c.junction.d);
    assert_eq!((c.junction.v_end, c.junction.j_start), (v.len(), v.len() + 3));
}

#[test]
fn inverted_d_found() {
    let r = human();
    let v = ungapped(r, "IGHV3-23*01");
    let d = revcomp(&ungapped(r, "IGHD3-3*01"));
    let read = [&v[..], b"CCTA", &d, b"GTCA", &ungapped(r, "IGHJ4*02")].concat();
    let c = calls(r, &read);
    let h = c.junction.d.expect("D call");
    assert_eq!(r.alleles[h.primary].name, "IGHD3-3*01");
    assert_eq!(h.ties[0], h.primary);
    assert!(h.inverted);
    let start = v.len() + 4;
    assert_eq!((h.start, h.end, h.g_start, h.g_end), (start, start + d.len(), 0, d.len()));
    assert_eq!((c.junction.v_end, c.junction.j_start), (v.len(), start + d.len() + 4));
}

#[test]
fn vj_locus_has_no_d() {
    let r = human();
    let c = calls(r, &[&ungapped(r, "IGKV3-20*01")[..], b"TTCCA", &ungapped(r, "IGKJ1*01")].concat());
    assert!(c.junction.d.is_none());
    // A heavy read with a full D, resolved as a locus without D: V end and J start only.
    let v = ungapped(r, "IGHV3-23*01");
    let read = [&v[..], b"CC", &ungapped(r, "IGHD3-3*01"), b"GG", &ungapped(r, "IGHJ4*02")].concat();
    let c = calls(r, &read);
    assert!(c.junction.d.is_some());
    let j = resolve_junction(&c.read, &c.v, &c.j, &r.alleles, &r.d, false);
    assert!(j.d.is_none());
    assert!(j.v_end >= v.len() && j.v_end <= j.j_start);
}

#[test]
fn c_call_igg1() {
    let r = human();
    let g1 = ungapped(r, "IGHG1*01");
    let vdj = [&ungapped(r, "IGHV3-23*01")[..], b"GCGAAAGATCGGGGATATAGTAGTGGCTGGTAC", &ungapped(r, "IGHJ4*02")].concat();
    let read = [&vdj[..], &g1[..30]].concat();
    let c = calls(r, &read);
    let h = c.c.expect("C call");
    // IGH C alleles identical to IGHG1*01 over its first 30 nt.
    let expected: Vec<usize> = (0..r.alleles.len())
        .filter(|&i| {
            let a = &r.alleles[i];
            a.segment == Segment::C && a.locus == Locus::IGH && a.sequence().as_bytes().get(..30) == Some(&g1[..30])
        })
        .collect();
    assert!(r.names(&expected).contains("IGHG1*01") && expected.len() > 1, "{:?}", r.names(&expected));
    assert_eq!(r.names(&h.ties), r.names(&expected));
    assert_eq!(h.primary, primary_by_tiebreak(&expected, &r.alleles));
    assert_eq!(h.ties[0], h.primary);
    assert_eq!((h.start, h.end, h.g_end), (vdj.len(), read.len(), 30));
    // Fewer than 10 nt after J: no C call.
    let c = calls(r, &[&vdj[..], &g1[..9]].concat());
    assert!(assign_c(&c.read, &c.j, &r.alleles, &r.c).is_none());
    assert!(calls(r, &[&vdj[..], &g1[..10]].concat()).c.is_some());
}

/// Recoverable reads: V, J and D are each absent or have the true allele in their valid set,
/// and V and J are present.
fn recoverable(g: &Generated, i: usize) -> bool {
    let bad = |s: &str| g.bool(&format!("unrecoverable_{s}"), i).unwrap();
    !g.list("v_call_valid", i).is_empty() && !g.list("j_call_valid", i).is_empty() && !bad("v") && !bad("d") && !bad("j")
}

#[derive(Default)]
struct Tally {
    n: usize,
    v_end: usize,
    d_start: usize,
    d_end: usize,
    j_start: usize,
    /// The J alignment's own start in `j_sequence_start_valid` (before the junction step).
    j_aln: usize,
    all: usize,
    /// The D call is correct: its primary is in `d_call_valid`, or it is null and
    /// `d_call_null_valid`.
    d_call: usize,
    /// The D call is placed right as far as the D call can show: null and
    /// `d_call_null_valid`, or in `d_call_valid` in the true orientation with a true D run of
    /// at least 5 nt (a shorter run puts every allele containing its few bases in
    /// `d_call_valid`, wherever the call is placed).
    d_right: usize,
    /// All four boundaries correct on `d_right` reads.
    all_d_right: usize,
    /// (correct, graded): reads where the annotator and truth both have a D.
    d_inverted: (usize, usize),
    /// (called inverted, graded): reads whose true D is inverted and at least 5 nt.
    inverted_recall: (usize, usize),
    /// (correct, graded): reads with C bases.
    c_call: (usize, usize),
}

/// Grades one read: each boundary in its `_valid` list (a null D's boundaries, like the null D
/// call itself, are correct iff `d_call_null_valid`), the D and C calls and `d_inverted`.
fn grade_read(r: &Ref, g: &Generated, i: usize, c: Option<&Calls>, t: &mut Tally, missed: &mut Vec<String>) {
    t.n += 1;
    let Some(c) = c else { return };
    let ints = |name: &str| -> BTreeSet<i64> {
        use arrow_array::{Array, Int64Array, ListArray};
        let a = g.0.column_by_name(name).unwrap().as_any().downcast_ref::<ListArray>().unwrap().value(i);
        let a = a.as_any().downcast_ref::<Int64Array>().unwrap();
        (0..a.len()).map(|x| a.value(x)).collect()
    };
    let null_ok = g.bool("d_call_null_valid", i).unwrap();
    let jn = &c.junction;
    let v_end = ints("v_sequence_end_valid").contains(&(jn.v_end as i64));
    let j_start = ints("j_sequence_start_valid").contains(&(jn.j_start as i64 + 1));
    let (d_start, d_end, d_call) = match &jn.d {
        Some(d) => (
            ints("d_sequence_start_valid").contains(&(d.start as i64 + 1)),
            ints("d_sequence_end_valid").contains(&(d.end as i64)),
            g.list("d_call_valid", i).contains(&r.alleles[d.primary].name),
        ),
        None => (null_ok, null_ok, null_ok),
    };
    t.v_end += v_end as usize;
    t.j_start += j_start as usize;
    t.j_aln += ints("j_sequence_start_valid").contains(&(c.j.aln.q_start as i64 + 1)) as usize;
    t.d_start += d_start as usize;
    t.d_end += d_end as usize;
    t.d_call += d_call as usize;
    let d_right = match &jn.d {
        Some(d) => d_call && !null_ok && g.bool("d_inverted", i) == Some(d.inverted),
        None => null_ok,
    };
    t.d_right += d_right as usize;
    let all = v_end && j_start && d_start && d_end;
    t.all += all as usize;
    t.all_d_right += (all && d_right) as usize;
    if d_right && !all && missed.len() < 10 {
        missed.push(format!(
            "{}: v_end {} {:?} j_start {} {:?} d {:?} d_start {:?} d_end {:?} null_ok {null_ok}",
            g.str("sequence_id", i).unwrap(),
            jn.v_end,
            ints("v_sequence_end_valid"),
            jn.j_start + 1,
            ints("j_sequence_start_valid"),
            jn.d.as_ref().map(|d| (r.alleles[d.primary].name.clone(), d.start + 1, d.end, d.inverted)),
            ints("d_sequence_start_valid"),
            ints("d_sequence_end_valid"),
        ));
    }
    if let (Some(d), Some(inv)) = (&jn.d, g.bool("d_inverted", i)) {
        t.d_inverted.1 += 1;
        t.d_inverted.0 += (d.inverted == inv) as usize;
    }
    if g.bool("d_inverted", i) == Some(true) && !null_ok {
        t.inverted_recall.1 += 1;
        t.inverted_recall.0 += jn.d.as_ref().is_some_and(|d| d.inverted) as usize;
    }
    let c_valid = g.list("c_call_valid", i);
    if !c_valid.is_empty() {
        t.c_call.1 += 1;
        t.c_call.0 += c.c.as_ref().is_some_and(|h| c_valid.contains(&r.alleles[h.primary].name)) as usize;
    }
}

/// Grades every recoverable read of `g` and times the junction step (one thread, 200 reads);
/// printed. Returns the tally and up to 10 reads with a correct D call but a wrong boundary.
fn grade(label: &str, g: &Generated) -> (Tally, Vec<String>) {
    let r = human();
    let rows: Vec<usize> = (0..g.len()).filter(|&i| recoverable(g, i)).collect();
    let all: Vec<Option<Calls>> = rows.par_iter().map(|&i| r.call_all(g.sequence(i).as_bytes())).collect();
    let (mut t, mut missed) = (Tally::default(), Vec::new());
    for (&i, c) in rows.iter().zip(&all) {
        grade_read(r, g, i, c.as_ref(), &mut t, &mut missed);
    }
    let s = Scoring::default();
    let vj: Vec<_> = rows
        .iter()
        .take(200)
        .filter_map(|&i| {
            let (o, v) = r.call(g.sequence(i).as_bytes())?;
            let j = assign_j(&o, &v, &r.alleles, &r.j, &s)?;
            Some((o, v, j))
        })
        .collect();
    let start = Instant::now();
    for (o, v, j) in &vj {
        resolve_junction(o, v, j, &r.alleles, &r.d, true);
    }
    let per_read = start.elapsed().as_secs_f64() / vj.len().max(1) as f64;
    let pct = |x: usize, n: usize| 100.0 * x as f64 / n.max(1) as f64;
    let n = t.n;
    println!(
        "junction {label}: {n} reads; v_end {:.2}% d_start {:.2}% d_end {:.2}% j_start {:.2}% (J alignment start {:.2}%) all four {:.2}%; \
         D call {:.2}%; D placed right {:.2}%, all four then {:.2}%; d_inverted {}/{} = {:.2}%; inverted D recall {}/{} = {:.2}%; C call {}/{} = {:.2}%; junction step {:.1} us/read ({})",
        pct(t.v_end, n),
        pct(t.d_start, n),
        pct(t.d_end, n),
        pct(t.j_start, n),
        pct(t.j_aln, n),
        pct(t.all, n),
        pct(t.d_call, n),
        pct(t.d_right, n),
        pct(t.all_d_right, t.d_right),
        t.d_inverted.0,
        t.d_inverted.1,
        pct(t.d_inverted.0, t.d_inverted.1),
        t.inverted_recall.0,
        t.inverted_recall.1,
        pct(t.inverted_recall.0, t.inverted_recall.1),
        t.c_call.0,
        t.c_call.1,
        pct(t.c_call.0, t.c_call.1),
        per_read * 1e6,
        if cfg!(debug_assertions) { "debug" } else { "release" }
    );
    (t, missed)
}

#[test]
fn boundary_accuracy_naive() {
    // V end and J start are graded on every recoverable read; all four boundaries on reads whose
    // D is placed right (`Tally::d_right`). A short D is often outscored by a chance D placement
    // that truth, which knows the D-derived bases, cannot accept (record 0030, ruling A26), so
    // the D call itself is reported (about 96%), not asserted.
    let (t, missed) = grade("naive, seed 1", &generated(1, "naive", 2000, |_| {}));
    let rate = |x: usize, n: usize| x as f64 / n as f64;
    assert!(rate(t.v_end, t.n) >= 0.99 && rate(t.j_start, t.n) >= 0.99, "v_end {} j_start {}", rate(t.v_end, t.n), rate(t.j_start, t.n));
    assert!(rate(t.all_d_right, t.d_right) >= 0.99, "boundaries when the D is placed right {}: {missed:#?}", rate(t.all_d_right, t.d_right));
}

#[test]
fn inverted_d_recall() {
    // Record 0031: with the inverted-D margin, most inverted Ds (true run >= 5 nt) are still
    // called inverted on naive reads from a model that inverts half its Ds.
    let (t, _) = grade("naive, d_inversion 0.5, seed 3", &generated_inverting(3, "naive", 2000, 0.5, |_| {}));
    let recall = t.inverted_recall.0 as f64 / t.inverted_recall.1 as f64;
    assert!(t.inverted_recall.1 > 500 && recall >= 0.9, "inverted D recall {recall}");
}

#[test]
fn junction_accuracy_report() {
    // Printed only (not asserted): the Task 5 sets, and a naive set with 40 nt of C.
    let sets = [
        ("mutation 0-30%, seed 1", generated(1, "naive", 2000, |p| p.mutation = Spec::Uniform(0.0, 0.3))),
        ("memory + 20% indels, seed 11", generated(11, "memory", 3000, |p| p.indel_rate = 0.2)),
        ("bnab, seed 1", generated(1, "bnab", 2000, |_| {})),
        ("naive + 40 nt C, seed 2", generated(2, "naive", 2000, |p| p.constant_len = Some(40))),
        ("mutation 0-30%, d_inversion 0.5, seed 3", generated_inverting(3, "naive", 2000, 0.5, |p| p.mutation = Spec::Uniform(0.0, 0.3))),
    ];
    for (label, g) in &sets {
        grade(label, g);
    }
}
