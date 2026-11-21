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
    // An IGK V followed by an IGH J: only IGK J alleles are searched. None fits the IGH J; the
    // Wave 3 low path may still call a chance IGK match (accepted false positive), but never an
    // IGH allele.
    let r = human();
    let read = [&ungapped(r, "IGKV3-20*01")[..], NP, &ungapped(r, "IGHJ4*02")].concat();
    let (o, v) = r.call(&read).expect("V hit");
    assert_eq!(r.alleles[v.primary].locus, Locus::IGK);
    if let Some(h) = assign_j(&o, &v, &r.alleles, &r.j, &Scoring::default()) {
        assert!(h.ties.iter().all(|&i| r.alleles[i].locus == Locus::IGK), "{:?}", r.names(&h.ties));
    }
    // An IGK J after it is called, from IGK.
    let read = [&ungapped(r, "IGKV3-20*01")[..], NP, &ungapped(r, "IGKJ1*01")].concat();
    let (o, v) = r.call(&read).expect("V hit");
    let h = assign_j(&o, &v, &r.alleles, &r.j, &Scoring::default()).expect("J hit");
    let loci: BTreeSet<Locus> = h.ties.iter().map(|&i| r.alleles[i].locus).collect();
    assert_eq!(loci, [Locus::IGK].into(), "{:?}", r.names(&h.ties));
    assert_eq!(r.alleles[h.primary].name, "IGKJ1*01");
}

#[test]
fn missing_j_is_none() {
    let r = human();
    let v = ungapped(r, "IGHV3-23*01");
    assert!(r.call_j(&v).is_none(), "V only");
    // A read that ends in N or D bases (no J base) gets a J exactly when its tail, laid on some
    // IGH J allele's 5' end, scores at least the low path's floor of 6 on the ±2 scale (3 net
    // matches): such a tail cannot be told from a truncated J (record 0029, Wave 3). On the NP
    // stretch, 27 truncation points stay below it and get no J; 6 reach it.
    let js: Vec<Vec<u8>> = r.alleles.iter().filter(|a| a.segment == Segment::J && a.locus == Locus::IGH).map(|a| a.sequence().into_bytes()).collect();
    let five = |q: &[u8]| {
        let score = |t: &[u8], g: &[u8]| t.iter().zip(g).map(|(a, b)| if a == b { 2 } else { -2 }).sum::<i32>();
        js.iter().flat_map(|g| (1..=q.len().min(g.len())).map(|l| score(&q[q.len() - l..], &g[..l]))).max().unwrap_or(0)
    };
    let mut below = 0;
    for k in 1..=NP.len() {
        let called = r.call_j(&[&v[..], &NP[..k]].concat()).is_some();
        assert_eq!(called, five(&NP[..k]) >= 6, "read ends {k} nt into N/D, 5' tail score {}", five(&NP[..k]));
        below += usize::from(!called);
    }
    assert_eq!(below, 27);
    // V followed by random bases of 1-120 nt. The random-tail rate is reported, not a bound: the
    // user accepts false J calls on junk after a V (record 0029, Wave 3). 309 of these 1,000 get a
    // J (11 in Wave 2, 64 before the low path); the cap only catches a gross change.
    let called: Vec<String> = (0..1000u64)
        .into_par_iter()
        .filter_map(|s| {
            let read = [&v[..], &rand_seq(1 + s as usize % 120, s + 1)].concat();
            r.call_j(&read).map(|h| format!("seed {s}: {} score {}", r.alleles[h.primary].name, h.aln.score))
        })
        .collect();
    assert!(called.len() <= 400, "random tails given a J: {called:?}");
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
    // Ruling B1 as amended in Waves 2 and 3: a 3' J fragment is called when its matches minus
    // mismatches reach 8 (IGH, IGK, IGL; 10 in Wave 2). An exact 8-nt fragment is called (Wave 2
    // left it uncalled), and a 12-nt fragment with 2 mismatches (score 8) is called. A 7-nt
    // fragment is below B1 and goes to the low path (`low_path_calls_a_short_fragment`).
    let r = human();
    let (v, j) = (ungapped(r, "IGHV3-23*01"), ungapped(r, "IGHJ4*02"));
    let h = r.call_j(&[&v[..], NP, &j[..8]].concat()).expect("J hit");
    assert_eq!(r.alleles[h.primary].name, "IGHJ4*02", "ties {:?}", r.names(&h.ties));
    assert_eq!((h.aln.q_end - h.aln.q_start, h.aln.g_start, h.mismatches), (8, 0, 0));
    let mut two = j[..12].to_vec();
    two[5] = transition(two[5]);
    two[10] = transition(two[10]);
    let h = r.call_j(&[&v[..], NP, &two].concat()).expect("J hit");
    assert_eq!(r.alleles[h.primary].name, "IGHJ4*02", "ties {:?}", r.names(&h.ties));
    assert_eq!((h.aln.q_end - h.aln.q_start, h.aln.g_start, h.mismatches), (12, 0, 2));
}

#[test]
fn low_path_calls_a_short_fragment() {
    // Record 0029, Wave 3: an exact 7-nt IGHJ4*02 fragment is below B1 (score 7 < 8) and holds no
    // Trp118 for the weak J; the low path (gapless cut scoring >= 6, junction bounded) calls it,
    // as the read's own 7 J bases.
    let r = human();
    let (v, j) = (ungapped(r, "IGHV3-23*01"), ungapped(r, "IGHJ4*02"));
    let read = [&v[..], NP, &j[..7]].concat();
    let h = r.call_j(&read).expect("J hit");
    assert!(r.names(&h.ties).contains("IGHJ4*02"), "{:?}", r.names(&h.ties));
    assert_eq!(r.alleles[h.primary].name, "IGHJ4*02", "ties {:?}", r.names(&h.ties));
    assert_eq!((h.aln.q_start, h.aln.q_end, h.aln.g_start, h.aln.g_end), (read.len() - 7, read.len(), 0, 7));
    assert_eq!(h.aln.ops, vec![(Op::Match, 7)]);
}

#[test]
fn short_fragment_needs_read_end_and_a_bounded_junction() {
    // Case 1, rejected by score: a 12-nt J fragment followed by 10 bases that mismatch the J's
    // continuation at every base. Every alignment runs to the read's or the J's end, so the cut
    // carries those 10 mismatches: B1 needs the fragment at the read's end, the weak J cannot
    // reach Trp118, and the low path's cut scores 24 - 20 = 4 < 6.
    let r = human();
    let (v, j) = (ungapped(r, "IGHV3-23*01"), ungapped(r, "IGHJ4*02"));
    let other: Vec<u8> = j[12..22].iter().map(|&b| transversion(b)).collect();
    assert!(r.call_j(&[&v[..], NP, &j[..12], &other].concat()).is_none(), "fragment not at the read end");
    // Case 2, rejected by the junction bound alone: an exact 14-nt fragment at the read's end
    // scores enough for B1 (14) and the low path (28), and is called after NP; after 200 random
    // bases its implied junction is over 180 nt and every path refuses it.
    assert!(r.call_j(&[&v[..], NP, &j[..14]].concat()).is_some(), "control: the fragment's score suffices");
    assert!(r.call_j(&[&v[..], &rand_seq(200, 3), &j[..14]].concat()).is_none(), "implied junction over 180 nt");
}

/// Benchmark igh_truncated g1554 (IGHV5-10-1*03, the read ends in N/D: no true J base).
const TRUNC_G1554: &str = "AGGTGAAAGAGCCCGGGGAGTCTCTGAGGATCTCCTGTAAGGGTTCTGGATACAGCTTCACCAGCTACTGGATCCGGTGGGTGTGCCAGATGTCCGGGAAAGGCCTGGAGTGGATGGGGAGGATTGGTCCTAGTGACTCTTACACCAACTACAGCCCGTCCTTCCAAGGCCACGTCACCATCTCAGCTGACTAGTCAATCAGCACTGTCTACCTACAGTGGAGTAGCCTGAAGGCCTCGGACACCGCCATGTATCACTGTGTGAGAGGGATCCGATATAATTGGACCGACCTTTGGGGGGGGACCAG";

#[test]
fn low_path_needs_a_gapless_cut() {
    // Record 0029, Wave 3: the low path takes only gapless cuts. This read's best normal cut is
    // IGHJ5*02 `13M5I3M` (a chance match in N/D bases, truth has no J); without the gapless guard
    // it would be called. Measured on igh_truncated, the guard changes 6 calls (3 right, 3 wrong).
    let r = human();
    assert!(r.call_j(TRUNC_G1554.as_bytes()).is_none());
}

#[test]
fn v_only_reads_get_no_j() {
    // Wave 3 fix round 1: the J search starts 10 nt before V's end, so a fallback could call a J
    // from V's own 3' bases. Every fallback's evidence must now lie past V's end. Every human V
    // allele, whole and cut back 1-12 nt, with nothing after it: no J (178 + 2,983 of these reads
    // got one before the fix).
    let r = human();
    let reads: Vec<Vec<u8>> = r
        .alleles
        .iter()
        .filter(|a| a.segment == Segment::V)
        .flat_map(|a| {
            let s = a.sequence().into_bytes();
            (0..=12usize).filter(|&c| s.len() > c + 60).map(|c| s[..s.len() - c].to_vec()).collect::<Vec<_>>()
        })
        .collect();
    assert!(reads.len() > 7000, "{}", reads.len());
    let called: Vec<String> = reads
        .par_iter()
        .filter_map(|read| r.call_j(read).map(|h| format!("{} len {}: {} {:?}", r.alleles[r.call(read).unwrap().1.primary].name, read.len(), r.alleles[h.primary].name, h.aln)))
        .collect();
    assert!(called.is_empty(), "{} V-only reads got a J: {:?}", called.len(), &called[..called.len().min(5)]);
}

#[test]
fn heavily_mutated_j_called() {
    // Record 0029, Wave 2: a full-length IGHJ1*01 with 20 of 52 bases mutated (38%, the bnab
    // median), its Trp118 codon intact, scores below the J floor of 30 (no call before Wave 2)
    // but is called on the weak path, with the anchor codon at its true read position.
    let r = human();
    let j = ungapped(r, "IGHJ1*01");
    let w = r.alleles[r.idx("IGHJ1*01")].anchor.unwrap();
    let m: Vec<u8> = j.iter().enumerate().map(|(i, &b)| if i % 5 % 2 == 1 && !(w..w + 3).contains(&i) { transition(b) } else { b }).collect();
    assert_eq!(m.iter().zip(&j).filter(|(a, b)| a != b).count(), 20);
    let v = ungapped(r, "IGHV3-23*01");
    let start = v.len() + NP.len();
    let h = r.call_j(&[&v[..], NP, &m].concat()).expect("J hit");
    assert_eq!(r.alleles[h.primary].name, "IGHJ1*01", "ties {:?}", r.names(&h.ties));
    assert!(h.aln.g_start <= w && w + 3 <= h.aln.g_end, "anchor codon in the cut: {:?}", h.aln);
    assert_eq!(h.aln.q_start - h.aln.g_start, start, "J placed at its true offset");
}

fn transition(b: u8) -> u8 {
    match b {
        b'A' => b'G',
        b'G' => b'A',
        b'C' => b'T',
        _ => b'C',
    }
}

fn transversion(b: u8) -> u8 {
    match b {
        b'A' => b'C',
        b'C' => b'A',
        b'G' => b'T',
        _ => b'G',
    }
}

#[test]
fn later_j_start_breaks_iglj2_iglj3_tie() {
    // B4: N bases TG match IGLJ3*02's 5' end (TTGGGTG...) by chance; the J itself is IGLJ2*01 from
    // its 4th base (GGTATTC...), where IGLJ3*02 has a mismatch (GGTGTTC...).
    let r = human();
    let j2 = ungapped(r, "IGLJ2*01");
    let read = [&ungapped(r, "IGLV2-14*01")[..], b"CTGCATG", &j2[3..]].concat();
    let h = r.call_j(&read).expect("J hit");
    assert_eq!(r.alleles[h.primary].name, "IGLJ2*01", "ties {:?}", r.names(&h.ties));
}

/// Benchmark igh_bnab_c g3609 (IGHV3-7*01, IGHJ4*02 mutated 18 of 45 nt, then 60 nt of C).
const BNAB_C_G3609: &str = "GGGACTCACATGGTGAAATCTGGGGGAAGTTTAGTCCACACTGGTAACCCCCGGGGACCCGCCTCTGGCGTCTCTGATTTCACGACTAATTACGATTGGGTGAATTGGACTGGTCAGGGTCCAGGGACGGAGCTTGAGTGGGTGGTCACTATGATGTATGGAGGCAGTGAACAATATGAGGACTCTGGGCCGGGTCATTATACTAAGTCCAGCAGTCATGTCGCGAATTCAGTGTTTCCTCAGATGATCAATCTGAGAGGCAAGGACCCGCCACTGGCTCATTGTGGGCGAGTCCTAATAAAGGGAATGATACATGACATCCAACGTCCCCAGAAATCTCCTCATCTGATCCTTGAGTGGTGGACTCGTGCATCTCTGGTCACTGCCCCCTCCAGGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGATA";
/// Benchmark tra g0440 (TRAV10*02, an exact 47-nt TRAJ8*01 tail from germline base 14).
const TRA_G0440: &str = "AAAAACCAAGTGGAGCAGAGTCCTCAGTCCCTGATCATCCTGGAGGGAAAGAACTGCACTCTTCAATGCAATTATACAGTGAGCCCCTTCAGCAACTTAAGGTGGTATAAGCAAGATACGGGGAGAGGTCCTGTTTCCCTGACAATCATGACTTTCAGTGAGAACACAAAGTCGAACGGAAGATATACAGCAACTCTGGATGCAGACACAAAGCAAAGCTCTCTGCACATCACAGCCTCCCAGCTCAGCGATTCAGCCTCCTACATCTGTGTGGTGAGCGCTGGCTTCCCCTCAGAAACTTGTATTTGGAACTGGCACCCGACTTCTGGTCAGTCCAA";

#[test]
fn weak_j_guard_holds_after_the_record_0008_pick() {
    // Wave 2, fix round 1: the weak J's cut held IGHJ4*02's Trp118, but the re-cut on the picked
    // primary (IGHJ2*01) started past its anchor (366S25N28M60S), so the reported J had no anchor.
    // The guard is now re-applied to the reported alignment; failing it, B1 is tried, then no J.
    let r = human();
    if let Some(h) = r.call_j(BNAB_C_G3609.as_bytes()) {
        let w = r.alleles[h.primary].anchor.unwrap();
        assert!(h.aln.g_start <= w && w + 3 <= h.aln.g_end, "{}: {:?}", r.alleles[h.primary].name, h.aln);
    }
}

#[test]
fn exact_3p_fragment_beats_a_weaker_weak_cut() {
    // Wave 2, fix round 1: an exact 47-nt TRAJ8*01 tail (B1 score 47, 94 on the ±2 scale) is called
    // over a weak cut scoring less, and is reported as the exact gapless fragment.
    let r = human();
    let read = TRA_G0440.as_bytes();
    let h = r.call_j(read).expect("J hit");
    assert_eq!(r.alleles[h.primary].name, "TRAJ8*01", "ties {:?}", r.names(&h.ties));
    assert_eq!((h.aln.q_start, h.aln.q_end, h.aln.g_start, h.aln.g_end, h.mismatches), (291, 338, 13, 60, 0));
    assert_eq!(h.aln.ops, vec![(Op::Match, 47)]);
}
