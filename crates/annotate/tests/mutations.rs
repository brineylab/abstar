//! Per-segment mutations, indels and masks (records 0033, 0034) on the built-in human database,
//! graded against generated truth by the grader's shared-span rule (`python/abstar/_grade.py`).
mod common;
use abstar_annotate::align::Scoring;
use abstar_annotate::assign::{assign_c, assign_j, resolve_junction, CHit, JHit, Junction, Oriented, VHit};
use abstar_annotate::derive::{masks, mutations, regions, Masks, SegmentMutations};
use abstar_generate::Spec;
use common::{generated, human, revcomp, Generated, Ref};
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet};

const IGHV: &str = "IGHV1-2*02";
const IGHD: &str = "IGHD3-3*01";
const IGHJ: &str = "IGHJ4*02";

fn ungapped(r: &Ref, name: &str) -> Vec<u8> {
    r.alleles[r.idx(name)].sequence().into_bytes()
}

/// `v` + CC + IGHD3-3*01 + GG + IGHJ4*02.
fn vdj(r: &Ref, v: &str) -> Vec<u8> {
    [&ungapped(r, v)[..], b"CC", &ungapped(r, IGHD), b"GG", &ungapped(r, IGHJ)].concat()
}

/// Ungapped index of IMGT position `p` (1-based) of a V allele.
fn at_imgt(r: &Ref, v: &str, p: usize) -> usize {
    let g = r.alleles[r.idx(v)].sequence_gapped.as_bytes();
    assert_ne!(g[p - 1], b'.');
    g[..p - 1].iter().filter(|&&b| b != b'.').count()
}

/// Every call for one read; J (and so the junction and C) may be missing.
struct Called {
    read: Oriented,
    v: VHit,
    j: Option<(JHit, Junction)>,
    c: Option<CHit>,
}

impl Called {
    fn j(&self) -> Option<(&JHit, &Junction)> {
        self.j.as_ref().map(|(j, jn)| (j, jn))
    }
    fn mutations(&self, r: &Ref) -> [Option<SegmentMutations>; 4] {
        mutations(&self.read, &self.v, self.j(), self.c.as_ref(), &r.alleles)
    }
    fn masks(&self, r: &Ref) -> Masks {
        let g = regions(&self.read, &self.v, self.j(), self.c.as_ref(), &r.alleles);
        masks(&self.read, &self.v, self.j(), &r.alleles, &g)
    }
}

fn call(r: &Ref, read: &[u8]) -> Option<Called> {
    let (o, v) = r.call(read)?;
    let j = assign_j(&o, &v, &r.alleles, &r.j, &Scoring::default()).map(|j| {
        let has_d = r.d.loci.contains(&r.alleles[v.primary].locus);
        let jn = resolve_junction(&o, &v, &j, &r.alleles, &r.d, has_d);
        (j, jn)
    });
    let c = j.as_ref().and_then(|(j, _)| assign_c(&o, j, &r.alleles, &r.c));
    Some(Called { read: o, v, j, c })
}

fn forward(r: &Ref, read: &[u8]) -> Called {
    let c = call(r, read).expect("V hit");
    assert!(!c.read.rev_comp);
    c
}

fn strs(x: &[&str]) -> Vec<String> {
    x.iter().map(|s| s.to_string()).collect()
}

#[test]
fn v_mutation_imgt_numbering() {
    let r = human();
    // IGHV3-64*01 has A at IMGT 153 (codon 51 GAA, Glu); A>G keeps Glu.
    let v = "IGHV3-64*01";
    let mut read = vdj(r, v);
    read[at_imgt(r, v, 153)] = b'G';
    let c = forward(r, &read);
    assert_eq!(r.alleles[c.v.primary].name, v);
    let m = c.mutations(r)[0].clone().expect("V mutations");
    assert_eq!(m.nt, strs(&["IMGT153:A>G"]));
    assert_eq!((m.aa, m.count, m.count_aa), (Some(vec![]), 1, Some(0)));
    assert!(m.insertions.is_empty() && m.deletions.is_empty());
    // IGHV1-2*02 codon 51 is GAG (Glu); A>G at IMGT 152 makes it GGG (Gly).
    let mut read = vdj(r, IGHV);
    read[at_imgt(r, IGHV, 152)] = b'G';
    let c = forward(r, &read);
    let m = c.mutations(r)[0].clone().expect("V mutations");
    assert_eq!(m.nt, strs(&["IMGT152:A>G"]));
    assert_eq!(m.aa, Some(strs(&["IMGT51:E>G"])));
    assert_eq!((m.count, m.count_aa, m.germline), (1, Some(1), (0, 296)));
    assert_eq!(m.frequency, 1.0 / 296.0);
    // 98 IMGT codons hold three aligned bases (gap codons and the partial last codon don't count)
    assert_eq!(m.frequency_aa, Some(1.0 / 98.0));
}

#[test]
fn d_and_j_mutations() {
    let r = human();
    let mut read = vdj(r, IGHV);
    // D base 10 (T) to C: D codon 3 (ATT, Ile, in the V frame) becomes ACT (Thr)
    read[296 + 2 + 9] = b'C';
    // J base 31 (T) to A: J codon 10 (CTG, Leu, frame from Trp118) becomes CAG (Gln)
    let j0 = 296 + 2 + 31 + 2;
    read[j0 + 30] = b'A';
    let c = forward(r, &read);
    let [v, d, j, cm] = c.mutations(r);
    assert!(cm.is_none());
    let v = v.unwrap();
    assert_eq!((v.nt.len(), v.aa.as_ref().map(Vec::len)), (0, Some(0)));
    let d = d.expect("D mutations");
    assert!(!c.j.as_ref().unwrap().1.d.as_ref().unwrap().inverted);
    assert_eq!(d.nt, strs(&["10:T>C"]));
    assert_eq!(d.aa, Some(strs(&["3:I>T"])));
    assert_eq!((d.count, d.frequency, d.germline), (1, 1.0 / 31.0, (0, 31)));
    // D codons ATT ACG ATT TTT GGA GTG GTT ATT ATA from base 3: 9 counted
    assert_eq!(d.frequency_aa, Some(1.0 / 9.0));
    let j = j.expect("J mutations");
    assert_eq!(j.nt, strs(&["31:T>A"]));
    assert_eq!(j.aa, Some(strs(&["10:L>Q"])));
    assert_eq!((j.frequency, j.frequency_aa, j.germline), (1.0 / 48.0, Some(1.0 / 15.0), (0, 48)));

    // An inverted D: positions are on the reverse-complemented allele (record 0031).
    let v = ungapped(r, "IGHV3-23*01");
    let mut drc = revcomp(&ungapped(r, IGHD));
    assert_eq!(drc[4], b'T');
    drc[4] = b'C';
    let read = [&v[..], b"CCTA", &drc, b"GTCA", &ungapped(r, IGHJ)].concat();
    let c = forward(r, &read);
    let h = c.j.as_ref().unwrap().1.d.clone().expect("D call");
    assert!(h.inverted && r.alleles[h.primary].name == IGHD, "{h:?}");
    let d = c.mutations(r)[1].clone().expect("D mutations");
    assert_eq!(d.nt, strs(&["5:T>C"]));
    assert_eq!(d.germline, (0, 31));
}

#[test]
fn frequency_denominator_truncated() {
    let r = human();
    let mut full = vdj(r, IGHV);
    for p in [100, 150] {
        full[p] = if full[p] == b'A' { b'C' } else { b'A' };
    }
    let m = forward(r, &full).mutations(r)[0].clone().unwrap();
    assert_eq!((m.count, m.frequency), (2, 2.0 / 296.0));
    // starting at V base 86: 211 aligned V positions
    let m = forward(r, &full[85..]).mutations(r)[0].clone().unwrap();
    assert_eq!((m.count, m.frequency, m.germline), (2, 2.0 / 211.0, (85, 296)));
}

#[test]
fn indel_strings() {
    let r = human();
    let v = ungapped(r, IGHV);
    let j = ungapped(r, IGHJ);
    // V: bases 121-123 (CCT, IMGT 136-138) deleted, TTT inserted after base 200 (IMGT 224);
    // J: bases 34-36 (TCA) deleted. Each placement is the only equivalent one.
    let vread = [&v[..120], &v[123..200], b"TTT", &v[200..]].concat();
    let jread = [&j[..33], &j[36..]].concat();
    let read = [&vread[..], b"CC", &ungapped(r, IGHD), b"GG", &jread].concat();
    let c = forward(r, &read);
    let [vm, _, jm, _] = c.mutations(r);
    let (vm, jm) = (vm.unwrap(), jm.unwrap());
    assert_eq!(vm.deletions, strs(&["IMGT136:del:CCT"]));
    assert_eq!(vm.insertions, strs(&["IMGT224:ins:TTT"]));
    assert_eq!((vm.nt.len(), vm.frequency), (0, 0.0));
    assert_eq!(jm.deletions, strs(&["34:del:TCA"]));
    assert!(jm.insertions.is_empty() && jm.nt.is_empty());
    // J codons from base 3: 15 complete; the two holding bases 34-36 lose deleted bases
    assert_eq!((jm.aa, jm.frequency_aa), (Some(vec![]), Some(0.0)));
}

#[test]
fn masks_hand_built() {
    let r = human();
    let c = forward(r, &vdj(r, IGHV));
    let m = c.masks(r);
    let rep = |parts: &[(&str, usize)]| parts.iter().map(|(s, n)| s.repeat(*n)).collect::<String>();
    assert_eq!(m.segment, rep(&[("V", 296), ("N", 2), ("D", 31), ("N", 2), ("J", 48)]));
    assert_eq!(m.cdr.as_deref(), Some(rep(&[("0", 75), ("1", 24), ("0", 51), ("2", 24), ("0", 114), ("3", 57), ("0", 34)]).as_str()));
    assert_eq!(m.nongermline, rep(&[("0", 296), ("1", 2), ("0", 31), ("1", 2), ("0", 48)]));
    assert_eq!(m.segment_aa, rep(&[("V", 98), ("N", 1), ("D", 10), ("N", 1), ("J", 16)]));
    assert_eq!(m.cdr_aa.as_deref(), Some(rep(&[("0", 25), ("1", 8), ("0", 17), ("2", 8), ("0", 38), ("3", 19), ("0", 11)]).as_str()));
    assert_eq!(m.nongermline_aa, rep(&[("0", 98), ("1", 2), ("0", 9), ("1", 2), ("0", 15)]));
    // a mismatch is nongermline; a silent one leaves the codon germline
    let mut read = vdj(r, IGHV);
    read[at_imgt(r, IGHV, 152)] = b'G';
    read[at_imgt(r, IGHV, 153)] = b'A';
    let m = forward(r, &read).masks(r);
    let p = at_imgt(r, IGHV, 152);
    assert_eq!(&m.nongermline[p - 1..p + 3], "0110");
    assert_eq!(&m.nongermline_aa[p / 3 - 1..p / 3 + 2], "010");
}

// ---- the grader's shared-span rule (python/abstar/_grade.py `mut_fn`, `mask_fn`), ported ----

const SEGS: [&str; 4] = ["v", "d", "j", "c"];
const PARTS: [&str; 8] =
    ["mutations", "mutations_aa", "insertions", "deletions", "mutation_count", "mutation_count_aa", "mutation_frequency", "mutation_frequency_aa"];
const MASKS: [&str; 6] = ["cdr_mask", "gene_segment_mask", "nongermline_mask", "cdr_mask_aa", "gene_segment_mask_aa", "nongermline_mask_aa"];
const BIG: i64 = 1_000_000_000;

/// The annotator's fields for one read, as the grader reads them (AIRR 1-based closed).
#[derive(Default)]
struct Row {
    call: [Option<String>; 4],
    germline: [Option<(i64, i64)>; 4],
    /// v_sequence_end, d_sequence_start, d_sequence_end, j_sequence_start
    cuts: BTreeMap<&'static str, i64>,
    v_start: i64,
    muts: [Option<SegmentMutations>; 4],
    masks: Option<Masks>,
}

fn row(r: &Ref, c: &Called) -> Row {
    let muts = c.mutations(r);
    let mut out = Row { germline: std::array::from_fn(|k| muts[k].as_ref().map(|m| (m.germline.0 as i64 + 1, m.germline.1 as i64))), ..Row::default() };
    out.call[0] = Some(r.alleles[c.v.primary].name.clone());
    out.v_start = c.v.aln.q_start as i64 + 1;
    if let Some((j, jn)) = &c.j {
        out.call[2] = Some(r.alleles[j.primary].name.clone());
        out.cuts.insert("v_sequence_end", jn.v_end as i64);
        out.cuts.insert("j_sequence_start", jn.j_start as i64 + 1);
        if let Some(d) = &jn.d {
            out.call[1] = Some(r.alleles[d.primary].name.clone());
            out.cuts.insert("d_sequence_start", d.start as i64 + 1);
            out.cuts.insert("d_sequence_end", d.end as i64);
        }
    } else {
        out.cuts.insert("v_sequence_end", c.v.aln.q_end as i64);
    }
    out.call[3] = c.c.as_ref().map(|x| r.alleles[x.primary].name.clone());
    out.muts = muts;
    out.masks = Some(c.masks(r));
    out
}

#[derive(Clone, Debug, PartialEq)]
enum Val {
    List(Vec<String>),
    Int(i64),
    Float(Option<f64>),
}

fn val_eq(a: &Val, b: &Val) -> bool {
    match (a, b) {
        (Val::Float(Some(x)), Val::Float(Some(y))) => (x - y).abs() <= 1e-9,
        _ => a == b,
    }
}

fn pos(s: &str) -> i64 {
    s.trim_start_matches("IMGT").split(':').next().unwrap().parse().unwrap()
}

fn ceil3(x: i64) -> i64 {
    -(-x).div_euclid(3)
}

/// One segment's four lists (nt, aa, insertions, deletions), each `None` when null.
type Lists = [Option<Vec<String>>; 4];

struct Grader<'a> {
    r: &'a Ref,
    t: &'a Generated,
}

impl Grader<'_> {
    fn allele(&self, name: &str) -> &abstar_germline::Allele {
        &self.r.alleles[self.r.idx(name)]
    }
    /// IMGT position (1-based) of each ungapped base of a V allele.
    fn gmap(&self, name: &str) -> Vec<i64> {
        self.allele(name).sequence_gapped.bytes().enumerate().filter(|&(_, b)| b != b'.').map(|(g, _)| g as i64 + 1).collect()
    }
    fn tlists(&self, i: usize, seg: &str) -> Lists {
        ["mutations", "mutations_aa", "insertions", "deletions"].map(|p| self.t.opt_list(&format!("{seg}_{p}"), i))
    }
    fn clists(m: &Option<SegmentMutations>) -> Lists {
        match m {
            None => [None, None, None, None],
            Some(m) => [Some(m.nt.clone()), m.aa.clone(), Some(m.insertions.clone()), Some(m.deletions.clone())],
        }
    }
    fn tgerm(&self, i: usize, k: usize) -> Option<(i64, i64)> {
        let s = SEGS[k];
        self.t.int(&format!("{s}_germline_start"), i).zip(self.t.int(&format!("{s}_germline_end"), i))
    }
    fn frame(&self, i: usize, seg: &str) -> Option<i64> {
        let t = self.t;
        match seg {
            "v" => None,
            "d" => {
                let phase = (t.int("d_sequence_start", i)? - t.int("v_sequence_end", i)? + t.int("v_germline_end", i)? - 1).rem_euclid(3);
                Some((t.int("d_germline_start", i)? - 1 - phase).rem_euclid(3))
            }
            _ => {
                let a = self.allele(t.str("j_call", i)?);
                let b = a.anchor? as i64;
                Some(if seg == "j" { b % 3 } else { (b - a.sequence().len() as i64).rem_euclid(3) })
            }
        }
    }
    /// (ulo, uhi, lo, hi, nlo, nhi)
    fn mkspan(&self, i: usize, seg: &str, ulo: i64, uhi: i64) -> [i64; 6] {
        if seg == "v" {
            let m = self.gmap(self.t.str("v_call", i).unwrap());
            if uhi > m.len() as i64 {
                return [ulo, uhi, 0, BIG, 0, BIG];
            }
            let (lo, hi) = if ulo <= uhi { (m[ulo as usize - 1], m[uhi as usize - 1]) } else { (1, 0) };
            return [ulo, uhi, lo, hi, ceil3(lo + 2), hi.div_euclid(3)];
        }
        match self.frame(i, seg) {
            None => [ulo, uhi, ulo, uhi, 0, BIG],
            Some(s0) => [ulo, uhi, ulo, uhi, 1.max(ceil3(ulo - s0 + 2)), (uhi - s0).div_euclid(3)],
        }
    }
    fn aa_denom(&self, i: usize, seg: &str, dels: &[String], sp: [i64; 6]) -> i64 {
        let (nlo, nhi) = (sp[4], sp[5]);
        if nhi >= BIG {
            return BIG;
        }
        let mut codons: BTreeSet<i64> = (nlo.max(1)..=nhi).collect();
        if seg == "v" {
            let name = self.t.str("v_call", i).unwrap();
            let gapped = self.allele(name).sequence_gapped.as_bytes();
            codons.retain(|&n| !gapped.get((3 * n - 3) as usize..((3 * n) as usize).min(gapped.len())).is_some_and(|c| c.contains(&b'.')));
            let m = self.gmap(name);
            for x in dels {
                let Some(u) = m.iter().position(|&p| p == pos(x)).map(|u| u as i64 + 1) else { continue };
                for k in 0..x.split(':').nth(2).unwrap().len() as i64 {
                    if ((u - 1 + k) as usize) < m.len() {
                        codons.remove(&ceil3(m[(u - 1 + k) as usize]));
                    }
                }
            }
        } else {
            let s0 = self.frame(i, seg).unwrap();
            for x in dels {
                for p in pos(x)..pos(x) + x.split(':').nth(2).unwrap().len() as i64 {
                    codons.remove(&((p - s0 - 1).div_euclid(3) + 1));
                }
            }
        }
        codons.len() as i64
    }
    fn val(&self, i: usize, seg: &str, part: &str, l: &Lists, sp: [i64; 6]) -> Val {
        let [ulo, uhi, lo, hi, nlo, nhi] = sp;
        let aa = part.ends_with("_aa");
        let mut muts: Vec<String> =
            l[aa as usize].clone().unwrap().into_iter().filter(|m| if aa { (nlo..=nhi).contains(&pos(m)) } else { (lo..=hi).contains(&pos(m)) }).collect();
        muts.sort();
        let ranged = |k: usize| {
            let mut x: Vec<String> = l[k].clone().unwrap().into_iter().filter(|m| (lo..=hi).contains(&pos(m))).collect();
            x.sort();
            x
        };
        match part {
            "mutations" | "mutations_aa" => Val::List(muts),
            "insertions" => Val::List(ranged(2)),
            "deletions" => Val::List(ranged(3)),
            "mutation_count" | "mutation_count_aa" => Val::Int(muts.len() as i64),
            _ => {
                let denom = if aa {
                    self.aa_denom(i, seg, l[3].as_ref().unwrap(), sp)
                } else {
                    let deleted: i64 = ranged(3).iter().map(|x| x.split(':').nth(2).unwrap().len() as i64).sum();
                    0.max(uhi - ulo + 1) - deleted
                };
                Val::Float((denom > 0).then(|| muts.len() as f64 / denom as f64))
            }
        }
    }
    fn span_consistent(&self, i: usize, c: &Row, k: usize, tg: (i64, i64), cg: (i64, i64)) -> bool {
        let cuts: [Option<&str>; 2] = match k {
            0 => [None, Some("v_sequence_end")],
            1 => [Some("d_sequence_start"), Some("d_sequence_end")],
            2 => [Some("j_sequence_start"), None],
            _ => [None, None],
        };
        [(tg.0, cg.0), (tg.1, cg.1)].iter().zip(cuts).all(|(&(t, cv), col)| match col {
            None => t == cv,
            Some(col) => {
                let (sc, st) = (c.cuts.get(col).copied(), self.t.int(col, i));
                let (Some(sc), Some(st)) = (sc, st) else { return false };
                (self.t.ints(&format!("{col}_valid"), i).contains(&sc) || sc == st) && cv - t == sc - st
            }
        })
    }
    fn recoverable(&self, i: usize, seg: &str) -> bool {
        let col = format!("unrecoverable_{seg}");
        self.t.0.column_by_name(&col).is_none() || self.t.bool(&col, i) != Some(true)
    }
    /// The grader's verdict on `{seg}_{part}`; `None` when the read isn't graded for it.
    fn mutation(&self, i: usize, c: &Row, k: usize, part: &str) -> Option<bool> {
        let seg = SEGS[k];
        let tcall = self.t.str(&format!("{seg}_call"), i);
        if !self.recoverable(i, seg) || tcall.is_none() || c.call[k].as_deref() != tcall {
            return None;
        }
        let aa = part.ends_with("_aa");
        let needed: Vec<usize> = match part {
            "insertions" => vec![2],
            "deletions" => vec![3],
            "mutation_frequency" => vec![0, 3],
            _ => vec![aa as usize],
        };
        let (tl, cl) = (self.tlists(i, seg), Self::clists(&c.muts[k]));
        let (tg, cg) = (self.tgerm(i, k), c.germline[k]);
        let derived = part.starts_with("mutation_");
        let reported = |m: &SegmentMutations| match part {
            "mutation_count" => Val::Int(m.count),
            "mutation_count_aa" => m.count_aa.map_or(Val::Float(None), Val::Int),
            "mutation_frequency" => Val::Float(Some(m.frequency)),
            _ => Val::Float(m.frequency_aa),
        };
        let truth_field = || -> Option<Val> {
            let f = format!("{seg}_{part}");
            match part {
                "mutation_count" | "mutation_count_aa" => self.t.int(&f, i).map(Val::Int),
                _ => self.t.float(&f, i).map(|x| Val::Float(Some(x))),
            }
        };
        let (Some(tg), Some(cg)) = (tg, cg) else {
            return Some(tg.is_none() && cg.is_none() && needed.iter().all(|&x| tl[x].is_none() && cl[x].is_none()));
        };
        if needed.iter().any(|&x| tl[x].is_none() || cl[x].is_none()) {
            return Some(needed.iter().all(|&x| tl[x].is_none() && cl[x].is_none()));
        }
        if !self.span_consistent(i, c, k, tg, cg) {
            return Some(false);
        }
        let tspan = self.mkspan(i, seg, tg.0, tg.1);
        let cspan = self.mkspan(i, seg, cg.0, cg.1);
        let shared = self.mkspan(i, seg, tg.0.max(cg.0), tg.1.min(cg.1));
        if !derived {
            return Some(val_eq(&self.val(i, seg, part, &tl, shared), &self.val(i, seg, part, &cl, shared)));
        }
        let rep = reported(c.muts[k].as_ref().unwrap());
        if tspan == cspan {
            let t = truth_field().unwrap_or_else(|| self.val(i, seg, part, &tl, tspan));
            return Some(val_eq(&rep, &t));
        }
        Some(val_eq(&rep, &self.val(i, seg, part, &cl, cspan)) && val_eq(&self.val(i, seg, part, &tl, shared), &self.val(i, seg, part, &cl, shared)))
    }
    /// Read positions (1-based) whose segment the _valid lists (plus the true cut) leave open.
    fn uncertain(&self, i: usize) -> BTreeSet<i64> {
        let mut out = BTreeSet::new();
        for (col, lo, hi) in [("v_sequence_end", 1, 0), ("d_sequence_start", 0, -1), ("d_sequence_end", 1, 0), ("j_sequence_start", 0, -1)] {
            let mut cuts = self.t.ints(&format!("{col}_valid"), i);
            cuts.extend(self.t.int(col, i));
            if let (Some(a), Some(b)) = (cuts.iter().min(), cuts.iter().max()) {
                out.extend(a + lo..=b + hi);
            }
        }
        out
    }
    fn codon_phase(&self, i: usize) -> i64 {
        let m = self.gmap(self.t.str("v_call", i).unwrap());
        match self.t.int("v_germline_start", i) {
            Some(s) if 1 <= s && s as usize <= m.len() => (3 - (m[s as usize - 1] - 1) % 3) % 3,
            _ => 0,
        }
    }
    /// The grader's verdict on mask `col`; `None` when the read isn't graded for it.
    fn mask(&self, i: usize, c: &Row, col: &str) -> Option<bool> {
        let segs: &[usize] = if self.t.str("d_call", i).is_some() { &[0, 2, 1] } else { &[0, 2] };
        if !segs.iter().all(|&k| self.recoverable(i, SEGS[k]) && self.t.str(&format!("{}_call", SEGS[k]), i).is_some() && c.call[k].as_deref() == self.t.str(&format!("{}_call", SEGS[k]), i)) {
            return None;
        }
        let m = c.masks.as_ref().unwrap();
        let get = |col: &str| -> Option<String> {
            match col {
                "cdr_mask" => m.cdr.clone(),
                "cdr_mask_aa" => m.cdr_aa.clone(),
                "gene_segment_mask" => Some(m.segment.clone()),
                "gene_segment_mask_aa" => Some(m.segment_aa.clone()),
                "nongermline_mask" => Some(m.nongermline.clone()),
                _ => Some(m.nongermline_aa.clone()),
            }
        };
        let aa = col.ends_with("_aa");
        let gcol = if aa { "gene_segment_mask_aa" } else { "gene_segment_mask" };
        let (Some(tm), Some(cm)) = (self.t.str(col, i), get(col)) else { return Some(self.t.str(col, i).is_none() && get(col).is_none()) };
        if tm.len() != cm.len() {
            return Some(false);
        }
        let (tg, cg) = (self.t.str(gcol, i).unwrap_or(""), get(gcol).unwrap());
        if tg.len() != tm.len() || cg.len() != cm.len() {
            return Some(false);
        }
        let (tm, cm, tg, cg) = (tm.as_bytes(), cm.as_bytes(), tg.as_bytes(), cg.as_bytes());
        if col == gcol {
            let (amb, start) = (self.uncertain(i), self.t.int("v_sequence_start", i).unwrap());
            let (step, off) = if aa { (3, self.codon_phase(i) + 2) } else { (1, 0) };
            return Some((0..tm.len()).all(|k| tm[k] == cm[k] || amb.contains(&(start + step * k as i64 + off))));
        }
        Some((0..tm.len()).all(|k| tg[k] != cg[k] || tm[k] == cm[k]))
    }
    /// The annotator's mutation fields equal truth's exactly (spans included): the derivation
    /// given truth's cuts and calls. `None` unless the call and both germline spans match truth's.
    fn exact(&self, i: usize, c: &Row, k: usize) -> Option<bool> {
        let seg = SEGS[k];
        let tcall = self.t.str(&format!("{seg}_call"), i);
        let tg = self.tgerm(i, k)?;
        if tcall.is_none() || c.call[k].as_deref() != tcall || c.germline[k] != Some(tg) {
            return None;
        }
        let cuts_true = match k {
            0 => vec!["v_sequence_end"],
            1 => vec!["d_sequence_start", "d_sequence_end"],
            2 => vec!["j_sequence_start"],
            _ => vec![],
        }
        .iter()
        .all(|col| c.cuts.get(col).copied() == self.t.int(col, i));
        if !cuts_true {
            return None;
        }
        let m = c.muts[k].as_ref()?;
        let f = |p: &str| format!("{seg}_{p}");
        Some(
            self.tlists(i, seg) == Self::clists(&c.muts[k])
                && self.t.int(&f("mutation_count"), i) == Some(m.count)
                && self.t.int(&f("mutation_count_aa"), i) == m.count_aa
                && self.t.float(&f("mutation_frequency"), i).is_some_and(|x| (x - m.frequency).abs() < 1e-9)
                && match (self.t.float(&f("mutation_frequency_aa"), i), m.frequency_aa) {
                    (Some(x), Some(y)) => (x - y).abs() < 1e-9,
                    (x, y) => x.is_none() && y.is_none(),
                },
        )
    }
}

impl Grader<'_> {
    /// All six masks equal truth's exactly, on reads whose calls, V start, cuts and J end are
    /// truth's: the derivation given truth's cuts and calls. `None` on other reads.
    fn masks_exact(&self, i: usize, c: &Row) -> Option<bool> {
        let t = self.t;
        let calls = (0..3).all(|k| c.call[k].as_deref() == t.str(&format!("{}_call", SEGS[k]), i));
        let cuts = ["v_sequence_end", "d_sequence_start", "d_sequence_end", "j_sequence_start"].iter().all(|x| c.cuts.get(x).copied() == t.int(x, i));
        let ends = Some(c.v_start) == t.int("v_sequence_start", i) && c.germline[2].map(|g| g.1) == t.int("j_germline_end", i);
        if !(calls && cuts && ends) {
            return None;
        }
        let m = c.masks.as_ref()?;
        let ours = [m.cdr.as_deref(), Some(m.segment.as_str()), Some(m.nongermline.as_str()), m.cdr_aa.as_deref(), Some(m.segment_aa.as_str()), Some(m.nongermline_aa.as_str())];
        Some(MASKS.iter().zip(ours).all(|(col, x)| t.str(col, i) == x))
    }
}

/// (correct, graded) per field.
type Tally = BTreeMap<String, (usize, usize)>;

fn bump(t: &mut Tally, k: String, ok: Option<bool>) {
    if let Some(ok) = ok {
        let e = t.entry(k).or_default();
        e.0 += ok as usize;
        e.1 += 1;
    }
}

/// Grades every mutation field and mask of `g` (graded reads only), printed; returns the tally
/// and up to 10 V/J list misses.
fn grade(label: &str, g: &Generated) -> (Tally, Vec<String>) {
    let r = human();
    let rows: Vec<Option<Row>> = (0..g.len()).into_par_iter().map(|i| call(r, g.sequence(i).as_bytes()).map(|c| row(r, &c))).collect();
    let gr = Grader { r, t: g };
    let (mut t, mut missed) = (Tally::new(), Vec::new());
    let empty = Row::default();
    for (i, c) in rows.iter().enumerate() {
        let c = c.as_ref().unwrap_or(&empty);
        for (k, seg) in SEGS.iter().enumerate() {
            for part in PARTS {
                let ok = gr.mutation(i, c, k, part);
                if ok == Some(false) && part == "mutations" && k != 1 && missed.len() < 10 {
                    missed.push(format!(
                        "{} {seg}: truth {:?} {:?} annotator {:?} {:?}",
                        g.str("sequence_id", i).unwrap(),
                        gr.tgerm(i, k),
                        g.opt_list(&format!("{seg}_mutations"), i),
                        c.germline[k],
                        c.muts[k].as_ref().map(|m| &m.nt)
                    ));
                }
                bump(&mut t, format!("{seg}_{part}"), ok);
            }
            bump(&mut t, format!("{seg} exact given truth's cuts"), gr.exact(i, c, k));
        }
        if c.masks.is_some() {
            for col in MASKS {
                bump(&mut t, col.to_string(), gr.mask(i, c, col));
            }
            bump(&mut t, "masks exact given truth's cuts".to_string(), gr.masks_exact(i, c));
        }
    }
    println!("mutations {label} ({} reads):", g.len());
    for (k, (ok, n)) in &t {
        println!("  {k:<40} {ok:>5}/{n:<5} {:.2}%", 100.0 * *ok as f64 / (*n).max(1) as f64);
    }
    (t, missed)
}

fn rate(t: &Tally, k: &str) -> f64 {
    let (ok, n) = t[k];
    ok as f64 / n as f64
}

#[test]
fn mutation_accuracy_memory() {
    // the memory preset with a 40 nt constant region, so C's lists are graded too
    let g = generated(12, "memory", 2000, |p| p.constant_len = Some(40));
    let (t, missed) = grade("memory + 40 nt C, seed 12", &g);
    // C is graded only where the C primary is the true allele; 40 nt rarely separate C alleles
    for (seg, min) in [("v", 1500), ("d", 1000), ("j", 1500), ("c", 150)] {
        let k = format!("{seg}_mutations");
        assert!(t[&k].1 > min, "{k} graded on {} reads", t[&k].1);
        assert!(rate(&t, &k) >= 0.99, "{k} {:?}: {missed:#?}", t[&k]);
    }
}

#[test]
fn mutation_accuracy_report() {
    // Printed only (not asserted).
    let sets = [
        ("naive, seed 1", generated(1, "naive", 2000, |_| {})),
        ("mutation 0-30%, seed 1", generated(1, "naive", 2000, |p| p.mutation = Spec::Uniform(0.0, 0.3))),
        ("memory + 20% indels, seed 11", generated(11, "memory", 3000, |p| p.indel_rate = 0.2)),
        ("bnab, seed 1", generated(1, "bnab", 2000, |_| {})),
    ];
    for (label, g) in &sets {
        let (_, missed) = grade(label, g);
        for m in missed.iter().take(3) {
            println!("  {m}");
        }
    }
}
