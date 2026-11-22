//! IMGT regions, junction, translation, alignment strings, productivity and assembled
//! sequences (record 0032) on the built-in human database.
mod common;
use abstar_annotate::align::Scoring;
use abstar_annotate::assign::{assign_c, assign_j, resolve_junction, CHit, JHit, Junction, Oriented, VHit};
use abstar_annotate::derive::{alignments, productivity, region_aa, regions, translate, Productivity, Regions};
use abstar_generate::Spec;
use common::{generated, human, rand_seq, Generated, Ref};
use rayon::prelude::*;

const IGHV: &str = "IGHV1-2*02";
const N1: &[u8] = b"CC";
const IGHD: &[u8] = b"GTATTACGATTTTTGGAGTGGTTATTATACC"; // IGHD3-3*01
const N2: &[u8] = b"GG";

fn ungapped(r: &Ref, name: &str) -> Vec<u8> {
    r.alleles[r.idx(name)].sequence().into_bytes()
}

/// IGHV1-2*02 + CC + IGHD3-3*01 + GG + IGHJ4*02: in frame, productive.
fn igh_read(r: &Ref) -> Vec<u8> {
    [&ungapped(r, IGHV)[..], N1, IGHD, N2, &ungapped(r, "IGHJ4*02")].concat()
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
    fn regions(&self, r: &Ref) -> Regions {
        regions(&self.read, &self.v, self.j(), self.c.as_ref(), &r.alleles)
    }
}

fn call(r: &Ref, read: &[u8]) -> Option<Called> {
    let (o, v) = r.call(read)?;
    let mut j = assign_j(&o, &v, &r.alleles, &r.j, &Scoring::default());
    let c = j.as_mut().and_then(|j| assign_c(&o, j, &r.alleles, &r.c, &Scoring::default()));
    let j = j.map(|j| {
        let has_d = r.d.loci.contains(&r.alleles[v.primary].locus);
        let jn = resolve_junction(&o, &v, &j, &r.alleles, &r.d, has_d);
        (j, jn)
    });
    Some(Called { read: o, v, j, c })
}

fn forward(r: &Ref, read: &[u8]) -> (Called, Regions) {
    let c = call(r, read).expect("V hit");
    assert!(!c.read.rev_comp);
    let g = c.regions(r);
    (c, g)
}

fn text(seq: &[u8], x: Option<(usize, usize)>) -> Option<String> {
    x.map(|(a, b)| String::from_utf8(seq[a..b].to_vec()).unwrap())
}

#[test]
fn regions_hand_built_igh() {
    let r = human();
    let read = igh_read(r);
    let (c, g) = forward(r, &read);
    assert_eq!(r.alleles[c.v.primary].name, IGHV);
    let expected = [(0, 75), (75, 99), (99, 150), (150, 174), (174, 288), (288, 345), (345, 379)];
    assert_eq!(g.all().map(|x| x.1), expected.map(Some));
    assert_eq!(g.junction, Some((285, 348)));
    assert_eq!((g.frame, g.vdj, g.vdjc), (0, (0, 379), None));
    let aa: Vec<String> = g.all().iter().map(|x| region_aa(&read, g.frame, x.1.unwrap())).collect();
    assert_eq!(aa, ["QVQLVQSGAEVKKPGASVKVSCKAS", "GYTFTGYY", "MHWVRQAPGQGLEWMGW", "INPNSGGT", "NYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC", "ARDRITIFGVVIIPDYFDY", "WGQGTLVTVSS"]);
    assert_eq!(text(&read, g.cdr3).unwrap(), "GCGAGAGACCGTATTACGATTTTTGGAGTGGTTATTATACCGGACTACTTTGACTAC");
    assert_eq!(text(&read, g.junction).unwrap(), "TGTGCGAGAGACCGTATTACGATTTTTGGAGTGGTTATTATACCGGACTACTTTGACTACTGG");
    assert_eq!(translate(&read[285..348], 0), "CARDRITIFGVVIIPDYFDYW");
}

#[test]
fn regions_igk() {
    let r = human();
    let read = [&ungapped(r, "IGKV3-20*01")[..], &ungapped(r, "IGKJ1*01")].concat();
    let (_, g) = forward(r, &read);
    let expected = [(0, 78), (78, 99), (99, 150), (150, 159), (159, 267), (267, 297), (297, 328)];
    assert_eq!(g.all().map(|x| x.1), expected.map(Some));
    assert_eq!(g.junction, Some((264, 300)));
    assert_eq!(translate(&read[264..300], 0), "CQQYGSSPPWTF");
    assert_eq!(region_aa(&read, g.frame, g.cdr3.unwrap()), "QQYGSSPPWT");
    assert_eq!(region_aa(&read, g.frame, g.fwr4.unwrap()), "FGQGTKVEIK");
}

#[test]
fn truncated_into_cdr1() {
    // Review Focus 4: the read starts at IGHV1-2*02 base 85 (IMGT 89, in CDR1).
    let r = human();
    let read = igh_read(r)[85..].to_vec();
    let (_, g) = forward(r, &read);
    assert_eq!(g.fwr1, None);
    let expected = [(0, 14), (14, 65), (65, 89), (89, 203), (203, 260), (260, 294)];
    assert_eq!([g.cdr1, g.fwr2, g.cdr2, g.fwr3, g.cdr3, g.fwr4], expected.map(Some));
    assert_eq!(g.junction, Some((200, 263)));
    // IMGT 89 is the 2nd base of codon 30, so the first full codon starts at read 2.
    assert_eq!(g.frame, 2);
    assert_eq!(region_aa(&read, g.frame, g.cdr1.unwrap()), "TGYY");
    assert_eq!(region_aa(&read, g.frame, g.fwr3.unwrap()), "NYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC");
}

#[test]
fn truncated_into_cdr3() {
    // The read ends 10 nt past V, in bases no J allele matches at the Wave 3 low-path floor: no J
    // hit (status `no_j_hit`, record 0035), no junction, CDR3
    // the covered part (ruling A4), sequence_vdj through the read end (ruling A8).
    let r = human();
    let v = ungapped(r, IGHV);
    let read = [&v[..], b"AACCAACCAA"].concat();
    let (c, g) = forward(r, &read);
    assert!(c.j.is_none());
    assert_eq!(g.junction, None);
    assert_eq!((g.fwr3, g.cdr3, g.fwr4), (Some((174, 288)), Some((288, read.len())), None));
    assert_eq!(g.vdj, (0, read.len()));
    let p = productivity(&c.read, &g);
    assert_eq!(p, Productivity { productive: false, vj_in_frame: false, stop_codon: false });
}

#[test]
fn large_indel_regions() {
    // Review Focus 5: 27 nt inserted after IGHV1-2*02 base 87 (CDR1) and bases 200..236 (FWR3)
    // deleted leave the IMGT regions and the junction in place.
    let r = human();
    let v = ungapped(r, IGHV);
    let ins = (0..).map(|seed| rand_seq(27, seed)).find(|x| !translate(x, 0).contains('*')).unwrap();
    let vread = [&v[..87], &ins, &v[87..200], &v[236..]].concat();
    let read = [&vread[..], N1, IGHD, N2, &ungapped(r, "IGHJ4*02")].concat();
    let (c, g) = forward(r, &read);
    assert_eq!(r.alleles[c.v.primary].gene, "IGHV1-2");
    let expected = [(0, 75), (75, 126), (126, 177), (177, 201), (201, 279), (279, 336), (336, 370)];
    assert_eq!(g.all().map(|x| x.1), expected.map(Some));
    assert_eq!(g.junction, Some((276, 339)));
    assert_eq!(translate(&read[276..339], 0), "CARDRITIFGVVIIPDYFDYW");
    let p = productivity(&c.read, &g);
    assert!(p.productive && p.vj_in_frame && !p.stop_codon, "{p:?}");
    let a = alignments(&c.read, &c.v, c.j(), &r.alleles, g.frame);
    assert_eq!(a.sequence.len(), a.germline.len());
    assert_eq!(a.sequence.bytes().filter(|&b| b == b'-').count(), 36);
    assert_eq!(a.germline.bytes().filter(|&b| b == b'-').count(), 27 + N1.len() + N2.len());
}

#[test]
fn alignment_strings() {
    let r = human();
    let read = igh_read(r);
    let (c, g) = forward(r, &read);
    let a = alignments(&c.read, &c.v, c.j(), &r.alleles, g.frame);
    let vg = &r.alleles[r.idx(IGHV)].sequence_gapped;
    // IMGT gaps in both strings over V; N bases align to germline gaps; D and J align as called.
    assert!(a.sequence.starts_with(vg.as_str()) && a.germline.starts_with(vg.as_str()));
    let (vn, jl) = (vg.len(), ungapped(r, "IGHJ4*02").len());
    let jn = c.j.as_ref().unwrap().1.clone();
    let d = jn.d.as_ref().expect("D call");
    let tail_s = &a.sequence[vn..];
    let tail_g = &a.germline[vn..];
    assert_eq!(tail_s.as_bytes(), &read[jn.v_end..]);
    let n1 = d.start - jn.v_end;
    let n2 = jn.j_start - d.end;
    let dg = &r.alleles[d.primary].sequence()[d.g_start..d.g_end];
    let jg = &ungapped(r, "IGHJ4*02")[jl - (read.len() - jn.j_start)..];
    let germ = ["-".repeat(n1), dg.to_string(), "-".repeat(n2), String::from_utf8(jg.to_vec()).unwrap()].concat();
    assert_eq!(tail_g, germ);
    assert_eq!(a.sequence.replace('.', "").as_bytes(), &read[..]);
    assert_eq!(a.sequence_aa.replace('.', ""), translate(&read, 0));
    assert_eq!(a.germline_aa.len(), a.sequence_aa.len());
    assert!(a.sequence_aa.contains("...."), "{}", a.sequence_aa);
    // N columns: a codon with germline gaps only is '-', one mixing gaps and bases is 'X'.
    assert!(a.germline_aa.contains('X'), "{}", a.germline_aa);
}

#[test]
fn productivity_cases() {
    let r = human();
    let base = igh_read(r);
    let p = |read: &[u8]| {
        let (c, g) = forward(r, read);
        productivity(&c.read, &g)
    };
    let yes = Productivity { productive: true, vj_in_frame: true, stop_codon: false };
    assert_eq!(p(&base), yes);
    // one more N base: out of frame
    let shifted = [&ungapped(r, IGHV)[..], b"CCA", IGHD, N2, &ungapped(r, "IGHJ4*02")].concat();
    let x = p(&shifted);
    assert!(!x.productive && !x.vj_in_frame, "{x:?}");
    // FWR3 Gln (CAG at 183) to TAG
    let mut stop = base.clone();
    stop[183] = b'T';
    assert_eq!(p(&stop), Productivity { productive: false, vj_in_frame: true, stop_codon: true });
    // Cys104 TGT to TCT
    let mut anchor = base.clone();
    anchor[286] = b'C';
    assert_eq!(p(&anchor), Productivity { productive: false, vj_in_frame: true, stop_codon: false });
}

#[test]
fn translate_rules() {
    assert_eq!(translate(b"ATGGCC", 0), "MA");
    assert_eq!(translate(b"ATGGC", 0), "M");
    assert_eq!(translate(b"CATGGCC", 1), "MA");
    assert_eq!(translate(b"ATNTAATGA", 0), "X**");
    assert_eq!(translate(b"atg", 0), "M");
    assert_eq!(translate(b"AT", 0), "");
    assert_eq!(translate(b"AT", 5), "");
}

const NAMES: [&str; 7] = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"];

/// Per-field agreement with truth over the graded reads.
#[derive(Default)]
struct Tally {
    n: usize,
    regions: [usize; 7],
    all_regions: usize,
    /// (all regions right, graded): reads whose V start, V end and J start are truth's.
    true_cuts: (usize, usize),
    junction: usize,
    junction_aa: usize,
    productive: usize,
    vj_in_frame: usize,
    stop_codon: usize,
    vdj: usize,
    vdjc: usize,
    v_start: usize,
}

/// (nt, aa, 1-based start, 1-based end) of a region; empty regions have no start or end.
type Fields = (Option<String>, Option<String>, Option<i64>, Option<i64>);

fn annotated(seq: &[u8], frame: usize, x: Option<(usize, usize)>) -> Fields {
    match x {
        None => (None, None, None, None),
        Some((a, b)) => {
            let pos = (a < b).then_some((a as i64 + 1, b as i64));
            (text(seq, x), Some(region_aa(seq, frame, (a, b))), pos.map(|p| p.0), pos.map(|p| p.1))
        }
    }
}

fn truth(g: &Generated, i: usize, k: &str) -> Fields {
    (
        g.str(k, i).map(String::from),
        g.str(&format!("{k}_aa"), i).map(String::from),
        g.int(&format!("{k}_start"), i),
        g.int(&format!("{k}_end"), i),
    )
}

/// Grades the reads of `g` whose V and J calls are the true alleles; printed. Returns the tally and up to
/// 10 reads with a wrong region.
fn grade(label: &str, g: &Generated) -> (Tally, Vec<String>) {
    let r = human();
    let calls: Vec<Option<Called>> = (0..g.len()).into_par_iter().map(|i| call(r, g.sequence(i).as_bytes())).collect();
    let (mut t, mut missed) = (Tally::default(), Vec::new());
    for (i, c) in calls.iter().enumerate() {
        // the grader's condition for region fields: the true V allele, and the true J allele (or
        // no J where the read has no J base)
        let Some(c) = c.as_ref().filter(|c| {
            let j_ok = match &c.j {
                Some((j, _)) => Some(r.alleles[j.primary].name.as_str()) == g.str("j_call", i),
                None => g.list("j_call_valid", i).is_empty(),
            };
            j_ok && Some(r.alleles[c.v.primary].name.as_str()) == g.str("v_call", i)
        }) else {
            continue;
        };
        t.n += 1;
        let seq = &c.read.seq;
        let rg = c.regions(r);
        let (v_end, j_start) = c.j.as_ref().map_or((c.v.aln.q_end, None), |(_, jn)| (jn.v_end, Some(jn.j_start as i64 + 1)));
        let cuts = [Some(c.v.aln.q_start as i64 + 1), Some(v_end as i64), j_start];
        let cuts_true = cuts == [g.int("v_sequence_start", i), g.int("v_sequence_end", i), g.int("j_sequence_start", i)];
        let mut all = true;
        for (k, (name, x)) in rg.all().iter().enumerate() {
            let ok = annotated(seq, rg.frame, *x) == truth(g, i, name);
            t.regions[k] += ok as usize;
            all &= ok;
            if !ok && missed.len() < 10 && cuts_true {
                missed.push(format!("{} {name}: {:?} truth {:?}", g.str("sequence_id", i).unwrap(), annotated(seq, rg.frame, *x), truth(g, i, name)));
            }
        }
        t.all_regions += all as usize;
        if cuts_true {
            t.true_cuts.1 += 1;
            t.true_cuts.0 += all as usize;
        }
        let junction = text(seq, rg.junction);
        t.junction += (junction.as_deref() == g.str("junction", i)) as usize;
        t.junction_aa += (rg.junction.map(|(a, b)| translate(&seq[a..b], 0)).as_deref() == g.str("junction_aa", i)) as usize;
        let p = productivity(&c.read, &rg);
        t.productive += (Some(p.productive) == g.bool("productive", i)) as usize;
        t.vj_in_frame += (Some(p.vj_in_frame) == g.bool("vj_in_frame", i)) as usize;
        t.stop_codon += (Some(p.stop_codon) == g.bool("stop_codon", i)) as usize;
        t.vdj += (text(seq, Some(rg.vdj)).as_deref() == g.str("sequence_vdj", i)) as usize;
        t.vdjc += (text(seq, rg.vdjc).as_deref() == g.str("sequence_vdjc", i)) as usize;
        t.v_start += (Some(c.v.aln.q_start as i64 + 1) == g.int("v_sequence_start", i)) as usize;
    }
    let pct = |x: usize| 100.0 * x as f64 / t.n.max(1) as f64;
    let per: Vec<String> = NAMES.iter().zip(t.regions).map(|(n, x)| format!("{n} {:.2}%", pct(x))).collect();
    println!(
        "derive {label}: {} of {} reads with the true V and J alleles; {}; all regions {:.2}% (with truth's V start, V end and J start {}/{} = {:.2}%); junction {:.2}% junction_aa {:.2}%; \
         productive {:.2}% vj_in_frame {:.2}% stop_codon {:.2}%; sequence_vdj {:.2}% sequence_vdjc {:.2}%; v_sequence_start {:.2}%",
        t.n,
        g.len(),
        per.join(" "),
        pct(t.all_regions),
        t.true_cuts.0,
        t.true_cuts.1,
        100.0 * t.true_cuts.0 as f64 / t.true_cuts.1.max(1) as f64,
        pct(t.junction),
        pct(t.junction_aa),
        pct(t.productive),
        pct(t.vj_in_frame),
        pct(t.stop_codon),
        pct(t.vdj),
        pct(t.vdjc),
        pct(t.v_start),
    );
    (t, missed)
}

#[test]
fn region_accuracy_naive() {
    // Truth places the regions from each base's origin. Where a chance N match lets the read
    // support another V end or J start (record 0025's valid cuts) and that cut decides a region
    // edge or whether an anchor base is in its segment, the annotator's (valid) cut can give
    // other regions than truth's; with truth's cuts the derivation must match exactly (record
    // 0032). The overall rate (about 98.3%) is asserted at 98%.
    let (t, missed) = grade("naive, seed 1", &generated(1, "naive", 2000, |_| {}));
    let rate = t.all_regions as f64 / t.n as f64;
    let exact = t.true_cuts.0 as f64 / t.true_cuts.1 as f64;
    assert!(t.n > 1900 && t.true_cuts.1 > 800, "graded {} with truth's cuts {}", t.n, t.true_cuts.1);
    assert!(exact >= 0.995, "all regions with truth's cuts {exact}: {missed:#?}");
    assert!(rate >= 0.98, "all regions {rate} over {}", t.n);
}

#[test]
fn derive_accuracy_report() {
    // Printed only (not asserted).
    let sets = [
        ("mutation 0-30%, seed 1", generated(1, "naive", 2000, |p| p.mutation = Spec::Uniform(0.0, 0.3))),
        ("memory + 20% indels, seed 11", generated(11, "memory", 3000, |p| p.indel_rate = 0.2)),
        ("bnab, seed 1", generated(1, "bnab", 2000, |_| {})),
        ("naive + 40 nt C, seed 2", generated(2, "naive", 2000, |p| p.constant_len = Some(40))),
        (
            "naive, truncated 0-120 nt 5' and 0-60 nt 3', seed 4",
            generated(4, "naive", 2000, |p| {
                p.truncate_5p = Spec::Uniform(0.0, 120.0);
                p.truncate_3p = Spec::Uniform(0.0, 60.0);
            }),
        ),
    ];
    for (label, g) in &sets {
        let (_, missed) = grade(label, g);
        for m in missed.iter().take(3) {
            println!("  {m}");
        }
    }
}

