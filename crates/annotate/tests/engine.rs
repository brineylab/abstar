//! The assembled annotator end to end (record 0035): `Engine::annotate` on hand-built reads with
//! exact AIRR values, the status rules, the built-in mixed-species databases, and
//! determinism through the run pipeline.
mod common;
use abstar_annotate::engine::{Annotate, Engine, EngineOptions};
use abstar_annotate::input::Input;
use abstar_annotate::pipeline::{run_to_writer, RunOptions};
use abstar_germline::{Database, Locus, Segment};
use abstar_generate::Spec;
use common::{generated, human_db, package_root, rand_seq, revcomp};
use std::sync::OnceLock;

fn load(name: &str) -> Database {
    abstar_germline::load(&package_root().join("germline_dbs").join(name)).unwrap()
}

fn human() -> &'static Database {
    static D: OnceLock<Database> = OnceLock::new();
    D.get_or_init(human_db)
}

fn engine() -> &'static Engine<'static> {
    static E: OnceLock<Engine<'static>> = OnceLock::new();
    E.get_or_init(|| Engine::new(human(), &EngineOptions::default()).unwrap())
}

fn ungapped(db: &Database, name: &str) -> String {
    db.alleles.iter().find(|a| a.name == name).unwrap_or_else(|| panic!("no allele {name}")).sequence()
}

/// Ungapped index of IMGT position `p` (1-based) of a V allele.
fn at_imgt(db: &Database, v: &str, p: usize) -> usize {
    let a = db.alleles.iter().find(|a| a.name == v).unwrap();
    let g = a.sequence_gapped.as_bytes();
    assert_ne!(g[p - 1], b'.');
    g[..p - 1].iter().filter(|&&b| b != b'.').count()
}

/// IGHV1-2*02 (A>G at IMGT 152) + CC + IGHD3-3*01 + GG + IGHJ4*02: in frame, productive.
fn igh_read() -> String {
    let db = human();
    let mut v = ungapped(db, "IGHV1-2*02").into_bytes();
    v[at_imgt(db, "IGHV1-2*02", 152)] = b'G';
    [String::from_utf8(v).unwrap(), "CC".into(), ungapped(db, "IGHD3-3*01"), "GG".into(), ungapped(db, "IGHJ4*02")].concat()
}

fn s(x: &str) -> Option<String> {
    Some(x.to_string())
}

fn strs(x: &[&str]) -> Option<Vec<String>> {
    Some(x.iter().map(|s| s.to_string()).collect())
}

#[test]
fn engine_end_to_end_igh() {
    let read = igh_read();
    assert_eq!(read.len(), 379);
    let r = engine().annotate("r1", &read);
    assert_eq!((r.annotation_status.as_str(), r.failure_reason.as_deref()), ("annotated", None));
    assert_eq!((r.sequence_id.as_str(), r.sequence.as_str()), ("r1", read.as_str()));
    assert_eq!((r.rev_comp, r.locus.as_deref()), (Some(false), Some("IGH")));
    // calls
    assert_eq!((r.v_call.as_deref(), r.v_gene.as_deref(), r.v_species.as_deref()), (Some("IGHV1-2*02"), Some("IGHV1-2"), Some("human")));
    assert_eq!((r.d_call.as_deref(), r.d_gene.as_deref(), r.d_inverted), (Some("IGHD3-3*01"), Some("IGHD3-3"), Some(false)));
    assert_eq!((r.j_call.as_deref(), r.j_gene.as_deref()), (Some("IGHJ4*02"), Some("IGHJ4")));
    assert_eq!(r.c_call, None);
    for (calls, species) in [(&r.v_call_ties, &r.v_species_ties), (&r.d_call_ties, &r.d_species_ties), (&r.j_call_ties, &r.j_species_ties)] {
        let (calls, species) = (calls.as_ref().unwrap(), species.as_ref().unwrap());
        assert_eq!(calls.len(), species.len());
        assert!(species.iter().all(|s| s == "human"));
    }
    assert_eq!(r.v_call_ties.as_ref().unwrap()[0], "IGHV1-2*02");
    assert_eq!(r.j_call_ties.as_ref().unwrap()[0], "IGHJ4*02");
    // coordinates: AIRR 1-based closed
    assert_eq!((r.v_sequence_start, r.v_sequence_end, r.v_germline_start, r.v_germline_end), (Some(1), Some(296), Some(1), Some(296)));
    assert_eq!((r.d_sequence_start, r.d_sequence_end, r.d_germline_start, r.d_germline_end), (Some(299), Some(329), Some(1), Some(31)));
    assert_eq!((r.j_sequence_start, r.j_sequence_end, r.j_germline_start, r.j_germline_end), (Some(332), Some(379), Some(1), Some(48)));
    assert_eq!((r.v_cigar.as_deref(), r.d_cigar.as_deref(), r.j_cigar.as_deref()), (Some("296M83S"), Some("298S31M50S"), Some("331S48M")));
    // scores: record 0029 scoring over the reported columns; identity = matches / aligned
    assert_eq!((r.v_score, r.v_identity), (Some(295.0 * 2.0 - 3.0), Some(295.0 / 296.0)));
    assert_eq!((r.d_score, r.d_identity, r.j_score, r.j_identity), (Some(62.0), Some(1.0), Some(96.0), Some(1.0)));
    assert_eq!((r.np1.as_deref(), r.np1_length, r.np2.as_deref(), r.np2_length), (Some("CC"), Some(2), Some("GG"), Some(2)));
    // junction and regions
    assert_eq!(r.junction, s("TGTGCGAGAGACCGTATTACGATTTTTGGAGTGGTTATTATACCGGACTACTTTGACTACTGG"));
    assert_eq!((r.junction_aa.as_deref(), r.junction_length), (Some("CARDRITIFGVVIIPDYFDYW"), Some(63)));
    let regions = [
        (&r.fwr1_start, &r.fwr1_end, &r.fwr1_aa, &r.fwr1_length, (1, 75), "QVQLVQSGAEVKKPGASVKVSCKAS"),
        (&r.cdr1_start, &r.cdr1_end, &r.cdr1_aa, &r.cdr1_length, (76, 99), "GYTFTGYY"),
        (&r.fwr2_start, &r.fwr2_end, &r.fwr2_aa, &r.fwr2_length, (100, 150), "MHWVRQAPGQGLGWMGW"),
        (&r.cdr2_start, &r.cdr2_end, &r.cdr2_aa, &r.cdr2_length, (151, 174), "INPNSGGT"),
        (&r.fwr3_start, &r.fwr3_end, &r.fwr3_aa, &r.fwr3_length, (175, 288), "NYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC"),
        (&r.cdr3_start, &r.cdr3_end, &r.cdr3_aa, &r.cdr3_length, (289, 345), "ARDRITIFGVVIIPDYFDY"),
        (&r.fwr4_start, &r.fwr4_end, &r.fwr4_aa, &r.fwr4_length, (346, 379), "WGQGTLVTVSS"),
    ];
    for (start, end, aa, len, (a, b), want) in regions {
        assert_eq!((*start, *end, aa.as_deref(), *len), (Some(a), Some(b), Some(want), Some(want.len() as i64)));
    }
    assert_eq!(r.cdr3, s(&read[288..345]));
    assert_eq!(r.fwr1, s(&read[..75]));
    // productivity and assembled sequences
    assert_eq!((r.productive, r.vj_in_frame, r.stop_codon), (Some(true), Some(true), Some(false)));
    assert_eq!(r.sequence_vdj, s(&read));
    assert_eq!((r.sequence_vdjc.as_deref(), r.sequence_vdjc_aa.as_deref()), (None, None));
    let aa = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLGWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARDRITIFGVVIIPDYFDYWGQGTLVTVSS";
    assert_eq!((r.sequence_aa.as_deref(), r.sequence_vdj_aa.as_deref()), (Some(aa), Some(aa)));
    assert_eq!(r.sequence_alignment_aa.as_deref().map(|x| x.replace('.', "")), Some(aa.to_string()));
    assert_eq!(r.sequence_alignment.as_deref().map(|x| x.replace(['.', '-'], "")), Some(read.clone()));
    // mutations
    assert_eq!((r.v_mutations.clone(), r.v_mutations_aa.clone()), (strs(&["IMGT152:A>G"]), strs(&["IMGT51:E>G"])));
    assert_eq!((r.v_mutation_count, r.v_mutation_count_aa, r.v_mutation_frequency), (Some(1), Some(1), Some(1.0 / 296.0)));
    assert_eq!((r.v_insertions.clone(), r.v_deletions.clone()), (strs(&[]), strs(&[])));
    assert_eq!((r.d_mutations.clone(), r.d_mutation_count, r.j_mutations.clone(), r.j_mutation_count), (strs(&[]), Some(0), strs(&[]), Some(0)));
    assert_eq!((r.c_mutations.clone(), r.c_mutation_count), (None, None));
    // masks
    let cdr = ["0".repeat(75), "1".repeat(24), "0".repeat(51), "2".repeat(24), "0".repeat(114), "3".repeat(57), "0".repeat(34)].concat();
    assert_eq!(r.cdr_mask, Some(cdr));
    assert_eq!(r.gene_segment_mask, Some(["V".repeat(296), "NN".into(), "D".repeat(31), "NN".into(), "J".repeat(48)].concat()));
    let mut ng = vec![b'0'; 379];
    for i in [at_imgt(human(), "IGHV1-2*02", 152), 296, 297, 329, 330] {
        ng[i] = b'1';
    }
    assert_eq!(r.nongermline_mask, Some(String::from_utf8(ng).unwrap()));
    assert_eq!(r.cdr_mask_aa.as_ref().map(|m| m.len()), Some(aa.len()));
    // a read and its reverse complement annotate identically apart from `sequence` and `rev_comp`
    let mut rc = engine().annotate("r1", &String::from_utf8(revcomp(read.as_bytes())).unwrap());
    assert_eq!(rc.rev_comp, Some(true));
    (rc.sequence, rc.rev_comp) = (r.sequence.clone(), r.rev_comp);
    assert_eq!(rc, r);
}

/// `igh_read` + IGHG1*01 2-61, as real transcripts show it (the C allele's base 1 is the J's last
/// base), with one C substitution: the C fields in AIRR coordinates, the C codons in the J frame.
#[test]
fn engine_c_fields() {
    let db = human();
    let g1 = ungapped(db, "IGHG1*01");
    assert_eq!(&g1[..12], "GCCTCCACCAAG");
    let mut tail = g1[1..61].to_string().into_bytes();
    tail[3] = b'A'; // allele base 5: C>A, codon 2 TCC>TAC (S>Y)
    let read = igh_read() + std::str::from_utf8(&tail).unwrap();
    let r = engine().annotate("r1", &read);
    assert_eq!(r.c_gene.as_deref(), Some("IGHG1"));
    assert!(r.c_call_ties.as_ref().unwrap().contains(&"IGHG1*01".to_string()));
    // J keeps its last base; C is read 380-439 = allele 2-61
    assert_eq!((r.j_sequence_end, r.j_germline_end), (Some(379), Some(48)));
    assert_eq!((r.c_sequence_start, r.c_sequence_end, r.c_germline_start, r.c_germline_end), (Some(380), Some(439), Some(2), Some(61)));
    assert_eq!(r.c_identity, Some(59.0 / 60.0));
    assert_eq!((r.c_mutations.clone(), r.c_mutations_aa.clone()), (strs(&["5:C>A"]), strs(&["2:S>Y"])));
    assert_eq!((r.c_mutation_count, r.c_mutation_count_aa), (Some(1), Some(1)));
    assert_eq!((r.c_mutation_frequency, r.c_mutation_frequency_aa), (Some(1.0 / 60.0), Some(1.0 / 19.0)));
    assert_eq!(r.sequence_vdjc, s(&read));
    assert!(r.sequence_vdjc_aa.as_deref().unwrap().ends_with("WGQGTLVTVSSAYTKGPSVFPLAPSSKSTSG"), "{:?}", r.sequence_vdjc_aa);
    // the masks stop at J's end
    assert_eq!(r.gene_segment_mask.as_ref().map(|m| m.len()), Some(379));
}

/// An ambiguous junction base makes the read nonproductive (record 0032), though it stays in frame
/// with no stop codon.
#[test]
fn ambiguous_junction_base_is_nonproductive() {
    let mut read = igh_read().into_bytes();
    read[310] = b'N'; // inside the D, in CDR3
    let r = engine().annotate("r1", std::str::from_utf8(&read).unwrap());
    assert_eq!(r.junction.as_deref().map(|j| j.contains('N')), Some(true));
    assert_eq!((r.productive, r.vj_in_frame, r.stop_codon), (Some(false), Some(true), Some(false)));
    // an N outside the junction leaves it productive
    let mut read = igh_read().into_bytes();
    read[100] = b'N';
    let r = engine().annotate("r1", std::str::from_utf8(&read).unwrap());
    assert_eq!((r.productive, r.vj_in_frame, r.stop_codon), (Some(true), Some(true), Some(false)));
}

#[test]
fn no_j_hit_keeps_v_fields() {
    let read = [ungapped(human(), "IGHV1-2*02"), "CCGTATTACG".into()].concat();
    let r = engine().annotate("r", &read);
    assert_eq!((r.annotation_status.as_str(), r.failure_reason.as_deref()), ("unassigned", Some("no_j_hit")));
    assert_eq!((r.v_call.as_deref(), r.v_species.as_deref(), r.rev_comp), (Some("IGHV1-2*02"), Some("human"), Some(false)));
    assert!(r.v_call_ties.as_ref().is_some_and(|t| t[0] == "IGHV1-2*02"));
    assert_eq!((r.v_sequence_start, r.v_sequence_end, r.v_germline_start, r.v_germline_end), (Some(1), Some(296), Some(1), Some(296)));
    let got = [(r.fwr1_start, r.fwr1_end), (r.cdr1_start, r.cdr1_end), (r.fwr2_start, r.fwr2_end), (r.cdr2_start, r.cdr2_end), (r.fwr3_start, r.fwr3_end)];
    assert_eq!(got, [(1, 75), (76, 99), (100, 150), (151, 174), (175, 288)].map(|(a, b)| (Some(a), Some(b))));
    assert_eq!(r.fwr3_aa.as_deref(), Some("NYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC"));
    assert_eq!((r.v_mutations.clone(), r.v_mutation_count), (strs(&[]), Some(0)));
    assert_eq!((r.j_call.clone(), r.d_call.clone(), r.junction.clone(), r.junction_length), (None, None, None, None));
    assert_eq!((r.j_sequence_start, r.j_mutations.clone(), r.np1.clone(), r.fwr4.clone()), (None, None, None, None));
}

#[test]
fn junction_undefined_keeps_calls() {
    // J truncated inside the Trp118 codon: V and J are called, the junction isn't
    let j = ungapped(human(), "IGHJ4*02");
    let read = [ungapped(human(), "IGHV1-2*02"), "CC".into(), ungapped(human(), "IGHD3-3*01"), "GG".into(), j[..16].to_string()].concat();
    let r = engine().annotate("r", &read);
    assert_eq!((r.annotation_status.as_str(), r.failure_reason.as_deref()), ("unassigned", Some("junction_undefined")));
    assert_eq!((r.v_call.as_deref(), r.d_call.as_deref()), (Some("IGHV1-2*02"), Some("IGHD3-3*01")));
    assert!(r.j_call.as_deref().is_some_and(|c| c.starts_with("IGHJ")));
    // the regions the read covers: CDR3 whole, FWR4 its first two bases (ruling A4)
    assert_eq!((r.junction.clone(), r.junction_length, r.fwr4.as_deref()), (None, None, Some("TG")));
    assert_eq!((r.fwr3_start, r.fwr3_end, r.cdr3_start, r.cdr3_end), (Some(175), Some(288), Some(289), Some(345)));
}

#[test]
fn unassignable_reads() {
    let e = engine();
    let status = |seq: &str| {
        let r = e.annotate("x", seq);
        (r.annotation_status, r.failure_reason)
    };
    let no_v = ("unassigned".to_string(), Some("no_v_hit".to_string()));
    assert_eq!(status(&String::from_utf8(rand_seq(400, 7)).unwrap()), no_v);
    assert_eq!(status("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"), no_v);
    assert_eq!(status("ACGT"), no_v);
    assert_eq!(status(""), ("invalid_input".into(), Some("empty_sequence".into())));
    assert_eq!(status("ACGTX"), ("invalid_input".into(), Some("invalid_characters".into())));
    // N and IUPAC bases never crash
    let mut read = igh_read().into_bytes();
    for (i, b) in [(5, b'N'), (100, b'R'), (200, b'n'), (300, b'Y'), (350, b'N')] {
        read[i] = b;
    }
    let r = e.annotate("n", std::str::from_utf8(&read).unwrap());
    assert_eq!((r.annotation_status.as_str(), r.v_call.as_deref()), ("annotated", Some("IGHV1-2*02")));
    // a TCR read against an IG-only database
    let tcr = [ungapped(human(), "TRBV7-2*01"), "GGGACAGG".into(), ungapped(human(), "TRBJ2-1*01")].concat();
    assert_eq!(e.annotate("t", &tcr).locus.as_deref(), Some("TRB"));
    let mouse = load("c57bl6");
    let ig_only = Engine::new(&mouse, &EngineOptions::default()).unwrap();
    let r = ig_only.annotate("t", &tcr);
    assert_eq!((r.annotation_status, r.failure_reason), no_v);
}

#[test]
fn mixed_species_fields() {
    let db = load("human+c57bl6");
    let e = Engine::new(&db, &EngineOptions::default()).unwrap();
    let find = |species: &str, seg: Segment| {
        db.alleles.iter().find(|a| a.species == species && a.segment == seg && a.locus == Locus::IGH && a.anchor.is_some()).unwrap()
    };
    let hv = db.alleles.iter().find(|a| a.name == "IGHV1-2*02" && a.species == "human").unwrap();
    for (v, j) in [(hv, find("human", Segment::J)), (find("mouse", Segment::V), find("mouse", Segment::J))] {
        let read = [v.sequence(), "CCGTATTACGATTTTTGG".into(), j.sequence()].concat();
        let r = e.annotate("x", &read);
        assert_eq!(r.annotation_status, "annotated", "{}: {:?}", v.name, r.failure_reason);
        assert!(v.id.contains('|'), "mixed database ids carry a species prefix: {}", v.id);
        assert_eq!(r.v_species.as_deref(), Some(v.species.as_str()));
        assert_eq!(r.j_species.as_deref(), Some(j.species.as_str()));
        assert!(r.j_call_ties.as_ref().unwrap().contains(&j.name));
        let (ties, species) = (r.v_call_ties.unwrap(), r.v_species_ties.unwrap());
        assert!(ties.contains(&v.name), "{}: {ties:?}", v.id);
        assert_eq!(r.v_call.as_ref(), Some(&ties[0]));
        assert_eq!(ties.len(), species.len());
        assert!(ties.iter().all(|t| !t.contains('|')), "{ties:?}");
        assert!(species.iter().all(|s| *s == v.species), "{species:?}");
    }
}

#[test]
fn loci_option_restricts_the_search() {
    let tcr = [ungapped(human(), "TRBV7-2*01"), "GGGACAGG".into(), ungapped(human(), "TRBJ2-1*01")].concat();
    let ig = Engine::new(human(), &EngineOptions { loci: Some(vec![Locus::IGH, Locus::IGK, Locus::IGL]) }).unwrap();
    assert_eq!(ig.annotate("t", &tcr).failure_reason.as_deref(), Some("no_v_hit"));
    assert_eq!(ig.annotate("r", &igh_read()).annotation_status, "annotated");
}

#[test]
fn deterministic_across_threads_and_batches() {
    let g = generated(5, "memory", 300, |p| p.mutation = Spec::Uniform(0.0, 0.2));
    let mut recs: Vec<(String, String)> = (0..g.len()).map(|i| (g.str("sequence_id", i).unwrap().to_string(), g.sequence(i).to_string())).collect();
    // reverse complements, junk and invalid records in the same stream
    for i in 0..50 {
        recs.push((format!("rc{i}"), String::from_utf8(revcomp(g.sequence(i).as_bytes())).unwrap()));
    }
    recs.push(("junk".into(), String::from_utf8(rand_seq(300, 3)).unwrap()));
    recs.push(("empty".into(), String::new()));
    let inputs = [Input::Records(recs)];
    let tsv = |threads, batch_size| {
        let mut out = Vec::new();
        let o = RunOptions { threads: Some(threads), batch_size, ..Default::default() };
        let counts = run_to_writer(&inputs, human(), &o, &mut out).unwrap();
        (out, counts)
    };
    let (a, ca) = tsv(1, 7);
    let (b, cb) = tsv(4, 4096);
    let (c, _) = tsv(3, 64);
    assert_eq!(ca, cb);
    assert!(a == b && b == c, "output depends on threads or batch size");
    assert!(ca.annotated >= 330, "{ca:?}");
    assert_eq!((ca.invalid_input, ca.error), (1, 0));
}

/// IGHV1-2*02 with the `n` V bases before Cys104 all mutated and the Cys codon replaced by
/// `cys_codon`, then CC + IGHD3-3*01 + GG + IGHJ4*02; and the Cys104 codon's read position.
fn mutated_v_end_read(n: usize, cys_codon: &[u8]) -> (String, usize) {
    let db = human();
    let mut v = ungapped(db, "IGHV1-2*02").into_bytes();
    let cys = at_imgt(db, "IGHV1-2*02", 310);
    assert_eq!(&v[cys..cys + 3], b"TGT");
    for b in &mut v[cys - n..cys] {
        *b = match *b { b'A' => b'C', b'C' => b'G', b'G' => b'T', _ => b'A' };
    }
    v[cys..cys + 3].copy_from_slice(cys_codon);
    ([String::from_utf8(v).unwrap(), "CC".into(), ungapped(db, "IGHD3-3*01"), "GG".into(), ungapped(db, "IGHJ4*02")].concat(), cys)
}

#[test]
fn cys_projected_past_a_mutated_v_end() {
    // Ruling B2: the 24 V bases before Cys104 all mutated, so V's alignment stops short of it; the
    // intact TGT is projected along V's diagonal and the junction starts there.
    let (read, cys) = mutated_v_end_read(24, b"TGT");
    let r = engine().annotate("r1", &read);
    assert!(r.v_sequence_end.unwrap() < cys as i64, "V reaches Cys: {:?}", r.v_sequence_end);
    assert_eq!(r.junction.as_deref().map(|j| &j[..3]), Some("TGT"));
    assert_eq!(r.cdr3_start, Some(cys as i64 + 4));
    assert_eq!(r.annotation_status, "annotated");
}

#[test]
fn cys_not_projected_unless_cys_and_near() {
    // B2 projects only a TGT/TGC codon, at most 50 nt past V's last base.
    let (read, cys) = mutated_v_end_read(24, b"TCT");
    let r = engine().annotate("r1", &read);
    assert!(r.v_sequence_end.unwrap() < cys as i64);
    assert_eq!((r.junction, r.failure_reason.as_deref()), (None, Some("junction_undefined")));
    let (read, cys) = mutated_v_end_read(60, b"TGT");
    let r = engine().annotate("r2", &read);
    assert!(r.v_sequence_end.unwrap() + 50 < cys as i64, "V end {:?}, Cys {cys}", r.v_sequence_end);
    assert_eq!(r.junction, None);
}

#[test]
fn no_v_no_short_j() {
    // B1 needs a V call: random bases then the first 14 nt of IGHJ4*02 get no call at all.
    let j = ungapped(human(), "IGHJ4*02");
    let read = [String::from_utf8(rand_seq(300, 5)).unwrap(), j[..14].to_string()].concat();
    let r = engine().annotate("r1", &read);
    assert_eq!((r.v_call, r.j_call, r.failure_reason.as_deref()), (None, None, Some("no_v_hit")));
}

#[test]
fn d_margin() {
    // Ruling B3: an exact 5-nt D gains only 5 over no D and is not called; 10 nt are.
    let db = human();
    let d = ungapped(db, "IGHD3-3*01");
    let read = |run: &str| [ungapped(db, "IGHV1-2*02"), "CCCCC".into(), run.to_string(), "CCCCC".into(), ungapped(db, "IGHJ4*02")].concat();
    let r = engine().annotate("r5", &read(&d[6..11]));
    assert_eq!(r.d_call, None, "5-nt D called");
    let r = engine().annotate("r10", &read(&d[6..16]));
    assert!(r.d_call_ties.as_ref().is_some_and(|t| t.iter().any(|x| x == "IGHD3-3*01")), "{:?}", r.d_call_ties);
}

fn transition(b: u8) -> u8 {
    match b {
        b'A' => b'G',
        b'G' => b'A',
        b'C' => b'T',
        _ => b'C',
    }
}

#[test]
fn weak_j_junction_ends_at_the_j_anchor() {
    // Record 0029, Wave 2: a full-length IGHJ1*01 with 20 of 52 bases mutated (Trp118 codon kept)
    // is called on the weak-J path, and the junction runs from the V's Cys104 to the end of the J's
    // Trp118 codon at its true read position.
    let db = human();
    let (v, j) = (ungapped(db, "IGHV3-23*01"), ungapped(db, "IGHJ1*01").into_bytes());
    let anchor = |name: &str| db.alleles.iter().find(|a| a.name == name).unwrap().anchor.unwrap();
    let (va, w) = (anchor("IGHV3-23*01"), anchor("IGHJ1*01"));
    let m: Vec<u8> = j.iter().enumerate().map(|(i, &b)| if i % 5 % 2 == 1 && !(w..w + 3).contains(&i) { transition(b) } else { b }).collect();
    let np = "GCGAAAGATCGGGGATATAGTAGTGGCTGGTAC";
    let read = [v.as_str(), np, std::str::from_utf8(&m).unwrap()].concat();
    let r = engine().annotate("weak", &read);
    assert_eq!(r.j_call.as_deref(), Some("IGHJ1*01"));
    let end = v.len() + np.len() + w + 3;
    assert_eq!(&read[end - 3..end], "TGG");
    assert_eq!(r.junction.as_deref(), Some(&read[va..end]));
}

#[test]
fn mutated_j_start_cannot_pass_a_held_anchor() {
    // Wave 2, fix round 1: benchmark igh_bnab g0327 (IGHJ4*02 mutated at its 5' end). Its J
    // alignment holds Trp118, but the joint step used to start J at germline 31, past the anchor,
    // leaving no junction. J now starts no later than the anchor and the junction is truth's.
    let read = "GCGGTACAACTCGTGGAGTCCAGGGCAGGAGTCCTATGTCCTGGGGGGTCCCCGCGACACTCCTATTCCGTCTCCAGAGACGCAGCAGATGATTGTGACCTGAGGTCGGTCCGCCAACTTGGAGGGACAGGGTCACTGTGGGTCTCTCGTCGCTCTTCGAACCGTGGGAGTGTAGGCTATGGCGATTCCGAAAATAATCACTTCAATATCCTCAGAGACAATTCCAAGACCTCTCTGGGTCTGCGGATGGTCAGTCCGAGAGTCGAAACAACCGCCGCATTTAACTGTGAGAGAGAGAAGTATTGGGTGACCACTTATTGGCGACGCATACTACCAAAGGGTCTCGATTGTTGGAGCCAACAGAACCTGGTCACCGCCCCCGTCG";
    let r = engine().annotate("g0327", read);
    assert_eq!(r.annotation_status, "annotated", "{:?}", r.failure_reason);
    assert_eq!(r.junction.as_deref(), Some("TGTGAGAGAGAGAAGTATTGGGTGACCACTTATTGGCGACGCATACTACCAAAGGGTCTCGATTGTTGG"));
}
