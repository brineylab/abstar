//! `AaEngine` end to end (spec §Fields, §Input validation) and amino-acid mode through the run
//! pipeline.
mod common;
use abstar_annotate::aa::AaEngine;
use abstar_annotate::engine::{Annotate, EngineOptions};
use abstar_annotate::input::Input;
use abstar_annotate::pipeline::{run_in_memory, RunOptions, SequenceType};
use abstar_annotate::schema::Record;
use abstar_germline::Database;
use arrow_array::cast::AsArray;
use arrow_array::Array;
use std::sync::OnceLock;

const IGH_AA: &str = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLGWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARDRITIFGVVIIPDYFDYWGQGTLVTVSS";
const REGIONS: [&str; 7] = [
    "QVQLVQSGAEVKKPGASVKVSCKAS",
    "GYTFTGYY",
    "MHWVRQAPGQGLGWMGW",
    "INPNSGGT",
    "NYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC",
    "ARDRITIFGVVIIPDYFDY",
    "WGQGTLVTVSS",
];
const JUNCTION: &str = "CARDRITIFGVVIIPDYFDYW";

fn human() -> &'static Database {
    static D: OnceLock<Database> = OnceLock::new();
    D.get_or_init(common::human_db)
}

fn engine() -> &'static AaEngine<'static> {
    static E: OnceLock<AaEngine<'static>> = OnceLock::new();
    E.get_or_init(|| AaEngine::new(human(), &EngineOptions::default()).unwrap())
}

fn ann(seq: &str) -> Record {
    engine().annotate("r", seq)
}

fn status(r: &Record) -> (&str, Option<&str>) {
    (r.annotation_status.as_str(), r.failure_reason.as_deref())
}

/// A region's `_aa`, `_length`, `_start` and `_end`.
type Region = (Option<String>, Option<i64>, Option<i64>, Option<i64>);

/// Every region's fields.
fn regions(r: &Record) -> Vec<Region> {
    vec![
        (r.fwr1_aa.clone(), r.fwr1_length, r.fwr1_start, r.fwr1_end),
        (r.cdr1_aa.clone(), r.cdr1_length, r.cdr1_start, r.cdr1_end),
        (r.fwr2_aa.clone(), r.fwr2_length, r.fwr2_start, r.fwr2_end),
        (r.cdr2_aa.clone(), r.cdr2_length, r.cdr2_start, r.cdr2_end),
        (r.fwr3_aa.clone(), r.fwr3_length, r.fwr3_start, r.fwr3_end),
        (r.cdr3_aa.clone(), r.cdr3_length, r.cdr3_start, r.cdr3_end),
        (r.fwr4_aa.clone(), r.fwr4_length, r.fwr4_start, r.fwr4_end),
    ]
}

fn igh_regions() -> Vec<Region> {
    REGIONS.iter().map(|s| (Some(s.to_string()), Some(s.len() as i64), None, None)).collect()
}

#[test]
fn aa_end_to_end_igh() {
    let r = ann(IGH_AA);
    assert_eq!(status(&r), ("annotated", None));
    assert_eq!(r.sequence_id, "r");
    assert_eq!(r.sequence, None);
    assert_eq!(r.sequence_aa.as_deref(), Some(IGH_AA));
    assert_eq!(r.locus.as_deref(), Some("IGH"));
    assert_eq!(r.v_gene.as_deref(), Some("IGHV1-2"));
    assert_eq!(r.j_gene.as_deref(), Some("IGHJ4"));
    for (call, ties, species) in [(&r.v_call, &r.v_call_ties, &r.v_species_ties), (&r.j_call, &r.j_call_ties, &r.j_species_ties)] {
        let (ties, species) = (ties.as_ref().unwrap(), species.as_ref().unwrap());
        assert_eq!(ties.len(), species.len());
        assert_eq!(Some(&ties[0]), call.as_ref());
    }
    assert_eq!(r.v_species.as_deref(), Some("human"));
    assert!(r.v_score.unwrap() > 0.0 && r.j_score.unwrap() > 0.0);
    // V covers residues 0..98 with one substitution; J 111..126 exactly
    assert_eq!(r.v_identity, Some(97.0 / 98.0));
    assert_eq!(r.j_identity, Some(1.0));
    assert_eq!(regions(&r), igh_regions());
    assert_eq!(r.junction_aa.as_deref(), Some(JUNCTION));
    assert_eq!(r.v_mutations_aa, Some(vec!["IMGT51:E>G".to_string()]));
    assert_eq!(r.v_mutation_count_aa, Some(1));
    assert_eq!(r.v_insertions_aa, Some(vec![]));
    assert_eq!(r.v_deletions_aa, Some(vec![]));
    assert_eq!((r.j_mutations_aa, r.j_mutation_count_aa, r.j_insertions_aa), (Some(vec![]), Some(0), Some(vec![])));
    assert_eq!((r.productive, r.stop_codon, r.vj_in_frame), (Some(true), Some(false), None));
    assert_eq!(r.sequence_vdj_aa.as_deref(), Some(IGH_AA));
    assert_eq!(r.gene_segment_mask_aa, Some(format!("{}{}{}", "V".repeat(98), "N".repeat(13), "J".repeat(15))));
    assert_eq!(r.cdr_mask_aa.as_ref().map(String::len), Some(IGH_AA.len()));
    assert_eq!(r.nongermline_mask_aa.as_ref().map(|m| m.matches('1').count()), Some(14));
    let (s, g) = (r.sequence_alignment_aa.unwrap(), r.germline_alignment_aa.unwrap());
    assert_eq!((s.replace('.', ""), s.len()), (IGH_AA.to_string(), g.len()));
}

#[test]
fn lowercase_input_is_uppercased() {
    let r = ann(&IGH_AA.to_lowercase());
    assert_eq!(status(&r), ("annotated", None));
    assert_eq!(r.sequence_aa.as_deref(), Some(IGH_AA));
    assert_eq!(r.junction_aa.as_deref(), Some(JUNCTION));
}

macro_rules! assert_null {
    ($r:expr, $($f:ident),+ $(,)?) => { $( assert!($r.$f.is_none(), concat!(stringify!($f), " is not null")); )+ };
}

#[test]
fn aa_null_contract() {
    let r = ann(IGH_AA);
    assert_null!(r, sequence, rev_comp, vj_in_frame, v_cigar, d_cigar, j_cigar);
    assert_null!(
        r,
        v_sequence_start, v_sequence_end, v_germline_start, v_germline_end,
        d_sequence_start, d_sequence_end, d_germline_start, d_germline_end,
        j_sequence_start, j_sequence_end, j_germline_start, j_germline_end,
        c_sequence_start, c_sequence_end, c_germline_start, c_germline_end,
    );
    assert_null!(r, junction, junction_length, np1, np1_length, np2, np2_length);
    assert_null!(r, sequence_vdj, sequence_vdjc, sequence_vdjc_aa, sequence_alignment, germline_alignment);
    assert_null!(r, fwr1, fwr1_start, fwr1_end, cdr1, cdr1_start, cdr1_end, fwr2, fwr2_start, fwr2_end);
    assert_null!(r, cdr2, cdr2_start, cdr2_end, fwr3, fwr3_start, fwr3_end, cdr3, cdr3_start, cdr3_end);
    assert_null!(r, fwr4, fwr4_start, fwr4_end);
    assert_null!(r, cdr_mask, gene_segment_mask, nongermline_mask);
    assert_null!(r, v_mutations, v_insertions, v_deletions, v_mutation_count, v_mutation_frequency);
    assert_null!(r, j_mutations, j_insertions, j_deletions, j_mutation_count, j_mutation_frequency);
    assert_null!(r, d_call, d_gene, d_species, d_call_ties, d_species_ties, d_inverted, d_score, d_identity);
    assert_null!(r, d_mutations, d_mutations_aa, d_insertions, d_deletions, d_mutation_count, d_mutation_count_aa);
    assert_null!(r, d_mutation_frequency, d_mutation_frequency_aa);
    assert_null!(r, c_call, c_gene, c_species, c_call_ties, c_species_ties, c_score, c_identity);
    assert_null!(r, c_mutations, c_mutations_aa, c_insertions, c_deletions, c_mutation_count, c_mutation_count_aa);
    assert_null!(r, c_mutation_frequency, c_mutation_frequency_aa);
}

#[test]
fn aa_statuses() {
    let r = ann("");
    assert_eq!(status(&r), ("invalid_input", Some("empty_sequence")));
    let r = ann("QVQ-LV");
    assert_eq!(status(&r), ("invalid_input", Some("invalid_characters")));
    assert_eq!((r.sequence_aa.as_deref(), r.sequence.as_deref()), (Some("QVQ-LV"), None));
    let ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
    let r = ann(ubiquitin);
    assert_eq!(status(&r), ("unassigned", Some("no_v_hit")));
    assert_eq!((r.sequence_aa.as_deref(), r.sequence.as_deref()), (Some(ubiquitin), None));
    let r = ann(&IGH_AA[..96]);
    assert_eq!(status(&r), ("unassigned", Some("no_j_hit")));
    assert_eq!(r.v_gene.as_deref(), Some("IGHV1-2"));
    let mut q = IGH_AA.to_string();
    q.drain(90..98);
    let r = ann(&q);
    assert_eq!(status(&r), ("unassigned", Some("junction_undefined")));
    assert_eq!((r.junction_aa.as_deref(), r.cdr3_aa.as_deref(), r.productive), (None, None, Some(false)));
}

#[test]
fn x_in_junction_is_nonproductive() {
    let mut q = IGH_AA.as_bytes().to_vec();
    q[101] = b'X'; // CDR3 `ARDRI[T]IF...`
    let r = ann(std::str::from_utf8(&q).unwrap());
    assert_eq!(status(&r), ("annotated", None));
    assert_eq!(r.junction_aa.as_deref(), Some("CARDRIXIFGVVIIPDYFDYW"));
    assert_eq!(r.productive, Some(false));
}

#[test]
fn trailing_stop_is_productive() {
    let r = ann(&format!("{IGH_AA}*"));
    assert_eq!(status(&r), ("annotated", None));
    assert_eq!(r.sequence_vdj_aa.as_deref(), Some(IGH_AA));
    assert_eq!((r.productive, r.stop_codon), (Some(true), Some(false)));
}

#[test]
fn leader_is_excluded() {
    let r = ann(&format!("MDWTWRILFLVAAATGAHS{IGH_AA}"));
    assert_eq!(status(&r), ("annotated", None));
    assert_eq!(regions(&r), igh_regions());
    assert_eq!(r.sequence_vdj_aa.as_deref(), Some(IGH_AA));
}

#[test]
fn fab_with_ch1() {
    let ch1 = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSC";
    let r = ann(&format!("{IGH_AA}{ch1}"));
    assert_eq!(status(&r), ("annotated", None));
    assert_eq!(regions(&r), igh_regions());
    assert_eq!(r.junction_aa.as_deref(), Some(JUNCTION));
    assert_eq!(r.productive, Some(true));
}

#[test]
fn scfv_one_row() {
    let vk = "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK";
    let r = ann(&format!("{IGH_AA}GGGGSGGGGSGGGGS{vk}"));
    assert_eq!(status(&r), ("annotated", None));
}

fn ungapped(name: &str) -> String {
    human().alleles.iter().find(|a| a.name == name).unwrap_or_else(|| panic!("no allele {name}")).sequence()
}

#[test]
fn nucleotide_in_aa_mode() {
    // `igh_read()` of tests/engine.rs: IGHV1-2*02 (A>G at IMGT 152) + CC + IGHD3-3*01 + GG + IGHJ4*02
    let a = human().alleles.iter().find(|a| a.name == "IGHV1-2*02").unwrap();
    let g = a.sequence_gapped.as_bytes();
    let mut v = ungapped("IGHV1-2*02").into_bytes();
    v[g[..151].iter().filter(|&&b| b != b'.').count()] = b'G';
    let read = [String::from_utf8(v).unwrap(), "CC".into(), ungapped("IGHD3-3*01"), "GG".into(), ungapped("IGHJ4*02")].concat();
    assert_eq!(read.len(), 379);
    assert_eq!(status(&ann(&read)), ("unassigned", Some("no_v_hit")));
}

#[test]
fn pipeline_aa_mode() {
    let inputs = [Input::Records(vec![("a".into(), IGH_AA.into()), ("b".into(), "QV\tQ".into())])];
    let o = RunOptions { sequence_type: SequenceType::Aa, ..Default::default() };
    let b = run_in_memory(&inputs, human(), &o).unwrap();
    assert_eq!(b.len(), 1);
    let b = &b[0];
    let col = |n: &str| b.column_by_name(n).unwrap().as_string::<i32>().clone();
    let (seq, aa, st, fr) = (col("sequence"), col("sequence_aa"), col("annotation_status"), col("failure_reason"));
    assert_eq!((st.value(0), seq.is_null(0), aa.value(0)), ("annotated", true, IGH_AA));
    assert_eq!((st.value(1), fr.value(1), seq.is_null(1)), ("invalid_input", "invalid_characters", true));
    assert_eq!(aa.value(1), "QV Q");
}

#[test]
fn sequence_type_parses() {
    assert_eq!("aa".parse::<SequenceType>(), Ok(SequenceType::Aa));
    assert_eq!("nt".parse::<SequenceType>(), Ok(SequenceType::Nt));
    assert_eq!(SequenceType::default(), SequenceType::Nt);
    assert_eq!(SequenceType::Aa.to_string(), "aa");
    assert_eq!("protein".parse::<SequenceType>(), Err("sequence_type must be 'nt' or 'aa', got \"protein\"".to_string()));
}
