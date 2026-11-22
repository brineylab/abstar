use abstar_germline::gap::{gap_v, GapResult};
use abstar_germline::{load, Allele};
use std::path::PathBuf;

fn alleles(db: &str) -> Vec<Allele> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../python/abstar/germline_dbs")
        .join(db);
    load(&dir).unwrap().alleles
}

fn human() -> Vec<Allele> {
    alleles("human")
}

fn get<'a>(db: &'a [Allele], id: &str) -> &'a Allele {
    db.iter().find(|a| a.id == id).unwrap()
}

fn without(db: &[Allele], id: &str) -> Vec<Allele> {
    db.iter().filter(|a| a.id != id).cloned().collect()
}

fn strip(s: &str) -> String {
    s.replace('.', "")
}

fn gap_one(name: &str, seq: &str, templates: &[Allele]) -> GapResult {
    gap_v(&[(name.to_string(), seq.to_string())], templates)
        .pop()
        .unwrap()
}

/// Regap `gapped` (stripped) as `id`, with human templates minus `id`.
fn regap(id: &str, gapped: &str) -> GapResult {
    gap_one(id, &strip(gapped), &without(&human(), id))
}

#[test]
fn regaps_a_human_allele_from_human_templates() {
    let db = human();
    let shipped = &get(&db, "IGHV1-2*02").sequence_gapped;
    let r = regap("IGHV1-2*02", shipped);
    assert_eq!(r.failure, None);
    assert_eq!(r.gapped.as_deref(), Some(shipped.as_str()));
    assert_eq!(r.cdr1_length, Some(8));
}

#[test]
fn macaque_style_extra_gaps_are_ignored() {
    let db = human();
    let shipped = &get(&db, "IGHV1-2*02").sequence_gapped;
    let kimdb = format!(
        "{}...{}...{}",
        &shipped[..27],
        &shipped[27..120],
        &shipped[120..]
    );
    let r = gap_one("IGHV1-2*02", &kimdb, &without(&db, "IGHV1-2*02"));
    assert_eq!(r.gapped.as_deref(), Some(shipped.as_str()));
}

#[test]
fn five_prime_partial() {
    let db = human();
    let shipped = &get(&db, "IGHV3-23*01").sequence_gapped;
    assert!(!shipped[..16].contains('.'));
    let r = regap("IGHV3-23*01", &shipped[16..]);
    assert_eq!(
        r.gapped,
        Some(format!("{}{}", ".".repeat(16), &shipped[16..]))
    );
}

#[test]
fn trailing_partial_codon_kept() {
    let db = human();
    let shipped = &get(&db, "IGKV1-39*01").sequence_gapped;
    let r = regap("IGKV1-39*01", &format!("{shipped}GA"));
    assert_eq!(r.failure, None);
    assert!(r.gapped.unwrap().ends_with(&format!("{shipped}GA")));
}

#[test]
fn iupac_bases_are_gapped() {
    let db = human();
    let shipped = &get(&db, "IGHV1-2*02").sequence_gapped;
    let idx = 3 * (80 - 1) + 1; // middle base of FR3 position 80
    let mut mutated = shipped.clone().into_bytes();
    mutated[idx] = b'N';
    let mutated = String::from_utf8(mutated).unwrap();
    let r = regap("IGHV1-2*02", &mutated.to_lowercase());
    assert_eq!(r.failure, None);
    let gapped = r.gapped.unwrap();
    assert_eq!(gapped.as_bytes()[idx], b'N');
    assert_eq!(gapped, mutated);
}

#[test]
fn unknown_locus_name_uses_all_templates() {
    let db = human();
    let shipped = &get(&db, "IGHV1-2*02").sequence_gapped;
    let r = gap_one("seq1", &strip(shipped), &without(&db, "IGHV1-2*02"));
    assert_eq!(r.gapped.as_deref(), Some(shipped.as_str()));
}

#[test]
fn missing_locus_falls_back() {
    let db = human();
    let query = strip(&get(&db, "TRBV7-2*01").sequence_gapped);
    let mouse = alleles("c57bl6");
    let r = gap_one("TRBV7-2*01", &query, &mouse);
    match r.failure.as_deref() {
        None => {
            assert!(r.warnings.iter().any(|w| w == "low_identity"), "{r:?}");
            let t = r.template.as_deref().unwrap();
            assert!(mouse.iter().any(|a| a.id == t));
            assert_eq!(r.template_species.as_deref(), Some("mouse"));
        }
        Some(reason) => assert!(
            reason == "framework_insertion" || reason == "cdr_too_long",
            "{r:?}"
        ),
    }
}

#[test]
fn stop_codon_is_a_warning() {
    let db = human();
    let shipped = &get(&db, "IGHV1-2*02").sequence_gapped;
    let mutated = format!("{}TAG{}", &shipped[..147], &shipped[150..]);
    let r = regap("IGHV1-2*02", &mutated);
    assert_eq!(r.gapped.as_deref(), Some(mutated.as_str()));
    assert!(r.warnings.iter().any(|w| w == "stop_codon"));
}

#[test]
fn prefers_cys104() {
    let db = human();
    let shipped = &get(&db, "IGHV1-2*02").sequence_gapped;
    let mutated = format!("{}TCT{}", &shipped[..309], &shipped[312..]);
    let r = regap("IGHV1-2*02", &mutated);
    assert_eq!(r.failure, None);
    assert!(r.warnings.iter().any(|w| w == "no_cys104"));
    assert_eq!(&r.gapped.unwrap()[309..312], "TCT");
}

#[test]
fn framework_insertion_fails() {
    let db = human();
    let shipped = &get(&db, "IGHV1-2*02").sequence_gapped;
    let mutated = format!("{}GGCAGCGGC{}", &shipped[..135], &shipped[135..]);
    let r = regap("IGHV1-2*02", &mutated);
    assert_eq!(r.failure.as_deref(), Some("framework_insertion"));
    assert_eq!(r.gapped, None);
}

#[test]
fn cdr_too_long_fails() {
    let db = human();
    let shipped = &get(&db, "IGHV1-2*02").sequence_gapped;
    // CDR1 is 27-30 and 35-38; insert 5 codons in its dotted middle (31-34).
    let mutated = format!("{}GGTAGTGGTAGTGGT{}", &shipped[..90], &shipped[90..]);
    let r = regap("IGHV1-2*02", &mutated);
    assert_eq!(r.failure.as_deref(), Some("cdr_too_long"));
}

#[test]
fn non_v_fails() {
    let db = human();
    let j = strip(&get(&db, "IGHJ4*02").sequence_gapped);
    let mut state: u64 = 0x2545_F491_4F6C_DD1D;
    let random: String = (0..300)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state % 4) as usize] as char
        })
        .collect();
    let queries = [
        ("IGHJ4*02".to_string(), j),
        ("IGHVrandom".to_string(), random),
        ("IGHVempty".to_string(), String::new()),
    ];
    let r = gap_v(&queries, &db);
    assert_eq!(r[0].failure.as_deref(), Some("no_alignment"));
    assert_eq!(r[1].failure.as_deref(), Some("no_alignment"));
    assert_eq!(r[2].failure.as_deref(), Some("empty_sequence"));
}

#[test]
fn order_preserved() {
    let db = human();
    let queries: Vec<(String, String)> = db
        .iter()
        .filter(|a| a.id.starts_with("IGKV"))
        .take(5)
        .map(|a| (a.id.clone(), a.sequence()))
        .collect();
    let r = gap_v(&queries, &db);
    let names: Vec<_> = r.iter().map(|g| g.name.clone()).collect();
    let expected: Vec<_> = queries.iter().map(|q| q.0.clone()).collect();
    assert_eq!(names, expected);
}
