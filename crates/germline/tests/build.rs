use abstar_germline::*;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

/// `v_gene_delineations[0].aligned_sequence` of IGHV4-28*05 in the OGRDB fixture.
const IGHV4_28_05_GAPPED: &str = "CAGGTGCAGCTGCAGGAGTCGGGCCCA...GGACTGGTGAAGCCTTCGGACACCCTGTCCCTCACCTGCGCTGTCTCTGGTTACTCCATCAGC.........AGTAGTAACTGGTGGGGCTGGATCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGTACATCTATTATAGT.........GGGAGCATCTACTACAACCCGTCCCTCAAG...AGTCGAGTCACCATGTCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGTGGACACGGCCGTGTATTACTGTGCGAGAAA";

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name)
}

fn opts() -> BuildOptions {
    BuildOptions { overwrite: false, abstar_version: "1.0.0-test".into() }
}

fn put(dir: &Path, name: &str, text: &str) -> PathBuf {
    let p = dir.join(name);
    fs::write(&p, text).unwrap();
    p
}

fn allele<'a>(db: &'a Database, id: &str) -> &'a Allele {
    db.alleles.iter().find(|a| a.id == id).unwrap_or_else(|| panic!("no allele {id}"))
}

fn no_dbs() -> BTreeMap<String, PathBuf> {
    BTreeMap::new()
}

#[test]
fn v_anchor_ighv4_28_05() {
    assert_eq!(find_v_anchor(IGHV4_28_05_GAPPED), Some(285));
}

#[test]
fn j_anchor_ighj6_03() {
    assert_eq!(
        find_j_anchor("ATTACTACTACTACTACTACATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCCTCAG"),
        Some(29)
    );
}

#[test]
fn j_anchor_ighj4_02() {
    assert_eq!(find_j_anchor("ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"), Some(14));
}

#[test]
fn j_anchor_none_without_motif() {
    assert_eq!(find_j_anchor("ACGTACGTACGTACGTACGT"), None);
}

#[test]
fn anchors_do_not_match_ambiguous_bases() {
    // [ACGT] stays strict: an IUPAC code never satisfies the motif or the Cys codon
    let j = "ACTACTTTGACTACTGGGGCCAGGGAACC";
    assert_eq!(find_j_anchor(j), Some(14));
    assert_eq!(find_j_anchor(&j.replace("TGGGGCCAGGG", "TGGGGCSAGGG")), None);
    assert_eq!(find_j_anchor(&j.replace("TGGGGCCAGGG", "TGGGGCCAGGS")), None);
    let mut v = IGHV4_28_05_GAPPED.to_string();
    assert_eq!(&v[309..312], "TGT");
    v.replace_range(311..312, "Y"); // TGY could be either Cys codon, but is not called
    assert_eq!(find_v_anchor(&v), None);
}

#[test]
fn ogrdb_build() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("db");
    let r = build(&data("ogrdb.toml"), &out, &no_dbs(), &opts()).unwrap();
    assert_eq!(r.name, "ogrdb-test");
    assert_eq!(r.path, out);
    assert_eq!(r.n_alleles, 4);
    assert_eq!(r.warnings, vec![]);
    let counts: BTreeMap<String, usize> =
        [("IGH/V/human", 1), ("IGH/J/human", 2), ("IGH/C/human", 1)]
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect();
    assert_eq!(r.counts, counts);
    assert!(r.restricted.is_empty() && r.passed_through.is_empty());

    let db = load(&out).unwrap();
    assert_eq!(db.alleles.len(), 4);
    assert!(db.manifest.warnings.is_empty());
    let j6 = allele(&db, "IGHJ6*03");
    assert_eq!(j6.anchor, Some(29));
    assert_eq!(j6.functionality, Functionality::F);
    assert_eq!(allele(&db, "IGHJ4*02").anchor, Some(14));
    let v = allele(&db, "IGHV4-28*05");
    assert_eq!(v.gene, "IGHV4-28");
    assert_eq!(v.name, "IGHV4-28*05");
    assert_eq!(v.sequence_gapped, IGHV4_28_05_GAPPED);
    assert_eq!(v.anchor, Some(285));
    assert_eq!(v.source, "Homo_sapiens_IGH_VDJ_rev_9_ex_trimmed.json");
    let c = allele(&db, "IGHG1*01");
    assert_eq!((c.segment, c.anchor, c.functionality), (Segment::C, None, Functionality::Unknown));
    assert_eq!(c.source, "c.fasta");

    let m = &db.manifest;
    assert_eq!((m.species.clone(), m.loci.clone()), (vec!["human".to_string()], vec![Locus::IGH]));
    assert_eq!(m.abstar_version, "1.0.0-test");
    assert_eq!(m.built.len(), 10);
    assert_eq!(m.provenance.len(), 1);
    let p = &m.provenance[0];
    assert_eq!(p.files, ["Homo_sapiens_IGH_VDJ_rev_9_ex_trimmed.json", "c.fasta"]);
    assert_eq!(p.sha256.len(), 2);
    assert!(p.description.contains("IGH_VDJ_extended"), "{}", p.description);
    assert!(p.description.contains("2024-10-12"), "{}", p.description);
}

#[test]
fn ogrdb_anchor_disagreement_warns_and_keeps_ogrdb() {
    let tmp = tempfile::tempdir().unwrap();
    let json = fs::read_to_string(data("Homo_sapiens_IGH_VDJ_rev_9_ex_trimmed.json")).unwrap();
    assert_eq!(json.matches("\"j_cdr3_end\": 30,").count(), 1);
    put(tmp.path(), "set.json", &json.replace("\"j_cdr3_end\": 30,", "\"j_cdr3_end\": 31,"));
    let recipe = put(
        tmp.path(),
        "r.toml",
        "name = \"x\"\n[[source]]\nspecies = \"human\"\nogrdb_json = [\"set.json\"]\n",
    );
    let out = tmp.path().join("db");
    let r = build(&recipe, &out, &no_dbs(), &opts()).unwrap();
    assert_eq!(
        r.warnings,
        vec![BuildWarning {
            allele: "IGHJ6*03".into(),
            message: "motif anchor 29 disagrees with OGRDB j_cdr3_end 31".into(),
        }]
    );
    let db = load(&out).unwrap();
    assert_eq!(allele(&db, "IGHJ6*03").anchor, Some(30));
    assert_eq!(db.manifest.warnings, r.warnings);
}

#[test]
fn missing_anchors_warn_once_per_allele() {
    let tmp = tempfile::tempdir().unwrap();
    put(tmp.path(), "v.fasta", ">IGHV1-2*01\nCAGGTGCAGCTGGTGCAG\n");
    put(tmp.path(), "j.fasta", ">IGHJ1*01\nACGTACGTACGTACGTACGT\n");
    let recipe = put(
        tmp.path(),
        "r.toml",
        "name = \"x\"\n[[source]]\nspecies = \"human\"\nv_fasta = [\"v.fasta\"]\nj_fasta = [\"j.fasta\"]\n",
    );
    let r = build(&recipe, &tmp.path().join("db"), &no_dbs(), &opts()).unwrap();
    assert_eq!(r.warnings.len(), 2, "{:?}", r.warnings);
    let v = r.warnings.iter().find(|w| w.allele == "IGHV1-2*01").unwrap();
    assert!(v.message.starts_with("no Cys104 codon at IMGT positions 310-312"), "{}", v.message);
    assert!(v.message.contains("not IMGT-gapped"), "{}", v.message);
    let j = r.warnings.iter().find(|w| w.allele == "IGHJ1*01").unwrap();
    assert!(j.message.contains("motif"), "{}", j.message);
}

fn mixed_setup(dir: &Path, prefixes: &str) -> PathBuf {
    put(dir, "human_v.fasta", ">IGHV1-2*01\nCAGGTGCAGCTGGTGCAGTCTGGG\n");
    put(dir, "mouse_v.fasta", ">IGHV1-2*01\nGAGGTCCAGCTGCAACAGTCTGGA\n");
    put(
        dir,
        "r.toml",
        &format!(
            "name = \"mix\"\n{prefixes}\n\
             [[source]]\nspecies = \"human\"\nv_fasta = [\"human_v.fasta\"]\n\
             [[source]]\nspecies = \"mouse\"\nv_fasta = [\"mouse_v.fasta\"]\n"
        ),
    )
}

#[test]
fn mixed_species_prefixes() {
    let tmp = tempfile::tempdir().unwrap();
    let recipe = mixed_setup(tmp.path(), "[prefixes]\nhuman = \"hs\"\nmouse = \"mm\"\n");
    let out = tmp.path().join("db");
    build(&recipe, &out, &no_dbs(), &opts()).unwrap();
    let db = load(&out).unwrap();
    let hs = allele(&db, "hs|IGHV1-2*01");
    let mm = allele(&db, "mm|IGHV1-2*01");
    assert_eq!((hs.name.as_str(), hs.species.as_str()), ("IGHV1-2*01", "human"));
    assert_eq!((mm.name.as_str(), mm.species.as_str()), ("IGHV1-2*01", "mouse"));
    assert_eq!(db.manifest.species, ["human", "mouse"]);
    assert_eq!(db.manifest.prefixes.get("mouse").map(String::as_str), Some("mm"));
    assert_eq!(db.manifest.provenance.len(), 2);

    // a database source is re-prefixed from plain names, never stacked
    let again = put(
        tmp.path(),
        "again.toml",
        "name = \"again\"\n[prefixes]\nhuman = \"h\"\nmouse = \"m\"\n[[source]]\ndatabase = \"db\"\n",
    );
    let out2 = tmp.path().join("again");
    build(&again, &out2, &no_dbs(), &opts()).unwrap();
    let ids: Vec<String> = load(&out2).unwrap().alleles.into_iter().map(|a| a.id).collect();
    assert_eq!(ids, ["h|IGHV1-2*01", "m|IGHV1-2*01"]);
}

#[test]
fn mixed_without_prefix_fails() {
    let tmp = tempfile::tempdir().unwrap();
    let recipe = mixed_setup(tmp.path(), "");
    let out = tmp.path().join("db");
    let e = build(&recipe, &out, &no_dbs(), &opts()).unwrap_err();
    assert_eq!(e.code(), "RECIPE_INVALID");
    assert!(e.to_string().contains("human") && e.to_string().contains("mouse"), "{e}");
    assert!(!out.exists());
}

/// Base database: IGHV x3, IGHJ x2, IGKV x1.
fn keep_base(dir: &Path) -> PathBuf {
    put(
        dir,
        "v.fasta",
        ">IGHV1-2*01\nCAGGTGCAGCTG\n>IGHV1-3*01\nCAGGTCCAGCTT\n>IGHV1-8*01\nCAGGTGCAGCTA\n\
         >IGKV1-5*01\nGACATCCAGATG\n",
    );
    put(
        dir,
        "j.fasta",
        ">IGHJ1*01\nGCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTC\n\
         >IGHJ2*01\nCTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTC\n",
    );
    let recipe = put(
        dir,
        "base.toml",
        "name = \"base\"\n[[source]]\nspecies = \"human\"\nv_fasta = [\"v.fasta\"]\nj_fasta = [\"j.fasta\"]\n",
    );
    let out = dir.join("base");
    build(&recipe, &out, &no_dbs(), &opts()).unwrap();
    out
}

#[test]
fn keep_list_restricts_only_named_pairs() {
    let tmp = tempfile::tempdir().unwrap();
    let base = keep_base(tmp.path());
    put(tmp.path(), "genotype.txt", "# donor 1\nIGHV1-3*01  # comment\n\n");
    let recipe = put(
        tmp.path(),
        "donor.toml",
        "name = \"donor\"\nkeep = \"genotype.txt\"\n[[source]]\ndatabase = \"base\"\n",
    );
    let dbs = BTreeMap::from([("base".to_string(), base)]);
    let out = tmp.path().join("donor");
    let r = build(&recipe, &out, &dbs, &opts()).unwrap();
    assert_eq!(r.restricted, ["IGH/V"]);
    assert_eq!(r.passed_through, ["IGH/J", "IGK/V"]);
    assert_eq!(r.counts.get("IGH/V/human"), Some(&1));
    assert_eq!(r.counts.get("IGH/J/human"), Some(&2));
    assert_eq!(r.counts.get("IGK/V/human"), Some(&1));
    let db = load(&out).unwrap();
    assert_eq!(db.alleles.len(), 4);
    let v = allele(&db, "IGHV1-3*01");
    assert_eq!(v.source, "database:base");
    assert_eq!(db.manifest.provenance[0].files, ["database:base"]);

    // without a map entry, the value is used as a directory path
    let by_path = put(
        tmp.path(),
        "donor2.toml",
        &format!("name = \"donor2\"\n[[source]]\ndatabase = {:?}\n", tmp.path().join("base")),
    );
    let r2 = build(&by_path, &tmp.path().join("donor2"), &no_dbs(), &opts()).unwrap();
    assert_eq!(r2.n_alleles, 6);

    let missing = put(tmp.path(), "bad.toml", "name = \"bad\"\n[[source]]\ndatabase = \"nope\"\n");
    let e = build(&missing, &tmp.path().join("bad"), &no_dbs(), &opts()).unwrap_err();
    assert_eq!(e.code(), "RECIPE_INVALID");
    assert!(e.to_string().contains("nope"), "{e}");
}

#[test]
fn keep_list_unknown_entry_fails() {
    let tmp = tempfile::tempdir().unwrap();
    let base = keep_base(tmp.path());
    put(tmp.path(), "genotype.txt", "IGHV1-2*01\nIGHV9-9*01\nIGHV7-7*01\n");
    let recipe = put(
        tmp.path(),
        "donor.toml",
        "name = \"donor\"\nkeep = \"genotype.txt\"\n[[source]]\ndatabase = \"base\"\n",
    );
    let dbs = BTreeMap::from([("base".to_string(), base)]);
    let out = tmp.path().join("donor");
    let e = build(&recipe, &out, &dbs, &opts()).unwrap_err();
    assert_eq!(e.code(), "RECIPE_INVALID");
    let msg = e.to_string();
    assert!(msg.contains("IGHV9-9*01") && msg.contains("IGHV7-7*01"), "{msg}");
    assert!(!msg.contains("IGHV1-2*01"), "{msg}");
    assert!(!out.exists());
}

fn one_v(dir: &Path, fasta: &str) -> PathBuf {
    put(dir, "v.fasta", fasta);
    put(dir, "r.toml", "name = \"x\"\n[[source]]\nspecies = \"human\"\nv_fasta = [\"v.fasta\"]\n")
}

#[test]
fn fasta_header_imgt_style() {
    let tmp = tempfile::tempdir().unwrap();
    let recipe = one_v(
        tmp.path(),
        ">X62106|IGHV1-2*01|Homo sapiens|(F)|V-REGION|...\nCAGGTGCAGCTG\n",
    );
    let out = tmp.path().join("db");
    build(&recipe, &out, &no_dbs(), &opts()).unwrap();
    let db = load(&out).unwrap();
    let a = allele(&db, "IGHV1-2*01");
    assert_eq!(a.name, "IGHV1-2*01");
    assert_eq!(a.gene, "IGHV1-2");
    assert_eq!(a.functionality, Functionality::F);
}

#[test]
fn unknown_locus_fails() {
    let tmp = tempfile::tempdir().unwrap();
    let recipe = one_v(tmp.path(), ">XYZV1*01\nCAGGTGCAGCTG\n");
    let out = tmp.path().join("db");
    let e = build(&recipe, &out, &no_dbs(), &opts()).unwrap_err();
    assert_eq!(e.code(), "RECIPE_INVALID");
    assert!(e.to_string().contains("XYZV1*01"), "{e}");
    assert!(!out.exists());
}

#[test]
fn failed_build_leaves_no_directory() {
    let tmp = tempfile::tempdir().unwrap();
    let recipe = one_v(tmp.path(), ">IGHV1-2*01\nCAGGTGCAGCTG\n>IGHV1-2*01\nCAGGTGCAGCTA\n");
    let out = tmp.path().join("db");
    let e = build(&recipe, &out, &no_dbs(), &opts()).unwrap_err();
    assert!(e.to_string().contains("IGHV1-2*01"), "{e}");
    assert!(!out.exists());
    let mut left: Vec<_> = fs::read_dir(tmp.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .collect();
    left.sort();
    assert_eq!(left, ["r.toml", "v.fasta"]);
}
