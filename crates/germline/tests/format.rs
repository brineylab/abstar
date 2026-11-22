use abstar_germline::*;
use sha2::Digest;
use std::collections::BTreeMap;
use std::fs;
use std::os::unix::fs::PermissionsExt;

fn allele(id: &str, locus: Locus, segment: Segment, seq: &str, anchor: Option<usize>) -> Allele {
    Allele {
        id: id.into(),
        name: format!("{id}*01"),
        gene: id.into(),
        locus,
        segment,
        species: "human".into(),
        functionality: Functionality::F,
        sequence_gapped: seq.into(),
        anchor,
        source: "test".into(),
        usage: None,
    }
}

/// A gapped V sequence with Cys104 (TGT) at gapped 309..312, then `tail`.
fn cys_seq(tail: &str) -> String {
    format!("A...{}TGT{tail}", "A".repeat(305))
}

fn db() -> Database {
    Database {
        manifest: Manifest {
            format_version: 1,
            name: "t".into(),
            description: "test db".into(),
            abstar_version: "1.0.0".into(),
            built: "2026-11-20".into(),
            alleles_sha256: String::new(),
            species: vec!["human".into()],
            loci: vec![Locus::IGH],
            prefixes: BTreeMap::new(),
            provenance: vec![Provenance {
                species: "human".into(),
                files: vec!["a.fasta".into()],
                sha256: vec!["abc".into()],
                description: "d".into(),
            }],
            warnings: vec![],
        },
        // deliberately unsorted
        alleles: vec![
            allele("C1", Locus::IGH, Segment::C, "ACGT", None),
            allele("J1", Locus::IGH, Segment::J, "ACGTACGT", Some(2)),
            allele("D1", Locus::IGH, Segment::D, "AC.GT", None),
            {
                let mut v = allele("V1", Locus::IGH, Segment::V, &cys_seq("ACGTAC..GTAC"), Some(3));
                v.usage = Some(0.25);
                v.functionality = Functionality::Unknown;
                v
            },
        ],
    }
}

const HEADER: &str =
    "id\tname\tgene\tlocus\tsegment\tspecies\tfunctionality\tsequence_gapped\tanchor\tsource\tusage";

#[test]
fn round_trip() {
    let t = tempfile::tempdir().unwrap();
    let dir = t.path().join("db");
    let d = db();
    write(&d, &dir, false).unwrap();
    let tsv = fs::read_to_string(dir.join("alleles.tsv")).unwrap();
    let lines: Vec<&str> = tsv.lines().collect();
    assert_eq!(lines[0], HEADER);
    let ids: Vec<&str> = lines[1..].iter().map(|l| l.split('\t').next().unwrap()).collect();
    assert_eq!(ids, ["V1", "D1", "J1", "C1"]);
    let loaded = load(&dir).unwrap();
    assert_eq!(loaded.manifest.built, "2026-11-20");
    assert_eq!(loaded.manifest.alleles_sha256.len(), 64);
    let mut expect = d;
    expect.manifest.alleles_sha256 = loaded.manifest.alleles_sha256.clone();
    expect.alleles.sort_by_key(|a| ["V", "D", "J", "C"].iter().position(|s| *s == a.segment.to_string()));
    assert_eq!(loaded, expect);
}

#[test]
fn hash_mismatch_is_invalid() {
    let t = tempfile::tempdir().unwrap();
    let dir = t.path().join("db");
    write(&db(), &dir, false).unwrap();
    let p = dir.join("alleles.tsv");
    let s = fs::read_to_string(&p).unwrap().replacen("test", "tesu", 1);
    fs::write(&p, s).unwrap();
    assert_eq!(load(&dir).unwrap_err().code(), "DATABASE_INVALID");
}

#[test]
fn missing_dir_not_found() {
    let t = tempfile::tempdir().unwrap();
    assert_eq!(load(&t.path().join("nope")).unwrap_err().code(), "DATABASE_NOT_FOUND");
}

fn load_hash(p: &std::path::Path) -> String {
    let t = fs::read_to_string(p).unwrap();
    let l = t.lines().find(|l| l.starts_with("alleles_sha256")).unwrap();
    l.split('"').nth(1).unwrap().to_string()
}

#[test]
fn bad_cell_message_has_row() {
    let t = tempfile::tempdir().unwrap();
    let dir = t.path().join("db");
    write(&db(), &dir, false).unwrap();
    let p = dir.join("alleles.tsv");
    // corrupt a locus on data row 2 and fix up the hash so the cell error is what fires
    let s = fs::read_to_string(&p).unwrap();
    let mut lines: Vec<String> = s.lines().map(String::from).collect();
    lines[2] = lines[2].replacen("IGH", "XYZ", 1);
    let new = lines.join("\n") + "\n";
    fs::write(&p, &new).unwrap();
    // keep the manifest hash consistent so the cell error is what fires
    let mp = dir.join("manifest.toml");
    let old_hash = load_hash(&mp);
    let new_hash: String = sha2::Sha256::digest(new.as_bytes()).iter().map(|b| format!("{b:02x}")).collect();
    fs::write(&mp, fs::read_to_string(&mp).unwrap().replace(&old_hash, &new_hash)).unwrap();
    let e = load(&dir).unwrap_err();
    assert_eq!(e.code(), "DATABASE_INVALID");
    assert!(e.to_string().contains("row 2"), "{e}");
}

#[test]
fn duplicate_id_rejected() {
    let a = allele("X", Locus::IGH, Segment::V, "ACGT", Some(0));
    let e = validate(&[a.clone(), a]).unwrap_err();
    assert_eq!(e.code(), "DATABASE_INVALID");
}

#[test]
fn bad_character_rejected() {
    for bad in ["ACGX", "ACG-T", "ACG*", "ACG1", "ACGU"] {
        let a = allele("X", Locus::IGH, Segment::V, bad, Some(0));
        assert!(validate(&[a]).is_err(), "{bad}");
    }
    let ok = allele("Y", Locus::IGH, Segment::V, "acgn.t", Some(0));
    assert!(validate(&[ok]).is_ok());
}

#[test]
fn empty_and_anchor_out_of_range_rejected() {
    assert!(validate(&[allele("X", Locus::IGH, Segment::D, "..", None)]).is_err());
    assert!(validate(&[allele("X", Locus::IGH, Segment::V, "ACGT", Some(4))]).is_err());
}

#[test]
fn missing_v_anchor_warns() {
    let w = validate(&[allele("V9", Locus::IGH, Segment::V, "ACGT", None)]).unwrap();
    assert_eq!(w.len(), 1);
    assert_eq!(w[0].allele, "V9");
    let w = validate(&[allele("J9", Locus::IGH, Segment::J, "ACGT", None)]).unwrap();
    assert_eq!(w.len(), 1);
}

#[test]
fn existing_dir_without_overwrite() {
    let t = tempfile::tempdir().unwrap();
    let dir = t.path().join("db");
    write(&db(), &dir, false).unwrap();
    assert_eq!(write(&db(), &dir, false).unwrap_err().code(), "DATABASE_EXISTS");
    write(&db(), &dir, true).unwrap();
    load(&dir).unwrap();
    let names: Vec<String> = fs::read_dir(t.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    assert_eq!(names, ["db"]);
}

#[test]
fn failed_write_leaves_nothing() {
    if std::env::var("USER").as_deref() == Ok("root") {
        return; // read-only dirs don't bind root
    }
    let t = tempfile::tempdir().unwrap();
    let parent = t.path().join("ro");
    fs::create_dir(&parent).unwrap();
    fs::set_permissions(&parent, fs::Permissions::from_mode(0o555)).unwrap();
    let e = write(&db(), &parent.join("db"), false).unwrap_err();
    fs::set_permissions(&parent, fs::Permissions::from_mode(0o755)).unwrap();
    assert_eq!(e.code(), "IO_ERROR");
    assert_eq!(fs::read_dir(&parent).unwrap().count(), 0);
}

fn retarget_hash(dir: &std::path::Path, new_tsv: &str) {
    let mp = dir.join("manifest.toml");
    let old = load_hash(&mp);
    let new: String = sha2::Sha256::digest(new_tsv.as_bytes()).iter().map(|b| format!("{b:02x}")).collect();
    fs::write(dir.join("alleles.tsv"), new_tsv).unwrap();
    fs::write(&mp, fs::read_to_string(&mp).unwrap().replace(&old, &new)).unwrap();
}

#[test]
fn anchor_out_of_range_on_load_has_row() {
    let t = tempfile::tempdir().unwrap();
    let dir = t.path().join("db");
    write(&db(), &dir, false).unwrap();
    let s = fs::read_to_string(dir.join("alleles.tsv")).unwrap();
    let mut lines: Vec<String> = s.lines().map(String::from).collect();
    // row 3 is J1 (anchor 2); set to 99
    let mut c: Vec<&str> = lines[3].split('\t').collect();
    c[8] = "99";
    lines[3] = c.join("\t");
    retarget_hash(&dir, &(lines.join("\n") + "\n"));
    let e = load(&dir).unwrap_err();
    assert_eq!(e.code(), "DATABASE_INVALID");
    assert!(e.to_string().contains("row 3"), "{e}");
}

#[test]
fn missing_files_are_invalid() {
    for name in ["manifest.toml", "alleles.tsv"] {
        let t = tempfile::tempdir().unwrap();
        let dir = t.path().join("db");
        write(&db(), &dir, false).unwrap();
        fs::remove_file(dir.join(name)).unwrap();
        let e = load(&dir).unwrap_err();
        assert_eq!(e.code(), "DATABASE_INVALID");
        assert!(e.to_string().contains(name), "{e}");
    }
}

#[test]
fn control_characters_rejected() {
    let t = tempfile::tempdir().unwrap();
    let dir = t.path().join("db");
    for c in ["\t", "\n", "\r"] {
        let mut d = db();
        d.alleles[0].source = format!("x{c}y");
        let e = validate(&d.alleles).unwrap_err();
        assert!(e.to_string().contains("source") && e.to_string().contains("C1"), "{e}");
        assert_eq!(write(&d, &dir, false).unwrap_err().code(), "DATABASE_INVALID");
        assert!(!dir.exists());
    }
}

#[test]
fn iupac_ambiguity_codes_accepted_and_round_trip() {
    // real OGRDB sets carry IUPAC codes, e.g. C57BL/6 IGHV0-KAOG*00 has an S
    let t = tempfile::tempdir().unwrap();
    let dir = t.path().join("db");
    let mut d = db();
    d.alleles.push(allele("V2", Locus::IGH, Segment::V, &cys_seq("ACGTS..sRYWKMBDHVnryswkmbdhv"), Some(4)));
    assert!(validate(&d.alleles).is_ok());
    write(&d, &dir, false).unwrap();
    let loaded = load(&dir).unwrap();
    let v2 = loaded.alleles.iter().find(|a| a.id == "V2").unwrap();
    assert_eq!(v2.sequence_gapped, cys_seq("ACGTS..sRYWKMBDHVnryswkmbdhv"));
}

#[test]
fn overwrite_refuses_a_directory_that_is_not_a_database() {
    let t = tempfile::tempdir().unwrap();
    let dir = t.path().join("victim");
    fs::create_dir(&dir).unwrap();
    fs::write(dir.join("thesis.docx"), "precious").unwrap();
    let e = write(&db(), &dir, true).unwrap_err();
    assert_eq!(e.code(), "DATABASE_EXISTS");
    assert!(e.to_string().contains("not a germline database"), "{e}");
    assert_eq!(fs::read_to_string(dir.join("thesis.docx")).unwrap(), "precious");
    assert_eq!(fs::read_dir(t.path()).unwrap().count(), 1);
}

#[test]
fn load_errors_name_the_directory() {
    let t = tempfile::tempdir().unwrap();
    let dir = t.path().join("db");
    write(&db(), &dir, false).unwrap();
    fs::write(dir.join("manifest.toml"), "not = [toml").unwrap();
    let e = load(&dir).unwrap_err();
    assert_eq!(e.code(), "DATABASE_INVALID");
    assert!(e.to_string().contains(&dir.display().to_string()), "{e}");
}

fn db_with_v(vs: Vec<Allele>) -> Database {
    let mut d = db();
    d.alleles = vs;
    d
}

/// Gapped V whose gapped 309..312 is `codon`.
fn gapped_v(id: &str, codon: &str) -> Allele {
    allele(id, Locus::IGH, Segment::V, &cys_seq("ACGT").replace("TGT", codon), None)
}

#[test]
fn ungapped_v_database_refused() {
    let t = tempfile::tempdir().unwrap();
    let d = db_with_v(vec![allele("V1", Locus::IGH, Segment::V, "ACGTACGT", None)]);
    let e = write(&d, &t.path().join("db"), false).unwrap_err().to_string();
    assert!(e.contains("database t") && e.contains("1 of 1") && e.contains("V1*01") && e.contains("IMGT-gapped"), "{e}");
}

#[test]
fn gapped_v_database_needs_95_percent_cys104() {
    let t = tempfile::tempdir().unwrap();
    let mut vs: Vec<Allele> = (0..19).map(|i| gapped_v(&format!("V{i}"), "TGT")).collect();
    vs.push(gapped_v("V19", "AAA"));
    // 19 of 20 = 95%: accepted, and loads back
    write(&db_with_v(vs.clone()), &t.path().join("ok"), false).unwrap();
    load(&t.path().join("ok")).unwrap();
    vs[18] = gapped_v("V18", "AAA");
    let e = write(&db_with_v(vs), &t.path().join("bad"), false).unwrap_err().to_string();
    assert!(e.contains("2 of 20") && e.contains("V18*01") && e.contains("V19*01"), "{e}");
}

#[test]
fn database_without_v_alleles_skips_the_gap_check() {
    let t = tempfile::tempdir().unwrap();
    let d = db_with_v(vec![allele("J1", Locus::IGH, Segment::J, "ACGT", None)]);
    write(&d, &t.path().join("db"), false).unwrap();
}

#[test]
fn examples_name_ungapped_alleles_first() {
    let t = tempfile::tempdir().unwrap();
    let mut vs: Vec<Allele> = (0..7).map(|i| gapped_v(&format!("A{i}"), "AAA")).collect();
    vs.push(allele("U1", Locus::IGH, Segment::V, "ACGTACGT", None));
    let e = write(&db_with_v(vs), &t.path().join("db"), false).unwrap_err().to_string();
    assert!(e.contains("U1*01"), "{e}");
}

#[test]
fn mixed_database_checked_per_species() {
    let t = tempfile::tempdir().unwrap();
    let mut vs: Vec<Allele> = (0..20).map(|i| gapped_v(&format!("H{i}"), "TGT")).collect();
    for i in 0..2 {
        let mut a = gapped_v(&format!("M{i}"), "AAA");
        a.species = "macaque".into();
        vs.push(a);
    }
    let e = write(&db_with_v(vs), &t.path().join("db"), false).unwrap_err().to_string();
    assert!(e.contains("2 of 2 macaque V alleles") && e.contains("M0*01"), "{e}");
}

#[test]
fn load_refuses_an_invalid_database_directory() {
    let t = tempfile::tempdir().unwrap();
    let dir = t.path().join("db");
    write(&db_with_v(vec![gapped_v("V1", "TGT")]), &dir, false).unwrap();
    // swap in an ungapped V, keeping the manifest hash consistent
    let tsv = fs::read_to_string(dir.join("alleles.tsv")).unwrap();
    let bad = tsv.replace(&cys_seq("ACGT"), "ACGTACGT");
    assert_ne!(tsv, bad);
    let (old, new) = (hex(&tsv), hex(&bad));
    fs::write(dir.join("alleles.tsv"), &bad).unwrap();
    let m = fs::read_to_string(dir.join("manifest.toml")).unwrap();
    assert!(m.contains(&old));
    fs::write(dir.join("manifest.toml"), m.replace(&old, &new)).unwrap();
    let e = load(&dir).unwrap_err().to_string();
    assert!(e.contains("not IMGT-gapped"), "{e}");
}

fn hex(s: &str) -> String {
    sha2::Sha256::digest(s.as_bytes()).iter().map(|b| format!("{b:02x}")).collect()
}
