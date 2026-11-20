use crate::{validate, Allele, Database, GermlineError, Manifest};
use sha2::{Digest, Sha256};
use std::fs;
use std::path::{Path, PathBuf};

const HEADER: &str =
    "id\tname\tgene\tlocus\tsegment\tspecies\tfunctionality\tsequence_gapped\tanchor\tsource\tusage";

pub(crate) fn hex_sha256(bytes: &[u8]) -> String {
    Sha256::digest(bytes).iter().map(|b| format!("{b:02x}")).collect()
}

fn opt<T: ToString>(v: &Option<T>) -> String {
    v.as_ref().map(T::to_string).unwrap_or_default()
}

fn render(alleles: &[Allele]) -> String {
    let mut rows: Vec<&Allele> = alleles.iter().collect();
    rows.sort_by(|a, b| {
        (a.locus, a.segment, a.id.as_bytes()).cmp(&(b.locus, b.segment, b.id.as_bytes()))
    });
    let mut out = format!("{HEADER}\n");
    for a in rows {
        let cells = [
            a.id.clone(), a.name.clone(), a.gene.clone(), a.locus.to_string(),
            a.segment.to_string(), a.species.clone(), a.functionality.to_string(),
            a.sequence_gapped.clone(), opt(&a.anchor), a.source.clone(), opt(&a.usage),
        ];
        out.push_str(&cells.join("\t"));
        out.push('\n');
    }
    out
}

fn parse(text: &str) -> Result<Vec<Allele>, GermlineError> {
    let mut lines = text.lines();
    if lines.next() != Some(HEADER) {
        return Err(GermlineError::Invalid("alleles.tsv: bad header".into()));
    }
    let mut out = Vec::new();
    for (i, line) in lines.enumerate() {
        let row = i + 1;
        let err = |m: String| GermlineError::Invalid(format!("alleles.tsv row {row}: {m}"));
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() != 11 {
            return Err(err(format!("expected 11 cells, found {}", c.len())));
        }
        fn some(s: &str) -> Option<&str> {
            (!s.is_empty()).then_some(s)
        }
        let a = Allele {
            id: c[0].into(),
            name: c[1].into(),
            gene: c[2].into(),
            locus: c[3].parse().map_err(|e: GermlineError| err(e.to_string()))?,
            segment: c[4].parse().map_err(|e: GermlineError| err(e.to_string()))?,
            species: c[5].into(),
            functionality: c[6].parse().map_err(|e: GermlineError| err(e.to_string()))?,
            sequence_gapped: c[7].into(),
            anchor: some(c[8])
                .map(|s| s.parse().map_err(|_| err(format!("bad anchor {s:?}"))))
                .transpose()?,
            source: c[9].into(),
            usage: some(c[10])
                .map(|s| s.parse().map_err(|_| err(format!("bad usage {s:?}"))))
                .transpose()?,
        };
        if let Some(p) = a.anchor {
            let len = a.sequence().len();
            if p >= len {
                return Err(err(format!("anchor {p} out of range for length {len}")));
            }
        }
        out.push(a);
    }
    Ok(out)
}

/// Loads and validates a database; every error message names `dir`.
pub fn load(dir: &Path) -> Result<Database, GermlineError> {
    if !dir.is_dir() {
        return Err(GermlineError::NotFound(dir.into()));
    }
    load_dir(dir).map_err(|e| match e {
        GermlineError::Invalid(m) => GermlineError::Invalid(format!("{}: {m}", dir.display())),
        GermlineError::Io(e) => {
            GermlineError::Io(std::io::Error::new(e.kind(), format!("{}: {e}", dir.display())))
        }
        e => e,
    })
}

fn load_dir(dir: &Path) -> Result<Database, GermlineError> {
    let read = |name: &str| {
        fs::read(dir.join(name)).map_err(|e| match e.kind() {
            std::io::ErrorKind::NotFound => GermlineError::Invalid(format!("missing {name}")),
            _ => e.into(),
        })
    };
    let manifest_text = read("manifest.toml")?;
    let manifest: Manifest = toml::from_str(&String::from_utf8_lossy(&manifest_text))
        .map_err(|e| GermlineError::Invalid(format!("manifest.toml: {e}")))?;
    let bytes = read("alleles.tsv")?;
    if hex_sha256(&bytes) != manifest.alleles_sha256 {
        return Err(GermlineError::Invalid("alleles.tsv: sha256 does not match manifest".into()));
    }
    let text = std::str::from_utf8(&bytes)
        .map_err(|e| GermlineError::Invalid(format!("alleles.tsv: {e}")))?;
    let alleles = parse(text)?;
    validate(&alleles)?;
    Ok(Database { manifest, alleles })
}

fn sibling(dir: &Path, suffix: &str) -> PathBuf {
    let mut n = dir.file_name().unwrap_or_default().to_os_string();
    n.push(format!(".{suffix}-{}", std::process::id()));
    dir.with_file_name(n)
}

pub fn write(db: &Database, dir: &Path, overwrite: bool) -> Result<(), GermlineError> {
    if dir.exists() {
        if !overwrite {
            return Err(GermlineError::Exists(dir.into()));
        }
        // overwrite replaces a database, never an arbitrary directory
        if !dir.join("manifest.toml").is_file() {
            return Err(GermlineError::NotADatabase(dir.into()));
        }
    }
    validate(&db.alleles)?; // never write a file load would reject
    let staging = sibling(dir, "staging");
    let res = stage_and_swap(db, dir, &staging);
    if res.is_err() {
        let _ = fs::remove_dir_all(&staging);
    }
    res
}

fn stage_and_swap(db: &Database, dir: &Path, staging: &Path) -> Result<(), GermlineError> {
    let _ = fs::remove_dir_all(staging);
    let tsv = render(&db.alleles);
    let mut manifest = db.manifest.clone();
    manifest.alleles_sha256 = hex_sha256(tsv.as_bytes());
    let toml = toml::to_string(&manifest)
        .map_err(|e| GermlineError::Invalid(format!("manifest: {e}")))?;
    if let Some(parent) = staging.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::create_dir(staging)?;
    fs::write(staging.join("alleles.tsv"), tsv)?;
    fs::write(staging.join("manifest.toml"), toml)?;
    if dir.exists() {
        let old = sibling(dir, "old");
        fs::rename(dir, &old)?;
        if let Err(e) = fs::rename(staging, dir) {
            let _ = fs::rename(&old, dir);
            return Err(e.into());
        }
        let _ = fs::remove_dir_all(&old);
    } else {
        fs::rename(staging, dir)?;
    }
    Ok(())
}
