use crate::io::hex_sha256;
use crate::{
    anchors, gap, load, recipe, sources, validate, write, Allele, BuildWarning, Database, GermlineError,
    Locus, Manifest, Provenance, Segment,
};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

pub struct BuildOptions {
    pub overwrite: bool,
    pub abstar_version: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BuildReport {
    pub name: String,
    pub path: PathBuf,
    pub n_alleles: usize,
    /// key "LOCUS/SEGMENT/species", e.g. "IGH/V/human"
    pub counts: BTreeMap<String, usize>,
    /// "LOCUS/SEGMENT" keys restricted by the keep-list / present but untouched by it
    pub restricted: Vec<String>,
    pub passed_through: Vec<String>,
    pub warnings: Vec<BuildWarning>,
}

fn invalid(m: String) -> GermlineError {
    GermlineError::RecipeInvalid(m)
}

/// Reads `base/rel`, returning (basename, bytes); a missing file names its path.
fn read(base: &Path, rel: &Path) -> Result<(String, Vec<u8>), GermlineError> {
    let path = base.join(rel);
    let bytes = fs::read(&path).map_err(|e| invalid(format!("{}: {e}", path.display())))?;
    let name = path.file_name().unwrap_or_default().to_string_lossy().into_owned();
    Ok((name, bytes))
}

/// An allele plus its anchor-search warning.
type Entry = (Allele, Option<String>);

pub fn build(
    recipe_path: &Path,
    output_dir: &Path,
    database_paths: &BTreeMap<String, PathBuf>,
    opts: &BuildOptions,
) -> Result<BuildReport, GermlineError> {
    let recipe = recipe::read(recipe_path)?;
    let base = recipe_path.parent().unwrap_or(Path::new("."));
    let mut entries: Vec<Entry> = Vec::new();
    let mut provenance = Vec::new();
    // gap_v outcomes as (species, allele name, message); ids are final only after the prefix step
    let mut gap_notes: Vec<(String, String, String)> = Vec::new();

    for src in &recipe.sources {
        if let Some(name) = &src.database {
            let path = match database_paths.get(name) {
                Some(p) => p.clone(),
                None if base.join(name).is_dir() => base.join(name),
                None => return Err(invalid(format!("database {name:?} not found"))),
            };
            let db = load(&path)?;
            let label = format!("database:{}", db.manifest.name);
            for a in db.alleles {
                // anchors come from the database as built; no new search
                entries.push((Allele { source: label.clone(), ..a }, None));
            }
            provenance.push(Provenance {
                species: db.manifest.species.join("+"),
                files: vec![label],
                sha256: vec![db.manifest.alleles_sha256],
                description: if src.description.is_empty() {
                    db.manifest.description
                } else {
                    src.description.clone()
                },
            });
            continue;
        }
        let species = src.species.clone().unwrap_or_default(); // recipe::read guarantees it
        let mut p = Provenance {
            species: species.clone(),
            files: vec![],
            sha256: vec![],
            description: src.description.clone(),
        };
        let mut descriptions: Vec<String> = Vec::new();
        let mut templates: Vec<Allele> = Vec::new();
        if src.gap_v {
            let names: Vec<String> = src
                .gap_templates
                .clone()
                .unwrap_or_else(|| gap::DEFAULT_TEMPLATES.iter().map(|s| s.to_string()).collect());
            for name in names {
                let path = match database_paths.get(&name) {
                    Some(p) => p.clone(),
                    None if base.join(&name).is_dir() => base.join(&name),
                    None => return Err(invalid(format!("gap template database {name:?} not found"))),
                };
                let db = load(&path)?;
                p.files.push(format!("gap-template:{}", db.manifest.name));
                p.sha256.push(db.manifest.alleles_sha256);
                templates.extend(db.alleles);
            }
        }
        let mut add = |alleles: Vec<Allele>, file: &str| {
            let alleles = if src.gap_v {
                let v: Vec<(String, String)> = alleles
                    .iter()
                    .filter(|a| a.segment == Segment::V)
                    .map(|a| (a.name.clone(), a.sequence()))
                    .collect();
                let mut results = gap::gap_v(&v, &templates).into_iter();
                let mut out = Vec::new();
                for mut a in alleles {
                    if a.segment == Segment::V {
                        let r = results.next().expect("one result per V allele");
                        match (r.gapped, r.failure) {
                            (Some(g), _) => {
                                a.sequence_gapped = g;
                                if !r.warnings.is_empty() {
                                    let m = format!("gapped: {}", r.warnings.join(", "));
                                    gap_notes.push((species.clone(), a.name.clone(), m));
                                }
                            }
                            (None, f) => {
                                let m = format!("not gapped: {}", f.unwrap_or_default());
                                gap_notes.push((species.clone(), a.name.clone(), m));
                                continue;
                            }
                        }
                    }
                    out.push(a);
                }
                out
            } else {
                alleles
            };
            for mut a in alleles {
                a.species = species.clone();
                a.source = file.into();
                let w = anchors::assign(&mut a);
                entries.push((a, w));
            }
        };
        for rel in &src.ogrdb_json {
            let (file, bytes) = read(base, rel)?;
            let (alleles, sets) = sources::ogrdb(&bytes, &file)?;
            add(alleles, &file);
            descriptions.extend(sets);
            p.sha256.push(hex_sha256(&bytes));
            p.files.push(file);
        }
        for (segment, rel) in src.fasta() {
            let (file, bytes) = read(base, rel)?;
            add(sources::fasta(&bytes, &file, segment)?, &file);
            p.sha256.push(hex_sha256(&bytes));
            p.files.push(file);
        }
        if !descriptions.is_empty() {
            let mut parts: Vec<String> = Some(p.description).into_iter().filter(|d| !d.is_empty()).collect();
            parts.extend(descriptions);
            p.description = parts.join("; ");
        }
        provenance.push(p);
    }

    // ids: plain names for one species, "{prefix}|{name}" otherwise (re-prefixing database alleles)
    let all_species: BTreeSet<String> = entries.iter().map(|(a, _)| a.species.clone()).collect();
    let mut prefixes = BTreeMap::new();
    if all_species.len() > 1 {
        let missing: Vec<&str> =
            all_species.iter().filter(|s| !recipe.prefixes.contains_key(*s)).map(String::as_str).collect();
        if !missing.is_empty() {
            return Err(invalid(format!(
                "database mixes species {}; [prefixes] has no entry for {}",
                all_species.iter().cloned().collect::<Vec<_>>().join(", "),
                missing.join(", ")
            )));
        }
        for s in &all_species {
            prefixes.insert(s.clone(), recipe.prefixes[s].clone());
        }
        let distinct: HashSet<&String> = prefixes.values().collect();
        if distinct.len() != prefixes.len() {
            return Err(invalid(format!("[prefixes] values must be distinct: {prefixes:?}")));
        }
    }
    for (a, _) in &mut entries {
        a.id = match prefixes.get(&a.species) {
            Some(p) => format!("{p}|{}", a.name),
            None => a.name.clone(),
        };
    }

    let (mut restricted, mut passed_through) = (Vec::new(), Vec::new());
    if let Some(rel) = &recipe.keep {
        let restrict = apply_keep(base, rel, &mut entries)?;
        let present: BTreeSet<(Locus, Segment)> = entries.iter().map(|(a, _)| (a.locus, a.segment)).collect();
        for &(l, s) in &present {
            let key = format!("{l}/{s}");
            if restrict.contains(&(l, s)) { restricted.push(key) } else { passed_through.push(key) }
        }
    }

    let mut warnings: Vec<BuildWarning> = entries
        .iter()
        .filter_map(|(a, w)| w.as_ref().map(|m| BuildWarning { allele: a.id.clone(), message: m.clone() }))
        .collect();
    warnings.extend(gap_notes.into_iter().map(|(species, name, message)| BuildWarning {
        allele: match prefixes.get(&species) {
            Some(p) => format!("{p}|{name}"),
            None => name,
        },
        message,
    }));
    let alleles: Vec<Allele> = entries.into_iter().map(|(a, _)| a).collect();
    // validate's generic "no anchor" warning is dropped for alleles that already have
    // a specific anchor-search warning: one message per allele
    let searched: HashSet<String> = warnings.iter().map(|w| w.allele.clone()).collect();
    warnings.extend(validate(&alleles)?.into_iter().filter(|w| !searched.contains(&w.allele)));

    let db = Database {
        manifest: Manifest {
            format_version: 1,
            name: recipe.name.clone(),
            description: recipe.description,
            abstar_version: opts.abstar_version.clone(),
            built: date(now_days()),
            alleles_sha256: String::new(), // set by write
            species: alleles.iter().map(|a| a.species.clone()).collect::<BTreeSet<_>>().into_iter().collect(),
            loci: alleles.iter().map(|a| a.locus).collect::<BTreeSet<_>>().into_iter().collect(),
            prefixes,
            provenance,
            warnings: warnings.clone(),
        },
        alleles,
    };
    write(&db, output_dir, opts.overwrite)?;
    Ok(BuildReport {
        name: recipe.name,
        path: output_dir.to_path_buf(),
        n_alleles: db.alleles.len(),
        counts: db.counts(),
        restricted,
        passed_through,
        warnings,
    })
}

/// Applies the keep-list: drops unkept alleles in every (locus, segment) pair it names.
/// Returns those pairs.
fn apply_keep(base: &Path, rel: &Path, entries: &mut Vec<Entry>) -> Result<BTreeSet<(Locus, Segment)>, GermlineError> {
    let (_, bytes) = read(base, rel)?;
    let text = String::from_utf8_lossy(&bytes);
    let mut kept = HashSet::new();
    let mut offenders = Vec::new();
    for line in text.lines() {
        let entry = line.split('#').next().unwrap_or("").trim();
        if entry.is_empty() {
            continue;
        }
        if let Some(i) = entries.iter().position(|(a, _)| a.id == entry) {
            kept.insert(i);
            continue;
        }
        let by_name: Vec<usize> = (0..entries.len()).filter(|&i| entries[i].0.name == entry).collect();
        match by_name.as_slice() {
            [i] => {
                kept.insert(*i);
            }
            [] => offenders.push(format!("{entry} (unknown)")),
            many => offenders.push(format!(
                "{entry} (ambiguous: {})",
                many.iter().map(|&i| entries[i].0.id.as_str()).collect::<Vec<_>>().join(", ")
            )),
        }
    }
    if !offenders.is_empty() {
        return Err(invalid(format!("keep-list {}: {}", base.join(rel).display(), offenders.join("; "))));
    }
    let restrict: BTreeSet<(Locus, Segment)> =
        kept.iter().map(|&i| (entries[i].0.locus, entries[i].0.segment)).collect();
    let mut i = 0;
    entries.retain(|(a, _)| {
        let keep = kept.contains(&i) || !restrict.contains(&(a.locus, a.segment));
        i += 1;
        keep
    });
    Ok(restrict)
}

fn now_days() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs() / 86_400)
}

/// Days since 1970-01-01 to a UTC YYYY-MM-DD (Hinnant's civil_from_days).
fn date(days: u64) -> String {
    let z = days as i64 + 719_468;
    let era = z / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!("{y:04}-{m:02}-{d:02}")
}

#[cfg(test)]
mod tests {
    #[test]
    fn date_from_days() {
        assert_eq!(super::date(0), "1970-01-01");
        assert_eq!(super::date(11_016), "2000-02-29");
        assert_eq!(super::date(20_777), "2026-11-20");
    }
}
