//! Context-dependent somatic hypermutation (record 0024).
//!
//! A [`ContextModel`] gives each ACGT 5-mer a mutability and a substitution row for
//! its centre base. [`motif_model`] is abstar's own heuristic default; [`load_s5f_dir`]
//! reads S5F tables a user downloaded from the Kleinstein lab (never bundled).

use crate::GenerateError;
use rand::Rng;
use std::path::Path;
use std::sync::OnceLock;

const BASES: [u8; 4] = *b"ACGT";

/// 5-mer (ACGT only) → mutability, and 5-mer → [A, C, G, T] substitution probabilities (self = 0).
#[derive(Debug)]
pub struct ContextModel {
    mutability: Vec<f64>,
    substitution: Vec<[f64; 4]>,
    /// Used for 5-mers containing N or an IUPAC code (including the padding at sequence ends).
    mean_mutability: f64,
    /// Mean substitution row for each centre base, for the same 5-mers.
    mean_substitution: [[f64; 4]; 4],
}

/// One hit: 0-based read position, base before and after, at the time of the hit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mutation {
    pub pos: usize,
    pub from: u8,
    pub to: u8,
}

fn code(b: u8) -> Option<usize> {
    BASES.iter().position(|&x| x == b.to_ascii_uppercase())
}

/// Index of an ACGT 5-mer (case-insensitive) in 0..1024.
fn index(kmer: &[u8]) -> Option<usize> {
    if kmer.len() != 5 {
        return None;
    }
    kmer.iter().try_fold(0, |acc, &b| Some(acc * 4 + code(b)?))
}

impl ContextModel {
    fn new(mutability: Vec<f64>, substitution: Vec<[f64; 4]>) -> Self {
        let mean_mutability = mutability.iter().sum::<f64>() / 1024.0;
        let mut mean_substitution = [[0.0; 4]; 4];
        for (i, row) in substitution.iter().enumerate() {
            let centre = (i >> 4) & 3;
            for (m, p) in mean_substitution[centre].iter_mut().zip(row) {
                *m += p / 256.0;
            }
        }
        Self { mutability, substitution, mean_mutability, mean_substitution }
    }

    /// Mutability of an ACGT 5-mer; `None` for anything else.
    pub fn mutability(&self, kmer: &[u8]) -> Option<f64> {
        index(kmer).map(|i| self.mutability[i])
    }

    /// Substitution probabilities [A, C, G, T] for the centre base of an ACGT 5-mer.
    pub fn substitution(&self, kmer: &[u8]) -> Option<[f64; 4]> {
        index(kmer).map(|i| self.substitution[i])
    }

    /// Mutability and substitution row at `pos`, padding the ends with N.
    /// `None` when the centre base is not ACGT (never mutated).
    fn at(&self, seq: &[u8], pos: usize) -> Option<(f64, [f64; 4])> {
        let centre = code(seq[pos])?;
        let kmer: Vec<u8> = (pos as isize - 2..=pos as isize + 2)
            .map(|j| usize::try_from(j).ok().and_then(|j| seq.get(j)).copied().unwrap_or(b'N'))
            .collect();
        Some(match index(&kmer) {
            Some(i) => (self.mutability[i], self.substitution[i]),
            None => (self.mean_mutability, self.mean_substitution[centre]),
        })
    }
}

/// abstar's heuristic motif weight for an ACGT 5-mer (record 0024: not fit to data).
fn motif_weight(k: &[u8]) -> f64 {
    let w = |b: u8| b == b'A' || b == b'T';
    let r = |b: u8| b == b'A' || b == b'G';
    let y = |b: u8| b == b'C' || b == b'T';
    let s = |b: u8| b == b'C' || b == b'G';
    // hot and cold motifs on the same centre are mutually exclusive, so the first
    // match in strongest-first order is the strongest match
    match k[2] {
        b'C' if w(k[0]) && r(k[1]) && y(k[3]) => 6.0, // WRCY
        b'G' if r(k[1]) && y(k[3]) && w(k[4]) => 6.0, // RGYW
        b'C' if w(k[0]) && r(k[1]) => 4.0,            // WRC
        b'G' if y(k[3]) && w(k[4]) => 4.0,            // GYW
        b'A' if w(k[1]) => 2.5,                       // WA
        b'T' if w(k[3]) => 2.5,                       // TW
        b'C' if s(k[0]) && y(k[1]) => 0.3,            // SYC
        b'G' if r(k[3]) && s(k[4]) => 0.3,            // GRS
        _ => 1.0,
    }
}

/// Transition-biased substitution row [A, C, G, T] for a centre base.
fn motif_substitution(centre: u8) -> [f64; 4] {
    match centre {
        b'A' => [0.0, 0.25, 0.5, 0.25],
        b'C' => [0.2, 0.0, 0.2, 0.6],
        b'G' => [0.6, 0.2, 0.0, 0.2],
        _ => [0.25, 0.5, 0.25, 0.0],
    }
}

/// abstar's own default context model, built in code once.
pub fn motif_model() -> &'static ContextModel {
    static MODEL: OnceLock<ContextModel> = OnceLock::new();
    MODEL.get_or_init(|| {
        let kmers: Vec<[u8; 5]> = (0..1024).map(|i| std::array::from_fn(|j| BASES[(i >> (2 * (4 - j))) & 3])).collect();
        ContextModel::new(
            kmers.iter().map(|k| motif_weight(k)).collect(),
            kmers.iter().map(|k| motif_substitution(k[2])).collect(),
        )
    })
}

fn invalid<T>(file: &Path, msg: impl std::fmt::Display) -> Result<T, GenerateError> {
    Err(GenerateError::InvalidArgument(format!("S5F table {}: {msg}", file.display())))
}

/// Rows of a Kleinstein-lab table: whitespace- or comma-separated, values optionally
/// quoted, header first. Returns, per data row, the values of `columns` (looked up by name).
fn read_table(file: &Path, columns: &[&str]) -> Result<Vec<Vec<String>>, GenerateError> {
    let text = match std::fs::read_to_string(file) {
        Ok(t) => t,
        Err(e) => return invalid(file, e),
    };
    let split = |line: &str| -> Vec<String> {
        line.split(|c: char| c == ',' || c.is_whitespace())
            .filter(|f| !f.is_empty())
            .map(|f| f.trim_matches('"').to_string())
            .collect()
    };
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header = split(lines.next().unwrap_or(""));
    let idx = columns
        .iter()
        .map(|c| header.iter().position(|h| h == c).map_or_else(|| invalid(file, format!("no '{c}' column")), Ok))
        .collect::<Result<Vec<_>, _>>()?;
    lines
        .enumerate()
        .map(|(n, line)| {
            let f = split(line);
            idx.iter()
                .map(|&i| f.get(i).cloned().map_or_else(|| invalid(file, format!("row {} is short", n + 2)), Ok))
                .collect()
        })
        .collect()
}

/// Fills a 1024-entry table from rows keyed by ACGT 5-mer; skips N-containing 5-mers and
/// rejects bad, duplicate or missing ones.
fn fill<T: Copy>(
    file: &Path,
    rows: Vec<Vec<String>>,
    parse: impl Fn(&[String], usize) -> Result<T, String>,
) -> Result<Vec<T>, GenerateError> {
    let mut table: Vec<Option<T>> = vec![None; 1024];
    for row in rows {
        let Some(i) = index(row[0].as_bytes()) else {
            // N-containing rows (some exports have them) are ignored: abstar uses the mean
            if row[0].len() == 5 && row[0].bytes().all(|b| b"ACGTNacgtn".contains(&b)) {
                continue;
            }
            return invalid(file, format!("bad 5-mer '{}'", row[0]));
        };
        if table[i].is_some() {
            return invalid(file, format!("duplicate 5-mer {}", row[0]));
        }
        table[i] = Some(parse(&row[1..], i).or_else(|e| invalid(file, format!("{}: {e}", row[0])))?);
    }
    let n = table.iter().filter(|x| x.is_some()).count();
    if n < 1024 {
        return invalid(file, format!("incomplete: {n} of 1024 ACGT 5-mers"));
    }
    Ok(table.into_iter().flatten().collect())
}

fn prob(s: &str) -> Result<f64, String> {
    match s.parse::<f64>() {
        Ok(p) if p.is_finite() && p >= 0.0 => Ok(p),
        _ => Err(format!("'{s}' is not a non-negative number")),
    }
}

/// Loads S5F tables the user downloaded: `DIR/Mutability.csv` (`Fivemer`, `Mutability`)
/// and `DIR/Substitution.csv` (`Fivemer`, `A`, `C`, `G`, `T`), columns found by header
/// name. Anything missing, malformed or incomplete is `INVALID_ARGUMENT` naming the file.
pub fn load_s5f_dir(dir: &Path) -> Result<ContextModel, GenerateError> {
    let mu_file = dir.join("Mutability.csv");
    let mutability = fill(&mu_file, read_table(&mu_file, &["Fivemer", "Mutability"])?, |v, _| prob(&v[0]))?;
    if mutability.iter().all(|&m| m == 0.0) {
        return invalid(&mu_file, "all mutabilities are 0");
    }
    let sub_file = dir.join("Substitution.csv");
    let rows = read_table(&sub_file, &["Fivemer", "A", "C", "G", "T"])?;
    let substitution = fill(&sub_file, rows, |v, i| {
        let row = [prob(&v[0])?, prob(&v[1])?, prob(&v[2])?, prob(&v[3])?];
        if row[(i >> 4) & 3] != 0.0 || (row.iter().sum::<f64>() - 1.0).abs() > 1e-6 {
            return Err("row must have self = 0 and sum to 1".into());
        }
        Ok(row)
    })?;
    Ok(ContextModel::new(mutability, substitution))
}

/// Samples an index ∝ `weights` (assumes a positive total).
fn pick(weights: &[f64], rng: &mut impl Rng) -> usize {
    let mut u = rng.random::<f64>() * weights.iter().sum::<f64>();
    let last = weights.iter().rposition(|&w| w > 0.0).unwrap_or(0);
    weights.iter().position(|&w| w > 0.0 && { u -= w; u < 0.0 }).unwrap_or(last)
}

/// Mutates `seq` in place, one hit at a time, until the realized divergence of the
/// `v_positions` from their original bases reaches `round(target_v_freq × v_positions.len())`.
///
/// - Only positions with `mutable[pos]` and an A/C/G/T base (either case) are hit; a
///   lowercase base is replaced by a lowercase base.
/// - Positions are picked ∝ the mutability of their 5-mer, recomputed after every hit;
///   repeat hits and back-mutations are allowed.
/// - Safety cap: at most 10 × target hits. A capped run returns with the V divergence
///   still below target; the caller, which records realized divergence anyway, detects
///   that and emits the warning (this function has no warning channel).
///
/// Panics if `mutable.len() != seq.len()` or a V position is out of range (caller bugs).
pub fn mutate(
    seq: &mut [u8],
    mutable: &[bool],
    v_positions: &[usize],
    target_v_freq: f64,
    model: &ContextModel,
    rng: &mut impl Rng,
) -> Vec<Mutation> {
    assert_eq!(mutable.len(), seq.len(), "mutable mask length");
    let target = (target_v_freq * v_positions.len() as f64).round() as usize;
    let mut v_orig: Vec<Option<u8>> = vec![None; seq.len()];
    for &p in v_positions {
        v_orig[p] = Some(seq[p]);
    }
    let weight = |seq: &[u8], i: usize| if mutable[i] { model.at(seq, i).map_or(0.0, |(m, _)| m) } else { 0.0 };
    let mut weights: Vec<f64> = (0..seq.len()).map(|i| weight(seq, i)).collect();
    let (mut hits, mut diverged) = (Vec::new(), 0);
    // ponytail: O(len) scan per hit; a Fenwick tree if profiling ever shows it
    while diverged < target && hits.len() < 10 * target && weights.iter().any(|&w| w > 0.0) {
        let pos = pick(&weights, rng);
        let (_, subs) = model.at(seq, pos).expect("positive weight implies an ACGT centre");
        let from = seq[pos];
        let mut to = BASES[pick(&subs, rng)];
        if from.is_ascii_lowercase() {
            to = to.to_ascii_lowercase();
        }
        if let Some(o) = v_orig[pos] {
            diverged = diverged + usize::from(to != o) - usize::from(from != o);
        }
        seq[pos] = to;
        hits.push(Mutation { pos, from, to });
        for (i, w) in weights.iter_mut().enumerate().take(pos + 3).skip(pos.saturating_sub(2)) {
            *w = weight(seq, i);
        }
    }
    hits
}
