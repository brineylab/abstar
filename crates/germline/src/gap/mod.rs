//! IMGT gapping of V germline sequences from IMGT-gapped protein templates.
pub mod align;
pub mod imgt;

use crate::{Allele, Locus, Segment};
use align::{align, blosum62, Col, Scoring};
use imgt::{cdr_positions, region, template_residues, translate, Cdr, Region};
use rayon::prelude::*;
use std::cmp::Reverse;

pub const DEFAULT_TEMPLATES: [&str; 3] = ["human", "c57bl6", "balbc"];
pub const TOP_N: usize = 10;
pub const LOW_IDENTITY: f64 = 0.60;
pub const TIE_MARGIN: i32 = 10;

/// Smallest number of query residues a (frame, template) pair must align to be tried. If the
/// best-scoring pair aligns fewer, the query is not a V (`no_alignment`).
const MIN_ALIGNED: usize = 50;

#[derive(Debug, Clone, PartialEq)]
pub struct GapResult {
    pub name: String,
    /// `Some` exactly when `failure` is `None`.
    pub gapped: Option<String>,
    /// The template allele's `id`.
    pub template: Option<String>,
    pub template_species: Option<String>,
    pub identity: Option<f64>,
    pub cdr1_length: Option<usize>,
    pub cdr2_length: Option<usize>,
    pub warnings: Vec<String>,
    pub failure: Option<String>,
}

struct Template<'a> {
    allele: &'a Allele,
    /// IMGT position of each residue in `aa`.
    pos: Vec<u16>,
    aa: Vec<u8>,
    cdr: Vec<bool>,
}

/// One (frame, template) pair that passed placement and the required checks.
struct Built {
    score: i32,
    template: usize,
    gapped: String,
    cys104: bool,
    identity: f64,
    warnings: Vec<&'static str>,
    cdr1: usize,
    cdr2: usize,
}

/// `Err(Some(reason))` for a named failure, `Err(None)` for a pair that failed the required checks.
type Attempt = Result<Built, Option<&'static str>>;

/// Gap each `(name, sequence)` against the anchored V alleles of `templates`, in input order.
pub fn gap_v(queries: &[(String, String)], templates: &[Allele]) -> Vec<GapResult> {
    let templates: Vec<Template> = templates
        .iter()
        .filter(|a| a.segment == Segment::V && a.anchor.is_some())
        .map(|a| {
            let res = template_residues(&a.sequence_gapped);
            Template {
                allele: a,
                pos: res.iter().map(|r| r.pos).collect(),
                aa: res.iter().map(|r| r.aa).collect(),
                cdr: res.iter().map(|r| is_cdr(region(r.pos))).collect(),
            }
        })
        .collect();
    queries
        .par_iter()
        .map(|(name, seq)| gap_one(name, seq, &templates))
        .collect()
}

fn gap_one(name: &str, seq: &str, all: &[Template]) -> GapResult {
    let mut out = GapResult {
        name: name.to_string(),
        gapped: None,
        template: None,
        template_species: None,
        identity: None,
        cdr1_length: None,
        cdr2_length: None,
        warnings: Vec::new(),
        failure: None,
    };
    let cleaned: Vec<u8> = seq
        .bytes()
        .filter(|b| *b != b'.' && !b.is_ascii_whitespace())
        .map(|b| b.to_ascii_uppercase())
        .collect();
    if cleaned.is_empty() {
        out.failure = Some("empty_sequence".into());
        return out;
    }

    let locus: Option<Locus> = name.get(..3).and_then(|l| l.parse().ok());
    let mut cands: Vec<&Template> = all
        .iter()
        .filter(|t| Some(t.allele.locus) == locus)
        .collect();
    if cands.is_empty() {
        cands = all.iter().collect();
    }
    let occ = occupancy(&cands);

    let frames: Vec<Vec<u8>> = (0..3)
        .map(|f| {
            cleaned
                .get(f..)
                .unwrap_or_default()
                .as_chunks::<3>()
                .0
                .iter()
                .map(|c| translate(c))
                .collect()
        })
        .collect();
    let mut ranked: Vec<(usize, usize, align::ProteinAlignment)> = cands
        .iter()
        .enumerate()
        .flat_map(|(ti, t)| {
            let frames = &frames;
            (0..3).map(move |f| (ti, f, align(&frames[f], &t.aa, &t.cdr, &Scoring::DEFAULT)))
        })
        .collect();
    ranked.sort_by_key(|(ti, f, a)| (Reverse(a.score), *ti, *f));

    let long = |a: &align::ProteinAlignment| {
        a.score > 0 && a.cols.iter().filter(|c| matches!(c, Col::Pair(..))).count() >= MIN_ALIGNED
    };
    if !ranked.first().is_some_and(|(_, _, a)| long(a)) {
        out.failure = Some("no_alignment".into());
        return out;
    }
    // Lower pairs are held to the same floor, so a short wrong-frame hit is never written in
    // place of a refusal.
    ranked.retain(|(_, _, a)| long(a));

    let attempts: Vec<Attempt> = ranked
        .iter()
        .take(TOP_N)
        .map(|(ti, f, a)| try_pair(&cleaned, *f, &frames[*f], *ti, cands[*ti], a, &occ))
        .collect();
    let passing: Vec<&Built> = attempts.iter().filter_map(|a| a.as_ref().ok()).collect();
    let Some(win) = passing.iter().find(|b| b.cys104).or(passing.first()) else {
        let reason = attempts.iter().find_map(|a| *a.as_ref().err()?);
        // Placement passes the required checks by construction, so some pair has a named reason.
        debug_assert!(
            reason.is_some(),
            "{name}: every top pair failed the required checks"
        );
        out.failure = Some(reason.unwrap_or("no_alignment").into());
        return out;
    };

    let mut warnings = Vec::new();
    if !win.cys104 {
        warnings.push("no_cys104");
    }
    warnings.extend(&win.warnings);
    if win.identity < LOW_IDENTITY {
        warnings.push("low_identity");
    }
    if passing
        .iter()
        .any(|b| (b.score - win.score).abs() <= TIE_MARGIN && b.gapped != win.gapped)
    {
        warnings.push("template_disagreement");
    }
    let t = cands[win.template].allele;
    out.gapped = Some(win.gapped.clone());
    out.template = Some(t.id.clone());
    out.template_species = Some(t.species.clone());
    out.identity = Some(win.identity);
    out.cdr1_length = Some(win.cdr1);
    out.cdr2_length = Some(win.cdr2);
    out.warnings = warnings.into_iter().map(String::from).collect();
    out
}

/// Per IMGT position, how many candidate templates have a residue there.
fn occupancy(cands: &[&Template]) -> Vec<u32> {
    let max = cands
        .iter()
        .flat_map(|t| t.pos.iter())
        .max()
        .copied()
        .unwrap_or(0);
    let mut occ = vec![0; max as usize + 1];
    for t in cands {
        for &p in &t.pos {
            occ[p as usize] += 1;
        }
    }
    occ
}

fn is_cdr(r: Region) -> bool {
    matches!(r, Region::Cdr1 | Region::Cdr2)
}

fn try_pair(
    cleaned: &[u8],
    frame: usize,
    q: &[u8],
    ti: usize,
    t: &Template,
    aln: &align::ProteinAlignment,
    occ: &[u32],
) -> Attempt {
    let mut cols = aln.cols.clone();
    slide_deletions(&mut cols, q, t, occ);
    let mut pairs: Vec<(usize, usize)> = cols
        .iter()
        .filter_map(|c| match *c {
            Col::Pair(i, j) => Some((i, j)),
            _ => None,
        })
        .collect();
    // Positions after 104 are contiguous whatever the template's tail does: drop the tail
    // pairs, so the residues after the last one at <= 104 run on from 105 as an overhang.
    let reaches_tail = pairs.last().is_some_and(|&(_, j)| t.pos[j] > 104);
    pairs.retain(|&(_, j)| t.pos[j] <= 104);
    let (Some(&(qf, tf)), Some(&(ql, tl))) = (pairs.first(), pairs.last()) else {
        return Err(Some("unplaceable_overhang"));
    };

    let n = q.len();
    let mut pos: Vec<Option<u16>> = vec![None; n];
    for &(i, j) in &pairs {
        pos[i] = Some(t.pos[j]);
    }

    // Extra query residues inside one framework fill its most-occupied empty positions.
    for w in pairs.windows(2) {
        let ((qa, ta), (qb, tb)) = (w[0], w[1]);
        let r = qb - qa - 1;
        let (pa, pb) = (t.pos[ta], t.pos[tb]);
        if r == 0 || region(pa) != region(pb) || is_cdr(region(pa)) {
            continue;
        }
        let mut free: Vec<u16> = (pa + 1..pb).collect();
        if r > free.len() {
            return Err(Some("framework_insertion"));
        }
        free.sort_by_key(|&p| Reverse(occ[p as usize])); // stable: ties stay leftmost
        free.truncate(r);
        free.sort_unstable();
        for (k, p) in free.into_iter().enumerate() {
            pos[qa + 1 + k] = Some(p);
        }
    }

    // Overhangs.
    let p0 = t.pos[tf] as usize;
    if qf >= p0 {
        return Err(Some("unplaceable_overhang"));
    }
    for (i, slot) in pos.iter_mut().enumerate().take(qf) {
        *slot = Some((p0 - qf + i) as u16);
    }
    if ql + 1 < n {
        if t.pos[tl] < 104 && !reaches_tail {
            return Err(Some("unplaceable_overhang"));
        }
        for (i, slot) in pos.iter_mut().enumerate().skip(ql + 1) {
            *slot = Some((104 + i - ql) as u16);
        }
    }

    // CDRs by the length rule, between the flanking framework residues.
    for (cdr, fr_end, fr_start) in [(Cdr::One, 26u16, 39u16), (Cdr::Two, 55, 66)] {
        let Some(r) = pos.iter().position(|p| p.is_some_and(|p| p >= fr_start)) else {
            continue;
        };
        let l = pos[..r]
            .iter()
            .rposition(|p| p.is_some_and(|p| p <= fr_end));
        let first = l.map_or(0, |l| l + 1);
        let k = r - first;
        let places = match l {
            Some(_) => cdr_positions(cdr, k),
            // 5'-partial start inside the CDR: its full length is unknown, so the template's
            // positions stand; right-align on the last positions only if a residue has none.
            None if pos[first..r].iter().all(Option::is_some) => continue,
            None => (k < (fr_start - fr_end) as usize)
                .then(|| (fr_start - k as u16..fr_start).collect()),
        }
        .ok_or(Some("cdr_too_long"))?;
        for (slot, p) in pos[first..r].iter_mut().zip(places) {
            *slot = Some(p);
        }
    }
    let pos: Vec<u16> = pos
        .into_iter()
        .collect::<Option<_>>()
        .ok_or(Some("framework_insertion"))?;

    // Nucleotides.
    let mut g = vec![b'.'; 3 * pos[n - 1] as usize];
    for (i, &p) in pos.iter().enumerate() {
        let s = 3 * (p as usize - 1);
        g[s..s + 3].copy_from_slice(&cleaned[frame + 3 * i..frame + 3 * i + 3]);
    }
    let start = 3 * (pos[0] as usize - 1);
    if frame > start {
        return Err(Some("unplaceable_overhang"));
    }
    g[start - frame..start].copy_from_slice(&cleaned[..frame]);
    g.extend_from_slice(&cleaned[frame + 3 * n..]);
    let gapped = String::from_utf8(g).map_err(|_| None)?;

    // Required checks.
    let mut runs = Vec::new(); // (start, len) of dot runs
    for (i, b) in gapped.bytes().enumerate() {
        match runs.last_mut() {
            Some((s, l)) if b == b'.' && *s + *l == i => *l += 1,
            _ if b == b'.' => runs.push((i, 1)),
            _ => {}
        }
    }
    let codon_runs = runs
        .iter()
        .all(|&(s, l)| s == 0 || (s % 3 == 0 && l % 3 == 0));
    if !codon_runs || gapped.replace('.', "").as_bytes() != cleaned {
        return Err(None);
    }

    let aa_at = |p: u16| pos.iter().position(|&x| x == p).map(|i| q[i]);
    let mut warnings = Vec::new();
    if q.contains(&b'*') {
        warnings.push("stop_codon");
    }
    if aa_at(23) != Some(b'C') {
        warnings.push("no_cys23");
    }
    if aa_at(41) != Some(b'W') {
        warnings.push("no_trp41");
    }
    let same = pairs.iter().filter(|&&(i, j)| q[i] == t.aa[j]).count();
    let in_cdr = |lo, hi| pos.iter().filter(|&&p| (lo..=hi).contains(&p)).count();
    Ok(Built {
        score: aln.score,
        template: ti,
        cys104: matches!(gapped.get(309..312), Some("TGT" | "TGC")),
        gapped,
        identity: same as f64 / pairs.len() as f64,
        warnings,
        cdr1: in_cdr(27, 38),
        cdr2: in_cdr(56, 65),
    })
}

/// Move each framework `TemplateOnly` run, among the placements with the same alignment score,
/// to the one whose emptied positions the candidate templates fill least (ties leftmost).
fn slide_deletions(cols: &mut [Col], q: &[u8], t: &Template, occ: &[u32]) {
    let is_pair = |c: &Col| matches!(c, Col::Pair(..));
    let mut i = 0;
    while i < cols.len() {
        if !matches!(cols[i], Col::TemplateOnly(_)) {
            i += 1;
            continue;
        }
        let start = i;
        while matches!(cols[i], Col::TemplateOnly(_)) {
            i += 1; // cols end in a Pair, so this stays in bounds
        }
        let end = i;
        if !is_pair(&cols[start - 1]) || !is_pair(&cols[end]) {
            continue;
        }
        let mut lo = start;
        while lo > 0 && is_pair(&cols[lo - 1]) {
            lo -= 1;
        }
        let mut hi = end;
        while hi < cols.len() && is_pair(&cols[hi]) {
            hi += 1;
        }
        let Col::Pair(q0, t0) = cols[lo] else {
            unreachable!()
        };
        let k = end - start;
        let npairs = hi - lo - k;
        let g0 = start - lo;
        // With g pairs before the run: pair u is (q0+u, t0+u) for u < g, (q0+u, t0+u+k) after.
        let score = |g: usize| -> i32 {
            (0..npairs)
                .map(|u| blosum62(q[q0 + u], t.aa[t0 + u + if u < g { 0 } else { k }]))
                .sum()
        };
        let emptied = |g: usize| (t0 + g..t0 + g + k).map(|j| t.pos[j]);
        let reg = region(t.pos[t0 + g0]);
        if is_cdr(reg) || emptied(g0).any(|p| region(p) != reg) {
            continue;
        }
        let s0 = score(g0);
        let ok = |g: usize| score(g) == s0 && emptied(g).all(|p| region(p) == reg);
        let (mut a, mut b) = (g0, g0);
        while a > 1 && ok(a - 1) {
            a -= 1;
        }
        while b + 1 < npairs && ok(b + 1) {
            b += 1;
        }
        let cost = |g: usize| emptied(g).map(|p| occ[p as usize]).sum::<u32>();
        let g = (a..=b).min_by_key(|&g| cost(g)).unwrap(); // first minimum: leftmost
        for u in 0..g {
            cols[lo + u] = Col::Pair(q0 + u, t0 + u);
        }
        for v in 0..k {
            cols[lo + g + v] = Col::TemplateOnly(t0 + g + v);
        }
        for u in g..npairs {
            cols[lo + k + u] = Col::Pair(q0 + u, t0 + u + k);
        }
        i = lo + g + k;
    }
}
