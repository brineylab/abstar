//! One synthetic read per call, with every event recorded and the origin of every base.
//!
//! **Coordinates.** All positions are 0-based. `Mutation::pos`, `Indel::read_pos` and
//! `Indel::equiv` address the *untruncated* read (V(D)J after indels, plus the C prefix);
//! `read` is that sequence with `truncate_5` and `truncate_3` bases removed, so a final
//! read position is the untruncated position minus `truncate_5`.
//!
//! **D inversion.** When `d_inverted`, the D used is the reverse complement of the allele:
//! `Origin::D(i)`, `d_trim_5` and `d_trim_3` all refer to that reverse-complemented
//! sequence (its 5' end is the allele's 3' end).

use std::collections::BTreeSet;
use std::ops::Range;

use abstar_germline::Allele;
use rand::Rng;

use crate::shm::pick;
use crate::{load_s5f_dir, motif_model, mutate, Bound, ContextModel, Dist, GenerateError, GenerateParams, Mutation};

const BASES: [u8; 4] = *b"ACGT";

/// Where a read base came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Origin {
    /// IMGT-gapped index in the true V allele.
    V(usize),
    /// Ungapped index in the true D allele, in the orientation used (see module docs).
    D(usize),
    /// Ungapped index in the true J allele.
    J(usize),
    /// N bases between V and D (or V and J).
    N1,
    /// N bases between D and J.
    N2,
    /// Palindromic P-nucleotides at any segment end.
    P,
    /// Ungapped index in the constant allele.
    C(usize),
    /// Bases of an inserted indel.
    Inserted,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndelKind {
    Ins,
    Del,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Indel {
    pub kind: IndelKind,
    /// Untruncated-read position of the first inserted base, or of the first base after
    /// the deleted span.
    pub read_pos: usize,
    /// Ungapped index in the V or J allele of the first deleted base, or of the germline
    /// base right after the insertion.
    pub germline_pos: usize,
    pub len: usize,
    /// Inserted bases, or the deleted germline bases.
    pub bases: String,
    /// Leftmost and rightmost `read_pos` at which the same indel, applied to the
    /// unmutated sequence, gives the same result (VCF-style normalization).
    pub equiv: (usize, usize),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Simulated {
    pub read: Vec<u8>,
    /// `origin[i]` is the origin of `read[i]`.
    pub origin: Vec<Origin>,
    pub v: Allele,
    pub d: Option<Allele>,
    pub d_inverted: bool,
    pub j: Allele,
    pub c: Option<Allele>,
    /// Number of C bases appended (before truncation).
    pub c_len: usize,
    /// Negative = that many P-nucleotides added.
    pub v_trim_3: i64,
    pub d_trim_5: i64,
    pub d_trim_3: i64,
    pub j_trim_5: i64,
    /// N bases as generated (before mutation); P-nucleotides are implied by negative trims.
    /// VJ loci use `np1` only.
    pub np1: String,
    pub np2: String,
    /// Hits in order. Only V(D)J bases that reach the final read are mutated, so every
    /// position lies in `truncate_5 .. truncate_5 + read.len()`.
    pub mutations: Vec<Mutation>,
    pub indels: Vec<Indel>,
    pub truncate_5: usize,
    pub truncate_3: usize,
    /// Realized divergence of the read's V-derived bases from the true V allele.
    pub v_mutation_frequency: f64,
    pub warnings: Vec<String>,
}

/// AIRR productivity of a read, from its origin map.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Productivity {
    /// Read positions from the first base of the Cys104 codon through the last base of the
    /// Trp/Phe118 codon; `None` when either codon is not in the read.
    pub junction: Option<Range<usize>>,
    /// The junction length is a multiple of 3 and the J anchor is in the V frame.
    pub vj_in_frame: bool,
    /// A stop codon in the V-frame translation of the V(D)J bases (C excluded).
    pub stop_codon: bool,
    /// `vj_in_frame`, no stop codon, and the read keeps Cys104 and Trp/Phe118.
    pub productive: bool,
}

pub(crate) fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

/// P-nucleotides for a negative trim: the reverse complement of the end's first/last k bases.
fn p_nt(seg: &[u8], trim: i64, three_prime: bool) -> Vec<u8> {
    let k = (-trim).max(0) as usize;
    revcomp(if three_prime { &seg[seg.len() - k..] } else { &seg[..k] })
}

/// IMGT-gapped index of each ungapped base.
fn gapped_index(a: &Allele) -> Vec<usize> {
    a.sequence_gapped.bytes().enumerate().filter(|&(_, c)| c != b'.').map(|(i, _)| i).collect()
}

fn pick_from<'b, T>(items: &'b [(T, f64)], rng: &mut impl Rng) -> &'b T {
    let w: Vec<f64> = items.iter().map(|x| x.1).collect();
    &items[pick(&w, rng)].0
}

fn draw(d: &Dist, rng: &mut impl Rng) -> i64 {
    d.values[pick(&d.probs, rng)]
}

fn invalid(msg: String) -> GenerateError {
    GenerateError::InvalidArgument(msg)
}

pub fn productivity(read: &[u8], origin: &[Origin], v: &Allele, j: &Allele) -> Productivity {
    let codon = |i: usize| read.get(i..i + 3).map(|c| c.to_ascii_uppercase()).unwrap_or_default();
    // codons start at read positions ≡ i − g (mod 3) for the first V base (IMGT phase)
    let frame = origin.iter().enumerate().find_map(|(i, o)| match *o {
        Origin::V(g) => Some((i as i64 - g as i64).rem_euclid(3) as usize),
        _ => None,
    });
    let vdj_end = origin.iter().position(|o| matches!(o, Origin::C(_))).unwrap_or(read.len());
    let stop_codon = frame.is_some_and(|f| {
        read[f.min(vdj_end)..vdj_end]
            .as_chunks::<3>()
            .0
            .iter()
            .any(|c| matches!(&c.to_ascii_uppercase()[..], b"TAA" | b"TAG" | b"TGA"))
    });
    let start = v.anchor.and_then(|a| {
        let g = gapped_index(v)[a];
        origin.iter().position(|&o| o == Origin::V(g))
    });
    let end = j.anchor.and_then(|b| origin.iter().position(|&o| o == Origin::J(b)));
    let junction = match (start, end) {
        (Some(s), Some(e)) if s + 3 <= e && e + 3 <= read.len() => Some(s..e + 3),
        _ => None,
    };
    let vj_in_frame = junction.as_ref().is_some_and(|r| r.len() % 3 == 0 && frame == Some((r.end - 3) % 3));
    let conserved = junction.as_ref().is_some_and(|r| {
        matches!(&codon(r.start)[..], b"TGT" | b"TGC") && matches!(&codon(r.end - 3)[..], b"TGG" | b"TTT" | b"TTC")
    });
    Productivity { productive: vj_in_frame && !stop_codon && conserved, junction, vj_in_frame, stop_codon }
}

/// A bound model and run parameters, resolved and validated once per run.
pub struct Simulator<'a> {
    bound: &'a Bound,
    params: &'a GenerateParams,
    /// User S5F tables (`params.shm_table`), loaded once.
    table: Option<ContextModel>,
    mutates: bool,
    /// C alleles to draw from when `constant_len` is set.
    c: Vec<(Allele, f64)>,
}

impl<'a> Simulator<'a> {
    /// Fails with `INVALID_ARGUMENT` for out-of-range rates, an unknown `shm_model`, an
    /// unreadable `shm_table`, or a constant region the database cannot supply.
    pub fn new(bound: &'a Bound, params: &'a GenerateParams) -> Result<Self, GenerateError> {
        let p = params;
        let unit = |name: &str, x: f64| {
            if (0.0..=1.0).contains(&x) {
                Ok(())
            } else {
                Err(invalid(format!("{name}: {x} is outside [0, 1]")))
            }
        };
        unit("indel_rate", p.indel_rate)?;
        unit("long_cdr3_weight", p.long_cdr3_weight)?;
        let (lo, hi) = p.mutation.range();
        unit("mutation", lo)?;
        unit("mutation", hi)?;
        for (name, t) in [("truncate_5p", p.truncate_5p), ("truncate_3p", p.truncate_3p)] {
            if t.range().0 < 0.0 {
                return Err(invalid(format!("{name}: truncation cannot be negative")));
            }
        }
        let mutates = match p.shm_model.as_deref().unwrap_or(&bound.model.shm_model) {
            "motif" => true,
            "none" => false,
            x => return Err(invalid(format!("unknown shm_model '{x}'; valid: motif, none"))),
        };
        let table = p.shm_table.as_deref().map(load_s5f_dir).transpose()?;
        if !mutates && table.is_none() && hi > 0.0 {
            return Err(invalid(format!(
                "mutation up to {hi} was requested, but model '{}' does not mutate (shm_model = \"none\"); \
                 set shm_model = \"motif\" or supply shm_table to mutate it",
                bound.model.name
            )));
        }

        let mut c = Vec::new();
        if p.constant_len.is_some() {
            let all = &bound.candidates.c;
            if all.is_empty() {
                return Err(invalid(format!(
                    "constant_len: the germline database has no C alleles for {}",
                    bound.model.locus
                )));
            }
            if p.isotype_weights.is_empty() {
                c = all.iter().map(|a| (a.clone(), 1.0)).collect();
            }
            let genes: BTreeSet<&str> = all.iter().map(|a| a.gene.as_str()).collect();
            for (gene, &w) in &p.isotype_weights {
                let alleles: Vec<&Allele> = all.iter().filter(|a| &a.gene == gene).collect();
                if alleles.is_empty() || !(w.is_finite() && w >= 0.0) {
                    return Err(invalid(format!(
                        "isotype_weights: '{gene}' = {w} is not a non-negative weight for a C gene in the \
                         database; C genes: {}",
                        genes.iter().copied().collect::<Vec<_>>().join(", ")
                    )));
                }
                c.extend(alleles.iter().map(|a| ((*a).clone(), w / alleles.len() as f64)));
            }
            if c.iter().all(|x| x.1 == 0.0) {
                return Err(invalid("isotype_weights: all weights are 0".into()));
            }
        }
        Ok(Self { bound, params, table, mutates, c })
    }

    /// User S5F tables when given, else the motif model if the model mutates.
    fn shm(&self) -> Option<&ContextModel> {
        self.table.as_ref().or(self.mutates.then(motif_model))
    }

    /// N bases for junction `key` ("vd", "dj" or "vj").
    fn n_bases(&self, key: &str, rng: &mut impl Rng) -> Vec<u8> {
        let m = &self.bound.model;
        let w = self.params.long_cdr3_weight;
        let long = m.insertions_long_tail.get(key).filter(|_| w > 0.0 && rng.random_bool(w));
        let len = draw(long.unwrap_or(&m.insertions[key]), rng).max(0) as usize;
        let mut prev = "start";
        (0..len)
            .map(|_| {
                let i = pick(&m.n_markov[prev], rng);
                prev = ["A", "C", "G", "T"][i];
                BASES[i]
            })
            .collect()
    }

    pub fn simulate(&self, rng: &mut impl Rng) -> Simulated {
        let (b, m, p) = (self.bound, &self.bound.model, self.params);
        let trim = |t: &crate::TrimSet, gene: &str, rng: &mut _| draw(t.genes.get(gene).unwrap_or(&t.default), rng);

        // alleles: J first, then D | J
        let v = pick_from(&b.v, rng).clone();
        let j = pick_from(&b.j, rng).clone();
        let d = b.d_given_j.get(&j.gene).map(|ds| pick_from(ds, rng).clone());
        let d_inverted = d.is_some() && m.d_inversion > 0.0 && rng.random_bool(m.d_inversion);

        // trims (clamped to the segment) and P-nucleotides
        let (vs, js) = (v.sequence().into_bytes(), j.sequence().into_bytes());
        let ds = d.as_ref().map(|d| if d_inverted { revcomp(d.sequence().as_bytes()) } else { d.sequence().into_bytes() });
        let clamp = |t: i64, len: usize| t.clamp(-(len as i64), len as i64);
        let v_trim_3 = clamp(trim(&m.trim.v_3, &v.gene, rng), vs.len());
        let (mut d_trim_5, mut d_trim_3) = (0, 0);
        if let (Some(d), Some(ds), Some(t5), Some(t3)) = (&d, &ds, &m.trim.d_5, &m.trim.d_3) {
            d_trim_5 = clamp(trim(t5, &d.gene, rng), ds.len());
            d_trim_3 = clamp(trim(t3, &d.gene, rng), ds.len() - d_trim_5.max(0) as usize);
        }
        let j_trim_5 = clamp(trim(&m.trim.j_5, &j.gene, rng), js.len());

        // assemble V [P N1 P D P N2] or [P N1] P J
        let (mut seq, mut org) = (Vec::new(), Vec::new());
        let mut push = |bases: &[u8], o: &mut dyn FnMut(usize) -> Origin| {
            for (i, &x) in bases.iter().enumerate() {
                seq.push(x);
                org.push(o(i));
            }
        };
        let vg = gapped_index(&v);
        push(&vs[..vs.len() - v_trim_3.max(0) as usize], &mut |u| Origin::V(vg[u]));
        push(&p_nt(&vs, v_trim_3, true), &mut |_| Origin::P);
        let np1 = self.n_bases(if d.is_some() { "vd" } else { "vj" }, rng);
        push(&np1, &mut |_| Origin::N1);
        let mut np2 = Vec::new();
        if let Some(ds) = &ds {
            let (t5, t3) = (d_trim_5.max(0) as usize, d_trim_3.max(0) as usize);
            push(&p_nt(ds, d_trim_5, false), &mut |_| Origin::P);
            push(&ds[t5..ds.len() - t3], &mut |i| Origin::D(t5 + i));
            push(&p_nt(ds, d_trim_3, true), &mut |_| Origin::P);
            np2 = self.n_bases("dj", rng);
            push(&np2, &mut |_| Origin::N2);
        }
        let jt = j_trim_5.max(0) as usize;
        push(&p_nt(&js, j_trim_5, false), &mut |_| Origin::P);
        push(&js[jt..], &mut |u| Origin::J(jt + u));

        // C prefix, never mutated
        let vdj_len = seq.len();
        let c = (p.constant_len.is_some()).then(|| pick_from(&self.c, rng).clone());
        let mut c_len = 0;
        if let (Some(c), Some(n)) = (&c, p.constant_len) {
            let cs = c.sequence().into_bytes();
            c_len = n.min(cs.len());
            seq.extend(&cs[..c_len]);
            org.extend((0..c_len).map(Origin::C));
        }

        // The indel placement and truncation are drawn before mutation, so mutation can
        // target the V bases of the final read; the indel itself is applied after it.
        let placed = (p.indel_rate > 0.0 && rng.random_bool(p.indel_rate))
            .then(|| self.place_indel(&org, &v, &j, rng))
            .flatten();
        let n = match placed {
            Some((IndelKind::Ins, _, l)) => seq.len() + l,
            Some((IndelKind::Del, _, l)) => seq.len() - l,
            None => seq.len(),
        };
        let nt = |s: crate::Spec, rng: &mut _| s.sample(rng).round().max(0.0) as usize;
        let truncate_5 = nt(p.truncate_5p, rng).min(n);
        let truncate_3 = nt(p.truncate_3p, rng).min(n - truncate_5);
        // untruncated-read position of pre-indel position `i`; `None` if deleted
        let after = |i: usize| match placed {
            Some((IndelKind::Ins, r, l)) if i >= r => Some(i + l),
            Some((IndelKind::Del, r, l)) if i >= r + l => Some(i - l),
            Some((IndelKind::Del, r, _)) if i >= r => None,
            _ => Some(i),
        };
        let kept = |i: usize| after(i).is_some_and(|x| (truncate_5..n - truncate_3).contains(&x));

        // mutate the V(D)J bases that reach the final read; `germ` keeps the unmutated sequence
        let germ = seq.clone();
        let mut mutations = Vec::new();
        let mut warnings = Vec::new();
        let target = p.mutation.sample(rng).clamp(0.0, 1.0);
        if let (Some(model), true) = (self.shm(), target > 0.0) {
            let mutable: Vec<bool> = (0..vdj_len).map(kept).collect();
            let v_pos: Vec<usize> = (0..vdj_len).filter(|&i| mutable[i] && matches!(org[i], Origin::V(_))).collect();
            // productive_only: purifying selection against stop codons in the V frame
            // (record 0024). Codons are those of the pre-indel V(D)J, from V's IMGT phase; a
            // codon an indel later changes is not covered, and a stop there resamples the read.
            let stop_frame = p.productive_only.then(|| {
                org.iter().enumerate().find_map(|(i, o)| match *o {
                    Origin::V(g) => Some((i as i64 - g as i64).rem_euclid(3) as usize),
                    _ => None,
                })
            });
            mutations = mutate(&mut seq[..vdj_len], &mutable, &v_pos, target, model, stop_frame.flatten(), rng);
            let want = (target * v_pos.len() as f64).round() as usize;
            let got = v_pos.iter().filter(|&&i| seq[i] != germ[i]).count();
            if got < want {
                warnings.push(format!(
                    "SHM safety cap reached: realized V divergence {:.4} is below the target {target:.4}",
                    got as f64 / v_pos.len() as f64
                ));
            }
        }

        let mut indels = Vec::new();
        if let Some((kind, r, len)) = placed {
            indels.push(self.apply_indel(kind, r, len, &mut seq, &germ, &mut org, rng));
            mutations.retain_mut(|x| after(x.pos).map(|y| x.pos = y).is_some());
        }

        debug_assert_eq!(seq.len(), n);
        let read = seq[truncate_5..n - truncate_3].to_vec();
        let origin = org[truncate_5..n - truncate_3].to_vec();

        let vgap = v.sequence_gapped.as_bytes();
        let (nv, diff) = read.iter().zip(&origin).fold((0, 0), |(nv, diff), (&x, o)| match *o {
            Origin::V(g) => (nv + 1, diff + usize::from(!x.eq_ignore_ascii_case(&vgap[g]))),
            _ => (nv, diff),
        });
        let v_mutation_frequency = if nv == 0 { 0.0 } else { diff as f64 / nv as f64 };

        Simulated {
            read,
            origin,
            v,
            d,
            d_inverted,
            j,
            c,
            c_len,
            v_trim_3,
            d_trim_5,
            d_trim_3,
            j_trim_5,
            np1: String::from_utf8(np1).expect("ACGT"),
            np2: String::from_utf8(np2).expect("ACGT"),
            mutations,
            indels,
            truncate_5,
            truncate_3,
            v_mutation_frequency,
            warnings,
        }
    }

    /// Draws one indel placement `(kind, position, length)` in V upstream of Cys104 or in J
    /// downstream of Trp/Phe118, at least 3 nt from either codon and with germline on both
    /// sides (so never in the junction, where inserted bases would be indistinguishable
    /// from N). `position` is the first deleted base or the base after the insertion.
    /// `None` if nothing fits.
    fn place_indel(&self, org: &[Origin], v: &Allele, j: &Allele, rng: &mut impl Rng) -> Option<(IndelKind, usize, usize)> {
        let (a, b) = (v.anchor?, j.anchor?);
        let lens: Vec<usize> =
            if self.params.frameshift_indels { (1..=36).collect() } else { (3..=36).step_by(3).collect() };
        let w: Vec<f64> = lens.iter().map(|&l| 1.0 / l as f64).collect();
        let len = lens[pick(&w, rng)];
        let kind = if rng.random_bool(0.5) { IndelKind::Del } else { IndelKind::Ins };
        // CDR1 (IMGT aa 27-38) and CDR2 (56-65) with 3-codon flanks, as gapped nt ranges
        let cdr = |g: usize| (69..123).contains(&g) || (156..204).contains(&g);

        // V bases sit at their ungapped index (V is never 5'-trimmed), so for V `end` is
        // also an ungapped index
        let (mut at, mut weights) = (Vec::new(), Vec::new());
        for i in 1..org.len() {
            let end = if kind == IndelKind::Del { i + len } else { i };
            let wt = match (org[i - 1], org[i], org.get(end)) {
                (Origin::V(_), Origin::V(g), Some(Origin::V(_))) if end + 3 <= a => if cdr(g) { 3.0 } else { 1.0 },
                (Origin::J(_), Origin::J(u), Some(Origin::J(_))) if u >= b + 6 => 1.0,
                _ => continue,
            };
            at.push(i);
            weights.push(wt);
        }
        (!at.is_empty()).then(|| (kind, at[pick(&weights, rng)], len))
    }

    /// Applies a placed indel to the (mutated) `seq`; equivalence is computed on `germ`.
    #[allow(clippy::too_many_arguments)]
    fn apply_indel(
        &self,
        kind: IndelKind,
        r: usize,
        len: usize,
        seq: &mut Vec<u8>,
        germ: &[u8],
        org: &mut Vec<Origin>,
        rng: &mut impl Rng,
    ) -> Indel {
        let p = self.params;
        let germline_pos = match org[r] {
            Origin::J(u) => u,
            _ => r,
        };

        let (mut lo, mut hi) = (r, r);
        let bases = match kind {
            IndelKind::Del => {
                while lo > 0 && germ[lo - 1] == germ[lo + len - 1] {
                    lo -= 1;
                }
                while hi + len < germ.len() && germ[hi] == germ[hi + len] {
                    hi += 1;
                }
                seq.drain(r..r + len);
                org.drain(r..r + len);
                germ[r..r + len].to_vec()
            }
            IndelKind::Ins => {
                // duplication copies the neighbouring bases, upstream when there are enough
                let bases: Vec<u8> = if p.indel_duplication && r >= len {
                    seq[r - len..r].to_vec()
                } else if p.indel_duplication && r + len <= seq.len() {
                    seq[r..r + len].to_vec()
                } else {
                    (0..len).map(|_| BASES[rng.random_range(0..4)]).collect()
                };
                let mut x = bases.clone();
                while lo > 0 && germ[lo - 1] == x[len - 1] {
                    x.rotate_right(1);
                    lo -= 1;
                }
                let mut x = bases.clone();
                while hi < germ.len() && germ[hi] == x[0] {
                    x.rotate_left(1);
                    hi += 1;
                }
                seq.splice(r..r, bases.iter().copied());
                org.splice(r..r, std::iter::repeat_n(Origin::Inserted, len));
                bases
            }
        };
        Indel { kind, read_pos: r, germline_pos, len, bases: String::from_utf8(bases).expect("ACGT"), equiv: (lo, hi) }
    }
}

/// Rebuilds the read from the alleles and recorded events alone (never from `origin`):
/// assemble with trims, P and N, append C, apply indels, then mutations, then truncate.
/// An event inconsistent with the sequence it applies to (a mutation's `from`, a
/// deletion's `bases`) makes the result end in `!`, so it cannot equal a real read.
pub fn replay(s: &Simulated) -> Vec<u8> {
    let (v, j) = (s.v.sequence().into_bytes(), s.j.sequence().into_bytes());
    let mut out = v[..v.len() - s.v_trim_3.max(0) as usize].to_vec();
    out.extend(p_nt(&v, s.v_trim_3, true));
    out.extend(s.np1.bytes());
    if let Some(d) = &s.d {
        let d = if s.d_inverted { revcomp(d.sequence().as_bytes()) } else { d.sequence().into_bytes() };
        out.extend(p_nt(&d, s.d_trim_5, false));
        out.extend(&d[s.d_trim_5.max(0) as usize..d.len() - s.d_trim_3.max(0) as usize]);
        out.extend(p_nt(&d, s.d_trim_3, true));
        out.extend(s.np2.bytes());
    }
    out.extend(p_nt(&j, s.j_trim_5, false));
    out.extend(&j[s.j_trim_5.max(0) as usize..]);
    if let Some(c) = &s.c {
        out.extend(&c.sequence().as_bytes()[..s.c_len]);
    }
    let mut ok = true;
    for i in &s.indels {
        let r = i.read_pos;
        match i.kind {
            IndelKind::Ins => drop(out.splice(r..r, i.bases.bytes())),
            IndelKind::Del => ok &= out.drain(r..r + i.len).eq(i.bases.bytes()),
        }
    }
    for m in &s.mutations {
        ok &= out[m.pos] == m.from;
        out[m.pos] = m.to;
    }
    let mut read = out[s.truncate_5..out.len() - s.truncate_3].to_vec();
    if !ok {
        read.push(b'!');
    }
    read
}
