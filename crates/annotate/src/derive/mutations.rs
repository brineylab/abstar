use super::{imgt_of, translate, Col, Layout};
use crate::assign::{reverse_complement, CHit, JHit, Junction, Oriented, VHit};
use abstar_germline::Allele;

/// One segment's substitutions, indels and their counts against its called allele (record 0033).
#[derive(Clone, Debug, PartialEq)]
pub struct SegmentMutations {
    /// `IMGT{pos}:REF>ALT` (V) or `{pos}:REF>ALT`, in germline order.
    pub nt: Vec<String>,
    /// `IMGT{codon}:REF>ALT` (V) or `{codon}:REF>ALT` over the germline codons whose three
    /// bases are aligned; `None` when the segment has no codon frame (a J allele without an
    /// anchor, or C after such a J).
    pub aa: Option<Vec<String>>,
    pub insertions: Vec<String>,
    pub deletions: Vec<String>,
    pub count: i64,
    pub count_aa: Option<i64>,
    /// `count` over the aligned germline positions.
    pub frequency: f64,
    /// `count_aa` over the counted codons; `None` when `aa` is or no codon is counted.
    pub frequency_aa: Option<f64>,
    /// The first and last aligned germline positions, 0-based half-open on the allele as used
    /// (the reverse complement for an inverted D).
    pub germline: (usize, usize),
}

/// Each segment's (V, D, J, C) alignment columns and germline as used: uppercase, the D
/// reverse-complemented when inverted. V and J follow the junction's cuts (`Layout`), D and C
/// are gapless.
pub(super) fn segments(
    read: &Oriented,
    v: &VHit,
    j: Option<(&JHit, &Junction)>,
    c: Option<&CHit>,
    alleles: &[Allele],
) -> [Option<(Vec<Col>, Vec<u8>)>; 4] {
    let Layout { v: vcols, j: jcols, .. } = Layout::new(read, v, j);
    let up = |k: usize| alleles[k].sequence().to_ascii_uppercase().into_bytes();
    let gapless = |r: std::ops::Range<usize>, g0: usize| -> Vec<Col> { r.enumerate().map(|(n, i)| (Some(i), Some(g0 + n))).collect() };
    let d = j.and_then(|(_, jn)| jn.d.as_ref()).map(|d| {
        let a = up(d.primary);
        (gapless(d.start..d.end, d.g_start), if d.inverted { reverse_complement(&a) } else { a })
    });
    [Some((vcols, up(v.primary))), d, j.map(|(jh, _)| (jcols, up(jh.primary))), c.map(|c| (gapless(c.start..c.end, 0), up(c.primary)))]
}

/// Per-segment mutations of a read (V, D, J, C; record 0033), from the same alignment the regions
/// use: V and J as the aligner placed them (leftmost minimum-mismatch indels, records 0026/0029),
/// cut at the junction step's V end and J start; the D as the junction step placed it; C over
/// the covered prefix. A segment with no aligned base is `None` (ruling A6), and so is V when its
/// allele isn't IMGT-gapped (ruling A10).
///
/// Codon frames: V by IMGT codon; D by the V frame (the phase of the read's last V base); J by
/// Trp/Phe118 (`anchor % 3`); C continuing the J frame past the J allele's end.
pub fn mutations(read: &Oriented, v: &VHit, j: Option<(&JHit, &Junction)>, c: Option<&CHit>, alleles: &[Allele]) -> [Option<SegmentMutations>; 4] {
    let q = &read.seq;
    let segs = segments(read, v, j, c, alleles);
    let va = &alleles[v.primary];
    let (vg, imgt) = (va.sequence_gapped.as_bytes(), imgt_of(va));
    // the D frame: the codon phase of the D's first base, counted from the last V base
    let d_frame = j.and_then(|(_, jn)| jn.d.as_ref()).and_then(|d| {
        let (il, ul) = segs[0].as_ref()?.0.iter().rev().find_map(|&(i, u)| i.zip(u))?;
        let phase = (d.start as i64 - il as i64 + imgt[ul] as i64).rem_euclid(3);
        Some((d.g_start as i64 - phase).rem_euclid(3) as usize)
    });
    let ja = j.map(|(jh, _)| &alleles[jh.primary]);
    let j_frame = ja.and_then(|a| a.anchor).map(|b| b % 3);
    let c_frame = ja.and_then(|a| a.anchor.map(|b| (b as i64 - a.sequence().len() as i64).rem_euclid(3) as usize));
    let frames = [None, d_frame, j_frame, c_frame];
    std::array::from_fn(|k| {
        let (cols, g) = segs[k].as_ref()?;
        let codons: Option<Vec<(String, usize)>> = if k == 0 {
            if !vg.contains(&b'.') {
                return None;
            }
            // (IMGT codon, ungapped index of its first base) of each codon without an IMGT gap
            let mut u = 0;
            let mut out = Vec::new();
            for (n, cod) in vg.as_chunks::<3>().0.iter().enumerate() {
                if !cod.contains(&b'.') {
                    out.push((format!("IMGT{}", n + 1), u));
                }
                u += cod.iter().filter(|&&b| b != b'.').count();
            }
            Some(out)
        } else {
            frames[k].map(|s0| (s0..).step_by(3).take_while(|u| u + 3 <= g.len()).map(|u| (((u - s0) / 3 + 1).to_string(), u)).collect())
        };
        let label = |u: usize| if k == 0 { format!("IMGT{}", imgt[u] + 1) } else { (u + 1).to_string() };
        segment(q, cols, g, label, codons)
    })
}

fn segment(q: &[u8], cols: &[Col], g: &[u8], label: impl Fn(usize) -> String, codons: Option<Vec<(String, usize)>>) -> Option<SegmentMutations> {
    // read position aligned to each germline position
    let mut at = vec![None; g.len()];
    for &(i, u) in cols {
        if let (Some(i), Some(u)) = (i, u) {
            at[u] = Some(i);
        }
    }
    let first = at.iter().position(Option::is_some)?;
    let last = at.iter().rposition(Option::is_some)?;
    let n = at.iter().flatten().count();
    let text = |b: &[u8]| String::from_utf8_lossy(b).into_owned();
    let nt: Vec<String> = (0..g.len())
        .filter_map(|u| {
            let i = at[u].filter(|&i| !q[i].eq_ignore_ascii_case(&g[u]))?;
            Some(format!("{}:{}>{}", label(u), g[u] as char, q[i] as char))
        })
        .collect();
    // an indel is a run of inserted or deleted columns with an aligned column on each side
    let kind = |c: &Col| (c.0.is_some(), c.1.is_some());
    let runs: Vec<&[Col]> = cols.chunk_by(|a, b| kind(a) == kind(b)).collect();
    let (mut insertions, mut deletions) = (Vec::new(), Vec::new());
    for w in runs.windows(3) {
        let (before, gap, after) = (w[0], w[1], w[2]);
        if kind(&before[0]) != (true, true) || kind(&after[0]) != (true, true) {
            continue;
        }
        match kind(&gap[0]) {
            (true, false) => {
                let bases: Vec<u8> = gap.iter().map(|c| q[c.0.unwrap()]).collect();
                insertions.push(format!("{}:ins:{}", label(before[before.len() - 1].1.unwrap()), text(&bases)));
            }
            (false, true) => {
                let u = gap[0].1.unwrap();
                deletions.push(format!("{}:del:{}", label(u), text(&g[u..u + gap.len()])));
            }
            _ => {}
        }
    }
    let aa = codons.map(|cs| {
        let (mut out, mut counted) = (Vec::new(), 0usize);
        for (lab, u) in cs {
            let Some(r) = (u..u + 3).map(|x| at[x].map(|i| q[i])).collect::<Option<Vec<u8>>>() else { continue };
            counted += 1;
            let (from, to) = (translate(&g[u..u + 3], 0), translate(&r, 0));
            if from != to {
                out.push(format!("{lab}:{from}>{to}"));
            }
        }
        (out, counted)
    });
    let count = nt.len() as i64;
    Some(SegmentMutations {
        count_aa: aa.as_ref().map(|a| a.0.len() as i64),
        frequency_aa: aa.as_ref().filter(|a| a.1 > 0).map(|a| a.0.len() as f64 / a.1 as f64),
        aa: aa.map(|a| a.0),
        nt,
        insertions,
        deletions,
        count,
        frequency: count as f64 / n as f64,
        germline: (first, last + 1),
    })
}
