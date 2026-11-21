//! Truth for the annotator's IMGT regions, amino-acid fields, masks and per-segment
//! mutations (decision 0026). Derived from the origin map and the true alleles alone,
//! independently of `abstar_annotate`.
//!
//! The reading frame of every translation is the V's: codons start at read positions in
//! phase with the first V base's IMGT-gapped position.

use crate::pipeline::translate;
use crate::recombine::revcomp;
use crate::schema::TruthRecord;
use crate::{IndelKind, Origin, Simulated};

/// Exclusive IMGT-gapped ends of FWR1, CDR1, FWR2, CDR2 and FWR3 (which ends with Cys104).
const V_ENDS: [usize; 5] = [78, 114, 165, 195, 312];
/// Region indices: FWR1, CDR1, FWR2, CDR2, FWR3, CDR3, FWR4.
const FWR3: usize = 4;
const CDR3: usize = 5;
const FWR4: usize = 6;

/// `R`, `R_aa`, `R_start`, `R_end`, `R_length`.
type Region = (Option<String>, Option<String>, Option<i64>, Option<i64>, Option<i64>);
/// `S_mutations`, `S_mutations_aa`, `S_insertions`, `S_deletions`, `S_mutation_count`,
/// `S_mutation_count_aa`, `S_mutation_frequency`, `S_mutation_frequency_aa`.
type Segment = (
    Option<Vec<String>>,
    Option<Vec<String>>,
    Option<Vec<String>>,
    Option<Vec<String>>,
    Option<i64>,
    Option<i64>,
    Option<f64>,
    Option<f64>,
);

fn text(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

/// Fills the region, sequence, mask and per-segment mutation columns of `r`.
pub(crate) fn fill(r: &mut TruthRecord, s: &Simulated) {
    use Origin::*;
    let (read, o) = (&s.read, &s.origin);
    let vg = s.v.sequence_gapped.as_bytes();
    let gapped = vg.contains(&b'.');
    // ungapped → gapped (`gi`) and gapped → ungapped (`ug`) V indices
    let gi: Vec<usize> = (0..vg.len()).filter(|&g| vg[g] != b'.').collect();
    let mut ug = vec![0; vg.len()];
    for (u, &g) in gi.iter().enumerate() {
        ug[g] = u;
    }
    // germline as used: V, D (reverse-complemented when inverted), J, C
    let germ: [Option<Vec<u8>>; 4] = [
        Some(s.v.sequence().into_bytes()),
        s.d.as_ref().map(|d| if s.d_inverted { revcomp(d.sequence().as_bytes()) } else { d.sequence().into_bytes() }),
        Some(s.j.sequence().into_bytes()),
        s.c.as_ref().map(|c| c.sequence().into_bytes()),
    ];
    // (segment, ungapped germline index) of a germline-derived base
    let seg_of = |x: Origin| match x {
        V(g) => Some((0, ug[g])),
        D(u) => Some((1, u)),
        J(u) => Some((2, u)),
        C(u) => Some((3, u)),
        _ => None,
    };
    let gbase = |x: Origin| seg_of(x).map(|(k, u)| germ[k].as_ref().expect("segment allele")[u]);
    let label = |k: usize, u: usize| if k == 0 { format!("IMGT{}", gi[u] + 1) } else { (u + 1).to_string() };
    // The Cys104 codon is the 3 read bases from the V anchor's first base, whatever their
    // origin (trimming may replace its tail with N bases); CDR3 runs from there to the J's
    // Trp/Phe118 base. Its ends are known when both anchors' first bases survive trimming;
    // truncation may still cut it, leaving the covered part. `cys`: read index of the Cys
    // codon's first base (before the read when truncated), placed from the last V base.
    let v_kept = s.v.sequence().len() as i64 - s.v_trim_3.max(0);
    let cys = s.v.anchor.filter(|&a| (a as i64) < v_kept).and_then(|a| {
        let il = o.iter().rposition(|x| matches!(x, V(_)))?;
        let V(gl) = o[il] else { unreachable!() };
        Some(il as i64 + a as i64 - ug[gl] as i64)
    });
    let cdr3_known = cys.is_some() && s.j.anchor.is_some_and(|b| s.j_trim_5 <= b as i64);

    if let (true, Some(i0)) = (gapped, o.iter().position(|x| matches!(x, V(_)))) {
        let V(g0) = o[i0] else { unreachable!() };
        // the annotated span: first V base to the last base that is not a dangling insertion
        let end = o.iter().rposition(|&x| x != Inserted).expect("a V base") + 1;
        let f0 = i0 + (3 - g0 % 3) % 3;
        let codons = |end: usize| (f0..).step_by(3).take_while(move |p| p + 3 <= end);
        let mut reg = vec![None; end];
        for i in i0..end {
            reg[i] = match o[i] {
                V(g) => Some(V_ENDS.iter().position(|&e| g < e).unwrap_or(CDR3)),
                C(_) => None,
                Inserted => reg[i - 1],
                _ if cys.is_some_and(|c| (i as i64) < c + 3) => Some(FWR3),
                J(u) => s.j.anchor.map(|b| if u < b { CDR3 } else { FWR4 }),
                N1 | N2 | P | D(_) => Some(CDR3),
            };
        }
        let regions: Vec<Region> = (0..7)
            .map(|k| {
                let mut idx = (i0..end).filter(|&i| reg[i] == Some(k));
                let (first, last) = (idx.next(), idx.next_back());
                if k == CDR3 && !cdr3_known {
                    (None, None, None, None, None)
                } else if let Some(a) = first {
                    let z = last.unwrap_or(a);
                    let aa: String = codons(end).filter(|p| reg[p + 2] == Some(k)).map(|p| translate(&read[p..p + 3])).collect();
                    let n = aa.len() as i64;
                    (Some(text(&read[a..=z])), Some(aa), Some(a as i64 + 1), Some(z as i64 + 1), Some(n))
                } else if reg.iter().any(|x| x.is_some_and(|x| x < k)) && reg.iter().any(|x| x.is_some_and(|x| x > k)) {
                    // inside the read, every base deleted
                    (Some(String::new()), Some(String::new()), None, None, Some(0))
                } else {
                    (None, None, None, None, None)
                }
            })
            .collect();
        let [f1, c1, f2, c2, f3, c3, f4]: [Region; 7] = regions.try_into().expect("7 regions");
        (r.fwr1, r.fwr1_aa, r.fwr1_start, r.fwr1_end, r.fwr1_length) = f1;
        (r.cdr1, r.cdr1_aa, r.cdr1_start, r.cdr1_end, r.cdr1_length) = c1;
        (r.fwr2, r.fwr2_aa, r.fwr2_start, r.fwr2_end, r.fwr2_length) = f2;
        (r.cdr2, r.cdr2_aa, r.cdr2_start, r.cdr2_end, r.cdr2_length) = c2;
        (r.fwr3, r.fwr3_aa, r.fwr3_start, r.fwr3_end, r.fwr3_length) = f3;
        (r.cdr3, r.cdr3_aa, r.cdr3_start, r.cdr3_end, r.cdr3_length) = c3;
        (r.fwr4, r.fwr4_aa, r.fwr4_start, r.fwr4_end, r.fwr4_length) = f4;

        r.sequence_aa = Some(translate(&read[f0 % 3..]));
        let last = |f: fn(&Origin) -> bool| o.iter().rposition(f);
        // through the J end, or without a J base through the last V/N/D-derived base
        let vdj = last(|x| matches!(x, J(_))).map_or(end, |z| z + 1);
        let vdjc = last(|x| matches!(x, C(_))).map(|z| z + 1);
        r.sequence_vdj = Some(text(&read[i0..vdj]));
        r.sequence_vdj_aa = Some(translate(&read[f0.min(vdj)..vdj]));
        r.sequence_vdjc = vdjc.map(|z| text(&read[i0..z]));
        r.sequence_vdjc_aa = vdjc.map(|z| translate(&read[f0.min(z)..z]));

        {
            // the masks cover the V(D)J region only (decision 0026): they stop before any constant-region base
            let m = o.iter().position(|x| matches!(x, C(_))).map_or(vdj, |c| c.min(vdj));
            let mut seg = Vec::with_capacity(m - i0);
            for &x in &o[i0..m] {
                seg.push(match x {
                    V(_) => b'V',
                    D(_) => b'D',
                    J(_) => b'J',
                    C(_) => unreachable!("masks stop before the constant region"),
                    N1 | N2 | P => b'N',
                    Inserted => seg[seg.len() - 1],
                });
            }
            let ng: Vec<u8> =
                (i0..m).map(|i| if gbase(o[i]).is_some_and(|g| g.eq_ignore_ascii_case(&read[i])) { b'0' } else { b'1' }).collect();
            let cdr: Option<Vec<u8>> = cdr3_known.then(|| {
                (i0..m)
                    .map(|i| match reg[i] {
                        Some(1) => b'1',
                        Some(3) => b'2',
                        Some(CDR3) => b'3',
                        _ => b'0',
                    })
                    .collect()
            });
            let third = |mask: &[u8]| text(&codons(m).map(|p| mask[p + 2 - i0]).collect::<Vec<_>>());
            let ng_aa: Vec<u8> = codons(m)
                .map(|p| {
                    let g: Option<Vec<u8>> = (p..p + 3).map(|i| gbase(o[i])).collect();
                    match g {
                        Some(g) if translate(&g) == translate(&read[p..p + 3]) => b'0',
                        _ => b'1',
                    }
                })
                .collect();
            r.gene_segment_mask_aa = Some(third(&seg));
            r.cdr_mask_aa = cdr.as_deref().map(third);
            r.cdr_mask = cdr.map(|x| text(&x));
            r.gene_segment_mask = Some(text(&seg));
            r.nongermline_mask = Some(text(&ng));
            r.nongermline_mask_aa = Some(text(&ng_aa));
        }
    } else if let Some(b) = s.j.anchor {
        // no V frame: FWR4 alone, phased by the J's Trp/Phe118 codon
        let end = o.iter().rposition(|&x| x != Inserted).map_or(0, |e| e + 1);
        let mut in4 = vec![false; end];
        for i in 0..end {
            in4[i] = match o[i] {
                J(u) => u >= b,
                Inserted => i > 0 && in4[i - 1],
                _ => false,
            };
        }
        if let (Some(a), Some(z)) = (in4.iter().position(|&x| x), in4.iter().rposition(|&x| x)) {
            let J(u) = o[a] else { unreachable!() };
            let p0 = a + (3 - (u - b) % 3) % 3;
            let codons = (p0..).step_by(3).take_while(|p| p + 3 <= end);
            let aa: String = codons.filter(|p| in4[p + 2]).map(|p| translate(&read[p..p + 3])).collect();
            r.fwr4_length = Some(aa.len() as i64);
            (r.fwr4, r.fwr4_aa) = (Some(text(&read[a..=z])), Some(aa));
            (r.fwr4_start, r.fwr4_end) = (Some(a as i64 + 1), Some(z as i64 + 1));
        }
    }

    // per-segment indels at their leftmost equivalent placement (VCF-style, against the
    // unmutated germline: the left end of `equiv`), kept where the base to their left is
    // still of the same segment and in the read. Listed only with a read base on both sides.
    let (mut ins, mut del): ([Vec<String>; 4], [Vec<String>; 4]) = Default::default();
    // (ungapped germline index, length) of each deletion at its left-normalized placement
    let mut del_at: [Vec<(usize, usize)>; 4] = Default::default();
    let t5 = s.truncate_5;
    for x in &s.indels {
        let is_ins = x.kind == IndelKind::Ins;
        let after = x.read_pos + if is_ins { x.len } else { 0 }; // first base after it
        if x.read_pos <= t5 || after >= t5 + read.len() {
            continue;
        }
        let Some((k, _)) = seg_of(o[after - t5]) else { continue };
        let first = t5 + o.iter().position(|&y| seg_of(y).is_some_and(|y| y.0 == k)).expect("segment base");
        let shift = x.read_pos - x.equiv.0.max(first + 1);
        let g = germ[k].as_ref().expect("segment allele");
        if is_ins {
            let mut bases = x.bases.clone().into_bytes();
            bases.rotate_right(shift % x.len);
            ins[k].push(format!("{}:ins:{}", label(k, x.germline_pos - 1 - shift), text(&bases)));
        } else {
            let u = x.germline_pos - shift;
            del_at[k].push((u, x.len));
            del[k].push(format!("{}:del:{}", label(k, u), text(&g[u..u + x.len])));
        }
    }

    // codon frames: V by IMGT codon; D by the V frame at the read's last V base; J by
    // Trp/Phe118; C continuing the J frame across the J end
    let jlen = s.j.sequence().len();
    let last_v = o.iter().rposition(|x| matches!(x, V(_)));
    let d_frame = last_v.zip(o.iter().position(|x| matches!(x, D(_)))).map(|(il, id)| {
        let (V(gl), D(ud)) = (o[il], o[id]) else { unreachable!() };
        let phase = (id as i64 - il as i64 + gl as i64).rem_euclid(3);
        (ud as i64 - phase).rem_euclid(3) as usize
    });
    let c_frame = s.j.anchor.map(|b| (b as i64 - jlen as i64).rem_euclid(3) as usize);
    let frames = [None, d_frame, s.j.anchor.map(|b| b % 3), c_frame];

    let segs: Vec<Segment> = (0..4)
        .map(|k| {
            let none = (None, None, None, None, None, None, None, None);
            let Some(g) = germ[k].as_ref().filter(|_| k > 0 || gapped) else { return none };
            let mut at = vec![None; g.len()];
            for (i, &x) in o.iter().enumerate() {
                if let Some((_, u)) = seg_of(x).filter(|y| y.0 == k) {
                    at[u] = Some(i);
                }
            }
            let n = at.iter().flatten().count();
            if n == 0 {
                return none;
            }
            let muts: Vec<String> = (0..g.len())
                .filter_map(|u| {
                    let i = at[u].filter(|&i| !read[i].eq_ignore_ascii_case(&g[u]))?;
                    Some(format!("{}:{}>{}", label(k, u), g[u] as char, read[i] as char))
                })
                .collect();
            // (label, first ungapped base) of each germline codon
            let codons: Option<Vec<(String, usize)>> = if k == 0 {
                Some(
                    (0..vg.len() / 3)
                        .filter(|c| !vg[3 * c..3 * c + 3].contains(&b'.'))
                        .map(|c| (format!("IMGT{}", c + 1), ug[3 * c]))
                        .collect(),
                )
            } else {
                frames[k].map(|s0| {
                    (s0..).step_by(3).take_while(|u| u + 3 <= g.len()).map(|u| (((u - s0) / 3 + 1).to_string(), u)).collect()
                })
            };
            // which germline positions count as aligned follows the left-normalized deletions (as the
            // deletion strings do): the read bases keep their order and fill the positions not deleted
            let aligned: Vec<bool> = if del_at[k].is_empty() {
                at.iter().map(Option::is_some).collect()
            } else {
                let deleted = |u: usize| del_at[k].iter().any(|&(d, len)| (d..d + len).contains(&u));
                let (first, n) = (at.iter().position(Option::is_some).expect("a base of the segment"), at.iter().flatten().count());
                let mut ok = vec![false; g.len()];
                for u in (first..g.len()).filter(|&u| !deleted(u)).take(n) {
                    ok[u] = true;
                }
                ok
            };
            let aa = codons.map(|cs| {
                let (mut out, mut nc) = (Vec::new(), 0usize);
                for (lab, u) in cs {
                    if !(u..u + 3).all(|x| aligned[x]) {
                        continue;
                    }
                    nc += 1;
                    // the amino acid comes from the generated placement
                    let Some(pos) = (u..u + 3).map(|x| at[x]).collect::<Option<Vec<_>>>() else { continue };
                    let (from, to) = (translate(&g[u..u + 3]), translate(&pos.iter().map(|&i| read[i]).collect::<Vec<_>>()));
                    if from != to {
                        out.push(format!("{lab}:{from}>{to}"));
                    }
                }
                (out, nc)
            });
            let count = muts.len();
            (
                Some(muts),
                aa.as_ref().map(|a| a.0.clone()),
                Some(std::mem::take(&mut ins[k])),
                Some(std::mem::take(&mut del[k])),
                Some(count as i64),
                aa.as_ref().map(|a| a.0.len() as i64),
                Some(count as f64 / n as f64),
                aa.as_ref().filter(|a| a.1 > 0).map(|a| a.0.len() as f64 / a.1 as f64),
            )
        })
        .collect();
    let [v, d, j, c]: [Segment; 4] = segs.try_into().expect("4 segments");
    (r.v_mutations, r.v_mutations_aa, r.v_insertions, r.v_deletions) = (v.0, v.1, v.2, v.3);
    (r.v_mutation_count, r.v_mutation_count_aa, r.v_mutation_frequency, r.v_mutation_frequency_aa) = (v.4, v.5, v.6, v.7);
    (r.d_mutations, r.d_mutations_aa, r.d_insertions, r.d_deletions) = (d.0, d.1, d.2, d.3);
    (r.d_mutation_count, r.d_mutation_count_aa, r.d_mutation_frequency, r.d_mutation_frequency_aa) = (d.4, d.5, d.6, d.7);
    (r.j_mutations, r.j_mutations_aa, r.j_insertions, r.j_deletions) = (j.0, j.1, j.2, j.3);
    (r.j_mutation_count, r.j_mutation_count_aa, r.j_mutation_frequency, r.j_mutation_frequency_aa) = (j.4, j.5, j.6, j.7);
    (r.c_mutations, r.c_mutations_aa, r.c_insertions, r.c_deletions) = (c.0, c.1, c.2, c.3);
    (r.c_mutation_count, r.c_mutation_count_aa, r.c_mutation_frequency, r.c_mutation_frequency_aa) = (c.4, c.5, c.6, c.7);
}
