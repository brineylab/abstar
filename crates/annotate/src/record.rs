//! The output `Record` of one read, built from its calls (record 0035).
//!
//! Internal coordinates are 0-based half-open on the oriented read. They become AIRR 1-based
//! closed positions here, once (record 0016): a half-open `(a, b)` is `start = a + 1`,
//! `end = b`. Germline positions convert the same way on the ungapped allele as used (the
//! reverse complement for an inverted D).
use crate::align::{cigar, Alignment, Op, Scoring};
use crate::assign::{CHit, JHit, Junction, Oriented, VHit};
use crate::derive::{alignments, masks, mutations, productivity, region_aa, regions, segments, translate, Col, SegmentMutations};
use crate::schema::Record;
use abstar_germline::Allele;

/// The `invalid_input` failure reason of a sequence, if it has one.
pub fn invalid_input(seq: &str) -> Option<&'static str> {
    if seq.is_empty() {
        Some("empty_sequence")
    } else if !seq.bytes().all(|b| b"ACGTNRYSWKMBDHV".contains(&b.to_ascii_uppercase())) {
        Some("invalid_characters")
    } else {
        None
    }
}

/// A record with only its identity and status.
pub fn failed(id: &str, seq: &str, status: &str, reason: &str) -> Record {
    Record {
        sequence_id: id.into(),
        sequence: seq.into(),
        annotation_status: status.into(),
        failure_reason: Some(reason.into()),
        ..Default::default()
    }
}

/// A V allele whose gapped sequence follows IMGT numbering. One without `.` gaps does not; nor,
/// until #137 regaps them, does a macaque V allele without an anchor (KIMDB's gapping puts no
/// Cys104 at IMGT 310–312).
// ponytail: species-specific rule; drop the macaque clause once #137 regaps the macaque V genes.
pub fn imgt_gapped(a: &Allele) -> bool {
    a.sequence_gapped.contains('.') && !(a.species == "macaque" && a.anchor.is_none())
}

type Call = (Option<String>, Option<String>, Option<String>, Option<Vec<String>>, Option<Vec<String>>);

/// `<seg>_call`, `_gene`, `_species`, `_call_ties`, `_species_ties` (plain names, record 0007).
fn call(alleles: &[Allele], primary: usize, ties: &[usize]) -> Call {
    let a = &alleles[primary];
    (
        Some(a.name.clone()),
        Some(a.gene.clone()),
        Some(a.species.clone()),
        Some(ties.iter().map(|&t| alleles[t].name.clone()).collect()),
        Some(ties.iter().map(|&t| alleles[t].species.clone()).collect()),
    )
}

type Span = (Option<i64>, Option<i64>, Option<i64>, Option<i64>, Option<String>, Option<f64>, Option<f64>);

/// `<seg>_sequence_start/end`, `_germline_start/end`, `_cigar`, `_score`, `_identity` from the
/// segment's reported columns: the same columns its mutations use. The score is under `s`
/// (record 0029's scoring; a non-ACGT read base mismatches) over those columns, and the
/// identity is matches over aligned (non-gap) columns.
fn span(q: &[u8], seg: Option<&(Vec<Col>, Vec<u8>)>, s: &Scoring) -> Span {
    let Some((cols, g)) = seg else { return (None, None, None, None, None, None, None) };
    let mut ops: Vec<(Op, u32)> = Vec::new();
    for &c in cols {
        let op = match c {
            (Some(i), Some(u)) if matches!(q[i], b'A' | b'C' | b'G' | b'T') && q[i] == g[u] => Op::Match,
            (Some(_), Some(_)) => Op::Mismatch,
            (Some(_), None) => Op::Ins,
            _ => Op::Del,
        };
        match ops.last_mut() {
            Some((o, n)) if *o == op => *n += 1,
            _ => ops.push((op, 1)),
        }
    }
    let aligned: Vec<(usize, usize)> = cols.iter().filter_map(|&(i, u)| i.zip(u)).collect();
    let (Some(&(_, g_start)), Some(&(_, g_last))) = (aligned.first(), aligned.last()) else { return (None, None, None, None, None, None, None) };
    let (q_start, q_last) = (cols.iter().find_map(|c| c.0).expect("aligned"), cols.iter().rev().find_map(|c| c.0).expect("aligned"));
    let n = |op: Op| ops.iter().filter(|o| o.0 == op).map(|o| o.1 as i32).sum::<i32>();
    let gaps: i32 = ops.iter().filter(|o| matches!(o.0, Op::Ins | Op::Del)).map(|o| s.gap_open as i32 + (o.1 as i32 - 1) * s.gap_extend as i32).sum();
    let (matches, mismatches) = (n(Op::Match), n(Op::Mismatch));
    let score = matches * s.matched as i32 + mismatches * s.mismatch as i32 + gaps;
    let aln = Alignment { score, q_start, q_end: q_last + 1, g_start, g_end: g_last + 1, ops };
    (
        Some(q_start as i64 + 1),
        Some(q_last as i64 + 1),
        Some(g_start as i64 + 1),
        Some(g_last as i64 + 1),
        Some(cigar(&aln, q.len())),
        Some(score as f64),
        Some(matches as f64 / aligned.len() as f64),
    )
}

type Muts = (
    Option<Vec<String>>,
    Option<Vec<String>>,
    Option<Vec<String>>,
    Option<Vec<String>>,
    Option<i64>,
    Option<i64>,
    Option<f64>,
    Option<f64>,
);

/// `<seg>_mutations`, `_mutations_aa`, `_insertions`, `_deletions`, `_mutation_count`,
/// `_mutation_count_aa`, `_mutation_frequency`, `_mutation_frequency_aa`.
fn muts(m: Option<SegmentMutations>) -> Muts {
    match m {
        None => (None, None, None, None, None, None, None, None),
        Some(m) => (Some(m.nt), m.aa, Some(m.insertions), Some(m.deletions), Some(m.count), m.count_aa, Some(m.frequency), m.frequency_aa),
    }
}

fn text(q: &[u8], (a, b): (usize, usize)) -> String {
    String::from_utf8_lossy(&q[a..b]).into_owned()
}

type Region = (Option<String>, Option<String>, Option<i64>, Option<i64>, Option<i64>);

/// A region's nt, aa, start, end and aa length; start and end are null when it is empty.
fn region(q: &[u8], frame: usize, x: Option<(usize, usize)>) -> Region {
    let Some((a, b)) = x else { return (None, None, None, None, None) };
    let aa = region_aa(q, frame, (a, b));
    let len = aa.len() as i64;
    let (start, end) = if a < b { (Some(a as i64 + 1), Some(b as i64)) } else { (None, None) };
    (Some(text(q, (a, b))), Some(aa), start, end, Some(len))
}

/// The fields every record with a V call carries: identity, strand, locus and the V call.
fn with_v(id: &str, seq: &str, read: &Oriented, v: &VHit, alleles: &[Allele]) -> Record {
    let mut r = Record {
        sequence_id: id.into(),
        sequence: seq.into(),
        annotation_status: "annotated".into(),
        rev_comp: Some(read.rev_comp),
        locus: Some(alleles[v.primary].locus.to_string()),
        ..Default::default()
    };
    (r.v_call, r.v_gene, r.v_species, r.v_call_ties, r.v_species_ties) = call(alleles, v.primary, &v.ties);
    r
}

/// A read whose V allele isn't IMGT-gapped (`imgt_gapped`): the V call and its alignment only,
/// `unassigned` / `macaque_not_imgt_gapped`. Regions, mutations and masks need IMGT numbering.
pub fn v_only(id: &str, seq: &str, read: &Oriented, v: &VHit, alleles: &[Allele], s: &Scoring) -> Record {
    let mut r = with_v(id, seq, read, v, alleles);
    let segs = segments(read, v, None, None, alleles);
    (r.v_sequence_start, r.v_sequence_end, r.v_germline_start, r.v_germline_end, r.v_cigar, r.v_score, r.v_identity) =
        span(&read.seq, segs[0].as_ref(), s);
    r.annotation_status = "unassigned".into();
    r.failure_reason = Some("macaque_not_imgt_gapped".into());
    r
}

/// The record of a read with a V call on an IMGT-gapped allele, and J (with its junction) and C
/// when found. Status: `no_j_hit` without J, `junction_undefined` without a junction (an anchor
/// missing or truncated), else `annotated`. Every field the calls support is filled either way.
#[allow(clippy::too_many_arguments)]
pub fn assemble(
    id: &str,
    seq: &str,
    read: &Oriented,
    v: &VHit,
    j: Option<(&JHit, &Junction)>,
    c: Option<&CHit>,
    alleles: &[Allele],
    s: &Scoring,
) -> Record {
    let q = &read.seq;
    let mut r = with_v(id, seq, read, v, alleles);
    let d = j.and_then(|(_, jn)| jn.d.as_ref());
    if let Some(d) = d {
        (r.d_call, r.d_gene, r.d_species, r.d_call_ties, r.d_species_ties) = call(alleles, d.primary, &d.ties);
        r.d_inverted = Some(d.inverted);
    }
    if let Some((jh, jn)) = j {
        (r.j_call, r.j_gene, r.j_species, r.j_call_ties, r.j_species_ties) = call(alleles, jh.primary, &jh.ties);
        let np1 = text(q, (jn.v_end, d.map_or(jn.j_start, |d| d.start)));
        (r.np1_length, r.np1) = (Some(np1.len() as i64), Some(np1));
        if let Some(d) = d {
            let np2 = text(q, (d.end, jn.j_start));
            (r.np2_length, r.np2) = (Some(np2.len() as i64), Some(np2));
        }
    }
    if let Some(c) = c {
        (r.c_call, r.c_gene, r.c_species, r.c_call_ties, r.c_species_ties) = call(alleles, c.primary, &c.ties);
    }

    let segs = segments(read, v, j, c, alleles);
    (r.v_sequence_start, r.v_sequence_end, r.v_germline_start, r.v_germline_end, r.v_cigar, r.v_score, r.v_identity) = span(q, segs[0].as_ref(), s);
    (r.d_sequence_start, r.d_sequence_end, r.d_germline_start, r.d_germline_end, r.d_cigar, r.d_score, r.d_identity) = span(q, segs[1].as_ref(), s);
    (r.j_sequence_start, r.j_sequence_end, r.j_germline_start, r.j_germline_end, r.j_cigar, r.j_score, r.j_identity) = span(q, segs[2].as_ref(), s);
    (r.c_sequence_start, r.c_sequence_end, r.c_germline_start, r.c_germline_end, _, r.c_score, r.c_identity) = span(q, segs[3].as_ref(), s);
    let [mv, md, mj, mc] = mutations(read, v, j, c, alleles);
    (r.v_mutations, r.v_mutations_aa, r.v_insertions, r.v_deletions, r.v_mutation_count, r.v_mutation_count_aa, r.v_mutation_frequency, r.v_mutation_frequency_aa) = muts(mv);
    (r.d_mutations, r.d_mutations_aa, r.d_insertions, r.d_deletions, r.d_mutation_count, r.d_mutation_count_aa, r.d_mutation_frequency, r.d_mutation_frequency_aa) = muts(md);
    (r.j_mutations, r.j_mutations_aa, r.j_insertions, r.j_deletions, r.j_mutation_count, r.j_mutation_count_aa, r.j_mutation_frequency, r.j_mutation_frequency_aa) = muts(mj);
    (r.c_mutations, r.c_mutations_aa, r.c_insertions, r.c_deletions, r.c_mutation_count, r.c_mutation_count_aa, r.c_mutation_frequency, r.c_mutation_frequency_aa) = muts(mc);

    let g = regions(read, v, j, c, alleles);
    let f = g.frame;
    (r.fwr1, r.fwr1_aa, r.fwr1_start, r.fwr1_end, r.fwr1_length) = region(q, f, g.fwr1);
    (r.cdr1, r.cdr1_aa, r.cdr1_start, r.cdr1_end, r.cdr1_length) = region(q, f, g.cdr1);
    (r.fwr2, r.fwr2_aa, r.fwr2_start, r.fwr2_end, r.fwr2_length) = region(q, f, g.fwr2);
    (r.cdr2, r.cdr2_aa, r.cdr2_start, r.cdr2_end, r.cdr2_length) = region(q, f, g.cdr2);
    (r.fwr3, r.fwr3_aa, r.fwr3_start, r.fwr3_end, r.fwr3_length) = region(q, f, g.fwr3);
    (r.cdr3, r.cdr3_aa, r.cdr3_start, r.cdr3_end, r.cdr3_length) = region(q, f, g.cdr3);
    (r.fwr4, r.fwr4_aa, r.fwr4_start, r.fwr4_end, r.fwr4_length) = region(q, f, g.fwr4);
    if let Some(x) = g.junction {
        let nt = text(q, x);
        (r.junction_aa, r.junction_length, r.junction) = (Some(translate(nt.as_bytes(), 0)), Some(nt.len() as i64), Some(nt));
    }
    r.sequence_aa = Some(translate(q, f % 3));
    (r.sequence_vdj, r.sequence_vdj_aa) = (Some(text(q, g.vdj)), Some(translate(&q[..g.vdj.1], f)));
    if let Some(x) = g.vdjc {
        (r.sequence_vdjc, r.sequence_vdjc_aa) = (Some(text(q, x)), Some(translate(&q[..x.1], f)));
    }
    let a = alignments(read, v, j, alleles, f);
    (r.sequence_alignment, r.germline_alignment, r.sequence_alignment_aa, r.germline_alignment_aa) =
        (Some(a.sequence), Some(a.germline), Some(a.sequence_aa), Some(a.germline_aa));
    let p = productivity(read, &g);
    (r.productive, r.stop_codon, r.vj_in_frame) = (Some(p.productive), Some(p.stop_codon), g.junction.map(|_| p.vj_in_frame));
    let m = masks(read, v, j, alleles, &g);
    (r.cdr_mask, r.cdr_mask_aa, r.gene_segment_mask, r.gene_segment_mask_aa, r.nongermline_mask, r.nongermline_mask_aa) =
        (m.cdr, m.cdr_aa, Some(m.segment), Some(m.segment_aa), Some(m.nongermline), Some(m.nongermline_aa));

    let reason = if j.is_none() {
        Some("no_j_hit")
    } else if g.junction.is_none() {
        Some("junction_undefined")
    } else {
        None
    };
    if let Some(reason) = reason {
        (r.annotation_status, r.failure_reason) = ("unassigned".into(), Some(reason.into()));
    }
    r
}
