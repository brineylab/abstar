//! Amino-acid (protein) input: germline translation, V/J assignment, derived fields and the
//! `AaEngine` that assembles them into records (spec 2026-11-22-amino-acid-input).
pub mod assign;
pub mod derive;
pub mod germline;

use self::assign::assign_aa;
use self::derive::{aa_alignments, aa_masks, aa_mutations, aa_regions};
use self::germline::AaSegment;
use crate::align::{blosum62, ProteinScoring};
use crate::derive::{columns, Col};
use crate::engine::{Annotate, EngineOptions};
use crate::error::RunError;
use crate::record::{call, failed};
use crate::schema::Record;
use abstar_germline::{Allele, Database, Segment};

/// The `invalid_input` failure reason of an uppercased protein sequence, if it has one.
fn invalid_aa(seq: &str) -> Option<&'static str> {
    if seq.is_empty() {
        Some("empty_sequence")
    } else if !seq.bytes().all(|b| b"ACDEFGHIKLMNPQRSTVWYXBZJUO*".contains(&b)) {
        Some("invalid_characters")
    } else {
        None
    }
}

/// `_score` (BLOSUM62, affine gaps) and `_identity` (identical over aligned residues) of a
/// segment's reported columns; `None` without an aligned residue.
fn score(q: &[u8], g: &[u8], cols: &[Col], s: &ProteinScoring) -> (Option<f64>, Option<f64>) {
    let (mut score, mut same, mut aligned) = (0, 0, 0);
    let mut prev = None; // the previous column's gap kind: Some(true) inserted, Some(false) deleted
    for &(i, u) in cols {
        let gap = match (i, u) {
            (Some(i), Some(u)) => {
                score += blosum62(q[i], g[u]);
                same += (q[i] == g[u]) as i32;
                aligned += 1;
                None
            }
            (Some(_), None) => Some(true),
            _ => Some(false),
        };
        if gap.is_some() {
            score += if gap == prev { s.gap_extend } else { s.gap_open };
        }
        prev = gap;
    }
    if aligned == 0 {
        return (None, None);
    }
    (Some(score as f64), Some(same as f64 / aligned as f64))
}

/// The database's alleles and their V and J translations, built once per run.
pub struct AaEngine<'a> {
    alleles: &'a [Allele],
    v: AaSegment,
    j: AaSegment,
    scoring: ProteinScoring,
}

impl<'a> AaEngine<'a> {
    pub fn new(db: &'a Database, opts: &EngineOptions) -> Result<AaEngine<'a>, RunError> {
        let seg = |s| AaSegment::new(&db.alleles, s, opts.loci.as_deref(), opts.include_nonfunctional);
        Ok(AaEngine { alleles: &db.alleles, v: seg(Segment::V), j: seg(Segment::J), scoring: ProteinScoring::default() })
    }
}

impl Annotate for AaEngine<'_> {
    /// V and J, then every field of spec §Fields; the rest stay null. Status as `record::assemble`.
    fn annotate(&self, id: &str, seq: &str) -> Record {
        let up = seq.to_ascii_uppercase();
        if let Some(reason) = invalid_aa(&up) {
            return self.failed(id, &up, "invalid_input", reason);
        }
        let (q, a, v, j, s) = (up.as_bytes(), self.alleles, &self.v, &self.j, &self.scoring);
        let Some(c) = assign_aa(q, v, j, a, s) else {
            return self.failed(id, &up, "unassigned", "no_v_hit");
        };
        let mut r = Record {
            sequence_id: id.into(),
            sequence_aa: Some(up.clone()),
            annotation_status: "annotated".into(),
            locus: Some(a[c.v.primary].locus.to_string()),
            ..Default::default()
        };
        (r.v_call, r.v_gene, r.v_species, r.v_call_ties, r.v_species_ties) = call(a, c.v.primary, &c.v.ties);
        (r.v_score, r.v_identity) = score(q, &v.groups[c.v.group].seq, &columns(&c.v.aln, c.v.aln.q_start..c.v_end), s);
        if let Some(jh) = &c.j {
            (r.j_call, r.j_gene, r.j_species, r.j_call_ties, r.j_species_ties) = call(a, jh.primary, &jh.ties);
            (r.j_score, r.j_identity) = score(q, &j.groups[jh.group].seq, &columns(&jh.aln, c.j_start..jh.aln.q_end), s);
        }
        let [mv, mj] = aa_mutations(q, &c, v, j);
        if let Some(m) = mv {
            (r.v_mutations_aa, r.v_insertions_aa, r.v_deletions_aa, r.v_mutation_count_aa, r.v_mutation_frequency_aa) =
                (Some(m.nt), Some(m.insertions), Some(m.deletions), Some(m.count), Some(m.frequency));
        }
        if let Some(m) = mj {
            (r.j_mutations_aa, r.j_insertions_aa, r.j_deletions_aa, r.j_mutation_count_aa, r.j_mutation_frequency_aa) =
                (Some(m.nt), Some(m.insertions), Some(m.deletions), Some(m.count), Some(m.frequency));
        }

        let g = aa_regions(q, &c, v, j);
        let text = |(x, y): (usize, usize)| up[x..y].to_string();
        let region = |x: Option<(usize, usize)>| x.map(|x| (text(x), (x.1 - x.0) as i64)).unzip();
        (r.fwr1_aa, r.fwr1_length) = region(g.fwr1);
        (r.cdr1_aa, r.cdr1_length) = region(g.cdr1);
        (r.fwr2_aa, r.fwr2_length) = region(g.fwr2);
        (r.cdr2_aa, r.cdr2_length) = region(g.cdr2);
        (r.fwr3_aa, r.fwr3_length) = region(g.fwr3);
        (r.cdr3_aa, r.cdr3_length) = region(g.cdr3);
        (r.fwr4_aa, r.fwr4_length) = region(g.fwr4);
        r.junction_aa = g.junction.map(text);
        let (sa, ga) = aa_alignments(q, &c, v, j);
        (r.sequence_alignment_aa, r.germline_alignment_aa) = (Some(sa), Some(ga));
        let m = aa_masks(q, &c, v, j, &g);
        (r.cdr_mask_aa, r.gene_segment_mask_aa, r.nongermline_mask_aa) = (m.cdr, Some(m.segment), Some(m.nongermline));
        let vdj = text(g.vdj);
        let stop = vdj.contains('*');
        let junction_ok = |x: &str| x.starts_with('C') && (x.ends_with('W') || x.ends_with('F')) && !x.contains('X');
        (r.productive, r.stop_codon) = (Some(!stop && r.junction_aa.as_deref().is_some_and(junction_ok)), Some(stop));
        r.sequence_vdj_aa = Some(vdj);

        let reason = if c.j.is_none() {
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

    /// The input, uppercased, goes in `sequence_aa` (also for the pipeline's `invalid_encoding`
    /// and panic records); `sequence` stays null.
    fn failed(&self, id: &str, seq: &str, status: &str, reason: &str) -> Record {
        Record { sequence: None, sequence_aa: Some(seq.to_ascii_uppercase()), ..failed(id, seq, status, reason) }
    }
}
