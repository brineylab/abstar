//! Run-level annotation engine: built once per run, shared by every worker thread.
use crate::align::Scoring;
use crate::assign::{assign_c, assign_j, call_v, resolve_junction, SegmentIndex};
use crate::error::RunError;
use crate::record::{assemble, failed, imgt_gapped, invalid_input, v_only};
use crate::schema::Record;
use crate::stub::annotate_stub;
use abstar_germline::{Allele, Database, Locus, Segment};

/// Annotates one sequence. Implementations are shared across threads.
pub trait Annotate: Sync {
    fn annotate(&self, id: &str, seq: &str) -> Record;
}

#[derive(Clone, Debug, Default)]
pub struct EngineOptions {
    /// Search only these loci (all of the database's when `None`).
    pub loci: Option<Vec<Locus>>,
}

/// The database's alleles and their V, D, J and C indexes, built once per run.
pub struct Engine<'a> {
    alleles: &'a [Allele],
    v: SegmentIndex,
    d: SegmentIndex,
    j: SegmentIndex,
    c: SegmentIndex,
    scoring: Scoring,
}

impl<'a> Engine<'a> {
    pub fn new(db: &'a Database, opts: &EngineOptions) -> Result<Engine<'a>, RunError> {
        let index = |s| SegmentIndex::for_loci(&db.alleles, s, opts.loci.as_deref());
        Ok(Engine {
            alleles: &db.alleles,
            v: index(Segment::V),
            d: index(Segment::D),
            j: index(Segment::J),
            c: index(Segment::C),
            scoring: Scoring::default(),
        })
    }
}

impl Annotate for Engine<'_> {
    /// Orient and V, J (V's locus), the junction with D, C, then the record (record 0035).
    fn annotate(&self, id: &str, seq: &str) -> Record {
        if let Some(reason) = invalid_input(seq) {
            return failed(id, seq, "invalid_input", reason);
        }
        let (a, s) = (self.alleles, &self.scoring);
        let Some((read, v)) = call_v(seq.as_bytes(), a, &self.v, &self.j, s) else {
            return failed(id, seq, "unassigned", "no_v_hit");
        };
        if !imgt_gapped(&a[v.primary]) {
            return v_only(id, seq, &read, &v, a);
        }
        let j = assign_j(&read, &v, a, &self.j, s).map(|j| {
            let has_d = self.d.loci.contains(&a[v.primary].locus);
            let jn = resolve_junction(&read, &v, &j, a, &self.d, has_d);
            (j, jn)
        });
        let c = j.as_ref().and_then(|(j, _)| assign_c(&read, j, a, &self.c));
        assemble(id, seq, &read, &v, j.as_ref().map(|(j, jn)| (j, jn)), c.as_ref(), a)
    }
}

/// The placeholder annotator as an `Annotate` (seam tests only).
pub struct Stub<'a>(pub &'a Database);

impl Annotate for Stub<'_> {
    fn annotate(&self, id: &str, seq: &str) -> Record {
        annotate_stub(id, seq, self.0)
    }
}
