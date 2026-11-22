use crate::{Allele, Segment};

/// V Cys104 anchor: IMGT positions 310-312 (gapped `309..312`) must be `TGT`/`TGC`;
/// returns the ungapped offset of that codon's first base. A sequence without any `.`
/// is not IMGT-gapped, so its offsets are not IMGT positions: never an anchor.
pub fn find_v_anchor(gapped: &str) -> Option<usize> {
    if !gapped.contains('.') {
        return None;
    }
    let codon = gapped.get(309..312)?.to_ascii_uppercase();
    matches!(codon.as_str(), "TGT" | "TGC")
        .then(|| gapped.as_bytes()[..309].iter().filter(|&&b| b != b'.').count())
}

/// J Trp118/Phe118 anchor: first offset (any frame) matching
/// `(TGG|TTT|TTC)GG[ACGT]{4}GG[ACGT]`.
pub fn find_j_anchor(ungapped: &str) -> Option<usize> {
    let nt = |b: &u8| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T');
    let is = |w: &[u8], s: &[u8]| w.eq_ignore_ascii_case(s);
    ungapped.as_bytes().windows(12).position(|w| {
        (is(&w[..3], b"TGG") || is(&w[..3], b"TTT") || is(&w[..3], b"TTC"))
            && is(&w[3..5], b"GG")
            && w[5..9].iter().all(nt)
            && is(&w[9..11], b"GG")
            && nt(&w[11])
    })
}

/// Sets `a.anchor` by the build rules and returns the allele's anchor warning, if any.
/// On entry `a.anchor` holds the OGRDB J anchor (`j_cdr3_end - 1`) when OGRDB supplied one.
pub(crate) fn assign(a: &mut Allele) -> Option<String> {
    match a.segment {
        Segment::V => {
            a.anchor = find_v_anchor(&a.sequence_gapped);
            a.anchor.is_none().then(|| {
                if a.sequence_gapped.contains('.') {
                    "no Cys104 codon at IMGT positions 310-312".into()
                } else {
                    "V sequence is not IMGT-gapped: no Cys104 anchor assigned".into()
                }
            })
        }
        Segment::J => {
            let motif = find_j_anchor(&a.sequence());
            match a.anchor {
                Some(ogrdb) if motif != Some(ogrdb) => Some(format!(
                    "motif anchor {} disagrees with OGRDB j_cdr3_end {}",
                    motif.map_or("none".to_string(), |m| m.to_string()),
                    ogrdb + 1
                )),
                Some(_) => None,
                None => {
                    a.anchor = motif;
                    motif.is_none().then(|| "no J Trp/Phe118 motif (TGG|TTT|TTC)GGNNNNGGN".into())
                }
            }
        }
        Segment::D | Segment::C => None,
    }
}
