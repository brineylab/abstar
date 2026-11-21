use super::{imgt_of, translate, Layout};
use crate::assign::{reverse_complement, JHit, Junction, Oriented, VHit};
use abstar_germline::Allele;

/// `sequence_alignment`, `germline_alignment` and their translations (record 0032).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignments {
    pub sequence: String,
    pub germline: String,
    pub sequence_aa: String,
    pub germline_aa: String,
}

fn upper(a: &Allele) -> Vec<u8> {
    a.sequence().to_ascii_uppercase().into_bytes()
}

/// Codons of alignment columns from `start`: all `.` (IMGT gap) is `.`, all `-` is `-`, a mix
/// of gaps and bases is `X`; a partial terminal codon is dropped.
fn translate_columns(row: &[u8], start: usize) -> String {
    row[start.min(row.len())..]
        .as_chunks::<3>()
        .0
        .iter()
        .map(|c| match c {
            b"..." => ".".to_string(),
            b"---" => "-".to_string(),
            _ if c.iter().any(|&b| b == b'.' || b == b'-') => "X".to_string(),
            _ => translate(c, 0),
        })
        .collect()
}

/// The V(D)J region (`Regions::vdj`) aligned to the called alleles: V with the primary's IMGT
/// gaps (`.` in both strings), N bases against `-` in the germline, D (in its called
/// orientation) and J as the junction step placed them; an inserted base is `-` in the
/// germline, a deleted one `-` in the sequence. The translations group the same columns in
/// threes from the first V-frame codon (`frame`, `Regions::frame`).
pub fn alignments(read: &Oriented, v: &VHit, j: Option<(&JHit, &Junction)>, alleles: &[Allele], frame: usize) -> Alignments {
    let q = &read.seq;
    let l = Layout::new(read, v, j);
    // (sequence, germline) columns
    let mut cols: Vec<(u8, u8)> = Vec::new();
    let base = |r: Option<usize>| r.map_or(b'-', |i| q[i]);
    let (vseq, imgt) = (upper(&alleles[v.primary]), imgt_of(&alleles[v.primary]));
    let (mut prev, mut start) = (None, 0);
    for &(i, u) in &l.v {
        if let Some(u) = u {
            let gaps = prev.map_or(0, |p: usize| imgt[u] - p - 1);
            cols.extend(std::iter::repeat_n((b'.', b'.'), gaps));
            prev = Some(imgt[u]);
        }
        if i == Some(frame) {
            start = cols.len();
        }
        cols.push((base(i), u.map_or(b'-', |u| vseq[u])));
    }
    let n = |cols: &mut Vec<(u8, u8)>, r: std::ops::Range<usize>| cols.extend(r.map(|i| (q[i], b'-')));
    let d = j.and_then(|(_, jn)| jn.d.as_ref());
    let j_start = l.j_start.min(l.end);
    n(&mut cols, l.v_end..d.map_or(j_start, |d| d.start));
    if let Some(d) = d {
        let a = upper(&alleles[d.primary]);
        let dseq = if d.inverted { reverse_complement(&a) } else { a };
        cols.extend((d.start..d.end).zip(&dseq[d.g_start..d.g_end]).map(|(i, &b)| (q[i], b)));
        n(&mut cols, d.end..j_start);
    }
    if let Some((jh, _)) = j {
        let jseq = upper(&alleles[jh.primary]);
        cols.extend(l.j.iter().map(|&(i, u)| (base(i), u.map_or(b'-', |u| jseq[u]))));
    }
    let (s, g): (Vec<u8>, Vec<u8>) = cols.into_iter().unzip();
    let text = |x: &[u8]| String::from_utf8_lossy(x).into_owned();
    Alignments { sequence_aa: translate_columns(&s, start), germline_aa: translate_columns(&g, start), sequence: text(&s), germline: text(&g) }
}
