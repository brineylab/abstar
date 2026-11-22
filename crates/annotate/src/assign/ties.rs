use abstar_germline::{Allele, Functionality};
use std::cmp::Ordering;

fn functionality_rank(f: Functionality) -> u8 {
    match f {
        Functionality::F => 0,
        Functionality::ORF => 1,
        Functionality::P => 2,
        Functionality::NF => 3,
        Functionality::Unknown => 4,
    }
}

/// Natural order: digit runs compare as numbers, everything else bytewise; names equal under
/// that order fall back to plain byte order.
fn natural_cmp(a: &str, b: &str) -> Ordering {
    let (mut x, mut y) = (a.as_bytes(), b.as_bytes());
    let digits = |s: &[u8]| s.iter().take_while(|c| c.is_ascii_digit()).count();
    let number = |s: &[u8]| {
        let z = s.iter().take_while(|&&c| c == b'0').count();
        s[z..].to_vec()
    };
    loop {
        match (x.first(), y.first()) {
            (None, None) => return a.cmp(b),
            (None, _) => return Ordering::Less,
            (_, None) => return Ordering::Greater,
            (Some(p), Some(q)) if p.is_ascii_digit() && q.is_ascii_digit() => {
                let (n, m) = (digits(x), digits(y));
                let (u, v) = (number(&x[..n]), number(&y[..m]));
                let o = u.len().cmp(&v.len()).then_with(|| u.cmp(&v));
                if o != Ordering::Equal {
                    return o;
                }
                (x, y) = (&x[n..], &y[m..]);
            }
            (Some(p), Some(q)) => {
                if p != q {
                    return p.cmp(q);
                }
                (x, y) = (&x[1..], &y[1..]);
            }
        }
    }
}

/// Record 0008 order among alleles already tied on score: F > ORF > P > NF > unknown, then
/// higher usage (a missing usage after any given one), then natural sort of names.
pub(crate) fn tiebreak_cmp(a: &Allele, b: &Allele) -> Ordering {
    functionality_rank(a.functionality)
        .cmp(&functionality_rank(b.functionality))
        .then_with(|| match (a.usage, b.usage) {
            (Some(x), Some(y)) => y.total_cmp(&x),
            (Some(_), None) => Ordering::Less,
            (None, Some(_)) => Ordering::Greater,
            (None, None) => Ordering::Equal,
        })
        .then_with(|| natural_cmp(&a.name, &b.name))
}

/// The primary call among `tied` (indices into `alleles`, all tied on fewest mismatches), by
/// record 0008; equal alleles fall back to the lower index. Panics on an empty `tied`.
pub fn primary_by_tiebreak(tied: &[usize], alleles: &[Allele]) -> usize {
    *tied
        .iter()
        .min_by(|&&x, &&y| tiebreak_cmp(&alleles[x], &alleles[y]).then(x.cmp(&y)))
        .expect("primary_by_tiebreak: empty tie set")
}
