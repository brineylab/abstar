/// K-mer prefilter: k-mer -> alleles containing it (CSR, each allele listed once per k-mer).
pub struct KmerIndex {
    k: usize,
    offsets: Vec<u32>,
    postings: Vec<u32>,
    n: usize,
}

fn code(b: u8) -> Option<u32> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Calls `f(kmer)` for every k-mer of `s` that has only ACGT bases.
fn for_kmers(s: &[u8], k: usize, mut f: impl FnMut(usize)) {
    let mask = (1u32 << (2 * k)) - 1;
    let (mut km, mut valid) = (0u32, 0usize);
    for &b in s {
        match code(b) {
            Some(c) => {
                km = ((km << 2) | c) & mask;
                valid += 1;
            }
            None => valid = 0,
        }
        if valid >= k {
            f(km as usize);
        }
    }
}

impl KmerIndex {
    pub fn new(seqs: &[&[u8]], k: usize) -> KmerIndex {
        assert!((1..=12).contains(&k), "k must be 1..=12");
        let mut pairs: Vec<(u32, u32)> = Vec::new();
        for (a, s) in seqs.iter().enumerate() {
            for_kmers(s, k, |km| pairs.push((km as u32, a as u32)));
        }
        pairs.sort_unstable();
        pairs.dedup();
        let mut offsets = vec![0u32; (1 << (2 * k)) + 1];
        for &(km, _) in &pairs {
            offsets[km as usize + 1] += 1;
        }
        for i in 1..offsets.len() {
            offsets[i] += offsets[i - 1];
        }
        KmerIndex { k, offsets, postings: pairs.into_iter().map(|p| p.1).collect(), n: seqs.len() }
    }

    /// (allele, number of query k-mer positions whose k-mer the allele contains), for every
    /// allele sharing at least one k-mer, sorted by allele index.
    pub fn hits(&self, query: &[u8]) -> Vec<(usize, u32)> {
        let mut counts = vec![0u32; self.n];
        for_kmers(query, self.k, |km| {
            for &a in &self.postings[self.offsets[km] as usize..self.offsets[km + 1] as usize] {
                counts[a as usize] += 1;
            }
        });
        counts.into_iter().enumerate().filter(|&(_, c)| c > 0).collect()
    }
}
