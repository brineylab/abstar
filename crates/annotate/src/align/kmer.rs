/// K-mer prefilter: k-mer -> alleles containing it (CSR, each allele listed once per k-mer).
pub struct KmerIndex {
    k: usize,
    offsets: Vec<u32>,
    postings: Vec<u32>,
    n: usize,
    /// Every occurrence (allele, position) per k-mer, CSR, for `diagonal_hits`.
    occ_offsets: Vec<u32>,
    occ: Vec<(u32, u32)>,
    max_len: usize,
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

/// Calls `f(kmer, start)` for every k-mer of `s` that has only ACGT bases.
fn for_kmers(s: &[u8], k: usize, mut f: impl FnMut(usize, usize)) {
    let mask = (1u32 << (2 * k)) - 1;
    let (mut km, mut valid) = (0u32, 0usize);
    for (i, &b) in s.iter().enumerate() {
        match code(b) {
            Some(c) => {
                km = ((km << 2) | c) & mask;
                valid += 1;
            }
            None => valid = 0,
        }
        if valid >= k {
            f(km as usize, i + 1 - k);
        }
    }
}

impl KmerIndex {
    pub fn new(seqs: &[&[u8]], k: usize) -> KmerIndex {
        assert!((1..=12).contains(&k), "k must be 1..=12");
        let mut pairs: Vec<(u32, u32)> = Vec::new();
        let mut occs: Vec<(u32, u32, u32)> = Vec::new();
        for (a, s) in seqs.iter().enumerate() {
            for_kmers(s, k, |km, pos| {
                pairs.push((km as u32, a as u32));
                occs.push((km as u32, a as u32, pos as u32));
            });
        }
        occs.sort_unstable();
        let mut occ_offsets = vec![0u32; (1 << (2 * k)) + 1];
        for &(km, _, _) in &occs {
            occ_offsets[km as usize + 1] += 1;
        }
        for i in 1..occ_offsets.len() {
            occ_offsets[i] += occ_offsets[i - 1];
        }
        let occ = occs.into_iter().map(|(_, a, p)| (a, p)).collect();
        let max_len = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
        pairs.sort_unstable();
        pairs.dedup();
        let mut offsets = vec![0u32; (1 << (2 * k)) + 1];
        for &(km, _) in &pairs {
            offsets[km as usize + 1] += 1;
        }
        for i in 1..offsets.len() {
            offsets[i] += offsets[i - 1];
        }
        KmerIndex { k, offsets, postings: pairs.into_iter().map(|p| p.1).collect(), n: seqs.len(), occ_offsets, occ, max_len }
    }

    /// (allele, the most query k-mer positions shared on one diagonal band), for every allele
    /// sharing at least one k-mer, sorted by allele index. A shared k-mer at query position i
    /// and allele position j lies on diagonal i - j; the band is two adjacent bins of `bin`
    /// diagonals, so it tolerates small indels. Unlike `hits`, chance matches scattered over
    /// many diagonals don't add up.
    pub fn diagonal_hits(&self, query: &[u8], bin: usize) -> Vec<(usize, u32)> {
        let bin = bin.max(1);
        let nbins = (query.len() + self.max_len) / bin + 2;
        let mut counts = vec![0u16; self.n * nbins];
        for_kmers(query, self.k, |km, qi| {
            for &(a, p) in &self.occ[self.occ_offsets[km] as usize..self.occ_offsets[km + 1] as usize] {
                let c = &mut counts[a as usize * nbins + (qi + self.max_len - p as usize) / bin];
                *c = c.saturating_add(1);
            }
        });
        (0..self.n)
            .filter_map(|a| {
                let row = &counts[a * nbins..(a + 1) * nbins];
                let best = row.windows(2).map(|w| w[0] as u32 + w[1] as u32).max().unwrap_or(0);
                (best > 0).then_some((a, best))
            })
            .collect()
    }

    /// (allele, number of query k-mer positions whose k-mer the allele contains), for every
    /// allele sharing at least one k-mer, sorted by allele index.
    pub fn hits(&self, query: &[u8]) -> Vec<(usize, u32)> {
        let mut counts = vec![0u32; self.n];
        for_kmers(query, self.k, |km, _| {
            for &a in &self.postings[self.offsets[km] as usize..self.offsets[km + 1] as usize] {
                counts[a as usize] += 1;
            }
        });
        counts.into_iter().enumerate().filter(|&(_, c)| c > 0).collect()
    }
}
