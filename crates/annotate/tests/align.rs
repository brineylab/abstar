use abstar_annotate::align::*;

const S: Scoring = Scoring { matched: 2, mismatch: -3, gap_open: -10, gap_extend: -1 };

/// Deterministic pseudo-random ACGT (xorshift64).
fn rand_seq(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x >> 33) as usize % 4]
        })
        .collect()
}

fn cat(parts: &[&[u8]]) -> Vec<u8> {
    parts.concat()
}

fn ops_score(ops: &[(Op, u32)]) -> i32 {
    ops.iter()
        .map(|&(op, n)| match op {
            Op::Match => n as i32 * S.matched as i32,
            Op::Mismatch => n as i32 * S.mismatch as i32,
            _ => S.gap_open as i32 + (n as i32 - 1) * S.gap_extend as i32,
        })
        .sum()
}

/// An alignment scored from its ops.
fn mk(q_start: usize, q_end: usize, g_start: usize, g_end: usize, ops: Vec<(Op, u32)>) -> Alignment {
    Alignment { score: ops_score(&ops), q_start, q_end, g_start, g_end, ops }
}

/// Ops and coordinates agree with each other and with the sequences; the score is the ops' score
/// plus any dropped terminal gaps' cost (so never above the ops' score).
fn check(a: &Alignment, q: &[u8], g: &[u8]) {
    let (mut qi, mut gi) = (a.q_start, a.g_start);
    for &(op, n) in &a.ops {
        let n = n as usize;
        assert!(n > 0);
        match op {
            Op::Match | Op::Mismatch => {
                for k in 0..n {
                    let m = q[qi + k] == g[gi + k] && b"ACGT".contains(&q[qi + k]);
                    assert_eq!(m, op == Op::Match, "op {op:?} at q{} g{}", qi + k, gi + k);
                }
                qi += n;
                gi += n;
            }
            Op::Ins => qi += n,
            Op::Del => gi += n,
        }
    }
    assert_eq!((qi, gi), (a.q_end, a.g_end));
    assert!(a.score <= ops_score(&a.ops));
    for w in a.ops.windows(2) {
        assert_ne!(w[0].0, w[1].0, "adjacent ops of one kind are merged");
    }
    if let (Some(f), Some(l)) = (a.ops.first(), a.ops.last()) {
        assert!(matches!(f.0, Op::Match | Op::Mismatch) && matches!(l.0, Op::Match | Op::Mismatch), "no terminal gaps");
    }
}

fn gaps(a: &Alignment) -> Vec<(Op, u32)> {
    a.ops.iter().copied().filter(|(op, _)| matches!(op, Op::Ins | Op::Del)).collect()
}

/// (query, germline) position where the first gap op starts.
fn gap_at(a: &Alignment) -> (usize, usize) {
    let (mut qi, mut gi) = (a.q_start, a.g_start);
    for &(op, n) in &a.ops {
        if matches!(op, Op::Ins | Op::Del) {
            return (qi, gi);
        }
        qi += n as usize;
        gi += n as usize;
    }
    panic!("no gap");
}

#[test]
fn kmer_hits_rank_true_allele() {
    // 10 unrelated alleles plus, for each, two siblings differing by 3 SNPs.
    let mut alleles = Vec::new();
    for i in 0..10u64 {
        let base = rand_seq(300, i + 1);
        for s in 0..3usize {
            let mut a = base.clone();
            for p in [60, 150, 240] {
                if s > 0 {
                    a[p + s] = if a[p + s] == b'A' { b'C' } else { b'A' };
                }
            }
            alleles.push(a);
        }
    }
    let refs: Vec<&[u8]> = alleles.iter().map(|a| a.as_slice()).collect();
    let idx = KmerIndex::new(&refs, 7);
    for x in [0, 4, 17, 29] {
        let read = &alleles[x][20..280];
        let hits = idx.hits(read);
        assert!(hits.windows(2).all(|w| w[0].0 < w[1].0), "sorted by allele index");
        let best = hits.iter().max_by_key(|&&(a, c)| (c, std::cmp::Reverse(a))).unwrap();
        assert_eq!(best.0, x);
        assert_eq!(best.1, (read.len() - 6) as u32);
    }
    // Non-ACGT bases break k-mers; nothing shared means no hits.
    assert!(idx.hits(b"NNNNNNNNNNNN").is_empty());
}

#[test]
fn align_exact() {
    let g = rand_seq(300, 7);
    let q = &g[50..250];
    let a = align_semiglobal(q, &g, &S);
    check(&a, q, &g);
    assert_eq!((a.q_start, a.q_end, a.g_start, a.g_end), (0, 200, 50, 250));
    assert_eq!(a.ops, vec![(Op::Match, 200)]);
    assert_eq!(a.score, 400);
}

#[test]
fn align_mismatches() {
    let g = rand_seq(300, 8);
    let mut q = g.clone();
    for p in [40, 41, 150, 260] {
        q[p] = if q[p] == b'G' { b'T' } else { b'G' };
    }
    let a = align_semiglobal(&q, &g, &S);
    check(&a, &q, &g);
    assert_eq!((a.q_start, a.q_end, a.g_start, a.g_end), (0, 300, 0, 300));
    assert_eq!(
        a.ops,
        vec![
            (Op::Match, 40),
            (Op::Mismatch, 2),
            (Op::Match, 108),
            (Op::Mismatch, 1),
            (Op::Match, 109),
            (Op::Mismatch, 1),
            (Op::Match, 39)
        ]
    );
}

#[test]
fn align_deletion_36nt() {
    let v = rand_seq(300, 9);
    let (lead, tail) = (rand_seq(20, 90), rand_seq(30, 91));
    let read = cat(&[&lead, &v[..150], &v[186..], &tail]);
    let a = align_semiglobal(&read, &v, &S);
    check(&a, &read, &v);
    assert_eq!(gaps(&a), vec![(Op::Del, 36)]);
    // leftmost of the equivalent placements
    let mut p = 150;
    while v[p - 1] == v[p + 35] {
        p -= 1;
    }
    assert_eq!(gap_at(&a), (20 + p, p));
    assert_eq!((a.q_start, a.q_end, a.g_start, a.g_end), (20, 20 + 264, 0, 300));
}

#[test]
fn align_insertion_27nt() {
    let v = rand_seq(300, 10);
    let (lead, ins, tail) = (rand_seq(20, 100), rand_seq(27, 101), rand_seq(30, 102));
    let read = cat(&[&lead, &v[..150], &ins, &v[150..], &tail]);
    let a = align_semiglobal(&read, &v, &S);
    check(&a, &read, &v);
    assert_eq!(gaps(&a), vec![(Op::Ins, 27)]);
    let mut s = 170;
    while read[s - 1] == read[s + 26] {
        s -= 1;
    }
    assert_eq!(gap_at(&a), (s, s - 20));
    assert_eq!((a.q_start, a.q_end, a.g_start, a.g_end), (20, 20 + 327, 0, 300));
}

#[test]
fn align_truncated_5p() {
    let v = rand_seq(300, 11);
    let tail = rand_seq(40, 110);
    let read = cat(&[&v[100..], &tail]);
    let a = align_semiglobal(&read, &v, &S);
    check(&a, &read, &v);
    assert_eq!((a.q_start, a.q_end, a.g_start, a.g_end), (0, 200, 100, 300));
    assert_eq!(a.ops, vec![(Op::Match, 200)]);
    // and 3'-truncated: the germline 3' end is free too
    let read = &v[..180];
    let a = align_semiglobal(read, &v, &S);
    assert_eq!((a.q_start, a.q_end, a.g_start, a.g_end), (0, 180, 0, 180));
}

#[test]
fn n_never_matches() {
    let mut g = rand_seq(60, 12);
    g[45] = b'N';
    let mut q = g.clone();
    q[10] = b'N';
    q[30] = b'n';
    let a = align_semiglobal(&q, &g, &S);
    assert_eq!((a.q_start, a.q_end, a.g_start, a.g_end), (0, 60, 0, 60));
    assert_eq!(
        a.ops,
        vec![
            (Op::Match, 10),
            (Op::Mismatch, 1),
            (Op::Match, 19),
            (Op::Mismatch, 1),
            (Op::Match, 14),
            (Op::Mismatch, 1), // N in the query against N in the germline
            (Op::Match, 14)
        ]
    );
    assert_eq!(a.score, 57 * 2 - 3 * 3);
}

#[test]
fn cigar_format() {
    let a = Alignment {
        score: 0,
        q_start: 5,
        q_end: 305,
        g_start: 0,
        g_end: 303,
        ops: vec![(Op::Match, 200), (Op::Mismatch, 1), (Op::Match, 89), (Op::Del, 3), (Op::Match, 10)],
    };
    assert_eq!(cigar(&a, 305), "5S290M3D10M");
    assert_eq!(cigar(&a, 310), "5S290M3D10M5S");
    let b = mk(0, 12, 7, 17, vec![(Op::Match, 5), (Op::Ins, 2), (Op::Match, 5)]);
    assert_eq!(cigar(&b, 12), "7N5M2I5M");
}

#[test]
fn deterministic() {
    let v = rand_seq(300, 13);
    let mut read = cat(&[&rand_seq(15, 130), &v[..120], &v[130..], &rand_seq(25, 131)]);
    for p in (0..read.len()).step_by(9) {
        read[p] = b"ACGT"[(p / 9) % 4]; // heavy mutation, many ties
    }
    let first = align_semiglobal(&read, &v, &S);
    let short = rand_seq(90, 14);
    for _ in 0..20 {
        assert_eq!(align_semiglobal(&read, &v, &S), first);
        align_semiglobal(&short, &v[..100], &S); // other sizes reuse the scratch buffers
    }
    let threads: Vec<_> = (0..4)
        .map(|_| {
            let (r, g) = (read.clone(), v.clone());
            std::thread::spawn(move || align_semiglobal(&r, &g, &S))
        })
        .collect();
    for t in threads {
        assert_eq!(t.join().unwrap(), first);
    }
}

#[test]
fn place_gaps_homopolymer_leftmost() {
    // germline with an A homopolymer; read lacks one A. Start with the deletion on the last A.
    let g = b"CGTACGGTCAAAAAATCGGATCCA".to_vec();
    let q = cat(&[&g[..10], &g[11..]]); // any A of g[9..15] is equivalent
    let a = mk(0, 23, 0, 24, vec![(Op::Match, 14), (Op::Del, 1), (Op::Match, 9)]);
    let p = place_gaps(&a, &q, &g, None, &S);
    check(&p, &q, &g);
    assert_eq!(p.ops, vec![(Op::Match, 9), (Op::Del, 1), (Op::Match, 14)]);
    // insertion of an extra A, placed at the end of the run, moves to its start
    let q = cat(&[&g[..15], b"A", &g[15..]]);
    let a = mk(0, 25, 0, 24, vec![(Op::Match, 15), (Op::Ins, 1), (Op::Match, 9)]);
    let p = place_gaps(&a, &q, &g, None, &S);
    check(&p, &q, &g);
    assert_eq!(p.ops, vec![(Op::Match, 9), (Op::Ins, 1), (Op::Match, 15)]);
}

#[test]
fn place_gaps_mutation_moves_deletion_left() {
    // germline TTGATGATTATG|G|C|ATGAGCT: C13 deleted and G12 mutated to C. Read as aligned with
    // the deletion on C13 has a G12>C mismatch; deleting G12 (or the equal G11) removes it.
    let g = b"TTGATGATTATGGCATGAGCT".to_vec();
    let q = cat(&[&g[..12], b"C", &g[14..]]);
    let a = mk(0, 20, 0, 21, vec![(Op::Match, 12), (Op::Mismatch, 1), (Op::Del, 1), (Op::Match, 7)]);
    let p = place_gaps(&a, &q, &g, None, &S);
    check(&p, &q, &g);
    assert_eq!(p.ops, vec![(Op::Match, 11), (Op::Del, 1), (Op::Match, 9)]);
    assert_eq!(p.score, 20 * 2 - 10);
    // the aligner itself agrees
    let b = align_semiglobal(&q, &g, &S);
    assert_eq!(b, p);
}

#[test]
fn place_gaps_respects_protect() {
    let g = b"TTGATGATTATGGCATGAGCT".to_vec();
    let q = cat(&[&g[..12], b"C", &g[14..]]);
    let a = mk(0, 20, 0, 21, vec![(Op::Match, 12), (Op::Mismatch, 1), (Op::Del, 1), (Op::Match, 7)]);
    // flanks of a deletion at 11 are germline 10 and 12; 10 is protected, so 12 is the best left
    let p = place_gaps(&a, &q, &g, Some(0..11), &S);
    check(&p, &q, &g);
    assert_eq!(p.ops, vec![(Op::Match, 12), (Op::Del, 1), (Op::Match, 8)]);
    // protecting the right side keeps every placement at or left of 10
    let p = place_gaps(&a, &q, &g, Some(12..21), &S);
    check(&p, &q, &g);
    assert!(gap_at(&p).1 + 1 < 12);
}

/// Best overlap-alignment score (scalar Gotoh): starts on the first base of the query or the
/// germline, ends on the last base of one; gaps cost as in `Scoring`; non-ACGT never matches.
fn overlap_dp(q: &[u8], g: &[u8]) -> i32 {
    let (n, m) = (q.len(), g.len());
    let neg = i32::MIN / 4;
    let (o, e) = (S.gap_open as i32, S.gap_extend as i32);
    // h: best ending in a match column; x: ending in Ins (query only); y: ending in Del
    let mut h = vec![vec![neg; m + 1]; n + 1];
    let mut x = vec![vec![neg; m + 1]; n + 1];
    let mut y = vec![vec![neg; m + 1]; n + 1];
    for row in h.iter_mut() {
        row[0] = 0;
    }
    h[0].fill(0);
    for i in 1..=n {
        for j in 1..=m {
            let best = |i: usize, j: usize| h[i][j].max(x[i][j]).max(y[i][j]);
            let s = if q[i - 1] == g[j - 1] && b"ACGT".contains(&q[i - 1]) { S.matched } else { S.mismatch } as i32;
            h[i][j] = best(i - 1, j - 1) + s;
            x[i][j] = (h[i - 1][j].max(y[i - 1][j]) + o).max(x[i - 1][j] + e);
            y[i][j] = (h[i][j - 1].max(x[i][j - 1]) + o).max(y[i][j - 1] + e);
        }
    }
    let best = |i: usize, j: usize| h[i][j].max(x[i][j]).max(y[i][j]);
    (0..=m).map(|j| best(n, j)).chain((0..=n).map(|i| best(i, m))).max().unwrap()
}

#[test]
fn matches_scalar_dp() {
    let mut dropped = 0;
    for t in 0..240u64 {
        let v = rand_seq(150 + (t as usize % 200), 1000 + t);
        let cut = (t as usize * 7) % 60;
        let mut read = match t % 4 {
            0 => [rand_seq(t as usize % 40, 1500 + t), v[cut..].to_vec()].concat(), // 5' lead, 5' truncation
            1 => v[..v.len() - 20 - cut].to_vec(),                                    // 3' germline truncation
            2 => [&v[..v.len() / 2], &rand_seq(1 + t as usize % 30, 1700 + t)[..], &v[v.len() / 2..]].concat(), // insertion
            _ => [&v[..90], &v[91 + t as usize % 13..]].concat(),                     // deletion
        };
        let mut x = t;
        for b in read.iter_mut() {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            if (x >> 60) < 3 {
                *b = b"ACGTN"[(x >> 40) as usize % 5]; // ~19% hit rate
            }
        }
        if t % 8 == 1 {
            // read ends inside the germline on junk: the optimum often ends in a trailing insertion
            read.truncate(read.len() - 30);
            read.extend(rand_seq(25, 1900 + t));
        } else {
            read.extend(rand_seq(t as usize % 50, 2000 + t));
        }
        let a = align_semiglobal(&read, &v, &S);
        check(&a, &read, &v);
        let dp = overlap_dp(&read, &v);
        // the documented rule: score is the overlap optimum; terminal gaps are dropped from the ops
        assert_eq!(a.score, dp, "case {t}");
        let (score, q_end, g_end) = score_semiglobal(&read, &v, &S);
        assert_eq!(score, dp, "case {t}");
        if ops_score(&a.ops) == a.score {
            assert_eq!((q_end, g_end), (a.q_end, a.g_end), "case {t}");
        } else {
            dropped += 1;
        }
    }
    assert!(dropped > 0, "covers a dropped terminal gap");
}

#[test]
fn empty_when_nothing_aligns() {
    let a = align_semiglobal(b"AAAAAAAA", b"CCCCCCCCCCCC", &S);
    assert_eq!(a, Alignment { score: 0, q_start: 0, q_end: 0, g_start: 0, g_end: 0, ops: vec![] });
    assert_eq!(score_semiglobal(b"AAAAAAAA", b"CCCCCCCCCCCC", &S), (0, 0, 0));
    assert_eq!(align_semiglobal(b"", b"ACGT", &S).ops, vec![]);
}

/// Gap list (is insertion, start, length) of an alignment: insertions in query coordinates,
/// deletions in germline coordinates.
fn gap_list(a: &Alignment) -> Vec<(bool, usize, usize)> {
    let (mut qi, mut gi, mut out) = (a.q_start, a.g_start, Vec::new());
    for &(op, n) in &a.ops {
        let n = n as usize;
        match op {
            Op::Ins => out.push((true, qi, n)),
            Op::Del => out.push((false, gi, n)),
            _ => {}
        }
        if op != Op::Del {
            qi += n;
        }
        if op != Op::Ins {
            gi += n;
        }
    }
    out
}

/// The generator's rule written out by brute force (crates/generate/src/derive.rs, ruling A17):
/// each gap in turn tries every placement, rebuilding the whole column map.
fn brute_place(a: &Alignment, q: &[u8], g: &[u8], protect: Option<std::ops::Range<usize>>) -> Vec<(bool, usize, usize)> {
    let mut ids = gap_list(a);
    let align = |ids: &[(bool, usize, usize)]| {
        let mut bases = (a.q_start..a.q_end).filter(|&i| !ids.iter().any(|&(ins, p, n)| ins && (p..p + n).contains(&i)));
        let mut at = vec![None; g.len()];
        for u in (a.g_start..a.g_end).filter(|&u| !ids.iter().any(|&(ins, p, n)| !ins && (p..p + n).contains(&u))) {
            at[u] = bases.next();
        }
        at
    };
    for j in 0..ids.len() {
        let (ins, p0, len) = ids[j];
        let mut best = (usize::MAX, p0);
        for p in 1..if ins { q.len() } else { g.len() } {
            if ids.iter().enumerate().any(|(k, &(i, x, n))| k != j && i == ins && x < p + len && p < x + n) {
                continue;
            }
            ids[j] = (ins, p, len);
            let at = align(&ids);
            let flanks = if ins {
                at.iter().position(|&i| i == Some(p - 1)).zip(at.iter().position(|&i| i == Some(p + len)))
            } else {
                (p + len < g.len() && at[p - 1].is_some() && at[p + len].is_some()).then_some((p - 1, p + len))
            };
            let Some((l, r)) = flanks else { continue };
            if protect.as_ref().is_some_and(|pr| l < pr.end && pr.start <= r) {
                continue;
            }
            let mm = (0..g.len()).filter(|&u| at[u].is_some_and(|i| q[i] != g[u])).count();
            best = best.min((mm, p));
        }
        ids[j] = (ins, best.1, len);
    }
    ids
}

#[test]
fn place_gaps_matches_brute_force() {
    for t in 0..200u64 {
        // a low-complexity germline so equivalent and near-equivalent placements abound
        let g: Vec<u8> = rand_seq(90, 3000 + t).iter().map(|&b| if b == b'T' { b'A' } else { b }).collect();
        let (dp, dl) = (10 + (t as usize * 7) % 30, 1 + t as usize % 4);
        let (ip, il) = (50 + (t as usize * 11) % 25, 1 + (t as usize / 3) % 4);
        // read: g[..dp] + g[dp+dl..ip] + inserted + g[ip..], with mutations
        let mut q = cat(&[&g[..dp], &g[dp + dl..ip], &rand_seq(il, 4000 + t), &g[ip..]]);
        for (k, b) in q.iter_mut().enumerate() {
            if (k as u64 * 31 + t).is_multiple_of(9) {
                *b = b"ACG"[(k + t as usize) % 3];
            }
        }
        let iq = ip - dl; // query position of the insertion
        // only gap positions matter to place_gaps, so the aligned runs need not be labelled exactly
        let ops = vec![(Op::Match, dp as u32), (Op::Del, dl as u32), (Op::Match, (iq - dp) as u32), (Op::Ins, il as u32), (Op::Match, (g.len() - ip) as u32)];
        let a = mk(0, q.len(), 0, g.len(), ops);
        for protect in [None, Some(0..30), Some(60..90)] {
            let p = place_gaps(&a, &q, &g, protect.clone(), &S);
            check(&p, &q, &g);
            assert_eq!(gap_list(&p), brute_place(&a, &q, &g, protect.clone()), "case {t} protect {protect:?}");
        }
    }
}

/// Ops labelled Match/Mismatch from the sequences for a layout of (kind, length) runs, where a
/// `Match` run means "aligned columns".
fn labelled(layout: &[(Op, u32)], q: &[u8], g: &[u8]) -> Vec<(Op, u32)> {
    let (mut qi, mut gi, mut ops) = (0, 0, Vec::<(Op, u32)>::new());
    let push = |ops: &mut Vec<(Op, u32)>, op: Op, n: u32| match ops.last_mut() {
        Some((l, m)) if *l == op => *m += n,
        _ => ops.push((op, n)),
    };
    for &(op, n) in layout {
        match op {
            Op::Ins => {
                push(&mut ops, Op::Ins, n);
                qi += n as usize;
            }
            Op::Del => {
                push(&mut ops, Op::Del, n);
                gi += n as usize;
            }
            _ => {
                for _ in 0..n {
                    push(&mut ops, if q[qi] == g[gi] { Op::Match } else { Op::Mismatch }, 1);
                    qi += 1;
                    gi += 1;
                }
            }
        }
    }
    ops
}

#[test]
fn collapse_spurious_pair() {
    // An insertion and a deletion of 3, 5 columns apart, on a read identical to the germline:
    // the diagonal is all matches and replaces the pair.
    let g = rand_seq(120, 41);
    let ops = labelled(&[(Op::Match, 40), (Op::Ins, 3), (Op::Match, 5), (Op::Del, 3), (Op::Match, 72)], &g, &g);
    let a = mk(0, 120, 0, 120, ops);
    check(&a, &g, &g);
    let c = collapse_gap_pairs(&a, &g, &g, &S);
    check(&c, &g, &g);
    assert_eq!(c.ops, vec![(Op::Match, 120)]);
    assert_eq!(c.score, 240);
    // Unbalanced or farther apart than MAX_PAIR_SPAN: kept.
    let far = labelled(&[(Op::Match, 10), (Op::Ins, 3), (Op::Match, MAX_PAIR_SPAN as u32 + 1), (Op::Del, 3), (Op::Match, 120 - 14 - MAX_PAIR_SPAN as u32)], &g, &g);
    let a = mk(0, 120, 0, 120, far);
    assert_eq!(collapse_gap_pairs(&a, &g, &g, &S), a);
}

#[test]
fn collapse_keeps_real_indels() {
    // A real 9-nt deletion, and a real insertion+deletion pair whose diagonal is far worse.
    let g = rand_seq(200, 43);
    let q = cat(&[&g[..80], &g[89..]]);
    let a = align_semiglobal(&q, &g, &S);
    assert_eq!(gaps(&a), vec![(Op::Del, 9)]);
    assert_eq!(collapse_gap_pairs(&a, &q, &g, &S), a);
    let q = cat(&[&g[..60], b"ACGTAC", &g[60..90], &g[96..]]);
    let a = align_semiglobal(&q, &g, &S);
    assert_eq!(gaps(&a), vec![(Op::Ins, 6), (Op::Del, 6)]);
    assert_eq!(collapse_gap_pairs(&a, &q, &g, &S), a);
}

#[test]
fn diagonal_hits_count_one_band() {
    // The source allele's k-mers all lie on one diagonal; a read of shuffled 20-mers of it
    // shares as many k-mers but spread over many diagonals.
    let a = rand_seq(300, 45);
    let b = rand_seq(300, 46);
    let index = KmerIndex::new(&[&a, &b], 7);
    let d = index.diagonal_hits(&a, 8);
    assert_eq!(d.iter().find(|h| h.0 == 0).unwrap().1, 294);
    let shuffled: Vec<u8> = (0..15).rev().flat_map(|i| a[i * 20..i * 20 + 20].to_vec()).collect();
    let plain = index.hits(&shuffled).iter().find(|h| h.0 == 0).unwrap().1;
    let banded = index.diagonal_hits(&shuffled, 8).iter().find(|h| h.0 == 0).unwrap().1;
    assert!(plain >= 15 * 14 && banded <= 2 * 14, "plain {plain} banded {banded}");
}
