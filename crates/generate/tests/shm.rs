use abstar_generate::{load_s5f_dir, motif_model, mutate, ContextModel, Selection};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::fs;

const BASES: &[u8; 4] = b"ACGT";

fn kmer(i: usize) -> Vec<u8> {
    (0..5).map(|j| BASES[(i >> (2 * (4 - j))) & 3]).collect()
}

fn random_seq(n: usize, seed: u64) -> Vec<u8> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    (0..n).map(|_| BASES[rng.random_range(0..4)]).collect()
}

fn complete(m: &ContextModel) {
    for i in 0..1024 {
        let k = kmer(i);
        let mu = m.mutability(&k).unwrap();
        assert!(mu.is_finite() && mu > 0.0, "{}", String::from_utf8_lossy(&k));
        let s = m.substitution(&k).unwrap();
        assert!((s.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        assert_eq!(s[BASES.iter().position(|&b| b == k[2]).unwrap()], 0.0);
    }
}

#[test]
fn motif_model_complete() {
    complete(motif_model());
}

#[test]
fn motif_weights() {
    let m = motif_model();
    assert_eq!(m.mutability(b"AGCTA"), Some(6.0)); // WRCY
    assert_eq!(m.mutability(b"TAGCT"), Some(6.0)); // RGYW
    assert_eq!(m.mutability(b"AGCGA"), Some(4.0)); // WRC, not WRCY
    assert_eq!(m.mutability(b"ACGCA"), Some(4.0)); // GYW, not RGYW
    assert_eq!(m.mutability(b"CTACC"), Some(2.5)); // WA
    assert_eq!(m.mutability(b"CCTAC"), Some(2.5)); // TW
    assert_eq!(m.mutability(b"GTCAA"), Some(0.3)); // SYC
    assert_eq!(m.mutability(b"AAGAC"), Some(0.3)); // GRS
    assert_eq!(m.mutability(b"ACATA"), Some(1.0));
    assert_eq!(m.mutability(b"ACNTA"), None);
    assert_eq!(m.substitution(b"AACAA"), Some([0.2, 0.0, 0.2, 0.6]));
    assert_eq!(m.substitution(b"AAGAA"), Some([0.6, 0.2, 0.0, 0.2]));
    assert_eq!(m.substitution(b"AAAAA"), Some([0.0, 0.25, 0.5, 0.25]));
    assert_eq!(m.substitution(b"AATAA"), Some([0.25, 0.5, 0.25, 0.0]));
}

fn v_divergence(orig: &[u8], seq: &[u8], v: &[usize]) -> usize {
    v.iter().filter(|&&p| orig[p] != seq[p]).count()
}

#[test]
fn realized_v_divergence_hits_target() {
    // 300 nt V followed by 60 nt of non-V junction, all mutable
    let orig = random_seq(360, 1);
    let v: Vec<usize> = (0..300).collect();
    let mutable = vec![true; 360];
    for (i, &target) in [0.0, 0.05, 0.20, 0.35, 0.50].iter().enumerate() {
        let mut seq = orig.clone();
        let mut rng = ChaCha8Rng::seed_from_u64(100 + i as u64);
        let muts = mutate(&mut seq, &mutable, &v, target, motif_model(), None, &mut rng);
        let frac = v_divergence(&orig, &seq, &v) as f64 / 300.0;
        assert!(frac >= target - 1e-12 && frac < target + 2.0 / 300.0, "target {target}: {frac}");
        assert_eq!(muts.is_empty(), target == 0.0);
        for m in &muts {
            assert_ne!(m.from, m.to);
        }
    }
}

#[test]
fn hotspots_mutate_more() {
    let orig = random_seq(300, 7);
    let v: Vec<usize> = (0..300).collect();
    let mutable = vec![true; 300];
    let w = |b: u8| b == b'A' || b == b'T';
    let r = |b: u8| b == b'A' || b == b'G';
    let s = |b: u8| b == b'C' || b == b'G';
    let y = |b: u8| b == b'C' || b == b'T';
    let wrc: Vec<usize> = (2..300).filter(|&i| orig[i] == b'C' && r(orig[i - 1]) && w(orig[i - 2])).collect();
    let syc: Vec<usize> = (2..300).filter(|&i| orig[i] == b'C' && y(orig[i - 1]) && s(orig[i - 2])).collect();
    assert!(wrc.len() >= 5 && syc.len() >= 5, "{} WRC, {} SYC", wrc.len(), syc.len());

    let mut hits = vec![0u64; 300];
    let mut rng = ChaCha8Rng::seed_from_u64(42);
    for _ in 0..2000 {
        let mut seq = orig.clone();
        for m in mutate(&mut seq, &mutable, &v, 0.05, motif_model(), None, &mut rng) {
            hits[m.pos] += 1;
        }
    }
    let mean = |ps: &[usize]| ps.iter().map(|&p| hits[p] as f64).sum::<f64>() / ps.len() as f64;
    let (hot, cold) = (mean(&wrc), mean(&syc));
    assert!(hot > 0.0 && hot >= 2.0 * cold, "WRC {hot} vs SYC {cold}");
}

#[test]
fn non_mutable_positions_untouched() {
    let orig = random_seq(330, 3);
    let v: Vec<usize> = (0..300).collect();
    let mutable: Vec<bool> = (0..330).map(|i| i % 3 != 0).collect();
    let mut seq = orig.clone();
    let mut rng = ChaCha8Rng::seed_from_u64(5);
    let muts = mutate(&mut seq, &mutable, &v, 0.2, motif_model(), None, &mut rng);
    assert!(!muts.is_empty());
    assert!(muts.iter().all(|m| mutable[m.pos]));
    for i in (0..330).step_by(3) {
        assert_eq!(seq[i], orig[i], "position {i}");
    }
}

#[test]
fn non_acgt_bases_untouched_and_case_kept() {
    let mut orig = random_seq(300, 9);
    orig[10] = b'N';
    orig[20] = b'R';
    orig[100..200].make_ascii_lowercase();
    let v: Vec<usize> = (0..300).collect();
    let mut seq = orig.clone();
    let mut rng = ChaCha8Rng::seed_from_u64(9);
    let muts = mutate(&mut seq, &[true; 300], &v, 0.5, motif_model(), None, &mut rng);
    assert_eq!((seq[10], seq[20]), (b'N', b'R'));
    for (i, b) in seq.iter().enumerate() {
        assert_eq!(b.is_ascii_lowercase(), (100..200).contains(&i), "position {i}");
    }
    assert!(muts.iter().any(|m| m.from.is_ascii_lowercase()));
}

#[test]
fn safety_cap_stops_unreachable_target() {
    // V positions are not mutable, so the target can never be reached
    let orig = random_seq(330, 4);
    let v: Vec<usize> = (0..300).collect();
    let mutable: Vec<bool> = (0..330).map(|i| i >= 300).collect();
    let mut seq = orig.clone();
    let mut rng = ChaCha8Rng::seed_from_u64(4);
    let muts = mutate(&mut seq, &mutable, &v, 0.1, motif_model(), None, &mut rng);
    assert_eq!(muts.len(), 10 * 30);
    assert_eq!(&seq[..300], &orig[..300]);
    // nothing mutable at all: stops at once
    let muts = mutate(&mut seq, &[false; 330], &v, 0.1, motif_model(), None, &mut rng);
    assert!(muts.is_empty());
}

#[test]
fn deterministic_for_seed() {
    let orig = random_seq(360, 11);
    let v: Vec<usize> = (0..300).collect();
    let mutable = vec![true; 360];
    let run = |seed: u64| {
        let mut seq = orig.clone();
        let muts = mutate(&mut seq, &mutable, &v, 0.2, motif_model(), None, &mut ChaCha8Rng::seed_from_u64(seed));
        (seq, muts.iter().map(|m| (m.pos, m.from, m.to)).collect::<Vec<_>>())
    };
    assert_eq!(run(1), run(1));
    assert_ne!(run(1), run(2));
}

/// Synthetic tables in the Kleinstein-lab layout (space-separated, quoted values).
fn write_tables(dir: &std::path::Path, drop_last_substitution: bool) {
    let mut mu = String::from("Fivemer Mutability Source lower25 upper25\n");
    // columns deliberately out of A/C/G/T order: they must be read by name
    let mut sub = String::from("Fivemer T G C A Source\n");
    let n_sub = if drop_last_substitution { 1023 } else { 1024 };
    // an N-containing row, as some exports have, is ignored
    mu.push_str("\"NNCNN\" \"9\" \"Inferred\" \"0\" \"1\"\n");
    for i in 0..1024 {
        let k = String::from_utf8(kmer(i)).unwrap();
        mu.push_str(&format!("\"{k}\" \"{}\" \"Inferred\" \"0\" \"1\"\n", (i + 1) as f64 / 1000.0));
        if i < n_sub {
            let row = synthetic_sub(k.as_bytes()[2]);
            sub.push_str(&format!("\"{k}\" \"{}\" \"{}\" \"{}\" \"{}\" \"Measured\"\n", row[3], row[2], row[1], row[0]));
        }
    }
    fs::write(dir.join("Mutability.csv"), mu).unwrap();
    fs::write(dir.join("Substitution.csv"), sub).unwrap();
}

/// [A, C, G, T] with self = 0 and the others 0.5, 0.3, 0.2 in A/C/G/T order.
fn synthetic_sub(centre: u8) -> [f64; 4] {
    let mut row = [0.0; 4];
    let mut vals = [0.5, 0.3, 0.2].into_iter();
    for (j, &b) in BASES.iter().enumerate() {
        if b != centre {
            row[j] = vals.next().unwrap();
        }
    }
    row
}

#[test]
fn s5f_dir_round_trip() {
    let dir = tempfile::tempdir().unwrap();
    write_tables(dir.path(), false);
    let m = load_s5f_dir(dir.path()).unwrap();
    complete(&m);
    for i in [0, 37, 512, 1023] {
        let k = kmer(i);
        assert_eq!(m.mutability(&k), Some((i + 1) as f64 / 1000.0));
        assert_eq!(m.substitution(&k), Some(synthetic_sub(k[2])));
    }

    let missing = tempfile::tempdir().unwrap();
    let e = load_s5f_dir(missing.path()).unwrap_err();
    assert_eq!(e.code(), "INVALID_ARGUMENT");
    assert!(e.to_string().contains("Mutability.csv"), "{e}");

    let incomplete = tempfile::tempdir().unwrap();
    write_tables(incomplete.path(), true);
    let e = load_s5f_dir(incomplete.path()).unwrap_err();
    assert_eq!(e.code(), "INVALID_ARGUMENT");
    assert!(e.to_string().contains("Substitution.csv"), "{e}");
}

#[test]
fn selection_keeps_frame_stop_free_and_anchors_conserved() {
    let in_frame_stops = |s: &[u8], f: usize| {
        s[f..].as_chunks::<3>().0.iter().filter(|c| matches!(&c.to_ascii_uppercase()[..], b"TAA" | b"TAG" | b"TGA")).count()
    };
    let cys_ok = |s: &[u8], c: usize| matches!(&s[c..c + 3], b"TGT" | b"TGC");
    let trp_ok = |s: &[u8], c: usize| matches!(&s[c..c + 3], b"TGG" | b"TTT" | b"TTC");
    let (mut neutral_stops, mut neutral_anchor_loss, mut synonymous, mut rejected_any) = (0, 0, 0, false);
    for seed in 0..50u64 {
        let f = (seed % 3) as usize;
        // 300 nt V (Cys codon ending it) + 60 nt junction with a Trp codon, stop-free in frame f
        let mut orig = random_seq(360, 1000 + seed);
        for c in (f..358).step_by(3) {
            if in_frame_stops(&orig[c..c + 3], 0) > 0 {
                orig[c] = b'C';
            }
        }
        let (cys, trp) = (f + 3 * 98, f + 3 * 112);
        orig[cys..cys + 3].copy_from_slice(b"TGT");
        orig[trp..trp + 3].copy_from_slice(b"TGG");
        let sel = Selection { frame: f, cys: Some(cys), trp: Some(trp) };
        let v: Vec<usize> = (0..300).collect();
        let mutable = vec![true; 360];
        let target = 0.15 + 0.2 * (seed as f64 / 50.0);
        let mut seq = orig.clone();
        let muts = mutate(&mut seq, &mutable, &v, target, motif_model(), Some(&sel), &mut ChaCha8Rng::seed_from_u64(seed));
        assert_eq!(in_frame_stops(&seq, f), 0, "seed {seed}");
        assert!(cys_ok(&seq, cys) && trp_ok(&seq, trp), "seed {seed}: anchors {:?}", (&seq[cys..cys + 3], &seq[trp..trp + 3]));
        synonymous += usize::from(seq[cys..cys + 3] != orig[cys..cys + 3]);
        let frac = v_divergence(&orig, &seq, &v) as f64 / 300.0;
        assert!(frac >= target - 0.5 / 300.0 && frac < target + 0.5 / 300.0, "seed {seed}: {frac} vs {target}");
        // same draws without selection: diverge at the first rejected draw
        let mut free = orig.clone();
        let n = mutate(&mut free, &mutable, &v, target, motif_model(), None, &mut ChaCha8Rng::seed_from_u64(seed));
        rejected_any |= n != muts;
        neutral_stops += in_frame_stops(&free, f);
        neutral_anchor_loss += usize::from(!cys_ok(&free, cys)) + usize::from(!trp_ok(&free, trp));
    }
    eprintln!("neutral: {neutral_stops} stops, {neutral_anchor_loss} anchors lost; selected: {synonymous} TGT->TGC");
    assert!(rejected_any && neutral_stops > 50 && neutral_anchor_loss > 10, "selection never engaged");
    assert!(synonymous > 0, "synonymous anchor changes must stay allowed");
}
