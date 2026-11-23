use abstar_annotate::align::{align_protein, blosum62, Op, ProteinScoring};
const S: ProteinScoring = ProteinScoring { gap_open: -11, gap_extend: -1 };

#[test] fn blosum62_values() {
    assert_eq!((blosum62(b'W', b'W'), blosum62(b'A', b'A'), blosum62(b'C', b'C')), (11, 4, 9));
    assert_eq!((blosum62(b'a', b'R'), blosum62(b'*', b'*'), blosum62(b'J', b'A')), (-1, 1, blosum62(b'X', b'A')));
}
#[test] fn identical_is_one_match_run() {
    let a = align_protein(b"QVQLVQSGAEV", b"QVQLVQSGAEV", &S);
    assert_eq!((a.q_start, a.q_end, a.g_start, a.g_end, a.ops.clone()), (0, 11, 0, 11, vec![(Op::Match, 11)]));
}
#[test] fn free_ends_on_both_sides() {
    // leader before the germline, germline running past the query's end
    let a = align_protein(b"MGWSQVQLVQSGAEV", b"QVQLVQSGAEVKKPG", &S);
    assert_eq!((a.q_start, a.q_end, a.g_start, a.g_end), (4, 15, 0, 11));
}
#[test] fn deletion_in_a_repeat_goes_leftmost() {
    let a = align_protein(b"MKVLLAAWGQGTLV", b"MKVLLLAAWGQGTLV", &S);
    assert_eq!(a.ops, vec![(Op::Match, 3), (Op::Del, 1), (Op::Match, 11)]);
}
#[test] fn insertion_in_a_repeat_goes_leftmost() {
    let a = align_protein(b"MKVLLLAAWGQGTLV", b"MKVLLAAWGQGTLV", &S);
    assert_eq!(a.ops, vec![(Op::Match, 3), (Op::Ins, 1), (Op::Match, 11)]);
}
#[test] fn x_never_matches() {
    let a = align_protein(b"QVQXVQSGAEV", b"QVQLVQSGAEV", &S);
    assert_eq!(a.ops, vec![(Op::Match, 3), (Op::Mismatch, 1), (Op::Match, 7)]);
}
#[test] fn unrelated_is_empty() {
    assert_eq!(align_protein(b"WWWW", b"PPPP", &S).ops, vec![]);
}
