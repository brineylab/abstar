# #!usr/env/python
# # filename: test_junction.py

# import sys

# from . import REFERENCE_BNAB_HC_ANTIBODIES, REFERENCE_BNAB_LC_ANTIBODIES
# from ..utils.junction import get_junction

# if sys.version_info[0] > 2:
#     STR = str
# else:
#     STR = basestring


# def test_junction_bnab_hcs():
#     compare_junction_to_reference(REFERENCE_BNAB_HC_ANTIBODIES)


# def test_junction_bnab_lcs():
#     compare_junction_to_reference(REFERENCE_BNAB_LC_ANTIBODIES)


# def compare_junction_to_reference(antibodies):
#     errors = []
#     for ab in antibodies:
#         ab.j.imgt_aa_positions = [p if not isinstance(p, STR) else None for p in ab.j.imgt_aa_positions]
#         ab.j.imgt_nt_positions = [p if not isinstance(p, STR) else None for p in ab.j.imgt_nt_positions]
#         # if sys.version_info[0] == 2:
#         #     ab.oriented_input.sequence = ab.oriented_input.sequence.encode('ascii')
#         #     ab.v.imgt_germline.gapped_nt_sequence = ab.v.imgt_germline.gapped_nt_sequence.encode('ascii')
#         #     ab.j.imgt_germline.gapped_nt_sequence = ab.j.imgt_germline.gapped_nt_sequence.encode('ascii')
#         ref = ab.junction
#         j = get_junction(ab)
#         if ref.in_frame != j.in_frame:
#             e = '{}: '.format(ab.id)
#             e += 'reference in_frame ({}) '.format(ref.in_frame)
#             e += 'did not match calculated in_frame ({})'.format(j.in_frame)
#             errors.append(e)
#         if ref.junction_nt_start != j.junction_nt_start:
#             e = '{}: '.format(ab.id)
#             e += 'reference junction_nt_start ({}) '.format(ref.junction_nt_start)
#             e += 'did not match calculated junction_nt_start ({})'.format(j.junction_nt_start)
#             errors.append(e)
#         if ref.junction_nt_end != j.junction_nt_end:
#             e = '{}: '.format(ab.id)
#             e += 'reference junction_nt_end ({}) '.format(ref.junction_nt_end)
#             e += 'did not match calculated junction_nt_end ({})'.format(j.junction_nt_end)
#             errors.append(e)
#         if ref.junction_nt != j.junction_nt:
#             e = '{}: '.format(ab.id)
#             e += 'reference junction_nt ({}) '.format(ref.junction_nt)
#             e += 'did not match calculated junction_nt ({})'.format(j.junction_nt)
#             errors.append(e)
#         if ref.junction_aa != j.junction_aa:
#             e = '{}: '.format(ab.id)
#             e += 'reference junction_aa ({}) '.format(ref.junction_aa)
#             e += 'did not match calculated junction_aa ({})'.format(j.junction_aa)
#             errors.append(e)
#         if ref.cdr3_nt != j.cdr3_nt:
#             e = '{}: '.format(ab.id)
#             e += 'reference cdr3_nt ({}) '.format(ref.cdr3_nt)
#             e += 'did not match calculated cdr3_nt ({})'.format(j.cdr3_nt)
#             errors.append(e)
#         if ref.cdr3_aa != j.cdr3_aa:
#             e = '{}: '.format(ab.id)
#             e += 'reference cdr3_aa ({}) '.format(ref.cdr3_aa)
#             e += 'did not match calculated cdr3_aa ({})'.format(j.cdr3_aa)
#             errors.append(e)
#         if ref.v_nt != j.v_nt:
#             e = '{}: '.format(ab.id)
#             e += 'reference v_nt ({}) '.format(ref.v_nt)
#             e += 'did not match calculated v_nt ({})'.format(j.v_nt)
#             errors.append(e)
#         if ref.d_nt != j.d_nt:
#             e = '{}: '.format(ab.id)
#             e += 'reference d_nt ({}) '.format(ref.d_nt)
#             e += 'did not match calculated d_nt ({})'.format(j.d_nt)
#             errors.append(e)
#         if ref.j_nt != j.j_nt:
#             e = '{}: '.format(ab.id)
#             e += 'reference j_nt ({}) '.format(ref.j_nt)
#             e += 'did not match calculated j_nt ({})'.format(j.j_nt)
#             errors.append(e)
#         if ab.d is not None:
#             if ref.n1_nt != j.n1_nt:
#                 e = '{}: '.format(ab.id)
#                 e += 'reference n1_nt ({}) '.format(ref.n1_nt)
#                 e += 'did not match calculated n1_nt ({})'.format(j.n1_nt)
#                 errors.append(e)
#             if ref.n2_nt != j.n2_nt:
#                 e = '{}: '.format(ab.id)
#                 e += 'reference n2_nt ({}) '.format(ref.n2_nt)
#                 e += 'did not match calculated n2_nt ({})'.format(j.n2_nt)
#                 errors.append(e)
#         else:
#             if ref.n_nt != j.n_nt:
#                 e = '{}: '.format(ab.id)
#                 e += 'reference n_nt ({}) '.format(ref.n_nt)
#                 e += 'did not match calculated n_nt ({})'.format(j.n_nt)
#                 errors.append(e)
#         if ref.cdr3_imgt_nt_numbering != j.cdr3_imgt_nt_numbering:
#             e = '{}: '.format(ab.id)
#             e += 'reference cdr3_imgt_nt_numbering ({}) '.format(ref.cdr3_imgt_nt_numbering)
#             e += 'did not match calculated cdr3_imgt_nt_numbering ({})'.format(j.cdr3_imgt_nt_numbering)
#             errors.append(e)
#         if ref.cdr3_imgt_aa_numbering != j.cdr3_imgt_aa_numbering:
#             e = '{}: '.format(ab.id)
#             e += 'reference cdr3_imgt_aa_numbering ({}) '.format(ref.cdr3_imgt_aa_numbering)
#             e += 'did not match calculated cdr3_imgt_aa_numbering ({})'.format(j.cdr3_imgt_aa_numbering)
#             errors.append(e)
#         if ref.junction_imgt_nt_numbering != j.junction_imgt_nt_numbering:
#             e = '{}: '.format(ab.id)
#             e += 'reference junction_imgt_nt_numbering ({}) '.format(ref.junction_imgt_nt_numbering)
#             e += 'did not match calculated junction_imgt_nt_numbering ({})'.format(j.junction_imgt_nt_numbering)
#             errors.append(e)
#         if ref.junction_imgt_aa_numbering != j.junction_imgt_aa_numbering:
#             e = '{}: '.format(ab.id)
#             e += 'reference junction_imgt_aa_numbering ({}) '.format(ref.junction_imgt_aa_numbering)
#             e += 'did not match calculated junction_imgt_aa_numbering ({})'.format(j.junction_imgt_aa_numbering)
#             errors.append(e)
#     assert len(errors) == 0, '\n'.join(errors)
