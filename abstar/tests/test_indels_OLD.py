# #!usr/env/python
# # filename: test_indels.py

# import sys

# from . import REFERENCE_BNAB_HC_ANTIBODIES, REFERENCE_BNAB_LC_ANTIBODIES
# from .mock_classes import MockAntibody
# from ..core.germline import GermlineSegment


# def test_var_ins_bnab_hcs():
#     compare_v_ins_to_reference(REFERENCE_BNAB_HC_ANTIBODIES)


# def test_var_del_bnab_hcs():
#     compare_v_del_to_reference(REFERENCE_BNAB_HC_ANTIBODIES)


# def test_var_ins_bnab_lcs():
#     compare_v_ins_to_reference(REFERENCE_BNAB_LC_ANTIBODIES)


# def test_var_del_bnab_lcs():
#     compare_v_del_to_reference(REFERENCE_BNAB_LC_ANTIBODIES)


# def compare_v_ins_to_reference(antibodies):
#     errors = []
#     for antibody in antibodies:
#         ref = antibody.v
#         v = GermlineSegment(ref.full, ref.species, ref.db_name, 'bcr', initialize_log=False)
#         v.query_alignment = ref.query_alignment
#         v.germline_alignment = ref.germline_alignment
#         v.alignment_midline = ref.alignment_midline
#         v._imgt_position_from_raw = ref._imgt_position_from_raw
#         v.query_start = ref.query_start
#         ab = MockAntibody()
#         ab.oriented_input = antibody.oriented_input
#         v._find_indels(ab)
#         v._calculate_imgt_indel_positions()
#         if ref._indel_check() != v._indel_check():
#             e = '{} {}: '.format(antibody.id, antibody.chain)
#             e += 'reference _indel_check ({}) '.format(ref._indel_check())
#             e += 'did not match calculated _indel_check ({})'.format(v._indel_check())
#             errors.append(e)
#         if len(ref.insertions) != len(v.insertions):
#             e = '{} {}: '.format(antibody.id, antibody.chain)
#             e += 'reference number of insertions ({}) '.format(len(ref.insertions))
#             e += 'did not match calculated number of insertions ({})'.format(len(v.insertions))
#             errors.append(e)
#         for ins in v.insertions:
#             try:
#                 ref_ins = [i for i in ref.insertions if i.raw_position == ins.raw_position][0]
#             except IndexError:
#                 e = '{} {}: '.format(antibody.id, antibody.chain)
#                 e += 'no reference matches were found for insertion {}'.format(ins.abstar_formatted)
#                 errors.append(e)
#                 continue
#             if ref_ins.length != ins.length:
#                 e = '{} {} {}: '.format(antibody.id, antibody.chain, ins.abstar_formatted)
#                 e += 'reference insertion length ({}) '.format(ref_ins.length)
#                 e += 'did not match calculated insertion length ({})'.format(ins.length)
#                 errors.append(e)
#             if ref_ins.sequence != ins.sequence:
#                 e = '{} {} {}: '.format(antibody.id, antibody.chain, ins.abstar_formatted)
#                 e += 'reference insertion sequence ({}) '.format(ref_ins.sequence)
#                 e += 'did not match calculated insertion sequence ({})'.format(ins.sequence)
#                 errors.append(e)
#             if ref_ins.imgt_position != ins.imgt_position:
#                 e = '{} {} {}: '.format(antibody.id, antibody.chain, ins.abstar_formatted)
#                 e += 'reference insertion imgt_position ({}) '.format(ref_ins.imgt_position)
#                 e += 'did not match calculated insertion imgt_position ({})'.format(ins.imgt_position)
#                 errors.append(e)
#             if ref_ins.imgt_codon != ins.imgt_codon:
#                 e = '{} {} {}: '.format(antibody.id, antibody.chain, ins.abstar_formatted)
#                 e += 'reference insertion imgt_codon ({}) '.format(ref_ins.imgt_codon)
#                 e += 'did not match calculated insertion imgt_codon ({})'.format(ins.imgt_codon)
#                 errors.append(e)
#     assert len(errors) == 0, '\n'.join(errors)


# def compare_v_del_to_reference(antibodies):
#     errors = []
#     for antibody in antibodies:
#         ref = antibody.v
#         v = GermlineSegment(ref.full, ref.species, ref.db_name, 'bcr', initialize_log=False)
#         v.query_alignment = ref.query_alignment
#         v.germline_alignment = ref.germline_alignment
#         v.alignment_midline = ref.alignment_midline
#         v._imgt_position_from_raw = ref._imgt_position_from_raw
#         v.query_start = ref.query_start
#         ab = MockAntibody()
#         ab.oriented_input = antibody.oriented_input
#         v._find_indels(ab)
#         v._calculate_imgt_indel_positions()
#         if ref._indel_check() != v._indel_check():
#             e = '{} {}: '.format(antibody.id, antibody.chain)
#             e += 'reference _indel_check ({}) '.format(ref._indel_check())
#             e += 'did not match calculated _indel_check ({})'.format(v._indel_check())
#             errors.append(e)
#         if len(ref.deletions) != len(v.deletions):
#             e = '{} {}: '.format(antibody.id, antibody.chain)
#             e += 'reference number of deletions ({}) '.format(len(ref.deletions))
#             e += 'did not match calculated number of deletions ({})'.format(len(v.deletions))
#             errors.append(e)
#         for ins in v.deletions:
#             try:
#                 ref_ins = [i for i in ref.deletions if i.raw_position == ins.raw_position][0]
#             except IndexError:
#                 e = '{} {}: '.format(antibody.id, antibody.chain)
#                 e += 'no reference matches were found for insertion {}'.format(ins.abstar_formatted)
#                 errors.append(e)
#                 continue
#             if ref_ins.length != ins.length:
#                 e = '{} {} {}: '.format(antibody.id, antibody.chain, ins.abstar_formatted)
#                 e += 'reference insertion length ({}) '.format(ref_ins.length)
#                 e += 'did not match calculated insertion length ({})'.format(ins.length)
#                 errors.append(e)
#             if ref_ins.sequence != ins.sequence:
#                 e = '{} {} {}: '.format(antibody.id, antibody.chain, ins.abstar_formatted)
#                 e += 'reference insertion sequence ({}) '.format(ref_ins.sequence)
#                 e += 'did not match calculated insertion sequence ({})'.format(ins.sequence)
#                 errors.append(e)
#             if ref_ins.imgt_position != ins.imgt_position:
#                 e = '{} {} {}: '.format(antibody.id, antibody.chain, ins.abstar_formatted)
#                 e += 'reference insertion imgt_position ({}) '.format(ref_ins.imgt_position)
#                 e += 'did not match calculated insertion imgt_position ({})'.format(ins.imgt_position)
#                 errors.append(e)
#             if ref_ins.imgt_codon != ins.imgt_codon:
#                 e = '{} {} {}: '.format(antibody.id, antibody.chain, ins.abstar_formatted)
#                 e += 'reference insertion imgt_codon ({}) '.format(ref_ins.imgt_codon)
#                 e += 'did not match calculated insertion imgt_codon ({})'.format(ins.imgt_codon)
#                 errors.append(e)
#     assert len(errors) == 0, '\n'.join(errors)
