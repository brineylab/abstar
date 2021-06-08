#!usr/env/python
# filename: test_mutations.py

import sys

from . import REFERENCE_BNAB_HC_ANTIBODIES, REFERENCE_BNAB_LC_ANTIBODIES
from .mock_classes import MockAntibody
from ..core.germline import GermlineSegment
from ..utils.mutations import aa_mutations, nt_mutations


def test_v_nt_mutations_bnab_hcs():
    compare_v_nt_mutations_to_reference(REFERENCE_BNAB_HC_ANTIBODIES)

def test_v_nt_mutations_bnab_lcs():
    compare_v_nt_mutations_to_reference(REFERENCE_BNAB_LC_ANTIBODIES)

def test_j_nt_mutations_bnab_hcs():
    compare_j_nt_mutations_to_reference(REFERENCE_BNAB_HC_ANTIBODIES)

def test_j_nt_mutations_bnab_lcs():
    compare_j_nt_mutations_to_reference(REFERENCE_BNAB_LC_ANTIBODIES)


def test_v_aa_mutations_bnab_hcs():
    compare_v_aa_mutations_to_reference(REFERENCE_BNAB_HC_ANTIBODIES)

def test_v_aa_mutations_bnab_lcs():
    compare_v_aa_mutations_to_reference(REFERENCE_BNAB_LC_ANTIBODIES)

def test_j_aa_mutations_bnab_hcs():
    compare_j_aa_mutations_to_reference(REFERENCE_BNAB_HC_ANTIBODIES)

def test_j_aa_mutations_bnab_lcs():
    compare_j_aa_mutations_to_reference(REFERENCE_BNAB_LC_ANTIBODIES)


def compare_v_nt_mutations_to_reference(antibodies):
    errors = []
    for antibody in antibodies:
        ref = antibody.v.nt_mutations
        ab = MockAntibody()
        v = GermlineSegment(antibody.v.full, antibody.v.species, antibody.v.db_name, initialize_log=False)
        v.gene_type = antibody.v.gene_type
        v.query_alignment = antibody.v.query_alignment
        v.germline_alignment = antibody.v.germline_alignment
        v._imgt_position_from_raw = antibody.v._imgt_position_from_raw
        v.query_start = antibody.v.query_start
        ab.v = v
        ab.j = None
        m = nt_mutations(ab)
        if len(ref.mutations) != len(m.mutations):
            e = '{} V-gene: '.format(antibody.id)
            e += 'reference number of nt_mutations ({}) '.format(len(ref.mutations))
            e += 'did not match calculated number of nt_mutations ({})'.format(len(m.mutations))
            errors.append(e)
        for mut in m.mutations:
            try:
                ref_mut = [_m for _m in ref.mutations if _m.raw_position == mut.raw_position][0]
            except IndexError:
                e = '{} V-gene: '.format(antibody.id)
                e += 'no reference matches were found for mutation {}'.format(mut.abstar_formatted)
                errors.append(e)
                continue
            if ref_mut.was != mut.was:
                e = '{} V-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference was nucleotide ({}) '.format(ref_mut.was)
                e += 'did not match calculated was nucleotide ({})'.format(mut.was)
                errors.append(e)
            if ref_mut.now != mut.now:
                e = '{} V-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference now nucleotide ({}) '.format(ref_mut.now)
                e += 'did not match calculated now nucleotide ({})'.format(mut.now)
                errors.append(e)
            if ref_mut.imgt_position != mut.imgt_position:
                e = '{} V-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference imgt_position ({}) '.format(ref_mut.imgt_position)
                e += 'did not match calculated imgt_position ({})'.format(mut.imgt_position)
                errors.append(e)
            if ref_mut.imgt_codon != mut.imgt_codon:
                e = '{} V-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference imgt_codon ({}) '.format(len(ref_mut.imgt_codon))
                e += 'did not match calculated imgt_codon ({})'.format(mut.imgt_codon)
                errors.append(e)
    assert len(errors) == 0, '\n'.join(errors)


def compare_j_nt_mutations_to_reference(antibodies):
    errors = []
    for antibody in antibodies:
        ref = antibody.j.nt_mutations
        ab = MockAntibody()
        j = GermlineSegment(antibody.j.full, antibody.j.species, antibody.j.db_name, initialize_log=False)
        j.gene_type = antibody.j.gene_type
        j.query_alignment = antibody.j.query_alignment
        j.germline_alignment = antibody.j.germline_alignment
        j._imgt_position_from_raw = antibody.j._imgt_position_from_raw
        j.query_start = antibody.j.query_start
        j._correct_imgt_nt_position_from_imgt = antibody.j._correct_imgt_nt_position_from_imgt
        j._initial_correct_imgt_nt_position_from_imgt = antibody.j._initial_correct_imgt_nt_position_from_imgt
        ab.j = j
        ab.v = None
        m = nt_mutations(ab)
        if len(ref.mutations) != len(m.mutations):
            e = '{} J-gene: '.format(antibody.id)
            e += 'reference number of nt_mutations ({}) '.format(len(ref.mutations))
            e += 'did not match calculated number of nt_mutations ({})'.format(len(m.mutations))
            errors.append(e)
        for mut in m.mutations:
            try:
                ref_mut = [_m for _m in ref.mutations if _m.raw_position == mut.raw_position][0]
            except IndexError:
                e = '{} J-gene: '.format(antibody.id)
                e += 'no reference matches were found for mutation {}'.format(mut.abstar_formatted)
                errors.append(e)
                continue
            if ref_mut.was != mut.was:
                e = '{} J-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference was nucleotide ({}) '.format(ref_mut.was)
                e += 'did not match calculated was nucleotide ({})'.format(mut.was)
                errors.append(e)
            if ref_mut.now != mut.now:
                e = '{} J-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference now nucleotide ({}) '.format(ref_mut.now)
                e += 'did not match calculated now nucleotide ({})'.format(mut.now)
                errors.append(e)
            if ref_mut.imgt_position != mut.imgt_position:
                e = '{} J-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference imgt_position ({}) '.format(ref_mut.imgt_position)
                e += 'did not match calculated imgt_position ({})'.format(mut.imgt_position)
                errors.append(e)
            if ref_mut.imgt_codon != mut.imgt_codon:
                e = '{} J-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference imgt_codon ({}) '.format(ref_mut.imgt_codon)
                e += 'did not match calculated imgt_codon ({})'.format(mut.imgt_codon)
                errors.append(e)
    assert len(errors) == 0, '\n'.join(errors)


def compare_v_aa_mutations_to_reference(antibodies):
    errors = []
    for antibody in antibodies:
        ref = antibody.v.aa_mutations
        ab = MockAntibody()
        v = GermlineSegment(antibody.v.full, antibody.v.species, antibody.v.db_name, initialize_log=False)
        v.gene_type = antibody.v.gene_type
        v.query_alignment = antibody.v.query_alignment
        v.germline_alignment = antibody.v.germline_alignment
        v._imgt_position_from_raw = antibody.v._imgt_position_from_raw
        v.query_start = antibody.v.query_start
        v.aa_sequence = antibody.v.aa_sequence
        ab.v = v
        ab.j = None
        m = aa_mutations(ab)

        if ab.exceptions:
            errors.append('\n'.join(ab.exceptions))

        if len(ref.mutations) != len(m.mutations):
            e = '{} V-gene: '.format(antibody.id)
            e += 'reference number of aa_mutations ({}) '.format(len(ref.mutations))
            e += 'did not match calculated number of aa_mutations ({})'.format(len(m.mutations))
            errors.append(e)
        for mut in m.mutations:
            try:
                ref_mut = [_m for _m in ref.mutations if _m.imgt_position == mut.imgt_position][0]
            except IndexError:
                e = '{} V-gene: '.format(antibody.id)
                e += 'no reference matches were found for mutation {}'.format(mut.abstar_formatted)
                errors.append(e)
                continue
            if ref_mut.was != mut.was:
                e = '{} V-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference was nucleotide ({}) '.format(ref_mut.was)
                e += 'did not match calculated was nucleotide ({})'.format(mut.was)
                errors.append(e)
            if ref_mut.now != mut.now:
                e = '{} V-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference now nucleotide ({}) '.format(ref_mut.now)
                e += 'did not match calculated now nucleotide ({})'.format(mut.now)
                errors.append(e)
            if ref_mut.imgt_position != mut.imgt_position:
                e = '{} V-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference imgt_position ({}) '.format(ref_mut.imgt_position)
                e += 'did not match calculated imgt_position ({})'.format(mut.imgt_position)
                errors.append(e)
            if ref_mut.imgt_codon != mut.imgt_codon:
                e = '{} V-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference imgt_codon ({}) '.format(ref_mut.imgt_codon)
                e += 'did not match calculated imgt_codon ({})'.format(mut.imgt_codon)
                errors.append(e)
    assert len(errors) == 0, '\n'.join(errors)


def compare_j_aa_mutations_to_reference(antibodies):
    errors = []
    for antibody in antibodies:
        ref = antibody.j.aa_mutations
        ab = MockAntibody()
        j = GermlineSegment(antibody.j.full, antibody.j.species, antibody.j.db_name, initialize_log=False)
        j.gene_type = antibody.j.gene_type
        j.query_alignment = antibody.j.query_alignment
        j.germline_alignment = antibody.j.germline_alignment
        j._imgt_position_from_raw = antibody.j._imgt_position_from_raw
        j.query_start = antibody.j.query_start
        j._correct_imgt_aa_position_from_imgt = antibody.j._correct_imgt_aa_position_from_imgt
        j._initial_correct_imgt_aa_position_from_imgt = antibody.j._initial_correct_imgt_aa_position_from_imgt
        j.aa_sequence = antibody.j.aa_sequence
        ab.j = j
        ab.v = None
        m = aa_mutations(ab)

        if ab.exceptions:
            errors.append('\n'.join(ab.exceptions))

        if len(ref.mutations) != len(m.mutations):
            e = '{} J-gene: '.format(antibody.id)
            e += 'reference number of aa_mutations ({}) '.format(len(ref.mutations))
            e += 'did not match calculated number of aa_mutations ({})'.format(len(m.mutations))
            errors.append(e)
        for mut in m.mutations:
            try:
                ref_mut = [_m for _m in ref.mutations if _m.imgt_position == mut.imgt_position][0]
            except IndexError:
                e = '{} J-gene: '.format(antibody.id)
                e += 'no reference matches were found for mutation {}'.format(mut.abstar_formatted)
                errors.append(e)
                continue
            if ref_mut.was != mut.was:
                e = '{} J-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference was nucleotide ({}) '.format(ref_mut.was)
                e += 'did not match calculated was nucleotide ({})'.format(mut.was)
                errors.append(e)
            if ref_mut.now != mut.now:
                e = '{} J-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference now nucleotide ({}) '.format(ref_mut.now)
                e += 'did not match calculated now nucleotide ({})'.format(mut.now)
                errors.append(e)
            if ref_mut.imgt_position != mut.imgt_position:
                e = '{} J-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference imgt_position ({}) '.format(ref_mut.imgt_position)
                e += 'did not match calculated imgt_position ({})'.format(mut.imgt_position)
                errors.append(e)
            if ref_mut.imgt_codon != mut.imgt_codon:
                e = '{} J-gene {}: '.format(antibody.id, ref_mut.abstar_formatted)
                e += 'reference imgt_codon ({}) '.format(ref_mut.imgt_codon)
                e += 'did not match calculated imgt_codon ({})'.format(mut.imgt_codon)
                errors.append(e)
    assert len(errors) == 0, '\n'.join(errors)