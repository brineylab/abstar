#!usr/env/python
# filename: test_regions.py

import sys

from . import REFERENCE_BNAB_HC_ANTIBODIES, REFERENCE_BNAB_LC_ANTIBODIES
from .mock_classes import MockAntibody
from ..core.germline import GermlineSegment
from ..utils.regions import get_joining_regions, get_variable_regions


def test_var_regions_bnab_hcs():
    compare_v_regions_to_reference(REFERENCE_BNAB_HC_ANTIBODIES)


def test_join_regions_bnab_hcs():
    compare_j_regions_to_reference(REFERENCE_BNAB_HC_ANTIBODIES)


def test_var_regions_bnab_lcs():
    compare_v_regions_to_reference(REFERENCE_BNAB_LC_ANTIBODIES)


def test_join_regions_bnab_lcs():
    compare_j_regions_to_reference(REFERENCE_BNAB_LC_ANTIBODIES)


def compare_v_regions_to_reference(antibodies):
    errors = []
    for antibody in antibodies:
        ref = antibody.v.regions
        ab = MockAntibody()
        ab.oriented_input = antibody.oriented_input
        ab.v_rf_offset = antibody.v_rf_offset
        v = GermlineSegment(antibody.v.full, antibody.v.species, antibody.v.db_name, initialize_log=False)
        v.query_start = antibody.v.query_start
        v.gene_type = antibody.v.gene_type
        v._imgt_position_from_raw = antibody.v._imgt_position_from_raw
        v._raw_position_from_imgt = antibody.v._raw_position_from_imgt
        ab.v = v
        regions = get_variable_regions(ab)
        for r in regions.region_names:
            if ref.raw_nt_positions[r] != regions.raw_nt_positions[r]:
                e = '{} {}: '.format(antibody.id, r)
                e += 'reference raw_nt_positions ({}) '.format(ref.raw_nt_positions[r]())
                e += 'did not match calculated raw_nt_positions ({})'.format(regions.raw_nt_positions[r])
                errors.append(e)
            if ref.nt_seqs[r] != regions.nt_seqs[r]:
                e = '{} {}: '.format(antibody.id, r)
                e += 'reference nt_seq ({}) '.format(ref.nt_seqs[r]())
                e += 'did not match calculated nt_seq ({})'.format(regions.nt_seqs[r])
                errors.append(e)
            if ref.aa_seqs[r] != regions.aa_seqs[r]:
                e = '{} {}: '.format(antibody.id, r)
                e += 'reference aa_seq ({}) '.format(ref.aa_seqs[r]())
                e += 'did not match calculated aa_seq ({})'.format(regions.aa_seqs[r])
                errors.append(e)
    assert len(errors) == 0, '\n'.join(errors)


def compare_j_regions_to_reference(antibodies):
    errors = []
    for antibody in antibodies:
        ref = antibody.j.regions
        ab = MockAntibody()
        ab.oriented_input = antibody.oriented_input
        ab.v_rf_offset = antibody.v_rf_offset
        j = GermlineSegment(antibody.j.full, antibody.j.species, antibody.j.db_name, initialize_log=False)
        j.query_start = antibody.j.query_start
        j.gene_type = antibody.j.gene_type
        j._imgt_position_from_raw = antibody.j._imgt_position_from_raw
        j._raw_position_from_imgt = antibody.j._raw_position_from_imgt
        ab.j = j
        regions = get_joining_regions(ab)
        for r in regions.region_names:
            if ref.raw_nt_positions[r] != regions.raw_nt_positions[r]:
                e = '{} {}: '.format(antibody.id, r)
                e += 'reference raw_nt_positions ({}) '.format(ref.raw_nt_positions[r]())
                e += 'did not match calculated raw_nt_positions ({})'.format(regions.raw_nt_positions[r])
                errors.append(e)
            if ref.nt_seqs[r] != regions.nt_seqs[r]:
                e = '{} {}: '.format(antibody.id, r)
                e += 'reference nt_seq ({}) '.format(ref.nt_seqs[r]())
                e += 'did not match calculated nt_seq ({})'.format(regions.nt_seqs[r])
                errors.append(e)
            if ref.aa_seqs[r] != regions.aa_seqs[r]:
                e = '{} {}: '.format(antibody.id, r)
                e += 'reference aa_seq ({}) '.format(ref.aa_seqs[r]())
                e += 'did not match calculated aa_seq ({})'.format(regions.aa_seqs[r])
                errors.append(e)
    assert len(errors) == 0, '\n'.join(errors)
