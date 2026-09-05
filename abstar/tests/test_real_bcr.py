"""Fixture integrity and loader contracts; these tests never annotate reads."""

import copy
from dataclasses import FrozenInstanceError
import hashlib
import json

import pytest

from abstar.tests import corpus


PILOT_IDS = (
    'ACGATACCAGGTTTCA-1_contig_2', 'CAAGATCAGAGCTTCT-1_contig_2',
    'CCATGTCCAGTCTTCC-1_contig_1', 'CTAAGACAGCAATCTC-1_contig_2',
    'CTGTTTACAGGTGCCT-1_contig_1', 'GCTGCTTAGAAACGAG-1_contig_2',
    'GGTGCGTAGGGAAACA-1_contig_1', 'GTAACTGAGTATTGGA-1_contig_2',
)


@pytest.fixture(autouse=True)
def no_external_commands(monkeypatch):
    # Applies before the shared case fixtures load, as well as inside test bodies.
    import subprocess
    def forbidden(*args, **kwargs):
        raise AssertionError('fixture integrity must never launch an external command')
    monkeypatch.setattr(subprocess, 'Popen', forbidden)


@pytest.fixture
def schema_files(tmp_path):
    case = {
        'dataset': '00123', 'sequence_id': '10E8',
        'sequence_sha256': hashlib.sha256(b'TGTAAATTT').hexdigest(),
        'selection_reasons': ['concordant_IGK'],
        'source': {
            'publication': {'record': '19617633', 'version': 3, 'license': 'MIT'},
            'cellranger': {'v_gene': 'IGKV1-5', 'cdr3': 'CKF'},
            'alignment': {'v': {'genes': ['IGKV1-5'], 'notes': ['direct comparison']}},
        },
        'expected': {
            'locus': 'IGK', 'rev_comp': False, 'status': 'annotated',
            'v_call': ['IGKV1-5'], 'j_call': 'IGKJ1*01', 'd_call': None,
            'v_sequence_start': 0, 'v_sequence_end': 3,
            'j_sequence_start': 6, 'j_sequence_end': 9,
            'junction_start': 0, 'junction_end': 9,
            'junction': 'TGTAAATTT', 'junction_aa': 'CKF',
            'cdr3': 'AAA', 'cdr3_aa': 'K',
            'productive': True, 'productivity_issues': [],
        },
        'evidence': ['Direct V/J comparison and conserved C/F delimit CKF.'],
    }
    def write(cases=None, fasta='>10E8\nTGTAAATTT\n'):
        (tmp_path / 'cases.json').write_text(json.dumps([case] if cases is None else cases))
        (tmp_path / 'sequences.fasta').write_text(fasta)
        return tmp_path
    return case, write


def test_loader_preserves_external_identity_and_deeply_freezes(schema_files):
    raw, write = schema_files
    root = write()
    cases = corpus.load_real_bcr_cases(root)
    case = cases[0]
    assert (case.dataset, case.sequence_id, case.sequence) == ('00123', '10E8', 'TGTAAATTT')
    assert case.expected['v_call'] == ('IGKV1-5',)
    assert case.source['alignment']['v']['notes'] == ('direct comparison',)
    with pytest.raises(FrozenInstanceError):
        case.sequence = 'AAA'
    for mapping in (case.expected, case.source, case.source['alignment']['v']):
        with pytest.raises(TypeError):
            mapping['new'] = 'mutation'
    assert corpus.load_real_bcr_cases(root) is not cases
    assert corpus.load_real_bcr_cases(root)[0] is not case
    sequence = case.as_sequence()
    sequence.id = 'changed'
    assert case.as_sequence().id == '10E8'
    assert case.as_sequence().sequence == 'TGTAAATTT'


def test_loader_allows_same_contig_in_distinct_datasets(schema_files):
    raw, write = schema_files
    other = copy.deepcopy(raw)
    other['dataset'] = '00234'
    cases = corpus.load_real_bcr_cases(write([raw, other], '>10E8\nTGTAAATTT\n>10E8\nTGTAAATTT\n'))
    assert [(c.dataset, c.sequence_id) for c in cases] == [('00123', '10E8'), ('00234', '10E8')]


@pytest.mark.parametrize('mutation', [
    'wrong_id', 'missing_record', 'extra_record', 'wrong_hash', 'duplicate',
    'unknown_expected', 'no_evidence', 'blank_evidence', 'unsorted_call',
    'duplicate_call', 'empty_call', 'numeric_call', 'alleles_in_allowed_set',
    'missing_source', 'invalid_boolean', 'invalid_coordinate', 'reversed_coordinate',
    'missing_required', 'wrong_junction_slice', 'invalid_reason', 'null_v',
])
def test_loader_rejects_corrupt_or_unsupported_cases(schema_files, mutation):
    raw, write = schema_files
    cases = [raw]
    fasta = '>10E8\nTGTAAATTT\n'
    if mutation == 'wrong_id': fasta = '>other\nTGTAAATTT\n'
    elif mutation == 'missing_record': fasta = ''
    elif mutation == 'extra_record': fasta += '>other\nAAA\n'
    elif mutation == 'wrong_hash': raw['sequence_sha256'] = '0' * 64
    elif mutation == 'duplicate': cases.append(copy.deepcopy(raw)); fasta += fasta
    elif mutation == 'unknown_expected': raw['expected']['invented'] = 1
    elif mutation == 'no_evidence': raw['evidence'] = []
    elif mutation == 'blank_evidence': raw['evidence'] = [' ']
    elif mutation == 'unsorted_call': raw['expected']['v_call'] = ['IGKV3-20', 'IGKV1-5']
    elif mutation == 'duplicate_call': raw['expected']['v_call'] = ['IGKV1-5', 'IGKV1-5']
    elif mutation == 'empty_call': raw['expected']['v_call'] = []
    elif mutation == 'numeric_call': raw['expected']['v_call'] = 5
    elif mutation == 'alleles_in_allowed_set': raw['expected']['v_call'] = ['IGKV1-5*01']
    elif mutation == 'missing_source': raw['source'] = {}
    elif mutation == 'invalid_boolean': raw['expected']['productive'] = 'true'
    elif mutation == 'invalid_coordinate': raw['expected']['v_sequence_start'] = True
    elif mutation == 'reversed_coordinate': raw['expected']['j_sequence_end'] = 5
    elif mutation == 'missing_required': del raw['expected']['junction']
    elif mutation == 'wrong_junction_slice': raw['expected']['junction'] = 'TGTCCCTTT'
    elif mutation == 'invalid_reason': raw['expected']['productivity_issues'] = ['']
    elif mutation == 'null_v': raw['expected']['v_call'] = None
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write(cases, fasta))


def test_committed_fixture_integrity(real_bcr_cases, pilot_loss_cases):
    cases = corpus.load_real_bcr_cases()
    assert 25 <= len(cases) <= 40
    assert len({(c.dataset, c.sequence_id) for c in cases}) == len(cases)
    assert {c.sequence_id for c in pilot_loss_cases} == set(PILOT_IDS)
    assert all(c.dataset == '1279068' for c in pilot_loss_cases)
    assert {c.expected['locus'] for c in cases} == {'IGH', 'IGK', 'IGL'}
    assert isinstance(real_bcr_cases, tuple) and isinstance(pilot_loss_cases, tuple)
    assert cases is not real_bcr_cases
    for case in cases:
        assert hashlib.sha256(case.sequence.encode('ascii')).hexdigest() == case.sequence_sha256
        assert case.source['publication']['record'] == '19617633'
        assert case.source['publication']['version'] == 3
        assert case.source['publication']['license'] == 'MIT'
        assert case.evidence and case.selection_reasons


def test_constructor_copies_nested_caller_owned_mappings(schema_files):
    raw, write = schema_files
    case = corpus.RealBCRCase(sequence='TGTAAATTT', **raw)
    raw['source']['alignment']['v']['notes'].append('mutated')
    raw['expected']['v_call'].append('IGKV3-20')
    assert case.source['alignment']['v']['notes'] == ('direct comparison',)
    assert case.expected['v_call'] == ('IGKV1-5',)


@pytest.mark.parametrize('field,value', [
    ('v_insertions', 'arbitrary'),
    ('v_deletions', [{'query_start': -1, 'query_end': -1, 'sequence': 'A'}]),
])
def test_loader_rejects_malformed_indel_evidence(schema_files, field, value):
    raw, write = schema_files
    raw['expected'][field] = value
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write())


def test_loader_checks_translation_in_addition_to_matching_slices(schema_files):
    raw, write = schema_files
    raw['expected']['junction_aa'] = 'CRF'
    raw['expected']['cdr3_aa'] = 'R'
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write())


def test_loader_uses_oriented_full_input_for_reverse_case(schema_files):
    raw, write = schema_files
    raw['expected']['rev_comp'] = True
    raw['sequence_sha256'] = hashlib.sha256(b'AAATTTACA').hexdigest()
    case = corpus.load_real_bcr_cases(write(fasta='>10E8\nAAATTTACA\n'))[0]
    assert case.expected['junction'] == 'TGTAAATTT'
    assert case.sequence == 'AAATTTACA'


def test_committed_evidence_and_translations_match_original_reads(real_bcr_cases):
    from Bio.Seq import Seq
    for case in real_bcr_cases:
        expected = case.expected
        assert str(Seq(expected['junction']).translate()) == expected['junction_aa']
        for segment in ('v', 'j'):
            alignment = case.source['alignment'][segment]
            assert alignment['query_aligned'].replace('-', '') == case.sequence[
                alignment['query_start']:alignment['query_end']]
        for insertion in expected.get('v_insertions', ()):
            assert case.sequence[insertion['query_start']:insertion['query_end']] == insertion['sequence']
    failed_anchors = {case.sequence_id: (case.expected['junction_aa'], case.expected['productivity_issues'])
                      for case in real_bcr_cases if not case.expected['productive']}
    assert failed_anchors == {
        'ATCATCTTCAGCAACT-1_contig_1': ('GARDEGWSGCSESYCSSYRIDFDYW', ('junction does not start with conserved C',)),
        'CTAAGACAGCAATCTC-1_contig_2': ('CQQYGHTSSI', ('junction does not end with conserved F',)),
        'CTGTTTACAGGTGCCT-1_contig_1': ('CTRDIKEV', ('junction does not end with conserved W',)),
    }


def test_constructor_freezes_non_dict_mapping_containers(schema_files):
    from collections import UserDict
    raw, write = schema_files
    notes = UserDict({'notes': ['original']})
    raw['source'] = UserDict({'nested': notes})
    case = corpus.RealBCRCase(sequence='TGTAAATTT', **raw)
    notes['notes'].append('changed')
    assert case.source['nested']['notes'] == ('original',)
    with pytest.raises(TypeError):
        case.source['nested']['new'] = 'changed'


def test_loader_accepts_exact_constant_alleles(schema_files):
    call = 'IGKC*01'
    raw, write = schema_files
    raw['expected']['c_call'] = call
    case = corpus.load_real_bcr_cases(write())[0]
    assert case.expected['c_call'] == call
