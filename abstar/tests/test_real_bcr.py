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
            'cellranger': {'chain': 'IGK', 'v_gene': 'IGKV1-5', 'j_gene': 'IGKJ1',
                           'd_gene': '', 'cdr3': 'CKF', 'cdr3_nt': 'TGTAAATTT', 'productive': 'true'},
            'abstar_evidence': {'status': 'annotated',
                'raw_calls': {'abstar': {'v': 'IGKV1-5*01', 'j': 'IGKJ1*01', 'd': None}},
                'junction': {'abstar_nt': 'TGTAAATTT', 'abstar_aa': 'CKF'},
                'productivity': {'abstar': True, 'cellranger': True}},
            'alignment': {
                'v': {'reference': 'IGKV1-5*01', 'query_start': 0, 'query_end': 3,
                      'germline_start': 0, 'germline_end': 3, 'query_aligned': 'TGT',
                      'germline_aligned': 'TGT', 'indels': [], 'notes': ['direct comparison']},
                'j': {'reference': 'IGKJ1*01', 'query_start': 6, 'query_end': 9,
                      'germline_start': 0, 'germline_end': 3, 'query_aligned': 'TTT',
                      'germline_aligned': 'TTT', 'indels': []},
                'v_imgt104_query_start': 0, 'v_imgt104_ungapped_offset': 0,
                'j_anchor_query_start': 6, 'j_anchor_germline_offset': 0,
                'j_anchor_germline_codon': 'TTT', 'coding_start': 0, 'coding_end': 9,
                'coding_scope': 'through_primary_j', 'coding_translation': 'CKF',
            },
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


@pytest.mark.parametrize('locus,field,gene', [
    ('IGH', 'v_call', 'IGKV1-5'), ('IGH', 'v_call', 'IGHJ4'),
    ('IGH', 'd_call', 'IGHD'), ('IGH', 'd_call', 'IGKV1-5'),
    ('IGH', 'j_call', 'IGKJ1'), ('IGH', 'j_call', 'IGHV3-23'),
    ('IGH', 'c_call', 'IGKC'), ('IGH', 'c_call', 'IGHD3-10'),
    ('IGK', 'v_call', 'IGLV2-14'), ('IGK', 'v_call', 'IGKJ1'),
    ('IGK', 'd_call', 'IGHD3-10'), ('IGK', 'd_call', 'IGKD1-1'),
    ('IGK', 'j_call', 'IGHJ4'), ('IGK', 'j_call', 'IGKV1-5'),
    ('IGK', 'c_call', 'IGHM'), ('IGK', 'c_call', 'IGKJ1'),
    ('IGL', 'v_call', 'IGHV3-23'), ('IGL', 'v_call', 'IGLJ1'),
    ('IGL', 'd_call', 'IGHD3-10'), ('IGL', 'd_call', 'IGLD1-1'),
    ('IGL', 'j_call', 'IGKJ1'), ('IGL', 'j_call', 'IGLV2-14'),
    ('IGL', 'c_call', 'IGKC'), ('IGL', 'c_call', 'IGLV2-14'),
    ('IGK', 'v_call', 'IGKVjunk'), ('IGK', 'v_call', 'IGKV1--5'),
    ('IGH', 'd_call', 'IGHD0-0'), ('IGH', 'j_call', 'IGHJjunk'),
    ('IGH', 'c_call', 'IGHG9'), ('IGL', 'c_call', 'IGLCjunk'),
])
@pytest.mark.parametrize('exact', [False, True])
def test_loader_rejects_incompatible_locus_segment_and_names(schema_files, locus, field, gene, exact):
    raw, write = schema_files
    fasta = set_schema_locus(raw, locus)
    corpus.load_real_bcr_cases(write(fasta=fasta))
    raw['expected'][field] = gene + '*01' if exact else [gene]
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write(fasta=fasta))


def set_schema_locus(raw, locus):
    """Make a consistent tiny schema example before testing one corrupt call."""
    v = {'IGH': 'IGHV3-23', 'IGK': 'IGKV1-5', 'IGL': 'IGLV2-14'}[locus]
    j = locus + 'J1'
    codon = 'TGG' if locus == 'IGH' else 'TTT'
    sequence = 'TGTAAA' + codon
    aa = 'CKW' if locus == 'IGH' else 'CKF'
    expected, source = raw['expected'], raw['source']
    expected.update(locus=locus, v_call=[v], j_call=j+'*01', junction=sequence, junction_aa=aa)
    source['cellranger'].update(chain=locus, v_gene=v, j_gene=j, cdr3=aa, cdr3_nt=sequence)
    source['abstar_evidence']['raw_calls']['abstar'].update(v=v+'*01', j=j+'*01')
    source['abstar_evidence']['junction'].update(abstar_nt=sequence, abstar_aa=aa)
    source['alignment']['v']['reference'] = v+'*01'
    source['alignment']['j'].update(reference=j+'*01', query_aligned=codon, germline_aligned=codon)
    source['alignment'].update(j_anchor_germline_codon=codon, coding_translation=aa)
    raw['selection_reasons'] = ['concordant_'+locus]
    raw['sequence_sha256'] = hashlib.sha256(sequence.encode()).hexdigest()
    return '>10E8\n'+sequence+'\n'


@pytest.mark.parametrize('locus,field,call', [
    ('IGH', 'v_call', 'IGHV3-23*01'), ('IGH', 'd_call', 'IGHD3-10*01'),
    ('IGH', 'j_call', 'IGHJ1*01'), ('IGH', 'c_call', 'IGHD*01'),
    ('IGH', 'c_call', 'IGHM*01'), ('IGH', 'c_call', 'IGHG4A*01'),
    ('IGH', 'c_call', ['IGHA1', 'IGHA2']), ('IGH', 'c_call', ['IGHE']),
    ('IGK', 'v_call', 'IGKV1-5*01'), ('IGK', 'c_call', 'IGKC*01'),
    ('IGL', 'v_call', 'IGLV2-14*01'), ('IGL', 'c_call', 'IGLC1*01'),
    ('IGL', 'c_call', ['IGLC2', 'IGLC3']),
])
def test_loader_retains_compatible_call_naming(schema_files, locus, field, call):
    raw, write = schema_files
    fasta = set_schema_locus(raw, locus)
    raw['expected'][field] = call
    case = corpus.load_real_bcr_cases(write(fasta=fasta))[0]
    assert case.expected[field] == (tuple(call) if isinstance(call, list) else call)


@pytest.fixture
def original_case_file(tmp_path):
    """Copy one committed record to a temporary mutation boundary."""
    from Bio import SeqIO
    cases = json.loads((corpus.REAL_BCR_DIRECTORY / 'cases.json').read_text())
    with (corpus.REAL_BCR_DIRECTORY / 'sequences.fasta').open() as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
    def choose(reason):
        raw, record = next((copy.deepcopy(c), r) for c, r in zip(cases, records)
                           if reason in c['selection_reasons'])
        def write():
            (tmp_path / 'cases.json').write_text(json.dumps([raw]))
            (tmp_path / 'sequences.fasta').write_text('>'+record.id+'\n'+str(record.seq)+'\n')
            return tmp_path
        return raw, write
    return choose


@pytest.mark.parametrize('mutation', [
    'v_start', 'v_end', 'j_start', 'j_end', 'j_before_v', 'trace_query',
    'trace_ref_length', 'v_anchor_query', 'v_anchor_germline', 'j_anchor_query',
    'j_anchor_germline', 'coding_start', 'coding_end', 'coding_translation',
    'coding_scope',
])
def test_loader_rejects_biology_outside_retained_evidence(original_case_file, mutation):
    raw, write = original_case_file('concordant_IGK')
    corpus.load_real_bcr_cases(write())
    expected, alignment = raw['expected'], raw['source']['alignment']
    if mutation in ('v_start', 'v_end', 'j_start', 'j_end'):
        segment, boundary = mutation.split('_')
        expected[segment+'_sequence_'+boundary] += 1
    elif mutation == 'j_before_v':
        expected['j_sequence_start'], expected['j_sequence_end'] = 0, 20
    elif mutation == 'trace_query': alignment['v']['query_aligned'] = 'A'+alignment['v']['query_aligned'][1:]
    elif mutation == 'trace_ref_length': alignment['v']['germline_aligned'] += 'A'
    elif mutation == 'v_anchor_query': alignment['v_imgt104_query_start'] += 3
    elif mutation == 'v_anchor_germline': alignment['v_imgt104_ungapped_offset'] += 3
    elif mutation == 'j_anchor_query': alignment['j_anchor_query_start'] += 3
    elif mutation == 'j_anchor_germline': alignment['j_anchor_germline_offset'] += 3
    elif mutation == 'coding_start': alignment['coding_start'] += 1
    elif mutation == 'coding_end': alignment['coding_end'] -= 3
    elif mutation == 'coding_translation': alignment['coding_translation'] = 'A'+alignment['coding_translation'][1:]
    elif mutation == 'coding_scope': alignment['coding_scope'] = 'through_secondary_j_repeat'
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write())


@pytest.mark.parametrize('mutation', ['boundary', 'bases', 'fabricated_source_indel', 'nonzero_width'])
def test_loader_rejects_deletions_not_supported_by_reference_trace(original_case_file, mutation):
    raw, write = original_case_file('deletion')
    corpus.load_real_bcr_cases(write())
    deletion = raw['expected']['v_deletions'][0]
    if mutation in ('boundary', 'fabricated_source_indel'):
        deletion['query_start'] += 3
        deletion['query_end'] += 3
    elif mutation == 'bases': deletion['sequence'] = 'AAAAAA'
    elif mutation == 'nonzero_width': deletion['query_end'] += 1
    if mutation == 'fabricated_source_indel':
        raw['source']['alignment']['v']['indels'][0].update(deletion)
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write())


@pytest.mark.parametrize('label', [
    'invented_bucket', 'concordant_IGL', 'productivity_disagreement_IGK',
    'tied_call_IGK', 'insertion', 'deletion', 'no_d_IGH', 'pilot_loss',
    'pilot_junction_disagreement',
])
def test_loader_rejects_labels_without_nomination_evidence(original_case_file, label):
    raw, write = original_case_file('concordant_IGK')
    corpus.load_real_bcr_cases(write())
    raw['selection_reasons'] = [label]
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write())


@pytest.mark.parametrize('reason,mutation', [
    ('concordant_IGH', 'different_gene'), ('concordant_IGH', 'different_junction'),
    ('productivity_disagreement_IGH', 'same_productivity'),
    ('tied_call_IGH', 'single_call'), ('insertion', 'missing_reported_indel'),
    ('deletion', 'missing_reported_indel'), ('no_d_IGH', 'reported_d'),
    ('no_d_IGH', 'expected_d'), ('pilot_loss', 'successful_report'),
    ('pilot_junction_disagreement', 'same_junction'),
])
def test_loader_rejects_broken_nomination_evidence(original_case_file, reason, mutation):
    raw, write = original_case_file(reason)
    corpus.load_real_bcr_cases(write())
    evidence = raw['source']['abstar_evidence']
    if mutation == 'different_gene': raw['source']['cellranger']['v_gene'] = 'IGHV2-5'
    elif mutation == 'different_junction': evidence['junction']['abstar_nt'] = 'TGTAAATGG'
    elif mutation == 'same_productivity': evidence['productivity']['abstar'] = True
    elif mutation == 'single_call':
        for segment in ('v', 'j'):
            evidence['raw_calls']['abstar'][segment] = evidence['raw_calls']['abstar'][segment].split(',')[0]
    elif mutation == 'missing_reported_indel': evidence['indels']['v_'+reason+'s'] = None
    elif mutation == 'reported_d': evidence['raw_calls']['abstar']['d'] = 'IGHD3-10*01'
    elif mutation == 'expected_d': raw['expected']['d_call'] = ['IGHD3-10']
    elif mutation == 'successful_report': evidence['status'] = 'annotated'
    elif mutation == 'same_junction': raw['source']['cellranger']['cdr3_nt'] = evidence['junction']
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write())


REQUIRED_BUCKETS = {
    'pilot_loss': 8, 'pilot_junction_disagreement': 2,
    'concordant_IGH': 2, 'concordant_IGK': 2, 'concordant_IGL': 2,
    'productivity_disagreement_IGH': 2, 'productivity_disagreement_IGK': 2,
    'productivity_disagreement_IGL': 2, 'tied_call_IGH': 2, 'tied_call_IGK': 2,
    'tied_call_IGL': 2, 'insertion': 2, 'deletion': 2, 'no_d_IGH': 2,
    'shortest_junction': 1, 'longest_junction': 1,
}


def test_complete_required_bucket_contract(real_bcr_cases):
    from collections import Counter
    assert Counter(reason for case in real_bcr_cases for reason in case.selection_reasons) == REQUIRED_BUCKETS
    assert {(case.dataset, case.sequence_id) for case in real_bcr_cases
            if 'pilot_junction_disagreement' in case.selection_reasons} == {
        ('1279068', 'ATCATCTTCAGCAACT-1_contig_1'),
        ('1279068', 'GTTACAGCACATAACC-1_contig_2'),
    }
    assert {reason: len(case.expected['junction']) for case in real_bcr_cases
            for reason in case.selection_reasons if reason in ('shortest_junction', 'longest_junction')} == {
        'shortest_junction': 18, 'longest_junction': 105,
    }
    corpus.validate_real_bcr_cohort(real_bcr_cases)


@pytest.mark.parametrize('bucket', tuple(REQUIRED_BUCKETS))
def test_cohort_rejects_each_missing_required_bucket(real_bcr_cases, bucket):
    remaining = tuple(case for case in real_bcr_cases if bucket not in case.selection_reasons)
    with pytest.raises(ValueError):
        corpus.validate_real_bcr_cohort(remaining)


@pytest.mark.parametrize('field,value', [
    ('abstar_evidence', None), ('raw_calls', []), ('junction', None), ('productivity', []),
])
def test_loader_rejects_malformed_nomination_containers(original_case_file, field, value):
    raw, write = original_case_file('concordant_IGH')
    if field == 'abstar_evidence': raw['source'][field] = value
    else: raw['source']['abstar_evidence'][field] = value
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write())


@pytest.mark.parametrize('reason,productive,issues', [
    ('productivity_disagreement_IGK', False, ['stop codon(s)']),
    ('pilot_junction_disagreement', False, ['stop codon(s)']),
])
def test_loader_rejects_productivity_reasons_contradicted_by_trace(original_case_file, reason, productive, issues):
    raw, write = original_case_file(reason)
    corpus.load_real_bcr_cases(write())
    raw['expected'].update(productive=productive, productivity_issues=issues)
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write())
