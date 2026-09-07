"""Fixture integrity, loader contracts, and curated real-BCR regressions."""

import copy
from dataclasses import FrozenInstanceError
import hashlib
import json

import pytest

from abstar.tests import corpus
from abstar.tests.derived import load_derived_bcr_cases


def expected_for_selected_reference(case, field, call):
    # The authenticated primary tract is a J4 trace, while J5 is explicitly
    # accepted. Direct J5*02 comparison supports only its exact 32-base suffix:
    # [455:487] / reference [17:49], score 64 with reviewed match=2. Extending
    # to the J4 start 452 would add TAT versus CCC, three mismatches. Keep the
    # source fixture intact and make this separately reviewed allele boundary
    # explicit instead of applying J4 coordinates to J5.
    if (case.sequence_id == 'CCATGTCCAGTCTTCC-1_contig_1'
            and field == 'j_sequence_start' and call == 'IGHJ5*02'):
        return 455
    return case.expected[field]


def annotate_review_probe(sequence, tmp_path):
    import abstar
    import polars as pl
    from abutils import Sequence
    from abstar.annotation.antibody import Antibody
    from abstar.annotation.annotator import annotate_single_sequence

    abstar.run([Sequence(sequence, id='review-probe')], project_path=str(tmp_path),
               output_format='parquet', n_processes=1, mmseqs_threads=1, debug=True)
    rows = pl.read_parquet(tmp_path / 'parquet' / 'sequences.parquet').to_dicts()
    assert len(rows) == 1
    assert rows[0]['sequence_id'] == 'review-probe'
    assignment = pl.read_parquet(tmp_path / 'tmp' / 'chunk_0.parquet').row(0, named=True)
    receptor = assignment.pop('receptor_type')
    ab = Antibody(**assignment)
    ab.receptor_type = receptor
    ab = annotate_single_sequence(ab, germline_database='human')
    return rows[0], ab


@pytest.mark.e2e
def test_compensating_v_indels_use_one_retained_alignment(tmp_path):
    parent = next(c for c in corpus.load_real_bcr_cases()
                  if c.sequence_id == 'GCTGCGAGTCCTGCTT-1_contig_1')
    sequence = parent.sequence[:210] + 'A' + parent.sequence[210:218] + parent.sequence[219:]
    row, ab = annotate_review_probe(sequence, tmp_path)
    assert row['annotation_status'] == 'annotated'
    assert ab.v_insertions == '105:1>A!'
    assert ab.v_deletions == '114:1>T!'
    assert 'out-of-frame indel(s)' not in ab.productivity_issues.split('|')
    query, reference = reconstruct_segment_evidence(ab, 'v')
    assert query == sequence[120:219] + '-' + sequence[219:414]
    source = named_segment_reference(ab, 'v')
    assert reference == source[:90] + '-' + source[90:294]
    assert_emitted_evidence_reconstructs(row, ab)


@pytest.mark.e2e
@pytest.mark.parametrize('reference_start', [0, 1, 2])
def test_full_constant_insertion_preserves_terminal_alignment_columns(tmp_path, reference_start):
    from Bio import SeqIO
    parent = next(c for c in corpus.load_real_bcr_cases()
                  if c.sequence_id == 'GCTGCGAGTCCTGCTT-1_contig_1')
    path = corpus.REAL_BCR_DIRECTORY.parents[2] / 'germline_dbs/bcr/human/ungapped/c.fasta'
    # Read the named source reference independently of annotation's lookup.
    with path.open() as handle:
        constant = next(str(r.seq) for r in SeqIO.parse(handle, 'fasta')
                        if r.id == 'IGHM*01')
    inserted = constant[reference_start:60] + 'AAA' + constant[60:]
    row, ab = annotate_review_probe(parent.sequence[:481] + inserted, tmp_path)
    assert (ab.c_sequence_start, ab.c_sequence_end) == (481, 1840 - reference_start)
    assert ab.c_germline_start == reference_start
    assert ab.c_germline_end == 1356
    assert (len(ab.c_sequence), len(ab.c_germline)) == (1359 - reference_start, 1356 - reference_start)
    assert ab.c_sequence == inserted
    assert ab.c_germline == constant[reference_start:]
    assert row['c_sequence_gapped'].replace('-', '') == inserted
    assert row['c_germline_gapped'].replace('-', '') == constant[reference_start:]
    assert len(row['c_sequence_gapped']) == len(row['c_germline_gapped']) == 1359 - reference_start
    assert row['c_sequence_gapped'][-3:] == constant[-3:]
    assert_emitted_evidence_reconstructs(row, ab)


@pytest.fixture(scope='module')
def constant_evidence_probe(tmp_path_factory):
    from Bio import SeqIO
    parent = next(c for c in corpus.load_real_bcr_cases()
                  if c.sequence_id == 'GCTGCGAGTCCTGCTT-1_contig_1')
    path = corpus.REAL_BCR_DIRECTORY.parents[2] / 'germline_dbs/bcr/human/ungapped/c.fasta'
    with path.open() as handle:
        reference = next(str(r.seq) for r in SeqIO.parse(handle, 'fasta')
                         if r.id == 'IGHM*01')
    # Source-coordinate edits: insert AAA at 60, change codon 30 to a
    # different amino acid, delete codon 50. No result fields define the edits.
    replacement = 'GCT' if reference[90:93] != 'GCT' else 'TGG'
    edited = reference[:60] + 'AAA' + reference[60:90] + replacement + reference[93:150] + reference[153:]
    return annotate_review_probe(parent.sequence[:481] + edited,
                                 tmp_path_factory.mktemp('constant-evidence'))


@pytest.mark.e2e
@pytest.mark.parametrize('field', [
    'c_insertions', 'c_deletions', 'c_mutations', 'c_mutations_aa',
    'c_identity_aa', 'c_sequence_aa', 'c_frame',
    'c_sequence_gapped_aa', 'c_germline_gapped_aa',
    'v_mutations_aa_reference', 'v_mutations_aa_alternate',
    'c_mutations_aa_reference', 'c_mutations_aa_alternate',
    'serialized_v_mutations',
])
def test_evidence_gate_rejects_missing_or_corrupt_evidence(constant_evidence_probe, field):
    row, original = constant_evidence_probe
    row, ab = dict(row), copy.deepcopy(original)
    assert all(getattr(ab, name) for name in (
        'c_insertions', 'c_deletions', 'c_mutations', 'c_mutations_aa', 'v_mutations_aa'))
    assert_emitted_evidence_reconstructs(row, ab)
    if field.startswith('serialized_'):
        row[field.removeprefix('serialized_')] = '1:A>T'
    else:
        key = field.removesuffix('_reference').removesuffix('_alternate')
        if field.endswith(('_reference', '_alternate')):
            events = getattr(ab, key).split('|')
            position, change = events[0].split(':')
            reference, alternate = change.split('>')
            if field.endswith('_reference'):
                reference = 'W' if reference != 'W' else 'A'
            else:
                alternate = 'W' if alternate != 'W' else 'A'
            events[0] = f'{position}:{reference}>{alternate}'
            value = '|'.join(events)
        elif field == 'c_identity_aa':
            value = 0.123456789
        elif field in ('c_sequence_aa', 'c_sequence_gapped_aa', 'c_germline_gapped_aa'):
            current = getattr(ab, field)
            value = ('W' if current[0] != 'W' else 'A') + current[1:]
        elif field == 'c_frame':
            value = ab.c_frame % 3 + 1
        else:
            value = ''
            if field.startswith('c_mutations'):
                count = 'c_mutation_count' + ('_aa' if field.endswith('_aa') else '')
                setattr(ab, count, 0)
                row[count] = 0
        setattr(ab, key, value)
        # Mutate both surfaces so this tests independent reconstruction, not
        # merely the internal/public equality check.
        row[key] = value
    with pytest.raises(AssertionError):
        assert_emitted_evidence_reconstructs(row, ab)


@pytest.mark.parametrize('kind,offset,affected,replacement,expected', [
    ('identity', 0, '', '', 'AACCGG'),
    ('reverse_complement', 0, '', '', 'CCGGTT'),
    ('truncate_5prime', 0, 'AA', '', 'CCGG'),
    ('truncate_3prime', 4, 'GG', '', 'AACC'),
    ('substitute', 2, 'C', 'T', 'AATCGG'),
    ('ambiguity', 2, 'C', 'N', 'AANCGG'),
    ('insert', 3, '', 'GGA', 'AACGGACGG'),
    # The brief's AAG typo is corrected: indices [1:4] are ACC, leaving A + GG.
    ('delete', 1, 'ACC', '', 'AGG'),
    ('insert', 0, '', 'T', 'TAACCGG'),
    ('insert', 6, '', 'T', 'AACCGGT'),
    ('delete', 0, 'A', '', 'ACCGG'),
    ('delete', 5, 'G', '', 'AACCG'),
])
def test_derived_operations_use_exact_input_slices(kind, offset, affected, replacement, expected):
    from abstar.tests.derived import DerivedOperation, derive_sequence
    operation = DerivedOperation(kind, offset, affected, replacement)
    assert derive_sequence('AACCGG', operation) == expected
    with pytest.raises(FrozenInstanceError):
        operation.offset = 1
    assert not hasattr(operation, '__dict__')


@pytest.mark.parametrize('kind,offset,affected,replacement', [
    ('delete', -1, 'G', ''), ('insert', 7, '', 'A'),
    ('delete', 5, 'GG', ''), ('delete', 1, 'CCC', ''),
    ('insert', True, '', 'A'), ('insert', 1.0, '', 'A'),
    ('insert', 1, 'A', 'C'), ('insert', 1, '', ''),
    ('delete', 1, '', ''), ('delete', 1, 'A', 'C'),
    ('delete', 0, 'AACCGG', ''),
    ('substitute', 1, 'A', 'A'), ('substitute', 1, 'AC', 'TT'),
    ('substitute', 1, 'A', 'N'), ('ambiguity', 1, 'A', 'C'),
    ('ambiguity', 1, 'A', 'NN'),
    ('truncate_5prime', 1, 'A', ''), ('truncate_3prime', 1, 'A', ''),
    ('truncate_5prime', 0, '', ''), ('truncate_3prime', 4, 'GG', 'A'),
    ('reverse_complement', 1, '', ''), ('reverse_complement', 0, 'A', ''),
    ('identity', 0, '', 'A'), ('unsupported', 0, '', ''),
    ('insert', 1, '', 'X'), ('insert', 1, '', None),
])
def test_derived_operations_reject_invalid_spans_and_payloads(kind, offset, affected, replacement):
    from abstar.tests.derived import DerivedOperation, derive_sequence
    with pytest.raises(ValueError):
        derive_sequence('AACCGG', DerivedOperation(kind, offset, affected, replacement))


@pytest.mark.parametrize('parent', ['', 'aacg', 'AC-G', 'ACUG', None, 123])
def test_derived_operations_reject_invalid_parent(parent):
    from abstar.tests.derived import DerivedOperation, derive_sequence
    with pytest.raises(ValueError):
        derive_sequence(parent, DerivedOperation('identity', 0, '', ''))


def test_derived_reverse_complement_preserves_iupac_meaning():
    from abstar.tests.derived import DerivedOperation, derive_sequence
    assert derive_sequence('ACGTRYSWKMBDHVN', DerivedOperation('reverse_complement', 0, '', '')) == 'NBDHVKMWSRYACGT'


def test_derived_matrix_covers_all_operations_and_indel_lengths():
    from abstar.tests.derived import load_derived_bcr_cases
    cases = load_derived_bcr_cases()
    assert len(cases) == 60
    assert len({c.case_id for c in cases}) == 60
    for locus in ('IGH', 'IGK', 'IGL'):
        subset = [c for c in cases if c.expected['locus'] == locus]
        assert len(subset) == 20
        assert {c.operation.kind for c in subset} == {
            'identity', 'reverse_complement', 'truncate_5prime', 'truncate_3prime',
            'substitute', 'ambiguity', 'insert', 'delete',
        }
        for context in ('v', 'junction'):
            for kind in ('insert', 'delete'):
                assert {max(len(c.operation.affected_bases), len(c.operation.replacement_bases))
                        for c in subset if c.context == context and c.operation.kind == kind} == {1, 2, 3}
            assert {c.operation.kind for c in subset if c.context == context} >= {'substitute', 'ambiguity'}


def test_derived_loader_checks_hashes_freezes_records_and_preserves_parents():
    from abstar.tests.derived import derive_sequence, load_derived_bcr_cases
    paths = [corpus.REAL_BCR_DIRECTORY / p for p in ('cases.json', 'sequences.fasta')]
    before = [p.read_bytes() for p in paths]
    parents = {(c.dataset, c.sequence_id): c for c in corpus.load_real_bcr_cases()}
    cases = load_derived_bcr_cases()
    for case in cases:
        parent = parents[case.parent['dataset'], case.parent['sequence_id']]
        assert case.sequence == derive_sequence(parent.sequence, case.operation)
        assert hashlib.sha256(case.sequence.encode('ascii')).hexdigest() == case.sequence_sha256
        assert case.parent['sequence_sha256'] == parent.sequence_sha256
        assert not hasattr(case, '__dict__')
        with pytest.raises(FrozenInstanceError):
            case.sequence = 'AAA'
        for mapping in (case.parent, case.expected, case.expected['homologous_junction']):
            with pytest.raises(TypeError):
                mapping['new'] = 'changed'
        assert isinstance(case.evidence, tuple)
    assert load_derived_bcr_cases()[0] is not cases[0]
    assert [p.read_bytes() for p in paths] == before


def test_derived_literal_anchor_projection_and_frame_expectations():
    from abstar.tests.derived import load_derived_bcr_cases
    cases = {c.case_id: c for c in load_derived_bcr_cases()}
    # IGH parent: input length 641, V [120:414], homologous junction [405:450].
    junction = 'TGTGCGAGATATCACCCGGTATTGCGGAATGGTTTTGATGTCTGG'
    forward = cases['IGH-forward'].expected
    reverse = cases['IGH-reverse'].expected
    assert forward['homologous_junction'] == {
        'oriented_start': 405, 'oriented_end': 450, 'input_start': 405,
        'input_end': 450, 'sequence': junction, 'length_mod3': 0,
    }
    assert reverse['rev_comp'] is True
    assert reverse['homologous_junction']['input_start'] == 191
    assert reverse['homologous_junction']['input_end'] == 236
    assert reverse['homologous_junction']['sequence'] == junction
    assert cases['IGH-truncate-5prime'].expected['homologous_junction']['oriented_start'] == 285
    assert cases['IGH-truncate-3prime'].expected['sequence_length'] == 481
    assert cases['IGH-v-insert-1'].expected['coding_frame_delta_mod3'] == 1
    assert cases['IGH-v-delete-1'].expected['coding_frame_delta_mod3'] == 2
    assert cases['IGH-junction-insert-3'].expected['homologous_junction']['sequence'] == 'TGTGCGGGAAGATATCACCCGGTATTGCGGAATGGTTTTGATGTCTGG'
    assert cases['IGH-junction-delete-3'].expected['homologous_junction']['sequence'] == 'TGTGCGTATCACCCGGTATTGCGGAATGGTTTTGATGTCTGG'
    assert cases['IGH-v-ambiguity'].expected['ambiguous_base_count'] == 1


@pytest.fixture
def derived_case_file(tmp_path):
    raw = json.loads((corpus.REAL_BCR_DIRECTORY / 'derived_cases.json').read_text())
    def write():
        path = tmp_path / 'derived_cases.json'
        path.write_text(json.dumps(raw))
        return path
    return raw, write


@pytest.mark.parametrize('mutation', [
    'hash', 'parent_hash', 'unknown_parent', 'numeric_parent', 'offset', 'affected',
    'duplicate', 'unknown_operation', 'extra_operation_field', 'missing_operation_field',
    'extra_case_field', 'sequence_copy', 'no_evidence', 'blank_evidence',
    'expected_length', 'expected_frame', 'expected_junction', 'expected_rev_comp_type',
    'context', 'unclean_parent', 'empty_cases', 'unknown_schema',
])
def test_derived_loader_rejects_corruption(derived_case_file, mutation):
    from abstar.tests.derived import load_derived_bcr_cases
    raw, write = derived_case_file
    case = next(c for c in raw['cases'] if c['case_id'] == 'IGH-v-delete-1')
    if mutation == 'hash':
        case['sequence_sha256'] = '0' * 64
    elif mutation == 'parent_hash':
        case['parent']['sequence_sha256'] = '0' * 64
    elif mutation == 'unknown_parent':
        case['parent']['sequence_id'] = 'unknown'
    elif mutation == 'numeric_parent':
        case['parent']['dataset'] = int(case['parent']['dataset'])
    elif mutation == 'offset':
        case['operation']['offset'] += 1
    elif mutation == 'affected':
        case['operation']['affected_bases'] = 'A' if case['operation']['affected_bases'] != 'A' else 'C'
    elif mutation == 'duplicate':
        raw['cases'].append(copy.deepcopy(case))
    elif mutation == 'unknown_operation':
        case['operation']['kind'] = 'invented'
    elif mutation == 'extra_operation_field':
        case['operation']['unused'] = 0
    elif mutation == 'missing_operation_field':
        del case['operation']['offset']
    elif mutation in ('extra_case_field', 'sequence_copy'):
        case['sequence' if mutation == 'sequence_copy' else 'unused'] = 'ACGT'
    elif mutation in ('no_evidence', 'blank_evidence'):
        case['evidence'] = [] if mutation == 'no_evidence' else [' ']
    elif mutation == 'expected_length':
        case['expected']['sequence_length'] += 1
    elif mutation == 'expected_frame':
        case['expected']['coding_frame_delta_mod3'] = 0
    elif mutation == 'expected_junction':
        case['expected']['homologous_junction']['input_start'] += 1
    elif mutation == 'expected_rev_comp_type':
        case['expected']['rev_comp'] = 0
    elif mutation == 'context':
        case['context'] = 'junction'
    elif mutation == 'unclean_parent':
        parent = corpus.load_real_bcr_cases()[0]
        case['parent'] = dict(dataset=parent.dataset, sequence_id=parent.sequence_id,
                              sequence_sha256=parent.sequence_sha256)
    elif mutation == 'empty_cases':
        raw['cases'] = []
    else:
        raw['schema_version'] = 99
    with pytest.raises(ValueError):
        load_derived_bcr_cases(write())


@pytest.mark.parametrize('location', ['before_v', 'v_anchor', 'j_anchor'])
def test_derived_loader_rejects_rehashed_edits_outside_supported_coding_context(derived_case_file, location):
    """A valid result hash must not turn a flank/anchor edit into an internal coding edit."""
    from abstar.tests.derived import load_derived_bcr_cases
    raw, write = derived_case_file
    case = next(c for c in raw['cases'] if c['case_id'] == 'IGH-v-insert-1')
    parent = next(p for p in corpus.load_real_bcr_cases()
                  if p.sequence_id == case['parent']['sequence_id'])
    offset = {'before_v': 120, 'v_anchor': 405, 'j_anchor': 447}[location]
    case['operation']['offset'] = offset
    if location == 'j_anchor':
        case['context'] = 'junction'
    sequence = parent.sequence[:offset] + 'G' + parent.sequence[offset:]
    case['sequence_sha256'] = hashlib.sha256(sequence.encode('ascii')).hexdigest()
    with pytest.raises(ValueError, match='context'):
        load_derived_bcr_cases(write())


PILOT_IDS = (
    'ACGATACCAGGTTTCA-1_contig_2', 'CAAGATCAGAGCTTCT-1_contig_2',
    'CCATGTCCAGTCTTCC-1_contig_1', 'CTAAGACAGCAATCTC-1_contig_2',
    'CTGTTTACAGGTGCCT-1_contig_1', 'GCTGCTTAGAAACGAG-1_contig_2',
    'GGTGCGTAGGGAAACA-1_contig_1', 'GTAACTGAGTATTGGA-1_contig_2',
)


@pytest.fixture(scope='module')
def original_annotations(tmp_path_factory):
    """Run the public API, retaining real assignment rows for internal checks."""
    import abstar
    import polars as pl
    from abstar.annotation.antibody import Antibody
    from abstar.annotation.annotator import annotate_single_sequence

    cases = corpus.load_real_bcr_cases()
    project = tmp_path_factory.mktemp('original-goldens')
    abstar.run([case.as_sequence() for case in cases], project_path=str(project),
               output_format='parquet', n_processes=1, mmseqs_threads=1, debug=True)
    rows = pl.read_parquet(project / 'parquet' / 'sequences.parquet').to_dicts()
    assert [row['sequence_id'] for row in rows] == [case.sequence_id for case in cases]
    assignments = pl.read_parquet(project / 'tmp' / 'chunk_0.parquet').to_dicts()
    assert [row['sequence_id'] for row in assignments] == [case.sequence_id for case in cases]
    internal = []
    for row, assignment in zip(rows, assignments):
        for segment in 'vjc':
            assert row[f'{segment}_support'] == assignment[f'{segment}_support']
            assert row[f'{segment}_call'] == assignment[f'{segment}_call'].replace('__homo_sapiens', '')
        receptor = assignment.pop('receptor_type')
        ab = Antibody(**assignment)
        ab.receptor_type = receptor
        internal.append(annotate_single_sequence(ab, germline_database='human'))
    return tuple(zip(cases, rows, internal))


@pytest.mark.e2e
@pytest.mark.parametrize('index', range(36))
def test_original_biological_productivity_goldens(original_annotations, index):
    case, row, ab = original_annotations[index]
    assert row['sequence_id'] == ab.sequence_id == case.sequence_id
    assert row['sequence_input'] == row['sequence_oriented'] == case.sequence
    assert row['annotation_status'] == case.expected['status']
    assert row['failure_reason'] is None
    assert 'row_id' not in row
    assert ab.receptor_type == case.source['alignment']['germline_receptor']
    assert row['germline_database'] == ab.germline_database == 'human'
    disagreements = {}
    for field, expected in case.expected.items():
        if field in ('status', 'v_insertions', 'v_deletions'):
            continue
        actual = getattr(ab, field)
        expected = expected_for_selected_reference(case, field, ab.j_call)
        if field.endswith('_call') and isinstance(expected, tuple):
            actual = corpus.normalize_gene(actual)
            if not actual or not set(actual) <= set(expected):
                disagreements[field] = (actual, expected)
        elif field == 'productivity_issues':
            actual = tuple(actual.split('|')) if actual else ()
            if actual != expected:
                disagreements[field] = (actual, expected)
        elif actual != expected:
            disagreements[field] = (actual, expected)
    assert not disagreements, (case.sequence_id, disagreements)
    for segment in 'vj':
        trace = case.source['alignment'][segment]
        if getattr(ab, f'{segment}_call') == trace['reference']:
            expected_score = trace['score']
            if segment == 'j' and case.expected['j_sequence_start'] > trace['query_start']:
                # Two authenticated lambda cases explicitly give the shared
                # V/J base to V. Their raw J trace includes one additional exact
                # match, worth two points in the reviewed boundary score.
                assert case.expected['j_sequence_start'] - trace['query_start'] == 1
                assert trace['query_aligned'][0] == trace['germline_aligned'][0]
                expected_score -= 2
            assert getattr(ab, f'{segment}_score') == expected_score


def named_segment_reference(ab, segment, *, imgt=False):
    """Read a named source FASTA; do not use annotation's lookup or row template."""
    from Bio import SeqIO

    kind = 'imgt_gapped' if imgt else 'ungapped'
    path = (corpus.REAL_BCR_DIRECTORY.parents[2] / 'germline_dbs' /
            ab.receptor_type / ab.germline_database / kind / f'{segment}.fasta')
    call = getattr(ab, f'{segment}_call').split(',')[0].split('__')[0]
    with path.open() as handle:
        matches = [str(r.seq) for r in SeqIO.parse(handle, 'fasta')
                   if r.id.split('__')[0] == call]
    assert len(matches) == 1, (call, path)
    return matches[0]


def reconstruct_segment_evidence(ab, segment):
    """Invert every NT event into one retained pair of alignment rows."""
    from abstar.tests.test_properties import _parse_indels

    # V events use IMGT positions. C events use complete ungapped reference
    # positions, independently of how the retained C pair is serialized.
    template = named_segment_reference(ab, segment, imgt=segment == 'v')
    numbered = [(i, base) for i, base in enumerate(template, 1) if base != '.']
    start, end = (getattr(ab, f'{segment}_germline_{edge}') for edge in ('start', 'end'))
    retained = numbered[start:end]
    mutations = {}
    for mutation in (getattr(ab, f'{segment}_mutations') or '').split('|'):
        if mutation:
            position, change = mutation.split(':')
            reference, alternate = change.split('>')
            position = int(position)
            assert template[position - 1] == reference != alternate
            assert len(reference) == len(alternate) == 1
            assert position not in mutations
            mutations[position] = alternate
    assert len(mutations) == getattr(ab, f'{segment}_mutation_count')
    deleted = set()
    for first, last, bases in _parse_indels(getattr(ab, f'{segment}_deletions')):
        residues = [(i, base) for i, base in retained if first <= i <= last]
        assert ''.join(base for _, base in residues) == bases
        assert (residues[0][0], residues[-1][0]) == (first, last)
        assert deleted.isdisjoint(i for i, _ in residues)
        deleted.update(i for i, _ in residues)
    events = _parse_indels(getattr(ab, f'{segment}_insertions'))
    inserted = {first: bases for first, last, bases in events}
    assert len(inserted) == len(events)
    preceding = numbered[start - 1][0] if start else 0
    query = inserted.pop(preceding, '')
    germline = '-' * len(query)
    for i, base in retained:
        query += '-' if i in deleted else mutations.pop(i, base)
        germline += base
        bases = inserted.pop(i, '')
        query += bases
        germline += '-' * len(bases)
    assert not mutations and not inserted
    assert query.replace('-', '') == getattr(ab, f'{segment}_sequence')
    assert germline.replace('-', '') == getattr(ab, f'{segment}_germline')
    assert query == getattr(ab, f'{segment}_sequence_gapped').replace('.', '')
    if segment == 'c':
        assert germline == ab.c_germline_gapped.replace('.', '')
    assert getattr(ab, f'{segment}_identity') == sum(
        q == g and q != '-' for q, g in zip(query, germline)) / len(query)
    return query, germline


def reconstruct_v_evidence(ab):
    return reconstruct_segment_evidence(ab, 'v')[0].replace('-', '')


def translate_gapped_codons(sequence, frame):
    from Bio.Data import CodonTable

    table = CodonTable.unambiguous_dna_by_id[1]
    codons = {**table.forward_table, **dict.fromkeys(table.stop_codons, '*'),
              '---': '-', '...': '.'}
    return ''.join(codons.get(sequence[i:i + 3], 'X')
                   for i in range(frame - 1, len(sequence) - 2, 3))


def assert_translated_segment_evidence(ab, segment):
    """Reproduce translation and all AA substitutions using Bio's aligner.

    Enumerate optimal global paths under the documented AA score (3/-2,
    affine gaps 35/1); do not assume a unique placement in repeats. Parasail's
    standard 20-AA matrix assigns zero to unknown X/stop characters, so the
    independent Bio matrix explicitly uses that same scoring convention.
    """
    from Bio.Seq import Seq
    from Bio.Align import PairwiseAligner, substitution_matrices

    start = getattr(ab, f'{segment}_germline_start')
    frame = (-start) % 3 + 1
    assert getattr(ab, f'{segment}_frame') == frame
    query_nt = getattr(ab, f'{segment}_sequence')[frame - 1:]
    reference_nt = getattr(ab, f'{segment}_germline')[frame - 1:]
    query = str(Seq(query_nt[:len(query_nt) // 3 * 3]).translate())
    reference = str(Seq(reference_nt[:len(reference_nt) // 3 * 3]).translate())
    assert getattr(ab, f'{segment}_sequence_aa') == query
    assert getattr(ab, f'{segment}_germline_aa') == reference
    if segment == 'c':
        for name in ('c_sequence_gapped', 'c_germline_gapped'):
            assert getattr(ab, name + '_aa') == translate_gapped_codons(getattr(ab, name), frame)
    template = named_segment_reference(ab, segment, imgt=segment == 'v')
    full_aa = ''.join('.' if template[i:i + 3] == '...' else str(Seq(
        template[i:i + 3]).translate()) for i in range(0, len(template) - 2, 3))
    numbered = [(i, residue) for i, residue in enumerate(full_aa, 1) if residue != '.']
    origin = (start + frame - 1) // 3
    assert ''.join(aa for _, aa in numbered[origin:origin + len(reference)]) == reference
    actual = tuple(filter(None, (getattr(ab, f'{segment}_mutations_aa') or '').split('|')))
    assert getattr(ab, f'{segment}_mutation_count_aa') == len(actual)
    alphabet = 'ACDEFGHIKLMNPQRSTVWYX*'
    matrix = substitution_matrices.Array(alphabet, dims=2)
    for a in alphabet[:20]:
        for b in alphabet[:20]:
            matrix[a, b] = 3 if a == b else -2
    aligner = PairwiseAligner(mode='global', substitution_matrix=matrix,
                             open_gap_score=-35, extend_gap_score=-1)
    for alignment in aligner.align(query, reference):
        reference_position = origin
        mutations = []
        for q, g in zip(alignment[0], alignment[1]):
            if g == '-':
                continue
            position, residue = numbered[reference_position]
            assert residue == g
            if q != '-' and q != g:
                mutations.append(f'{position}:{g}>{q}')
            reference_position += 1
        identity = sum(q == g and q != '-' for q, g in zip(
            alignment[0], alignment[1])) / alignment.shape[1]
        if tuple(mutations) == actual and identity == getattr(ab, f'{segment}_identity_aa'):
            break
    else:
        raise AssertionError(f'{ab.sequence_id}: {segment} AA mutations/identity match no optimal source alignment')


def assert_emitted_evidence_reconstructs(row, ab):
    import math
    from abstar.tests.helpers import expected_airr_amino_acids

    assert row['annotation_status'] == 'annotated'
    assert row['failure_reason'] is None
    assert 'row_id' not in row
    official = {'sequence': ab.sequence_input, **expected_airr_amino_acids(vars(ab))}
    for field, value in row.items():
        if hasattr(ab, field):
            expected = official[field] if field in official else getattr(ab, field)
            assert expected == value, (ab.sequence_id, field)
    for suffix in ('', '_aa'):
        # File sequence/AA fields now use official meanings. The published
        # masks still index the same internal assembled V(D)J and translation.
        sequence = getattr(ab, f'sequence{suffix}')
        assert getattr(ab, f'sequence_alignment{suffix}').replace('-', '') == sequence
        assert len(row[f'sequence_alignment{suffix}']) == len(row[f'germline_alignment{suffix}'])
        assert len(row[f'gene_segment_mask{suffix}']) == len(sequence)
        assert len(row[f'nongermline_mask{suffix}']) == len(sequence)
    for segment in 'vdjc':
        call = getattr(ab, f'{segment}_call')
        if call is None:
            assert getattr(ab, f'{segment}_sequence') is None
            assert getattr(ab, f'{segment}_score') is None
            assert getattr(ab, f'{segment}_support') is None
            continue
        start, end = (getattr(ab, f'{segment}_sequence_{edge}') for edge in ('start', 'end'))
        assert 0 <= start < end <= len(ab.sequence_oriented)
        sequence = getattr(ab, f'{segment}_sequence')
        assert ab.sequence_oriented[start:end] == sequence
        reference = named_segment_reference(ab, segment)
        gl_start, gl_end = (getattr(ab, f'{segment}_germline_{edge}') for edge in ('start', 'end'))
        assert 0 <= gl_start < gl_end <= len(reference)
        germline = getattr(ab, f'{segment}_germline')
        assert reference[gl_start:gl_end] == germline
        score = getattr(ab, f'{segment}_score')
        assert math.isfinite(score) and score > 0
        support = getattr(ab, f'{segment}_support')
        if segment in 'vjc':
            assert support is not None
        if support is not None:
            assert math.isfinite(support) and support >= 0
        assert 0 <= getattr(ab, f'{segment}_identity') <= 1
        if segment == 'c':
            query_gapped, germline_gapped = ab.c_sequence_gapped, ab.c_germline_gapped
            assert len(query_gapped) == len(germline_gapped)
            assert len(ab.c_sequence_gapped_aa) == len(ab.c_germline_gapped_aa)
            assert query_gapped.replace('.', '').replace('-', '') == sequence
            assert germline_gapped.replace('.', '').replace('-', '') == germline
            assert ab.c_identity == sum(q == g and q not in '.-' for q, g in zip(
                query_gapped, germline_gapped)) / sum(q != '.' or g != '.' for q, g in zip(
                    query_gapped, germline_gapped))
            reconstruct_segment_evidence(ab, 'c')
            assert_translated_segment_evidence(ab, 'c')
    reconstruct_v_evidence(ab)
    for suffix in ('', '_aa'):
        regions = [getattr(ab, region + suffix)
                   for region in ('fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3')]
        assert all(regions[:-1])
        if not regions[-1]:
            # A truncated local V can omit the complete FWR3 endpoint; these
            # fixtures do not adjudicate partial-region extraction behavior.
            assert ab.v_germline_end < len(ab.v_germline_gapped[:312].replace('.', ''))
        assert getattr(ab, 'v_sequence' + suffix).startswith(''.join(regions)), ab.sequence_id
    assert_translated_segment_evidence(ab, 'v')

@pytest.mark.e2e
def test_partial_v_codon_keeps_full_imgt_amino_acid_positions(original_annotations):
    case, row, ab = next(item for item in original_annotations
                         if item[0].sequence_id == 'TTTGGTTAGGTTACCT-1_contig_2')
    assert (ab.v_sequence_start, ab.v_germline_start, ab.v_frame) == (138, 2, 2)
    assert case.source['alignment']['coding_start'] == 139
    assert ab.v_sequence_aa.startswith('VQLLESGGGL')
    assert '29:T>I' in ab.v_mutations_aa.split('|')
    assert '28:T>I' not in ab.v_mutations_aa.split('|')


@pytest.mark.e2e
@pytest.mark.parametrize('index', range(36))
def test_original_emitted_segment_evidence_reconstructs(original_annotations, index):
    case, row, ab = original_annotations[index]
    assert_emitted_evidence_reconstructs(row, ab)
    # Authenticate explicitly adjudicated indel payloads and query positions
    # by decoding the serialized IMGT events against the retained source trace.
    from abstar.tests.test_properties import _parse_indels
    for field, is_insertion in (('v_insertions', True), ('v_deletions', False)):
        if field not in case.expected:
            continue
        trace = case.source['alignment']['v']
        numbered = [i for i, base in enumerate(ab.v_germline_gapped, 1) if base != '.']
        query_position = trace['query_start']
        germline_position = trace['germline_start']
        query_at_reference = {}
        for q, g in zip(trace['query_aligned'], trace['germline_aligned']):
            query_at_reference.setdefault(germline_position, query_position)
            query_position += q != '-'
            germline_position += g != '-'
        decoded = []
        for start, end, bases in _parse_indels(getattr(ab, field)):
            raw = numbered.index(start) + (1 if is_insertion else 0)
            query_start = query_at_reference[raw]
            decoded.append(dict(query_start=query_start,
                                query_end=query_start + (len(bases) if is_insertion else 0),
                                sequence=bases))
        assert decoded == list(case.expected[field])


@pytest.fixture(scope='module')
def derived_annotations(tmp_path_factory):
    import abstar
    import polars as pl
    from abutils import Sequence
    from abstar.tests.derived import load_derived_bcr_cases
    from abstar.annotation.antibody import Antibody
    from abstar.annotation.annotator import annotate_single_sequence

    cases = load_derived_bcr_cases()
    project = tmp_path_factory.mktemp('derived-goldens')
    abstar.run([Sequence(case.sequence, id=case.case_id) for case in cases],
               project_path=str(project), output_format='parquet',
               n_processes=1, mmseqs_threads=1, debug=True)
    rows = pl.read_parquet(project / 'parquet' / 'sequences.parquet').to_dicts()
    assert [row['sequence_id'] for row in rows] == [case.case_id for case in cases]
    assignments = pl.read_parquet(project / 'tmp' / 'chunk_0.parquet').to_dicts()
    assert [row['sequence_id'] for row in assignments] == [case.case_id for case in cases]
    internal = []
    for assignment in assignments:
        receptor = assignment.pop('receptor_type')
        ab = Antibody(**assignment)
        ab.receptor_type = receptor
        internal.append(annotate_single_sequence(ab, germline_database='human'))
    return tuple(zip(cases, rows, internal))


@pytest.mark.e2e
@pytest.mark.parametrize('case_id', [case.case_id for case in load_derived_bcr_cases()])
def test_derived_named_invariants_and_changed_fields(derived_annotations, case_id):
    from Bio.Seq import Seq

    case, row, ab = next(item for item in derived_annotations if item[0].case_id == case_id)
    parents = {(c.dataset, c.sequence_id): c for c in corpus.load_real_bcr_cases()}
    parent = parents[case.parent['dataset'], case.parent['sequence_id']]
    expected, operation = case.expected, case.operation
    assert row['sequence_id'] == ab.sequence_id == case_id
    assert row['sequence_input'] == case.sequence
    assert len(row['sequence_input']) == expected['sequence_length']
    assert len(case.sequence) - len(parent.sequence) == expected['length_delta']
    assert row['rev_comp'] is expected['rev_comp']
    oriented = str(Seq(case.sequence).reverse_complement()) if expected['rev_comp'] else case.sequence
    assert row['sequence_oriented'] == oriented
    assert row['locus'] == expected['locus']
    assert ab.receptor_type == 'bcr'
    assert row['germline_database'] == 'human'
    assert sum(base not in 'ACGT' for base in oriented) == expected['ambiguous_base_count']
    homologous = expected['homologous_junction']
    assert oriented[homologous['oriented_start']:homologous['oriented_end']] == homologous['sequence']
    original_span = case.sequence[homologous['input_start']:homologous['input_end']]
    if expected['rev_comp']:
        original_span = str(Seq(original_span).reverse_complement())
    assert original_span == homologous['sequence']
    assert len(original_span) % 3 == homologous['length_mod3']
    if operation.kind not in ('identity', 'reverse_complement'):
        a, b = operation.offset, operation.offset + len(operation.affected_bases)
        assert parent.sequence[a:b] == operation.affected_bases
        assert case.sequence == parent.sequence[:a] + operation.replacement_bases + parent.sequence[b:]
    if case.context in ('v', 'junction'):
        assert (len(operation.replacement_bases) - len(operation.affected_bases)) % 3 == expected['coding_frame_delta_mod3']
    else:
        assert expected['coding_frame_delta_mod3'] == 0
    assert_emitted_evidence_reconstructs(row, ab)
    # Pure orientation/flank edits preserve every adjudicated biological value.
    if operation.kind in ('identity', 'reverse_complement', 'truncate_5prime', 'truncate_3prime'):
        for field in ('junction', 'junction_aa', 'cdr3', 'cdr3_aa', 'productive'):
            assert getattr(ab, field) == parent.expected[field], field
        assert not ab.productivity_issues
        for segment in 'vj':
            assert set(corpus.normalize_gene(getattr(ab, f'{segment}_call'))) <= set(parent.expected[f'{segment}_call'])
        assert (ab.junction_start, ab.junction_end) == (homologous['oriented_start'], homologous['oriented_end'])
    if operation.kind == 'ambiguity':
        assert ab.productive is False
        assert 'ambiguous nucleotide(s)' in ab.productivity_issues.split('|')
    if case.context == 'v':
        # Clean derived parents have gap-free authenticated traces. Project
        # the explicit operation onto that source reference, independently of
        # annotator positions/indel helpers and of the selected query endpoints.
        trace = parent.source['alignment']['v']
        references = corpus._packaged_bcr_references()
        template = references['imgt_gapped', 'v'][trace['reference']]
        numbered = [i for i, base in enumerate(template, 1) if base != '.']
        raw = trace['germline_start'] + operation.offset - trace['query_start']
        if operation.kind == 'insert':
            bases = operation.replacement_bases
            marker = '!' if len(bases) % 3 else ''
            assert ab.v_insertions == f'{numbered[raw - 1]}:{len(bases)}>{bases}{marker}'
            assert not ab.v_deletions
        elif operation.kind == 'delete':
            # Repeated bases admit several equivalent gap placements. Enumerate
            # the exact allowed representations from parent/result strings;
            # the annotator must reconstruct the edit, not invent a unique
            # molecular deletion position that the sequence cannot establish.
            width = len(operation.affected_bases)
            representations = set()
            for offset in range(trace['query_start'], parent.expected['junction_start'] - width + 1):
                if parent.sequence[:offset] + parent.sequence[offset + width:] != case.sequence:
                    continue
                first = trace['germline_start'] + offset - trace['query_start']
                start, end = numbered[first], numbered[first + width - 1]
                span = str(start) if start == end else f'{start}-{end}'
                bases = parent.sequence[offset:offset + width]
                marker = '!' if width % 3 else ''
                representations.add(f'{span}:{width}>{bases}{marker}')
            assert ab.v_deletions in representations
            assert not ab.v_insertions
        elif operation.kind in ('substitute', 'ambiguity'):
            reference = template[numbered[raw] - 1]
            expected_mutation = f'{numbered[raw]}:{reference}>{operation.replacement_bases}'
            assert expected_mutation in ab.v_mutations.split('|')


def assert_derived_annotation_consequences(case, row):
    expected = case.expected['annotation']
    assert row['annotation_status'] == 'annotated'
    for field in ('junction', 'junction_aa', 'productive', 'vj_in_frame', 'stop_codon'):
        assert row[field] == expected[field], (case.case_id, field, row[field], expected[field])
    assert (row['v_sequence_start'] + row['frame'] - 1) % 3 == expected['coding_origin_mod3']
    issues = set((row['productivity_issues'] or '').split('|')) - {''}
    conditional = set(expected['alignment_dependent_issues'])
    assert issues - conditional == set(expected['productivity_issues']), (case.case_id, issues, expected)
    assert ('out-of-frame indel(s)' in issues) is row['v_frameshift']
    if case.expected['locus'] in ('IGK', 'IGL'):
        assert row['d_call'] is None


@pytest.mark.e2e
@pytest.mark.parametrize('case_id', [case.case_id for case in load_derived_bcr_cases()])
def test_derived_biological_consequences_match_public_annotation(derived_annotations, case_id):
    case, row, ab = next(item for item in derived_annotations if item[0].case_id == case_id)
    assert_derived_annotation_consequences(case, row)
    homologous = case.expected['homologous_junction']
    assert (ab.junction_start, ab.junction_end) == (homologous['oriented_start'], homologous['oriented_end'])


@pytest.mark.parametrize('field', ['productive', 'vj_in_frame', 'stop_codon', 'junction', 'junction_aa',
                                 'productivity_issues', 'coding_origin_mod3', 'alignment_dependent_issues'])
def test_derived_loader_authenticates_annotation_consequences(derived_case_file, field):
    document, write = derived_case_file
    case = next(case for case in document['cases'] if case['case_id'] == 'IGK-junction-insert-2')
    expected = case['expected']['annotation']
    value = expected[field]
    expected[field] = not value if isinstance(value, bool) else 'corrupted'
    with pytest.raises(ValueError, match='expectations disagree'):
        load_derived_bcr_cases(write())


@pytest.mark.e2e
@pytest.mark.parametrize('field', ['productive', 'vj_in_frame', 'stop_codon', 'junction', 'junction_aa', 'productivity_issues', 'frame'])
def test_derived_consequence_gate_rejects_corrupted_public_output(derived_annotations, field):
    case, original, _ = next(item for item in derived_annotations if item[0].case_id == 'IGK-junction-insert-2')
    row = dict(original)
    assert_derived_annotation_consequences(case, row)
    if field == 'frame':
        row[field] = row[field] % 3 + 1
    elif isinstance(row[field], bool):
        row[field] = not row[field]
    else:
        row[field] = 'corrupted'
    with pytest.raises(AssertionError):
        assert_derived_annotation_consequences(case, row)


@pytest.fixture(autouse=True)
def no_external_commands(monkeypatch, request):
    if request.node.get_closest_marker('e2e') is not None:
        return
    # Applies before the shared case fixtures load, as well as inside test bodies.
    import subprocess
    def forbidden(*args, **kwargs):
        raise AssertionError('fixture integrity must never launch an external command')
    monkeypatch.setattr(subprocess, 'Popen', forbidden)


@pytest.mark.e2e
def test_pilot_loss_cohort_has_no_internal_failures(pilot_loss_cases):
    import abstar

    result = abstar.run(
        [case.as_sequence() for case in pilot_loss_cases],
        n_processes=1, mmseqs_threads=1,
    )
    assert [row.id for row in result] == [case.sequence_id for case in pilot_loss_cases]
    assert all(row['annotation_status'] in {'annotated', 'unassigned'} for row in result)
    for row, case in zip(result, pilot_loss_cases):
        assert_pilot_adjudicated_outcome(row.annotations, case)
        assert 'row_id' not in row.annotations
        for suffix in ('', '_aa'):
            sequence = row[f'sequence{suffix}']
            assert row[f'sequence_alignment{suffix}'].replace('-', '') == sequence
            assert len(row[f'gene_segment_mask{suffix}']) == len(sequence)
            assert len(row[f'nongermline_mask{suffix}']) == len(sequence)


@pytest.mark.e2e
@pytest.mark.parametrize('injected_reason', [
    'unexpected unrelated biological failure',
    'V/J junction is out of frame ',
    'v/j junction is out of frame',
])
def test_pilot_loss_golden_helper_rejects_unrelated_reason(pilot_loss_cases, injected_reason):
    import abstar

    result = abstar.run(
        [case.as_sequence() for case in pilot_loss_cases],
        n_processes=1, mmseqs_threads=1,
    )
    assert [row.id for row in result] == [case.sequence_id for case in pilot_loss_cases]
    for row, case in zip(result, pilot_loss_cases):
        assert_pilot_adjudicated_outcome(row.annotations, case)
        mutated = dict(row.annotations)
        mutated['productivity_issues'] = f"{row['productivity_issues']}|{injected_reason}"
        with pytest.raises(AssertionError, match='productivity_issues'):
            assert_pilot_adjudicated_outcome(mutated, case)


@pytest.mark.parametrize('sequence,germline', [('', 'A'), ('A', ''), ('', '')])
def test_pilot_loss_empty_d_clears_stale_evidence(sequence, germline):
    from abstar.annotation.antibody import Antibody
    from abstar.annotation.annotator import _clear_empty_d_alignment

    ab = Antibody(
        row_id='internal-row', sequence_id='10E8', locus='IGH', rev_comp=True,
        germline_database='human', sequence_oriented='AAACCGGTTT',
        v_sequence_end=3, j_sequence_start=7,
        d_call='IGHD1-14*01', d_gene='IGHD1-14', d_score=12,
        d_support=0.1, d_cigar='4M', d_identity=1.0, d_identity_aa=1.0,
        d_sequence=sequence, d_germline=germline, d_frame=3,
        cdr3_d='C', cdr3_d_aa='P', cdr3_n2='G', cdr3_n2_aa='G',
        np1='C', np1_length=1, np2='G', np2_length=1,
    )
    ab.receptor_type = 'bcr'
    ab.d_sequence_start, ab.d_sequence_end = 4, 6
    ab.d_germline_start, ab.d_germline_end = 7, 9

    result = _clear_empty_d_alignment(ab)

    assert result is ab
    assert (ab.np1, ab.np1_length, ab.np2, ab.np2_length) == ('CCGG', 4, None, None)
    for field in (
        'd_call', 'd_gene', 'd_score', 'd_support', 'd_cigar',
        'd_identity', 'd_identity_aa', 'd_sequence', 'd_germline', 'd_frame',
        'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end',
        'cdr3_d', 'cdr3_d_aa', 'cdr3_n2', 'cdr3_n2_aa',
    ):
        assert getattr(ab, field) is None, field
    assert (ab.row_id, ab.sequence_id, ab.locus, ab.rev_comp,
            ab.receptor_type, ab.germline_database) == (
        'internal-row', '10E8', 'IGH', True, 'bcr', 'human',
    )


@pytest.mark.parametrize('sequence,germline,frame,identities', [
    ('CCGG', 'CCGG', 3, (1.0, None)),
    ('GAAC', 'GAAC', 3, (1.0, None)),
    ('CGTC', 'CGTC', 3, (1.0, None)),
    ('AC', 'AT', 1, (0.5, None)),
    ('AAT', 'A', 1, (1 / 3, None)),
    ('A', 'AAT', 1, (1 / 3, None)),
    ('GCT', 'GCT', 1, (1.0, 1.0)),
])
def test_pilot_loss_empty_d_translation_preserves_nt_identity(sequence, germline, frame, identities):
    from abstar.annotation.annotator import _segment_identities

    assert _segment_identities(sequence, germline, frame) == identities


PILOT_SHORT_D = (
    ('ACGATACCAGGTTTCA-1_contig_2', 'IGHD1-14*01__homo_sapiens', 'CCGG', 441, 445, 7, 11),
    ('CAAGATCAGAGCTTCT-1_contig_2', 'IGHD1-20*01,IGHD1-7*01', 'GAAC', 430, 434, 10, 14),
    ('GGTGCGTAGGGAAACA-1_contig_1', 'IGHD1-14*01__homo_sapiens', 'CCGG', 407, 411, 7, 11),
)


def test_pilot_short_d_old_gt_interval_retains_real_nucleotide_evidence():
    """The old narrower V/J interval supplied a four-base fallback D match."""
    from abstar.annotation.antibody import Antibody
    from abstar.annotation.annotator import ALIGNMENT_PARAMS, _segment_identities
    from abstar.annotation.germline import reassign_dgene, process_dgene_alignment

    case = next(c for c in corpus.load_real_bcr_cases()
                if c.sequence_id == 'GTAACTGAGTATTGGA-1_contig_2')
    local = reassign_dgene(case.sequence[433:442], 'human', 'IGH', 'bcr', ALIGNMENT_PARAMS)
    ab = process_dgene_alignment(case.sequence, 433, local, Antibody())
    assert ab.d_call == 'IGHD6-6*01__homo_sapiens'
    assert (ab.d_sequence_start, ab.d_sequence_end) == (433, 437)
    assert (ab.d_germline_start, ab.d_germline_end) == (13, 17)
    assert ab.d_sequence == ab.d_germline == 'CGTC'
    assert ab.d_score == 12
    assert _segment_identities(ab.d_sequence, ab.d_germline, ab.d_frame) == (1.0, None)


def test_corrected_interval_d_alignment_has_a_unique_stronger_gt_candidate():
    """Score regression for a supplied interval; IGH D biological truth is unadjudicated."""
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
    from abstar.annotation.antibody import Antibody
    from abstar.annotation.annotator import ALIGNMENT_PARAMS, _segment_identities
    from abstar.annotation.germline import reassign_dgene, process_dgene_alignment

    case = next(c for c in corpus.load_real_bcr_cases()
                if c.sequence_id == 'GTAACTGAGTATTGGA-1_contig_2')
    start, end = case.expected['v_sequence_end'], case.expected['j_sequence_start']
    assert (start, end) == (431, 442)
    interval = case.sequence[start:end]
    assert interval == 'ACCGTCACCCG'
    path = corpus.REAL_BCR_DIRECTORY.parents[2] / 'germline_dbs/bcr/human/ungapped/d.fasta'
    with path.open() as handle:
        references = {r.id: str(r.seq) for r in SeqIO.parse(handle, 'fasta')}
    aligner = PairwiseAligner(mode='local', match_score=3, mismatch_score=-2,
                             open_gap_score=-35, extend_gap_score=-1)
    scores = {name: aligner.score(interval, reference) for name, reference in references.items()}
    # Independently score every packaged candidate, including both tied runners-up.
    assert {name: score for name, score in scores.items() if score >= 12} == {
        'IGHD1-14*01__homo_sapiens': 17,
        'IGHD4-4*01__homo_sapiens': 14,
        'IGHD4-17*01__homo_sapiens': 14,
        'IGHD6-6*01__homo_sapiens': 12,
    }
    best = list(aligner.align(interval, references['IGHD1-14*01__homo_sapiens']))
    assert len(best) == 1
    assert best[0].coordinates.tolist() == [[0, 9], [6, 15]]
    assert (best[0][0], best[0][1]) == ('ACCGTCACC', 'ACCGGAACC')
    local = reassign_dgene(interval, 'human', 'IGH', 'bcr', ALIGNMENT_PARAMS)
    ab = process_dgene_alignment(case.sequence, start, local, Antibody())
    assert ab.d_call == 'IGHD1-14*01__homo_sapiens'
    assert (ab.d_sequence_start, ab.d_sequence_end) == (431, 440)
    assert (ab.d_germline_start, ab.d_germline_end) == (6, 15)
    assert (ab.d_sequence, ab.d_germline, ab.d_score) == ('ACCGTCACC', 'ACCGGAACC', 17)
    assert _segment_identities(ab.d_sequence, ab.d_germline, ab.d_frame) == (7 / 9, 2 / 3)


@pytest.mark.e2e
def test_primary_j_boundary_rejects_downstream_j5_repeat(original_annotations):
    import abutils
    from abstar.annotation.germline import get_germline, VJ_BOUNDARY_PARAMS

    case, row, ab = next(item for item in original_annotations
                         if item[0].sequence_id == 'CCATGTCCAGTCTTCC-1_contig_1')
    germline = get_germline('IGHJ5*02', 'human', receptor='bcr', exact_match=True)
    downstream = abutils.tl.local_alignment(case.sequence[412:], germline, **VJ_BOUNDARY_PARAMS)
    assert (412 + downstream.query_begin, 413 + downstream.query_end) == (487, 528)
    assert ab.junction_start == 401
    assert ab.junction_end == 458
    assert 412 + downstream.query_begin >= ab.junction_end
    assert (ab.j_sequence_start, ab.j_sequence_end) == (455, 487)
    assert (ab.j_germline_start, ab.j_germline_end) == (17, 49)
    assert ab.j_sequence == ab.j_germline == 'TGGGGCCAGGGAACCCTGGTCACCGTCTCCTC'
    assert row['j_score'] == 64
    assert corpus.normalize_gene(row['j_call']) == ('IGHJ5',)


def assert_pilot_adjudicated_outcome(row, case):
    expected = case.expected
    assert row['sequence_id'] == case.sequence_id
    assert row['sequence_input'] == row['sequence_oriented'] == case.sequence
    assert row['germline_database'] == 'human'
    assert row['annotation_status'] == expected['status']
    assert row['failure_reason'] is None
    for field in ('locus', 'rev_comp', 'v_sequence_start', 'v_sequence_end',
                  'j_sequence_start', 'j_sequence_end',
                  'junction', 'junction_aa', 'cdr3', 'cdr3_aa'):
        assert row[field] == expected_for_selected_reference(case, field, row['j_call']), (case.sequence_id, field)
    for field in ('v_call', 'j_call'):
        calls = {call.split('*')[0] for call in row[field].split(',')}
        assert calls and calls <= set(expected[field]), (case.sequence_id, field)
    if 'd_call' in expected:
        assert row['d_call'] == expected['d_call']
    issues = row['productivity_issues'].split('|') if row['productivity_issues'] else []
    assert issues == list(expected['productivity_issues']), (case.sequence_id, 'productivity_issues')
    assert row['productive'] is expected['productive']


@pytest.mark.e2e
@pytest.mark.parametrize('sequence_id,call,sequence,start,end,gl_start,gl_end', PILOT_SHORT_D)
def test_pilot_loss_empty_d_translation_retains_nucleotide_evidence(
    pilot_loss_cases, tmp_path, sequence_id, call, sequence, start, end, gl_start, gl_end,
):
    import abstar
    import polars as pl
    from abstar.annotation.antibody import Antibody
    from abstar.annotation.annotator import annotate_single_sequence
    from abstar.annotation.schema import OUTPUT_SCHEMA
    from abstar.tests.helpers import assert_airr_matches_annotations

    case = next(case for case in pilot_loss_cases if case.sequence_id == sequence_id)
    result = abstar.run(
        [case.as_sequence()], project_path=str(tmp_path), n_processes=1,
        mmseqs_threads=1, debug=True, output_format=['airr', 'parquet'],
    )
    assert result is None  # An explicit project path writes public output files.
    rows = pl.read_parquet(tmp_path / 'parquet' / 'sequences.parquet').to_dicts()
    assert_airr_matches_annotations(rows, tmp_path / 'airr' / 'sequences.tsv', tuple(OUTPUT_SCHEMA))
    assert [row['sequence_id'] for row in rows] == [sequence_id]
    row = rows[0]
    assert_pilot_adjudicated_outcome(row, case)
    assert row['annotation_status'] == 'annotated'
    assert row['failure_reason'] is None
    assert row['d_call'] == call
    assert row['d_sequence'] == row['d_germline'] == sequence
    assert row['d_identity'] == 1.0
    assert row['d_identity_aa'] is None
    assert row['d_score'] == 12
    assert row['d_frame'] == 3
    assert case.sequence[start:end] == sequence
    assert 'row_id' not in row

    assert (row['d_sequence_start'], row['d_sequence_end']) == (start, end)
    assert (row['d_germline_start'], row['d_germline_end']) == (gl_start, gl_end)
    # Replay the assignment for the internal junction boundaries too.
    assignment = pl.read_parquet(tmp_path / 'tmp' / 'chunk_0.parquet').row(0, named=True)
    receptor = assignment.pop('receptor_type')
    ab = Antibody(**assignment)
    ab.receptor_type = receptor
    ab = annotate_single_sequence(ab, germline_database='human')
    assert (ab.d_sequence_start, ab.d_sequence_end) == (start, end)
    assert (ab.d_germline_start, ab.d_germline_end) == (gl_start, gl_end)
    assert (ab.junction_start, ab.junction_end) == (
        case.expected['junction_start'], case.expected['junction_end'],
    )
    assert ab.np1 == case.sequence[ab.v_sequence_end:start]
    assert ab.np2 == case.sequence[end:ab.j_sequence_start]


@pytest.mark.e2e
@pytest.mark.parametrize('sequence_id', [
    'CCATGTCCAGTCTTCC-1_contig_1', 'CTAAGACAGCAATCTC-1_contig_2',
    'CTGTTTACAGGTGCCT-1_contig_1', 'GCTGCTTAGAAACGAG-1_contig_2',
])
def test_pilot_loss_mask_preserves_adjudicated_outcome(pilot_loss_cases, tmp_path, sequence_id):
    import abstar
    import polars as pl
    from abstar.annotation.antibody import Antibody
    from abstar.annotation.annotator import annotate_single_sequence

    case = next(case for case in pilot_loss_cases if case.sequence_id == sequence_id)
    abstar.run(
        [case.as_sequence()], project_path=str(tmp_path), output_format='parquet',
        n_processes=1, mmseqs_threads=1, debug=True,
    )
    rows = pl.read_parquet(tmp_path / 'parquet' / 'sequences.parquet').to_dicts()
    assert len(rows) == 1
    row = rows[0]
    assert_pilot_adjudicated_outcome(row, case)
    assert 'row_id' not in row

    assignment = pl.read_parquet(tmp_path / 'tmp' / 'chunk_0.parquet').row(0, named=True)
    receptor = assignment.pop('receptor_type')
    ab = Antibody(**assignment)
    ab.receptor_type = receptor
    ab = annotate_single_sequence(ab, germline_database='human')
    assert_emitted_evidence_reconstructs(row, ab)
    assert ab.receptor_type == 'bcr'
    assert (ab.junction_start, ab.junction_end) == (
        case.expected['junction_start'], case.expected['junction_end'],
    )


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
                'germline_receptor': 'bcr', 'germline_database': 'human',
                'v': {'reference': 'IGKV1-5*01', 'query_start': 0, 'query_end': 3,
                      'germline_start': 261, 'germline_end': 264, 'query_aligned': 'TGT',
                      'germline_aligned': 'TGC', 'indels': [], 'notes': ['direct comparison']},
                'j': {'reference': 'IGKJ1*01', 'query_start': 6, 'query_end': 9,
                      'germline_start': 7, 'germline_end': 10, 'query_aligned': 'TTT',
                      'germline_aligned': 'TTC', 'indels': []},
                'v_imgt104_query_start': 0, 'v_imgt104_ungapped_offset': 261,
                'j_anchor_query_start': 6, 'j_anchor_germline_offset': 7,
                'j_anchor_germline_codon': 'TTC', 'coding_start': 0, 'coding_end': 9,
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
    """Tiny schema traces use literal packaged anchor slices, not gene assignment."""
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
    voffset, vcodon = {'IGH': (285, 'TGT'), 'IGK': (261, 'TGC'), 'IGL': (267, 'TGC')}[locus]
    joffset, jcodon = (18, 'TGG') if locus == 'IGH' else (7, 'TTC')
    source['alignment']['v'].update(reference=v+'*01', germline_start=voffset,
                                   germline_end=voffset+3, germline_aligned=vcodon)
    source['alignment']['j'].update(reference=j+'*01', query_aligned=codon,
                                   germline_aligned=jcodon, germline_start=joffset,
                                   germline_end=joffset+3)
    source['alignment'].update(v_imgt104_ungapped_offset=voffset,
                              j_anchor_germline_offset=joffset,
                              j_anchor_germline_codon=jcodon, coding_translation=aa)
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
    def choose(reason, sequence_id=None):
        raw, record = next((copy.deepcopy(c), r) for c, r in zip(cases, records)
                           if reason in c['selection_reasons']
                           and (sequence_id is None or c['sequence_id'] == sequence_id))
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


def test_loader_rejects_coordinated_shift_to_second_cysteine(original_case_file):
    raw, write = original_case_file('productivity_disagreement_IGL')
    assert (raw['dataset'], raw['sequence_id']) == ('1287179', 'CATATTCTCATACGGT-1_contig_1')
    corpus.load_real_bcr_cases(write())
    alignment, expected = raw['source']['alignment'], raw['expected']
    assert alignment['v']['reference'] == 'IGLV2-23*03'
    assert (alignment['v_imgt104_ungapped_offset'], expected['junction_start']) == (267, 362)
    alignment['v_imgt104_ungapped_offset'] = 270
    alignment['v_imgt104_query_start'] = expected['junction_start'] = 365
    expected['junction'] = expected['junction'][3:]
    expected['junction_aa'] = expected['junction_aa'][1:]
    expected['cdr3'] = expected['junction'][3:-3]
    expected['cdr3_aa'] = expected['junction_aa'][1:-1]
    assert expected['junction_aa'] == 'CSYAGSSTFVVF'
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write())


def test_loader_rejects_coordinated_edit_of_deleted_reference_bases(original_case_file):
    raw, write = original_case_file('deletion')
    corpus.load_real_bcr_cases(write())
    trace = raw['source']['alignment']['v']
    column = trace['query_aligned'].index('------')
    assert trace['germline_aligned'][column:column+6] == 'TAGTGG'
    trace['germline_aligned'] = trace['germline_aligned'][:column]+'AAAAAA'+trace['germline_aligned'][column+6:]
    trace['indels'][0]['sequence'] = 'AAAAAA'
    raw['expected']['v_deletions'][0]['sequence'] = 'AAAAAA'
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write())


@pytest.mark.parametrize('segment', ['v', 'j'])
def test_loader_authenticates_named_reference_trace(original_case_file, segment):
    raw, write = original_case_file('concordant_IGH')
    corpus.load_real_bcr_cases(write())
    trace = raw['source']['alignment'][segment]
    original = trace['germline_aligned'][0]
    replacement = 'A' if original != 'A' else 'C'
    trace['germline_aligned'] = replacement+trace['germline_aligned'][1:]
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write())


@pytest.mark.parametrize('gene', ['IGHG1A', 'IGHG2A', 'IGHG3A'])
@pytest.mark.parametrize('exact', [True, False])
def test_loader_rejects_invented_ighg_a_forms(schema_files, gene, exact):
    raw, write = schema_files
    fasta = set_schema_locus(raw, 'IGH')
    raw['expected']['c_call'] = gene+'*01' if exact else [gene]
    with pytest.raises(ValueError):
        corpus.load_real_bcr_cases(write(fasta=fasta))


def test_loader_accepts_supported_ighg4a_allowed_set(schema_files):
    raw, write = schema_files
    fasta = set_schema_locus(raw, 'IGH')
    raw['expected']['c_call'] = ['IGHG4A']
    assert corpus.load_real_bcr_cases(write(fasta=fasta))[0].expected['c_call'] == ('IGHG4A',)


@pytest.mark.parametrize('segment', ['v', 'j'])
def test_loader_requires_the_named_packaged_allele(original_case_file, segment):
    raw, write = original_case_file('concordant_IGH')
    corpus.load_real_bcr_cases(write())
    trace = raw['source']['alignment'][segment]
    trace['reference'] = trace['reference'].split('*')[0] + '*999'
    with pytest.raises(ValueError, match='named packaged germline'):
        corpus.load_real_bcr_cases(write())


def test_loader_authenticates_secondary_j_trace(original_case_file):
    raw, write = original_case_file('pilot_loss', 'CCATGTCCAGTCTTCC-1_contig_1')
    corpus.load_real_bcr_cases(write())
    trace = raw['source']['alignment']['j_secondary_repeat']
    assert trace['reference'] == 'IGHJ5*02'
    trace['germline_aligned'] = 'A' + trace['germline_aligned'][1:]
    with pytest.raises(ValueError, match='named packaged germline'):
        corpus.load_real_bcr_cases(write())


def test_loader_uses_packaged_resources_without_home_lookup(monkeypatch):
    from pathlib import Path
    def forbidden(*args, **kwargs):
        raise AssertionError('fixture germlines must not use HOME lookup')
    monkeypatch.setattr(Path, 'home', forbidden)
    assert len(corpus.load_real_bcr_cases()) == 36
