"""AIRR 2.0 serialization and biological contracts against the reference reader."""

import importlib
import importlib.util
from copy import deepcopy

import airr
from airr.schema import RearrangementSchema
import polars as pl
import pytest

from abstar.annotation.schema import OUTPUT_SCHEMA


def read_validated(path):
    """The official reader normalizes AIRR starts back to Python coordinates."""
    reader = airr.read_rearrangement(str(path), validate=True)
    try:
        return list(reader)
    finally:
        reader.close()


@pytest.fixture
def serialization():
    class Boundary:
        def __getattr__(self, name):
            assert importlib.util.find_spec('abstar.annotation.airr'), 'named AIRR boundary is missing'
            return getattr(importlib.import_module('abstar.annotation.airr'), name)
    return Boundary()


def test_exact_public_field_order(serialization):
    required = ('sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call',
                'd_call', 'j_call', 'sequence_alignment', 'germline_alignment',
                'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar')
    assert serialization.AIRR_SCHEMA_VERSION == '2.0'
    assert tuple(RearrangementSchema.required) == required
    assert tuple(serialization.AIRR_REQUIRED_FIELDS) == required
    assert tuple(serialization.AIRR_FIELDS) == required + tuple(
        field for field in OUTPUT_SCHEMA if field not in required and field != 'row_id')
    assert len(serialization.AIRR_FIELDS) == len(set(serialization.AIRR_FIELDS))
    assert 'row_id' not in serialization.AIRR_FIELDS


@pytest.mark.parametrize('start,end,expected', [
    (0, 1, (1, 1)), (137, 439, (138, 439)), (None, None, (None, None)),
    (None, 9, (None, None)), (4, None, (None, None)),
])
def test_half_open_to_closed(serialization, start, end, expected):
    assert serialization.to_airr_interval(start, end) == expected


@pytest.mark.parametrize('start,end', [(-1, 1), (0, 0), (10, 2), (1.5, 3), (True, 4)])
def test_invalid_intervals_raise(serialization, start, end):
    with pytest.raises(ValueError, match='half-open interval'):
        serialization.to_airr_interval(start, end)


def test_row_conversion_is_pure_and_converts_every_coordinate(serialization):
    row = {'row_id': 'private', 'sequence_id': '00123', 'sequence': 'assembled',
           'sequence_input': 'AACG', 'sequence_oriented': 'CGTT', 'rev_comp': True,
           'v_call': 'IGHV1-2*02,IGHV1-2*04', 'annotation_status': 'annotated'}
    prefixes = [f'{s}_{a}' for s in 'vdjc' for a in ('sequence', 'germline')]
    prefixes += ['fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4']
    for prefix in prefixes:
        row[prefix + '_start'], row[prefix + '_end'] = 1, 4
    before = deepcopy(row)
    result = serialization.to_airr_row(row)
    assert row == before
    assert result['sequence'] == 'AACG'
    assert result['sequence_id'] == '00123'
    assert result['v_call'] == row['v_call']
    assert result['rev_comp'] is True
    for prefix in prefixes:
        assert (result[prefix + '_start'], result[prefix + '_end']) == (2, 4)
    assert 'row_id' not in result


@pytest.mark.parametrize('query,germline,query_start,germline_start,expected', [
    ('AC-GT', 'ACCGT', 2, 1, '2S1N2M1D2M'),
    ('ACTGT', 'AC-GT', 0, 0, '2M1I2M'),
    ('ACGT', 'TGCA', 0, 0, '4M'),
    ('', '', 0, 0, ''),
])
def test_cigar_run_lengths(serialization, query, germline, query_start, germline_start, expected):
    assert serialization.build_cigar(query, germline, query_start=query_start,
                                     germline_start=germline_start) == expected


@pytest.mark.parametrize('query,germline,query_start,germline_start', [
    ('A', 'AA', 0, 0), ('A-', 'A-', 0, 0), ('A', 'A', -1, 0),
    ('A', 'A', 0, -1), ('A', 'A', 1.5, 0), ('A', 'A', False, 0),
])
def test_invalid_cigar_evidence_raises(serialization, query, germline, query_start, germline_start):
    with pytest.raises(ValueError):
        serialization.build_cigar(query, germline, query_start=query_start,
                                   germline_start=germline_start)


def test_writer_boolean_null_lf_and_reference_validation(serialization, tmp_path):
    frame = pl.DataFrame([
        {'sequence_id': '10E8', 'sequence_input': 'ACG', 'rev_comp': True, 'productive': False},
        {'sequence_id': '00123', 'sequence_input': 'N', 'rev_comp': False,
         'productive': None, 'annotation_status': 'unassigned',
         'failure_reason': 'no compatible V gene assignment'},
    ], schema=OUTPUT_SCHEMA)
    path = tmp_path / 'airr.tsv'
    serialization.write_airr_tsv(frame, path)
    raw = path.read_bytes()
    assert b'\r' not in raw and raw.endswith(b'\n')
    lines = raw.decode().splitlines()
    fields = lines[0].split('\t')
    rows = [dict(zip(fields, line.split('\t'))) for line in lines[1:]]
    assert rows[0]['rev_comp'] == 'T' and rows[0]['productive'] == 'F'
    assert rows[1]['rev_comp'] == 'F' and rows[1]['productive'] == ''
    assert rows[1]['v_call'] == ''
    validated = read_validated(path)
    assert [r['sequence_id'] for r in validated] == ['10E8', '00123']
    assert validated[1]['productive'] is None
    assert validated[1]['annotation_status'] == 'unassigned'


@pytest.mark.parametrize('delimiter', ['\t', '\n', '\r'])
def test_writer_rejects_delimiters(serialization, tmp_path, delimiter):
    frame = pl.DataFrame([{'sequence_id': 'bad' + delimiter + 'id'}], schema=OUTPUT_SCHEMA)
    with pytest.raises(ValueError, match='forbidden delimiter'):
        serialization.write_airr_tsv(frame, tmp_path / 'bad.tsv')


def test_writer_allows_literal_quotes(serialization, tmp_path):
    frame = pl.DataFrame([{'sequence_id': 'a"b', 'sequence_input': 'N'}], schema=OUTPUT_SCHEMA)
    path = tmp_path / 'quoted.tsv'
    serialization.write_airr_tsv(frame, path)
    assert read_validated(path)[0]['sequence_id'] == 'a"b'


def test_coordinate_fields_are_public_model_and_schema_fields():
    from abstar.annotation.antibody import Antibody
    fields = [f'd_{axis}_{end}' for axis in ('sequence', 'germline') for end in ('start', 'end')]
    fields += [f'{region}_{end}' for region in ('fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4')
               for end in ('start', 'end')]
    for field in fields:
        assert field in Antibody.__dataclass_fields__
        assert OUTPUT_SCHEMA[field] == pl.Int64
        assert RearrangementSchema.type(field) == 'integer'


@pytest.fixture(scope='module')
def shared_airr_run(tmp_path_factory):
    """One real public run supplies both formats for all boundary checks."""
    import abstar
    from abutils import Sequence
    from abstar.tests.corpus import load_real_bcr_cases
    from abstar.tests.derived import load_derived_bcr_cases

    originals = load_real_bcr_cases()
    selected = [next(case for case in originals if reason in case.selection_reasons)
                for reason in ('concordant_IGH', 'concordant_IGK', 'concordant_IGL',
                               'insertion', 'deletion')]
    selected.append(next(case for case in originals if not case.expected['productive']))
    variants = [case for case in load_derived_bcr_cases()
                if case.case_id in ('IGH-reverse', 'IGH-truncate-5prime', 'IGH-truncate-3prime')]
    assert len(variants) == 3
    # Opaque and duplicated public IDs must survive both formats in input order.
    ids = ['10E8', '00123', 'same', 'same', 'del', 'nonproductive']
    sequences = [Sequence(case.sequence, id=identifier) for case, identifier in zip(selected, ids)]
    sequences += [Sequence(case.sequence, id=case.case_id) for case in variants]
    sequences.append(Sequence('N', id='unassigned'))
    project = tmp_path_factory.mktemp('shared-airr')
    abstar.run(sequences, project_path=str(project), output_format=['airr', 'parquet'],
               n_processes=1, chunksize=3, mmseqs_threads=1)
    return project, sequences, selected


@pytest.mark.e2e
def test_shared_run_passes_official_validator_and_preserves_outcomes(shared_airr_run):
    from importlib.metadata import version
    project, sequences, selected = shared_airr_run
    assert version('airr') == '2.0.0'
    rows = read_validated(project / 'airr/sequences.tsv')
    assert len(rows) == len(sequences) == 10
    assert [r['sequence_id'] for r in rows] == [s.id for s in sequences]
    for row, sequence in zip(rows, sequences):
        assert row['sequence'] == sequence.sequence
    for row, case in zip(rows, selected):
        assert row['junction'] == case.expected['junction']
        assert row['cdr3'] == case.expected['cdr3']
        assert row['productive'] == case.expected['productive']
        for segment in ('v', 'j'):
            for axis in ('sequence', 'germline'):
                prefix = f'{segment}_{axis}'
                evidence_axis = 'query' if axis == 'sequence' else 'germline'
                trace = case.source['alignment'][segment]
                assert row[prefix + '_start'] == trace[evidence_axis + '_start']
                assert row[prefix + '_end'] == trace[evidence_axis + '_end']
    assert rows[-1]['annotation_status'] == 'unassigned'
    assert rows[-1]['failure_reason'] == 'no compatible V gene assignment'
    assert rows[-1]['productive'] is None and rows[-1]['v_cigar'] == ''
    assert any(r['rev_comp'] for r in rows)
    assert {r['productive'] for r in rows} == {True, False, None}


@pytest.mark.e2e
def test_normalized_tsv_and_parquet_agree_every_public_field(shared_airr_run):
    import csv
    project, _, _ = shared_airr_run
    parquet = pl.read_parquet(project / 'parquet/sequences.parquet')
    with (project / 'airr/sequences.tsv').open(newline='') as handle:
        tsv = list(csv.DictReader(handle, delimiter='\t', quoting=csv.QUOTE_NONE))
    assert len(tsv) == parquet.height
    assert set(tsv[0]) == set(OUTPUT_SCHEMA)
    coordinate_prefixes = {f'{s}_{a}' for s in 'vdjc' for a in ('sequence', 'germline')}
    coordinate_prefixes |= {'fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4'}
    for source, encoded in zip(parquet.iter_rows(named=True), tsv):
        for field, dtype in OUTPUT_SCHEMA.items():
            expected = source['sequence_input'] if field == 'sequence' else source[field]
            if field.endswith('_start') and field[:-6] in coordinate_prefixes and expected is not None:
                expected += 1
            value = encoded[field]
            if expected is None or expected == '':
                assert value == '', (source['sequence_id'], field)
                continue
            if dtype == pl.Boolean:
                assert value in ('T', 'F')
                value = value == 'T'
            elif dtype == pl.Int64:
                value = int(value)
            elif dtype == pl.Float64:
                value = float(value)
            assert value == expected, (source['sequence_id'], field, value, expected)


def replay_cigar(cigar, query, reference):
    """Independently replay serialized operations against original sequences."""
    import re
    operations = re.findall(r'([1-9][0-9]*)([SMIDN])', cigar)
    assert ''.join(n + op for n, op in operations) == cigar
    q, g = 0, 0
    aligned_query, aligned_germline = [], []
    query_start = germline_start = 0
    seen_alignment = False
    for length, operation in operations:
        count = int(length)
        if operation in 'SN':
            assert not seen_alignment
            if operation == 'S': q += count
            else: g += count
            query_start, germline_start = q, g
            continue
        seen_alignment = True
        if operation == 'M':
            aligned_query.append(query[q:q + count])
            aligned_germline.append(reference[g:g + count])
            q += count
            g += count
        elif operation == 'I':
            aligned_query.append(query[q:q + count])
            aligned_germline.append('-' * count)
            q += count
        elif operation == 'D':
            aligned_query.append('-' * count)
            aligned_germline.append(reference[g:g + count])
            g += count
    assert q <= len(query) and g <= len(reference)
    return ''.join(aligned_query), ''.join(aligned_germline), (query_start, q), (germline_start, g)


@pytest.mark.e2e
def test_cigar_replay_and_np_columns_preserve_retained_traces(shared_airr_run):
    from abstar.annotation.germline import get_germline
    project, _, selected = shared_airr_run
    rows = pl.read_parquet(project / 'parquet/sequences.parquet').to_dicts()
    for index, row in enumerate(rows[:-1]):
        query = row['sequence_oriented']
        pairs = {}
        for segment in 'vdjc':
            if row[f'{segment}_sequence'] is None:
                assert row[f'{segment}_cigar'] is None
                continue
            assert row[f'{segment}_cigar'], (row['sequence_id'], segment)
            reference = get_germline(row[f'{segment}_call'].split(',')[0].split('__')[0], 'human',
                                     receptor='bcr', exact_match=True,
                                     force_constant=segment == 'c').sequence
            aq, ag, qs, gs = replay_cigar(row[f'{segment}_cigar'], query, reference)
            assert qs == (row[f'{segment}_sequence_start'], row[f'{segment}_sequence_end'])
            assert gs == (row[f'{segment}_germline_start'], row[f'{segment}_germline_end'])
            assert aq.replace('-', '') == row[f'{segment}_sequence'] == query[slice(*qs)]
            assert ag.replace('-', '') == row[f'{segment}_germline'] == reference[slice(*gs)]
            assert len(aq) == len(ag)
            matches = sum(a == b and a != '-' for a, b in zip(aq, ag))
            assert row[f'{segment}_identity'] == matches / len(aq)
            pairs[segment] = (aq, ag)
            if index < len(selected) and segment in 'vj':
                evidence = selected[index].source['alignment'][segment]
                assert aq == evidence['query_aligned']
                assert ag == evidence['germline_aligned']
        expected_query = pairs['v'][0] + row['np1']
        expected_germline = pairs['v'][1] + '-' * len(row['np1'])
        if 'd' in pairs:
            expected_query += pairs['d'][0] + row['np2']
            expected_germline += pairs['d'][1] + '-' * len(row['np2'])
        expected_query += pairs['j'][0]
        expected_germline += pairs['j'][1]
        assert row['sequence_alignment'] == expected_query
        assert row['germline_alignment'] == expected_germline
        assert len(expected_query) == len(expected_germline)


@pytest.mark.e2e
def test_airr_region_coordinates_slice_the_oriented_query(shared_airr_run):
    project, _, _ = shared_airr_run
    rows = read_validated(project / 'airr/sequences.tsv')
    for row in rows:
        query = row['sequence_oriented']
        for region in ('fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4'):
            start, end = row[region + '_start'], row[region + '_end']
            if not row[region]:
                assert start is end is None
                continue
            assert start is not None and end is not None, (row['sequence_id'], region)
            assert 0 <= start < end <= len(query)
            assert query[start:end] == row[region], (row['sequence_id'], region)
        if row['junction']:
            assert row['cdr3'] == row['junction'][3:-3]
            assert query[row['cdr3_start'] - 3:row['cdr3_end'] + 3] == row['junction']


@pytest.mark.parametrize('aligned,start,end,origin,expected', [
    ('A-CGTT', 3, 5, 137, (139, 142)),
    ('AC-GT', 1, 3, 2, (3, 5)),
    ('AC--GT', 2, 3, 2, (None, None)),
    ('ACGT', None, None, 2, (None, None)),
])
def test_region_alignment_boundary_counts_query_residues(aligned, start, end, origin, expected):
    from abstar.annotation.annotator import region_alignment_to_query_interval
    assert region_alignment_to_query_interval(aligned, start, end, origin) == expected


@pytest.mark.e2e
def test_gap_free_region_coordinates_match_independent_imgt_projection(shared_airr_run):
    from abstar.annotation.germline import get_germline
    project, _, selected = shared_airr_run
    rows = read_validated(project / 'airr/sequences.tsv')
    # These are literal IMGT NT boundaries, inclusive, before removing spacers.
    regions = {'fwr1': (1, 78), 'cdr1': (79, 114), 'fwr2': (115, 165),
               'cdr2': (166, 195), 'fwr3': (196, 312)}
    for case, row in zip(selected[:3], rows[:3]):
        trace = case.source['alignment']['v']
        assert '-' not in trace['query_aligned'] + trace['germline_aligned']
        template = get_germline(trace['reference'], 'human', receptor='bcr',
                                exact_match=True, imgt_gapped=True).sequence
        for region, (first, last) in regions.items():
            reference_start = len(template[:first - 1].replace('.', ''))
            reference_end = len(template[:last].replace('.', ''))
            start = trace['query_start'] + reference_start - trace['germline_start']
            end = trace['query_start'] + reference_end - trace['germline_start']
            assert (row[region + '_start'], row[region + '_end']) == (start, end)
            assert row[region] == case.sequence[start:end]
        assert (row['cdr3_start'], row['cdr3_end']) == (
            case.expected['junction_start'] + 3, case.expected['junction_end'] - 3)
        assert (row['fwr4_start'], row['fwr4_end']) == (
            case.expected['junction_end'] - 3, case.expected['j_sequence_end'])
