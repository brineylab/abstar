"""Compare every public field from one real run writing both final formats."""

import csv
from copy import deepcopy
from importlib.metadata import version

import abstar
import airr
from abutils import Sequence
from Bio.Seq import Seq
import polars as pl
import pytest

from abstar.annotation.schema import OUTPUT_SCHEMA
from abstar.tests.derived import load_derived_bcr_cases
from abstar.tests.helpers import (
    assert_same_annotations,
    expected_airr_amino_acids,
    normalize_airr_row,
    normalize_parquet_row,
)


@pytest.fixture(scope='module', params=(
    'three-locus', 'reverse-complement', 'indel', 'unassigned', 'duplicate-ID',
))
def paired_output(request, public_bcr_cases, tmp_path_factory):
    scenario = request.param
    sequences = [case.as_sequence() for case in public_bcr_cases]
    if scenario == 'reverse-complement':
        sequences = [Sequence(str(Seq(case.sequence).reverse_complement()), id=case.sequence_id)
                     for case in public_bcr_cases]
    elif scenario == 'indel':
        variants = {case.case_id: case for case in load_derived_bcr_cases()}
        sequences = [Sequence(variants[key].sequence, id=key) for key in (
            'IGH-v-insert-1', 'IGH-v-delete-1', 'IGH-v-insert-3', 'IGH-v-delete-3',
        )]
    elif scenario == 'unassigned':
        sequences.insert(1, Sequence('N', id='unassigned'))
    elif scenario == 'duplicate-ID':
        cases = list(public_bcr_cases) + list(public_bcr_cases[:2])
        sequences = [Sequence(case.sequence, id=identifier) for case, identifier in zip(
            cases, ('10E8', '00123', 'same', 'same', '"unterminated'),
        )]
    project = tmp_path_factory.mktemp('parity-' + scenario)
    result = abstar.run(
        sequences, project_path=str(project), output_format=['airr', 'parquet'],
        n_processes=2 if scenario in ('reverse-complement', 'duplicate-ID') else 1,
        chunksize=2, mmseqs_threads=1,
    )
    assert result is None
    # RAW CSV parsing handles literal quotes but performs no type/coordinate
    # conversion. normalize_airr_row owns the one inverse start conversion.
    with (project / 'airr/sequences.tsv').open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        raw_rows = list(reader)
        fields = reader.fieldnames
    parquet = pl.read_parquet(project / 'parquet/sequences.parquet')
    assert len(fields) == len(set(fields)) == len(OUTPUT_SCHEMA) == 168
    assert set(fields) == set(parquet.columns) == set(OUTPUT_SCHEMA)
    assert parquet.schema == pl.Schema(OUTPUT_SCHEMA)
    assert not list((project / 'tmp').rglob('*.parquet'))
    return scenario, sequences, raw_rows, parquet, project


@pytest.mark.e2e
def test_every_shared_public_field_has_identical_meaning(paired_output, public_bcr_cases):
    scenario, inputs, raw_rows, parquet, project = paired_output
    assert len(raw_rows) == parquet.height == len(inputs)
    assert [row['sequence_id'] for row in raw_rows] == parquet['sequence_id'].to_list() == [
        record.id for record in inputs]
    normalized_airr = [normalize_airr_row(row) for row in raw_rows]
    normalized_parquet = [normalize_parquet_row(row) for row in parquet.to_dicts()]
    assert_same_annotations(normalized_airr, normalized_parquet, tuple(OUTPUT_SCHEMA))
    assert normalized_airr == normalized_parquet
    if scenario in ('three-locus', 'reverse-complement'):
        assert len(raw_rows) == 3
        assert parquet['locus'].to_list() == ['IGH', 'IGK', 'IGL']
        for row, case in zip(normalized_parquet, public_bcr_cases):
            assert row['rev_comp'] is (scenario == 'reverse-complement')
            for field in ('junction', 'junction_aa', 'cdr3', 'cdr3_aa', 'productive',
                          'v_sequence_start', 'v_sequence_end', 'j_sequence_start', 'j_sequence_end'):
                assert row[field] == case.expected[field], (case.sequence_id, field)
    elif scenario == 'indel':
        assert parquet['productive'].to_list() == [False, False, True, True]
        assert parquet['v_frameshift'].to_list() == [True, True, False, False]
    elif scenario == 'unassigned':
        row = normalized_parquet[1]
        assert row['annotation_status'] == 'unassigned'
        assert row['failure_reason'] == 'no compatible V gene assignment'
        assert row['productive'] is row['v_call'] is row['junction'] is None
    else:
        assert parquet['locus'].to_list() == ['IGH', 'IGK', 'IGL', 'IGH', 'IGK']
        assert normalized_parquet[2]['sequence'] != normalized_parquet[3]['sequence']
    assert version('airr') == '2.0.0'
    official_reader = airr.read_rearrangement(str(project / 'airr/sequences.tsv'), validate=True)
    try:
        official_rows = list(official_reader)
    finally:
        official_reader.close()
    assert len(official_rows) == len(inputs)
    assert [row['sequence_id'] for row in official_rows] == [record.id for record in inputs]
    # The official reader has ALREADY converted these starts. Never feed these
    # dictionaries to normalize_airr_row or subtract again.
    assert [row['v_sequence_start'] for row in official_rows] == parquet['v_sequence_start'].to_list()


@pytest.mark.e2e
@pytest.mark.parametrize('field', (
    'sequence', 'sequence_aa', 'sequence_alignment_aa', 'germline_alignment_aa',
))
def test_final_parquet_official_sequences_match_independent_source_evidence(paired_output, field):
    _, inputs, _, parquet, _ = paired_output
    for row, original in zip(parquet.to_dicts(), inputs):
        assert row['sequence_input'] == original.sequence
        expected = (original.sequence if field == 'sequence'
                    else expected_airr_amino_acids(row)[field])
        assert row[field] == expected, (original.id, field)


@pytest.mark.e2e
@pytest.mark.parametrize('as_dataframe', (False, True))
@pytest.mark.parametrize('input_count', (1, 3))
def test_no_project_api_keeps_internal_sequence_meanings_and_return_shapes(
    public_bcr_cases, as_dataframe, input_count,
):
    cases = public_bcr_cases[:input_count]
    result = abstar.run([case.as_sequence() for case in cases], as_dataframe=as_dataframe,
                        output_format=['airr', 'parquet'], n_processes=1, mmseqs_threads=1)
    if as_dataframe:
        assert isinstance(result, pl.DataFrame)
        assert result.schema == pl.Schema(OUTPUT_SCHEMA)
        rows = result.to_dicts()
    else:
        if input_count == 1:
            assert isinstance(result, Sequence)
            result = [result]
        assert isinstance(result, list) and all(isinstance(row, Sequence) for row in result)
        rows = [row.annotations for row in result]
    assert len(rows) == input_count
    for row, case in zip(rows, cases):
        assert set(row) == set(OUTPUT_SCHEMA)
        assert row['sequence_id'] == case.sequence_id
        assembled = case.sequence[case.expected['v_sequence_start']:case.expected['j_sequence_end']]
        assert row['sequence'] == assembled != case.sequence
        coding = assembled[row['frame'] - 1:]
        assert row['sequence_aa'] == str(Seq(coding[:len(coding) // 3 * 3]).translate())
        assert len(row['gene_segment_mask']) == len(assembled)
        assert len(row['gene_segment_mask_aa']) == len(row['sequence_aa'])
        assert row['sequence_aa'] != expected_airr_amino_acids(row)['sequence_aa']
        assert row['productive'] == case.expected['productive']


def test_raw_tsv_normalization_changes_only_representation():
    raw = dict.fromkeys(OUTPUT_SCHEMA, '')
    raw.update(sequence_id='00123', sequence='T', sequence_input='A',
               sequence_alignment='A--N', germline_alignment='-ACN',
               sequence_aa='F', sequence_alignment_aa='X-', germline_alignment_aa='-X',
               v_call='IGHV1-2*04,IGHV1-2*02', v_identity='0.875', v_score='137.0',
               v_support='1.23e-45', productive='F', rev_comp='T',
               v_sequence_start='1', v_sequence_end='9', cdr3_start='4', cdr3_end='6',
               annotation_status='unassigned', failure_reason='literal "reason"')
    before = deepcopy(raw)
    expected = dict.fromkeys(OUTPUT_SCHEMA)
    expected.update({field: value for field, value in raw.items() if value != ''})
    expected.update(v_identity=0.875, v_score=137.0, v_support=1.23e-45,
                    productive=False, rev_comp=True, v_sequence_start=0, v_sequence_end=9,
                    cdr3_start=3, cdr3_end=6)
    assert normalize_airr_row(raw) == expected
    assert raw == before
    assert normalize_parquet_row(expected) == expected
    assert normalize_parquet_row({**expected, 'failure_reason': ''})['failure_reason'] is None
    assert {field for field in OUTPUT_SCHEMA if field.endswith('_start')} == {
        *(f'{segment}_{axis}_start' for segment in 'vdjc' for axis in ('sequence', 'germline')),
        *(f'{region}_start' for region in ('fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4')),
    }


@pytest.mark.parametrize('field,value', (
    ('sequence', 'wrong sequence'), ('sequence_aa', 'wrong AA'),
    ('sequence_alignment_aa', 'X-'), ('germline_alignment_aa', '-X'),
    ('v_call', 'IGHV1-2*04,IGHV1-2*02'), ('v_score', 1.0),
    ('v_sequence_start', 1), ('sequence_id', '10E8'),
    ('annotation_status', 'unassigned'), ('failure_reason', 'changed reason'),
))
def test_parquet_normalization_does_not_hide_logical_mismatches(field, value):
    row = dict.fromkeys(OUTPUT_SCHEMA)
    changed = {**row, field: value}
    assert normalize_parquet_row(changed) != normalize_parquet_row(row)
    assert normalize_parquet_row(changed)[field] == value


@pytest.mark.parametrize('normalize,empty', ((normalize_airr_row, ''), (normalize_parquet_row, None)))
def test_normalizers_reject_missing_or_unexpected_public_fields(normalize, empty):
    row = dict.fromkeys(OUTPUT_SCHEMA, empty)
    for changed in ({key: value for key, value in row.items() if key != 'sequence_aa'},
                    {**row, 'new_public_field': empty}, {**row, 'row_id': 'private'}):
        with pytest.raises(AssertionError, match='public field set'):
            normalize(changed)
