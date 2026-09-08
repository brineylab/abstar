"""Continuous-query AA regions must partition the established NT intervals."""
import json
from pathlib import Path

from Bio.Seq import Seq
import pytest

from abstar.annotation.antibody import Antibody
from abstar.annotation.annotator import annotate_single_sequence

REGIONS = ('fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4')
CASES = json.loads((Path(__file__).parents[1] / 'test_data/aa_region_boundaries.json').read_text())['records']


@pytest.mark.parametrize('case', CASES, ids=lambda c: c['name'])
@pytest.mark.parametrize('reverse', [False, True])
def test_real_aa_regions_partition_continuous_query(case, reverse):
    assignment = dict(case['assignment'])
    if reverse:
        assignment['sequence_input'] = str(Seq(assignment['sequence_input']).reverse_complement())
        assignment['rev_comp'] = True
    ab = annotate_single_sequence(Antibody(**assignment), 'human')
    if case['name'] == 'frame3_control':
        assert ab.frame == 3
    for key, value in case['expected'].items():
        assert getattr(ab, key) == value, key
    row = ab.to_dict()
    for key, value in case['unchanged'].items():
        assert row[key] == value, key
    nt = ab.sequence_oriented[ab.v_sequence_start:ab.j_sequence_end]
    coding = nt[ab.frame - 1:]
    protein = str(Seq(coding[:len(coding) // 3 * 3]).translate())
    assert ''.join(getattr(ab, r) for r in REGIONS) == nt
    assert ''.join(getattr(ab, r + '_aa') for r in REGIONS) == protein
    assert len(ab.cdr_mask_aa) == len(ab.gene_segment_mask_aa) == len(ab.nongermline_mask_aa) == len(protein)
    assert ab.cdr3_length == len(ab.cdr3_aa)
    assert ab.rev_comp is reverse


@pytest.mark.e2e
@pytest.mark.parametrize('entrypoint', ['api', 'cli'])
def test_aa_region_partition_public_outputs(tmp_path, entrypoint):
    import abstar
    import csv
    import polars as pl
    from click.testing import CliRunner
    from abstar.scripts.abstar import cli

    source = tmp_path / 'aa.fasta'
    source.write_text(''.join(f">{c['assignment']['sequence_id']}\n{c['assignment']['sequence_input']}\n" for c in CASES))
    project = tmp_path / 'project'
    if entrypoint == 'api':
        abstar.run(str(source), project_path=str(project), output_format=['airr', 'parquet'],
                   n_processes=1, chunksize=1, mmseqs_threads=1, strict=True)
    else:
        result = CliRunner().invoke(cli, ['run', str(source), str(project), '-o', 'airr', '-o', 'parquet',
                                         '--n_processes', '2', '--chunksize', '3', '--mmseqs_threads', '1', '--strict', '--quiet'])
        assert result.exit_code == 0, (result.output, result.exception)
    rows = pl.read_parquet(project / 'parquet/aa.parquet').to_dicts()
    with (project / 'airr/aa.tsv').open() as handle:
        airr = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == len(airr) == len(CASES)
    for row, text, case in zip(rows, airr, CASES):
        assert row['sequence'] == case['assignment']['sequence_input']
        for key, expected in (case['expected'] | case['unchanged'] | case.get('public_unchanged', {})).items():
            assert row[key] == expected, key
        assert text['productive'] == ('T' if row['productive'] else 'F')
        assert text['junction_aa'] == row['junction_aa']
        for r in REGIONS:
            assert text[r + '_aa'] == row[r + '_aa']
        assert int(text['cdr3_length']) == row['cdr3_length'] == len(row['cdr3_aa'])
        assert len(row['cdr_mask_aa']) == sum(len(row[r + '_aa']) for r in REGIONS)
    with (project / 'logs/failures.tsv').open() as handle:
        assert list(csv.DictReader(handle, delimiter='\t')) == []


@pytest.mark.parametrize('removed,query_start,expected_fwr1', [
    (0, 0, 'SSELTQDPAVSVALGQTVRITCQGD'),
    (1, 2, 'SELTQDPAVSVALGQTVRITCQGD'),
    (2, 1, 'SELTQDPAVSVALGQTVRITCQGD'),
])
def test_partial_leading_codon_does_not_enter_aa_regions(removed, query_start, expected_fwr1):
    case = next(c for c in CASES if c['name'] == 'fwr3_2_True')
    assignment = dict(case['assignment'])
    assignment['sequence_input'] = assignment['sequence_input'][case['unchanged']['v_sequence_start'] + removed:]
    ab = annotate_single_sequence(Antibody(**assignment), 'human')
    assert ab.v_sequence_start == query_start
    assert ab.fwr1_aa == expected_fwr1
    assert ab.productive is True
    assert ab.junction_aa == case['unchanged']['junction_aa']
    assert ''.join(getattr(ab, r + '_aa') for r in REGIONS) == ab.sequence_aa
    assert len(ab.cdr_mask_aa) == len(ab.sequence_aa)


@pytest.mark.e2e
def test_no_project_python_return_uses_continuous_regions():
    import abstar
    import abutils

    case = CASES[0]
    result = abstar.run(abutils.Sequence(case['assignment']['sequence_input'], id=case['assignment']['sequence_id']),
                        n_processes=1, mmseqs_threads=1, strict=True)
    for key, value in case['expected'].items():
        assert result[key] == value
    assert result['sequence_aa'] == ''.join(result[r + '_aa'] for r in REGIONS)
    assert len(result['cdr_mask_aa']) == len(result['sequence_aa'])
