"""FWR3 and CDR3 share the established junction-start codon boundary."""
import json
from pathlib import Path

from Bio.Seq import Seq
import pytest

from abstar.annotation.antibody import Antibody
from abstar.annotation.annotator import annotate_single_sequence

REGIONS = ('fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4')
CASES = json.loads((Path(__file__).parents[1] / 'test_data/fwr3_cdr3_boundaries.json').read_text())['records']


def assert_partition(row):
    nt = row['sequence_oriented'][row['v_sequence_start']:row['j_sequence_end']]
    coding = nt[row['frame'] - 1:]
    protein = str(Seq(coding[:len(coding) // 3 * 3]).translate())
    assert row['fwr3_end'] == row['cdr3_start']
    assert ''.join(row[r] for r in REGIONS) == nt
    assert ''.join(row[r + '_aa'] for r in REGIONS) == protein
    assert row['sequence_alignment'].replace('-', '') == nt
    for mask in ('cdr_mask', 'gene_segment_mask', 'nongermline_mask'):
        assert len(row[mask]) == len(nt), mask
        assert len(row[mask + '_aa']) == len(protein), mask


@pytest.mark.parametrize('case', CASES, ids=lambda c: c['name'])
@pytest.mark.parametrize('reverse', [False, True])
def test_real_fwr3_cdr3_boundary_preserves_junction_and_v_evidence(case, reverse):
    assignment = dict(case['assignment'])
    if reverse:
        assignment['sequence_input'] = str(Seq(assignment['sequence_input']).reverse_complement())
        assignment['rev_comp'] = True
    ab = annotate_single_sequence(Antibody(**assignment), 'human')
    row = ab.to_dict()
    for key, value in (case['expected'] | case['unchanged']).items():
        assert row[key] == value, key
    assert ab.fwr3_end == ab.junction_start + 3
    assert ab.rev_comp is reverse
    assert_partition(row)


@pytest.mark.e2e
@pytest.mark.parametrize('entrypoint', ['api', 'cli'])
def test_shared_boundary_survives_public_outputs(tmp_path, entrypoint):
    import csv
    import abstar
    import polars as pl
    from click.testing import CliRunner
    from abstar.scripts.abstar import cli

    source = tmp_path / 'boundaries.fasta'
    source.write_text(''.join(f">{c['assignment']['sequence_id']}\n{c['assignment']['sequence_input']}\n" for c in CASES))
    project = tmp_path / 'project'
    if entrypoint == 'api':
        abstar.run(str(source), project_path=str(project), output_format=['airr', 'parquet'],
                   n_processes=1, chunksize=1, mmseqs_threads=1, strict=True)
    else:
        result = CliRunner().invoke(cli, ['run', str(source), str(project), '-o', 'airr', '-o', 'parquet',
                                         '--n_processes', '2', '--chunksize', '3', '--mmseqs_threads', '1', '--strict', '--quiet'])
        assert result.exit_code == 0, (result.output, result.exception)
    rows = pl.read_parquet(project / 'parquet/boundaries.parquet').to_dicts()
    with (project / 'airr/boundaries.tsv').open() as handle:
        airr = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == len(airr) == len(CASES)
    for row, text, case in zip(rows, airr, CASES):
        for key, expected in (case['expected'] | case['unchanged'] | case.get('public_unchanged', {})).items():
            assert row[key] == expected, key
        assert row['sequence'] == case['assignment']['sequence_input']
        assert int(text['fwr3_start']) == row['fwr3_start'] + 1
        assert int(text['fwr3_end']) + 1 == int(text['cdr3_start'])
        assert int(text['fwr3_end']) == row['fwr3_end']
        assert text['fwr3'] == row['fwr3']
        assert text['fwr3_aa'] == row['fwr3_aa']
        assert text['junction'] == row['junction']
        assert text['productive'] == ('T' if row['productive'] else 'F')
        assert_partition(row)
    with (project / 'logs/failures.tsv').open() as handle:
        assert list(csv.DictReader(handle, delimiter='\t')) == []


@pytest.mark.e2e
def test_shared_boundary_no_project_python_return():
    import abstar
    import abutils

    case = CASES[0]
    result = abstar.run(abutils.Sequence(case['assignment']['sequence_input'], id=case['assignment']['sequence_id']),
                        n_processes=1, mmseqs_threads=1, strict=True)
    for key, expected in case['expected'].items():
        assert result[key] == expected
    assert_partition(result)


@pytest.mark.parametrize('invalid_boundary', ['start_after_anchor', 'anchor_after_j'])
def test_inconsistent_boundaries_are_rejected_before_slicing(monkeypatch, invalid_boundary):
    import abstar.annotation.annotator as annotator
    from abstar.annotation.regions import RegionSequence

    case = next(c for c in CASES if c['name'] == 'post-anchor_insertion_ownership_-3_True')
    original = annotator.get_region_sequence

    def corrupted_boundary(region, *args, **kwargs):
        result = original(region, *args, **kwargs)
        if region == 'fwr3':
            if invalid_boundary == 'start_after_anchor':
                # Simulate an inconsistent mapped start in the post-anchor insertion.
                return RegionSequence(result.end, result.end, result.sequence[-1:])
            # Simulate a J endpoint that excludes the previously established anchor.
            ab = kwargs['ab']
            ab.j_sequence_end = ab.junction_start + 2
        return result

    monkeypatch.setattr(annotator, 'get_region_sequence', corrupted_boundary)
    with pytest.raises(ValueError, match='FWR3 and junction boundaries do not define a valid interval'):
        annotate_single_sequence(Antibody(**case['assignment']), 'human')
