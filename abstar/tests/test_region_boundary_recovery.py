"""Recover anatomy beyond retained V evidence only when its mapping is unique."""
import json
from pathlib import Path
from Bio.Seq import Seq
import pytest
from abstar.annotation.antibody import Antibody
from abstar.annotation.annotator import annotate_single_sequence

CASES = json.loads((Path(__file__).parents[1] / 'test_data/region_boundary_recovery.json').read_text())['records']
REGIONS = ('fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4')


@pytest.mark.parametrize('case', CASES, ids=lambda c: c['source']['dataset'])
@pytest.mark.parametrize('reverse', [False, True])
def test_recovery_requires_unique_boundaries_and_preserves_assignment(case, reverse):
    assignment = dict(case['assignment'])
    if reverse:
        assignment['sequence_input'] = str(Seq(assignment['sequence_input']).reverse_complement())
        assignment['rev_comp'] = True
    if 'error' in case:
        with pytest.raises(ValueError, match=case['error']):
            annotate_single_sequence(Antibody(**assignment), 'human')
        return
    ab = annotate_single_sequence(Antibody(**assignment), 'human')
    for key, expected in case['expected'].items():
        assert getattr(ab, key) == expected, key
    for key, expected in case['unchanged'].items():
        assert ab.to_dict()[key] == expected, key
    assert 'REGION BOUNDARY RECOVERY' in ab.format_log()
    nt = ab.sequence_oriented[ab.v_sequence_start:ab.j_sequence_end]
    assert ''.join(getattr(ab, r) for r in REGIONS) == nt == ab.sequence
    coding = nt[ab.frame - 1:]
    protein = str(Seq(coding[:len(coding) // 3 * 3]).translate())
    assert ''.join(getattr(ab, r + '_aa') for r in REGIONS) == protein
    assert len(ab.cdr_mask) == len(nt)
    assert len(ab.cdr_mask_aa) == len(protein)


def test_region_tie_is_not_hidden_by_a_unique_anchor():
    from abstar.annotation.junction import _recover_query_positions, recover_junction_anchor
    upstream = 'ACGATCCAGTACGCTGACCATGAGTACCGATGACCTG'
    reference = upstream + 'A' * 8 + 'GTCAGTGT'
    query = upstream + 'A' * 7 + 'GTCAGTGT'
    assert recover_junction_anchor(query, reference, query_origin=0).start == len(query) - 3
    with pytest.raises(ValueError, match='competing optimal boundaries'):
        _recover_query_positions(query, reference, query_origin=0,
                                 projection_positions=(len(upstream) + 3, len(reference) - 3))


@pytest.mark.parametrize('codon', ['TGT', 'TGG', 'TGA', 'TGN'])
def test_extended_regions_do_not_prefer_productive_anchor_codons(codon):
    from types import SimpleNamespace
    from abstar.annotation.junction import recover_v_region_boundaries
    reference = ('ACGT' * 78)[:309] + 'TGT'
    query = reference[:309] + codon
    ab = SimpleNamespace(v_germline_gapped=reference, v_germline_start=0, v_germline_end=150,
                         v_sequence_start=11, v_sequence_end=161, junction_start=None,
                         sequence_oriented='G' * 11 + query, j_sequence_start=323)
    recovered = recover_v_region_boundaries(ab, query[:150], reference[:150])
    assert recovered.intervals == {'fwr2': (125, 176), 'cdr2': (176, 206), 'fwr3': (206, 323)}
    assert (recovered.anchor.start, recovered.anchor.end) == (320, 323)
    assert ab.sequence_oriented[320:323] == codon


def test_retained_mapping_constrains_competing_repeat_placements():
    from abstar.annotation.junction import _recover_query_positions
    upstream = 'ACGATCCAGTACGCTGACCATGAGTACCGATGACCTG'
    reference = upstream + 'A' * 8 + 'GTCAGTGT'
    query = upstream + 'A' * 7 + 'GTCAGTGT'
    boundary = len(upstream) + 3
    _, positions = _recover_query_positions(
        query, reference, query_origin=0, projection_positions=(boundary, len(reference) - 3),
        fixed_positions={boundary: boundary - 1},
    )
    assert positions == (boundary - 1, len(query) - 3)


@pytest.mark.e2e
@pytest.mark.parametrize('entrypoint', ['api', 'cli'])
def test_public_recovery_preserves_ambiguous_failure_and_other_records(tmp_path, entrypoint):
    import csv
    import abstar
    import polars as pl
    from click.testing import CliRunner
    from abstar.scripts.abstar import cli

    source = tmp_path / 'recovery.fasta'
    source.write_text(''.join(f">{c['assignment']['sequence_id']}\n{c['assignment']['sequence_input']}\n" for c in CASES))
    project = tmp_path / 'project'
    if entrypoint == 'api':
        with pytest.warns(RuntimeWarning, match='1 sequence.*failed annotation'):
            abstar.run(str(source), project_path=str(project), output_format=['airr', 'parquet'],
                       n_processes=1, chunksize=1, mmseqs_threads=1)
    else:
        result = CliRunner().invoke(cli, ['run', str(source), str(project), '-o', 'airr', '-o', 'parquet',
                                         '--n_processes', '2', '--chunksize', '2', '--mmseqs_threads', '1', '--quiet'])
        assert result.exit_code == 0, (result.output, result.exception)
    rows = pl.read_parquet(project / 'parquet/recovery.parquet').to_dicts()
    with (project / 'airr/recovery.tsv').open() as handle:
        airr = list(csv.DictReader(handle, delimiter='\t'))
    recovered = [c for c in CASES if 'expected' in c]
    assert len(rows) == len(airr) == len(recovered) == 2
    for row, text, case in zip(rows, airr, recovered):
        for key, value in case['expected'].items():
            if key not in ('junction_start', 'junction_end'):
                assert row[key] == value, key
        for key, value in case['unchanged'].items():
            assert row[key] == value, key
        assert row['sequence'] == case['assignment']['sequence_input']
        assert text['junction'] == row['junction']
        assert text['productive'] == 'T'
        for r in REGIONS:
            assert text[r] == row[r]
            assert text[r + '_aa'] == row[r + '_aa']
        assert int(text['fwr3_start']) == row['fwr3_start'] + 1
        assert int(text['fwr3_end']) == row['fwr3_end']
    with (project / 'logs/failures.tsv').open() as handle:
        failures = list(csv.DictReader(handle, delimiter='\t'))
    assert len(failures) == 1
    assert failures[0]['sequence_id'] == CASES[1]['assignment']['sequence_id']
    diagnostic = (project / 'logs' / failures[0]['diagnostic_path']).read_text()
    assert CASES[1]['error'] in diagnostic
    assert '244, 265, 379, 380, 381' in diagnostic
    assert '247, 268, 379, 380, 381' in diagnostic
    assert len(rows) + len(failures) == len(CASES)


@pytest.mark.e2e
def test_strict_recovery_still_aborts_on_ambiguous_regions(tmp_path):
    import abstar
    source = tmp_path / 'ambiguous.fasta'
    case = CASES[1]
    source.write_text(f">{case['assignment']['sequence_id']}\n{case['assignment']['sequence_input']}\n")
    with pytest.raises(abstar.AnnotationRunError) as caught:
        abstar.run(str(source), project_path=str(tmp_path / 'project'), strict=True,
                   n_processes=1, mmseqs_threads=1)

    assert len(caught.value.failures) == 1
    assert caught.value.failures[0].sequence_id == case['assignment']['sequence_id']
    assert case['error'] in caught.value.failures[0].message


def _recovery_input(reference, query=None, retained_end=150, germline_start=0):
    from types import SimpleNamespace
    query = reference if query is None else query
    return SimpleNamespace(v_germline_gapped=reference, v_germline_start=germline_start,
                           v_germline_end=retained_end, v_sequence_start=11,
                           v_sequence_end=11 + retained_end - germline_start,
                           junction_start=None, sequence_oriented='G' * 11 + query,
                           j_sequence_start=11 + len(query))


@pytest.mark.parametrize('defect,message', [
    ('incomplete_anchor', 'complete IMGT anchor'),
    ('no_upstream_boundary', 'mapped upstream region boundary'),
    ('j_before_upstream', 'invalid upstream/J interval'),
    ('overlapping_reference_boundaries', 'overlapping reference boundaries'),
])
def test_missing_mapping_requires_valid_reference_and_upstream_support(defect, message):
    from abstar.annotation.junction import recover_v_region_boundaries
    reference = ('ACGT' * 78)[:309] + 'TGT'
    ab = _recovery_input(reference)
    query = target = reference[:150]
    if defect == 'incomplete_anchor':
        ab.v_germline_gapped = reference[:-1] + '.'
    elif defect == 'no_upstream_boundary':
        ab.v_germline_start = 125
        query = target = reference[125:150]
    elif defect == 'j_before_upstream':
        ab.j_sequence_start = 100
    else:
        ab.v_germline_gapped = reference[:165] + '.' * 30 + reference[195:]
    with pytest.raises(ValueError, match=message):
        recover_v_region_boundaries(ab, query, target)


def test_deleted_boundary_in_retained_evidence_is_not_fabricated():
    from abstar.annotation.junction import recover_v_region_boundaries
    reference = ('ACGT' * 78)[:309] + 'TGT'
    query = reference[:165] + reference[166:]
    ab = _recovery_input(reference, query, retained_end=170)
    ab.v_sequence_end -= 1
    retained_query = reference[:165] + '-' + reference[166:170]
    with pytest.raises(ValueError, match='deleted boundary or anchor base'):
        recover_v_region_boundaries(ab, retained_query, reference[:170])


def test_extended_mapping_rejects_interrupted_anchor_codon():
    from abstar.annotation.junction import recover_v_region_boundaries
    reference = ('ACGT' * 78)[:309] + 'TGT'
    query = reference[:310] + 'AC' + reference[310:]
    ab = _recovery_input(reference, query)
    with pytest.raises(ValueError, match='interrupted anchor codon'):
        recover_v_region_boundaries(ab, query[:150], reference[:150],
                                    mismatch=-20, gap_open=-1, gap_extend=-1)


@pytest.mark.parametrize('locus', ['TRA', 'TRB', 'TRD', 'TRG'])
def test_broader_recovery_uses_tcr_reference_boundaries(locus):
    from abstar.annotation.germline import get_germline
    from abstar.annotation.junction import recover_v_region_boundaries
    definitions = json.loads((Path(__file__).parent / 'data/tcr/cases.json').read_text())['cases']
    case = next(c for c in definitions if c['locus'] == locus)
    gapped = get_germline(case['source_alleles']['v']['allele'], 'human', receptor='tcr',
                         imgt_gapped=True, exact_match=True, truncate_species=False).sequence
    reference = gapped[:312].replace('.', '')
    # Retain through IMGT position 150, before the FWR3 boundary.
    retained_end = len(gapped[:150].replace('.', ''))
    ab = _recovery_input(reference, retained_end=retained_end)
    ab.v_germline_gapped = gapped
    recovered = recover_v_region_boundaries(ab, reference[:retained_end], reference[:retained_end])
    assert recovered.anchor.start == 11 + len(reference) - 3
    assert recovered.anchor.end == 11 + len(reference)
    assert recovered.intervals['fwr3'] == (11 + len(gapped[:195].replace('.', '')), 11 + len(reference))
    assert ab.sequence_oriented[recovered.anchor.start:recovered.anchor.end] == gapped[309:312]
