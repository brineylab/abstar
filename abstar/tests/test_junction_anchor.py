"""Exact biological anchors and conservative recovery from unmappable endpoints."""
import hashlib
import json
from pathlib import Path

import abutils
import pytest
from Bio.Seq import Seq

from abstar.annotation.antibody import Antibody
from abstar.annotation.annotator import annotate_single_sequence

CASES = json.loads((Path(__file__).parents[1] / 'test_data/lc_anchor_failures.json').read_text())['records']


@pytest.mark.parametrize('case', CASES, ids=lambda c: c['source_file'] + ':' + c['expected']['sequence_id'])
def test_lc_anchor_exact_saved_assignment(case):
    record = dict(case['assignment'])
    record.pop('receptor_type')
    ab = annotate_single_sequence(Antibody(**record), 'human')
    for field, expected in case['internal_expected'].items():
        assert getattr(ab, field) == expected
    for field, expected in case['expected'].items():
        if field != 'sequence':  # final files use original query; object retains assembled VDJ.
            assert getattr(ab, field) == expected, field


def test_reviewed_fixture_anchor_evidence():
    assert len(CASES) == 20
    for case in CASES:
        sequence = case['assignment']['sequence_input']
        expected = case['expected']
        assert hashlib.sha256(sequence.encode()).hexdigest() == case['source']['sequence_sha256']
        start, end = [case['internal_expected'][k] for k in ('junction_start', 'junction_end')]
        assert sequence[start:end] == expected['junction']
        assert sequence.count(expected['junction']) == 1
        assert str(Seq(expected['junction']).translate()) == expected['junction_aa']
        assert case['source']['cellranger']['cdr3_nt'] == expected['junction']


# Unique upstream sequence; reference ends at the homologous anchor codon.
UPSTREAM = 'ACGATCCAGTACGCTGACCATGAGTACCGATGACCTG'


@pytest.mark.parametrize('codon', ['TGT', 'TGC', 'TGG', 'TGA', 'TGN'])
def test_anchor_projection_preserves_mutations_and_ambiguity(codon):
    from abstar.annotation.junction import recover_junction_anchor
    anchor = recover_junction_anchor(
        UPSTREAM + codon + 'AAATGTTTT', UPSTREAM + 'TGT', query_origin=17,
    )
    assert (anchor.start, anchor.end) == (17 + len(UPSTREAM), 20 + len(UPSTREAM))


@pytest.mark.parametrize('prefix', ['', 'AC', 'GCTAGCTAGCTAGCTAGCTAG'])
def test_query_origin_preserves_coordinate_space(prefix):
    from abstar.annotation.junction import recover_junction_anchor
    anchor = recover_junction_anchor(UPSTREAM + 'TGTAAA', UPSTREAM + 'TGT', query_origin=len(prefix))
    assert anchor.start == len(prefix) + len(UPSTREAM)


@pytest.mark.parametrize('reference,query', [
    (UPSTREAM + 'TGT', UPSTREAM + 'TG'),
    (UPSTREAM + 'TGT', 'N' * (len(UPSTREAM) + 3)),
])
def test_missing_anchor_or_unsupported_upstream_fails_explicitly(reference, query):
    from abstar.annotation.junction import recover_junction_anchor, JunctionAnchorError
    with pytest.raises(JunctionAnchorError):
        recover_junction_anchor(query, reference, query_origin=0)


@pytest.mark.parametrize('edit', ['A', 'AT', 'ATC'])
def test_upstream_insertions_shift_anchor_without_imposing_frame(edit):
    from abstar.annotation.junction import recover_junction_anchor
    query = UPSTREAM[:10] + edit + UPSTREAM[10:] + 'TGTAAA'
    anchor = recover_junction_anchor(query, UPSTREAM + 'TGT', query_origin=0)
    assert anchor.start == len(UPSTREAM) + len(edit)


@pytest.mark.parametrize('size', [1, 2, 3])
def test_upstream_deletions_shift_anchor_without_imposing_frame(size):
    from abstar.annotation.junction import recover_junction_anchor
    query = UPSTREAM[:10] + UPSTREAM[10 + size:] + 'TGTAAA'
    anchor = recover_junction_anchor(query, UPSTREAM + 'TGT', query_origin=0)
    assert anchor.start == len(UPSTREAM) - size


def test_competing_optimal_anchor_projections_are_rejected():
    from abstar.annotation.junction import recover_junction_anchor, JunctionAnchorError
    with pytest.raises(JunctionAnchorError, match='competing optimal'):
        recover_junction_anchor(UPSTREAM + 'T' * 6, UPSTREAM + 'T' * 7, query_origin=0)


def test_equivalent_upstream_gap_placements_do_not_make_anchor_ambiguous():
    from abstar.annotation.junction import recover_junction_anchor
    reference = UPSTREAM[:10] + 'A' * 8 + UPSTREAM[10:] + 'TGT'
    query = UPSTREAM[:10] + 'A' * 7 + UPSTREAM[10:] + 'TGTAAA'
    anchor = recover_junction_anchor(query, reference, query_origin=0)
    assert anchor.start == len(UPSTREAM) + 7


@pytest.mark.parametrize('case', CASES, ids=lambda c: c['source_file'] + ':' + c['expected']['sequence_id'])
def test_recovered_anchor_reverse_complement_and_prefix(case):
    record = dict(case['assignment'])
    record.pop('receptor_type')
    prefix = 'GACTA'
    record['sequence_input'] = abutils.tl.reverse_complement(prefix + record['sequence_input'])
    record['rev_comp'] = True
    ab = annotate_single_sequence(Antibody(**record), 'human')
    assert ab.rev_comp is True
    assert ab.junction_start == case['internal_expected']['junction_start'] + len(prefix)
    assert ab.junction_end == case['internal_expected']['junction_end'] + len(prefix)
    assert ab.junction == case['expected']['junction']
    assert ab.junction_aa == case['expected']['junction_aa']
    assert ab.productive is True


def test_recovery_does_not_override_conflicting_retained_anchor():
    from types import SimpleNamespace
    from abstar.annotation.junction import recover_fwr3_anchor, JunctionAnchorError
    prefix = 'A' * 195 + UPSTREAM
    ab = SimpleNamespace(v_germline_gapped='A' * 195 + UPSTREAM + '.' * (114 - len(UPSTREAM)) + 'TGT',
                         v_sequence_start=0, v_germline_start=0,
                         sequence_oriented=prefix + 'TGTTGT', j_sequence_start=len(prefix) + 6)
    with pytest.raises(JunctionAnchorError, match='conflicts with retained'):
        recover_fwr3_anchor(ab, prefix + 'TGTTGT', prefix + '---TGT')


def test_reference_with_incomplete_imgt_anchor_is_rejected():
    from types import SimpleNamespace
    from abstar.annotation.junction import recover_fwr3_anchor, JunctionAnchorError
    ab = SimpleNamespace(v_germline_gapped='A' * 309 + 'TG.')
    with pytest.raises(JunctionAnchorError, match='complete IMGT anchor'):
        recover_fwr3_anchor(ab, '', '')


@pytest.mark.e2e
@pytest.mark.parametrize('workers,chunksize', [(1, 1), (2, 7)])
def test_real_public_api_recovers_all_cases_in_input_order(workers, chunksize):
    import abstar
    records = [abutils.Sequence(c['assignment']['sequence_input'], id='duplicate') for c in CASES]
    result = abstar.run(iter(records), n_processes=workers, chunksize=chunksize,
                        mmseqs_threads=1, strict=True)
    assert isinstance(result, list) and len(result) == 20
    for ab, case in zip(result, CASES):
        assert ab.id == 'duplicate'
        for field, expected in case['expected'].items():
            if field not in ('sequence', 'sequence_id'):
                assert ab[field] == expected, (case['source_file'], field)


@pytest.mark.e2e
def test_real_cli_recovers_cases_in_airr_and_parquet(tmp_path):
    import csv
    import polars as pl
    from click.testing import CliRunner
    from abstar.scripts.abstar import cli
    source = tmp_path / 'anchors.fasta'
    source.write_text(''.join(f">{c['expected']['sequence_id']}\n{c['assignment']['sequence_input']}\n" for c in CASES))
    project = tmp_path / 'project'
    result = CliRunner().invoke(cli, ['run', str(source), str(project), '-o', 'airr', '-o', 'parquet',
                                     '--n_processes', '2', '--chunksize', '3', '--mmseqs_threads', '1',
                                     '--strict', '--quiet'])
    assert result.exit_code == 0, (result.output, result.exception)
    parquet = pl.read_parquet(project / 'parquet/anchors.parquet').to_dicts()
    with (project / 'airr/anchors.tsv').open() as handle:
        airr = list(csv.DictReader(handle, delimiter='\t'))
    assert len(parquet) == len(airr) == 20
    for native, text, case in zip(parquet, airr, CASES):
        for field, expected in case['expected'].items():
            assert native[field] == expected, (case['source_file'], field)
        assert text['junction'] == case['expected']['junction']
        assert text['junction_aa'] == case['expected']['junction_aa']
        assert text['productive'] == 'T' and text['rev_comp'] == 'F'
        for field in ('v_sequence_start', 'j_sequence_start'):
            assert int(text[field]) == native[field] + 1
        for field in ('v_sequence_end', 'j_sequence_end'):
            assert int(text[field]) == native[field]
    with (project / 'logs/failures.tsv').open() as handle:
        assert list(csv.DictReader(handle, delimiter='\t')) == []


@pytest.mark.parametrize('codon,aa,issues', [
    ('TGA', '*AYATDGTLDF', 'stop codon(s)|junction does not start with conserved C'),
    ('TGG', 'WAYATDGTLDF', 'junction does not start with conserved C'),
    ('TGN', 'XAYATDGTLDF', 'ambiguous nucleotide(s)|junction does not start with conserved C|ambiguous nucleotide(s) in junction'),
])
def test_real_fallback_preserves_nonproductive_anchor(codon, aa, issues):
    record = dict(CASES[0]['assignment'])
    record.pop('receptor_type')
    record['sequence_input'] = record['sequence_input'][:365] + codon + record['sequence_input'][368:]
    ab = annotate_single_sequence(Antibody(**record), 'human')
    assert 'JUNCTION ANCHOR RECOVERY' in ab.format_log()
    assert (ab.junction_start, ab.junction_end) == (365, 398)
    assert ab.junction == codon + 'GCATATGCAACTGACGGCACTCTCGACTTC'
    assert ab.junction_aa == aa
    assert ab.productive is False
    assert ab.productivity_issues == issues


def test_deleted_anchor_does_not_select_a_downstream_cysteine():
    from abstar.annotation.junction import recover_junction_anchor
    # A genuine removed anchor leaves AAA followed by a tempting cysteine.
    # Preserve the homologous non-C boundary; do not shift to the later motif.
    anchor = recover_junction_anchor(UPSTREAM + 'AAATGTAAA', UPSTREAM + 'TGT', query_origin=0)
    assert anchor.start == len(UPSTREAM)


@pytest.mark.parametrize('locus', ['TRA', 'TRB', 'TRD', 'TRG'])
def test_fallback_projects_supported_tcr_anchors(monkeypatch, locus):
    from Bio import SeqIO
    root = Path(__file__).parent / 'data/tcr'
    definition = next(c for c in json.loads((root / 'cases.json').read_text())['cases'] if c['locus'] == locus)
    with (root / 'sequences.fasta').open() as handle:
        sequence = next(str(r.seq) for r in SeqIO.parse(handle, 'fasta') if r.id == definition['sequence_id'])
    # Simulate an unmappable upstream endpoint to exercise recovery independently
    # of whether this control happens to trigger the original aligner's edge case.
    monkeypatch.setattr('abstar.annotation.annotator.get_ungapped_position_from_aligned', lambda **kwargs: None)
    calls = {f'{segment}_call': definition['source_alleles'].get(segment, {}).get('allele') for segment in 'vdjc'}
    ab = annotate_single_sequence(Antibody(sequence_id=definition['sequence_id'], sequence_input=sequence, **calls), 'human')
    assert 'JUNCTION ANCHOR RECOVERY' in ab.format_log()
    assert (ab.junction_start, ab.junction_end) == (definition['junction_start'], definition['junction_end'])
    assert ab.junction == definition['junction']
    assert ab.productive is True and ab.productivity_issues == ''


def test_fallback_preserves_junction_frameshift(monkeypatch):
    record = dict(CASES[0]['assignment'])
    record.pop('receptor_type')
    record['sequence_input'] = record['sequence_input'][:380] + 'A' + record['sequence_input'][380:]
    monkeypatch.setattr('abstar.annotation.annotator.get_ungapped_position_from_aligned', lambda **kwargs: None)
    ab = annotate_single_sequence(Antibody(**record), 'human')
    assert 'JUNCTION ANCHOR RECOVERY' in ab.format_log()
    assert (ab.junction_start, ab.junction_end) == (365, 399)
    assert ab.junction == 'TGCGCATATGCAACTAGACGGCACTCTCGACTTC'
    assert ab.productive is False
    assert ab.productivity_issues == 'junction does not end with conserved F|junction length is not a multiple of 3'
