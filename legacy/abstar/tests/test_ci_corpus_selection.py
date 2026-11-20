"""Selection contracts: original identity, deterministic diversity, and safe writes."""
import gzip
import importlib.util
import json
from pathlib import Path

import polars as pl
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / 'scripts/build_ci_corpus.py'


@pytest.fixture
def builder():
    assert SCRIPT.exists(), 'The explicit-path corpus builder must exist'
    spec = importlib.util.spec_from_file_location('ci_corpus_builder', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def inputs(tmp_path):
    fastas = tmp_path / 'fasta'
    annotations = tmp_path / 'annotations'
    fastas.mkdir()
    annotations.mkdir()
    manifest = tmp_path / 'samples.csv'
    manifest.write_text('dataset,donor,flow_class\n001,01,naive\n10E8,02,switched\n')
    for dataset, prefix in [('001', 'AA'), ('10E8', 'CC')]:
        rows, fasta = [], []
        for i in range(24):
            sequence = prefix + ''.join('ACGT'[(i // 4 ** k) % 4] for k in range(4)) + 'ACGT' * 30
            sid = 'repeat' if i < 2 else f'{i:04d}'
            fasta.append(f'>{sid}\n{sequence}\n')
            rows.append({'sequence_id': sid, 'sequence': sequence, 'locus': 'IGH' if i % 2 else 'IGK',
                         'productive': i >= 12, 'v_identity': 0.8 if i < 6 else 1.0,
                         'v_germline_start': 0, 'cdr3_length': 12, 'v_call': 'IGHV1-2*01'})
        (fastas / f'{dataset}.fasta').write_text(''.join(fasta))
        pl.DataFrame(rows).write_parquet(annotations / f'{dataset}.parquet')
    return fastas, manifest, annotations


def build(builder, inputs, output, **kwargs):
    fasta, manifest, annotations = inputs
    return builder.build_corpus(fasta_dir=fasta, manifest=manifest, annotations=annotations,
                                output=output, size=20, difficult_fraction=0.2, mandatory=[], **kwargs)


def test_selection_is_reproducible_preserves_strings_and_stratifies(builder, inputs, tmp_path):
    first, second = tmp_path / 'first', tmp_path / 'second'
    build(builder, inputs, first)
    build(builder, inputs, second)
    assert (first / 'sequences.fasta.gz').read_bytes() == (second / 'sequences.fasta.gz').read_bytes()
    assert (first / 'records.parquet').read_bytes() == (second / 'records.parquet').read_bytes()
    rows = pl.read_parquet(first / 'records.parquet')
    assert rows.height == 20
    assert rows['corpus_ordinal'].to_list() == list(range(20))
    assert rows['sequence_sha256'].n_unique() == 20
    assert set(rows['donor']) == {'01', '02'}
    assert set(rows['source_file']) == {'001.fasta', '10E8.fasta'}
    assert rows.filter(pl.col('panel') == 'difficult').height == 4
    assert set(rows.filter(pl.col('panel') == 'background')['locus']) == {'IGH', 'IGK'}
    metadata = json.loads((first / 'manifest.json').read_text())
    assert metadata['expected_failures'] == []
    assert metadata['record_count'] == 20
    assert metadata['files']['sequences.fasta.gz'] == builder.sha256(first / 'sequences.fasta.gz')
    for row in rows.iter_rows(named=True):
        lines = (inputs[0] / row['source_file']).read_text().splitlines()
        assert lines[row['record_ordinal'] * 2] == '>' + row['sequence_id']


def test_duplicate_id_evidence_requires_sequence_disambiguation(builder):
    records = [{'sequence_id': 'same', 'locus': 'IGH'}, {'sequence_id': 'same', 'locus': 'IGK'}]
    with pytest.raises(ValueError, match='ambiguous.*same'):
        builder.annotation_index(pl.DataFrame(records))
    records[0]['sequence'] = 'AAAA'
    records[1]['sequence'] = 'CCCC'
    evidence = builder.annotation_index(pl.DataFrame(records))
    assert builder.match_annotation(evidence, 'same', 'CCCC')['locus'] == 'IGK'
    records[1]['sequence'] = 'AAAA'
    with pytest.raises(ValueError, match='ambiguous.*same'):
        builder.annotation_index(pl.DataFrame(records))


def test_difficulty_uses_explicit_mechanisms(builder):
    reasons = builder.difficulty_reasons({'productive': False, 'v_identity': .79,
        'v_germline_start': 24, 'cdr3_length': 32, 'locus': 'IGH',
        'v_insertions': '12:1>A!', 'productivity_issues': 'out-of-frame indel(s)',
        'v_call': 'IGHV1*01,IGHV1*02'}, 'ACNT')
    assert set(reasons) >= {'nonproductive', 'high_shm', 'v_truncation', 'long_cdr3',
                           'v_insertion', 'frameshift', 'ambiguous_bases', 'v_allele_tie'}
    assert builder.difficulty_reasons(None, 'ACGT') == ['missing_annotation']


def test_mandatory_original_and_cross_panel_exact_dedup(builder, inputs, tmp_path):
    fasta = inputs[0] / '001.fasta'
    lines = fasta.read_text().splitlines()
    sequence = lines[1]
    # Same exact original sequence appears under another external ID.
    fasta.write_text(fasta.read_text() + '>duplicate-sequence\n' + sequence + '\n')
    mandatory = [{'source_file': fasta.name, 'sequence_id': lines[0][1:],
                  'sequence_sha256': builder.sequence_hash(sequence), 'reason': 'reviewed_original'}]
    builder.build_corpus(fasta_dir=inputs[0], manifest=inputs[1], annotations=inputs[2],
                         output=tmp_path / 'selected', size=20, difficult_fraction=.2, mandatory=mandatory)
    rows = pl.read_parquet(tmp_path / 'selected/records.parquet')
    match = rows.filter(pl.col('sequence_sha256') == builder.sequence_hash(sequence))
    assert match.height == 1
    assert match['record_ordinal'].item() == 0
    assert match['panel'].item() == 'difficult'
    assert 'reviewed_original' in match['selection_reasons'].item().to_list()
    assert rows['sequence_sha256'].n_unique() == 20


@pytest.mark.parametrize('location', ['existing', 'inside_fasta', 'inside_annotations'])
def test_output_cannot_replace_or_enter_inputs(builder, inputs, tmp_path, location):
    output = {'existing': tmp_path, 'inside_fasta': inputs[0] / 'out',
              'inside_annotations': inputs[2] / 'out'}[location]
    with pytest.raises(ValueError, match='output'):
        build(builder, inputs, output)


def test_missing_mandatory_sequence_fails_without_partial_output(builder, inputs, tmp_path):
    output = tmp_path / 'selected'
    with pytest.raises(ValueError, match='mandatory'):
        builder.build_corpus(fasta_dir=inputs[0], manifest=inputs[1], annotations=inputs[2],
            output=output, size=20, difficult_fraction=.2,
            mandatory=[{'source_file': '001.fasta', 'sequence_id': 'not-original',
                        'sequence_sha256': '0' * 64, 'reason': 'reviewed'}])
    assert not output.exists()


def test_insufficient_unique_reads_fails_explicitly(builder, inputs, tmp_path):
    with pytest.raises(ValueError, match='insufficient'):
        builder.build_corpus(fasta_dir=inputs[0], manifest=inputs[1], annotations=inputs[2],
            output=tmp_path / 'out', size=100, difficult_fraction=.2, mandatory=[])


def test_duplicate_external_ids_remain_distinct_original_rows(builder, inputs, tmp_path):
    output = tmp_path / 'all'
    builder.build_corpus(fasta_dir=inputs[0], manifest=inputs[1], annotations=inputs[2],
        output=output, size=48, difficult_fraction=.5, mandatory=[])
    rows = pl.read_parquet(output / 'records.parquet')
    duplicates = rows.filter(pl.col('sequence_id') == 'repeat')
    assert duplicates.height == 4
    assert duplicates.select('source_file', 'record_ordinal').unique().height == 4
    assert duplicates['sequence_sha256'].n_unique() == 4
    with gzip.open(output / 'sequences.fasta.gz', 'rt') as handle:
        assert handle.read().count('>repeat\n') == 4


def test_manifest_order_does_not_change_selection(builder, inputs, tmp_path):
    build(builder, inputs, tmp_path / 'first')
    lines = inputs[1].read_text().splitlines()
    inputs[1].write_text('\n'.join([lines[0]] + list(reversed(lines[1:]))) + '\n')
    build(builder, inputs, tmp_path / 'second')
    assert (tmp_path / 'first/records.parquet').read_bytes() == (tmp_path / 'second/records.parquet').read_bytes()


def test_absent_assignment_probe_is_recorded_and_never_synthesized(builder, inputs, tmp_path):
    candidate = {'source_file': '001.fasta', 'sequence_id': 'synthetic-probe',
                 'sequence_sha256': '0' * 64, 'reason': 'original_fixture:probe', 'allow_absent': True}
    result = builder.build_corpus(fasta_dir=inputs[0], manifest=inputs[1], annotations=inputs[2],
        output=tmp_path / 'selected', size=20, difficult_fraction=.2, mandatory=[candidate])
    assert result['selection']['absent_fixture_candidates'] == [candidate]
    assert result['selection']['matched_mandatory_candidates'] == 0
    assert 'synthetic-probe' not in pl.read_parquet(tmp_path / 'selected/records.parquet')['sequence_id']


@pytest.mark.parametrize('sequence_field', ['sequence', 'sequence_input'])
def test_annotation_sequence_mismatch_cannot_label_another_original(builder, sequence_field):
    evidence = builder.annotation_index(pl.DataFrame([
        {'sequence_id': 'same', sequence_field: 'AAAA', 'locus': 'IGH'}]))
    with pytest.raises(ValueError, match='sequence mismatch'):
        builder.match_annotation(evidence, 'same', 'CCCC')


def test_mandatory_fixture_catalog_covers_all_boundary_mechanisms(builder):
    candidates, sources = builder.mandatory_fixtures()
    assert {Path(s['file']).stem for s in sources} >= {
        'region_boundary_recovery', 'lc_anchor_failures', 'missing_fwr3',
        'fwr3_cdr3_boundaries', 'aa_region_boundaries', 'fwr4_endpoints', 'cases'}
    assert {'original_fixture:fwr3_cdr3_boundaries', 'original_fixture:aa_region_boundaries',
            'original_fixture:fwr4_endpoints'} <= {c['reason'] for c in candidates}


def test_mandatory_duplicate_replacement_preserves_all_reviewed_reasons(builder, inputs, tmp_path):
    path = inputs[0] / '001.fasta'
    lines = path.read_text().splitlines()
    later_rank = builder.stable_rank(path.name, 24)
    original = next(i for i in range(24) if builder.stable_rank(path.name, i) > later_rank)
    sequence = lines[original * 2 + 1]
    path.write_text(path.read_text() + '>later-original\n' + sequence + '\n')
    requirements = [
        {'source_file': path.name, 'sequence_id': lines[original * 2][1:],
         'sequence_sha256': builder.sequence_hash(sequence), 'reason': 'first_reviewed_mechanism'},
        {'source_file': path.name, 'sequence_id': 'later-original',
         'sequence_sha256': builder.sequence_hash(sequence), 'reason': 'second_reviewed_mechanism'},
    ]
    builder.build_corpus(fasta_dir=inputs[0], manifest=inputs[1], annotations=inputs[2],
        output=tmp_path / 'selected', size=20, difficult_fraction=.2, mandatory=requirements)
    selected = pl.read_parquet(tmp_path / 'selected/records.parquet').filter(
        pl.col('sequence_sha256') == builder.sequence_hash(sequence))
    assert selected['sequence_id'].item() == 'later-original'
    assert {'first_reviewed_mechanism', 'second_reviewed_mechanism'} <= set(
        selected['selection_reasons'].item())
