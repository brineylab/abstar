#!/usr/bin/env python
"""Build a deterministic original-read CI corpus using saved sampling evidence.

This command never annotates. Existing native Parquets supply difficulty labels
only; a separate reviewed run must create the regression baseline. Source record
ordinals are zero-based and FASTA IDs remain strings, including duplicates.
"""
import argparse
from collections import Counter, defaultdict, deque
import csv
import gzip
import hashlib
import heapq
import json
import math
from pathlib import Path
import sys
import tempfile

from Bio.SeqIO.FastaIO import SimpleFastaParser
import polars as pl

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE_COLUMNS = (
    'sequence_id', 'sequence', 'sequence_input', 'locus', 'productive',
    'productivity_issues', 'v_identity', 'v_germline_start', 'cdr3_length',
    'v_call', 'j_call', 'v_insertions', 'v_deletions', 'c_insertions',
    'c_deletions', 'rev_comp', 'annotation_status',
)
RECORD_SCHEMA = {
    'corpus_ordinal': pl.Int64, 'source_file': pl.String,
    'record_ordinal': pl.Int64, 'sequence_id': pl.String,
    'sequence_sha256': pl.String, 'donor': pl.String, 'flow_class': pl.String,
    'locus': pl.String, 'panel': pl.String, 'selection_reasons': pl.List(pl.String),
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def stable_rank(source_file, ordinal):
    payload = f'bcr-ci-v1\0{source_file}\0{ordinal}'
    return int(hashlib.sha256(payload.encode()).hexdigest(), 16)


def fasta_records(path):
    with Path(path).open(encoding='utf-8') as handle:
        for ordinal, (title, sequence) in enumerate(SimpleFastaParser(handle)):
            if not title or not sequence:
                raise ValueError(f'empty FASTA identifier/sequence at {path}:{ordinal}')
            yield ordinal, title.split()[0], sequence


def annotation_index(frame):
    """Index sampling labels without silently selecting a duplicate ID's row.

    Native final Parquet sequence is the original query. If duplicate IDs have
    differing labels, exact original sequence disambiguation is mandatory. Old
    Parquets whose sequence has legacy assembled semantics cannot resolve these
    inputs and are rejected at the matching boundary.
    """
    if frame.schema.get('sequence_id') != pl.String:
        raise ValueError('annotation sequence_id must have string dtype')
    result = {}
    for row in frame.iter_rows(named=True):
        sid = row['sequence_id']
        if not isinstance(sid, str) or not sid:
            raise ValueError('annotation sequence_id must be a nonempty string')
        sequence = row.get('sequence_input') or row.get('sequence')
        digest = sequence_hash(sequence) if sequence else None
        evidence = {k: v for k, v in row.items() if k not in ('sequence', 'sequence_input')}
        bucket = result.setdefault(sid, {})
        if digest in bucket and bucket[digest] != evidence:
            raise ValueError(f'ambiguous annotation evidence for ID {sid!r}')
        bucket[digest] = evidence
        if None in bucket and len(bucket) > 1:
            raise ValueError(f'ambiguous annotation evidence for ID {sid!r}')
    return result


def match_annotation(index, sid, sequence):
    bucket = index.get(sid)
    if not bucket:
        return None
    exact = bucket.get(sequence_hash(sequence))
    if exact is not None:
        return exact
    if set(bucket) == {None}:
        return bucket[None]
    if len(bucket) == 1:
        raise ValueError(f'annotation sequence mismatch for ID {sid!r}')
    raise ValueError(f'ambiguous annotation evidence for ID {sid!r}: no exact original query match')


def difficulty_reasons(row, sequence):
    reasons = []
    if any(base not in 'ACGTacgt' for base in sequence):
        reasons.append('ambiguous_bases')
    if row is None:
        return sorted(reasons + ['missing_annotation'])
    if row.get('productive') is False:
        reasons.append('nonproductive')
    issues = row.get('productivity_issues') or ''
    if 'frame' in issues:
        reasons.append('frameshift')
    if 'stop' in issues:
        reasons.append('stop_codon')
    if row.get('v_identity') is not None and row['v_identity'] < .9:
        reasons.append('high_shm')
    if (row.get('v_germline_start') or 0) >= 15:
        reasons.append('v_truncation')
    length = row.get('cdr3_length')
    if length is not None:
        if length <= 5:
            reasons.append('short_cdr3')
        if length >= (25 if row.get('locus') == 'IGH' else 13):
            reasons.append('long_cdr3')
    for segment in ('v', 'c'):
        for event in ('insertion', 'deletion'):
            if row.get(f'{segment}_{event}s'):
                reasons.append(f'{segment}_{event}')
    for segment in ('v', 'j'):
        if ',' in (row.get(f'{segment}_call') or ''):
            reasons.append(f'{segment}_allele_tie')
    if row.get('rev_comp'):
        reasons.append('reverse_complement')
    if row.get('annotation_status') == 'unassigned':
        reasons.append('unassigned')
    return sorted(reasons)


def mandatory_fixtures():
    """Return candidate fixture originals; source matching authenticates them.

    Derived cases are excluded. Some assignment fixtures are edited probes:
    unmatched candidates are explicitly recorded as absent, never synthesized.
    """
    candidates, sources = [], []
    paths = [ROOT / 'abstar/test_data' / (name + '.json') for name in (
        'region_boundary_recovery', 'lc_anchor_failures', 'missing_fwr3',
        'fwr3_cdr3_boundaries', 'aa_region_boundaries', 'fwr4_endpoints',
    )]
    for path in paths:
        sources.append({'file': str(path.relative_to(ROOT)), 'sha256': sha256(path)})
        for case in json.loads(path.read_text())['records']:
            assignment = case['assignment']
            filename = case.get('source_file') or case['source']['dataset'] + '.fasta'
            candidates.append({'source_file': filename, 'sequence_id': assignment['sequence_id'],
                'sequence_sha256': sequence_hash(assignment['sequence_input']),
                'reason': 'original_fixture:' + path.stem, 'allow_absent': True})
    path = ROOT / 'abstar/tests/data/real_bcr/cases.json'
    sources.append({'file': str(path.relative_to(ROOT)), 'sha256': sha256(path)})
    for case in json.loads(path.read_text()):
        candidates.append({'source_file': case['dataset'] + '.fasta',
            'sequence_id': case['sequence_id'], 'sequence_sha256': case['sequence_sha256'],
            'reason': 'original_fixture:real_bcr', 'allow_absent': False})
    return candidates, sources


def _offer(pool, capacity, candidate):
    entry = (-candidate['_rank'], candidate['source_file'], candidate['record_ordinal'], candidate)
    if len(pool) < capacity:
        heapq.heappush(pool, entry)
    elif entry[:3] > pool[0][:3]:
        heapq.heapreplace(pool, entry)


def _choose(pools, target, seen, panel, selected):
    queues = {key: deque(entry[3] for entry in sorted(pool, key=lambda e: (-e[0], e[1], e[2])))
              for key, pool in pools.items()}
    active = deque(sorted(queues))
    count = 0
    while active and count < target:
        key = active.popleft()
        queue = queues[key]
        while queue and queue[0]['sequence_sha256'] in seen:
            queue.popleft()
        if queue:
            candidate = dict(queue.popleft())
            candidate['panel'] = panel
            selected.append(candidate)
            seen.add(candidate['sequence_sha256'])
            count += 1
        if queue:
            active.append(key)
    if count != target:
        raise ValueError(f'insufficient unique {panel} candidates: requested {target}, selected {count}; '
                         'increase --pool-factor or reduce corpus size')


def build_corpus(*, fasta_dir, manifest, annotations, output, size=50000,
                 difficult_fraction=.2, mandatory=None, pool_factor=4):
    fasta_dir, manifest, annotations = (Path(p).resolve(strict=True) for p in
                                         (fasta_dir, manifest, annotations))
    output = Path(output).resolve()
    if output.exists() or not output.parent.is_dir() or any(
            output.is_relative_to(p) for p in (fasta_dir, annotations, ROOT)) or output == manifest:
        raise ValueError('output must be new, have an existing parent, and be outside source/input trees')
    if not fasta_dir.is_dir() or not annotations.is_dir() or not manifest.is_file():
        raise ValueError('explicit input paths must be directories and a manifest file')
    if isinstance(size, bool) or size <= 0 or not 0 < difficult_fraction < 1 or pool_factor < 1:
        raise ValueError('size/pool-factor must be positive and difficult fraction between zero and one')
    with manifest.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        if not {'dataset', 'donor', 'flow_class'} <= set(reader.fieldnames or []):
            raise ValueError('manifest requires dataset, donor, flow_class')
        samples = sorted(reader, key=lambda row: row['dataset'])
    if not samples or len({s['dataset'] for s in samples}) != len(samples):
        raise ValueError('manifest datasets must be nonempty and unique')
    for sample in samples:
        if any(not sample[k] or '\0' in sample[k] for k in ('dataset', 'donor', 'flow_class')):
            raise ValueError('manifest fields must be nonempty strings')
        dataset = sample['dataset']
        if Path(dataset).name != dataset or dataset in ('.', '..'):
            raise ValueError('manifest dataset must be a single safe file component')
    fixture_sources = []
    if mandatory is None:
        mandatory, fixture_sources = mandatory_fixtures()
    requirements = defaultdict(list)
    for number, item in enumerate(mandatory):
        requirements[(item['source_file'], item['sequence_id'], item['sequence_sha256'])].append((number, item))
    found, mandatory_candidates = set(), {}
    background, difficult = defaultdict(list), defaultdict(list)
    difficult_count = round(size * difficult_fraction)
    capacity = max(32, math.ceil(size / len(samples)) * pool_factor)
    difficult_capacity = max(16, math.ceil(difficult_count / len(samples)) * pool_factor)
    sources, counts, difficulty_counts = [], Counter(), Counter()
    total = 0
    for sample in samples:
        filename = sample['dataset'] + '.fasta'
        fasta = (fasta_dir / filename).resolve(strict=True)
        parquet = (annotations / (sample['dataset'] + '.parquet')).resolve(strict=True)
        if not fasta.is_relative_to(fasta_dir) or not parquet.is_relative_to(annotations):
            raise ValueError('input files must stay inside their explicit roots')
        schema = pl.read_parquet_schema(parquet)
        evidence = annotation_index(pl.read_parquet(parquet, columns=[c for c in EVIDENCE_COLUMNS if c in schema]))
        source_count = 0
        for ordinal, sid, sequence in fasta_records(fasta):
            source_count += 1
            row = match_annotation(evidence, sid, sequence)
            reasons = difficulty_reasons(row, sequence)
            locus = (row or {}).get('locus') or 'unknown'
            digest = sequence_hash(sequence)
            candidate = {'source_file': filename, 'record_ordinal': ordinal, 'sequence_id': sid,
                'sequence_sha256': digest, 'donor': sample['donor'], 'flow_class': sample['flow_class'],
                'locus': locus, 'selection_reasons': reasons, '_rank': stable_rank(filename, ordinal)}
            stratum = (sample['donor'], sample['flow_class'], locus, sample['dataset'])
            counts[stratum] += 1
            _offer(background[stratum], capacity, candidate)
            for reason in reasons:
                difficulty_counts[reason] += 1
                _offer(difficult[(reason,) + stratum], difficult_capacity, candidate)
            matches = requirements.get((filename, sid, digest), [])
            if matches:
                key = (filename, digest)
                previous = mandatory_candidates.get(key)
                if previous is None or candidate['_rank'] < previous['_rank']:
                    mandatory_candidates[key] = dict(candidate)
                entry = mandatory_candidates[key]
                previous_reasons = set(previous['selection_reasons']) if previous else set()
                entry['selection_reasons'] = sorted(set(entry['selection_reasons']) | previous_reasons
                                                    | {item['reason'] for _, item in matches})
                found.update(number for number, _ in matches)
        if not source_count:
            raise ValueError(f'empty FASTA dataset: {filename}')
        total += source_count
        sources.append({'source_file': filename, 'fasta_sha256': sha256(fasta),
            'annotation_file': parquet.name, 'annotation_sha256': sha256(parquet), 'record_count': source_count})
        print(f'Indexed {filename}: {source_count} original reads', file=sys.stderr, flush=True)
    absent = [item for n, item in enumerate(mandatory) if n not in found]
    if any(not item.get('allow_absent', False) for item in absent):
        raise ValueError(f'mandatory original sequence not found: {absent}')
    selected, seen = [], set()
    for candidate in sorted(mandatory_candidates.values(), key=lambda c: c['_rank']):
        if candidate['sequence_sha256'] not in seen:
            candidate['panel'] = 'difficult'
            selected.append(candidate)
            seen.add(candidate['sequence_sha256'])
        else:
            retained = next(r for r in selected if r['sequence_sha256'] == candidate['sequence_sha256'])
            retained['selection_reasons'] = sorted(set(retained['selection_reasons']) | set(candidate['selection_reasons']))
    if len(selected) > difficult_count:
        raise ValueError('mandatory originals exceed difficult panel size')
    _choose(difficult, difficult_count - len(selected), seen, 'difficult', selected)
    _choose(background, size - difficult_count, seen, 'background', selected)
    selected.sort(key=lambda row: (row['source_file'], row['record_ordinal']))
    by_source = defaultdict(dict)
    for ordinal, row in enumerate(selected):
        row.pop('_rank')
        row['corpus_ordinal'] = ordinal
        by_source[row['source_file']][row['record_ordinal']] = row
    # Re-read original queries only for selected rows; never store all sequences.
    with tempfile.TemporaryDirectory(prefix='.bcr-ci-build-', dir=output.parent) as staging:
        staging = Path(staging)
        with (staging / 'sequences.fasta.gz').open('wb') as raw:
            with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as compressed:
                for filename in sorted(by_source):
                    pending = dict(by_source[filename])
                    for ordinal, sid, sequence in fasta_records(fasta_dir / filename):
                        if ordinal in pending:
                            row = pending.pop(ordinal)
                            if row['sequence_id'] != sid or row['sequence_sha256'] != sequence_hash(sequence):
                                raise ValueError('source changed during selection')
                            compressed.write(f'>{sid}\n{sequence}\n'.encode())
                    if pending:
                        raise ValueError('source changed during selection: missing records')
        pl.DataFrame(selected, schema=RECORD_SCHEMA).write_parquet(staging / 'records.parquet', compression='zstd')
        metadata = {'schema_version': 1, 'corpus_version': 'bcr-ci-v1', 'record_count': size,
            'files': {name: sha256(staging / name) for name in ('sequences.fasta.gz', 'records.parquet')},
            'parameters': {'receptor': 'bcr', 'germline_database': 'human', 'n_processes': 2,
                           'mmseqs_threads': 2, 'chunksize': 500},
            'expected_failures': [],
            'selection': {'size': size, 'difficult_fraction': difficult_fraction, 'pool_factor': pool_factor,
                'algorithm': 'sha256(bcr-ci-v1 NUL source_file NUL zero-based ordinal); smallest ranks in bounded per-stratum pools; round-robin sorted strata',
                'background_strata': ['donor', 'flow_class', 'locus', 'dataset'],
                'difficult_strata': ['mechanism', 'donor', 'flow_class', 'locus', 'dataset'],
                'deduplication': 'exact case-sensitive sequence SHA256 across both panels; mandatory first, then difficult, then background; first selected identity retained; mandatory reasons unioned',
                'background_sampling': 'all original reads eligible, including difficult reads; panel membership is exclusive',
                'difficulty_thresholds': {'high_shm_v_identity_below': .9, 'v_truncation_start_at_least': 15,
                    'short_cdr3_aa_at_most': 5, 'long_igh_cdr3_aa_at_least': 25, 'long_light_cdr3_aa_at_least': 13},
                'panel_counts': dict(Counter(r['panel'] for r in selected)),
                'reason_counts': dict(sorted(Counter(reason for r in selected for reason in r['selection_reasons']).items())),
                'source_difficulty_counts': dict(sorted(difficulty_counts.items())),
                'source_record_count': total, 'source_strata_count': len(counts),
                'selected_strata_count': len({(r['donor'], r['flow_class'], r['locus'], r['source_file']) for r in selected}),
                'mandatory_candidates': len(mandatory), 'matched_mandatory_candidates': len(found),
                'absent_fixture_candidates': absent,
                'fixture_policy': 'exact source file, ID, and original-sequence hash must match; absent assignment fixtures excluded; curated real_bcr originals required; derived fixtures never loaded'},
            'sources': {'manifest_sha256': sha256(manifest), 'files': sources, 'fixtures': fixture_sources,
                        'builder_sha256': sha256(__file__)}}
        (staging / 'manifest.json').write_text(json.dumps(metadata, indent=2, sort_keys=True) + '\n')
        staging.rename(output)
    return metadata


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta-dir', required=True, type=Path)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--annotations', required=True, type=Path,
                        help='Existing native Parquet directory; sampling labels only')
    parser.add_argument('--output', required=True, type=Path, help='New directory outside checkout and input trees')
    parser.add_argument('--size', type=int, default=50000)
    parser.add_argument('--difficult-fraction', type=float, default=.2)
    parser.add_argument('--pool-factor', type=int, default=4,
                        help='Bounded candidate pool multiplier; increase on deduplication exhaustion')
    args = parser.parse_args(argv)
    try:
        result = build_corpus(**vars(args))
    except (ValueError, OSError) as error:
        parser.exit(1, f'{error}\n')
    print(json.dumps({'record_count': result['record_count'], 'selection': result['selection']}, indent=2))


if __name__ == '__main__':
    main()
