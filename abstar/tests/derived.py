# Copyright (c) 2026 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Reproducible edits of adjudicated reads, without running annotation.

Offsets always address the ORIGINAL INPUT QUERY, zero-based and half-open.
Insertions occur before the base at offset (len(parent) appends). Affected
bases describe [offset:offset + len(affected_bases)]. A 5' truncation removes
that prefix at offset zero; a 3' truncation removes that suffix at its original
start. Reverse complement acts on the whole input, with offset zero and empty
payloads; an original span [a:b] becomes [n-b:n-a]. Operations are single edits,
never a pipeline with implicitly changing coordinate systems.
"""

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from abstar.tests import corpus


_DNA = frozenset('ACGTRYSWKMBDHVN')
_COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')


@dataclass(frozen=True, slots=True)
class DerivedOperation:
    kind: str
    offset: int
    affected_bases: str
    replacement_bases: str


def derive_sequence(parent: str, operation: DerivedOperation) -> str:
    """Apply one checked edit; reject invalid or empty resulting DNA with ValueError."""
    if not isinstance(parent, str) or not parent or set(parent) - _DNA:
        raise ValueError('parent must be nonempty uppercase IUPAC DNA')
    if not isinstance(operation, DerivedOperation):
        raise ValueError('operation must be a DerivedOperation')
    kind, offset, affected, replacement = (
        operation.kind, operation.offset, operation.affected_bases, operation.replacement_bases,
    )
    if type(offset) is not int or not 0 <= offset <= len(parent):
        raise ValueError('offset outside original input query')
    if any(not isinstance(s, str) or set(s) - _DNA for s in (affected, replacement)):
        raise ValueError('payloads must be uppercase IUPAC DNA strings')
    end = offset + len(affected)
    if end > len(parent) or parent[offset:end] != affected:
        raise ValueError('affected bases disagree with original input span')
    if kind in ('identity', 'reverse_complement'):
        valid = offset == 0 and not affected and not replacement
    elif kind == 'substitute':
        valid = (len(affected) == len(replacement) == 1 and affected != replacement
                 and replacement in 'ACGT' and affected in 'ACGT')
    elif kind == 'ambiguity':
        valid = (len(affected) == len(replacement) == 1 and affected in 'ACGT'
                 and replacement in 'RYSWKMBDHVN')
    elif kind == 'insert':
        valid = not affected and bool(replacement)
    elif kind in ('delete', 'truncate_5prime', 'truncate_3prime'):
        valid = bool(affected) and not replacement
        if kind == 'truncate_5prime':
            valid = valid and offset == 0
        elif kind == 'truncate_3prime':
            valid = valid and end == len(parent)
    else:
        raise ValueError('unsupported derived operation')
    if not valid:
        raise ValueError('invalid operation payload combination')
    if kind == 'reverse_complement':
        return parent.translate(_COMPLEMENT)[::-1]
    result = parent[:offset] + replacement + parent[end:]
    if not result:
        raise ValueError('operation removes the entire parent')
    return result


@dataclass(frozen=True, slots=True)
class DerivedBCRCase:
    case_id: str
    parent: Mapping[str, str]
    context: str
    operation: DerivedOperation
    sequence: str
    sequence_sha256: str
    expected: Mapping[str, object]
    evidence: tuple[str, ...]

    def __post_init__(self):
        for field in ('parent', 'expected', 'evidence'):
            object.__setattr__(self, field, corpus._freeze(getattr(self, field)))


def _validate_context(parent, operation, context):
    """Protect conserved anchors; edits have one explicitly supported context."""
    expected = parent.expected
    start, end = operation.offset, operation.offset + len(operation.affected_bases)
    if operation.kind in ('identity', 'reverse_complement'):
        valid = context == 'whole'
    elif operation.kind == 'truncate_5prime':
        valid = context == '5prime' and end <= expected['v_sequence_start']
    elif operation.kind == 'truncate_3prime':
        valid = context == '3prime' and start >= expected['j_sequence_end']
    elif context == 'v':
        # Insertion at the first V boundary could simply extend the 5' flank;
        # it does not establish an internal coding-frame displacement.
        coding_start = parent.source['alignment']['coding_start']
        valid = coding_start < start < expected['junction_start'] and end <= expected['junction_start']
    elif context == 'junction':
        valid = (expected['junction_start'] + 3 <= start < expected['junction_end'] - 3
                 and end <= expected['junction_end'] - 3)
    else:
        valid = False
    if not valid:
        raise ValueError('operation does not lie within the stated parent evidence context')


def _project_expectations(parent, operation, sequence, context):
    """Check literal expectations using only sequence arithmetic and parent anchors.

    These are homologous-anchor projections, NOT annotation output assertions.
    In particular, frameshifting indels do not justify claims about the selected
    alignment, translated junction, productivity, or exact gene/indel calls.
    coding_frame_delta_mod3 is the signed coding-length change modulo three,
    measured between the retained V coding origin and conserved J anchor.
    Truncations remove only flanks outside those retained coding traces.
    """
    delta = len(sequence) - len(parent.sequence)
    junction_start = parent.expected['junction_start']
    junction_end = parent.expected['junction_end']
    if operation.offset < junction_start:
        junction_start += delta
    if operation.offset < junction_end:
        junction_end += delta
    rev_comp = operation.kind == 'reverse_complement'
    oriented = sequence.translate(_COMPLEMENT)[::-1] if rev_comp else sequence
    input_start, input_end = junction_start, junction_end
    if rev_comp:
        input_start, input_end = len(sequence) - junction_end, len(sequence) - junction_start
    junction = oriented[junction_start:junction_end]
    return {
        'sequence_length': len(sequence), 'length_delta': delta,
        'rev_comp': rev_comp, 'locus': parent.expected['locus'],
        'homologous_junction': {
            'oriented_start': junction_start, 'oriented_end': junction_end,
            'input_start': input_start, 'input_end': input_end,
            'sequence': junction, 'length_mod3': len(junction) % 3,
        },
        'coding_frame_delta_mod3': delta % 3 if context in ('v', 'junction') else 0,
        'ambiguous_base_count': sum(base not in 'ACGT' for base in sequence),
    }


def load_derived_bcr_cases(path=None) -> tuple[DerivedBCRCase, ...]:
    """Load fresh immutable variants from checked committed original parents.

    The v1 schema accepts only productive, forward, ambiguity-free concordant
    parents with gap-free V/J traces and primary-J coding evidence. Parents are
    always loaded by the authenticated original fixture loader; a custom path
    selects only derived JSON, never an unvalidated replacement parent corpus.
    No full generated input sequences are stored in JSON. Every call rebuilds
    them and checks parent/result hashes, edit context and exact expectations.
    """
    path = corpus.REAL_BCR_DIRECTORY / 'derived_cases.json' if path is None else Path(path)
    document = json.loads(path.read_text(encoding='utf-8'))
    if (not isinstance(document, dict) or document.keys() != {'schema_version', 'cases'}
            or type(document['schema_version']) is not int or document['schema_version'] != 1
            or not isinstance(document['cases'], list) or not document['cases']):
        raise ValueError('invalid derived case document schema')
    parents = {(p.dataset, p.sequence_id): p for p in corpus.load_real_bcr_cases()}
    case_fields = set(DerivedBCRCase.__dataclass_fields__) - {'sequence'}
    operation_fields = set(DerivedOperation.__dataclass_fields__)
    cases, seen = [], set()
    for raw in document['cases']:
        if not isinstance(raw, dict) or raw.keys() != case_fields:
            raise ValueError('invalid derived case schema')
        corpus._validate_identity(raw['case_id'], 'case_id')
        if raw['case_id'] in seen:
            raise ValueError('duplicate derived case_id')
        seen.add(raw['case_id'])
        key = raw['parent']
        if not isinstance(key, dict) or key.keys() != {'dataset', 'sequence_id', 'sequence_sha256'}:
            raise ValueError('invalid parent key schema')
        for field in ('dataset', 'sequence_id'):
            corpus._validate_identity(key[field], field)
        parent = parents.get((key['dataset'], key['sequence_id']))
        if parent is None or key['sequence_sha256'] != parent.sequence_sha256:
            raise ValueError('unknown parent or parent SHA-256 mismatch')
        alignment = parent.source['alignment']
        if (parent.selection_reasons != ('concordant_' + parent.expected['locus'],)
                or parent.expected['rev_comp'] or not parent.expected['productive']
                or set(parent.sequence) - set('ACGT')
                or alignment['coding_scope'] != 'through_primary_j'
                or any(alignment[s]['indels'] for s in ('v', 'j'))):
            raise ValueError('derived cases require clean adjudicated forward parents')
        operation = raw['operation']
        if not isinstance(operation, dict) or operation.keys() != operation_fields:
            raise ValueError('invalid operation schema')
        operation = DerivedOperation(**operation)
        sequence = derive_sequence(parent.sequence, operation)
        if hashlib.sha256(sequence.encode('ascii')).hexdigest() != raw['sequence_sha256']:
            raise ValueError('derived sequence SHA-256 mismatch')
        _validate_context(parent, operation, raw['context'])
        expected = _project_expectations(parent, operation, sequence, raw['context'])
        # JSON comparison preserves Boolean vs integer types (False != 0 here),
        # exact nested keys, integer coordinates and explicit null semantics.
        if json.dumps(raw['expected'], sort_keys=True) != json.dumps(expected, sort_keys=True):
            raise ValueError('derived expectations disagree with transformation and parent evidence')
        corpus._nonempty_strings(raw['evidence'], 'evidence')
        cases.append(DerivedBCRCase(**{**raw, 'operation': operation, 'sequence': sequence}))
    return tuple(cases)
