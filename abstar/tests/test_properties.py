# Copyright (c) 2026 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Bounded properties for annotation and input-boundary contracts."""

from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import polars as pl
import pytest
from abutils import Sequence
from abutils.tl import reverse_complement
from hypothesis import HealthCheck, example, given, note, seed, settings, strategies as st

from abstar.annotation.antibody import Antibody
from abstar.annotation.germline import (
    process_cgene_alignment,
    process_dgene_alignment,
    process_jgene_alignment,
    process_vgene_alignment,
)
from abstar.annotation.indels import annotate_deletions, annotate_insertions
from abstar.annotation.positions import (
    get_gapped_position_from_raw,
    get_raw_position_from_gapped,
)
from abstar.assigners.mmseqs import select_best_hits
from abstar.core.abstar import _process_inputs
from abstar.tests.derived import DerivedOperation, derive_sequence


PROPERTY_SEED = 20260904
property_settings = settings(max_examples=200, deadline=None, print_blob=True)
dna = st.text(alphabet="ACGT", min_size=1, max_size=300)
external_ids = st.text(
    alphabet=st.characters(blacklist_categories=("Cs",), blacklist_characters="\t\n\r"),
    min_size=1,
    max_size=128,
)


@seed(PROPERTY_SEED)
@property_settings
@given(st.text('ACGTRYSWKMBDHVN', min_size=1, max_size=300))
def test_derived_reverse_complement_restores_original(sequence):
    """Catch loss of IUPAC information during orientation changes."""
    operation = DerivedOperation('reverse_complement', 0, '', '')
    assert derive_sequence(derive_sequence(sequence, operation), operation) == sequence


@seed(PROPERTY_SEED)
@property_settings
@given(dna, st.text('ACGT', min_size=1, max_size=12), st.data())
def test_derived_insertion_deletion_restores_original(sequence, payload, data):
    """Catch insertion boundary shifts, including prepend and append."""
    offset = data.draw(st.integers(0, len(sequence)))
    inserted = derive_sequence(sequence, DerivedOperation('insert', offset, '', payload))
    assert derive_sequence(inserted, DerivedOperation('delete', offset, payload, '')) == sequence


@seed(PROPERTY_SEED)
@property_settings
@given(dna, st.data())
def test_derived_substitution_changes_only_named_offset(sequence, data):
    """Catch a substitution that changes the wrong base or surrounding sequence."""
    offset = data.draw(st.integers(0, len(sequence) - 1))
    base = sequence[offset]
    replacement = data.draw(st.sampled_from(sorted(set('ACGT') - {base})))
    result = derive_sequence(sequence, DerivedOperation('substitute', offset, base, replacement))
    assert len(result) == len(sequence)
    assert [i for i, (a, b) in enumerate(zip(sequence, result)) if a != b] == [offset]
    assert result[offset] == replacement


@seed(PROPERTY_SEED)
@settings(max_examples=200, deadline=None, print_blob=True,
          suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(st.data())
def test_derived_operations_leave_loaded_parent_unchanged(real_bcr_cases, data):
    """Catch mutation of a shared real fixture or its Sequence adapters."""
    parent = data.draw(st.sampled_from(real_bcr_cases))
    original = (parent.sequence, parent.sequence_sha256, dict(parent.expected))
    offset = data.draw(st.integers(0, len(parent.sequence)))
    payload = data.draw(st.text('ACGT', min_size=1, max_size=3))
    result = derive_sequence(parent.sequence, DerivedOperation('insert', offset, '', payload))
    adapter = parent.as_sequence()
    adapter.sequence = result
    assert (parent.sequence, parent.sequence_sha256, dict(parent.expected)) == original
    assert parent.as_sequence().sequence == original[0]


@st.composite
def gapped_dna(draw, sequence_strategy=dna):
    """DNA with up to three IMGT dots at every boundary, including either end."""
    sequence = draw(sequence_strategy)
    gaps = draw(st.lists(st.integers(0, 3), min_size=len(sequence) + 1,
                         max_size=len(sequence) + 1))
    return "".join("." * gap + base for gap, base in zip(gaps, sequence)) + "." * gaps[-1]


@st.composite
def alignment_pairs(draw):
    """Matching DNA and indel blocks; no double gaps or substitutions.

    A shared anchor keeps both ungapped sequences nonempty. The anchor may be
    anywhere, so leading/trailing and adjacent insertion/deletion runs occur.
    """
    blocks = draw(st.lists(
        st.tuples(st.sampled_from("MID"), st.text("ACGT", min_size=1, max_size=9)),
        max_size=12,
    ))
    anchor = draw(st.integers(0, len(blocks)))
    blocks.insert(anchor, ("M", draw(st.text("ACGT", min_size=1, max_size=9))))
    query = "".join("-" * len(bases) if kind == "D" else bases for kind, bases in blocks)
    germline = "".join("-" * len(bases) if kind == "I" else bases for kind, bases in blocks)
    prefix = draw(st.text("ACGT", max_size=12))
    full_germline = prefix + germline.replace("-", "")
    gapped = draw(gapped_dna(st.just(full_germline)))
    return query, germline, gapped, len(prefix)


@seed(PROPERTY_SEED)
@property_settings
@given(dna)
@example("ACG")
def test_reverse_complement_is_involutive(sequence):
    """Catch a broken reverse-complement dependency used by sequence orientation."""
    note(f"Hypothesis seed: {PROPERTY_SEED}")
    assert reverse_complement(reverse_complement(sequence)) == sequence


@seed(PROPERTY_SEED)
@property_settings
@given(gapped_dna())
@example("...A..CG...T...")
def test_imgt_gapped_round_trip(gapped):
    """Catch off-by-one mapping and dropped leading/internal IMGT gaps."""
    note(f"Hypothesis seed: {PROPERTY_SEED}")
    sequence = gapped.replace(".", "")
    mapped = []
    for raw_position in range(1, len(sequence) + 1):
        imgt = get_gapped_position_from_raw(raw_position, gapped)
        assert 1 <= imgt <= len(gapped)
        assert gapped[imgt - 1] == sequence[raw_position - 1]
        assert get_raw_position_from_gapped(imgt, gapped) == raw_position
        mapped.append(imgt)
    assert mapped == sorted(set(mapped))


def _parse_indels(description):
    """Decode the public notation, without detecting indels from the alignment."""
    if not description:
        return []
    events = []
    for event in description.split("|"):
        coordinates, payload = event.split(":")
        length, bases = payload.split(">")
        frameshift = bases.endswith("!")
        bases = bases.removesuffix("!")
        start, _, end = coordinates.partition("-")
        assert int(length) == len(bases)
        assert frameshift == (len(bases) % 3 != 0)
        events.append((int(start), int(end or start), bases))
    return events


@seed(PROPERTY_SEED)
@property_settings
@given(alignment_pairs())
@example(("AG--TAC", "A-CGT-C", "..A.C..GTC.", 0))
def test_indel_descriptions_reconstruct_aligned_pair(pair):
    """Catch wrong IMGT endpoints, missing runs, payloads, or frame markers."""
    note(f"Hypothesis seed: {PROPERTY_SEED}")
    query, germline, gapped, prefix_length = pair
    insertions = _parse_indels(annotate_insertions(query, germline, gapped, prefix_length))
    deletions = _parse_indels(annotate_deletions(query, germline, gapped, prefix_length))

    # Decode by editing the numbered IMGT template; this is the inverse of
    # annotation, not a copy of its regex scan or raw-to-IMGT conversion.
    numbered = [(position, base) for position, base in enumerate(gapped, 1) if base != "."]
    numbered = numbered[prefix_length:]
    deleted = set()
    for start, end, bases in deletions:
        residues = [(position, base) for position, base in numbered if start <= position <= end]
        assert residues[0][0] == start
        assert residues[-1][0] == end
        assert "".join(base for _, base in residues) == bases
        positions = {position for position, _ in residues}
        assert deleted.isdisjoint(positions)
        deleted.update(positions)

    inserted = {}
    for start, end, bases in insertions:
        assert start == end
        assert start not in inserted
        assert start == 0 or gapped[start - 1] != "."
        inserted[start] = bases
    preceding = 0 if prefix_length == 0 else [i for i, base in enumerate(gapped, 1) if base != "."][prefix_length - 1]
    reconstructed_query = inserted.pop(preceding, "")
    reconstructed_germline = "-" * len(reconstructed_query)
    for position, base in numbered:
        added = inserted.pop(position, "")
        reconstructed_query += ("-" if position in deleted else base) + added
        reconstructed_germline += base + "-" * len(added)
    assert not inserted
    assert (reconstructed_query, reconstructed_germline) == (query, germline)


def _exact_alignment(query, target, query_begin=0):
    """An exact-match alignment record in its own query coordinate space."""
    return SimpleNamespace(
        query=query, target=target, query_begin=query_begin,
        query_end=len(query) - 1, target_begin=0, target_end=len(target) - 1,
        aligned_query=query, aligned_target="-" * query_begin + target,
        score=2 * len(target),
    )


@seed(PROPERTY_SEED)
@property_settings
@given(
    st.lists(st.text("ACGT", min_size=3, max_size=60), min_size=4, max_size=4),
    st.lists(st.text("ACGT", max_size=12), min_size=4, max_size=4),
)
def test_segment_intervals_are_ordered_and_reconstruct_query(segments, spacers):
    """Catch lost/doubled query offsets and inclusive/exclusive endpoint errors."""
    note(f"Hypothesis seed: {PROPERTY_SEED}")
    windows = [spacer + segment for spacer, segment in zip(spacers, segments)]
    query = "".join(windows)
    ab = Antibody(sequence_id="intervals", d_call="IGHD1-1*01")
    local = [_exact_alignment(segment, segment) for segment in segments]
    v_sg = _exact_alignment(windows[0], segments[0], len(spacers[0]))
    process_vgene_alignment(v_sg, local[0], ab)
    d_loc = _exact_alignment(windows[1], segments[1], len(spacers[1]))
    process_dgene_alignment(query, ab.v_sequence_end, d_loc, ab)
    j_window = windows[1] + windows[2]
    j_sg = _exact_alignment(j_window, segments[2], len(windows[1]) + len(spacers[2]))
    process_jgene_alignment(query, ab.v_sequence_end, j_sg, local[2], ab)
    c_sg = _exact_alignment(windows[3], segments[3], len(spacers[3]))
    process_cgene_alignment(query, ab.j_sequence_end, c_sg, local[3], ab)

    previous_end = 0
    reconstructed = ""
    for name, segment, spacer in zip("vdjc", segments, spacers):
        start = getattr(ab, f"{name}_sequence_start")
        end = getattr(ab, f"{name}_sequence_end")
        assert previous_end <= start < end <= len(query)
        assert query[previous_end:start] == spacer
        assert query[start:end] == getattr(ab, f"{name}_sequence") == segment
        reconstructed += spacer + getattr(ab, f"{name}_sequence")
        previous_end = end
    assert previous_end == len(query)
    assert reconstructed == query


@seed(PROPERTY_SEED)
@property_settings
@given(st.data())
def test_best_hits_are_permutation_invariant_with_sorted_allele_ties(data):
    """Catch row-order ranking, dropped/unsorted ties, and unstable detail selection."""
    note(f"Hypothesis seed: {PROPERTY_SEED}")
    rows = []
    expected_calls = {}
    expected_evidence = []
    for query_number in range(data.draw(st.integers(1, 3))):
        query = f"query-{query_number}"
        alleles = data.draw(st.lists(st.integers(1, 99), min_size=2, max_size=5, unique=True))
        calls = sorted(f"IGHV1-2*{allele:02}" for allele in alleles)
        expected_calls[query] = ",".join(calls)
        bits = data.draw(st.integers(100, 1000))
        support = data.draw(st.sampled_from([1e-30, 1e-15, 1e-5]))
        expected_evidence.append({
            "v_query": query, "v_bits": float(bits), "v_support": support,
            "v_fident": 0.9, "v_tcov": 0.8, "v_qcov": 0.7, "v_alnlen": 60,
            "v_qstart": 1, "v_qend": 60, "v_qseq": "A" * 100,
        })
        for index, call in enumerate(calls):
            rows.append(dict(v_query=query, v_call=call, v_bits=float(bits),
                             v_support=support, v_fident=0.9, v_tcov=0.8,
                             v_qcov=0.7, v_alnlen=60, v_qstart=index + 1,
                             v_qend=index + 60, v_qseq="A" * 100))
            shift = data.draw(st.integers(1, 20))
            # Same allele and all ranking metrics, different aligned interval.
            # Both intervals remain within the same full query sequence.
            rows.append(dict(rows[-1], v_qstart=index + 1 + shift,
                             v_qend=index + 60 + shift))
        # A weaker bit score loses even with stronger secondary evidence.
        rows.append(dict(rows[-1], v_call="IGHV9-9*01", v_bits=float(bits - 1),
                         v_support=support / 10, v_fident=1.0))
        rows.append(dict(rows[-2]))  # identical rows must not duplicate an allele
    frame = pl.DataFrame(rows)
    expected = select_best_hits(frame, "v").sort("v_query")
    assert dict(expected.select("v_query", "v_call").iter_rows()) == expected_calls
    assert expected.select(list(expected_evidence[0])).to_dicts() == expected_evidence
    permutations = [list(reversed(range(frame.height))),
                    data.draw(st.permutations(range(frame.height)))]
    for permutation in permutations:
        actual = select_best_hits(frame[list(permutation)], "v").sort("v_query")
        assert actual.equals(expected)
    assert select_best_hits(frame.lazy(), "v").sort("v_query").equals(expected)


class _SinglePassIterator:
    """Expose accidental repeated traversal of a consumable input source."""
    def __init__(self, records):
        self.records = iter(records)
        self.yielded = 0
        self.exhaustions = 0

    def __iter__(self):
        # Reacquiring an iterator does not consume it. Generator expressions
        # may call this more than once, depending on the Python version.
        return self

    def __next__(self):
        try:
            record = next(self.records)
        except StopIteration:
            self.exhaustions += 1
            assert self.exhaustions == 1, "exhausted input was consumed again"
            raise
        self.yielded += 1
        return record


def test_single_pass_guard_allows_repeated_iter_calls():
    source = _SinglePassIterator(["first", "second"])
    assert iter(source) is source
    assert iter(source) is source
    assert list(source) == ["first", "second"]
    assert source.yielded == 2
    assert source.exhaustions == 1


def test_single_pass_guard_rejects_actual_second_consumption():
    source = _SinglePassIterator(["first", "second"])
    assert list(source) == ["first", "second"]
    with pytest.raises(AssertionError, match="exhausted input was consumed again"):
        list(source)


@pytest.mark.parametrize("input_kind", ["list", "iterator", "generator"])
@seed(PROPERTY_SEED)
@settings(
    max_examples=200, deadline=None, print_blob=True,
    # Each example owns a fresh subdirectory below pytest's temporary root.
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(st.lists(st.tuples(external_ids, dna), min_size=1, max_size=8))
@example([("123", "A"), ("0001", "C"), ("10E8", "G"), ("β抗体", "T"),
          (" leading and trailing ", "AC"), ("x" * 128, "GT"), ("10E8", "CG")])
def test_process_inputs_preserves_headers_and_consumes_once(tmp_path, input_kind, records):
    """Catch identifier normalization, duplicate loss, and consuming generators twice."""
    note(f"Hypothesis seed: {PROPERTY_SEED}")
    sequences = [Sequence(sequence, id=identifier) for identifier, sequence in records]
    source = _SinglePassIterator(sequences)
    if input_kind == "list":
        inputs = sequences
    elif input_kind == "iterator":
        inputs = source
    else:
        inputs = (sequence for sequence in source)
    with TemporaryDirectory(prefix="abstar-properties-", dir=tmp_path) as temp_dir:
        paths = _process_inputs(inputs, temp_dir)
        assert len(paths) == 1
        payload = Path(paths[0]).read_bytes()
        expected = "\n".join(f">{identifier}\n{sequence}" for identifier, sequence in records)
        assert payload == expected.encode("utf-8")
        # Split only on the actual line delimiter: IDs may contain Unicode
        # whitespace that bytes/str.splitlines() treats as extra line breaks.
        lines = payload.split(b"\n")
        assert len(lines) == 2 * len(records)
        assert lines[::2] == [b">" + identifier.encode("utf-8") for identifier, _ in records]
        assert lines[1::2] == [sequence.encode("ascii") for _, sequence in records]
        if input_kind != "list":
            assert source.yielded == len(records)
            assert source.exhaustions == 1
    assert not Path(temp_dir).exists()
