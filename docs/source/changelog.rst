.. _changelog:

Changelog
=========

0.9.0 (unreleased)
------------------

Breaking changes
~~~~~~~~~~~~~~~~

Every public output (AIRR TSV, Parquet, and the dataframe or ``Sequence``
objects returned by ``abstar.run()``) now has the same fields, in the same
order, with the same meanings. See the breaking-change notice in
:ref:`output-formats` for a before/after example.

- Coordinates are AIRR 1-based closed intervals in every output. Every
  ``*_start`` value is one greater than in 0.8.x; every ``*_end`` value is
  unchanged. To compare with 0.8.x output, add one to its starts.
- ``sequence`` is the original input and ``sequence_aa`` translates the full
  oriented query. In 0.8.x both held the assembled V(D)J; use the new
  ``sequence_vdj`` / ``sequence_vdj_aa`` fields instead. A returned
  ``Sequence`` object's ``.sequence`` is now the input.
- ``sequence_alignment_aa`` and ``germline_alignment_aa`` translate the paired
  nucleotide alignment columns.
- AIRR TSV writes ``T``/``F`` booleans and empty null cells, and the
  AIRR-required fields lead the field order in every output.

Added
~~~~~

- ``sequence_vdj`` and ``sequence_vdj_aa``: the assembled V(D)J nucleotide
  sequence (no constant region) and its translation, for every locus. The
  gene-segment, CDR, and non-germline masks index these fields.
