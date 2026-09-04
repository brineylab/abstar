# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""Shared helpers for the abstar test suite."""

from pathlib import Path


def read_fasta_records(path: Path) -> dict[str, str]:
    """Read FASTA records, rejecting duplicate identifiers."""
    records = {}
    identifier = None
    sequence = []

    with path.open(encoding="utf-8") as fasta:
        for raw_line in fasta:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(sequence)
                identifier = line[1:].split()[0]
                if identifier in records:
                    raise ValueError(
                        f"duplicate FASTA identifier '{identifier}' in {path}"
                    )
                sequence = []
            else:
                sequence.append(line)

    if identifier is not None:
        records[identifier] = "".join(sequence)
    return records
