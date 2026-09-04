# Copyright (c) 2026 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import os
import subprocess

import pytest
from abutils import Sequence
from abutils.io import parse_fastx

from ..core.abstar import run
from ..preprocess import merging


def _touch_fastqs(tmp_path, names):
    paths = []
    for name in names:
        path = tmp_path / name
        path.write_text("@read\nACGT\n+\nIIII\n")
        paths.append(str(path))
    return paths


def test_fastq_suffix_removal_preserves_stem_characters():
    assert merging._strip_fastq_suffix("iraq.fastq.gz") == "iraq"
    assert merging._strip_fastq_suffix("staff.fq") == "staff"


def test_group_paired_fastqs_requires_one_r1_and_r2_per_lane(tmp_path):
    files = _touch_fastqs(
        tmp_path,
        ["sample_S1_L001_R1_001.fastq", "sample_S1_L001_R2_001.fastq"],
    )

    groups = merging.group_paired_fastqs(files)

    assert len(groups) == 1
    assert {fastq.read for fastq in groups[0].files} == {"R1", "R2"}


@pytest.mark.parametrize(
    "names",
    [
        ["sample_S1_L001_R1_001.fastq"],
        [
            "sample_S1_L001_R1_001.fastq",
            "sample_S1_L001_R1_002.fastq",
            "sample_S1_L001_R2_001.fastq",
        ],
    ],
)
def test_group_paired_fastqs_rejects_missing_or_duplicate_reads(tmp_path, names):
    files = _touch_fastqs(tmp_path, names)

    with pytest.raises(ValueError, match="exactly one R1 and one R2"):
        merging.group_paired_fastqs(files)


def test_fastp_uses_checked_argument_list(monkeypatch, tmp_path):
    forward, reverse = _touch_fastqs(tmp_path, ["R 1.fastq", "R 2.fastq"])
    merged = tmp_path / "merged reads.fastq"
    observed = {}

    def fake_run(command, **kwargs):
        observed["command"] = command
        observed["kwargs"] = kwargs
        return subprocess.CompletedProcess(command, 0, "stdout", "stderr")

    monkeypatch.setattr(merging.sp, "run", fake_run)

    merging.merge_fastqs_fastp(
        forward,
        reverse,
        str(merged),
        binary_path="/opt/fastp",
        log_directory=str(tmp_path / "logs with spaces"),
        additional_args="--thread 2 --report_title 'report; still one argument'",
    )

    command = observed["command"]
    assert isinstance(command, list)
    assert forward in command
    assert reverse in command
    assert str(merged) in command
    assert "report; still one argument" in command
    assert observed["kwargs"] == {
        "check": True,
        "capture_output": True,
        "text": True,
    }


def test_fastp_failure_includes_stderr(monkeypatch, tmp_path):
    forward, reverse = _touch_fastqs(tmp_path, ["R1.fastq", "R2.fastq"])

    def fail(command, **kwargs):
        raise subprocess.CalledProcessError(2, command, stderr="bad overlap")

    monkeypatch.setattr(merging.sp, "run", fail)

    with pytest.raises(ValueError, match="bad overlap"):
        merging.merge_fastqs_fastp(
            forward,
            reverse,
            str(tmp_path / "merged.fastq"),
            binary_path="fastp",
        )


@pytest.mark.integration
def test_fastp_merges_a_real_overlapping_pair(tmp_path):
    sequence = (
        "ACGTTGCAAGTCGATCGTACGATGCTAGCTACGTTAGCGATCGATGACCTGACTGATCGTAGCTAGTCGATG"
    )
    reverse_complement = Sequence(sequence).reverse_complement
    forward = tmp_path / "sample_R1.fastq"
    reverse = tmp_path / "sample_R2.fastq"
    merged = tmp_path / "sample.fastq"
    forward.write_text(f"@read\n{sequence}\n+\n{'I' * len(sequence)}\n")
    reverse.write_text(
        f"@read\n{reverse_complement}\n+\n{'I' * len(reverse_complement)}\n"
    )

    merging.merge_fastqs_fastp(
        str(forward),
        str(reverse),
        str(merged),
        minimum_overlap=30,
        trim_adapters=False,
        quality_trim=False,
    )

    merged_reads = list(parse_fastx(str(merged)))
    assert len(merged_reads) == 1
    assert merged_reads[0].sequence == sequence


def test_interleaved_merge_normalizes_text_and_cleans_temporary_file(
    monkeypatch, tmp_path
):
    interleaved = tmp_path / "iraq.fastq"
    interleaved.write_text(
        "@read/1\nACGT\n+\nIIII\n@read/2\nACGT\n+\nIIII\n"
    )
    observed = {}

    def fake_merge(*, forward, merged, **kwargs):
        observed["temporary"] = forward
        observed["content"] = open(forward).read()
        observed["merged"] = merged

    monkeypatch.setattr(merging, "merge_fastqs_fastp", fake_merge)

    outputs = merging.merge_fastqs(
        [str(interleaved)], str(tmp_path / "output"), interleaved=True
    )

    assert observed["content"].count("@read") == 2
    assert not os.path.exists(observed["temporary"])
    assert outputs == [str(tmp_path / "output" / "iraq.fastq")]


@pytest.mark.e2e
def test_run_merge_accepts_default_merge_kwargs(monkeypatch, tmp_path):
    sentinel = RuntimeError("merge reached")
    observed = {}

    def fake_merge(files, output_directory, **kwargs):
        observed["files"] = files
        observed["output_directory"] = output_directory
        observed["kwargs"] = kwargs
        raise sentinel

    monkeypatch.setattr("abstar.core.abstar.merge_fastqs", fake_merge)

    with pytest.raises(RuntimeError, match="merge reached"):
        run(
            Sequence("ACGT", id="read"),
            project_path=str(tmp_path / "project"),
            merge=True,
        )

    assert len(observed["files"]) == 1
    assert observed["output_directory"].endswith("merged")
    assert observed["kwargs"]["interleaved"] is False
