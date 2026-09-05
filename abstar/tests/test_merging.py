# Copyright (c) 2026 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

import os
import subprocess
from pathlib import Path

import pytest
from abstar import AnnotationRunError, pp, run
from abutils import Sequence
from abutils.bin import get_path as get_binary_path
from abutils.io import parse_fastx

from ..preprocess import merging


def _touch_fastqs(tmp_path, names):
    paths = []
    for name in names:
        path = tmp_path / name
        path.write_text("@read\nACGT\n+\nIIII\n")
        paths.append(str(path))
    return paths


PAIR_SEQUENCES = (
    "ACGTTGCAAGTCGATCGTACGATGCTAGCTACGTTAGCGATCGATGACCTGACTGATCGTAGCTAGTCGATGACGATCT",
    "TGCACTGATCGACGTTAGCATCGATGCATGCTAGCATCGTTACGATCGGATCCGATGCTAGCATCGATGATCGTACGTA",
)
PAIR_IDS = ("pair-A", "pair-B")


def _fastq_record(identifier, sequence, mate):
    return f"@{identifier}/{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n"


@pytest.fixture
def paired_fastq_inputs(tmp_path):
    input_directory = tmp_path / "paired"
    input_directory.mkdir()
    forward = input_directory / "sample_S1_L001_R1_001.fastq"
    reverse = input_directory / "sample_S1_L001_R2_001.fastq"
    forward.write_text("".join(
        _fastq_record(identifier, sequence, 1)
        for identifier, sequence in zip(PAIR_IDS, PAIR_SEQUENCES)
    ))
    reverse.write_text("".join(
        _fastq_record(identifier, Sequence(sequence).reverse_complement, 2)
        for identifier, sequence in zip(PAIR_IDS, PAIR_SEQUENCES)
    ))
    forward_reads = list(parse_fastx(str(forward)))
    reverse_reads = list(parse_fastx(str(reverse)))
    assert len(forward_reads) == len(reverse_reads) == 2
    assert [read.id.removesuffix("/1") for read in forward_reads] == list(PAIR_IDS)
    assert [read.id.removesuffix("/2") for read in reverse_reads] == list(PAIR_IDS)
    return input_directory, forward, reverse


@pytest.fixture
def interleaved_fastq_input(tmp_path):
    interleaved = tmp_path / "interleaved.fastq"
    interleaved.write_text("".join(
        _fastq_record(identifier, sequence, 1)
        + _fastq_record(identifier, Sequence(sequence).reverse_complement, 2)
        for identifier, sequence in zip(PAIR_IDS, PAIR_SEQUENCES)
    ))
    reads = list(parse_fastx(str(interleaved)))
    assert len(reads) == 4
    assert [read.id.removesuffix(f"/{mate}") for read, mate in zip(reads, (1, 2, 1, 2))] == [
        "pair-A", "pair-A", "pair-B", "pair-B",
    ]
    return interleaved


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
        additional_args=["--thread", "2", "--report_title", "report; $(touch never)"],
    )

    command = observed["command"]
    assert isinstance(command, list)
    assert forward in command
    assert reverse in command
    assert str(merged) in command
    assert "report; $(touch never)" in command
    assert observed["kwargs"] == {
        "check": True,
        "capture_output": True,
        "text": True,
    }


def test_fastp_cleans_implicit_report_directory(monkeypatch, tmp_path):
    forward, reverse = _touch_fastqs(tmp_path, ["R1.fastq", "R2.fastq"])
    report_directory = tmp_path / "implicit-fastp-reports"

    def make_report(command, **kwargs):
        Path(command[command.index("--html") + 1]).write_text("html")
        Path(command[command.index("--json") + 1]).write_text("json")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(
        merging.tempfile, "mkdtemp", lambda **kwargs: str(report_directory)
    )
    monkeypatch.setattr(merging.sp, "run", make_report)

    merging.merge_fastqs_fastp(
        forward,
        reverse,
        str(tmp_path / "merged.fastq"),
        binary_path="fastp",
    )

    assert not report_directory.exists()


def test_fastp_failure_includes_process_diagnostics(monkeypatch, tmp_path):
    forward, reverse = _touch_fastqs(tmp_path, ["R1.fastq", "R2.fastq"])

    def fail(command, **kwargs):
        Path(command[command.index("--merged_out") + 1]).write_text("partial output")
        raise subprocess.CalledProcessError(
            37, command, output="TASK18-STDOUT", stderr="TASK18-STDERR"
        )

    monkeypatch.setattr(merging.sp, "run", fail)

    with pytest.raises(merging.MergeExternalToolError) as captured:
        merging.merge_fastqs_fastp(
            forward,
            reverse,
            str(tmp_path / "merged.fastq"),
            binary_path="fastp",
        )
    error = captured.value
    assert error.returncode == 37
    assert error.stdout == "TASK18-STDOUT"
    assert error.stderr == "TASK18-STDERR"
    assert "EXIT STATUS: 37" in str(error)
    assert "TASK18-STDOUT" in str(error)
    assert "TASK18-STDERR" in str(error)
    assert not (tmp_path / "merged.fastq").exists()


@pytest.mark.integration
def test_fastp_paired_merge_conserves_records_and_order(paired_fastq_inputs, tmp_path):
    _, forward, reverse = paired_fastq_inputs
    output_directory = tmp_path / "paired-output"
    outputs = pp.merge_fastqs(
        [str(forward), str(reverse)],
        str(output_directory),
        binary_path=get_binary_path("fastp"),
        minimum_overlap=30,
        trim_adapters=False,
        quality_trim=False,
    )

    assert outputs == [str(output_directory / "sample.fastq")]
    merged_reads = list(parse_fastx(outputs[0]))
    assert len(merged_reads) == 2
    assert [read.id for read in merged_reads] == [f"{identifier}/1" for identifier in PAIR_IDS]
    assert [read.sequence for read in merged_reads] == list(PAIR_SEQUENCES)
    assert set(output_directory.iterdir()) == {output_directory / "sample.fastq"}


@pytest.mark.integration
def test_fastp_interleaved_merge_conserves_records_and_cleans_temporary_input(
    interleaved_fastq_input, tmp_path,
):
    output_directory = tmp_path / "interleaved-output"
    outputs = pp.merge_fastqs(
        [str(interleaved_fastq_input)],
        str(output_directory),
        interleaved=True,
        binary_path=get_binary_path("fastp"),
        minimum_overlap=30,
        trim_adapters=False,
        quality_trim=False,
    )

    assert outputs == [str(output_directory / "interleaved.fastq")]
    merged_reads = list(parse_fastx(outputs[0]))
    assert len(merged_reads) == 2
    assert [read.id for read in merged_reads] == [f"{identifier}/1" for identifier in PAIR_IDS]
    assert [read.sequence for read in merged_reads] == list(PAIR_SEQUENCES)
    assert set(output_directory.iterdir()) == {output_directory / "interleaved.fastq"}


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
        with open(forward) as temporary:
            observed["content"] = temporary.read()
        observed["merged"] = merged

    monkeypatch.setattr(merging, "merge_fastqs_fastp", fake_merge)

    outputs = merging.merge_fastqs(
        [str(interleaved)], str(tmp_path / "output"), interleaved=True
    )

    assert observed["content"].count("@read") == 2
    assert not os.path.exists(observed["temporary"])
    assert outputs == [str(tmp_path / "output" / "iraq.fastq")]


@pytest.mark.e2e
def test_run_reports_fastp_failure_as_structured_preprocess_error(
    monkeypatch, paired_fastq_inputs, tmp_path,
):
    input_directory, _, _ = paired_fastq_inputs
    project = tmp_path / "failed-project"

    def fail(command, **kwargs):
        assert isinstance(command, list)
        assert kwargs == {"check": True, "capture_output": True, "text": True}
        Path(command[command.index("--merged_out") + 1]).write_text("partial output")
        raise subprocess.CalledProcessError(
            41, command, output="PUBLIC-STDOUT", stderr="PUBLIC-STDERR"
        )

    monkeypatch.setattr(merging.sp, "run", fail)
    with pytest.raises(AnnotationRunError) as captured:
        run(
            str(input_directory),
            project_path=str(project),
            merge=True,
            merge_kwargs={"merge_args": ["--report_title", "x; $(false)"]},
            n_processes=1,
        )

    assert len(captured.value.failures) == 1
    failure = captured.value.failures[0]
    assert (failure.stage, failure.category) == ("preprocess", "external_tool")
    assert "EXIT STATUS: 41" in failure.message
    assert "PUBLIC-STDOUT" in failure.message
    assert "PUBLIC-STDERR" in failure.message
    assert captured.value.partial_output_paths == ()
    assert not list((project / "merged").glob("*.fastq"))
    failure_log = project / "logs" / "merge_fastqs.failed"
    assert failure_log.is_file()
    assert "EXIT STATUS: 41" in failure_log.read_text()
    assert "PUBLIC-STDOUT" in failure_log.read_text()
    assert "PUBLIC-STDERR" in failure_log.read_text()


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
