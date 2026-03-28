# Copyright (c) 2025 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT

"""
Tests for the preprocessing/merging module.

Tests filename parsing and grouping logic. Does not test actual merging
(requires fastp binary).
"""

import pytest

from abstar.preprocess.merging import (
    ElementFile,
    FASTQFile,
    IlluminaFile,
    MergeGroup,
    group_paired_fastqs,
)

# =============================================
#          FASTQFILE BASE CLASS TESTS
# =============================================


def test_fastqfile_basic_properties(tmp_path):
    """Test FASTQFile parses path components."""
    fq = tmp_path / "sample_S1_L001_R1_001.fastq"
    fq.touch()
    f = FASTQFile(fq)
    assert f.basename == "sample_S1_L001_R1_001.fastq"
    assert f.dir == tmp_path
    assert f.filename == "sample_S1_L001_R1_001"


def test_fastqfile_gz_suffix(tmp_path):
    """Test FASTQFile strips .gz and .fastq suffixes."""
    fq = tmp_path / "sample_S1_L001_R1_001.fastq.gz"
    fq.touch()
    f = FASTQFile(fq)
    assert f.filename == "sample_S1_L001_R1_001"


def test_fastqfile_fq_suffix(tmp_path):
    """Test FASTQFile handles .fq extension."""
    fq = tmp_path / "sample_S1_L001_R1_001.fq"
    fq.touch()
    f = FASTQFile(fq)
    assert f.filename == "sample_S1_L001_R1_001"


def test_fastqfile_fq_gz_suffix(tmp_path):
    """Test FASTQFile handles .fq.gz extension."""
    fq = tmp_path / "sample_S1_L001_R1_001.fq.gz"
    fq.touch()
    f = FASTQFile(fq)
    assert f.filename == "sample_S1_L001_R1_001"


# =============================================
#       ILLUMINA FILE PARSING TESTS
# =============================================


def test_illumina_file_name(tmp_path):
    """Test IlluminaFile extracts sample name."""
    fq = tmp_path / "MySample_S1_L001_R1_001.fastq"
    fq.touch()
    f = IlluminaFile(fq)
    assert f.name == "MySample"
    assert f.schema == "illumina"


def test_illumina_file_read(tmp_path):
    """Test IlluminaFile extracts read number."""
    r1 = tmp_path / "MySample_S1_L001_R1_001.fastq"
    r2 = tmp_path / "MySample_S1_L001_R2_001.fastq"
    r1.touch()
    r2.touch()
    assert IlluminaFile(r1).read == "R1"
    assert IlluminaFile(r2).read == "R2"


def test_illumina_file_lane(tmp_path):
    """Test IlluminaFile extracts lane."""
    fq = tmp_path / "MySample_S1_L001_R1_001.fastq"
    fq.touch()
    f = IlluminaFile(fq)
    assert f.lane == "L001"


def test_illumina_file_sample(tmp_path):
    """Test IlluminaFile extracts sample index."""
    fq = tmp_path / "MySample_S1_L001_R1_001.fastq"
    fq.touch()
    f = IlluminaFile(fq)
    assert f.sample == "S1"


def test_illumina_file_number(tmp_path):
    """Test IlluminaFile extracts file number."""
    fq = tmp_path / "MySample_S1_L001_R1_001.fastq"
    fq.touch()
    f = IlluminaFile(fq)
    assert f.number == "001"


def test_illumina_file_with_underscores_in_name(tmp_path):
    """Test IlluminaFile handles sample names with underscores."""
    fq = tmp_path / "My_Long_Sample_Name_S1_L001_R1_001.fastq"
    fq.touch()
    f = IlluminaFile(fq)
    assert f.name == "My_Long_Sample_Name"
    assert f.sample == "S1"
    assert f.lane == "L001"
    assert f.read == "R1"


def test_illumina_file_gz(tmp_path):
    """Test IlluminaFile works with .fastq.gz files."""
    fq = tmp_path / "MySample_S1_L001_R1_001.fastq.gz"
    fq.touch()
    f = IlluminaFile(fq)
    assert f.name == "MySample"
    assert f.read == "R1"
    assert f.lane == "L001"


def test_illumina_file_equality(tmp_path):
    """Test IlluminaFile equality is based on name."""
    r1 = tmp_path / "MySample_S1_L001_R1_001.fastq"
    r2 = tmp_path / "MySample_S1_L001_R2_001.fastq"
    r1.touch()
    r2.touch()
    assert IlluminaFile(r1) == IlluminaFile(r2)


# =============================================
#        ELEMENT FILE PARSING TESTS
# =============================================


def test_element_file_name(tmp_path):
    """Test ElementFile extracts sample name."""
    fq = tmp_path / "MySample_R1.fastq"
    fq.touch()
    f = ElementFile(fq)
    assert f.name == "MySample"
    assert f.schema == "element"


def test_element_file_read(tmp_path):
    """Test ElementFile extracts read number."""
    r1 = tmp_path / "MySample_R1.fastq"
    r2 = tmp_path / "MySample_R2.fastq"
    r1.touch()
    r2.touch()
    assert ElementFile(r1).read == "R1"
    assert ElementFile(r2).read == "R2"


def test_element_file_empty_properties(tmp_path):
    """Test ElementFile returns empty strings for lane/number/sample."""
    fq = tmp_path / "MySample_R1.fastq"
    fq.touch()
    f = ElementFile(fq)
    assert f.lane == ""
    assert f.number == ""
    assert f.sample == ""


def test_element_file_with_underscores_in_name(tmp_path):
    """Test ElementFile handles sample names with underscores."""
    fq = tmp_path / "My_Long_Sample_R1.fastq"
    fq.touch()
    f = ElementFile(fq)
    assert f.name == "My_Long_Sample"
    assert f.read == "R1"


def test_element_file_gz(tmp_path):
    """Test ElementFile works with .fastq.gz files."""
    fq = tmp_path / "MySample_R1.fastq.gz"
    fq.touch()
    f = ElementFile(fq)
    assert f.name == "MySample"
    assert f.read == "R1"


def test_element_file_equality(tmp_path):
    """Test ElementFile equality is based on name."""
    r1 = tmp_path / "MySample_R1.fastq"
    r2 = tmp_path / "MySample_R2.fastq"
    r1.touch()
    r2.touch()
    assert ElementFile(r1) == ElementFile(r2)


# =============================================
#       GROUP PAIRED FASTQS TESTS
# =============================================


def test_group_illumina_single_sample(tmp_path):
    """Test grouping a single Illumina sample pair."""
    (tmp_path / "Sample_S1_L001_R1_001.fastq").touch()
    (tmp_path / "Sample_S1_L001_R2_001.fastq").touch()

    groups = group_paired_fastqs(str(tmp_path), schema="illumina")
    assert len(groups) == 1
    assert groups[0].name == "Sample"
    assert len(groups[0].files) == 2


def test_group_illumina_multiple_samples(tmp_path):
    """Test grouping multiple Illumina samples."""
    for sample in ["SampleA", "SampleB"]:
        (tmp_path / f"{sample}_S1_L001_R1_001.fastq").touch()
        (tmp_path / f"{sample}_S1_L001_R2_001.fastq").touch()

    groups = group_paired_fastqs(str(tmp_path), schema="illumina")
    assert len(groups) == 2
    names = {g.name for g in groups}
    assert names == {"SampleA", "SampleB"}


def test_group_illumina_multiple_lanes(tmp_path):
    """Test grouping Illumina sample with multiple lanes."""
    for lane in ["L001", "L002"]:
        (tmp_path / f"Sample_S1_{lane}_R1_001.fastq").touch()
        (tmp_path / f"Sample_S1_{lane}_R2_001.fastq").touch()

    groups = group_paired_fastqs(str(tmp_path), schema="illumina")
    assert len(groups) == 1
    assert len(groups[0].files) == 4


def test_group_illumina_filters_index_reads(tmp_path):
    """Test that index reads (I1/I2) are filtered out."""
    (tmp_path / "Sample_S1_L001_R1_001.fastq").touch()
    (tmp_path / "Sample_S1_L001_R2_001.fastq").touch()
    (tmp_path / "Sample_S1_L001_I1_001.fastq").touch()
    (tmp_path / "Sample_S1_L001_I2_001.fastq").touch()

    groups = group_paired_fastqs(str(tmp_path), schema="illumina")
    assert len(groups) == 1
    assert len(groups[0].files) == 2
    reads = {f.read for f in groups[0].files}
    assert reads == {"R1", "R2"}


def test_group_element_single_sample(tmp_path):
    """Test grouping a single Element sample pair."""
    (tmp_path / "Sample_R1.fastq").touch()
    (tmp_path / "Sample_R2.fastq").touch()

    groups = group_paired_fastqs(str(tmp_path), schema="element")
    assert len(groups) == 1
    assert groups[0].name == "Sample"
    assert len(groups[0].files) == 2


def test_group_invalid_schema(tmp_path):
    """Test that invalid schema raises ValueError."""
    (tmp_path / "Sample_R1.fastq").touch()
    with pytest.raises(ValueError, match="Invalid schema"):
        group_paired_fastqs(str(tmp_path), schema="invalid")


def test_group_nonexistent_directory():
    """Test that nonexistent directory raises ValueError."""
    with pytest.raises(ValueError):
        group_paired_fastqs("/nonexistent/path", schema="illumina")


def test_group_file_list(tmp_path):
    """Test grouping from an explicit list of file paths."""
    r1 = tmp_path / "Sample_S1_L001_R1_001.fastq"
    r2 = tmp_path / "Sample_S1_L001_R2_001.fastq"
    r1.touch()
    r2.touch()

    groups = group_paired_fastqs([str(r1), str(r2)], schema="illumina")
    assert len(groups) == 1
    assert groups[0].name == "Sample"


# =============================================
#         MERGE GROUP TESTS
# =============================================


def test_merge_group_lane_grouping(tmp_path):
    """Test MergeGroup groups files by lane internally."""
    files = []
    for lane in ["L001", "L002"]:
        for read in ["R1", "R2"]:
            fq = tmp_path / f"Sample_S1_{lane}_{read}_001.fastq"
            fq.touch()
            files.append(IlluminaFile(fq))

    group = MergeGroup("Sample", files)
    lanes = group._group_by_lane()
    assert len(lanes) == 2
    assert len(lanes[0]) == 2
    assert len(lanes[1]) == 2
