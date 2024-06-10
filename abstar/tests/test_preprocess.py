import os
import tempfile
from subprocess import PIPE, Popen

import pytest
from Bio import SeqIO

from ..utils.preprocess import adapter_trim, fastqc, quality_trim


@pytest.fixture()
def tmp_path():
    return tempfile.mkdtemp()


@pytest.fixture()
def test_fastq(tmp_path):
    tmp_path = tempfile.mkdtemp()
    test_file = os.path.join(tmp_path, "test.fastq")
    with open(test_file, "w") as f:
        f.write("@test\nACGTCGCGTATA\n+\nIIIIIIIIIIII\n")
    return tmp_path


@pytest.fixture()
def adapter_file():
    tmp_path = tempfile.mkdtemp()
    adapter_file = os.path.join(tmp_path, "adapter.fasta")
    with open(adapter_file, "w") as f:
        f.write(">adapter\nACGT\n")
    return adapter_file


@pytest.fixture
# def input_directory(tmpdir):
def input_directory():
    # create temporary input directory
    input_dir = tempfile.mkdtemp()
    # create temporary input files
    sequence = "ACGT" * 100
    qual1 = "IIII" * 100
    qual2 = "IIII" * 90 + "!!!!" * 10
    input_file1 = os.path.join(input_dir, "input_file1.fastq")
    with open(input_file1, "w") as f:
        f.write(f"@seq1\n{sequence}\n+\n{qual1}\n")
    input_file2 = os.path.join(input_dir, "input_file2.fastq")
    with open(input_file2, "w") as f:
        f.write(f"@seq1\n{sequence}\n+\n{qual2}\n")
    # return path to input directory
    return str(input_dir)


@pytest.fixture
def output_directory():
    # create temporary output directory
    output_dir = tempfile.mkdtemp()
    # return path to output directory
    return str(output_dir)


# ----------------------------
#      quality trim
# ----------------------------


def test_quality_trim_compressed(input_directory, output_directory):
    # run quality_trim function
    quality_trim(
        input_directory=input_directory,
        output_directory=output_directory,
        quality_cutoff=20,
        length_cutoff=50,
        quality_type="sanger",
        compress_output=True,
        file_pairs=None,
        singles_directory=None,
        # nextseq=False,
        paired_reads=False,
        allow_5prime_trimming=False,
        print_debug=False,
    )
    # check that input files were created
    input_files = os.listdir(input_directory)
    assert len(input_files) == 2
    # check that output files were created
    output_files = os.listdir(output_directory)
    assert len(output_files) == 2
    assert os.path.exists(os.path.join(output_directory, "input_file1.fastq.gz"))
    assert os.path.exists(os.path.join(output_directory, "input_file2.fastq.gz"))


def test_quality_trim_uncompressed(input_directory, output_directory):
    # run quality_trim function
    quality_trim(
        input_directory=input_directory,
        output_directory=output_directory,
        quality_cutoff=20,
        length_cutoff=50,
        quality_type="sanger",
        compress_output=False,
        file_pairs=None,
        singles_directory=None,
        # nextseq=False,
        paired_reads=False,
        allow_5prime_trimming=False,
        print_debug=False,
    )
    # check that input files were created
    input_files = os.listdir(input_directory)
    assert len(input_files) == 2
    # check that output files were created
    output_files = os.listdir(output_directory)
    assert len(output_files) == 2
    assert os.path.exists(os.path.join(output_directory, "input_file1.fastq"))
    assert os.path.exists(os.path.join(output_directory, "input_file2.fastq"))


# ----------------------------
#       adapter trim
# ----------------------------


def test_adapter_trim(test_fastq, adapter_file):
    # run adapter_trim
    output_dir = adapter_trim(input_directory=test_fastq, adapter_both=adapter_file)
    # check that output files were created
    assert os.path.exists(os.path.join(output_dir, "test.fastq"))


def test_adapter_trim_compressed(test_fastq, adapter_file):
    # run adapter_trim
    output_dir = adapter_trim(
        input_directory=test_fastq,
        adapter_both=adapter_file,
        compress_output=True,
    )
    # check that output files were created
    assert os.path.exists(os.path.join(output_dir, "test.fastq.gz"))


def test_adapter_trim_output_directory(test_fastq, adapter_file):
    # run adapter_trim with a custom output directory
    custom_output_dir = tempfile.mkdtemp()
    output_dir = adapter_trim(
        input_directory=test_fastq,
        output_directory=custom_output_dir,
        adapter_both=adapter_file,
    )
    # check that output files were created in the custom output directory
    assert output_dir == custom_output_dir
    assert os.path.exists(os.path.join(custom_output_dir, "test.fastq"))


def test_adapter_trim_5prime(test_path, adapter_file):
    # run adapter_trim with a 5' adapter
    output_dir = adapter_trim(str(test_path), adapter_5prime=str(adapter_file))
    # check that output files were created
    assert os.path.exists(os.path.join(output_dir, "test.fastq"))


def test_adapter_trim_3prime(test_path, adapter_file):
    # run adapter_trim with a 3' adapter
    output_dir = adapter_trim(str(test_path), adapter_3prime=str(adapter_file))
    # check that output files were created
    assert os.path.exists(os.path.join(output_dir, "test.fastq"))


def test_adapter_trim_both(test_path, adapter_file):
    # run adapter_trim with a 5' and 3' adapter
    output_dir = adapter_trim(str(test_path), adapter_both=str(adapter_file))
    # check that output files were created
    assert os.path.exists(os.path.join(output_dir, "test.fastq"))


def test_adapter_trim_5prime_anchored(test_path, adapter_file):
    # run adapter_trim with a 5' anchored adapter
    output_dir = adapter_trim(str(test_path), adapter_5prime_anchored=str(adapter_file))
    # check that output files were created
    assert os.path.exists(os.path.join(output_dir, "test.fastq"))


def test_adapter_trim_3prime_anchored(test_path, adapter_file):
    # run adapter_trim with a 3' anchored adapter
    output_dir = adapter_trim(str(test_path), adapter_3prime_anchored=str(adapter_file))
    # check that output files were created
    assert os.path.exists(os.path.join(output_dir, "test.fastq"))


def test_adapter_trim_no_input_files(tmp_path):
    # run adapter_trim on an empty input directory
    with pytest.raises(ValueError):
        adapter_trim(str(tmp_path))


def test_adapter_trim_invalid_input_directory(tmp_path):
    # run adapter_trim on a non-existent input directory
    with pytest.raises(ValueError):
        adapter_trim("nonexistent_directory", str(tmp_path))


def test_adapter_trim_invalid_threads(test_path):
    # run adapter_trim with an invalid number of threads
    with pytest.raises(ValueError):
        adapter_trim(str(test_path), threads=-2)


def test_adapter_trim_command_line(test_path):
    # run adapter_trim using the command line
    adapter_string = "-g {}".format(str(SeqIO.read(str(adapter_file), "fasta").seq))
    cutadapt_cmd = "cutadapt -o {}/test.fastq.gz {} {}".format(
        str(test_path), adapter_string, str(test_path)
    )
    p = Popen(cutadapt_cmd, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = p.communicate()
    # run adapter_trim using the function
    output_dir = adapter_trim(str(test_path))
    # check that the output is the same
    assert stdout == b""
    assert stderr == b""
    assert os.path.exists(os.path.join(output_dir, "test.fastq.gz"))


# ----------------------------
#          fastqc
# ----------------------------


def test_fastqc(tmp_path):
    # create a test file
    test_file = tmp_path / "test.fastq"
    with open(test_file, "w") as f:
        f.write("@test\nACGT\n+\nIIII\n")
    # run fastqc
    output_dir = fastqc(str(tmp_path))
    # check that output files were created
    assert os.path.exists(os.path.join(output_dir, "test_fastqc.html"))
    assert os.path.exists(os.path.join(output_dir, "test_fastqc.zip"))


def test_fastqc_threads(tmp_path):
    # create a test file
    test_file = tmp_path / "test.fastq"
    with open(test_file, "w") as f:
        f.write("@test\nACGT\n+\nIIII\n")
    # run fastqc with 2 threads
    output_dir = fastqc(str(tmp_path), threads=2)
    # check that output files were created
    assert os.path.exists(os.path.join(output_dir, "test_fastqc.html"))
    assert os.path.exists(os.path.join(output_dir, "test_fastqc.zip"))


def test_fastqc_output_directory(tmp_path):
    # create a test file
    test_file = tmp_path / "test.fastq"
    with open(test_file, "w") as f:
        f.write("@test\nACGT\n+\nIIII\n")
    # run fastqc with a custom output directory
    custom_output_dir = tmp_path / "custom_output"
    output_dir = fastqc(str(tmp_path), str(custom_output_dir))
    # check that output files were created in the custom output directory
    assert os.path.exists(os.path.join(output_dir, "test_fastqc.html"))
    assert os.path.exists(os.path.join(output_dir, "test_fastqc.zip"))


def test_fastqc_no_input_files(tmp_path):
    # run fastqc on an empty input directory
    with pytest.raises(ValueError):
        fastqc(str(tmp_path))


def test_fastqc_invalid_input_directory(tmp_path):
    # run fastqc on a non-existent input directory
    with pytest.raises(ValueError):
        fastqc("nonexistent_directory", str(tmp_path))


def test_fastqc_invalid_threads(tmp_path):
    # run fastqc with an invalid number of threads
    with pytest.raises(ValueError):
        fastqc(str(tmp_path), threads=-2)


def test_fastqc_command_line(tmp_path):
    # create a test file
    test_file = tmp_path / "test.fastq"
    with open(test_file, "w") as f:
        f.write("@test\nACGT\n+\nIIII\n")
    # run fastqc using the command line
    fastqc_cmd = "fastqc --noextract -o={} {}".format(str(tmp_path), test_file)
    p = Popen(fastqc_cmd, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = p.communicate()
    # run fastqc using the function
    output_dir = fastqc(str(tmp_path))
    # check that the output is the same
    assert stdout == b""
    assert stderr == b""
    assert os.path.exists(os.path.join(output_dir, "test_fastqc.html"))
    assert os.path.exists(os.path.join(output_dir, "test_fastqc.zip"))
