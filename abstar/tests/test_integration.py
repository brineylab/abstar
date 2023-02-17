#!usr/env/python
# filename: test_integration.py

import os
import pytest
import sys
import tempfile

from Bio import SeqIO

from abutils.core.sequence import Sequence

from ..core import abstar


def test_api_integration_bcr_1k():
    test_data = os.path.abspath('abstar/test_data/test_1k.fasta')
    with open(test_data, 'r') as f:
        test_seqs = [Sequence(s) for s in SeqIO.parse(f, 'fasta')]
    seqs = abstar.run(*test_seqs)
    assert len(seqs) == 999


def test_api_integration_bcr_hiv_bnab_hcs():
    test_data = os.path.abspath('abstar/test_data/test_hiv_bnab_hcs.fasta')
    with open(test_data, 'r') as f:
        test_seqs = [Sequence(s) for s in SeqIO.parse(f, 'fasta')]
    seqs = abstar.run(*test_seqs)
    assert len(seqs) == 221


def test_api_integration_bcr_hiv_bnab_lcs():
    test_data = os.path.abspath('abstar/test_data/test_hiv_bnab_lcs.fasta')
    with open(test_data, 'r') as f:
        test_seqs = [Sequence(s) for s in SeqIO.parse(f, 'fasta')]
    seqs = abstar.run(*test_seqs)
    assert len(seqs) == 207


@pytest.mark.parametrize("num_cores", [1, 4])
def test_chunks(num_cores):
    """Test that a small file is correctly split into chunks, and returns correct set of json files"""
    temp_dir = tempfile.mkdtemp()
    out_dir = tempfile.mkdtemp()
    arg_list = ["-i", "abstar/test_data/test.fastq", "-o", out_dir, "-t", temp_dir, "-r", "tcr", "--num-cores", str(num_cores), "--chunksize", "6"]
    abstar.run_main(arg_list)
    assert len(os.listdir(os.path.join(temp_dir, "input"))) == len(os.listdir(os.path.join(temp_dir, "json")))
