#!usr/env/python
# filename: test_integration.py

import os
import sys

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
