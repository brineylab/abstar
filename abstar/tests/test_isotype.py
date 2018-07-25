#!usr/env/python
# filename: test_integration.py

import os
import sys

from Bio import SeqIO

from abutils.core.sequence import Sequence

from .mock_classes import MockAntibody
from ..core import abstar
from ..core.antibody import Antibody
from ..utils.isotype import get_isotype


def test_isotype_bcr():
    test_fasta = os.path.abspath('abstar/test_data/test_isotype.fasta')
    with open(test_fasta, 'r') as f:
        test_seqs = [Sequence(s) for s in SeqIO.parse(f, 'fasta')]
    vdj_nt = [s for s in test_seqs if s.id == 'vdj_nt'][0].sequence
    oriented_input = [s for s in test_seqs if s.id == 'oriented_input'][0].sequence
    antibody = MockAntibody()
    antibody.species = 'human'
    antibody.vdj_nt = vdj_nt
    antibody.oriented_input = oriented_input
    isotype = get_isotype(antibody)
    assert isotype.isotype == 'IgG1'
