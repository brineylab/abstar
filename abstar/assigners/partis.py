#!/usr/bin/env python
# filename: partis.py

#
# Copyright (c) 2016 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

import csv
import os
import platform
import subprocess as sp
from tempfile import NamedTemporaryFile
import traceback

from abutils.core.sequence import Sequence
from abutils.utils.alignment import local_alignment

from .assigner import BaseAssigner
from ..core.germline import GermlineSegment
from ..core.vdj import VDJ



class Partis(BaseAssigner):
    """
    docstring for Partis
    """

    def __init__(self, species):
        super(Partis, self).__init__(species)


    def __call__(self, sequence_file, file_format, locus='ig'):
        self.assigned, self.unassigned = run_partis(sequence_file, file_format, locus)


    def run_partis(self, sequence_file, file_format, locus):
        with open(sequence_file, 'r') as sequence_handle:
            seqs = [Sequence(s) for s in SeqIO.parse(sequence_handle, file_format)]
        seq_dict = {s.id.replace(':', 'c'): s for s in seqs}
        partis_out = NamedTemporaryFile(delete=False)
        germline_dir = os.path.join(self.germline_directory, 'partis/')
        locus = 'igh' if locus == 'ig' else 'tra'
        partis_cmd = ['partis', 'run-viterbi',
                      '--infname', sequence_file,
                      '--outfname', partis_out.name,
                      '--locus', locus,
                      '--initial-germline-dir', germline_dir]
        p = sp.Popen(partis_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        parsed_vdjs = self.parse_partis_output(partis_out.name, seq_dict)
        os.unlink(partis_out.name)
        return parsed_vdjs


    def parse_partis_output(self, output_file, seq_dict):
        assigned = []
        unassigned = []
        with open(output_file, 'rb') as f:
            csv_dict = list(csv.DictReader(f))
        for c in csv_dict:
            partis_id = c['unique_ids']
            seq = seq_dict[partis_id]
            seq_id = seq.id
            v = c['v_gene']
            if v[:3].upper() in ['IGH', 'TRA', 'TRD']:
                d = c['d_gene'] if c['d_gene'].strip() else None
            else:
                d = None
            j = c['j_gene']
            vdjs.append(VDJ(seq, v=v, d=d, j=j))
            if c['duplicates'].strip():
                for dup_id in c['duplicates'].strip().split(':'):
                    seq = seq_dict[dup_id]
                    vdjs.append(VDJ(seq, v=v, d=d, j=j))
        return vdjs


    def compute_alignment_positions(self):
        pass
