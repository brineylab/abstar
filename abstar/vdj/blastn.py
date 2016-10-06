#!/usr/bin/env python
# filename: blastn.py

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

import os
import platform
from tempfile import NamedTemporaryFile

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

from .assigner import AbstractAssigner
from .germline import Germline
from .vdj import VDJ


class BlastnAssigner(AbstractAssigner):
    """
    docstring for BlastnAssigner
    """

    def __init__(self):
        super(BlastnAssigner, self).__init__()


    def __call__(self, sequence_file, species):
        seqs = [Sequence(s) for s in SeqIO.parse(open(sequence_file, 'r'), 'fasta')]
        vdjs = []

        # assign V-genes
        vblast_records = self.blast(sequence_file, species, 'V')
        if not vblast_records:
            # TODO: log that there are no valid antibody sequences in the input file
            self.unassigned = seqs
            return
        jquery_seqs = []
        for seq, vbr in zip(seqs, vblast_records):
            try:
                germ = self.process_blast_record(vbr, species)
                vdj = VDJ(seq, v=germ)
                self.orient_query(vdj, vbr)
                jquery = self.get_jquery_sequence(vdj.oriented, vbr)
                if len(jquery) >= 10:
                    vdjs.append(vdj)
                    jquery_seqs.append(jquery)
                else:
                    # TODO: log that the jquery sequence was too short
                    self.unassigned.append(seq)
            except:
                # TODO: log the exception
                self.unassigned.append(seq)

        # assign J-genes
        _vdjs = []
        dquery_seqs = []
        jblast_infile = self.build_jblast_input(jquery_seqs)
        jblast_records = self.blast(jblast_infile, species, 'J')
        for vdj, jquery, jbr in zip(vdjs, jquery_seqs, jblast_records):
            germ = self.process_blast_record(jbr, species)
            vdj.j = germ
            if vdj.v.chain != vdj.j.chain:
                # TODO: log that the V and J chains don't match
                self.unassigned.append(vdj.sequence)
                continue
            dquery = self.get_dquery_sequence(jquery, jbr)
            dquery_seqs.append(dquery)
            _vdjs.append(vdj)
        os.unlink(jblast_infile)
        vdjs = _vdjs

        # assign D-genes
        _vdjs = []
        for vdj, dquery in zip(vdjs, dquery_seqs):
            if vdj.v.chain == 'heavy':
                germ = self.assign_dgene(dquery, species)
                vdj.d = germ
            _vdjs.append(vdj)
        self.assigned = _vdjs


    def blast(self, seq_file, species, segment):
        '''
        Runs BLASTn against an antibody germline database.

        Args:
        -----

            seq_file (str): Path to a FASTA-formatted file of input sequences.

            species (str): Species of origin of the antibody sequences in ``seq_file``.
                Options are: ``human``, ``macaque``, ``mouse`` and ``rabbit``.

            segment (str): Germline segment to query. Options are ``V`` and ``J``.
        '''
        blast_path = os.path.join(self.binary_directory, 'blastn_{}'.format(platform.system().lower()))
        blast_db_path = os.path.join(self.germline_directory, 'blast/{}_gl_{}'.format(species.lower(), segment.upper()))
        blastout = NamedTemporaryFile(delete=False)
        blastn_cmd = NcbiblastnCommandline(cmd=blast_path,
                                           db=blast_db_path,
                                           query=seq_file,
                                           out=blastout.name,
                                           outfmt=5,
                                           dust='no',
                                           word_size=self._word_size(segment),
                                           max_target_seqs=10,
                                           evalue=self._evalue(segment),
                                           reward=self._match_reward(segment),
                                           penalty=self._mismatch_penalty(segment),
                                           gapopen=self._gap_open(segment),
                                           gapextend=self._gap_extend(segment))
        stdout, stderr = blastn_cmd()
        blast_records = [br for br in NCBIXML.parse(blastout)]
        os.unlink(blastout.name)
        return blast_records


    def assign_dgene(self, seq, species):
        db_file = os.path.join(self.germline_directory, '/fasta/{}_D.fasta'.format(species.lower()))
        db_handle = open(db_file, 'r')
        germs = [Sequence(s) for s in SeqIO.parse(db_handle, 'fasta')]
        rc_germs = [Sequence(s.reverse_complement, id=s.id) for s in germs]
        germs.extend(rc_germs)
        db_handle.close()
        alignments = local_alignment(seq, targets=germs,
                                     gap_open=-20, gap_extend=-2)
        alignments.sort(key=lambda x: x.score, reverse=True)
        all_gls = [a.target.id for a in alignments]
        all_scores = [a.score for a in alignments]
        top_gl = all_gls[0]
        top_score = all_scores[0]
        others = [Germline(germ, species, score=score) for germ, score in zip(all_gls[1:6], all_scores[1:6])]
        return Germline(top_gl, species, score=top_score, others=others)


    @staticmethod
    def process_blast_record(blast_record, species):
        all_gls = [a.title.split()[0] for a in blast_record.alignments]
        all_scores = [a.hsps[0].bits for a in blast_record.alignments]
        top_gl = all_gls[0]
        top_score = all_scores[0]
        others = [Germline(germ, species, score=score) for germ, score in zip(all_gls[1:], all_scores[1:])]
        return Germline(top_gl, species, score=top_score, others=others[:5])


    @staticmethod
    def orient_query(vdj, vbr):
        hsp = vbr.alignments[0].hsps[0]
        if hsp.sbjct_start > hsp.sbjct_end:
            vdj.oriented = Sequence(vdj.sequence.reverse_complement, id=vdj.sequence.id)


    @staticmethod
    def get_jquery_sequence(seq, vbr):
        query_end = vbr.alignments[0].hsps[0].query_end
        return Sequence(seq[query_end:], id=seq.id)


    @staticmethod
    def get_dquery_sequence(seq, jbr):
        query_start = jbr.alignments[0].hsps[0].query_start - 1
        return Sequence(seq[:query_start], id=seq.id)


    @staticmethod
    def build_jblast_input(jseqs):
        jblast_input = NamedTemporaryFile(delete=False)
        jblast_input.write('\n'.join([s.fasta for s in jseqs]))
        jblast_input.close()
        return jblast_input.name


    @staticmethod
    def _word_size(segment):
        word_sizes = {'V': 11,
                      'D': 4,
                      'J': 7}
        return word_sizes[segment]

    @staticmethod
    def _gap_open(segment):
        gap_open = {'V': 5,
                    'D': 4,
                    'J': 5}
        return gap_open[segment]

    @staticmethod
    def _gap_extend(segment):
        gap_extend = {'V': 2,
                      'D': 2,
                      'J': 2}
        return gap_extend[segment]

    @staticmethod
    def _match_reward(segment):
        match = {'V': 1,
                 'D': 1,
                 'J': 1}
        return match[segment]

    @staticmethod
    def _mismatch_penalty(segment):
        mismatch = {'V': -1,
                    'D': -1,
                    'J': -1}
        return mismatch[segment]

    @staticmethod
    def _evalue(segment):
        evalue = {'V': 1,
                  'D': 100000,
                  'J': 1000}
        return evalue[segment]
