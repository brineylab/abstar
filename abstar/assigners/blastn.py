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
import traceback

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

from abutils.core.sequence import Sequence
from abutils.utils.alignment import local_alignment

from .assigner import BaseAssigner
from ..core.germline import GermlineSegment
from ..core.vdj import VDJ



class Blastn(BaseAssigner):
    """
    docstring for Blastn
    """

    def __init__(self, species):
        super(Blastn, self).__init__(species)


    def __call__(self, sequence_file, file_format):
        with open(sequence_file, 'r') as sequence_handle:
            seqs = [Sequence(s) for s in SeqIO.parse(sequence_handle, file_format)]
        vdjs = []

        # if the input file is FASTQ-formatted, need to convert it to FASTA for BLASTn to work
        if file_format == 'fastq':
            with open(sequence_file, 'w') as handle:
                handle.write('\n'.join(s.fasta for s in seqs))

        # assign V-genes
        vblast_records = self.blast(sequence_file, self.species, 'V')
        # if there aren't any vblast_records, that means that none of the
        # sequences in the input file contained sequences with a significant
        # match to any germline V-gene. These are likely all non-antibody sequences.
        if not vblast_records:
            vdjs = [VDJ(seq) for seq in seqs]
            for vdj in vdjs:
                vdj.log('V-GENE ASSIGNMENT ERROR:',
                        'No variable gene was found.',
                        'Query sequence does not appear to contain a rearranged antibody.')
            self.unassigned = vdjs
            return
        jquery_seqs = []
        for seq, vbr in zip(seqs, vblast_records):
            try:
                germ = self.process_blast_record(vbr, self.species)
                vdj = VDJ(seq, v=germ)
                self.orient_query(vdj, vbr)
                jquery = self.get_jquery_sequence(vdj.oriented, vbr)
                # only try to find J-genes if there's a minimum of 10 nucleotides
                # remaining after removal of the V-gene alignment
                if len(jquery) >= 10:
                    vdjs.append(vdj)
                    jquery_seqs.append(jquery)
                # abort VDJ assignment if the J-gene query sequence is too short
                else:
                    vdj = VDJ(seq)
                    vdj.log('J-GENE QUERY ERROR:', 'Query sequence for J-gene assignment is too short.')
                    vdj.log('')
                    query = vbr.alignments[0].hsps[0].query
                    subject = vbr.alignments[0].hsps[0].sbjct
                    vdj.log(' QUERY :', query)
                    vdj.log('        ', ''.join(['|' if q == s else ' ' for q, s in zip(query, subject)]))
                    vdj.log('SUBJECT:', subject)
                    vdj.log('')
                    vdj.log('J-QUERY SEQUENCE:', jquery.sequence)
                    self.unassigned.append(vdj)
            except:
                vdj = VDJ(seq)
                vdj.exception('V-GENE ASSIGNMENT ERROR', traceback.format_exc())
                self.unassigned.append(vdj)

        # assign J-genes
        _vdjs = []
        dquery_seqs = []
        jblast_infile = self.build_jblast_input(jquery_seqs)
        jblast_records = self.blast(jblast_infile, self.species, 'J')
        for vdj, jquery, jbr in zip(vdjs, jquery_seqs, jblast_records):
            try:
                germ = self.process_blast_record(jbr, self.species)
                vdj.j = germ
                # sanity check to make sure there's not an obvious problem with the V/J
                # assignments (likely due to poor germline matches to a non-antibody sequence)
                if vdj.v.chain != vdj.j.chain:
                    vdj.log('GERMLINE ASSIGNMENT ERROR:',
                            'V-gene ({}) and J-gene ({}) chains do not match'.format(vdj.v.chain, vdj.j.chain))
                    self.unassigned.append(vdj)
                    continue
                dquery = self.get_dquery_sequence(jquery, jbr)
                dquery_seqs.append(dquery)
                _vdjs.append(vdj)
            except:
                vdj.exception('V-GENE ASSIGNMENT ERROR', traceback.format_exc())
                self.unassigned.append(vdj)
        os.unlink(jblast_infile)
        vdjs = _vdjs

        # assign D-genes
        _vdjs = []
        for vdj, dquery in zip(vdjs, dquery_seqs):
            if all([vdj.v.chain == 'heavy', dquery]):
                try:
                    germ = self.assign_dgene(dquery, self.species)
                    vdj.d = germ
                except:
                    vdj.exception('D-GENE ASSIGNMENT ERROR:', traceback.format_exc())
                    self.unassigned.append(vdj)
                    continue
            _vdjs.append(vdj)
        self.assigned = _vdjs


    # @property
    # def name(self):
    #     return 'blastn'


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
        blast_db_path = os.path.join(self.germline_directory, 'blast/{}'.format(segment.lower()))
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
        with open(blastout.name, 'rb') as blastout_handle:
            blast_records = [br for br in NCBIXML.parse(blastout_handle)]
        os.unlink(blastout.name)
        return blast_records


    def assign_dgene(self, seq, species):
        db_file = os.path.join(self.germline_directory, 'ungapped/d.fasta')
        with open(db_file, 'r') as db_handle:
            germs = [Sequence(s) for s in SeqIO.parse(db_handle, 'fasta')]
            rc_germs = [Sequence(s.reverse_complement, id=s.id) for s in germs]
            germs.extend(rc_germs)
        alignments = local_alignment(seq, targets=germs,
                                     gap_open=-20, gap_extend=-2)
        alignments.sort(key=lambda x: x.score, reverse=True)
        all_gls = [a.target.id for a in alignments]
        all_scores = [a.score for a in alignments]
        if not all([all_gls, all_scores]):
            return None
        top_gl = all_gls[0]
        top_score = all_scores[0]
        others = [GermlineSegment(germ, species, score=score) for germ, score in zip(all_gls[1:6], all_scores[1:6])]
        return GermlineSegment(top_gl, species, score=top_score, others=others, assigner_name=self.name)


    def process_blast_record(self, blast_record, species):
        all_gls = [a.title.split()[0] for a in blast_record.alignments]
        all_scores = [a.hsps[0].bits for a in blast_record.alignments]
        top_gl = all_gls[0]
        top_score = all_scores[0]
        others = [GermlineSegment(germ, species, score=score) for germ, score in zip(all_gls[1:], all_scores[1:])]
        return GermlineSegment(top_gl, species, score=top_score, others=others[:5], assigner_name=self.name)


    @staticmethod
    def orient_query(vdj, vbr):
        hsp = vbr.alignments[0].hsps[0]
        # BLASTn always reverse complements the Subject sequence, never the query.
        # To determine whether the input sequence is the reverse complement, check
        # to see if the Subject sequence was reverse-complemented by BLASTn
        if hsp.sbjct_start > hsp.sbjct_end:
            vdj.oriented = Sequence(vdj.sequence.reverse_complement, id=vdj.sequence.id)


    @staticmethod
    def get_jquery_sequence(seq, vbr):
        hsp = vbr.alignments[0].hsps[0]
        # check to see if the raw input was reverse-complemented
        if hsp.sbjct_start > hsp.sbjct_end:
            # since the BLASTn alignment was done on the raw input
            # (which has since been reverse complemented), we need
            # to take the portion of the sequence that was 5' of the alignment
            # with the raw input (which is the 3' end of the correctly oriented sequence)
            return Sequence(seq[-hsp.query_start:], id=seq.id)
        else:
            return Sequence(seq[hsp.query_end:], id=seq.id)


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
