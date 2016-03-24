#!/usr/bin/python
# filename: blast.py

#
# Copyright (c) 2015 Bryan Briney
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


from __future__ import print_function

import os
import platform
import tempfile

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

from abtools.alignment import local_alignment, local_alignment_biopython


class BlastResult(object):
    '''
    Data structure for parsing and holding a BLASTn result.
    Input is a file handle for the XML-formatted BLASTn output file.
    '''
    def __init__(self, seq, blastout, species):
        super(BlastResult, self).__init__()
        self.seq = seq
        self.id = seq.id
        self.alignments = blastout.alignments
        self.input_sequence = seq.sequence
        self.species = species
        self.top_germline = self._get_top_germline()
        self.all_germlines = self._get_all_germlines()
        self.top_score = self._get_top_score()
        self.all_scores = self._get_all_scores()
        self.top_evalue = self._get_top_evalue()
        self.all_evalues = self._get_all_evalues()
        self.top_bitscore = self._get_top_bitscore()
        self.all_bitscores = self._get_all_bitscores()
        self.strand = 'plus'
        self.query_alignment = self._get_query_alignment()
        self.germline_alignment = self._get_germline_alignment()
        self.alignment_midline = self._get_alignment_midline()
        self.alignment_length = self._get_alignment_length()
        self.query_start = self._get_query_start()
        self.query_end = self._get_query_end()
        self.germline_start = self._get_germline_start()
        self.germline_end = self._get_germline_end()
        self.gene_type = self._gene_type()
        self.chain = self._chain()

    def annotate(self):
        self.fs_indel_adjustment = 0
        self.nfs_indel_adjustment = 0
        self._fix_ambigs()
        self._find_indels()
        self.regions = self._regions()
        self.nt_mutations = self._nt_mutations()
        self.aa_mutations = self._aa_mutations()


    def realign_variable(self, germline_gene):
        '''
        Due to restrictions on the available scoring parameters in BLASTn, incorrect truncation
        of the v-gene alignment can occur. This function re-aligns the query sequence with
        the identified germline variable gene using more appropriate alignment parameters.

        Input is the name of the germline variable gene (ex: 'IGHV1-2*02').
        '''
        self.germline_seq = self._get_germline_sequence_for_realignment(germline_gene, 'V')
        alignment = local_alignment(self.seq.sequence, self.germline_seq,
                                    match=3, mismatch=-2,
                                    gap_open_penalty=22, gap_extend_penalty=1)
        rc = self.seq.reverse_complement
        alignment_rc = local_alignment(rc, self.germline_seq,
                                       match=3, mismatch=-2,
                                       gap_open_penalty=22, gap_extend_penalty=1)
        if alignment.score > alignment_rc.score:
            self._process_realignment(alignment)
        else:
            self.strand = 'minus'
            self.input_sequence = rc
            self._process_realignment(alignment_rc)


    def realign_variable_biopython(self, germline_gene):
        '''
        Due to restrictions on the available scoring parameters in BLASTn, incorrect truncation
        of the v-gene alignment can occur. This function re-aligns the query sequence with
        the identified germline variable gene using more appropriate alignment parameters.

        Input is the name of the germline variable gene (ex: 'IGHV1-2*02').
        '''
        self.germline_seq = self._get_germline_sequence_for_realignment(germline_gene, 'V')
        alignment = local_alignment_biopython(self.seq.sequence, self.germline_seq,
                                    match=6, mismatch=-2,
                                    gap_open_penalty=-22, gap_extend_penalty=0)
        rc = self.seq.reverse_complement
        alignment_rc = local_alignment_biopython(rc, self.germline_seq,
                                       match=6, mismatch=-2,
                                       gap_open_penalty=-22, gap_extend_penalty=0)
        if alignment.score > alignment_rc.score:
            self._process_realignment(alignment)
        else:
            self.strand = 'minus'
            self.input_sequence = rc
            self._process_realignment(alignment_rc)


    def _get_germline_sequence_for_realignment(self, germ, gene):
        '''
        Identifies the appropriate germline variable gene from a database of all
        germline variable genes.

        Input is the name of the germline variable gene (ex: 'IGHV1-2*02') and
        the gene region ('V' or 'J').

        Output is the germline sequence.
        '''
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        db_file = os.path.join(mod_dir, 'ssw/dbs/{}_{}.fasta'.format(self.species.lower(), gene))
        for s in SeqIO.parse(open(db_file), 'fasta'):
            if s.id == germ:
                return str(s.seq)
        return None


    def _process_realignment(self, alignment):
        '''
        Processes the result of variable gene realignment and updates BlastResult
        attributes accordingly.

        Input is an AbTools SSWAlignment object.
        '''
        self.query_alignment = alignment.aligned_query
        self.germline_alignment = alignment.aligned_target
        self.alignment_midline = alignment.alignment_midline
        if self.gene_type == 'variable':
            self.query_start = alignment.query_begin
            self.germline_start = alignment.target_begin
        self.query_end = alignment.query_end
        self.germline_end = alignment.target_end


    def _fix_ambigs(self):
        '''
        Fixes ambiguous nucleotides by replacing them with the germline nucleotide.
        '''
        from abstar.utils import ambigs
        self.query_alignment = ambigs.fix_ambigs(self)


    def _find_indels(self):
        '''
        Identifies and annotates indels in the query sequence.
        '''
        self.insertions = []
        self.deletions = []
        if self._indel_check():
            from abstar.utils import indels
            self.insertions = indels.find_insertions(self)
            self.deletions = indels.find_deletions(self)

    def _indel_check(self):
        '''
        Checks the sequence for evidence of insertions or deletions.
        '''
        if '-' in self.query_alignment:
            return True
        elif '-' in self.germline_alignment:
            return True
        return False

    def _regions(self):
        '''
        Identifies and annotates variable/joining gene regions.
        '''
        from abstar.utils import regions
        return regions.regions(self)

    def _nt_mutations(self):
        '''
        Identifies and annotates nucleotide mutations.
        '''
        from abstar.utils import mutations
        return mutations.nt_mutations(self)

    def _aa_mutations(self):
        '''
        Identifies and annotates amino acid mutations.
        '''
        from abstar.utils import mutations
        return mutations.aa_mutations(self)

    def _get_top_germline(self):
        'Returns the top scoring germline gene'
        return self.alignments[0].title.split()[0]

    def _get_all_germlines(self):
        'Returns all germline genes'
        return [alignment.title.split()[0] for alignment in self.alignments]

    def _get_top_score(self):
        'Returns the score for the top scoring germline gene'
        top_alignment = self.alignments[0]
        return top_alignment.hsps[0].score

    def _get_all_scores(self):
        'Returns all germline gene scores'
        return [alignment.hsps[0].score for alignment in self.alignments]

    def _get_top_evalue(self):
        'Returns the e-value for the top scoring germline gene'
        top_alignment = self.alignments[0]
        return top_alignment.hsps[0].expect

    def _get_all_evalues(self):
        return [alignment.hsps[0].expect for alignment in self.alignments]

    def _get_top_bitscore(self):
        'Returns the bitscore for the top scoring germline gene'
        top_alignment = self.alignments[0]
        return top_alignment.hsps[0].bits

    def _get_all_bitscores(self):
        return [alignment.hsps[0].bits for alignment in self.alignments]

    def _get_strand(self):
        top_alignment = self.alignments[0]
        return top_alignment.hsps[0].strand

    def _get_query_alignment(self):
        '''Returns the query alignment string for the
        top scoring germline gene alignment'''
        top_alignment = self.alignments[0]
        return top_alignment.hsps[0].query

    def _get_germline_alignment(self):
        'Returns the top scoring germline alignment'
        top_alignment = self.alignments[0]
        return top_alignment.hsps[0].sbjct

    def _get_alignment_midline(self):
        '''Returns the alignment midline string for the
        top scoring germline gene alignment'''
        top_alignment = self.alignments[0]
        return top_alignment.hsps[0].match

    def _get_alignment_length(self):
        'Returns the alignment length for the top scoring germline gene alignment'
        return self.alignments[0].length

    def _get_query_start(self):
        '''Returns the start position of the query alignment with
        the top germline gene'''
        top_alignment = self.alignments[0]
        return top_alignment.hsps[0].query_start - 1

    def _get_query_end(self):
        '''Returns the start position of the query alignment with
        the top germline gene'''
        return self.query_start + self.alignment_length

    def _get_germline_start(self):
        'Returns the start position of the top germline alignment'
        top_alignment = self.alignments[0]
        return top_alignment.hsps[0].sbjct_start - 1

    def _get_germline_end(self):
        'Returns the end position of the top germline alignment'
        return self.germline_start + self.alignment_length

    def _gene_type(self):
        "Returns the gene type, either 'variable' or 'joining'"
        if self.top_germline[3] == 'V':
            return 'variable'
        elif self.top_germline[3] == 'J':
            return 'joining'
        return None

    def _chain(self):
        'Returns the chain type'
        if self.top_germline.startswith('IGH'):
            return 'heavy'
        elif self.top_germline.startswith('IGK'):
            return 'kappa'
        elif self.top_germline.startswith('IGL'):
            return 'lambda'
        # hack to accomodate IgBLAST's crappy germline database
        # to be removed when I curate the database
        elif self.top_germline.startswith('VH'):
            return 'heavy'
        return None


class DiversityResult(object):
    """
    Data structure for holding information about diversity germline gene assignments.
    Designed to have (mostly) the same attributes as BlastResult objects, so that
    DiversityResult objects and BlastResult objects can (mostly) be used
    interchangeably in downstream operations.

    Main differences between DiversityResult and BlastResult are that DiversityResults
    have empty bitscore and e-value attributes, and the alignments attribute contains
    an Alignment object instead of parsed BLASTn output.

    Input is a list of Alignment objects, representing the top-scoring diversity genes.
    """
    def __init__(self, seq, alignments):
        super(DiversityResult, self).__init__()
        self.id = seq.id
        self.input_sequence = seq.sequence
        self.alignments = alignments
        self.top_germline = self._get_top_germline()
        self.all_germlines = self._get_all_germlines()
        self.top_score = self._get_top_score()
        self.all_scores = self._get_all_scores()
        self.top_evalue = None
        self.all_evalues = []
        self.top_bitscore = None
        self.all_bitscores = []
        self.query_alignment = self._get_query_alignment()
        self.germline_alignment = self._get_germline_alignment()
        self.alignment_midline = self._get_alignment_midline()
        self.alignment_length = self._get_alignment_length()
        self.query_start = self._get_query_start()
        self.query_end = self._get_query_end()
        self.germline_start = self._get_germline_start()
        self.germline_end = self._get_germline_end()
        self.sequence = self._get_sequence()
        self.reading_frame = self._get_reading_frame()
        self.gene_type = 'diversity'
        self.chain = 'heavy'
        self.nt_mutations = self._nt_mutations()

    def _get_top_germline(self):
        '''
        Returns the top scoring germline gene. If no germline gene scores
        higher than 9, then the top_germline attribute is set to 'None'.
        The minimum for reaching a score of 9 would be a 6nt alignment
        region with at least 5 matching nucleotides and a single mismatch.
        '''
        top_alignment = self.alignments[0]
        if top_alignment.score > 9:
            return top_alignment.target.id
        return None

    def _get_all_germlines(self):
        'Returns all germline genes'
        return [a.target.id for a in self.alignments]

    def _get_top_score(self):
        'Returns the score for the top scoring germline gene'
        top_alignment = self.alignments[0]
        return top_alignment.score

    def _get_all_scores(self):
        'Returns all germline gene scores'
        return [a.score for a in self.alignments]

    def _get_query_alignment(self):
        '''Returns the query alignment string for the
        top scoring germline gene alignment'''
        top_alignment = self.alignments[0]
        return top_alignment.aligned_query

    def _get_germline_alignment(self):
        'Returns the top scoring germline alignment'
        top_alignment = self.alignments[0]
        return top_alignment.aligned_target

    def _get_alignment_midline(self):
        '''Returns the alignment midline string for the
        top scoring germline gene alignment'''
        query = self.query_alignment
        germ = self.germline_alignment
        midline = ['|' if q == g else ' ' for q, g in zip(query, germ)]
        return ''.join(midline)

    def _get_alignment_length(self):
        'Returns the alignment length for the top scoring germline gene alignment'
        return len(self.germline_alignment)

    def _get_query_start(self):
        '''Returns the start position of the query alignment with
        the top germline gene'''
        top_alignment = self.alignments[0]
        return top_alignment.query_begin

    def _get_query_end(self):
        '''Returns the start position of the query alignment with
        the top germline gene'''
        top_alignment = self.alignments[0]
        return top_alignment.query_end

    def _get_germline_start(self):
        'Returns the start position of the top germline alignment'
        top_alignment = self.alignments[0]
        return top_alignment.target_begin

    def _get_germline_end(self):
        'Returns the end position of the top germline alignment'
        top_alignment = self.alignments[0]
        return top_alignment.target_end

    def _get_reading_frame(self):
        '''
        Identifies the diverstiy gene reading frame
        '''
        rf = (self.germline_start % 3) + 1
        return rf

    def _get_sequence(self):
        return self.query_alignment[self.query_start:self.query_end]

    def _nt_mutations(self):
        '''
        Identifies and annotates nucleotide mutations.
        '''
        from abstar.utils import mutations
        return mutations.nt_mutations(self)


def build_j_blast_input(seqs, v_blast_results):
    j_fastas = []
    for seq, vbr in zip(seqs, v_blast_results):
        start = len(vbr.query_alignment) + vbr.query_start + vbr.fs_indel_adjustment + vbr.nfs_indel_adjustment
        j_fastas.append('>{}\n{}'.format(seq.id, seq.sequence[start:]))
    j_blastin = tempfile.NamedTemporaryFile(delete=False)
    j_blastin.write('\n'.join(j_fastas))
    j_blastin.close()
    return j_blastin


def blast(seq_file, species, segment):
    '''
    Runs BLASTn against an antibody germline database.

    Input is a FASTA file of sequences (the file path, not a handle), the species of origin
    of the sequences to be queried, and the gene segment (options are: 'V', 'D', or 'J')
    '''
    mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    blast_path = os.path.join(mod_dir, 'blast/blastn_{}'.format(platform.system().lower()))
    blast_db_path = os.path.join(mod_dir, 'blast/dbs/{}_gl_{}'.format(species.lower(), segment.upper()))
    blastout = tempfile.NamedTemporaryFile(delete=False)
    blastn_cmd = NcbiblastnCommandline(cmd=blast_path,
                                       db=blast_db_path,
                                       query=seq_file,
                                       out=blastout.name,
                                       outfmt=5,
                                       dust='no',
                                       word_size=_word_size(segment),
                                       max_target_seqs=10,
                                       evalue=_evalue(segment),
                                       reward=_match_reward(segment),
                                       penalty=_mismatch_penalty(segment),
                                       gapopen=_gap_open(segment),
                                       gapextend=_gap_extend(segment))
    stdout, stderr = blastn_cmd()
    blast_records = [br for br in NCBIXML.parse(blastout)]
    os.unlink(blastout.name)
    return blast_records


def _word_size(segment):
    'Returns BLASTn word size for the given gene segment'
    word_sizes = {'V': 11,
                  'D': 4,
                  'J': 7}
    return word_sizes[segment]


def _gap_open(segment):
    'Returns BLASTn gap-open penalty for the given gene segment'
    gap_open = {'V': 5,
                'D': 4,
                'J': 5}
    return gap_open[segment]


def _gap_extend(segment):
    'Returns BLASTn gap-extend penalty for the given gene segment'
    gap_extend = {'V': 2,
                  'D': 2,
                  'J': 2}
    return gap_extend[segment]


def _match_reward(segment):
    'Returns BLASTn match reward for the given gene segment'
    match = {'V': 1,
             'D': 1,
             'J': 1}
    return match[segment]


def _mismatch_penalty(segment):
    'Returns BLASTn mismatch penalty for the given gene segment'
    mismatch = {'V': -1,
                'D': -1,
                'J': -1}
    return mismatch[segment]


def _evalue(segment):
    'Returns minimum BLASTn e-value for the given gene segment'
    evalue = {'V': 1,
              'D': 100000,
              'J': 1000}
    return evalue[segment]
