#!/usr/bin/python
# filename: regions.py

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


import os
import re
import sys
import traceback

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from abtools import log


def regions(blast_result):
    '''
    Returns a VarRegions (or JoinRegions) object, containing a variety of information about
    Variable (or Joining) gene regions.  Both Regions objects contain equivalent attributes
    to ease downstream handling.

    Input is a BlastResult object for a variable gene.
    '''
    global logger
    logger = log.get_logger(__name__)
    try:
        if blast_result.gene_type == 'variable':
            return VarRegions(blast_result)
        if blast_result.gene_type == 'joining':
            return JoinRegions(blast_result)
    except Exception:
        logger.debug('REGIONS ERROR: {} {}'.format(blast_result.id, blast_result.seq))
        logger.debug(traceback.format_exc())



class VarRegions(object):
    '''
    Structure for information about antibody Variable gene regions.

    Input is a BlastResult object for a variable gene.
    '''
    def __init__(self, blast_result):
        self.raw_positions = self._get_region_positions(blast_result)
        self.adjusted_positions = self._adjusted_region_positions(blast_result)
        self.nt_seqs = self._get_region_nt_sequences(blast_result)
        self.nt_lengths = self._get_region_lengths(self.nt_seqs)
        self.aa_seqs = self._translate_regions(blast_result)
        self.aa_lengths = self._get_region_lengths(self.aa_seqs)
        self.germline_nt_seqs = self._get_region_nt_sequences(blast_result, germline=True)
        self.germline_aa_seqs = self._translate_regions(blast_result, germline=True)


    def _get_region_positions(self, br):
        '''
        Parses the region positions file to get the region positions for a
        specific germline gene.

        Input is a BlastResult object.

        Output is a list of region end positions, with the end position being
        'None' if the region is not present in the sequence (for example, if the
        sequence is truncated and doesn't contain FR1)
        '''
        region_positions = []
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        regions_file = os.path.join(mod_dir, 'utils/germline_data/v_region_info_{}'.format(br.species.lower()))
        with open(regions_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                sline = line.split()
                if sline[0] == br.top_germline:
                    region_positions = [int(i) for i in sline[1:]]
                    break
        if br.insertions or br.deletions:
            region_positions = self._recalibrate_region_positions(br, region_positions)
        return region_positions


    def _recalibrate_region_positions(self, br, positions):
        '''
        Recalibrates the region position numbering to account for indels.

        Input is a BlastResult object.

        Output is a list of recalibrated region end positions.
        '''
        if br.insertions:
            insertions = [(i['pos'], i['len']) for i in br.insertions]
            insertions.sort(key=lambda x: x[0])
            for i in insertions:
                istart = i[0]
                ilength = i[1]
                adjustments = []
                for p in positions:
                    adjustments.append(0 if p <= istart else ilength)
                positions = [p + a for p, a in zip(positions, adjustments)]
        if br.deletions:
            deletions = [(d['pos'], d['len']) for d in br.deletions]
            deletions.sort(key=lambda x: x[0])
            for d in deletions:
                dstart = d[0]
                dlength = d[1]
                adjustments = []
                for p in positions:
                    adjustments.append(0 if p < dstart else dlength)
                positions = [p - a for p, a in zip(positions, adjustments)]
        return positions


    def _adjusted_region_positions(self, br):
        '''
        Adjusts the region positions relative to the start of the query sequence.
        If the query sequence is sufficiently truncated that a region is not present,
        the end position for that region will be set to 'None'.

        Input is a BlastResult object.

        Output is a list of adjusted region end positions.
        '''
        adjusted_positions = [p - br.germline_start for p in self.raw_positions]
        return [p if p > 0 else None for p in adjusted_positions]


    def _get_region_nt_sequences(self, br, germline=False):
        '''
        Identifies the nucleotide sequence for each region of a variable gene.

        Input is a BlastResult object and, optionally, the germline flag. If set
        to True, region sequences for the top germline match are returned.

        Output is a dict of region sequences.
        '''
        regions = self._regions(br)
        region_nt_seqs = {}
        start = 0
        if germline:
            alignment = br.germline_alignment
        else:
            alignment = br.query_alignment
        for r in regions:
            reg = r[0]
            end = r[1]
            if end:
                region_nt_seqs[reg] = alignment[start:end]
                start = end
            else:
                region_nt_seqs[reg] = None
        return region_nt_seqs


    def _region_list(self):
        return ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3']


    def _get_region_lengths(self, regions):
        lengths = {}
        for r in regions:
            if regions[r]:
                lengths[r] = len(regions[r])
            else:
                lengths[r] = 0
        return lengths


    def _regions(self, region_positions):
        '''
        Produces a list of regions and region end positions for all regions.

        Input is a list of region positions.

        Output is a list of tuples, of the format:
        (region_name, region_end_position)
        '''
        all_regions = self._region_list()
        return zip(all_regions, self.adjusted_positions)


    def _translate_regions(self, br, germline=False):
        '''
        Translates the region sequence.

        Input is a BlastResult object and, optionally, the germline flag. If germline is set
        to True, translated regions for the top germline match will be returned.

        Output is a dict of region sequences.
        '''
        if germline:
            nt_seqs = self.germline_nt_seqs
        else:
            nt_seqs = self.nt_seqs
        aa_seqs = {}
        rf = br.germline_start % 3
        rf_offset = (rf * 2) % 3
        first_region = True
        regions = self._region_list()
        for i, r in enumerate(regions):
            if r in nt_seqs.keys():
                if not nt_seqs[r]:
                    aa_seqs[r] = None
                    continue
                if i <= len(regions) - 2:
                    nt_seqs = self._fix_region_spanning_indel(nt_seqs, r, regions[i + 1])
                if first_region:
                    aa_seqs[r] = str(Seq(nt_seqs[r][rf_offset:].replace('-', ''), generic_dna).translate())
                    if '-' in nt_seqs[r]:
                        aa_seqs[r] = self._fix_translated_indels(nt_seqs[r], aa_seqs[r])
                    first_region = False
                else:
                    trim = len(nt_seqs[r].replace('-', '')) % 3
                    seq_to_trans = nt_seqs[r].replace('-', '')[:-trim] if trim else nt_seqs[r].replace('-', '')
                    aa_seqs[r] = str(Seq(seq_to_trans, generic_dna).translate())
                    if '-' in nt_seqs[r]:
                        aa_seqs[r] = self._fix_translated_indels(nt_seqs[r], aa_seqs[r])
            else:
                aa_seqs[r] = None
        if germline:
            self.germline_nt_seqs = nt_seqs
        else:
            self.nt_seqs = nt_seqs
        return aa_seqs


    def _fix_region_spanning_indel(self, nt_seqs, r1, r2):
        r1_len = len(nt_seqs[r1]) - len(nt_seqs[r1].rstrip('-'))
        if r1_len % 3 == 0:
            return nt_seqs
        r2_len = len(nt_seqs[r2]) - len(nt_seqs[r2].lstrip('-'))
        if r1_len > r2_len:
            len_to_move = r2_len % 3
            chunk_to_move = nt_seqs[r1].lstrip('-')[:len_to_move]
            r1_gap = r1_len + len_to_move
            r2_gap = r2_len - len_to_move
            nt_seqs[r1] = '{}{}'.format(nt_seqs[r1].rstrip('-')[:-len_to_move], '-' * r1_gap)
            nt_seqs[r2] = '{}{}{}'.format('-' * r2_gap, chunk_to_move, nt_seqs[r2].lstrip('-'))
        else:
            len_to_move = r1_len % 3
            chunk_to_move = nt_seqs[r1].rstrip('-')[-len_to_move:]
            r1_gap = r1_len - len_to_move
            r2_gap = r2_len + len_to_move
            nt_seqs[r1] = '{}{}{}'.format(nt_seqs[r1].rstrip('-'), chunk_to_move, '-' * r1_gap)
            nt_seqs[r2] = '{}{}'.format('-' * r2_gap, nt_seqs[r1].lstrip('-')[len_to_move:])
        return nt_seqs


    def _fix_translated_indels(self, nt_seq, aa_seq):
        for indel in re.finditer('-+', nt_seq):
            s = indel.start() / 3
            e = indel.end() / 3
            i = '-' * (e - s)
            return aa_seq[:s] + i + aa_seq[s:]




class JoinRegions(object):
    '''
    Structure for information about antibody Joining gene regions.

    Input is a BlastResult object.
    '''
    def __init__(self, br):
        self.fix_v_overlap = False
        self.v_overlap_length = 0
        start = self._get_fr4_nt_start(br)
        end = self._get_fr4_nt_end(br)
        self.raw_positions = [start, end]
        self.adjusted_positions = self.raw_positions
        self.nt_seqs = self._get_fr4_nt_seq(br, start, end)
        self.nt_lengths = self._get_region_lengths(self.nt_seqs)
        self.aa_seqs = self._translate_regions(self.nt_seqs)
        self.aa_lengths = self._get_region_lengths(self.aa_seqs)
        self.germline_nt_seqs = self._get_fr4_germ_nt_seq(br, start, end)
        self.germline_aa_seqs = self._translate_regions(self.germline_nt_seqs)


    def _get_fr4_nt_start(self, br):
        '''
        Returns the start position (in the joining gene alignment) of FR4.

        Input is a BlastResult object.
        '''
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        regions_file = os.path.join(mod_dir, 'utils/germline_data/j_region_info_{}'.format(br.species.lower()))
        with open(regions_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                sline = line.split()
                if sline[0] == br.top_germline:
                    raw_fr4_start = int(sline[2])
                    break
        # in some rare cases where the V-gene alignment overlaps
        # the J-gene alignment, which results in the first nt or two
        # of FR4 being truncated (because they're part of the V-gene alignment)
        if br.chain != 'heavy' and raw_fr4_start - br.germline_start < 0:
            self.fix_v_overlap = True
            self.v_overlap_length = abs(raw_fr4_start - br.germline_start)
            logger.debug('V-OVERLAP FOUND: {}'.format(br.id))
            logger.debug('V-OVERLAP LENGTH: {}'.format(self.v_overlap_length))
        return max(0, raw_fr4_start - br.germline_start)


    def _get_fr4_nt_end(self, br):
        '''
        Returns the end position (in the joining gene alignment) of FR4.

        Input is a BlastResult object.
        '''
        return br.query_end


    def _get_fr4_nt_seq(self, br, start, end):
        '''
        Returns the nucleotide sequence of FR4.

        Input is a BlastResult object and the start and end positions of FR4.
        '''
        return {'FR4': br.query_alignment[start:end]}


    def _get_fr4_germ_nt_seq(self, br, start, end):
        '''
        Returns the germline nucleotide sequence of FR4.

        Input is a BlastResult object and the start and end positions of FR4.
        '''
        return {'FR4': br.germline_alignment[start:end]}


    def _translate_regions(self, regions):
        '''
        Returns the translated sequence for each region.

        Input is a dict of nucleotide sequences for each region.
        '''
        translated = {}
        for r in regions:

            trim = len(regions[r].replace('-', '')) % 3
            seq_for_trans = regions[r].replace('-', '')[:-trim] if trim else regions[r].replace('-', '')
            translated[r] = str(Seq(seq_for_trans, generic_dna).translate())
            if '-' in regions[r]:
                translated[r] = self._fix_translated_indels(regions[r], translated[r])
        return translated


    def _fix_translated_indels(self, nt_seq, aa_seq):
        for indel in re.finditer('-+', nt_seq):
            s = indel.start() / 3
            e = indel.end() / 3
            i = '-' * (e - s)
            return aa_seq[:s] + i + aa_seq[s:]


    def _get_region_lengths(self, regions):
        '''
        Returns the sequence lengths of each region.

        Input is a dict of nucleotide or amino acid sequences for each region.
        '''
        lengths = {}
        for r in regions:
            if regions[r]:
                lengths[r] = len(regions[r])
            else:
                lengths[r] = 0
        return lengths
