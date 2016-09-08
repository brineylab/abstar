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



def regions(antibody):
    global logger
    logger = log.get_logger(__name__)
    # V-gene regions
    self.v.regions = VarRegions(antibody)
    # J-gene regions
    self.j.regions = JoinRegions(antibody)


class VarRegions(object):
    '''
    Structure for information about antibody Variable gene regions.

    Input is a Antibody object containing an assigned variable gene.
    '''
    def __init__(self, antibody):
        segment = antibody.v
        self.raw_positions = self._get_region_positions(segment)
        self.adjusted_positions = self._adjusted_region_positions(segment)
        self.nt_seqs = self._get_region_nt_sequences(segment)
        self.nt_lengths = self._get_region_lengths(self.nt_seqs)
        self.aa_seqs = self._translate_regions(segment)
        self.aa_lengths = self._get_region_lengths(self.aa_seqs)
        self.germline_nt_seqs = self._get_region_nt_sequences(segment, germline=True)
        self.germline_aa_seqs = self._translate_regions(segment, germline=True)


    def _get_region_positions(self, segment):
        '''
        Parses the region positions file to get the region positions for a
        specific germline gene.

        Input is a Germline object corresponding to an assigned V-gene.

        Output is a list of region end positions, with the end position being
        'None' if the region is not present in the sequence (for example, if the
        sequence is truncated and doesn't contain FR1)
        '''
        region_positions = []
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        regions_file = os.path.join(mod_dir, 'utils/germline_data/v_region_info_{}'.format(segment.species.lower()))
        with open(regions_file) as f:
            for line in f:
                if line.lstrip().startswith('#'):
                    continue
                sline = line.split()
                if sline[0] == segment.full:
                    region_positions = [int(i) for i in sline[1:]]
                    break
        if segment.insertions:
            region_positions = self._recalibrate_region_positions(segment, region_positions)
        return region_positions


    def _recalibrate_region_positions(self, segment, positions):
        '''
        Recalibrates the region position numbering to account for insertions.

        Input is a BlastResult object and a list of region positions.

        Output is a list of recalibrated region end positions.
        '''
        if segment.insertions:
            insertions = [(i['pos'], i['len']) for i in segment.insertions]
            insertions.sort(key=lambda x: x[0])
            for i in insertions:
                istart = i[0]
                ilength = i[1]
                adjustments = []
                for p in positions:
                    adjustments.append(0 if p <= istart else ilength)
                positions = [p + a for p, a in zip(positions, adjustments)]
        return positions


    def _adjusted_region_positions(self, segment):
        '''
        Adjusts the region positions relative to the start of the query sequence.
        If the query sequence is sufficiently truncated that a region is not present,
        the end position for that region will be set to 'None'.

        Input is a BlastResult object.

        Output is a list of adjusted region end positions.
        '''
        adjusted_positions = [p - segment.germline_start for p in self.raw_positions]
        return [p if p > 0 else None for p in adjusted_positions]


    def _get_region_nt_sequences(self, segment, germline=False):
        '''
        Identifies the nucleotide sequence for each region of a variable gene.

        Input is a BlastResult object and, optionally, the germline flag. If set
        to True, region sequences for the top germline match are returned.

        Output is a dict of region sequences.
        '''
        # regions = self._regions(segment)
        regions = self._regions()
        region_nt_seqs = {}
        start = 0
        if germline:
            alignment = segment.germline_alignment
        else:
            alignment = segment.query_alignment
        for r in regions:
            reg = r[0]
            end = r[1]
            if end is not None:
                region_nt_seqs[reg] = alignment[start:end]
                if end > len(alignment) and reg == 'FR3':
                    logger.debug('TRUNCATED FR3: {}'.format(segment.id))
                    logger.debug('OLD: {}'.format(region_nt_seqs['FR3']))
                    region_nt_seqs[reg] = self._fix_truncated_fr3_alignment(segment, start, end, germline)
                    logger.debug('NEW: {}'.format(region_nt_seqs['FR3']))
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


    # def _regions(self, region_positions):
    def _regions(self):
        '''
        Produces a list of regions and region end positions for all regions.

        Input is a list of region positions.

        Output is a list of tuples, of the format:
        (region_name, region_end_position)
        '''
        all_regions = self._region_list()
        return zip(all_regions, self.adjusted_positions)


    def _translate_regions(self, segment, germline=False):
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
        rf = segment.germline_start % 3
        rf_offset = (rf * 2) % 3
        first_region = True
        regions = self._region_list()
        for i, r in enumerate(regions):
            if r in nt_seqs.keys():
                if not nt_seqs[r]:
                    aa_seqs[r] = None
                    continue
                if i <= len(regions) - 2:
                    nt_seqs = self._fix_region_spanning_indel(nt_seqs, r, regions[i + 1:])
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


    def _fix_region_spanning_indel(self, nt_seqs, r1, next_regions):
        r1_len = len(nt_seqs[r1].rstrip('-'))
        if r1_len % 3 == 0:
            return nt_seqs
        logger.debug(next_regions)
        reg_iter = iter(next_regions)
        r2 = reg_iter.next()
        while len(nt_seqs[r2].strip('-')) == 0:
            logger.debug(r2)
            r2 = reg_iter.next()
        r2_len = len(nt_seqs[r2].lstrip('-'))
        if r1_len % 3 > r2_len % 3:
            len_to_move = r2_len % 3
            chunk_to_move = nt_seqs[r2].lstrip('-')[:len_to_move]
            r1_gap = len(nt_seqs[r1]) - r1_len - len_to_move
            r2_gap = len(nt_seqs[r2]) - r2_len + len_to_move
            nt_seqs[r1] = '{}{}{}'.format(nt_seqs[r1].rstrip('-'), chunk_to_move, '-' * r1_gap)
            nt_seqs[r2] = '{}{}'.format('-' * r2_gap, nt_seqs[r2].lstrip('-')[len_to_move:])
        else:
            len_to_move = r1_len % 3
            chunk_to_move = nt_seqs[r1].rstrip('-')[-len_to_move:]
            r1_gap = len(nt_seqs[r1]) - r1_len + len_to_move
            r2_gap = len(nt_seqs[r2]) - r2_len - len_to_move
            nt_seqs[r1] = '{}{}'.format(nt_seqs[r1].rstrip('-')[:-len_to_move], '-' * r1_gap)
            nt_seqs[r2] = '{}{}{}'.format('-' * r2_gap, chunk_to_move, nt_seqs[r2].lstrip('-'))
        return nt_seqs


    def _fix_translated_indels(self, nt_seq, aa_seq):
        for indel in re.finditer('-+', nt_seq):
            s = indel.start() / 3
            e = indel.end() / 3
            i = '-' * (e - s)
            return aa_seq[:s] + i + aa_seq[s:]


    def _fix_truncated_fr3_alignment(self, segment, start, end, germline):
        if germline:
            alignment = segment.germline_alignment
            seq = segment.germline_seq[segment.germline_end + 1:]
        else:
            alignment = segment.query_alignment
            seq = segment.input_sequence[segment.query_end + 1:]
        trunc = end - len(alignment)
        return alignment[start:] + seq[:trunc]


class JoinRegions(object):
    '''
    Structure for information about antibody Joining gene regions.

    Input is a BlastResult object.
    '''
    def __init__(self, antibody):
        segment = antibody.j
        self.v_overlap_length = 0
        start = self._get_fr4_nt_start(segment)
        end = self._get_fr4_nt_end(segment)
        self.raw_positions = [start, end]
        self.adjusted_positions = self.raw_positions
        self.nt_seqs = self._get_fr4_nt_seq(antibody, start, end)
        self.nt_lengths = self._get_region_lengths(self.nt_seqs)
        self.aa_seqs = self._translate_regions(self.nt_seqs)
        self.aa_lengths = self._get_region_lengths(self.aa_seqs)
        self.germline_nt_seqs = self._get_fr4_germ_nt_seq(antibody, start, end)
        self.germline_aa_seqs = self._translate_regions(self.germline_nt_seqs)


    def _get_fr4_nt_start(self, segment):
        '''
        Returns the start position (in the joining gene alignment) of FR4.

        Input is a Germline object, corresponding to a J-gene assignment.
        '''
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        regions_file = os.path.join(mod_dir, 'utils/germline_data/j_region_info_{}'.format(segment.species))
        with open(regions_file) as f:
            for line in f:
                if line.lstrip().startswith('#'):
                    continue
                sline = line.split()
                if sline[0] == segment.full:
                    raw_fr4_start = int(sline[2])
                    break
        # in some rare cases where the V-gene alignment overlaps
        # the J-gene alignment, which results in the first nt or two
        # of FR4 being truncated (because they're part of the V-gene alignment)
        if segment.chain != 'heavy' and raw_fr4_start - segment.germline_start < 0:
            self.v_overlap_length = abs(raw_fr4_start - segment.germline_start)
            logger.debug('CHAIN: {}'.format(segment.chain))
            logger.debug('RAW GERMLINE FR4 START: {}'.format(raw_fr4_start))
            logger.debug('ALIGNED GERMLINE FR4 START: {}'.format(segment.germline_start))
            logger.debug('V-OVERLAP FOUND: {}'.format(segment.id))
            logger.debug('V-OVERLAP LENGTH: {}'.format(self.v_overlap_length))
        return max(0, raw_fr4_start - segment.germline_start)


    def _get_fr4_nt_end(self, segment):
        '''
        Returns the end position (in the joining gene alignment) of FR4.

        Input is a BlastResult object.
        '''
        return segment.query_end


    def _get_fr4_nt_seq(self, antibody, start, end):
        '''
        Returns the nucleotide sequence of FR4.

        Input is a BlastResult object and the start and end positions of FR4.
        '''
        if antibody.j.v_overlap_length:
            fr4 = antibody.v.query_alignment[-antibody.j.v_overlap_length:] + antibody.j.query_alignment[start:end]
        else:
            fr4 = antibody.j.query_alignment[start:end]
        return {'FR4': fr4}


    def _get_fr4_germ_nt_seq(self, antibody, start, end):
        '''
        Returns the germline nucleotide sequence of FR4.

        Input is a BlastResult object and the start and end positions of FR4.
        '''
        if antibody.j.v_overlap_length:
            fr4 = antibody.v.germline_alignment[-antibody.j.v_overlap_length:] + antibody.j.germline_alignment[start:end]
        else:
            fr4 = antibody.j.germline_alignment[start:end]
        return {'FR4': fr4}


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
