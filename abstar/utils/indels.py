#!/usr/bin/env python
# filename: indels.py

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


from __future__ import absolute_import, division, print_function, unicode_literals

import re
import traceback

from abutils.utils import log
from abutils.utils.decorators import lazy_property



class Indel(object):
    """Base class for insertions and deletions"""
    def __init__(self, indel):
        super(Indel, self).__init__()
        self.raw = indel
        self.length = indel.get('len', None)
        self.raw_position = indel.get('pos', None)
        self.sequence = indel.get('seq', None)
        self.fixed = indel.get('fixed', None)
        self.in_frame = self._get_frame()
        # IMGT numbering gets annotated later by the parent Germline object
        # once IMGT numbering for the germline segment is computed.
        self.imgt_position = None
        self.imgt_codon = None


    def __contains__(self, item):
        return item in self.raw.keys()


    def __getitem__(self, key):
        return self.raw.get(key, None)


    def _get_frame(self):
        if 'in frame' in self.raw:
            if self.raw['in frame'] in ['yes', 'no']:
                return self.raw['in frame']
            elif self.raw['in frame'] is True:
                return 'yes'
            elif self.raw['in frame'] is False:
                return 'no'
        return None



class Insertion(Indel):
    """docstring for Insertion"""
    def __init__(self, indel):
        super(Insertion, self).__init__(indel)
        self.type = 'insertion'


    @property
    def imgt_formatted(self):
        return '{}^{}>ins^{}{}'.format(self.imgt_position,
                                       self.imgt_position + 1,
                                       self.sequence.lower(),
                                       '' if self.in_frame else '#')


    @property
    def abstar_formatted(self):
        return '{}:{}{}>{}'.format(self.imgt_position,
                                   self.length,
                                   '' if self.in_frame else '!',
                                   self.sequence,)


    @property
    def json_formatted(self):
        j = {'in_frame': 'yes' if self.in_frame else 'no',
             'length': self.length,
             'sequence': self.sequence,
             'position': self.imgt_position,
             'codon': self.imgt_codon}
        return j



class Deletion(Indel):
    """docstring for Deletion"""
    def __init__(self, indel):
        super(Deletion, self).__init__(indel)
        self.type = 'deletion'


    @property
    def imgt_formatted(self):
        if self.length == 1:
            return '{}{}>del#'.format(self.sequence.lower(),
                                      self.imgt_position)
        else:
            return '{}{}-{}{}>del({}nt){}'.format(self.sequence[0].lower(),
                                                  self.imgt_position,
                                                  self.sequence[-1].lower(),
                                                  self.imgt_position + self.length - 1,
                                                  self.length,
                                                  '' if self.in_frame else '#')


    @property
    def abstar_formatted(self):
        if self.length == 1:
            return '{}:{}{}>{}'.format(self.imgt_position,
                                       self.length,
                                       '' if self.in_frame else '!',
                                       self.sequence,)
        else:
            return '{}-{}:{}{}>{}'.format(self.imgt_position,
                                          self.imgt_position + self.length - 1,
                                          self.length,
                                          '' if self.in_frame else '!',
                                          self.sequence,)

    @property
    def json_formatted(self):
        j = {'in_frame': 'yes' if self.in_frame else 'no',
             'length': self.length,
             'sequence': self.sequence,
             'position': self.imgt_position,
             'codon': self.imgt_codon}
        return j


# ------------
#  Insertions
# ------------

def find_insertions(antibody, segment):
    '''
    Identifies and annotates/fixes insertions. Frameshift insertions (those with a length
    not evenly divisible by 3) will be removed. Codon-length insertions will be annotated.

    Input is a BlastResult object.

    Output is a list of insertion annotations, or None if there were no codon-length
    insertions.
    '''
    try:
        insertions = []
        o = 0
        for i in re.finditer('-+', segment.germline_alignment):
            s = i.start() - o
            e = i.end() - o
            l = e - s
            ins_sequence = segment.query_alignment[s:e]
            if l % 3 == 0 or l > 3:
                insertions.append(_annotate_insertion(segment, s + segment.query_start, l, ins_sequence))
            else:
                insertions.append(_annotate_insertion(segment, s + segment.query_start, l, ins_sequence, fixed=True))
                _fix_frameshift_insertion(antibody, segment, s, e)
                o += l
        return insertions
    except:
        segment.exception('FIND INSERTIONS ERROR', traceback.format_exc(), sep='\n')


def _annotate_insertion(segment, start, length, sequence, fixed=False):
    '''
    Annotates codon-length (non-frameshift) insertions.

    Input is a BlastResult object, the starting postion of the insertion (s) and the
    ending position of the insertion (e).

    Output is a dict containing insertion start position, the insertion length, and
    the inserted sequence.
    '''
    in_frame = 'yes' if length % 3 == 0 else 'no'
    return Insertion({'pos': start, 'len': length, 'seq': sequence, 'in frame': in_frame, 'fixed': fixed})


def _fix_frameshift_insertion(antibody, segment, s, e):
    '''
    Fixes (removes) frameshift insertions.

    Input is a Germline object (segment), the starting postion of the insertion (s) and the
    ending position of the insertion (e).
    '''
    # remove the frameshift indel from the alignment
    segment.query_alignment = segment.query_alignment[:s] + segment.query_alignment[e:]
    segment.germline_alignment = segment.germline_alignment[:s] + segment.germline_alignment[e:]
    segment.alignment_midline = segment.alignment_midline[:s] + segment.alignment_midline[e:]

    # remove the frameshift indel from the oriented input
    oi_firsthalf = antibody.oriented_input.sequence[:s + segment.query_start]
    oi_secondhalf = antibody.oriented_input.sequence[e + segment.query_start:]
    antibody.oriented_input.sequence = oi_firsthalf + oi_secondhalf

    # adjust the alignment end position
    segment.query_end -= (e - s)


# ------------
#  Deletions
# ------------

def find_deletions(antibody, segment):
    '''
    Identifies and annotates/fixes deletions. Frameshift deletions (those with a length
    not evenly divisible by 3) will be removed. Codon-length deletions will be annotated.

    Input is a BlastResult object.

    Output is a list of deletion annotations, or None if there were no codon-length
    deletions.
    '''
    try:
        deletions = []
        o = 0
        for i in re.finditer('-+', segment.query_alignment):
            s = i.start() - o
            e = i.end() - o
            l = e - s
            del_sequence = segment.germline_alignment[s:e]
            if l % 3 == 0 or l > 3:
                deletions.append(_annotate_deletion(segment, s + segment.query_start, l, del_sequence))
            else:
                deletions.append(_annotate_deletion(segment, s + segment.query_start, l, del_sequence, fixed=True))
                _fix_frameshift_deletion(antibody, segment, s, e)
                o += l
        return deletions if deletions else []
    except:
        segment.exception('FIND DELETIONS ERROR', traceback.format_exc(), sep='\n')


def _annotate_deletion(segment, start, length, sequence, fixed=False):
    '''
    Annotates codon-length (non-frameshift) deletions.

    Input is a BlastResult object, the starting postion of the deletion (s) and the
    ending position of the deletion (e).

    Output is a dict containing deletion start position, the deletion length, and
    the deleted sequence.
    '''
    in_frame = 'yes' if length % 3 == 0 else 'no'
    return Deletion({'pos': start, 'len': length, 'seq': sequence, 'in frame': in_frame, 'fixed': fixed})


def _fix_frameshift_deletion(antibody, segment, s, e):
    '''
    Fixes (removes) frameshift deletions.

    Input is a BlastResult object, the starting postion of the deletion (s) and the
    ending position of the deletion (e).

    Output is a modified BlastResult object.
    '''
    # remove the frameshift indel from the alignment
    segment.query_alignment = segment.query_alignment[:s] + segment.germline_alignment[s:e] + segment.query_alignment[e:]
    segment.alignment_midline = segment.alignment_midline[:s] + ''.join(['|'] * (e - s)) + segment.alignment_midline[e:]

    # fix the frameshift deletion in the oriented input
    oi_firsthalf = antibody.oriented_input.sequence[:s + segment.query_start]
    middle = segment.germline_alignment[s:e]
    oi_secondhalf = antibody.oriented_input.sequence[s + segment.query_start:]
    antibody.oriented_input.sequence = oi_firsthalf + middle + oi_secondhalf

    # adjust the alignment end position
    segment.query_end += (e - s)
