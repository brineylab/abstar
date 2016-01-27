#!/usr/bin/env python
# filename: indels.py

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


import re
import traceback

from abtools import log


# ------------
#  Insertions
# ------------

def find_insertions(blast_result):
    '''
    Identifies and annotates/fixes insertions. Frameshift insertions (those with a length
    not evenly divisible by 3) will be removed. Codon-length insertions will be annotated.

    Input is a BlastResult object.

    Output is a list of insertion annotations, or None if there were no codon-length
    insertions.
    '''
    logger = log.get_logger(__name__)
    try:
        insertions = []
        o = 0
        for i in re.finditer('-+', blast_result.germline_alignment):
            s = i.start() - o
            e = i.end() - o
            l = e - s
            if l % 3 == 0:
                insertions.append(_nonframeshift_insertion(blast_result, s, e))
            else:
                blast_result = _fix_frameshift_insertion(blast_result, s, e)
                o += l
        return insertions if insertions else []
    except:
        logger.debug('FIND INSERTIONS ERROR: {}, {}'.format(blast_result.id,
                                                             blast_result.input_sequence))
        logger.debug(traceback.format_exc())


def _nonframeshift_insertion(br, s, e):
    '''
    Annotates codon-length (non-frameshift) insertions.

    Input is a BlastResult object, the starting postion of the insertion (s) and the
    ending position of the insertion (e).

    Output is a dict containing insertion start position, the insertion length, and
    the inserted sequence.
    '''
    ins_start = s + br.germline_start
    ins_length = e - s
    ins_seq = br.query_alignment[s:e]
    return {'pos': ins_start, 'len': ins_length, 'seq': ins_seq}


def _fix_frameshift_insertion(br, s, e):
    '''
    Fixes (removes) frameshift insertions.

    Input is a BlastResult object, the starting postion of the insertion (s) and the
    ending position of the insertion (e).

    Output is a modified BlastResult object.
    '''
    br.query_alignment = br.query_alignment[:s] + br.query_alignment[e:]
    br.germline_alignment = br.germline_alignment[:s] + br.germline_alignment[e:]
    br.alignment_midline = br.alignment_midline[:s] + br.alignment_midline[e:]
    br.fs_indel_adjustment += 1
    return br


# ------------
#  Deletions
# ------------

def find_deletions(blast_result):
    '''
    Identifies and annotates/fixes deletions. Frameshift deletions (those with a length
    not evenly divisible by 3) will be removed. Codon-length deletions will be annotated.

    Input is a BlastResult object.

    Output is a list of deletion annotations, or None if there were no codon-length
    deletions.
    '''
    logger = log.get_logger(__name__)
    try:
        deletions = []
        o = 0
        for i in re.finditer('-+', blast_result.query_alignment):
            s = i.start() - o
            e = i.end() - o
            l = e - s
            if l % 3 == 0:
                deletions.append(_nonframeshift_deletion(blast_result, s, e))
            else:
                blast_result = _fix_frameshift_deletion(blast_result, s, e)
                o += l
        return deletions if deletions else []
    except:
        logger.debug('FIND DELETIONS ERROR: {}, {}'.format(blast_result.id,
                                                             blast_result.input_sequence))
        logger.debug(traceback.format_exc())


def _nonframeshift_deletion(br, s, e):
    '''
    Annotates codon-length (non-frameshift) deletions.

    Input is a BlastResult object, the starting postion of the deletion (s) and the
    ending position of the deletion (e).

    Output is a dict containing deletion start position, the deletion length, and
    the deleted sequence.
    '''
    del_start = s + br.germline_start
    del_length = e - s
    del_seq = br.germline_alignment[s:e]
    br.nfs_indel_adjustment -= del_length
    return {'pos': del_start, 'len': del_length, 'seq': del_seq}


def _fix_frameshift_deletion(br, s, e):
    '''
    Fixes (removes) frameshift deletions.

    Input is a BlastResult object, the starting postion of the deletion (s) and the
    ending position of the deletion (e).

    Output is a modified BlastResult object.
    '''
    br.query_alignment = br.query_alignment[:s] + br.germline_alignment[s:e] + br.query_alignment[e:]
    br.fs_indel_adjustment -= 1
    return br
