#!/usr/bin/env python
# filename: vdj.py

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

from ..utils.mixins import LoggingMixin

from abutils.core.sequence import Sequence


class VDJ(LoggingMixin):
    """

    Args:
    -----

        sequence: The query sequence, in any format that ```abutils.core.sequence.Sequence``` can accept:

                1. a nucleotide sequence, as a string
                2. a Biopython ``SeqRecord`` object
                3. an abutils ``Sequence`` object
                4. a list/tuple of the format ``[seq_id, sequence]``

            Note that if a plain nucleotide sequence is provided, a random sequence seq_id
            will be generated.

        v (Germline): an AbStar Germline object representing the assigned Variable gene.

        d (Germline): an AbStar Germline object representing the assigned Diversity gene

        j (Germline): an AbStar Germline object representing the assigned Joining gene

    """
    def __init__(self, sequence, v=None, d=None, j=None):
        super(VDJ, self).__init__()
        LoggingMixin.__init__(self)
        self.sequence = Sequence(sequence)
        self.id = self.sequence.id
        self.oriented = self.sequence
        self.v = v
        self.d = d
        self.j = j
        self.initialize_log()


    def initialize_log(self):
        log = []
        log.append('=' * len(self.id))
        log.append(self.id)
        log.append('=' * len(self.id))
        log.append('')
        log.append(self.sequence.fasta)
        log.append('')
        if self.v is not None:
            log.append('V-GENE: {}'.format(self.v.full))
        if self.d is not None:
            log.append('D-GENE: {}'.format(self.d.full))
        if self.j is not None:
            log.append('J-GENE: {}'.format(self.j.full))
        self._log = log
