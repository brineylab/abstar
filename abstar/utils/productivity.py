#!/usr/bin/python
# filename: productivity.py

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


import traceback


def check_productivity(antibody):
    return Productivity(antibody)


class Productivity(object):
    """docstring for Productivity"""
    def __init__(self, antibody):
        super(Productivity, self).__init__()
        self.antibody = antibody
        self.productivity_issues = []
        self._is_productive = None


    @property
    def is_productive(self):
        if self._is_productive is None:
            try:
                problems = any([self.stop_codons(self.antibody),
                                self.ambig_codons(self.antibody),
                                self.vdj_disagreement(self.antibody),
                                self.missing_conserved_junc_residues(self.antibody),
                                self.junction_not_in_frame(self.antibody),
                                self.out_of_frame_indels(self.antibody),
                                ])
                if problems:
                    self._is_productive = False
                else:
                    self._is_productive = True
            except:
                self.antibody.exception('PRODUCTIVITY ERROR', traceback.format_exc())
                raise
        return self._is_productive


    def stop_codons(self, antibody):
        if '*' in antibody.vdj_aa:
            self.productivity_issues.append('Contains stop codon(s)')
            return True
        return False


    def ambig_codons(self, antibody):
        if 'X' in antibody.vdj_aa.upper():
            self.productivity_issues.append('Contains ambiguous codon(s)')
            return True
        return False


    def vdj_disagreement(self, antibody):
        if antibody.v.full[:3] != antibody.j.full[:3]:
            self.productivity_issues.append('V-gene ({}) and J-gene({}) disagree'.format(antibody.v.gene, antibody.j.gene))
            return True
        if antibody.d:
            if not antibody.d.gene:
                return False
            if antibody.v.full[:3] != antibody.d.full[:3]:
                self.productivity_issues.append('V-gene ({}) and D-gene({}) disagree'.format(antibody.v.gene, antibody.d.gene))
                return True
        return False


    def missing_conserved_junc_residues(self, antibody):
        if not antibody.junction:
            return True
        start = 'C'
        end = 'W' if antibody.chain == 'heavy' else 'F'
        if antibody.junction.junction_aa[0] != start or antibody.junction.junction_aa[-1] != end:
            self.productivity_issues.append('Junction ({}) lacks conserved start and/or end residue'.format(antibody.junction.junction_aa))
            return True
        return False


    def junction_not_in_frame(self, antibody):
        if not antibody.junction.in_frame:
            self.productivity_issues.append('Junction ({}) is out of frame'.format(antibody.junction.junction_nt))
            return True
        return False


    def out_of_frame_indels(self, antibody):
        problems = 0
        for indels in [antibody.v.insertions, antibody.j.insertions, antibody.v.deletions, antibody.j.deletions]:
            for i in indels:
                if i.in_frame == 'no':
                    self.productivity_issues.append('Contains a frameshift {}: {}'.format(i.type, i.abstar_formatted))
                    problems += 1
        return False if problems == 0 else True
