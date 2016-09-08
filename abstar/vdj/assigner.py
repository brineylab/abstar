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


import abc
import os

from abtools.log import get_logger


class AbstractAssigner(object):
    """
    docstring for AbstractAssigner
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        super(AbstractAssigner, self).__init__()
        self.logger = get_logger('assigner')
        self._assigned = None
        self._unassigned = None
        self._germline_directory = None
        self._binary_directory = None


    @abc.abstractmethod
    def __call__(self, sequence_file, species):
        pass


    @property
    def germline_directory(self):
        if self._germline_directory is None:
            mod_dir = os.path.dirname(os.path.abspath(__file__))
            self._germline_directory = os.path.join(mod_dir, 'germline_dbs')
        return self._germline_directory


    @property
    def binary_directory(self):
        if self._binary_directory is None:
            mod_dir = os.path.dirname(os.path.abspath(__file__))
            self._binary_directory = os.path.join(mod_dir, 'bin')
        return self._binary_directory

    @property
    def assigned(self):
        if self._assigned is None:
            self._assigned = []
        return self._assigned

    @assigned.setter
    def assigned(self, assigned):
        self._assigned = assigned

    @property
    def unassigned(self):
        if self._unassigned is None:
            self._unassigned = []
        return self._unassigned

    @unassigned.setter
    def unassigned(self, unassigned):
        self._unassigned = unassigned
