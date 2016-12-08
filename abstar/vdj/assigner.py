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


class BaseAssigner():
    """
    ``BaseAssigner`` provides an abstract base class for custom VDJ Assigners.

    Several different tools exist for inferring the germline components of recombined
    antibody sequences. Due to their differing design goals, these tools often prioritize
    certain aspects over others -- some are optimized for speed with a small (often
    negligible) accuracy penalty, while others define recombination junctions with high
    precision at the cost of speed. Because there isn't a single germline assignment tool
    that is optimal for every use case, the germline assignment component of AbStar has been
    designed to be modular and replaceable, allowing the relatively straightforward addition
    of new germline assignment tools. AbStar's default assigner identifies V- and J-genes
    using BLASTn and identifies D-genes using rapid Smith-Waterman alignment.

    A high-level overview of any custom Assigner class is that it accepts sequences and
    produces ``VDJ`` objects. ``VDJ`` objects package the input sequence together with
    ``GermlineSegment`` objects that contain information about V, D and/or J assignments.
    Following germline assignment by the assigner, addtional annotation will be performed
    by AbStar. This additional annotation is consistent regardless of the assigner (as is
    the output schema), and provides a unified format by which different germline assigners
    can be used across different projects while maintaining compatibility with downstream
    analysis tools.

    ``BaseAssigner`` provides the following attributes:

      - ``assigned``: a list designed to contain ``VDJ`` objects for which
                    successful V(D)J assignment has been completed. These
                    ``VDJ`` objects will be additionally annotated by AbStar
                    following assignment. Exception/log information in these
                    ``VDJ`` objects will only be written to the log file if AbStar
                    is run in debug mode.
      - ``unassigned``: a list designed to contain ``VDJ`` objects for which
                      V(D)J assigjment was unsuccessful. These ``VDJ`` objects
                      will not be additionally annotated by AbStar and any exception/log
                      information contained in these ``VDJ`` objects will be written
                      to the log file (even if AbStar was not run in debug mode).
      - ``germline_directory``: path to the directory containing germline databases.
                              Germline DBs are in separate folders, corresponding
                              to the Assigner for which they were designed (``blastn``,
                              for example). Within each folder, germline DB files are
                              named as ``{species}_gl_{segment}``, so the Blast V-gene DB for
                              human would be ``germline_directory/blast/human_gl_V``.
      - ``binary_directory``: path to the directory containing Assigner binaries. Binaries
                            are named as ``{binary}_{system}``, where ``system`` is the
                            lowercase output from platform.system(). For example, the
                            Blastn binary for OSX would be at ``binary_directory/blastn_darwin``.

    In order to build a custom assigner, you simply need to subclass BaseAssigner and
    implement both the ``__call__()`` method and the ``name`` property. The ``name`` property
    should be a string that describes the Assigner (for example, the built-in
    Blastn assigner is named ``'blastn'``). The ``__call__()`` method should perform V(D)J
    assignment.


    Implementing the ``__call__()`` method
    ----------------------------------
      - ``__call__()`` must accept wo arguments. The first (``seqs``) is a list of abtools.Sequence
        objects. Depending on the user-provided input (FASTA vs FASTQ), the Sequence objects may
        or may not contain quality information. If quoality information is not present, the ``qual``
        attribute of the Sequence object will be ``None``. If FASTQ files are required for your assigner,
        you should verify the presence of the ``qual`` attribute and raise an exception if necessary. The
        second argument passed to ``__call__()`` is ``species``, which is a string corresponding to the
        user-selected species from which the input sequences are derived. ``species`` may be used to
        select the appropriate germline database for germline assignment.
      - The result of ``__call__()`` should be the generation of a ``VDJ`` object for each input sequence
        (more on this below). The ``VDJ`` objects should be appended to either the ``assigned`` or ``unassigned``
        attribute, depending on the outcome of the V(D)J assignment. A high-level overview of how __call__()
        should work is as follows::

            def __call__(self, seqs, species):
                for seq in seqs:
                    vdj = self.assign_vdj(seq)
                    if self.is_properly_assigned(vdj):
                        self.assigned.append(vdj)
                    else:
                        self.unassigned.append(vdj)

        In this case, the output of `assign_vdj()` should be a `VDJ` object. `VDJ` objects must contain
        the input sequence and may also contain information about V, D and J assignments. While the
        sequence is required, the germline assignments are optional, to allow creation of `VDJ` objects for
        light chains (missing a D-gene) and/or for unsuccessful assignments. Because `VDJ` objects contain
        built-in logging, such objects should be created even in the case of unsuccessful assignments so
        that assignment issues can be logged and reported.


    ``VDJ`` objects
    ---------------


    """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        super(BaseAssigner, self).__init__()
        self._assigned = None
        self._unassigned = None
        self._germline_directory = None
        self._binary_directory = None


    @abc.abstractmethod
    def __call__(self, sequence_file, species):
        pass


    @abc.abstractproperty
    def name(self):
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
