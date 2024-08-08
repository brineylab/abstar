# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


import logging
import os
from typing import Optional

import abutils

from ..annotation.germline import get_germline_database_path


class AssignerBase:
    """
    docstring for AssignerBase
    """

    def __init__(
        self,
        output_directory: str,
        germdb_name: str,
        receptor: str,
        log_directory: str,
        logger: Optional[logging.Logger] = None,
        debug: bool = False,
    ):
        self.output_directory = os.path.abspath(output_directory)
        self.log_directory = os.path.abspath(log_directory)
        self.receptor = receptor
        self.germdb_name = germdb_name
        self.debug = debug
        self.logger = logger if logger is not None else abutils.log.null_logger()
        self.germdb_path = get_germline_database_path(germdb_name, receptor)
        self.to_delete = []  # files that should be deleted during cleanup

        abutils.io.make_dir(self.output_directory)

    def __call__(self, sequence_file: str) -> str:
        """
        All classes subclassing ``AssignerBase`` must implement this method.

        The only expected argument is `sequence_file`, which is a single FASTA or FASTQ-formatted
        file containing input sequences. Gzipped files are supported.

        The return value should be the path to a parquet-formatted file containing the following
        columns:
          * ``sequence_id``: a unique identifier for each input sequence.
          * ``sequence``: the input sequence.
          * ``quality``: the quality score for each input sequence.
          * ``rev_comp``: a boolean indicating whether the input sequence is in the reverse_complement orientation
          * ``v_call``: the V-gene call for each input sequence.
          * ``v_support``: the V-gene evalue for each input sequence.
          * ``d_call``: the D-gene call for each input sequence.
          * ``d_support``: the D-gene evalue for each input sequence.
          * ``j_call``: the J-gene call for each input sequence.
          * ``j_support``: the J-gene evalue for each input sequence.

        """
        raise NotImplementedError("Subclasses must implement this method")

    def cleanup(self):
        abutils.io.delete_files(self.to_delete)
        self.to_delete = []

    # @staticmethod
    # def get_germdb_path(germdb_name: str, receptor: str = "bcr") -> str:
    #     """
    #     Get the path to a germline database.

    #     Parameters
    #     ----------
    #     germdb_name : str
    #         The name of the germline database.

    #     receptor : str, default: "bcr"
    #         The receptor type. Options are "bcr" and "tcr".

    #     Returns
    #     -------
    #     str
    #         The path to the germline database.
    #     """
    #     germdb_name = germdb_name.lower()
    #     receptor = receptor.lower()
    #     if receptor not in ["bcr", "tcr"]:
    #         raise ValueError(f"Receptor type {receptor} not supported")

    #     # check the addon directory first
    #     addon_dir = os.path.expanduser(f"~/.abstar/germline_dbs/{receptor}")
    #     if os.path.isdir(addon_dir):
    #         if germdb_name.lower() in [
    #             os.path.basename(d[0]) for d in os.walk(addon_dir)
    #         ]:
    #             return os.path.join(addon_dir, f"{receptor}/{germdb_name}")

    #     # if a user-generated DB isn't found, use the built-in DB
    #     curr_dir = os.path.dirname(os.path.abspath(__file__))
    #     germdb_path = os.path.join(curr_dir, f"germline_dbs/{receptor}/{germdb_name}")
    #     if not os.path.exists(germdb_path):
    #         raise FileNotFoundError(
    #             f"Germline database {germdb_name} for receptor {receptor} not found"
    #         )
    #     return germdb_path


##
##
##
##
##
##


# class BaseAssigner:
#     """
#     ``BaseAssigner`` provides an abstract base class for custom VDJ Assigners.

#     Overview
#     --------

#     Several different tools exist for inferring the germline components of recombined
#     antibody sequences. Due to their differing design goals, these tools often have
#     different priorities -- some are optimized for speed with a small (often
#     negligible) accuracy penalty, while others define recombination junctions with high
#     precision at the cost of speed. Because there isn't a single germline assignment tool
#     that is optimal for every use case, the germline assignment component of AbStar has been
#     designed to be modular and replaceable, allowing the relatively straightforward addition
#     of new germline assignment tools and the selection of the desired Assigner at runtime.
#     AbStar's default Assigner (``blastn``) identifies V- and J-genes using BLASTn and identifies
#     D-genes using rapid Smith-Waterman local alignment.

#     The basic purpose of any custom Assigner class is to accept sequences and
#     produce ``VDJ`` objects. ``VDJ`` objects package the input sequence together with
#     ``GermlineSegment`` objects that contain information about V, D and/or J assignments.
#     Following germline assignment by the Assigner, addtional annotation will be performed
#     by AbStar. The ``VDJ`` objects provide a common interface by which various Assigners
#     can communicate with the additional sequence annotation components in AbStar. This additional
#     annotation is consistent regardless of the Assigner (as is the output schema), and provides
#     a unified method by which different germline Assigners can be used across different projects
#     while maintaining schematic compatibility with downstream analysis tools.

#     ``BaseAssigner`` provides the following attributes:

#       - ``assigned``: a list designed to contain ``VDJ`` objects for which
#                     successful V(D)J assignment has been completed. These
#                     ``VDJ`` objects will be additionally annotated by AbStar
#                     following assignment. Exception/log information in these
#                     ``VDJ`` objects will only be written to the log file if AbStar
#                     is run in debug mode.

#       - ``unassigned``: a list designed to contain ``VDJ`` objects for which
#                       V(D)J assigjment was unsuccessful. These ``VDJ`` objects
#                       will not be additionally annotated by AbStar and any exception/log
#                       information contained in these ``VDJ`` objects will be logged
#                       (even if AbStar was not run in debug mode).

#       - ``germline_directory``: path to the directory containing germline databases.
#                               Germline DBs are heirarchically grouped into separate folders
#                               first by species, then by the Assigner for which they were
#                               designed (``blastn``, for example). Within each folder, germline
#                               DB files are named as ``{segment}.fasta``, so the Blast V-gene
#                               DB for human would be ``human/blast/v.fasta``. Note that this will
#                               change as AbStar is updated to annotate TCRs in addition to antibodies

#       - ``binary_directory``: path to the directory containing Assigner binaries. Binaries
#                             are named as ``{binary}_{system}``, where ``system`` is the
#                             lowercase output from platform.system(). For example, the
#                             Blastn binary for OSX would be at ``binary_directory/blastn_darwin``.
#                             Some assigners require more than a simple binary for execution (Partis,
#                             for example). In this case, all required elements will be placed into
#                             a folder within the binary directory, named using the same convention
#                             as binaries. So the folder containing Linux-compiled Partis components
#                             will be ``binary_directory/partis_linux``. The main Partis executable
#                             would then be located at ``binary_directory/partis_linux/bin/partis``.

#     In order to build a custom assignment class, you simply need to subclass BaseAssigner, implement
#     your assigner, and register your assigner in abstar.assigners.registry.ASSIGNERS. Obviously, this
#     is in addition to adding any required binaries and specialized germline database formats. AbStar
#     prefers to incorporate assigner binaries if at all possible (as opposed to relying on system
#     installations) to smooth the process of installing AbStar and eliminate the need for users to
#     install several additional programs before having a fully operational AbStar installation.

#     There are two critical elements to the design of your assigner class. First, ensure that
#     you call the BaseAssigner's ``__init__()`` at the start of your class's ``__init__()``. Second,
#     you must implement the  ``__call__()`` method.


#     ``__init__()``
#     --------------

#     Your custom assigner class's ``__init__()`` must accept a single argument (``species``, in
#     addition to the required ``self``). You must then call the BaseAssigner's ``__init__()`` and pass
#     the ``species`` argument. The simplest was to do this is by using ``super()``, like this:
#     ``super(MyAssigner, self).__init__(species)``, where ``MyAssigner`` is your assigner class.


#     Implementing the ``__call__()`` method
#     --------------------------------------

#       - ``__call__()`` must accept two arguments. The first (``sequence_file``) is a FASTA- or
#         FASTQ-formatted file of input sequences, depending on the user-provided input. The second argument
#         (``file_format``) defines the format of ``sequence_file``, either ``'fasta'`` or ``'fastq'``. If
#         a FASTQ file is required by your assigner, ensure that ``__call__()`` raises the appropriate exception
#         (and provides a sufficiently clear description of the problem) if a FASTA-formatted file is provided.

#       - The result of ``__call__()`` should be the generation of a ``VDJ`` object for each input sequence
#         (more on this below). The ``VDJ`` objects should be appended to either the ``assigned`` or ``unassigned``
#         attribute, depending on the outcome of the V(D)J assignment. A high-level overview of how ``__call__()``
#         should work is as follows::

#             def __call__(self, seqs, species):
#                 for seq in seqs:
#                     vdj = self.assign_vdj(seq, species)
#                     if self.is_properly_assigned(vdj):
#                         self.assigned.append(vdj)
#                     else:
#                         self.unassigned.append(vdj)

#         In this case, the output of ``assign_vdj()`` should be a ``VDJ`` object.


#     V(D)J assignment
#     ----------------

#     AbStar defines ``VDJ`` objects, which are designed to package the input sequence and ``GermlineSegment``
#     objects that provide information about germline assignment for V, D and/or J genes. ``VDJ`` objects must
#     contain the input sequence and may also contain information about V, D and J assignments. While the
#     sequence is required, the germline assignments are optional, to allow creation of ``VDJ`` objects for
#     light chains (which lack a D-gene) and/or for unsuccessful assignments. Additionally, this allows
#     flexibility in designing the ``__call__()`` method. Assigners can either instantiate a ``VDJ`` object at
#     the start of ``__call__()`` and attach ``GermlineSegment``s as they are computed or they can instantiate
#     ``VDJ`` objects upon completion of V(D)J assignment using the input sequences and computed
#     ``GermlineSegment``s as instantiation arguments. Because ``VDJ`` objects contain built-in logging
#     capability, such objects should be created for each input sequence even in the case of unsuccessful
#     assignments so that assignment issues can be logged and reported.

#     ``VDJ`` objects include a logging method: ``log()``. Using ``log()`` is similar to using Python's
#     builtin ``print()`` function. Multiple arguments to ``log()`` will be space-delimited, and a ``sep``
#     keyword allows you to provide an alternate delimiter. Log messages are typically used to record useful
#     information or to note the occurance of errors/exceptions that are caught and handled by the Assigner.
#     It is also a useful means by which to note the reason a specific input sequence is unable to be fully
#     processed (because it isn't a valid antibody sequence, etc). By convention, AbStar log entries should
#     contain the name of the log event (in all caps) followed by a colon and additional log information.
#     An example ``log()`` call in the case of a V(D)J assignment that was aborted because the V-gene alignment
#     score didn't meet the alignment quality threshold::

#         for seq in seqs:
#             vdj = VDJ(seq)
#             v = assign_v_gene(seq)
#             if v.score < threshold:
#                 vdj.log('V-GENE ASSIGNMENT ERROR:', 'Alignment score ({}) was too low'.format(v.score))
#                 continue

#     ``VDJ`` objects also contain a mechanism for logging exceptions: ``exception()``. This is typically used
#     to log unexpected exceptions (ie, those that are not caught and handled by the Assigner). It is useful
#     to provide some information about when/where in the Assigner the exception occured, as well as a formatted
#     version of the exception traceback. An example of ``exceptions()`` being used to log an exception that
#     occurs during V-gene assignment::

#         vdj = VDJ(seq)
#         try:
#             v = assign_v_gene(seq)
#         except:
#             vdj.exception('V-GENE ASSIGNMENT EXCEPTION:', traceback.format_exc())

#     One of the primary benefits of recording log and exception information through the ``VDJ`` objects is
#     that it allows for much easier formatting of such events when multiple sequences are being processed
#     in parallel. Using Python's ``logging`` utility would result in log events being logged as they occur,
#     which means that multiple log entries for the same input sequence may be far apart in the resulting log
#     file. Packaging log and exception information with each ``VDJ`` object alows AbStar to create a much
#     more useful log file in which all information related to a single input sequence is found in a single
#     log location.


#     ``GermlineSegment`` objects
#     ---------------------------

#     ``GermlineSegment`` objects contain information about the assigned germline gene segment (V, D or J).
#     Instantiation of a ``GermlineSegment`` object requires only the name of the germline gene segment (in
#     IMGT format, like 'IGHV3-23*01') and the database name. Several other optional arguments can be included
#     at instantiation. ``score`` is the assignment score, typically an ``int`` or ``float``. ``strand``
#     indicates the strand orientation of the input sequence (``'+'`` or ``'-'``). ``others`` is a list
#     of additional high scoring ``GermlineSegment`` objects. ``assigner_name`` is the name of the custom
#     Assigner (as a ``str``), and will be converted to lowercase before recording in AbStar output.

#     ``GermlineSegment`` objects also provide a few additional public attributes that may be desirable
#     for certain custom Assigners, especially those that have been specifically designed to accurately
#     define recombination regions. ``query_start``, ``query_end``, ``germline_start`` and ``germline_end``
#     allow the Assigner to define the start and end points of the germline alignment for both the query
#     and germline sequences. Note that both ``query_start`` and ``query_start`` and ``query_end`` MUST
#     be relative to the complete query sequence in the ``'+'`` orientation (in other words, V-gene at
#     the 5' end and J-gene at the 3' end). AbStar realigns all assigned germline genes using alignment
#     parameters designed to accurately identify somatic mutations as well as somatic mutation-induced
#     indels. If the additional assignment position attributes are provided by the Assigner, AbStar will
#     force their use during realignment. If not, AbStar will use the assigned germline genes for realignment
#     and will define the start and end points of the germline alignment during realignment.


#     Example
#     -------

#     The following is an example of a custom Assigner class, named MyAssigner. It uses BioPython's SeqIO
#     to parse the input file, and because BioPython is a dependency of AbStar, you can assume that BioPython
#     will be present on any system running AbStar. All custom Assigners should be located in the ``assigners``
#     directory::

#         import os
#         import platform
#         import traceback

#         from abstar.assigners.assigner import BaseAssigner
#         from abstar.core.vdj import VDJ
#         from abstar.core.germline import GermlineSegment

#         from Bio import SeqIO

#         class MyAssigner(BaseAssigner):

#             def __init__(self, db_name):
#                 super(MyAssigner, self).__init__(db_name)
#                 self.binary = os.path.join(self.binary_directory, 'mybinary_{}'.format(platform.system()))

#             def __call__(self, sequence_file, file_format):
#                 for sequence in SeqIO.parse(open(sequence_file, 'r'), file_format):
#                     seq = Sequence(sequence)
#                     vdj = VDJ(seq)

#                     # V-gene assignment
#                     # The VDJ object (rather than just the sequence) is passed to the assignment
#                     # function for three reasons:
#                     #    1) it contains the sequence, so we don't need to pass it separately
#                     #    2) it allows logging during germline assignment, via vdj.log()
#                     #    3) subsequent assignments will have access to information about
#                     #       previous assignments (J-gene assignment operations have access
#                     #       to information about the V-gene assignment, etc.)
#                     try:
#                         vdj.v = self.assign_germline(vdj, 'V')
#                     except:
#                         vdj.exception('V-ASSIGNMENT:', traceback.format_exc())

#                     # J-gene assignment
#                     try:
#                         vdj.j = self.assign_germline(vdj, 'J')
#                     except:
#                         vdj.exception('J-ASSIGNMENT:', traceback.format_exc())

#                     # D-gene assignment
#                     try:
#                         vdj.d = self.assign_germline(vdj, 'D')
#                     except:
#                         vdj.exception('D-ASSIGNMENT:', traceback.format_exc())

#                     return vdj

#             def assign_germline(self, vdj, segment):
#                 germ_db = os.path.join(self.germline_directory, '{}/ungapped/{}.fasta'.format(self.db_name,
#                                                                                               segment.lower()))

#                 # do stuff to assign the germline gene, using the species-appropriate germline DB
#                 # ...
#                 # ...
#                 # ...
#                 # germs is a list of germline genes, ordered by score (highest first)

#                 if germs[0].score < threshold:
#                     vdj.log('{}-ASSIGNMENT ERROR:'.format(segment),
#                             'Score ({}) is too low'.format(germs[0].score))
#                     return None
#                 others = [GermlineSegment(germ.name, self.db_name, score=germ.score) for germ in germs[1:6]]
#                 return GermlineSegment(germs[0].name, self.db_name, score=germs[0].score, others=others)

#     """

#     def __init__(self, db_name, receptor):
#         super(BaseAssigner, self).__init__()
#         self.name = self.__class__.__name__.lower()
#         self.species = db_name
#         self.db_name = db_name
#         self.receptor = receptor.lower()
#         self._assigned = None
#         self._unassigned = None
#         self._germline_directory = None
#         self._binary_directory = None

#     def __call__(self, sequence_file: str, file_format: str):
#         # __call__ should be implemented by subclasses
#         pass

#     @property
#     def germline_directory(self):
#         """
#         returns a path to the directory containing the germline database for the species and receptor
#         """
#         if self._germline_directory is None:
#             self._germline_directory = get_germline_database_directory(
#                 self.db_name, self.receptor
#             )
#         return self._germline_directory

#     @germline_directory.setter
#     def germline_directory(self, directory):
#         self._germline_directory = directory

#     @property
#     def binary_directory(self):
#         """
#         returns a path to the directory containing the binary for the assigner
#         """
#         if self._binary_directory is None:
#             mod_dir = os.path.dirname(os.path.abspath(__file__))
#             self._binary_directory = os.path.join(mod_dir, "bin")
#         return self._binary_directory

#     @property
#     def assigned(self):
#         """
#         returns a list of VDJ objects that were successfully assigned
#         """
#         if self._assigned is None:
#             self._assigned = []
#         return self._assigned

#     @assigned.setter
#     def assigned(self, assigned):
#         self._assigned = assigned

#     @property
#     def unassigned(self):
#         """
#         returns a list of VDJ objects that were not successfully assigned
#         """
#         if self._unassigned is None:
#             self._unassigned = []
#         return self._unassigned

#     @unassigned.setter
#     def unassigned(self, unassigned):
#         self._unassigned = unassigned
