#!/usr/bin/python
# filename: vdj.py

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

import logging
import os
import sys
import tempfile
import traceback

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from celery.utils.log import get_task_logger

from abstar.utils import blast
from abstar.utils.queue.celery import celery

from abtools import log
from abtools.alignment import local_alignment
from abtools.sequence import Sequence


class VDJ(object):
    '''
    Data structure for managing germline gene assignments for a single antibody sequence.

    Input is a Sequence object, V and J BlastResult objects, and (optionally) a
    DiversityResult object.

    VDJ objects without an identified V or J gene will have 'None' for the unidentified gene
    segment and the 'rearrangement' attribute will be False.
    '''

    def __init__(self, seq, args, v, j, d=None):
        super(VDJ, self).__init__()
        self.id = seq.id
        self.raw_input = seq.sequence
        self.raw_query = seq._input_sequence
        self.args = args
        self.strand = v.strand
        self._uaid = None
        self._isotype = None
        self.v = v
        self.j = j
        self.d = d
        if not self._valid_assignments():
            self.rearrangement = False
        elif all([v is not None, j is not None]):
            try:
                self.rearrangement = True
                logger.debug('V-BLAST INPUT: {}'.format(v.seq.sequence))
                logger.debug('V-BLAST QUERY: {}'.format(v.query_alignment))
                logger.debug('J-BLAST INPUT: {}'.format(j.seq.sequence))
                logger.debug('J-BLAST QUERY: {}'.format(j.query_alignment))
                self._get_attributes()
            except:
                self.rearrangement = False
                logger.debug('VDJ ATTRIBUTE ERROR: {}'.format(seq.id))
                logger.debug(traceback.format_exc())
        else:
            if v is None:
                logger.debug('V ASSIGNMENT ERROR: {}'.format(seq.id))
            if j is None:
                logger.debug('J ASSIGNMENT ERROR: {}'.format(seq.id))
            self.rearrangement = False

    @property
    def uaid(self):
        if self._uaid is None:
            if self.args.uaid == 0:
                self._uaid = ''
            elif self.args.uaid > 0:
                self._uaid = self.raw_query[:self.args.uaid]
            elif self.args.uaid < 0:
                self._uaid = self.raw_query[self.args.uaid:]
            else:
                self._uaid = ''
        return self._uaid

    @property
    def isotype(self):
        return self._isotype


    def _valid_assignments(self):
        if any([self.v is None, self.j is None]):
            return False
        if self.v.chain != self.j.chain:
            logger.debug("ASSIGNMENT ERROR: {} - V and J chains don't match".format(self.id))
            return False
        return True

    def _get_attributes(self):
        'Adds VDJ attributes, only for VDJ objects with an identified V and J gene.'
        self.species = self.v.species
        self.chain = self._get_chain()
        self.query_reading_frame = self._query_reading_frame()
        self.v_start = self._get_v_start()
        self.v_end = self._get_v_end()
        self.j_start = self._get_j_start()
        self.j_end = self._get_j_end()
        if self.chain in ['kappa', 'lambda'] and self.j.regions.fix_v_overlap:
            self._fix_v_overlap()
        if self.d:
            self.d_start = self._get_d_start()
            self.d_end = self._get_d_end()
            self.n1_start = self._get_n1_start()
            self.n1_end = self._get_n1_end()
            self.n2_start = self._get_n2_start()
            self.n2_end = self._get_n2_end()
            self.n_start = None
            self.n_end = None
        else:
            self.d_start = None
            self.d_end = None
            self.n1_start = None
            self.n1_end = None
            self.n2_start = None
            self.n2_end = None
            self.n_start = self._get_n_start()
            self.n_end = self._get_n_end()

        self.gapped_vdj_nt = self._gapped_vdj_nt()
        self.vdj_nt = self._vdj_nt()
        self.gapped_vdj_region_string = self._get_gapped_vdj_region_string()
        self.vdj_region_string = self._get_vdj_region_string()
        self.vdj_aa = self._vdj_aa()
        self.junction = self._get_junction()
        self.gapped_vdj_germ_nt = self._gapped_vdj_germ_nt()
        self.vdj_germ_nt = self.gapped_vdj_germ_nt.replace('-', '')
        self.vdj_germ_aa = self._vdj_germ_aa()
        self.germ_junction = self._get_junction(germ=True)
        self.codons = self._get_codons()
        self.gapped_codons = self._get_codons(gapped=True)
        self.productive = self._check_productivity()
        if self.args.isotype:
            from abstar.utils.isotype import get_isotype
            self._isotype = get_isotype(self)
        else:
            self._isotype = ''
        if not self.junction:
            self.rearrangement = False

    def _get_chain(self):
        'Returns the antibody chain.'
        try:
            if self.v.top_germline.startswith('IGH'):
                return 'heavy'
            elif self.v.top_germline.startswith('IGK'):
                return 'kappa'
            elif self.v.top_germline.startswith('IGL'):
                return 'lambda'
        except:
            logger.debug('GET CHAIN ERROR: {}, {}'.format(self.id,
                                                          self.v.top_germline))
            return ''

    def _query_reading_frame(self):
        'Returns the reading frame of the query sequence.'
        try:
            return self.v.germline_start % 3
        except:
            logger.debug('QUERY READING FRAME ERROR: {}, {}'.format(self.id,
                                                                    self.v.germline_start))

    def _fix_v_overlap(self):
        logger.debug('FIXING V_OVERLAP: {}'.format(self.id))
        nt = self.j.regions.v_overlap_length
        # fix the query sequence
        new_fr4_nt = self.v.query_alignment[-nt:] + self.j.regions.nt_seqs['FR4']
        logger.debug('OVERLAP LENGTH: {}'.format(nt))
        logger.debug('OVERLAP SEQUENCE: {}'.format(self.v.query_alignment[-nt:]))
        logger.debug('RAW GERMLINE: {}'.format(self.j.germline_seq))
        logger.debug('OLD FR4: {}'.format(self.j.regions.nt_seqs['FR4']))
        logger.debug('NEW FR4: {}'.format(new_fr4_nt))
        self.j.regions.nt_seqs['FR4'] = new_fr4_nt
        self.j.regions.aa_seqs['FR4'] = str(Seq(new_fr4_nt, generic_dna).translate())
        # fix the germline sequence
        new_fr4_germ_nt = self.v.germline_alignment[-nt:] + self.j.regions.germline_nt_seqs['FR4']
        self.j.regions.germline_nt_seqs['FR4'] = new_fr4_germ_nt
        self.j.regions.germline_aa_seqs['FR4'] = str(Seq(new_fr4_germ_nt, generic_dna).translate())

    def _gapped_vdj_nt(self):
        try:
            gapped_vdj_nt = self.v.query_alignment + \
                self.j.input_sequence[:self.j.query_start] + \
                self.j.query_alignment
            return gapped_vdj_nt
        except:
            logger.debug('VDJ NT ERROR: {}, {}'.format(self.id, self.raw_query))

    def _vdj_nt(self):
        'Returns the nucleotide sequence of the VDJ region.'
        return self.gapped_vdj_nt.replace('-', '')

    def _vdj_aa(self):
        'Returns the amino acid sequence of the VDJ region.'
        offset = (self.query_reading_frame * 2) % 3
        trim = len(self.vdj_nt) - (len(self.vdj_nt[offset:]) % 3)
        translated_seq = Seq(self.vdj_nt[offset:trim], generic_dna).translate()
        return str(translated_seq)

    def _gapped_vdj_germ_nt(self):
        'Returns the nucleotide sequence of the VDJ region.'
        try:
            if self.d:
                germ_junction = self.junction.n1_nt + \
                    self.d.germline_alignment + \
                    self.junction.n2_nt
            else:
                germ_junction = self.junction.n_nt
            germ_vdj_nt = self.v.germline_alignment + \
                germ_junction + self.j.germline_alignment
            return germ_vdj_nt
        except:
            logger.debug('VDJ GERM NT ERROR: {}, {}'.format(self.id, self.raw_query))
            logger.debug(traceback.format_exc())

    def _vdj_germ_aa(self):
        'Returns the amino acid sequence of the VDJ region.'
        offset = (self.query_reading_frame * 2) % 3
        trim = len(self.vdj_germ_nt) - (len(self.vdj_germ_nt[offset:]) % 3)
        translated_seq = Seq(self.vdj_germ_nt[offset:trim], generic_dna).translate()
        return str(translated_seq)

    def _get_junction(self, germ=False):
        from abstar.utils import junction
        return junction.get_junction(self, germ=germ)

    def _get_v_start(self):
        return self.v.query_start

    def _get_v_end(self):
        return len(self.v.query_alignment)

    def _get_j_start(self):
        return self.v_end + self.j.query_start

    def _get_j_end(self):
        return self.j_start + len(self.j.query_alignment)

    def _get_n1_start(self):
        return self.v_end

    def _get_n1_end(self):
        return self.d_start

    def _get_d_start(self):
        return self.v_end + self.d.query_start

    def _get_d_end(self):
        return self.d_start + len(self.d.query_alignment)

    def _get_n2_start(self):
        return self.d_end

    def _get_n2_end(self):
        return self.j_start

    def _get_n_start(self):
        return self.v_end

    def _get_n_end(self):
        return self.j_start

    def _get_vdj_region_string(self):
        rstring = ''
        for s, r in zip(self.gapped_vdj_nt, self.gapped_vdj_region_string):
            if s == '-':
                continue
            rstring += r
        return rstring

    def _get_gapped_vdj_region_string(self):
        region_string = ''
        region_string += 'V' * self.v_end
        if self.d:
            region_string += 'N' * (self.d_start - self.v_end)
            region_string += 'D' * (self.d_end - self.d_start)
            region_string += 'N' * (self.j_start - self.d_end)
        else:
            region_string += 'N' * (self.j_start - self.v_end)
        region_string += 'J' * (self.j_end - self.j_start)
        return region_string

    def _get_codons(self, gapped=False):
        from abstar.utils import codons
        return codons.parse_codons(self, gapped=gapped)

    def _check_productivity(self):
        from abstar.utils import productivity
        return productivity.check_productivity(self)

    def _build_output(self, output_type):
        from abstar.utils import output
        return output.build_output(self, output_type)


@celery.task
def run(seq_file, output_dir, arg_dict):
    '''
    Wrapper function to multiprocess (or not) the assignment of V, D and J
    germline genes. Also writes the JSON-formatted output to file.

    Input is a a FASTA-formatted file of antibody sequences and the output directory.
    Optional input items include the species (supported species: 'human'); length of
    the unique antibody identifier (UAID); and debug mode (which forces single-threading
    and prints more verbose errors.)

    Output is the number of functional antibody sequences identified in the input file.
    '''
    try:
        from abstar.abstar import Args
        args = Args(**arg_dict)
        # global logger
        # if args.cluster:
        #     logger = get_task_logger(__name__)
        # else:
        #     logger = log.get_logger(__name__)
        output_filename = os.path.basename(seq_file)
        if args.output_type == 'json':
            output_file = os.path.join(output_dir, output_filename + '.json')
        elif args.output_type in ['imgt', 'hadoop']:
            output_file = os.path.join(output_dir, output_filename + '.txt')
        vdj_output = process_sequence_file(seq_file, args)
        if not vdj_output:
            return None
        clean_vdjs = [vdj for vdj in vdj_output if vdj.rearrangement]
        output_count = write_output(clean_vdjs, output_file, args.output_type, args.pretty, args.padding)
        return (output_file, output_count)
    except:
        # logger.debug(traceback.format_exc())
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))
        # run.retry(exc=exc, countdown=5)


def write_output(output, outfile, output_type, pretty, padding):
    from abstar.utils.output import build_output
    logger.debug("Padding - {}\t Pretty - {}\t".format(padding, pretty))
    output_data = build_output(output, output_type, pretty, padding)
    open(outfile, 'w').write('\n'.join(output_data))
    if output_type in ['json', 'hadoop']:
        return len(output_data)
    else:
        return len(output_data) - 1


def process_sequences(sequences, args):
    seq_file = tempfile.NamedTemporaryFile(delete=False)
    seq_file.write('\n'.join([s.fasta for s in sequences]))
    seq_file.close()
    vdjs = process_sequence_file(seq_file.name, args)
    os.unlink(seq_file.name)
    return vdjs


def process_sequence_file(seq_file, args):
    '''
    Runs BLASTn to identify germline V, D, and J genes.

    Input is a Sequence object.

    Output is a list of VDJ objects.
    '''

    global logger
    logger = logging.getLogger()

    # Variable gene assignment
    vs = []
    seqs = [Sequence(s) for s in SeqIO.parse(open(seq_file, 'r'), 'fasta')]
    v_blast_records = blast.blast(seq_file, args.species, 'V')
    for seq, vbr in zip(seqs, v_blast_records):
        try:
            v = assign_germline(seq, vbr, args.species, 'V')
            if v.strand == 'minus':
                seq.sequence = seq.reverse_complement
            logger.debug('ASSIGNED V-GENE: {}, {}'.format(seq.id, v.top_germline))
        except Exception:
            logger.debug('V-GENE ASSIGNMENT ERROR: {}'.format(seq.id))
            logger.debug(traceback.format_exc())
            v = None
        finally:
            if v:
                # v_end = v.query_end + v.query_start + 1
                v_end = v.query_end + 1
                if len(seq.sequence[v_end:]) <= 10:
                    logger.debug('REMANING REGION TOO SHORT AFTER V-ALIGNMENT REMOVAL: {}'.format(v.id))
                    # v = None
                    logger.debug('PERFORMING SECOND V-GENE REALIGNMENT: {}'.format(v.id))
                    v.realign_variable(v.top_germline,
                        match=30, mismatch=-20,
                        gap_open_penalty=220, gap_extend_penalty=1)
                    v.annotate()
                    # v_end = v.query_end + v.query_start + 1
                    v_end = v.query_end + 1
                    if len(seq.sequence[v_end:]) <= 10:
                        logger.debug('REMANING REGION TOO SHORT AFTER SECOND V-ALIGNMENT REMOVAL: {}'.format(v.id))
                        v = None
            vs.append((seq, v))
    v_blast_results = [v[1] for v in vs if v[1]]
    seqs = [v[0] for v in vs if v[1]]
    failed_seqs = [v[0] for v in vs if not v[1]]
    logger.debug('V-ASSIGNMENT RESULT: for {}, {} of {} sequences failed v-gene assignment'.format(
        os.path.basename(seq_file), len(failed_seqs), len(seqs) + len(failed_seqs)))
    for fs in failed_seqs:
        logger.debug('NO V-GENE ASSIGNMENT: {}'.format(fs.id))
    if not v_blast_results:
        seq_filename = os.path.basename(seq_file)
        logger.debug('NO VALID REARRANGEMENTS IN FILE: {}'.format(seq_filename))
        return None

    # Joining gene assignment
    js = []
    j_blastin = blast.build_j_blast_input(seqs, v_blast_results)
    j_blast_records = blast.blast(j_blastin.name, args.species, 'J')
    j_seqs = [Sequence(s) for s in SeqIO.parse(open(j_blastin.name, 'r'), 'fasta')]
    for i, jdata in enumerate(zip(j_seqs, j_blast_records)):
        j_seq, jbr = jdata
        try:
            j = assign_germline(j_seq, jbr, args.species, 'J')
            if not j:
                logger.debug('NO ASSIGNED J-GENE: {}'.format(j_seq.id))
            logger.debug('ASSIGNED J-GENE: {}, {}'.format(j_seq.id, j.top_germline))
        except:
            logger.debug('J-GENE ASSIGNMENT ERROR: {}'.format(j_seq.id))
            logger.debug(traceback.format_exc())
            try:
                vbr = v_blast_records[i]
                vseq = seqs[i]
                logger.debug('J-GENE ASSIGNMENT ERROR: {}\n{}'.format(j_seq.id,
                                                                          vseq.sequence, ))
                                                                          # vbr.query_alignment))
            except:
                logger.debug('J-GENE ASSIGNMENT ERROR: {}, could not print query info'.format(j_seq.id))
                logger.debug(traceback.format_exc())
            j = None
        finally:
            js.append(j)
    os.unlink(j_blastin.name)

    # Build VDJ objects (including optional D gene assignment)
    vdjs = []
    for seq, v, j in zip(seqs, v_blast_results, js):
        try:
            if not v or not j:
                continue
            if v.chain == 'heavy':
                junc_start = len(v.query_alignment) + v.query_start
                junc_end = junc_start + j.query_start
                junction = Sequence(seq.sequence[junc_start:junc_end], id=seq.id)
                if junction:
                    d = assign_d(junction, args.species)
                    if d is not None:
                        logger.debug('ASSIGNED D-GENE: {}, {}'.format(seq.id, d.top_germline))
                    vdjs.append(VDJ(seq, args, v, j, d))
                    continue
                vdjs.append(VDJ(seq, args, v, j))
            else:
                vdjs.append(VDJ(seq, args, v, j))
        except:
            logger.debug('VDJ ERROR: {}'.format(seq.id))
            logger.debug(traceback.format_exc())
    return vdjs


def assign_germline(seq, blast_record, species, segment):
    '''
    Identifies germline genes for a given antibody sequence (seq).

    Input is a Sequence object, the species of origin, the gene
    segment to be assigned (options are 'V' or 'J') and, optionally, the
    starting and ending points of the possible germline gene location.

    Output is a BlastResult object.
    '''
    blast_result = blast.BlastResult(seq, blast_record, species)
    if blast_result.gene_type == 'variable':
        blast_result.realign_variable(blast_result.top_germline)
    blast_result.annotate()
    return blast_result


def assign_d(seq, species):
    '''
    Identifies the germline diversity gene for a given sequence.
    Alignment is performed using the ssw_wrap.Aligner.align function.

    Input is a junction sequence (as a string) and the species of origin.

    Output is a DiversityResult object.
    '''
    mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    db_file = os.path.join(mod_dir, 'ssw/dbs/{}_D.fasta'.format(species.lower()))
    db_handle = open(db_file, 'r')
    germs = [Sequence(s) for s in SeqIO.parse(db_handle, 'fasta')]
    rc_germs = [Sequence(s.reverse_complement, id=s.id) for s in germs]
    germs.extend(rc_germs)
    db_handle.close()
    alignments = local_alignment(seq, targets=germs,
                                 gap_open_penalty=20, gap_extend_penalty=2)
    alignments.sort(key=lambda x: x.score, reverse=True)
    try:
        return blast.DiversityResult(seq, alignments[:5])
    except IndexError:
        return None
