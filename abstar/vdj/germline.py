#!/usr/bin/env python
# filename: germline.py

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


from Bio import SeqIO

from abtools.utils.properties import lazy_property


class Germline(object):
    """
    docstring for Germline

    Args:
    -----

        full (str): Full name of the germline gene, following IMGT naming conventions
            (eg: IGHV3-23*02)

        species (str): Species from which the germline gene was derived. Choices include
            'human', 'macaque', mouse', rabbit'.

        score (float): Score of the top germline alignment. Optional. If ``score`` is not
            provided, the realignment score will be used instead.

        strand (str): Strand of the alignment. Options are ``'+'`` for the positive strand
            (meaning the germline gene and the query sequence are in the same orientation)
            or ``'-'`` for the negative strand (query sequence is the reverse complement of
            the germline gene). Optional.

        others (list(Germline)): An optional list of additional high scoring germline genes.
            Can be a list of arbitrary length, with each member of the list being another
            ``Germline`` instance.
    """
    def __init__(self, full, species, score=None, strand=None, others=None):
        super(Germline, self).__init__()
        self.full = full
        self.species = species
        self.score = score
        self.strand = strand
        self.others = others
        self.gene_type = self.full[3].upper()
        self._family = None
        self._gene = None
        self._chain = None

        # properties that get assigned after re-alignment.
        # New assigners don't need to populate these.
        self.realignment = None
        self.raw_query = None
        self.raw_germline = None
        self.query_alignment = None
        self.germline_alignment = None
        self.alignment_midline = None
        self.alignment_length = None
        self.query_start = None
        self.query_end = None
        self.germline_start = None
        self.germline_end = None
        self.fs_indel_adjustment = 0
        self.nfs_indel_adjustment = 0
        self.has_insertion = 'no'
        self.has_deletion = 'no'
        self.insertions = None
        self.deletions = None
        self.regions = None
        self.nt_mutations = None
        self.aa_mutations = None


    @property
    def family(self):
        if self._family is None:
            sep = '-' if '-' in self.full else '*'
            self._family = self.full.split(sep)[0]
        return self._family

    @family.setter
    def family(self, family):
        self._family = family


    @property
    def gene(self):
        if self._gene is None:
            if all(['-' in self.full, '*' in self.full]):
                self._gene = self.full.split('*')[0]
        return self._gene

    @gene.setter
    def gene(self, gene):
        self._gene = gene


    @property
    def chain(self):
        if self._chain is None:
            c = {'H': 'heavy',
                 'K': 'kappa',
                 'L': 'lambda'}
            self._chain = c.get(self.full[2], None)
        return self._chain


    def process_realignment(self, aln):
        self.realignment = aln
        self.raw_query = aln.raw_query
        self.raw_germline = aln.raw_target
        self.query_alignment = aln.aligned_query
        self.germline_alignment = aln.aligned_target
        self.alignment_midline = ''.join(['|' if q == g else ' ' for q, g in zip(aln.aligned_query,
                                                                                 aln.aligned_target)])
        self.query_start = aln.query_begin
        self.query_end = aln.query_end
        self.germline_start = aln.target_begin
        self.germline_end = aln.target_end
        self._fix_ambigs()
        self._find_indels()


    def _get_germline_sequence_for_realignment(self):
        '''
        Identifies the appropriate germline variable gene from a database of all
        germline variable genes.

        Returns:
        --------

            str: Germline sequence, or ``None`` if the requested germline gene could not be found.
        '''
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        db_file = os.path.join(mod_dir, 'ssw/dbs/{}_{}.fasta'.format(self.species.lower(), self.gene))
        try:
            for s in SeqIO.parse(open(db_file), 'fasta'):
                if s.id == self.germ:
                    return str(s.seq)
            # TODO: log that the germline gene wasn't found in the database file
            return None
        except:
            # TODO: log that the germline database file couldn't be found
            return None


    def _fix_ambigs(self):
        '''
        Fixes ambiguous nucleotides by replacing them with the germline nucleotide.
        '''
        self.query_alignment = ''.join([q if q.upper() != 'N' else g for q, g in zip(self.query_alignment,
                                                                                     self.germline_alignment)])



    def _indel_check(self):
        if any(['-' in self.query_alignment, '-' in self.germline_alignment]):
            return True
        return False


    def _find_indels(self):
        '''
        Identifies and annotates indels in the query sequence.
        '''
        if self._indel_check():
            from abstar.utils import indels
            self.insertions = indels.find_insertions(self)
            if self.insertions:
                if 'yes' in [ins['in frame'] for ins in self.insertions]:
                    self.has_insertion = 'yes'
            self.deletions = indels.find_deletions(self)
            if self.deletions:
                if 'yes' in [deletion['in frame'] for deletion in self.deletions]:
                    self.has_deletion = 'yes'




def get_imgt_germlines(species, gene_type, gene=None):
    '''
    Returns one or more IMGTGermlineGene objects that each contain a single IMGT-gapped germline gene.

    Args:
    -----

        species (str): Species for which the germline genes should be obtained.

        gene_type (str): Options are 'V', 'D', and 'J'.

        gene (str): Full name of a germline gene (using IMGT-style names, like IGHV1-2*02).
                    If provided, a single ``IMGTGermlineGene`` object will be returned, or None if the
                    specified gene could not be found. If not provided, a list of ``IMGTGermlineGene``
                    objects for all germline genes matching the ``species`` and ``gene_type`` will be returned.

    Returns:
    --------

        IMGTGermlineGene: a single ``IMGTGermlineGene`` object (if ``gene`` is provided) or a list of
                          ``IMGTGermlineGene`` objects.
    '''
    mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    db_file = os.path.join(mod_dir, 'vdj/germline_dbs/imgt_gapped/{}_{}_imgt-gapped.fasta'.format(species, gene_type))
    try:
        germs = [IMGTGermlineGene(g) for g in SeqIO.parse(open(db_file, 'r'), 'fasta')]
    except:
        # TODO: log that the germline database file couldn't be found
        return None
    if gene is None:
        return germs
    try:
        return [g for g in germs if g.name == gene][0]
    except IndexError:
        return None




class IMGTGermlineGene(object):
    """docstring for IMGTGermlineGene"""

    species_lookup = {'homo sapiens': 'human'}

    def __init__(self, sequence, species=None):
        self.raw_sequence = Sequence(str(sequence.seq), id=sequence.description)
        self._species = species
        self.gapped_nt_sequence = self.raw_sequence.sequence
        self.ungapped_nt_sequence = self.gapped_nt_sequence.replace('.', '')


    @lazy_property
    def accession(self):
        return self.raw_sequence.id.split('|')[0].strip()

    @lazy_property
    def name(self):
        return self.raw_sequence.id.split('|')[1].strip()

    @property
    def species(self):
        if self._species is None:
            self._species = self.species_lookup[self.raw_sequence.id.split('|')[2].strip().lower()]
        return self._species

    @lazy_property
    def functionality(self):
        return self.raw_sequence.id.split('|')[3].strip()

    @lazy_property
    def gene_type(self):
        return self.raw_sequence.id.split('|')[4].strip()[0].upper()

    @lazy_property
    def coding_start(self):
        # NOTE: uses 1-based indexing, so need to adjust if using for slicing
        return int(self.raw_sequence.id.split('|')[7].strip())

    @lazy_property
    def nt_length(self):
        return int(self.raw_sequence.id.split('|')[12].strip().split('+')[0])

    @lazy_property
    def gap_length(self):
        return int(self.raw_sequence.id.split('|')[12].strip().split('+')[1].split('=')[0])

    @lazy_property
    def total_length(self):
        return int(self.raw_sequence.id.split('|')[12].strip().split('=')[1])

    @lazy_property
    def partial(self):
        p = []
        partial = self.raw_sequence.id.split('|')[13]
        if "3'" in partial:
            p.append("3'")
        if "5'" in partial:
            p.append("5'")
        return p

    @lazy_property
    def is_rev_comp(self):
        if self.raw_sequence.id.split('|')[14].strip() == 'rev-compl':
            is_rev_comp = True
        else:
            is_rev_comp = False
        return is_rev_comp

    @lazy_property
    def gapped_aa_sequence(self):
        res = []
        coding = self.gapped_nt_sequence[self.coding_start - 1:]
        for codon in (coding[pos:pos + 3] for pos in xrange(0, len(coding), 3)):
            if len(codon) != 3:
                continue
            if '.' in codon:
                res.append('.')
            else:
                res.append(codons.get(codon, 'X'))
        return ''.join(res)

    @lazy_property
    def ungapped_aa_sequence(self):
        return self.gapped_aa_sequence.replace('.', '')
