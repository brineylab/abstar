#!/usr/bin/env python
# filename: germline_dbs.py

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

from argparse import ArgumentParser
import os
import platform
import shutil
import subprocess as sp
import sys

from Bio import SeqIO

if sys.version_info[0] > 2:
    raw_input = input



def parse_arguments():
    parser = ArgumentParser("Creates AbStar germline databases using an IMGT-gapped FASTA files of germline sequences. \
        Properly formatted germline sequence files can be obtained from: http://www.imgt.org/genedb/")
    parser.add_argument('-v', '--variable', dest='v', required=True,
                        help="Path to an IMGT-gapped, FASTA-formatted file containing Variable gene sequences. \
                        Sequences for both heavy and light chains should be included in a single file. \
                        Required.")
    parser.add_argument('-d', '--diversity', dest='d', required=True,
                        help="Path to an IMGT-gapped, FASTA-formatted file containing Diversity gene sequences. \
                        Required.")
    parser.add_argument('-j', '--joining', dest='j', required=True,
                        help="Path to an IMGT-gapped, FASTA-formatted file containing Joining gene sequences. \
                        Sequences for both heavy and light chains should be included in a single file. \
                        Required.")
    parser.add_argument('-i', '--isotypes', dest='isotypes', default=None,
                        help="Path to a FASTA-formatted file containing isotype sequences. \
                        The name of the isotype in the FASTA file is what will be reported by AbStar, \
                        so the use of standard nomenclature (IgG, IgG1, etc) is encouraged. \
                        If an isotype file is not provided, isotypes will not be parsed by AbStar when \
                        using the germline database.")
    parser.add_argument('-m', '--manifest', dest='manifest', default=None,
                        help="Path to a plain-text file containing information about the germline database. \
                        Format is not important, but this is the place for optional inforamtion like the origin \
                        of the germline database, the date of download, etc. Optional. If provided, the contents \
                        of the file will be stored as 'manifest.txt' in the top level of the germline database \
                        (for user-provided databases, that would be '~/.abstar/germline_dbs/<species>/).")
    parser.add_argument('-n', '--name', dest='name', required=True,
                        help="Name of the resulting germline database. \
                        If an abstar germline database with the same already exists, it will be replaced. \
                        Germline database names are converted to lowercase, so 'Human' and 'human' are equivalent. \
                        User-added germline databases will persist even after abstar updates, so if you have added a \
                        'human' database and a new version of AbStar contains an updated 'human' database, the user-added \
                        database will still be used after the update.")
    parser.add_argument('-l', '--location', dest='db_location', default=None,
                        help="Location into which the new germline databases will be deposited. \
                        Default is '~/.abstar/'. \
                        Note that AbStar will only use user-generated databases found in '~/abstar/', \
                        so this option is provided primarily to test database creation without overwriting \
                        current databases of the same name. \
                        If the directory does not exist, it will be created.")
    parser.add_argument('--include-species-in-gene-name', default=True,
                        help="If `True`, the species is parsed from the IMGT header and appended to the \
                        germline gene name. If `False`, the species is not parsed, and the species name reported \
                        by abstar will be the database name ('human', for example). Default is `True`.")
    parser.add_argument('--allow-partials', dest='allow_partials', default=None, choices=("3'", "5'", "both"),
                        help="If set, will include 'partial' sequences (partial in 3', partial in 5' or both) to the \
                        AbStar germline database. 'Partial' is an annotation (field 14 in the IMGT/GENE-DB header) \
                        that indicates a germline sequence is incomplete, even if it is otherwise considered \
                        functional. Default is to exclude all partial sequences.")
    parser.add_argument('-D', '--debug', dest="debug", action='store_true', default=False,
                        help="More verbose logging if set.")
    args = parser.parse_args()
    return args


# -------------------------
#
#    FILES/DIRECTORIES
#
# -------------------------


def get_addon_directory(db_location):
    if db_location is not None:
        print('\n')
        print('NOTE: You have selected a non-default location for the germline directory.')
        string = 'AbStar only looks in the default location (~/.abstar/) for user-created germline databases, '
        string += 'so this database will not be used by AbStar. The custom database location option is primarily '
        string += 'provided so that users can test the database creation process without overwriting existing databases.'
        print(string)
        addon_dir = db_location
    else:
        addon_dir = os.path.expanduser('~/.abstar/germline_dbs')
    if not os.path.isdir(addon_dir):
        os.makedirs(addon_dir)
    return addon_dir


def get_binary_directory():
    mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    bin_dir = os.path.join(mod_dir, 'assigners/bin')
    return bin_dir


def check_for_existing_db(addon_dir, dbname):
    dbs = [os.path.basename(d[0]) for d in os.walk(addon_dir)]
    if dbname.lower() in dbs:
        print('\n')
        print('WARNING: A {} germline database already exists.'.format(dbname.lower()))
        print('Creating a new database with that name will overwrite the old one.')
        keep_going = raw_input('Do you want to continue? [y/N]: ')
        if keep_going.upper() not in ['Y', 'YES']:
            print('')
            print('Aborting germline database creation.')
            print('\n')
            sys.exit()


def make_db_directories(addon_dir, dbname, isotypes):
    dbname_dir = os.path.join(addon_dir, dbname.lower())
    if not os.path.isdir(dbname_dir):
        os.makedirs(dbname_dir)
    db_names = ['imgt_gapped', 'ungapped', 'blast']
    if isotypes is not None:
        db_names.append('isotypes')
    for db_name in db_names:
        db_dir = os.path.join(dbname_dir, db_name)
        if os.path.isdir(db_dir):
            shutil.rmtree(db_dir)
        os.makedirs(db_dir)


# -------------------------
#
#    DATABASE CREATION
#
# -------------------------


def make_blast_db(ungapped_germline_file, addon_directory, segment, dbname):
    print('  - BLASTn')
    bin_dir = get_binary_directory()
    mbd_binary = os.path.join(bin_dir, 'makeblastdb_{}'.format(platform.system().lower()))
    mbd_output = os.path.join(addon_directory, '{}/blast/{}'.format(dbname.lower(), segment.lower()))
    mbd_log = os.path.join(addon_directory, '{}/blast/{}.blastlog'.format(dbname.lower(), segment.lower()))
    mbd_cmd = '{} -in {} -out {} -parse_seqids -dbtype nucl -title {} -logfile {}'.format(mbd_binary,
                                                                                          ungapped_germline_file,
                                                                                          mbd_output,
                                                                                          segment.lower(),
                                                                                          mbd_log)
    p = sp.Popen(mbd_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    return mbd_output, stdout, stderr


def make_ungapped_db(ungapped_germline_file, addon_directory, segment, dbname, include_species_in_gene_name):
    print('  - ungapped FASTA')
    output_file = os.path.join(addon_directory, '{}/ungapped/{}.fasta'.format(dbname.lower(), segment.lower()))
    seqs = SeqIO.parse(open(ungapped_germline_file), 'fasta')
    if include_species_in_gene_name:
        fastas = []
        for s in seqs:
            gene = s.description.split('|')[1]
            species = s.description.split('|')[2].lower().replace(' ', '-')
            fastas.append('>{}__{}\n{}'.format(gene, species, str(s.seq).replace('.', '')))
    else:
        fastas = ['>{}\n{}'.format(s.description.split('|')[1], str(s.seq).replace('.', '')) for s in seqs]
    open(output_file, 'w').write('\n'.join(fastas))
    return output_file


def make_imgt_gapped_db(input_file, addon_directory, segment, dbname, allow_partials):
    print('  - IMGT-gapped FASTA')
    output_file = os.path.join(addon_directory, '{}/imgt_gapped/{}.fasta'.format(dbname.lower(), segment.lower()))
    seqs = sorted(list(SeqIO.parse(open(input_file), 'fasta')), key=lambda x: x.id)
    fastas = []
    for s in seqs:
        partial = s.description.split('|')[13].strip()
        if "3" in partial:
            if any([allow_partials == "3'", allow_partials == 'both']):
                fastas.append('>{}\n{}'.format(s.description, str(s.seq).upper()))
        elif "5" in partial:
            if any([allow_partials == "5'", allow_partials == 'both']):
                fastas.append('>{}\n{}'.format(s.description, str(s.seq).upper()))
        else:
            fastas.append('>{}\n{}'.format(s.description, str(s.seq).upper()))
    open(output_file, 'w').write('\n'.join(fastas))
    return output_file


def make_isotype_db(input_file, addon_directory, dbname):
    print('  - isotypes')
    output_file = os.path.join(addon_directory, '{}/isotypes/isotypes.fasta'.format(dbname.lower()))
    seqs = sorted(list(SeqIO.parse(open(input_file), 'fasta')), key=lambda x: x.id)
    fastas = ['>{}\n{}'.format(s.id, str(s.seq).upper()) for s in seqs]
    open(output_file, 'w').write('\n'.join(fastas))
    return output_file


def transfer_manifest_data(manifest, addon_directory, dbname):
    output_file = os.path.join(addon_directory, '{}/manifest.txt'.format(dbname.lower()))
    manifest_data = open(manifest, 'r').read()
    open(output_file, 'w').write(manifest_data)
    return output_file


# -------------------------
#
#        PRINTING
#
# -------------------------


def print_segment_info(segment, input_file):
    seqs = list(SeqIO.parse(open(input_file), 'fasta'))
    print('\n')
    seg_string = '  ' + segment.upper() + '  '
    print('-' * len(seg_string))
    print(seg_string)
    print('-' * len(seg_string))
    print(input_file)
    print('input file contains {} sequences'.format(len(seqs)))
    print('')
    print('Building germline databases:')


def print_manifest_info(manifest):
    seg_string = '  MANIFEST  '
    print('\n')
    print('-' * len(seg_string))
    print(seg_string)
    print('-' * len(seg_string))
    print(manifest)
    print('')
    print('Transferring manifest data...')


# -------------------------
#
#          MAIN
#
# -------------------------


def main():
    args = parse_arguments()
    addon_dir = get_addon_directory(args.db_location)
    check_for_existing_db(addon_dir, args.name)
    make_db_directories(addon_dir, args.name, args.isotypes)
    for segment, infile in [('Variable', args.v), ('Diversity', args.d), ('Joining', args.j)]:
        print_segment_info(segment, infile)
        imgt_gapped_file = make_imgt_gapped_db(infile,
                                               addon_dir,
                                               segment[0].lower(),
                                               args.name,
                                               args.allow_partials)
        ungapped_file = make_ungapped_db(imgt_gapped_file,
                                         addon_dir,
                                         segment[0].lower(),
                                         args.name,
                                         args.include_species_in_gene_name)
        blast_file, stdout, stderr = make_blast_db(ungapped_file, addon_dir, segment[0].lower(), args.name)
        if args.debug:
            print(stdout)
            print(stderr)
    if args.isotypes is not None:
        print_segment_info('ISOTYPES', args.isotypes)
        isotype_file = make_isotype_db(args.isotypes, addon_dir, args.name)
    if args.manifest is not None:
        print_manifest_info(args.manifest)
        transfer_manifest_data(args.manifest, addon_dir, args.name)
    print('\n')


if __name__ == '__main__':
    main()
