#!/usr/bin/python
# filename: mongoimport.py

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

import itertools
import multiprocessing as mp
import os
import subprocess
import sys
import time
import uuid

from abtools.utils import progbar
from abtools.utils.pipeline import make_dir


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser("Performs the mongoimport operation on all files in the input directory. \
                                      Determines the appropriate collection by parsing the filename.")
    parser.add_argument('-i', '--ip', dest='ip', default='localhost',
                        help="The IP address of the MongoDB server. Defaults to 'localhost'.")
    parser.add_argument('-P', '--port', dest='port', type=int, default=27017,
                        help="The MongoDB port. Defaults to 27017.")
    parser.add_argument('-u', '--user', dest='user', default=None,
                        help="Username for the MongoDB server. Not used if not provided.")
    parser.add_argument('-p', '--password', dest='password', default=None,
                        help="Password for the MongoDB server. Not used if not provided.")
    parser.add_argument('-f', '--input', dest='input', required=True,
                        help="A directory containing multiple JSON files for import to MongoDB, \
                        or a list of JSON file paths.")
    parser.add_argument('-d', '--db', dest='db', required=True,
                        help="The MongoDB database for import.")
    parser.add_argument('-l', '--log', dest='log', default='sys.stdout',
                        help="Log file for the mongoimport stdout.")
    parser.add_argument('-t', '--temp', dest='temp', default=None,
                        help="Temp directory if splitting JSON files prior to mongoimport. \
                        if not provided, a temp directory will be created in the --in directory.")
    parser.add_argument('-e', '--delim1', dest='delim1', default=None,
                        help="The first character delimiter used to split the filename to get the collection name. \
                        If splitting with a single delimiter, use this option to provide the delimiter. Required.")
    parser.add_argument('--delim2', dest='delim2', default=None,
                        help="If splitting with two different delimiters, use this option to provide the second delimiter. \
                        Required if splitting with two delimiters.")
    parser.add_argument('-s', '--split1', dest='split1_pos', default=1, type=int,
                        help="Builds the collection name by truncating at the <split> occurance of the <delim> character. \
                        If splitting with multiple delimiters, this option is used to specifiy the occurance of the first delimiter. \
                        Default is 1.")
    parser.add_argument('--split2', dest='split2_pos', default=1, type=int,
                        help="If splitting with multiple delimiters, this option is used to specify the occurance of the \
                        second delimiter at which to split. \
                        Required if splitting with two different delimiters.")
    parser.add_argument('--split-file', dest='split_file', action='store_true', default=False,
                        help="Splits each input file into subfiles containing --split_file_lines lines. \
                        Useful when performing mongoimport via an SSH tunnel, where for some reason MongoDB \
                        errors when importing files greater than 16MB (even if no individual documents are >16MB).")
    parser.add_argument('--split=file-lines', dest='split_file_lines', type=int, default=500,
                        help="Number of lines in each mongoimport subfile (each line should be a complete JSON document).")
    parser.add_argument('-x', '--split-only', dest='split_only', default=False, action='store_true',
                        help="Instead of truncating the filename to get the collection name, takes only the split for the collection. \
                        Default is False.")
    return parser.parse_args()


class Args(object):
    def __init__(self, ip='localhost', port=27017, user=None, password=None, db=None,
                 input=None, log=None, split_file=False, split_file_lines=500,
                 delim1=None, delim2=None,
                 split1_pos=1, split2_pos=1, split_only=False):
        super(Args, self).__init__()
        self.ip = ip
        self.port = int(port)
        self.user = user
        self.password = password
        self.input = input
        self.db = db
        self.log = log if log else sys.stdout
        self.split_file = split_file
        self.split_file_lines = int(split_file_lines)
        self.delim1 = delim1
        self.delim2 = delim2
        self.split1_pos = int(split1_pos)
        self.split2_pos = int(split2_pos)
        self.split_only = split_only


def mongo_import(json, db, coll, log, args):
    if args.split_file:
        jsons = split_file(json, args)
        multiprocess_mongoimport(jsons, db, coll, args)
    else:
        do_mongoimport(json, args.ip, args.port, db, coll, args.user, args.password)

    #     jsons = [json, ]
    # username = " -u {}".format(args.user) if args.user else ""
    # password = " -p {}".format(args.password) if args.password else ""
    # # user_password = "{}{} --authenticationDatabase admin".format(username, password)
    # for i, json in enumerate(jsons):
    #     mongo_cmd = "mongoimport --host {}:27017{}{} --db {} --collection {} --file {}".format(
    #         args.ip, username, password, db, coll, json)
    #     mongo = subprocess.Popen(mongo_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #     mongo.communicate()
    #     if args.split_file:
    #         progbar.progress_bar(i + 1, len(jsons))
    # if args.split_file:
    #     print('')
    #     remove_temp_files(args)


def multiprocess_mongoimport(jsons, db, coll, args):
    progbar.progress_bar(0, len(jsons))
    async_results = []
    p = mp.Pool()
    for j in jsons:
        async_results.append(p.apply_async(do_mongoimport, args=(j, args.ip, args.port, db, coll, args.user, args.password)))
    monitor_results(async_results)
    remove_temp_files(args)
    print('')


def monitor_results(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        ready = [ar for ar in results if ar.ready()]
        finished = len(ready)
        progbar.progress_bar(finished, jobs)
    progbar.progress_bar(finished, jobs)


def do_mongoimport(json, ip, port, db, coll, user, password):
    username = " -u {}".format(user) if user else ""
    password = " -p {}".format(password) if password else ""
    mongo_cmd = "mongoimport --host {}:{}{}{} --db {} --collection {} --file {} --numInsertionWorkers 12 --batchSize 100".format(
        ip, port, username, password, db, coll, json)
    mongo = subprocess.Popen(mongo_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    mongo.communicate()



def listdir_fullpath(d):
    if os.path.isdir(d):
        files = [os.path.join(d, f) for f in os.listdir(d) if f.split('.')[-1].lower() == 'json']
        if files:
            return sorted(files)
    else:
        if d.split('.')[-1].lower() == 'json':
            return [d, ]
    print("ERROR: no importable JSON files were found. Double-check your input file/directory")
    sys.exit(1)


def split_file(json, args):
    split_files = []
    temp_dir = args.temp if args.temp is not None else os.path.join(args.mongo_input_dir, 'temp')
    make_dir(temp_dir)
    with open(json) as f:
        for chunk in itertools.izip_longest(*[f] * args.split_file_lines):
            chunk = [c for c in chunk if c is not None]
            fname = os.path.join(temp_dir, str(uuid.uuid4()) + '.json')
            open(fname, 'w').write(''.join(chunk))
            split_files.append(fname)
    return split_files


def remove_temp_files(args):
    temp_dir = args.temp if args.temp is not None else os.path.join(args.mongo_input_dir, 'temp')
    files = listdir_fullpath(temp_dir)
    for f in files:
        os.unlink(f)



def get_collection(i, args):
    # split with two different delimiters
    if args.delim1 and args.delim2:
        bname = os.path.basename(i)
        split1 = args.delim1.join(bname.split(args.delim1)[args.split1_pos:])
        removed_delim2s = bname.count(args.delim2) - split1.count(args.delim2)
        return args.delim2.join(split1.split(args.delim2)[:args.split2_pos - removed_delim2s])
    # split with a single delimiter
    elif args.delim1:
        delim = str(args.delim1)
        split = args.split1_pos
        if args.split_only:
            return os.path.basename(i).split(delim)[split - 1]
        if split <= 1:
            return os.path.basename(i).split(delim)[0]
        else:
            pre_collection = os.path.basename(i).split(delim)
            return delim.join(pre_collection[:split])
    else:
        return os.path.basename(i)


def preflight_checks(args):
    if not args.db:
        err = "ERROR: a MongoDB database name must be provided."
        raise RuntimeError(err)
    if not args.input:
        err = "ERROR: a directory of JSON files must be provided."
        raise RuntimeError(err)
    if args.delim2 and not args.delim1:
        err = "ERROR: --delim1 must be provided with --delim2.\n"
        err += "If you only want to use one delimiter, use --delim1."
        raise RuntimeError(err)


def run(**kwargs):
    '''
    Imports one or more JSON files into a MongoDB database.

    Examples:

        To import a single JSON file into MyDatabase on a local MongoDB database::

            from abstar.utils import mongoimport

            mongoimport.run(input='/path/to/MySequences.json', db='MyDatabase')

        This will result in a collection named 'MySequences.json' being created in MyDatabase
        on your local MongoDB instance (if it doesn't already exist) and the data from
        MySequences.json being imported into that collection.

        Doing the same thing, but with a remote MongoDB server running on port 27017::

            mongoimport.run(ip='123.45.67.89',
                            user='my_username',
                            password='Secr3t',
                            input='/path/to/MySequences.json',
                            db='MyDatabase')

        But what if we want the collection name to be different than the file name? We can
        truncate the filename at the first occurance of any given pattern with ``delim1``::

            mongoimport.run(input='/path/to/MySequences.json,
                            db='MyDatabase',
                            delim1='.')

        In this case, the collection name is created by truncating the input file name
        at the first occurance of ``.``, so the collection name would be MySequences.
        We can also truncate the filename at the Nth occurance of any given pattern
        by using ``delim1`` with ``split1_pos``::

            mongoimport.run(input='/path/to/my_sequences_2016-01-01.json,
                            db='MyDatabase',
                            delim1='_',
                            split1_pos=2)

        which results in a collection name of my_sequences.

        If we have more complex filenames, we can use ``delim1`` in combination with
        ``delim2``. When ``delim1`` and ``delim2`` are used together, ``delim1`` becomes
        the pattern used to cut the filename on the left and ``delim2`` is used to cut
        the filename on the right. For example, if our filename is
        ``plate-2_SampleName-01_redo.json`` and we want the collection to be named
        ``SampleName``, we would set ``delim1`` to ``_`` and ``delim2`` to ``-``. We also
        need to specify that we want to cut at the second occurance of ``delim2``,
        which we can do with ``split2_pos``::

            mongoimport.run(input='/path/to/plate-2_SampleName-01_redo.json,
                            db='MyDatabase',
                            delim1='_',
                            delim2='-',
                            split2_pos=2)

        Trimming filenames this way is nice, but it becomes much more useful if you're
        importing more than one file at a time. ``mongoimport.run()`` will accept a list
        of file names, and will generate separate collection names for each input file::

            files = ['/path/to/A01-Sample01_2016-01-01',
                     '/path/to/A02-Sample02_2016-01-01',
                     '/path/to/A03-Sample03_2016-01-01]

            mongoimport.run(input=files,
                            db='MyDatabase',
                            delim1='-',
                            delim2='_')

        The three input files will be imported into collections ``Sample01``, ``Sample02``
        and ``Sample03``, respectively. Finally, you can pass the path to a directory
        containing oen or more JSON files, and all the JSON files will be imported::

            mongoimport.run(input='/path/to/output/directory',
                            db='MyDatabase',
                            delim1='-',
                            delim2='_')

    Args:

        input (str, list): Input is required and may be one of three things:

            1) A list/tuple of JSON file paths
            2) A path to a single JSON file
            3) A path to a directory containing one or more JSON files.

        ip (str): The IP address of the MongoDB server. Default is 'localhost'.

        port (int): MongoDB port. Default is 27017.

        user (str): Username with which to connect to the MongoDB database. If either
            of ``user`` or ``password`` is not provided, ``mongoimport.run()`` will attempt
            to connect to the MongoDB database without authentication.

        password (str): Password with which to connect to the MongoDB database. If either
            of ``user`` or ``password`` is not provided, ``mongoimport.run()`` will attempt
            to connect to the MongoDB database without authentication.

        db (str): Name of the MongoDB database for import. Required.

        log (str): Path to a logfile. If not provided log information will be written to stdout.

        delim1 (str): Pattern on which to split the input file to generate the collection name.
            Default is ``None``, which results in the file name being used as the collection name.

        split1_pos (int): Occurance of ``delim1`` on which to split the input file name.
            Default is 1.

        delim2 (str): Second pattern on which to split the input file name to generate the
            collection name. Default is None, which results in only ``delim1`` being used.

        split2_pos (int): Occurance of ``delim2`` on which to split the input file name.
            Default is 1.
    '''
    args = Args(**kwargs)
    main(args)


def main(args):
    preflight_checks(args)
    if type(args.input) in [list, tuple]:
        in_files = args.input
    elif os.path.isfile(args.input):
        in_files = [args.input, ]
    elif os.path.isdir(args.input):
        in_files = listdir_fullpath(args.input)
    else:
        err = 'ERROR: Input not recognized. Valid input can be one of three things:\n'
        err += '  1. a list of JSON file paths\n'
        err += '  2. the path to a single JSON file\n'
        err += '  3. the path to a directory of JSON files\n\n'
        err += 'You supplied: {}\n'.format(input)
        raise RuntimeError(err)
    log_handle = open(args.log, 'a')
    open(args.log, 'w').write('')
    log_handle.write('\nImporting {} files.\n'.format(len(in_files)))
    for i, f in enumerate(in_files):
        coll = get_collection(f, args)
        print("\n[ {} ] Importing {} into collection {}.".format(i + 1, os.path.basename(f), coll))
        log_handle.write('\n\n----------------------------------------')
        log_handle.write('File: {0}\nCollection: {1}'.format(f, coll))
        log_handle.write('----------------------------------------\n')
        mongo_import(f, args.db, coll, log_handle, args)
    print("\nDone. {0} files were imported into MongoDB.\n\n".format(len(in_files)))


if __name__ == '__main__':
    args = parse_arguments()
    main(args)
