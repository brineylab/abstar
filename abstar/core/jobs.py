#!/usr/bin/python
# filename: jobs.py

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


from __future__ import absolute_import, division, print_function, unicode_literals

from multiprocessing import Pool

from ..assigners.registry import ASSIGNERS
from ..utils.output import get_abstar_results, write_output
from ..utils.queue.celery import celery

from abutils.utils import log


@celery.task
def run_abstar(sequence_file, output_directory, args):
    """
    Wrapper function to multiprocess (or not) the assignment of V, D and J
    germline genes. Also writes the JSON-formatted output to file.

    Input is a a FASTA-formatted file of antibody sequences and the output directory.
    Optional input items include the species (supported species: 'human'); length of
    the unique antibody identifier (UAID); and debug mode (which forces single-threading
    and prints more verbose errors.)

    Output is the number of functional antibody sequences identified in the input file.
    """
    try:
        # setup logging
        global logger
        logger = log.get_logger(__name__)
        assigned_log = ""
        unassigned_log = ""
        # identify output file
        output_filename = os.path.basename(seq_file)
        if args.output_type == "json":
            output_file = os.path.join(output_dir, output_filename + ".json")
        elif args.output_type in ["imgt", "hadoop"]:
            output_file = os.path.join(output_dir, output_filename + ".txt")
        # start assignment
        assigner = ASSIGNERS[args.assigner]
        assigner(sequence_file, args.species)
        # process all of the successfully assigned sequences
        assigned = [Antibody(vdj, args.species) for vdj in assigner.assigned]
        for ab in assigned:
            ab.annotate()
            if args.debug:
                assigned_log += ab.format_log()
        results = get_abstar_results(
            assigned, pretty=args.pretty, padding=args.padding, raw=args.raw
        )
        write_output(results, output_file, args.output_type, args.parquet)
        # capture the log for all unsuccessful sequences
        for vdj in unassigned:
            unassigned_log += vdj.format_log()

        return (len(assigned), assigned_log, unassigned_log)

    #     vdj_output = process_sequence_file(seq_file, args)
    #     if not vdj_output:
    #         return None
    #     clean_vdjs = [vdj for vdj in vdj_output if vdj.rearrangement]
    #     output_count = write_output(clean_vdjs, output_file, args.output_type, args.pretty, args.padding)
    #     return (output_file, output_count)
    except:
        logger.debug(traceback.format_exc())
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def run_jobs(files, output_dir, args):
    sys.stdout.write("\nRunning VDJ...\n")
    if args.cluster:
        return _run_jobs_via_celery(files, output_dir, args)
    elif args.debug or args.chunksize == 0:
        return _run_jobs_singlethreaded(files, output_dir, args)
    else:
        return _run_jobs_via_multiprocessing(files, output_dir, args)


def _run_jobs_singlethreaded(files, output_dir, args):
    # from utils.vdj import run as run_abstar
    results = []
    for i, f in enumerate(files):
        try:
            results.append(run_abstar(f, output_dir, args))
            update_progress(i + 1, len(files))
        except:
            logger.debug("FILE-LEVEL EXCEPTION: {}".format(f))
            logging.debug(traceback.format_exc())
    logger.info("")
    return results


def _run_jobs_via_multiprocessing(files, output_dir, args):
    # from utils.vdj import run as run_abstar
    p = Pool(maxtasksperchild=50)
    async_results = []
    for f in files:
        async_results.append((f, p.apply_async(run_abstar, (f, output_dir, args))))
    monitor_mp_jobs([ar[1] for ar in async_results])
    results = []
    for a in async_results:
        try:
            results.append(a[1].get())
        except:
            logger.debug("FILE-LEVEL EXCEPTION: {}".format(a[0]))
            if args.debug:
                traceback.print_exc()
            logging.debug("".join(traceback.format_exc()))
            continue
    p.close()
    p.join()
    return results


def monitor_mp_jobs(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        ready = [ar for ar in results if ar.ready()]
        finished = len(ready)
        update_progress(finished, jobs)
    sys.stdout.write("\n\n")


def _run_jobs_via_celery(files, output_dir, args):
    # from utils.vdj import run as run_abstar
    async_results = []
    for f in files:
        async_results.append(run_abstar.delay(f, output_dir, args))
    succeeded, failed = monitor_celery_jobs(async_results)
    failed_files = [f for i, f in enumerate(files) if async_results[i].failed()]
    for ff in failed_files:
        logger.debug("FAILED FILE: {}".format(f))
    return [s.get() for s in succeeded]


def monitor_celery_jobs(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        succeeded = [ar for ar in results if ar.successful()]
        failed = [ar for ar in results if ar.failed()]
        finished = len(succeeded) + len(failed)
        update_progress(finished, jobs, failed=len(failed))
    sys.stdout.write("\n\n")
    return succeeded, failed


def update_progress(finished, jobs, failed=None):
    pct = int(100.0 * finished / jobs)
    ticks = pct / 2
    spaces = 50 - ticks
    if failed:
        prog_bar = "\r({}/{}) |{}{}|  {}% ({}, {})".format(
            finished, jobs, "|" * ticks, " " * spaces, pct, finished - failed, failed
        )
    else:
        prog_bar = "\r({}/{}) |{}{}|  {}%".format(
            finished, jobs, "|" * ticks, " " * spaces, pct
        )
    sys.stdout.write(prog_bar)
    sys.stdout.flush()
