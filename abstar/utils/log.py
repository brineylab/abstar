#!/usr/bin/env python
# filename: log.py


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


class LoggingMixin():
    """docstring for LoggingMixin"""
    def __init__(self):
        self._log = None
        self._exceptions = None


    def log(self, *args, **kwargs):
        sep = kwargs.get('sep', ' ')
        lstring = sep.join([str(a) for a in args])
        if self._log is None:
            self.log = [lstring, ]
        else:
            self._log.append(lstring)


    def exception(self, *args, **kwargs):
        sep = kwargs.get('sep', '\n')
        estring = sep.join([str(a) for a in args])
        if self._exceptions is None:
            self._exceptions = [estring, ]
        else:
            self._exceptions.append(estring)


    def format_log(self):
        '''
        Formats the antibody log.

        Log formatting is only performed on sequences that had an
        error during annotation, unless AbStar is run in debug
        mode. In debug mode, all sequences will be logged.

        Returns:
        --------

            str: Formatted log string.
        '''
        # self._log += ['', '']
        output = '\n'.join(self._log)
        if self._check_for_exceptions():
            output += self._format_exceptions()
        return output


    def initialize_vdj_log(self):
        log = ['', ]
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


    def initialize_antibody_log(self):
        log = []
        log.append('=' * len(self.id))
        log.append(self.id)
        log.append('=' * len(self.id))
        log.append('')
        log.append(self.raw_input.fasta)
        log.append('')
        log.append('RAW INPUT: {}'.format(self.raw_input.sequence))
        log.append('ORIENTED INPUT: {}'.format(self.oriented_input.sequence))
        log.append('CHAIN: {}'.format(self.chain))
        log.append('')
        if self.v is not None:
            log.append('V-GENE: {}'.format(self.v.full))
        if self.d is not None:
            log.append('D-GENE: {}'.format(self.d.full))
        if self.j is not None:
            log.append('J-GENE: {}'.format(self.j.full))
        self._log = log


    def _check_for_exceptions(self):
        if self._exceptions:
            return True
        if self.v is not None:
            if self.v._exceptions:
                return True
        if self.d is not None:
            if self.d._exceptions:
                return True
        if self.j is not None:
            if self.j._exceptions:
                return True
        return False


    def _format_exceptions(self):
        estring = 'EXCEPTIONS\n'
        estring += '----------\n\n'
        if self.v is not None:
            if self.v._exceptions is not None:
                self._exceptions += self.v._exceptions
        if self.d is not None:
            if self.d._exceptions is not None:
                self._exceptions += self.d._exceptions
        if self.j is not None:
            if self.j._exceptions is not None:
                self._exceptions += self.j._exceptions
        estring += '\n\n'.join([e for e in self._exceptions])
        return estring


# def setup_logging(logfile, debug=False):
#     fmt = '[%(levelname)s] %(name)s %(asctime)s %(message)s'
#     if debug:
#         logging.basicConfig(filename=logfile,
#                             filemode='w',
#                             format=fmt,
#                             level=logging.DEBUG)
#     else:
#         logging.basicConfig(filename=logfile,
#                             filemode='w',
#                             format=fmt,
#                             level=logging.INFO)
#     logger = logging.getLogger('log')
#     logger = add_stream_handler(logger)
#     logger.info('LOG LOCATION: {}'.format(logfile))


# def get_logger(name=None):
#     logger = logging.getLogger(name)
#     if len(logger.handlers) == 0:
#         logger = add_stream_handler(logger)
#     return logger


# def add_stream_handler(logger):
#     formatter = logging.Formatter("%(message)s")
#     ch = logging.StreamHandler()
#     ch.setFormatter(formatter)
#     ch.setLevel(logging.INFO)
#     logger.addHandler(ch)
#     return logger
