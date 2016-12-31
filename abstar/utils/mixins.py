#!/usr/bin/env python
# filename: mixins.py


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


class LoggingMixin():
    """docstring for LoggingMixin"""
    def __init__(self):
        self._log = None
        self._exceptions = None


    def log(self, *args, **kwargs):
        '''
        Records a log message
        '''
        sep = kwargs.get('sep', ' ')
        lstring = sep.join([str(a) for a in args])
        if self._log is None:
            self.log = [lstring, ]
        else:
            self._log.append(lstring)


    def exception(self, *args, **kwargs):
        '''
        Records an exception.
        '''
        sep = kwargs.get('sep', '\n')
        estring = sep.join([str(a) for a in args])
        if self._exceptions is None:
            self._exceptions = [estring, ]
        else:
            self._exceptions.append(estring)


    def format_log(self):
        '''
        Formats the log, including exceptions.

        Log formatting will only be performed on sequences that had an
        error during annotation, unless AbStar is run in debug
        mode. In debug mode, all sequences will be logged.

        Returns:
        --------

            str: Formatted log string.
        '''
        output = ''
        if self._log is not None:
            output += '\n'.join(self._log)
        if self._check_for_exceptions():
            output += '\n\n'
            output += self._format_exceptions()
        output += '\n\n'
        return output


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
        if self._exceptions is None:
            self._exceptions = []
        estring = '\n\nEXCEPTIONS\n'
        estring += '----------\n\n'
        if self.v is not None:
            self._exceptions += self.v._exceptions
        if self.d is not None:
            self._exceptions += self.d._exceptions
        if self.j is not None:
            self._exceptions += self.j._exceptions
        estring += '\n\n'.join([e for e in self._exceptions])
        return estring
