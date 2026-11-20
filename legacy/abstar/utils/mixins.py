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


class LoggingMixin:
    """
    Mixin for logging messages and exceptions.

    Methods:
    --------
        log(self, *args, separator: str = " ") -> None:
            Records a log message, with each argument joined by the separator.

        exception(self, *args, separator: str = "\n") -> None:
            Records an exception, with each argument joined by the separator.

        format_log(self, separator: str = "\n") -> str:
            Formats the log, including exceptions, as a string to be written to file or stdout.

    """

    def __init__(self):
        super().__init__()
        self.logs = []
        self.exceptions = []

    def log(self, *args, separator: str = " ") -> None:
        """
        Records a log message

        Parameters:
        ----------
            *args : str
                The log message to record.

            separator : str, default: " "
                The separator to join the arguments.

        """
        log_str = separator.join([str(a) for a in args])
        self.logs.append(log_str)

    def exception(self, *args, separator: str = "\n") -> None:
        """
        Records an exception

        Parameters:
        ----------
            *args : str
                The exception to record.

            separator : str, default: "\n"
                The separator to join the arguments.

        """
        exception_str = separator.join([str(a) for a in args])
        self.exceptions.append(exception_str)

    def format_log(self, separator="\n") -> str:
        """
        Formats the log, including exceptions.

        Log formatting will only be performed on sequences that had an
        error during annotation, unless AbStar is run in debug
        mode. In debug mode, all sequences will be logged.

        Parameters:
        ----------
            separator : str, default: "\n"
                The separator to join the log messages.

        Returns:
        --------
            str: Formatted log string.

        """
        output = ""
        output += separator.join(self.logs)
        if self.exceptions:
            output += "\n\n"
            output += self._format_exceptions()
        output += "\n\n"
        return output

    def _format_exceptions(self) -> str:
        """
        Formats the exceptions as a string to be written to file or stdout.

        Returns:
        --------
            str: Formatted exceptions string.

        """
        exception_str = "\n\nEXCEPTIONS\n"
        exception_str += "----------\n\n"
        exception_str += "\n\n".join([e for e in self.exceptions])
        return exception_str


# from __future__ import absolute_import, division, print_function, unicode_literals


# class LoggingMixin(object):
#     """docstring for LoggingMixin"""

#     def __init__(self):
#         self._log = None
#         self._exceptions = None

#     @property
#     def logs(self):
#         if self._log is not None:
#             return self._log
#         return []

#     @property
#     def exceptions(self):
#         if self._exceptions is not None:
#             return self._exceptions
#         return []

#     def log(self, *args, **kwargs):
#         """
#         Records a log message
#         """
#         sep = kwargs.get("sep", " ")
#         lstring = sep.join([str(a) for a in args])
#         if self._log is None:
#             self._log = [
#                 lstring,
#             ]
#         else:
#             self._log.append(lstring)

#     def exception(self, *args, **kwargs):
#         """
#         Records an exception.
#         """
#         sep = kwargs.get("sep", "\n")
#         estring = sep.join([str(a) for a in args])
#         if self._exceptions is None:
#             self._exceptions = [
#                 estring,
#             ]
#         else:
#             self._exceptions.append(estring)

#     def format_log(self):
#         """
#         Formats the log, including exceptions.

#         Log formatting will only be performed on sequences that had an
#         error during annotation, unless AbStar is run in debug
#         mode. In debug mode, all sequences will be logged.

#         Returns:
#         --------

#             str: Formatted log string.
#         """
#         output = ""
#         output += "\n".join(self.logs)
#         if self._check_for_exceptions():
#             output += "\n\n"
#             output += self._format_exceptions()
#         output += "\n\n"
#         return output

#     def _check_for_exceptions(self):
#         if self.exceptions:
#             return True
#         if self.v is not None:
#             if self.v._exceptions:
#                 return True
#         if self.d is not None:
#             if self.d._exceptions:
#                 return True
#         if self.j is not None:
#             if self.j._exceptions:
#                 return True
#         return False

#     def _format_exceptions(self):
#         exceptions = []
#         exceptions += self.exceptions
#         estring = "\n\nEXCEPTIONS\n"
#         estring += "----------\n\n"
#         if self.v is not None:
#             exceptions += self.v.exceptions
#         if self.d is not None:
#             exceptions += self.d.exceptions
#         if self.j is not None:
#             exceptions += self.j.exceptions
#         estring += "\n\n".join([e for e in exceptions])
#         return estring
