#!/usr/bin/env python
# filename: registry.py

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


# from .assigner import BaseAssigner
# from .blastn import Blastn

# import all assigner modules so that I can add them to the ASSIGNERS dict
from . import *

# ASSIGNERS = {'blastn': Blastn}
ASSIGNERS = {cls.__name__.lower(): cls for cls in vars()['BaseAssigner'].__subclasses__()}


# class Registration(type):
#     """docstring for Registration"""
#     def __init__(self, cls, name, bases, attrs):
#         ASSIGNERS[cls.__name__] = cls
#         super(Registration, cls).__init__(name, bases, attrs)


# def register(cls):
#     ASSIGNERS[cls.__name__.lower()] = cls
#     return cls
