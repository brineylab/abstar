#!/usr/bin/python
# filename: tasks.py


###########################################################################
#
# Copyright (c) 2014 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT) 
#
###########################################################################


from utils.vdj import run
from utils.celery.celery import celery


@celery.task
def run_vdj(*args):
	return run(*args)
		


