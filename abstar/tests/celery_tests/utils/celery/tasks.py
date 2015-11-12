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


from utils.celery.celery import celery


@celery.task
def add(x, y):
	return x + y

