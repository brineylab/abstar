#!/usr/bin/python
# filename: run_add_task.py


###########################################################################
#
# Copyright (c) 2014 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT) 
#
###########################################################################


import sys

from utils.celery.tasks import add


if __name__ == '__main__':
	res = []
	num = int(sys.argv[1]) if len(sys.argv) > 1 else 10
	for x, y in zip(range(num), range(num)):
		res.append(add.delay(x, y))
	results = [r.get() for r in res]
	print results