#!/usr/bin/python
# filename: celeryconfig.py


###########################################################################
#
# Copyright (c) 2014 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################



# config file for Celery Daemon

# RabbitMQ broker
broker_url = 'pyamqp://abcloud:abcloud@master:5672/abcloud_host'

# Redis broker
# BROKER_URL = 'redis://master:6379/0'

# RabbitMQ backend
result_backend = 'rpc://abcloud:abcloud@master:5672/abcloud_host'

# Redis backend
# CELERY_RESULT_BACKEND = 'redis://master:6379/0'

# Additional Redis-specific configs
broker_transport_options = {'fanout_prefix': True,
                            'fanout_patterns': True,
                            'visibility_timeout': 3600}

# Other configs
worker_max_tasks_per_child = 320
worker_prefetch_multiplier = 1
task_acks_late = True
