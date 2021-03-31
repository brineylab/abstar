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
BROKER_URL = 'amqp://abcloud:abcloud@master:5672/abcloud_host'

# Redis broker
# BROKER_URL = 'redis://master:6379/0'

# RabbitMQ backend
CELERY_RESULT_BACKEND = 'amqp://abcloud:abcloud@master:5672/abcloud_host'

# Redis backend
# CELERY_RESULT_BACKEND = 'redis://master:6379/0'

# Additional Redis-specific configs
BROKER_TRANSPORT_OPTIONS = {'fanout_prefix': True,
                            'fanout_patterns': True,
                            'visibility_timeout': 3600}

# Other configs
CELERYD_MAX_TASKS_PER_CHILD = 320
CELERYD_PREFETCH_MULTIPLIER = 1
CELERY_ACKS_LATE = True
