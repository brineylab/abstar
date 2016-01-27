#!/usr/bin/python
# filename: make_basespace_credfile.py

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


from __future__ import print_function, absolute_import

import json
import platform
import os
import sys


def main():
    creds = {}
    print('\nEnter your BaseSpace credentials')
    creds['client_id'] = raw_input('Client ID: ')
    if not creds['client_id']:
        print('ERROR: Client ID is a required field.')
        sys.exit(1)
    creds['client_secret'] = raw_input('Client Secret: ')
    if not creds['client_secret']:
        print('ERROR: Client Secret is a required field.')
        sys.exit(1)
    creds['access_token'] = raw_input('Access Token: ')
    if not creds['access_token']:
        print('ERROR: Access Token is a required field.')
        sys.exit(1)
    print("\nThe following options should only be provided if you know what you're doing.")
    print("If unsure, leave blank to use the default value.")
    creds['api_server'] = raw_input('API Server [https://api.basespace.illumina.com/]: ')
    if not creds['api_server'].strip():
        creds['api_server'] = 'https://api.basespace.illumina.com/'
    creds['version'] = raw_input('Version [v1pre3]: ')
    if not creds['version'].strip():
        creds['version'] = 'v1pre3'
    write_cred_file(creds)


def write_cred_file(creds):
    jcreds = json.dumps(creds)
    cred_file = os.path.expanduser('~/.abstar/basespace_credentials')
    if not os.path.exists(os.path.dirname(cred_file)):
        os.makedirs(os.path.dirname(cred_file))
    cred_handle = open(cred_file, 'w')
    cred_handle.write(jcreds)
    print('\nCredentials have been saved to {}\n'.format(cred_file))


if __name__ == '__main__':
    main()
