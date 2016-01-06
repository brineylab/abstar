#!/usr/bin/python
# filename: basespace.py

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
import os
import platform
import sys
import time

from BaseSpacePy.api.BaseSpaceAPI import BaseSpaceAPI
from BaseSpacePy.model.QueryParameters import QueryParameters as qp

from abtools.utils import log
logger = log.get_logger('basespace')


class BaseSpace(object):
	"""docstring for BaseSpace"""
	def __init__(self, project_id=None, project_name=None, undetermined=False):
		super(BaseSpace, self).__init__()
		# self.log = log
		# BaseSpace credentials
		creds = self._get_credentials()
		self.client_key = creds['client_id']
		self.client_secret = creds['client_secret']
		self.access_token = creds['access_token']
		self.version = creds['version']
		self.api_server = creds['api_server']
		self.api = BaseSpaceAPI(self.client_key, self.client_secret, self.api_server, self.version, AccessToken=self.access_token)
		self.params = qp(pars={'Limit': 1024, 'SortDir': 'Desc'})
		if project_id is not None:
			self.project_id = project_id
			self.project_name = None
		elif project_name is not None:
			self.project_name = project_name
			self.project_id = self._get_project_id_from_name(project_name)
		else:
			self.project_id, self.project_name = self._user_selected_project_id()


	def _get_credentials(self):
		# BaseSpace credentials file should be in JSON format
		cred_file = os.path.expanduser('~/.abstar/basespace_credentials')
		cred_handle = open(cred_file, 'r')
		return json.load(cred_handle)


	def _get_project_id_from_name(self):
		projects = self.api.getProjectByUser(queryPars=self.params)
		for project in projects:
			name = project.Name.encode('ascii', 'ignore')
			if name == self.project_name:
				return project.Id
		print('No projects matched the given project name ({})'.format(name))
		sys.exit(1)


	def _user_selected_project_id(self):
		projects = self.api.getProjectByUser(queryPars=self.params)
		self.print_basespace_project()
		offset = 0
		while True:
			for i, project in enumerate(projects[offset * 25:(offset * 25) + 25]):
				project_name = project.Name.encode('ascii', 'ignore')
				print('[ {} ] {}'.format(i + (offset * 25), project_name))
			print('')
			project_index = raw_input("Select the project number (or 'next' to see more projects): ")
			try:
				project_index = int(project_index)
				return projects[project_index].Id, projects[project_index].Name.encode('ascii', 'ignore')
			except:
				offset += 1
		return projects[project_index].Id, projects[project_index].Name.encode('ascii', 'ignore')


	def _get_projects(self, start=0):
		projects = self.api.getProjectByUser(queryPars=self.params)
		self.print_basespace_project()
		for i, project in enumerate(projects[:25]):
			project_name = project.Name.encode('ascii', 'ignore')
			print('[ {} ] {}'.format(i, project_name))
		print('')
		return projects


	def _get_samples(self, project_id):
		samples = []
		offset = 0
		while True:
			query_params = qp(pars={'Limit': 1024, 'SortDir': 'Asc', 'Offset': offset * 1024})
			s = self.api.getSamplesByProject(self.project_id, queryPars=query_params)
			if not s:
				break
			samples.extend(s)
			offset += 1
		return samples


	def _get_files(self):
		files = []
		samples = self._get_samples(self.project_id)
		for sample in samples:
			files.extend(self.api.getFilesBySample(sample.Id, queryPars=self.params))
		return files


	def download(self, direc):
		files = self._get_files()
		self.print_download_info(files)
		start = time.time()
		for i, f in enumerate(files):
			# self.log.write('[ {} ] {}\n'.format(i, str(f)))
			logger.info('[ {} ] {}\n'.format(i, str(f)))
			f.downloadFile(self.api, direc)
		end = time.time()
		self.print_completed_download_info(start, end)
		return len(files)


	def print_basespace_project(self):
		print('')
		print('')
		print('========================================')
		print('BaseSpace Project Selection')
		print('========================================')
		print('')


	def print_download_info(self, files):
		logger.info('')
		logger.info('')
		logger.info('========================================')
		logger.info('Downloading files from BaseSpace')
		logger.info('========================================')
		logger.info('')
		logger.info('Identified {0} files for download.'.format(len(files)))
		logger.info('')


	def print_completed_download_info(self, start, end):
		logger.info('')
		logger.info('Download completed in {0} seconds'.format(end - start))


# def _setup_logging():
# 	try:
# 		from abtools.utils import log
# 		global logger
# 		logger = log.get_logger('basespace')
# 	except ImportError:
# 		import logging
# 		fmt = '[%(levelname)s] %(name)s %(asctime)s %(message)s'
# 		logging.basicConfig()
# 		global logger
# 		logger = logging.getLogger('basespace')
# 		formatter = logging.Formatter("%(message)s")
# 		ch = logging.StreamHandler()
# 		ch.setFormatter(formatter)
# 		ch.setLevel(logging.INFO)
# 		logger.addHandler(ch)


def download(direc, project_id=None, project_name=None, undetermined=False):
	# _setup_logging()
	bs = BaseSpace(project_id, project_name, undetermined)
	return bs.download(direc)


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser("Parses the output of IgBLAST into something suitable for import into a MySQL database")
	parser.add_argument('-i', '--in',
						dest='input',
						required=True,
						help="The input file, to be split and processed in parallel. \
						If a directory is given, all files in the directory will be iteratively processed.")
	args = parser.parse_args()
	download(args.input)
