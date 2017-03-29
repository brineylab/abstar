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

import argparse
import json
import os
import platform
import sys
import time

from BaseSpacePy.api.BaseSpaceAPI import BaseSpaceAPI
from BaseSpacePy.model.QueryParameters import QueryParameters as qp

from abtools import log
from abtools.pipeline import make_dir

logger = log.get_logger('basespace')


class BaseSpace(object):
    def __init__(self, project_id=None, project_name=None, get_all_projects=False):
        super(BaseSpace, self).__init__()
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
            self.project_id = None
            self.project_name = None
            # self.project_id, self.project_name = self._user_selected_project_id()
        self._runs = None


    @property
    def runs(self):
        if self._runs is None:
            self._runs = self.api.getAccessibleRunsByUser(queryPars=self.params)
        return self._runs


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
        if all([self.project_id is None, self.project_name is None]):
            self.project_id, self.project_name = self._user_selected_project_id()
        files = self._get_files()
        self.print_download_info(files)
        start = time.time()
        for i, f in enumerate(files):
            # self.log.write('[ {} ] {}\n'.format(i, str(f)))
            logger.info('[ {} ] {}'.format(i, str(f)))
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


def parse_args():
    parser = argparse.ArgumentParser("Downloads sequencing data from BaseSpace, Illumina's cloud storage platform.")
    parser.add_argument('-d', '--download-directory',
                        dest='download_directory',
                        required=True,
                        help="Directory into which BaseSpace data will be downloaded.")
    parser.add_argument('--project-id',
                        default=None,
                        help='ID of the project to be downloaded. Optional.')
    parser.add_argument('--project-name',
                        default=None,
                        help='Name of the project to be downloaded. Optional.')
    args = parser.parse_args()
    return args


def download(download_directory, project_id=None, project_name=None):
    '''
    Downloads sequencing data from BaseSpace (Illumina's cloud storage platform).

    Before accessing BaseSpace through the AbStar API, you need to set up a
    credentials file:

    1. You need a BaseSpace access token. The easiest way to do this is to
       set up a BaseSpace developer account following
       `these instructions <https://support.basespace.illumina.com/knowledgebase/articles/403618-python-run-downloader>`_

    2. Make a BaseSpace credentials file using your developer credentials::

        $ make_basespace_credfile

    and follow the instructions.


    Examples:

        If you know the name of the project you'd like to download::

            from abstar.utils import basespace

            basespace.download('/path/to/download_directory', project_name='MyProject')

        If you know the ID of the project you'd like to download::

            basespace.download('/path/to/download_directory', project_id='ABC123')

        If neither ``project_id`` nor ``project_name`` is provided, a list of your available
        BaseSpace projects will be provided and you can select a project from that list::

            basespace.download('/path/to/download_directory')

    Args:

        download_directory (str): Directory into which the raw sequences files should
            be downloaded. If the directory does not exist, it will be created.

        project_id (str): ID of the project to be downloaded.

        project_name (str): Name of the project to be downloaded.

    Returns:

        int: The number of sequence files downloaded.
    '''
    make_dir(download_directory)
    bs = BaseSpace(project_id, project_name)
    return bs.download(download_directory)


if __name__ == '__main__':
    args = parse_args()
    download(args.download_directory,
             project_id=args.project_id,
             project_name=args.project_name)
