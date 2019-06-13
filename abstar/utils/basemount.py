#!/usr/bin/python
# filename: basemount.py

#
# Copyright (c) 2018 Bryan Briney
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

from datetime import datetime
from glob import glob
import os
import shutil
import stat
import sys

import arrow

from abutils.utils.decorators import lazy_property
from abutils.utils.pipeline import make_dir
from abutils.utils.progbar import progress_bar

if sys.version_info[0] > 2:
    raw_input = input


class BasemountProject():
    '''

    '''
    def __init__(self, parent_dir):
        self.parent_dir = os.path.abspath(parent_dir)
        self.name = os.path.basename(self.parent_dir)
    

    @lazy_property
    def date(self):
        app_sessions_dir = os.path.join(self.parent_dir, 'AppSessions')
        app_sessions = glob(os.path.join(app_sessions_dir, 'FASTQ Generation*'))
        try:
            app_session = sorted(app_sessions)[0]
            date = arrow.get(os.path.basename(app_session).lstrip('FASTQ Generation '))
        except:
            date = None
        return date

    @lazy_property
    def fastqs(self):
        fastq_dir = os.path.join(self.parent_dir, 'Samples')
        fastqs = []
        for subdir in glob(os.path.join(fastq_dir, '*')):
            files_subdir = os.path.join(subdir, 'Files')
            if os.path.exists(files_subdir):
                fastqs += glob(os.path.join(files_subdir, '*.fastq.gz'))
        return sorted(fastqs)


    @property
    def num_fastqs(self):
        return len(self.fastqs)

    
    def copy_fastqs(self, destination_dir, show_progress=True):
        make_dir(destination_dir)
        if show_progress:
            start = datetime.now()
            progress_bar(0, self.num_fastqs, start_time=start)
        for i, fastq_file in enumerate(self.fastqs, 1):
            shutil.copy(fastq_file, destination_dir)
            dest_file = os.path.join(destination_dir, os.path.basename(fastq_file))
            try:
                os.chmod(dest_file, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
            except PermissionError:
                continue
            if show_progress:
                is_complete = i == self.num_fastqs
                progress_bar(i, self.num_fastqs, start_time=start, complete=is_complete)



def copy_fastqs(basemount_dir, destination_dir, project_name=None):
    if project_name is not None:
        try:
            project = get_basemount_project(basemount_dir, project_name)
        except IndexError:
            print('The provided project name cannot be found. Please select a project:\n')
            projects = get_basemount_projects(basemount_dir)
            project = select_project(basemount_dir, projects)
    else:
        projects = get_basemount_projects(basemount_dir)
        project = select_project(basemount_dir, projects)
    if project is not None:
        project.copy_fastqs(destination_dir)
    else:
        error = "If you didn't see your project listed, refresh BaseMount "
        error += "(unmount and remount) and try again."
        print(error)
        print('\n')
        sys.exit(1)


def print_project_selection_header():
    print('')
    print('')
    print('========================================')
    print('BaseMount Project Selection')
    print('========================================')
    print('')


def select_project(basemount_dir, projects):
    print_project_selection_header()
    print('Retrieving project information (this may take a minute)...', end='')
    projects = sorted(projects, key=lambda x: x.date, reverse=True)
    print('\n')
    for chunk_start in range(0, len(projects), 25):
        chunk = projects[chunk_start:chunk_start + 25]
        for i, c in enumerate(chunk, chunk_start):
            print('[ {} ]  {}'.format(i, c.name))
        project_index = raw_input("Enter the project number (or 'next' to see more projects): ")
        try:
            project_index = int(project_index)
            return projects[project_index]
        except:
            print('')
            continue
    error = 'ERROR: There are no more projects.\n'
    error = "If you didn't see your project listed, refresh BaseMount "
    error += "(unmount and remount) and try again."
    print(error)
    print('\n')
    sys.exit(1)


def get_basemount_projects(basemount_dir):
    projects_dir = os.path.join(basemount_dir, 'Projects')
    projects = [BasemountProject(d) for d in glob(os.path.join(projects_dir, '*'))]
    projects = [p for p in projects if p.date is not None]
    return projects


def get_basemount_project(basemount_dir, project_name):
    projects_dir = os.path.join(basemount_dir, 'Projects')
    projects = [BasemountProject(d) for d in glob(os.path.join(projects_dir, '*'))]
    project = [p for p in projects if p.name.lower() == project_name.lower()]
    if project:
        return project[0]
    return None
