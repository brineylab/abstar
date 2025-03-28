# Copyright (c) 2024 brineylab @ scripps
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


import os

from setuptools import find_packages, setup

# the actual __version__ is read from version.py, not assigned directly
# this causes the linter to freak out when we mention __version__ in setup()
# to fix that, we fake assign it here
__version__ = None

# read version
version_file = os.path.join(os.path.dirname(__file__), "abstar", "version.py")
with open(version_file) as version:
    exec(version.read())

# read requirements
with open("requirements.txt") as reqs:
    requirements = reqs.read().splitlines()

# read long description
with open("README.md", encoding="utf-8") as readme:
    long_description = readme.read()

setup(
    name="abstar",
    version=__version__,
    author="Bryan Briney",
    author_email="briney@scripps.edu",
    description="Germline assignment and annotation of adaptive immune receptor repertoire (AIRR) data. Scalable from a single sequence to billions of sequences.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/brineylab/abstar",
    packages=find_packages(),
    include_package_data=True,
    # scripts=[
    #     "bin/abstar",
    #     # "bin/batch_mongoimport",
    #     # "bin/basespace_downloader",
    #     # "bin/build_abstar_germline_db",
    #     # "bin/make_basespace_credfile",
    # ],
    entry_points={
        "console_scripts": [
            "abstar=abstar.scripts.abstar:cli",
        ]
    },
    install_requires=requirements,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.10",
)


# import os
# import sys

# try:
#     from setuptools import setup
# except ImportError:
#     sys.exit('ERROR: setuptools is required.\n')

# try: # for pip >= 10
#     from pip._internal.req import parse_requirements
# except ImportError: # for pip <= 9.0.3
#     from pip.req import parse_requirements
# # try:
# #     from pip.req import parse_requirements
# # except ImportError:
# #     sys.exit('ERROR: pip is required.\n')


# if os.environ.get('READTHEDOCS', None):
#     # Set empty install_requires to get install to work on readthedocs
#     install_requires = []
# else:
#     req_file = 'requirements.txt'
#     try:
#         reqs = parse_requirements(req_file, session=False)
#     except TypeError:
#         reqs = parse_requirements(req_file)
#     try:
#         install_requires = [str(r.req) for r in reqs]
#     except AttributeError:
#         install_requires = [str(r.requirement) for r in reqs]

# # read version
# exec(open('abstar/version.py').read())

# # read long description
# with open("README.md", "r") as fh:
#     long_description = fh.read()


# config = {
#     'name': 'abstar',
#     'version': __version__,
#     'author': 'Bryan Briney',
#     'author_email': 'briney@scripps.edu',
#     'description': 'VDJ assignment and antibody sequence annotation. Scalable from a single sequence to billions of sequences.',
#     'long_description': long_description,
#     'long_description_content_type': 'text/markdown',
#     'url': 'https://www.github.com/briney/abstar',
#     'install_requires': install_requires,
#     'packages': ['abstar'],
#     'scripts': ['bin/abstar',
#                 'bin/batch_mongoimport',
#                 'bin/basespace_downloader',
#                 'bin/build_abstar_germline_db',
#                 'bin/make_basespace_credfile',],
#     'include_package_data': True,
#     'classifiers': ['License :: OSI Approved :: MIT License',
#                     'Programming Language :: Python :: 3.8',
#                     'Programming Language :: Python :: 3.9',
#                     'Topic :: Scientific/Engineering :: Bio-Informatics']
# }

# setup(**config)
