import os
import sys

try:
    from setuptools import setup
except ImportError:
    sys.exit('ERROR: setuptools is required.\n')

try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements
# try:
#     from pip.req import parse_requirements
# except ImportError:
#     sys.exit('ERROR: pip is required.\n')


if os.environ.get('READTHEDOCS', None):
    # Set empty install_requires to get install to work on readthedocs
    install_requires = []
else:
    if sys.version_info[0] > 2:
        req_file = 'requirements.txt'
    else:
        req_file = 'requirements2.txt'
    try:
        reqs = parse_requirements(req_file, session=False)
    except TypeError:
        reqs = parse_requirements(req_file)
    try:
        install_requires = [str(r.req) for r in reqs]
    except AttributeError:
        install_requires = [str(r.requirement) for r in reqs]

# read version
exec(open('abstar/version.py').read())

# read long description
with open("README.md", "r") as fh:
    long_description = fh.read()


config = {
    'name': 'abstar',
    'version': __version__,
    'author': 'Bryan Briney',
    'author_email': 'briney@scripps.edu',
    'description': 'VDJ assignment and antibody sequence annotation. Scalable from a single sequence to billions of sequences.',
    'long_description': long_description,
    'long_description_content_type': 'text/markdown',
    'url': 'https://www.github.com/briney/abstar',
    'install_requires': install_requires,
    'packages': ['abstar'],
    'scripts': ['bin/abstar',
                'bin/batch_mongoimport',
                'bin/basespace_downloader',
                'bin/build_abstar_germline_db',
                'bin/make_basespace_credfile',],
    'include_package_data': True,
    'classifiers': ['License :: OSI Approved :: MIT License',
                    'Programming Language :: Python :: 3.6',
                    'Programming Language :: Python :: 3.7',
                    'Programming Language :: Python :: 3.8',
                    'Topic :: Scientific/Engineering :: Bio-Informatics']
}

setup(**config)
