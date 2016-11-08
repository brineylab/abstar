import os

try:
    from setuptools import setup
except ImportError:
    sys.exit('ERROR: setuptools is required.\n')

try:
    from pip.req import parse_requirements
except ImportError:
    sys.exit('ERROR: pip is required.\n')


if os.environ.get('READTHEDOCS', None):
    # Set empty install_requires to get install to work on readthedocs
    install_requires = []
else:
    req_file = 'requirements.txt'
    try:
        reqs = parse_requirements(req_file, session=False)
    except TypeError:
        reqs = parse_requirements(req_file)
    install_requires = [str(r.req) for r in reqs]


config = {
    'description': 'VDJ assignment and antibody sequence annotation. Scalable from a single sequence to billions of sequences.',
    'author': 'Bryan Briney',
    'url': 'https://www.github.com/briney/abstar',
    # 'download_url': 'www.github.com/briney/abstar/',
    'author_email': 'briney@scripps.edu',
    'version': '0.1.10',
    'install_requires': install_requires,
    'packages': ['abstar'],
    'scripts': ['bin/abstar',
                'bin/batch_mongoimport',
                'bin/make_basespace_credfile'],
    'name': 'abstar',
    'include_package_data': True
}

setup(**config)
