try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'AbStar',
    'author': 'Bryan Briney',
    'url': 'www.github.com/briney/abstar/',
    'download_url': 'www.github.com/briney/abstar/',
    # 'author_email': 'My email.',
    'version': '0.1.0',
    'install_requires': ['nose', 'biopython', 'celery', 'nwalign'],
    'packages': ['abstar'],
    'scripts': [],
    'name': 'abstar',
    'include_package_data': True
}

setup(**config)
