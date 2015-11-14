from setuptools import setup

config = {
	'description': 'AbStar',
	'author': 'Bryan Briney',
	'url': 'www.github.com/briney/abstar/',
	# 'download_url': 'www.github.com/briney/abstar/',
	'author_email': 'briney@scripps.edu',
	'version': '0.1.0',
	'install_requires': ['nose', 'biopython', 'celery', 'nwalign', 'pymongo'],
	'packages': ['abstar'],
	'scripts': ['bin/abstar'],
	'name': 'abstar',
	'include_package_data': True
}

setup(**config)
