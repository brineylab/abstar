from setuptools import setup

config = {
	'description': 'VDJ assignment and antibody sequence annotation. Scalable from a single sequence to billions of sequences.',
	'author': 'Bryan Briney',
	'url': 'www.github.com/briney/abstar/',
	# 'download_url': 'www.github.com/briney/abstar/',
	'author_email': 'briney@scripps.edu',
	'version': '0.1.1',
	'install_requires': ['abtools',
						 'biopython',
						 'celery',
						 'nwalign',
						 'pymongo',
						 'scikit-bio'],
	'packages': ['abstar'],
	'scripts': ['bin/abstar',
				'bin/batch_mongoimport',
				'bin/make_basespace_credfile'],
	'name': 'abstar',
	'include_package_data': True
}

setup(**config)
