from setuptools import  find_packages
from distutils.core import setup

setup(
	name = 'DarkCapPy',

	packages = [
		'DarkCapPy',
		'DarkCapPy.Configure'
		   ],

	version = '2.0',

	description = 'Pyton functions for ArXiV:1509.07525v3',

	url = 'https://github.com/agree019/DarkCapPy',

	author = 'Adam Green',

	author_email = 'agree019@ucr.edu',

	include_package_data = True,

	license = 'MIT',

	install_requires = ['numpy','scipy','pandas'],

	download_url = 'https://github.com/agree019/DarkCapPy'
	)
