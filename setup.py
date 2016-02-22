#!/usr/bin/env python

from distutils.core import setup

setup(name='molniac',
      version='0.1',
      description='Manipulate macromolecular structures',
      author='Oliver Schillinger',
      author_email='oliver@schillinger-ttp.de',
      url='https://github.com/schilli/molniac',
      packages=['molniac'],
      package_data={'molniac': ['data/*']}
     )

