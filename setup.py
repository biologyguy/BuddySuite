#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, Command
import os


class CleanCommand(Command):
    """
    Custom clean command to tidy up the project root.
    http://bit.ly/2bw7xXb
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    @staticmethod
    def run():
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')

PACKAGES = [
    'buddysuite',
]

DEPENDENCIES = [
    'pytest',
    'pytest-xdist',
    'pytest-cov',
    'numpy',
    'scipy',
    'six',
    'dill',
    'suds-py3',
    'biopython',
    'dendropy',
    'ete3',
]

CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: Public Domain',
    'Natural Language :: English',
    'Operating System :: MacOS',
    'Operating System :: POSIX :: Linux',
    'Programming Language :: Python :: 3 :: only',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Medical Science Apps.',
    'Topic :: Software Development :: Libraries :: Python Modules'
]

KEYWORDS = [
    'computational biology',
    'computational genetics',
    'genetics',
    'genome',
    'phylogenetics',
    'biology',
    'bioinformatics'
]

ENTRY_POINTS = {
    'console_scripts': [
        'alignbuddy = buddysuite.AlignBuddy:main',
        'databasebuddy = buddysuite.DatabaseBuddy:main',
        'phylobuddy = buddysuite.PhyloBuddy:main',
        'seqbuddy = buddysuite.SeqBuddy:main',
        'buddysuite = buddysuite.BuddySuite:main'
    ]
}

setup(name='buddysuite',
      version='1.2b.7',
      description='BuddySuite is a collection of command line utilities written in Python for '
                  'working with biological data.',
      long_description=open(os.path.join(os.path.dirname(__file__), 'README.rst'), encoding="utf-8").read(),
      author='Stephen Bond',
      maintainer='Stephen Bond',
      author_email='steve.bond@gmail.com',
      url='https://github.com/biologyguy/BuddySuite',
      packages=PACKAGES,
      setup_requires=['numpy'],
      install_requires=DEPENDENCIES,
      entry_points=ENTRY_POINTS,
      license='Public Domain',
      keywords=KEYWORDS,
      zip_safe=False,
      cmdclass={'clean': CleanCommand},
      )
