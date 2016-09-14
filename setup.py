#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, Command
from pkg_resources import Requirement, resource_filename
import os
import sys
import shutil
import random
import string
from configparser import ConfigParser, NoOptionError
from tempfile import TemporaryDirectory


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
    'biopython',
    'ete3',
    'six',
    'dendropy',
    'pytest',
    'scipy',
    'numpy',
    'dill',
    'suds-py3'
]

DATA_FILES = [
    ('config', ['setup_files/config.ini']),
    ('buddy_data', ['setup_files/buddysuite_usage.json', 'setup_files/cmd_history'])
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
    'biology'
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
      version='1.2b.0',
      description='BuddySuite is a collection of command line utilities written in Python for '
                  'working with biological data.',
      author='Stephen Bond',
      author_email='steve.bond@gmail.com',
      url='https://github.com/biologyguy/BuddySuite',
      packages=PACKAGES,
      install_requires=DEPENDENCIES,
      entry_points=ENTRY_POINTS,
      license='Public Domain',
      keywords=KEYWORDS,
      data_files=DATA_FILES,
      zip_safe=False,
      cmdclass={'setup': SetupCommand, 'clean': CleanCommand, 'uninstall': UninstallCommand, 'remove': UninstallCommand}
      )
