#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

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
    'numpy'
]

DATA_FILES = [
    ('config', ['config.ini'])
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
        'seqbuddy = buddysuite.SeqBuddy:main'
    ]
}

setup(name='buddysuite',
      version='1.1.0',
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
      zip_safe=True)
