#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: dbbuddy_tests.py
version: 1, alpha
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Collection of PyTest unit tests for the DatabaseBuddy.py program
"""

import pytest
import os
import sys
from hashlib import md5

sys.path.insert(0, "./")
import DatabaseBuddy as Db
import buddy_resources as br
import MyFuncs


def result_to_hash(_dbbuddy, mode='hash'):
    if mode != "hash":
        return str(_dbbuddy)
    _hash = md5("{0}\n".format(str(_dbbuddy).rstrip()).encode()).hexdigest()
    return _hash

root_dir = os.getcwd()


def resource(file_name):
    return "{0}/unit_test_resources/{1}".format(root_dir, file_name)


accession_files = ["genbank.csv", "genbank.tsv"]
accession_types = ["genbank", "genbank"]

input_tuples = [(next_file, accession_types[indx]) for indx, next_file in enumerate(accession_files)]


@pytest.mark.parametrize("align_file,file_type", input_tuples)
def test_instantiate_dbbuddy_from_file(align_file, file_type):
    assert type(Db.DbBuddy(resource(align_file), _in_format=file_type)) == Db.DbBuddy
