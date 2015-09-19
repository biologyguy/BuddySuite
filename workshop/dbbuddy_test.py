#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
