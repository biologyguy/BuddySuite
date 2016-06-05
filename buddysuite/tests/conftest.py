#!/usr/bin/env python3
# coding=utf-8
""" Fixtures for py.test  """
import pytest
from hashlib import md5
from tests import Resources, RESOURCE_PATH


@pytest.fixture(scope="session")
def alb_resources():
    """
    Create a shared Resources object of alignment files and Alb objects
    """

    resources_obj = Resources()
    return resources_obj


@pytest.fixture(scope="session")
def alignment_bad_resources():
    """ A dict of invalid file resources """
    resource_list = {
        'dna': {
            'single': {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
                ('fasta', '{path}/gibberish.fa'), ('phylipss_recs', '{path}/malformed_phylip_records.physs'),
                ('phylipss_cols', '{path}/malformed_phylip_columns.physs')]}
        },
        'protein': {
            'single': {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
                ('phylip', '{path}/unrecognizable.phy')]}
        },
        'blank': "%s/blank.fa" % RESOURCE_PATH
    }
    return resource_list


@pytest.fixture(scope="session")
def helpers():
    class Helpers(object):
        def __init__(self):
            self.resource_path = RESOURCE_PATH

        @staticmethod
        def align_to_hash(alignbuddy=None, mode='hash'):
            if not alignbuddy:
                raise AttributeError("AlignBuddy object required")

            if mode != "hash":
                return str(alignbuddy)
            _hash = md5("{0}".format(str(alignbuddy)).encode("utf-8")).hexdigest()
            return _hash

        @staticmethod
        def string2hash(_input):
            return md5(_input.encode("utf-8")).hexdigest()

    helper_obj = Helpers()
    return helper_obj
