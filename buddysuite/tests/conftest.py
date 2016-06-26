#!/usr/bin/env python3
# coding=utf-8
""" Fixtures for py.test  """
import pytest
from hashlib import md5
import re

from __init__ import SbResources, AlbResources, PbResources, RESOURCE_PATH, br, MyFuncs


# #################################  -  SeqBuddy  -  ################################## #
@pytest.fixture(scope="session")
def sb_resources():
    """
    Create a shared Resources object of alignment files and Alb objects
    """
    resources_obj = SbResources()
    return resources_obj

"""
@pytest.fixture(scope="session")
def sb_bad_resources():
    # A dict of invalid file resources
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
"""


@pytest.fixture(scope="session")
def sb_helpers():
    class Helpers(object):
        def __init__(self):
            self.resource_path = RESOURCE_PATH
            self.write_file = MyFuncs.TempFile()

        def seqs_to_hash(self, _seqbuddy, mode='hash'):
            if _seqbuddy.out_format in ["gb", "genbank"]:
                for _rec in _seqbuddy.records:
                    try:
                        if re.search("(\. )+", _rec.annotations['organism']):
                            _rec.annotations['organism'] = "."
                    except KeyError:
                        pass

            if _seqbuddy.out_format == "phylipsr":
                self.write_file.write(br.phylip_sequential_out(_seqbuddy, relaxed=True, _type="seqbuddy"))
            elif _seqbuddy.out_format == "phylipss":
                self.write_file.write(br.phylip_sequential_out(_seqbuddy, relaxed=False, _type="seqbuddy"))
            else:
                _seqbuddy.write(self.write_file.path)

            seqs_string = "{0}\n".format(self.write_file.read().rstrip())
            self.write_file.clear()

            if mode != "hash":
                return seqs_string

            _hash = md5(seqs_string.encode()).hexdigest()
            return _hash

        @staticmethod
        def string2hash(_input):
            return md5(_input.encode("utf-8")).hexdigest()

    helper_obj = Helpers()
    return helper_obj


# #################################  -  AlignBuddy  -  ################################ #
@pytest.fixture(scope="session")
def alb_resources():
    """
    Create a shared Resources object of alignment files and Alb objects
    """
    resources_obj = AlbResources()
    return resources_obj


@pytest.fixture(scope="session")
def alb_bad_resources():
    """ A dict of invalid file resources """
    resource_list = {
        'dna': {
            'single': {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
                ('fasta', '{path}/gibberish.fa'),
                ('phylipss_recs', '{path}/malformed_phylip_records.physs'),
                ('phylipss_cols', '{path}/malformed_phylip_columns.physs'),
                ('ambiguous', '{path}/ambiguous_dna_alignment.fa')]}
        },
        'protein': {
            'single': {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
                ('phylip', '{path}/unrecognizable.phy')]}
        },
        'blank': "%s/blank.fa" % RESOURCE_PATH
    }
    return resource_list


@pytest.fixture(scope="session")
def alb_helpers():
    class Helpers(object):
        def __init__(self):
            self.resource_path = RESOURCE_PATH

        @staticmethod
        def align_to_hash(alignbuddy=None, mode='hash'):
            if not alignbuddy:
                raise AttributeError("AlignBuddy object required")

            if mode != "hash":
                return "{0}".format(str(alignbuddy))
            _hash = md5("{0}".format(str(alignbuddy)).encode("utf-8")).hexdigest()
            return _hash

        @staticmethod
        def string2hash(_input):
            return md5(_input.encode("utf-8")).hexdigest()

    helper_obj = Helpers()
    return helper_obj


# ################################  -  PhyloBuddy  -  ################################# #
@pytest.fixture(scope="session")
def pb_resources():
    """
    Create a shared Resources object of alignment files and Alb objects
    """
    resources_obj = PbResources()
    return resources_obj


@pytest.fixture(scope="session")
def pb_odd_resources():
    # A dict of invalid file resources
    resource_list = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
        ('blank', '{path}/blank.fa'),
        ('unrecognizable', '{path}/unrecognizable.phy'),
        ('figtree', '{path}/figtree.nexus'),
        ('compare', '{path}/compare_trees.newick')
    ]}
    return resource_list


@pytest.fixture(scope="session")
def pb_helpers():
    class Helpers(object):
        def __init__(self):
            self.resource_path = RESOURCE_PATH
            self.write_file = MyFuncs.TempFile()

        @staticmethod
        def phylo_to_hash(_phylobuddy, mode='hash'):
            if mode != "hash":
                return "{0}\n".format(str(_phylobuddy).rstrip())
            _hash = md5("{0}\n".format(str(_phylobuddy).rstrip()).encode('utf-8')).hexdigest()
            return _hash

        @staticmethod
        def string2hash(_input):
            return md5(_input.encode("utf-8")).hexdigest()

    helper_obj = Helpers()
    return helper_obj
# ###############################  -  DatabaseBuddy  -  ############################### #

