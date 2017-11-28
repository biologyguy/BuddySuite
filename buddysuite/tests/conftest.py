#!/usr/bin/env python3
# coding=utf-8
""" Fixtures for py.test  """
import pytest
import sys
import os

from .__init__ import RESOURCE_PATH
from . import __init__ as init
from .. import buddy_resources as br


# #################################  -  Helper functions  -  ################################## #
@pytest.fixture(scope="session")
def hf():
    """
    Collection of helper methods
    """
    return init.HelperMethods()


# #################################  -  SeqBuddy  -  ################################## #
@pytest.fixture(scope="session")
def sb_resources():
    """
    Create a shared Resources object of alignment files and Alb objects
    """
    resources_obj = init.SbResources()
    return resources_obj


@pytest.fixture(scope="session")
def sb_odd_resources():
    # A dict of invalid file resources
    resource_list = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
        ('blank', '{path}blank.fa'),
        ('circular', '{path}circular.gb'),
        ('circular_digest', '{path}circular_digest.gb'),
        ('figtree', '{path}figtree.nexus'),
        ('unrecognizable', '{path}unrecognizable.phy'),
        ('gibberish', '{path}gibberish.fa'),
        ('phylipss_cols', '{path}malformed_phylip_columns.physs'),
        ('duplicate', '{path}Duplicate_seqs.fa'),
        ('ambiguous_dna', '{path}ambiguous_dna.fa'),
        ('ambiguous_rna', '{path}ambiguous_rna.fa'),
        ('blastn', '{path}blast/Mnemiopsis_cds.n'),
        ('blastp', '{path}blast/Mnemiopsis_pep.p'),
        ('dummy_feats', '{path}Mnemiopsis_cds_dummy_features.gb'),
        ('cnidaria_pep', '{path}Cnidaria_pep.nexus'),
        ('mixed', '{path}mixed_alpha.fa')
    ]}
    return resource_list


# #################################  -  AlignBuddy  -  ################################ #
@pytest.fixture(scope="session")
def alb_resources():
    """
    Create a shared Resources object of alignment files and Alb objects
    """
    resources_obj = init.AlbResources()
    return resources_obj


@pytest.fixture(scope="session")
def alb_odd_resources():
    """ A dict of invalid file resources """
    resource_list = {
        'dna': {
            'single': {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
                ('fasta', '{path}gibberish.fa'),
                ('phylipss_recs', '{path}malformed_phylip_records.physs'),
                ('phylipss_cols', '{path}malformed_phylip_columns.physs'),
                ('ambiguous', '{path}ambiguous_dna_alignment.fa')]}
        },
        'protein': {
            'single': {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
                ('phylip', '{path}unrecognizable.phy')]}
        },
        'blank': "%sblank.fa" % RESOURCE_PATH
    }
    return resource_list


# ################################  -  PhyloBuddy  -  ################################# #
@pytest.fixture(scope="session")
def pb_resources():
    """
    Create a shared Resources object of alignment files and Alb objects
    """
    resources_obj = init.PbResources()
    return resources_obj


@pytest.fixture(scope="session")
def pb_odd_resources():
    # A dict of invalid file resources
    resource_list = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
        ('blank', '{path}blank.fa'),
        ('unrecognizable', '{path}unrecognizable.phy'),
        ('figtree', '{path}figtree.nexus'),
        ('compare', '{path}compare_trees.newick'),
        ('node_lables', '{path}tree_with_node_lables.nwk'),
        ('lengths', '{path}Mnemiopsis_pep.newick'),
        ('bootstraps', '{path}Mnemiopsis_pep_bootstraps.newick'),
        ('support', '{path}Mnemiopsis_pep_support.newick')
    ]}
    return resource_list
