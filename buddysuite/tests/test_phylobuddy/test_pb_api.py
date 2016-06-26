#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Tests PhyloBuddy API functions """
import pytest
try:
    from buddysuite import PhyloBuddy as Pb
    from buddysuite import AlignBuddy as Alb
    from buddysuite import MyFuncs
except ImportError:
    import PhyloBuddy as Pb
    import AlignBuddy as Alb
    import MyFuncs

# ###################### 'ct', '--consensus_tree' ###################### #
hashes = [('m k', 'acd3fb34cce867c37684244701f9f5bf'), ('m n', 'eede64c804e531cb1c99e4240589b04b'),
          ('m l', '73ac98a1656d1c4a52da16d3f096f8ce'), ('o k', '64f7df66253b104c300d13e344e2f216')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_consensus_tree(key, next_hash, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    tester = Pb.consensus_tree(tester)
    assert pb_helpers.phylo_to_hash(tester) == next_hash, tester.write("error_files/%s" % next_hash)

hashes = [('m k', 'baf9b2def2c0fa2ff97d3a16d24b4738'), ('m n', '3ab21ed1202222f17a44fff7b4051aa0'),
          ('m l', 'd680eece99ad907bbc8cd69ceaeb7b8a'), ('o k', '64f7df66253b104c300d13e344e2f216')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_consensus_tree_95(key, next_hash, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    tester = Pb.consensus_tree(tester, frequency=0.95)
    assert pb_helpers.phylo_to_hash(tester) == next_hash
