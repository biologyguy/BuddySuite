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

"""
# ###################### 'dt', '--display_trees' ###################### #
@pytest.mark.display
def test_display_trees(monkeypatch):
    show = mock.Mock(return_value=True)
    monkeypatch.setattr(ete3.TreeNode, "show", show)
    assert Pb.display_trees(pb_resources.get_one("o k"))


def test_display_trees_error():
    # noinspection PyUnresolvedReferences
    with mock.patch.dict('os.environ'):
        if 'DISPLAY' in os.environ:
            del os.environ['DISPLAY']
        with pytest.raises(SystemError):
            Pb.display_trees(Pb.make_copy(pb_objects[0]))
"""

# ###################### 'dis', '--distance' ###################### #
hashes = [('m k', 'aed02a04021dcbc99d41e34592ebeed0'), ('m n', 'cc23d0f2f6ae9b5a14425caef2d8ccb2'),
          ('m l', 'aed02a04021dcbc99d41e34592ebeed0')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_distance_wrf(key, next_hash, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    output = str(Pb.distance(tester, method='wrf'))
    assert pb_helpers.string2hash(output) == next_hash

hashes = [('m k', '49173bc33d89cc3d912a6af0fd51801d'), ('m n', '56574d4305bb5f094363a0fad351fd42'),
          ('m l', '49173bc33d89cc3d912a6af0fd51801d')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_distance_uwrf(key, next_hash, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    output = str(Pb.distance(tester, method='uwrf'))
    assert pb_helpers.string2hash(output) == next_hash, print(output)

hashes = [('m k', 'bae2c660250d42d6ba9bac7d311d6ffb'), ('m n', '0596ee4e2d4e76b27fe20b66a8fbea51'),
          ('m l', 'bae2c660250d42d6ba9bac7d311d6ffb')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_distance_ed(key, next_hash, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    output = str(Pb.distance(tester, method='ed'))
    assert pb_helpers.string2hash(output) == next_hash, print(output)


def test_distance_unknown_method(pb_resources):
    with pytest.raises(AttributeError) as err:
        Pb.distance(pb_resources.get_one("m k"), method='foo')
        assert "foo is an invalid comparison method." in err
