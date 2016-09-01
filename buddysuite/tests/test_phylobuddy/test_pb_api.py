#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Tests PhyloBuddy API functions """
import pytest
from ... import PhyloBuddy as Pb
from ... import buddy_resources as br

from unittest import mock
import ete3
import os
import shutil
import re


class MockPopen(object):
    def __init__(self, *args, **kwargs):
        return

    @staticmethod
    def wait():
        return True

    @staticmethod
    def communicate(output=""):
        return output


# ###################### 'cpt', '--collapse_polytomies' ###################### #
def test_collapse_polytomies(pb_odd_resources, pb_helpers):
    tester = Pb.PhyloBuddy(pb_odd_resources['support'])
    tester = Pb.collapse_polytomies(tester, 20)
    assert pb_helpers.phylo2hash(tester) == "1b0979265205b17ca7f34abbd02f6e26"

    tester = Pb.PhyloBuddy(pb_odd_resources['support'])
    tester = Pb.collapse_polytomies(tester, threshold=0.1, mode='length')
    assert pb_helpers.phylo2hash(tester) == "252572f7b9566c62df24d57065412240"

    with pytest.raises(NameError) as err:
        Pb.collapse_polytomies(tester, threshold=0.1, mode='foo')
        assert "Mode must be 'support' or 'length'" in str(err)

# ###################### 'ct', '--consensus_tree' ###################### #
hashes = [('m k', 'acd3fb34cce867c37684244701f9f5bf'), ('m n', 'eede64c804e531cb1c99e4240589b04b'),
          ('m l', '73ac98a1656d1c4a52da16d3f096f8ce'), ('o k', '64f7df66253b104c300d13e344e2f216')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_consensus_tree(key, next_hash, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    tester = Pb.consensus_tree(tester)
    assert pb_helpers.phylo2hash(tester) == next_hash, tester.write("error_files/%s" % next_hash)

hashes = [('m k', 'baf9b2def2c0fa2ff97d3a16d24b4738'), ('m n', '3ab21ed1202222f17a44fff7b4051aa0'),
          ('m l', 'd680eece99ad907bbc8cd69ceaeb7b8a'), ('o k', '64f7df66253b104c300d13e344e2f216')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_consensus_tree_95(key, next_hash, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    tester = Pb.consensus_tree(tester, frequency=0.95)
    assert pb_helpers.phylo2hash(tester) == next_hash


# ###################### 'dt', '--display_trees' ###################### #
def test_display_trees(monkeypatch, pb_resources):
    show = mock.Mock(return_value=True)
    monkeypatch.setattr(ete3.TreeNode, "show", show)
    try:
        assert Pb.display_trees(pb_resources.get_one("o k"))
    except SystemError as err:
        assert "This system is not graphical, so display_trees() will not work." in str(err)


def test_display_trees_error(pb_resources):
    # noinspection PyUnresolvedReferences
    with mock.patch.dict('os.environ'):
        if 'DISPLAY' in os.environ:
            del os.environ['DISPLAY']
        with pytest.raises(SystemError):
            Pb.display_trees(pb_resources.get_one("o k"))


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


# ######################  'gt', '--generate_trees' ###################### #
# ToDo: All of these tests need to be run on mock output. Actual 3rd party software is tested in test_alb_3rd_party.py
def test_raxml_inputs(monkeypatch, alb_resources, pb_helpers):
    tmp_dir = br.TempDir()
    if not os.path.isdir("%s/mock_resources/test_raxml_inputs" % alb_resources.res_path):
        raise NotADirectoryError("Unable to find mock resources")
    for root, dirs, files in os.walk("%s/mock_resources/test_raxml_inputs" % alb_resources.res_path):
        for _file in files:
            shutil.copyfile("%s/%s" % (root, _file), "%s/%s" % (tmp_dir.path, _file))

    tester = alb_resources.get_one("o d n")
    monkeypatch.setattr("buddysuite.buddy_resources.TempDir", lambda: tmp_dir)
    with mock.patch('buddysuite.PhyloBuddy.Popen', MockPopen):
        tester = Pb.generate_tree(tester, 'raxml', " -r_seed 12345")

    assert pb_helpers.phylo2hash(tester) == "1cede6c576bb88125e2387d850f813ab"


# PhyML version 20120412
def test_phyml_inputs(monkeypatch, alb_resources, pb_helpers):
    tmp_dir = br.TempDir()
    if not os.path.isdir("%s/mock_resources/test_phyml_inputs" % alb_resources.res_path):
        raise NotADirectoryError("Unable to find mock resources")
    for root, dirs, files in os.walk("%s/mock_resources/test_phyml_inputs" % alb_resources.res_path):
        for _file in files:
            shutil.copyfile("%s/%s" % (root, _file), "%s/%s" % (tmp_dir.path, _file))

    tester = alb_resources.get_one("o d n")
    monkeypatch.setattr("buddysuite.buddy_resources.TempDir", lambda: tmp_dir)
    with mock.patch('buddysuite.PhyloBuddy.Popen', MockPopen):
        tester = Pb.generate_tree(tester, 'phyml', " -r_seed 12345")

    assert pb_helpers.phylo2hash(tester) == "01699cec16b5dc41036d9c9f5e4894d2"


def test_fasttree_inputs(monkeypatch, alb_resources, pb_helpers):
    if not os.path.isdir("%s/mock_resources/test_phyml_inputs" % alb_resources.res_path):
        raise NotADirectoryError("Unable to find mock resources")
    tester = alb_resources.get_one("o d n")
    with open("%s/mock_resources/test_fasttree_inputs/fasttree_result.aln" % alb_resources.res_path) as treefile:
        mock_check_output = mock.Mock(return_value=treefile.read())
    monkeypatch.setattr(Pb, 'check_output', mock_check_output)
    monkeypatch.setattr(Pb, 'Popen', MockPopen)
    tester = Pb.generate_tree(tester, 'fasttree')
    assert pb_helpers.phylo2hash(tester) == "f6caf053c1cd4287ac2f80381a991d8c"


# ###################### 'hi', '--hash_ids' ###################### #
def test_hash_ids(pb_resources, pb_helpers):
    for phylobuddy in pb_resources.get_list("m o k n l"):
        orig_hash = pb_helpers.phylo2hash(phylobuddy)
        Pb.hash_ids(phylobuddy)
        assert pb_helpers.phylo2hash(phylobuddy) != orig_hash


def test_hash_ids_edges(pb_resources, pb_helpers, pb_odd_resources):
    with pytest.raises(TypeError) as e:
        Pb.hash_ids(Pb.PhyloBuddy, hash_length="foo")
    assert "Hash length argument must be an integer, not <class 'str'>" in str(e)

    with pytest.raises(ValueError) as e:
        Pb.hash_ids(Pb.PhyloBuddy, hash_length=0)
    assert "Hash length must be greater than 0" in str(e)

    with pytest.raises(ValueError) as e:
        Pb.hash_ids(pb_resources.get_one("m n"), hash_length=1)
    assert "Insufficient number of hashes available to cover all sequences." in str(e)

    tester = Pb.PhyloBuddy(pb_odd_resources['node_lables'])
    test_hash = pb_helpers.phylo2hash(tester)
    tester = Pb.hash_ids(tester, hash_length=5, nodes=True)
    assert pb_helpers.phylo2hash(tester) != test_hash

# ###################### 'li', '--list_ids' ###################### #
hashes = [('m k', '514675543e958d5177f248708405224d'), ('m n', '229e5d7cd8bb2bfc300fd45ec18e8424'),
          ('m l', '514675543e958d5177f248708405224d')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_list_ids(key, next_hash, pb_resources, pb_helpers):
    tester = str(Pb.list_ids(pb_resources.get_one(key)))
    assert pb_helpers.string2hash(tester) == next_hash

# ###################### 'ptr', '--print_trees' ###################### #
hashes = [('m k', '1f62f35b8fe64c5eaaf394c272febf4f'), ('m n', 'd5b502ae7e83f0c6ae2211b7d8281b90'),
          ('m l', '1f62f35b8fe64c5eaaf394c272febf4f'), ('o k', '8a5a1f3e409471184da4138ab45fc54e'),
          ('o n', '370b1c17c1d96ebb866277747f276add'), ('o l', '8a5a1f3e409471184da4138ab45fc54e')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_print_trees(key, next_hash, pb_resources, pb_helpers):
    tester = Pb.trees_to_ascii(pb_resources.get_one(key))
    assert pb_helpers.phylo2hash(tester) == next_hash

# ###################### 'pr', '--prune_taxa' ###################### #
pt_hashes = [('m k', '99635c6dbf708f94cf4dfdca87113c44'), ('m n', 'fc03b4f100f038277edf6a9f48913dd0'),
             ('m l', '001db76033cba463a0f187266855e8d5')]


@pytest.mark.parametrize("key, next_hash", pt_hashes)
def test_prune_taxa(key, next_hash, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    Pb.prune_taxa(tester, 'fir')
    assert pb_helpers.phylo2hash(tester) == next_hash

# ######################  'ri', '--rename_ids' ###################### #
hashes = [('m k', '6843a620b725a3a0e0940d4352f2036f'), ('m n', '543d2fc90ca1f391312d6b8fe896c59c'),
          ('m l', '6ce146e635c20ad62e21a1ed6fddbd3a'), ('o k', '4dfed97b2a23b8957ee5141bf4681fe4'),
          ('o n', '77d00fdc512fa09bd1146037d25eafa0'), ('o l', '9b1014be1b38d27f6b7ef73d17003dae')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_rename_ids(key, next_hash, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    tester = Pb.rename(tester, 'Mle', 'Phylo')
    assert pb_helpers.phylo2hash(tester) == next_hash


def test_rename_nodes(pb_odd_resources, pb_helpers):
    tester = Pb.PhyloBuddy(pb_odd_resources["node_lables"])
    Pb.rename(tester, "Equus|Ruminantiamorpha|Canis", "Family")
    assert pb_helpers.phylo2hash(tester) == "46ba6ce0f3a1859b4c7326a6f5e69263"


# ###################### 'rt', '--root' ###################### #
# Note that there's a weird behaviour with dendropy here. If running pytest with -n, it returns slightly different trees
hashes = [('m k', ['8a7fdd9421e0752c9cd58a1e073186c7', '25ea14c2e89530a0fb48163c0ef2a102']),
          ('m n', ['c16b19aaa11678595f6ed4e7c6b77955', 'e3711259d579cbf0511a5ded66dfd437']),
          ('m l', ['773e71730ccf270ea3a7cd37a7c0990d', '8135cec021240619c27d61288885d8e1']),
          ('o k', ['eacf232776eea70b5de156328e10ecc7', '05a83105f54340839dca64a62a22026e']),
          ('o n', ['53caffda3fed5b9004b79effc6d29c36', 'f0e26202274a191c9939835b25c1fae4']),
          ('o l', ['3137d568fe07d88620c08480a15006d3', 'e97c246dc7ebf4d80363f836beff4a81'])]


@pytest.mark.parametrize("key, next_hashes", hashes)
def test_root_middle(key, next_hashes, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    tester = Pb.root(tester)
    assert pb_helpers.phylo2hash(tester) in next_hashes

hashes = [('m k', 'f32bdc34bfe127bb0453a80cf7b01302'), ('m n', 'a7003478d75ad76ef61fcdc643ccdab8'),
          ('m l', 'c490e3a937b6ee2073c74119984a896e'), ('o k', 'eacf232776eea70b5de156328e10ecc7'),
          ('o n', '53caffda3fed5b9004b79effc6d29c36'), ('o l', '3137d568fe07d88620c08480a15006d3')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_root_leaf(key, next_hash, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    tester = Pb.root(tester, "firSA25a")
    assert pb_helpers.phylo2hash(tester) == next_hash

hashes = [('m k', 'edcc2400de6b0a4fb05c0a5159215ecd'), ('m n', '09e4c41d22f43b847677eec8be899a72'),
          ('m l', '89d62bd49d89daa14d2986fe8b826221'), ('o k', 'b6be77f1d16776554c5a61598ddb6899'),
          ('o n', '034f6fcb778284bd1bb9db34525f108e'), ('o l', '21d6e13646df1eb75cecb7b10e913eb0')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_root_mrca1(key, next_hash, pb_resources, pb_helpers):
    tester = pb_resources.get_one(key)
    tester = Pb.root(tester, "ovi47[ab]", "penIT11b")
    assert pb_helpers.phylo2hash(tester) == next_hash


def test_root_mrca2(pb_resources, pb_helpers):
    tester = pb_resources.get_one("o k")
    starting_hash = pb_helpers.string2hash(re.sub("\[&U\]", "[&R]", str(tester)))
    tester = Pb.root(tester, "NoMatchHere!")
    assert pb_helpers.phylo2hash(tester) == starting_hash

    # Step through a few nodes to find MRCA
    tester = pb_resources.get_one("m k")
    tester = Pb.root(tester, "ovi47[ab]", "penIT12b")
    assert pb_helpers.phylo2hash(tester) == "292f739318987701aedac62115be38c8"


# ######################  'su', '--show_unique' ###################### #
def test_show_unique(pb_odd_resources, pb_resources, pb_helpers):
    tester = Pb.PhyloBuddy(pb_odd_resources['compare'])
    Pb.show_unique(tester)
    assert pb_helpers.phylo2hash(tester) == "ea5b0d1fcd7f39cb556c0f5df96281cf"

    with pytest.raises(AssertionError):  # If only a single tree is present
        tester = Pb.PhyloBuddy(pb_resources.get_one("m k"))
        Pb.show_unique(tester)


def test_show_unique_unrooted(pb_odd_resources, pb_helpers):
    tester = Pb.PhyloBuddy(pb_odd_resources['compare'])
    Pb.unroot(tester)
    Pb.show_unique(tester)
    assert pb_helpers.phylo2hash(tester) == "2bba16e2c77102ba150adecc352407a9"


# ###################### 'sp', '--split_polytomies' ###################### #
def test_split_polytomies():
    tester = Pb.PhyloBuddy('(A,(B,C,D));')
    Pb.split_polytomies(tester)
    assert str(tester) in ['(A:1.0,(B:1.0,(C:1.0,D:1.0):1e-06):1.0):1.0;\n',
                           '(A:1.0,(B:1.0,(D:1.0,C:1.0):1e-06):1.0):1.0;\n',
                           '(A:1.0,(C:1.0,(B:1.0,D:1.0):1e-06):1.0):1.0;\n',
                           '(A:1.0,(C:1.0,(D:1.0,B:1.0):1e-06):1.0):1.0;\n',
                           '(A:1.0,(D:1.0,(B:1.0,C:1.0):1e-06):1.0):1.0;\n',
                           '(A:1.0,(D:1.0,(C:1.0,B:1.0):1e-06):1.0):1.0;\n']


# ###################### 'ur', '--unroot' ###################### #
def test_unroot(pb_odd_resources, pb_helpers):
    tester = Pb.PhyloBuddy(pb_odd_resources['figtree'])
    Pb.unroot(tester)
    assert pb_helpers.phylo2hash(tester) == "10e9024301b3178cdaed0b57ba33f615"
