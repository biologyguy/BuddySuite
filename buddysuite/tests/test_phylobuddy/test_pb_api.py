#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Tests PhyloBuddy API functions """
import pytest
import PhyloBuddy as Pb
import buddy_resources as br
from .. import __init__

from unittest import mock
import os
from os.path import join
import shutil
import re
import webbrowser
import pylab
from collections import OrderedDict
from dendropy import TreeList

RES_PATH = __init__.RESOURCE_PATH
HASH_MAP = OrderedDict([('mYSiElpW', 'Mle-Panxα9'), ('wPDsBSFF', 'Mle-Panxα7A'), ('eMBqjkZe', 'Mle-Panxα1'),
                        ('RSHOU1Si', 'Mle-Panxα3'), ('9aY8Cnau', 'Mle-Panxα12'), ('ap5o0Lez', 'Mle-Panxα11'),
                        ('SITWTd8S', 'Mle-Panxα4'), ('6pukL30B', 'Mle-Panxα8'), ('IqzLlJ2r', 'Mle-Panxα6'),
                        ('8rWWJ2h3', 'Mle-Panxα10B'), ('scPvTxj1', 'Mle-Panxα5'), ('i6SVQvFT', 'Mle-Panxα2'),
                        ('9V0yrhWi', 'Mle-Panxα10A'), ('5FHehw3fRm', 'Ate-PanxβG'), ('9J38W9UwOl', 'Ate-PanxβD'),
                        ('t7YmMoZMmc', 'Hvu-PanxβI'), ('iqwtNz79Pq', 'Che-PanxβF'), ('71B14dBaei', 'Ael-PanxβB'),
                        ('JPmIAQtH4Y', 'Ael-PanxβD'), ('ToKark9AhV', 'Hvu-PanxβL'), ('TwCggovLaz', 'Cla-PanxβD'),
                        ('QIxkMXuSco', 'Hec-PanxβD'), ('b0rsCIFOWv', 'Che-PanxβA'), ('RB5B3tt56a', 'Pph-PanxβE'),
                        ('nM2JtozKiL', 'Pph-PanxβD'), ('V6onrYnwqL', 'Nbi-PanxβF'), ('l7keL2rmHu', 'Hvu-PanxβE'),
                        ('zbHrcxNBJB', 'Ccr-PanxγA'), ('vsDsQ3V1ZY', 'Che-PanxβC'), ('PJnknRqrKD', 'Hvu-PanxβJ'),
                        ('h7zOT9CNlZ', 'Che-PanxβD'), ('a0BEtaUs4Y', 'Nbi-PanxβG'), ('FMLbM2xkqq', 'Ate-PanxβC')])


# ###################### 'cpt', '--collapse_polytomies' ###################### #
def test_add_branch(pb_odd_resources, hf):
    # Single taxon on a single sister taxon
    tester = Pb.PhyloBuddy(pb_odd_resources['lengths'])
    tester = Pb.add_branch(tester, "Foo", "Mle-Panxα8")
    assert hf.buddy2hash(tester) == "b83c4cfb07bd5df9f4bc41a7f9ebe718", print(tester)

    # Single taxon on multiple sister taxa
    tester = Pb.PhyloBuddy(pb_odd_resources['lengths'])
    tester = Pb.add_branch(tester, "Foo", ["Mle-Panxα2", "Mle-Panxα5"])
    assert hf.buddy2hash(tester) == "b74b7af2e3ebc3af1d49bfd508475c63", print(tester)

    # Subtree on a single sister taxon
    tester = Pb.PhyloBuddy(pb_odd_resources['lengths'])
    tester = Pb.add_branch(tester, Pb.PhyloBuddy("((Foo:0.78,Bar:0.34):1.1,Baz:0.2);"), "Mle-Panxα8")
    assert hf.buddy2hash(tester) == "1319f48d46800d28179655a2125fd2c3", print(tester)

    # Subtree on multiple sister taxa
    tester = Pb.PhyloBuddy(pb_odd_resources['lengths'])
    trees = TreeList()
    trees.read(data="((Foo:0.78,Bar:0.34):1.1,Baz:0.2);", schema="newick")
    tester = Pb.add_branch(tester, trees.pop(), ["Mle-Panxα2", "Mle-Panxα5"])
    assert hf.buddy2hash(tester) == "828a289b16f3a818230f54a0d4aa8c5d", print(tester)

    # Subtree on multiple sister taxa at root of tree
    tester = Pb.PhyloBuddy(pb_odd_resources['lengths'])
    trees = TreeList()
    trees.read(data="((Foo:0.78,Bar:0.34):1.1,Baz:0.2);", schema="newick")
    tester = Pb.add_branch(tester, trees.pop(), ["Mle-Panxα2", "Mle-Panxα9"])
    assert hf.buddy2hash(tester) == "4e8c323a5104154bc8be844742bca000", print(tester)

    # Error handling
    with pytest.raises(TypeError) as err:
        Pb.add_branch(tester, 123456, "Mle-Panxα8")
    assert "new_branch parameter must be a string, dentropy tree, or PhyloBuddy object." in str(err)

    with pytest.raises(AttributeError) as err:
        Pb.add_branch(tester, "Foo", "Mle-PanxαFoo")
    assert "Unable to identify any sister taxa in tree." in str(err)


# ###################### 'cpt', '--collapse_polytomies' ###################### #
def test_collapse_polytomies(pb_odd_resources, hf):
    tester = Pb.PhyloBuddy(pb_odd_resources['support'])
    tester = Pb.collapse_polytomies(tester, 20)
    assert hf.buddy2hash(tester) == "1b0979265205b17ca7f34abbd02f6e26"

    tester = Pb.PhyloBuddy(pb_odd_resources['support'])
    tester = Pb.collapse_polytomies(tester, threshold=0.1, mode='length')
    assert hf.buddy2hash(tester) == "252572f7b9566c62df24d57065412240"

    with pytest.raises(NameError) as err:
        Pb.collapse_polytomies(tester, threshold=0.1, mode='foo')
        assert "Mode must be 'support' or 'length'" in str(err)


# ###################### 'ct', '--consensus_tree' ###################### #
hashes = [('m k', 'acd3fb34cce867c37684244701f9f5bf'), ('m n', 'eede64c804e531cb1c99e4240589b04b'),
          ('m l', '73ac98a1656d1c4a52da16d3f096f8ce'), ('o k', '64f7df66253b104c300d13e344e2f216')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_consensus_tree(key, next_hash, pb_resources, hf):
    tester = pb_resources.get_one(key)
    tester = Pb.consensus_tree(tester)
    assert hf.buddy2hash(tester) == next_hash, tester.write(join("error_files", next_hash))


hashes = [('m k', 'baf9b2def2c0fa2ff97d3a16d24b4738'), ('m n', '3ab21ed1202222f17a44fff7b4051aa0'),
          ('m l', 'd680eece99ad907bbc8cd69ceaeb7b8a'), ('o k', '64f7df66253b104c300d13e344e2f216')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_consensus_tree_95(key, next_hash, pb_resources, hf):
    tester = pb_resources.get_one(key)
    tester = Pb.consensus_tree(tester, frequency=0.95)
    assert hf.buddy2hash(tester) == next_hash


# ###################### 'dt', '--display_trees' ###################### #
def test_display_trees(monkeypatch, pb_resources, capsys):
    monkeypatch.setattr(webbrowser, "open", lambda *_, **__: "")
    monkeypatch.setattr(Pb.Bio.Phylo, "draw", lambda *_, **__: True)
    monkeypatch.setattr(pylab, "axis", lambda *_: True)
    monkeypatch.setattr(pylab, "savefig", lambda *_, **__: True)
    monkeypatch.setattr("builtins.input", lambda *_: print("Mock input"))
    monkeypatch.setattr(Pb, "Popen", lambda *_, **__: True)

    with mock.patch.dict('os.environ'):
        os.environ["DISPLAY"] = ":0"
        phylobuddy = pb_resources.get_one("o k")
        tree = phylobuddy.trees[0]
        for node in tree:
            if node.taxon:
                node.annotations['!color'] = '#ff0000'
        assert Pb.display_trees(phylobuddy) is True
        out, err = capsys.readouterr()
        assert "Mock input" in out

        monkeypatch.setattr(Pb, "Popen", lambda *_, **__: True)
        monkeypatch.setattr(shutil, "which", lambda *_: True)
        assert Pb.display_trees(pb_resources.get_one("o k"), program="figtree") is True
        out, err = capsys.readouterr()
        assert "Mock input" in out

        monkeypatch.setattr(os, "name", "nt")
        assert Pb.display_trees(pb_resources.get_one("o k"), program="figtree") is False
        out, err = capsys.readouterr()
        assert "Sorry! I'm not sure how to launch FigTree on Windows" in out


def test_display_trees_pylab_setup(monkeypatch, pb_resources, capsys):
    def mock_raise_runtime1():
        raise RuntimeError("The Mac OS X backend will not be able to function correctly")

    tmp_dir = br.TempDir()
    mplibrc = os.path.join(tmp_dir.path, ".matplotlib", "matplotlibrc")

    def mock_expand_user(path):
        return re.sub("~", tmp_dir.path, path)

    monkeypatch.setattr(br, "dummy_func", mock_raise_runtime1)
    monkeypatch.setattr("builtins.input", lambda *_: print("Mock input"))
    monkeypatch.setattr(os.path, "expanduser", mock_expand_user)
    monkeypatch.setattr(br, "ask", lambda *_, **__: True)

    # First time this runs, matplotlibrc is not present in temp dir
    with mock.patch.dict('os.environ'):
        os.environ["DISPLAY"] = ":0"
        with pytest.raises(SystemExit):
            Pb.display_trees(pb_resources.get_one("o k"))

    out, err = capsys.readouterr()
    assert "%s amended.\nPlease re-run your command." % mplibrc in err
    assert os.path.isfile(mplibrc)
    with open(mplibrc, "r") as ifile:
        assert ifile.read() == "# Backend set by PhyloBuddy\nbackend: TkAgg\n"

    # Second time this runs, matplotlibrc is now present in temp dir
    with mock.patch.dict('os.environ'):
        os.environ["DISPLAY"] = ":0"
        with pytest.raises(SystemExit):
            Pb.display_trees(pb_resources.get_one("o k"))

    out, err = capsys.readouterr()
    assert "%s amended.\nPlease re-run your command." % mplibrc in err
    with open(mplibrc, "r") as ifile:
        output = ifile.read()
        assert output == ("# Backend set by PhyloBuddy\nbackend: TkAgg\n" * 2).strip(), print(output)

    # Do not agree to update matplotlibrc
    monkeypatch.setattr(br, "ask", lambda *_, **__: False)
    with open(mplibrc, "w") as ofile:
        ofile.write("Foo bar\n")

    with mock.patch.dict('os.environ'):
        os.environ["DISPLAY"] = ":0"
        with pytest.raises(SystemExit):
            Pb.display_trees(pb_resources.get_one("o k"))

    out, err = capsys.readouterr()
    assert "Stack Overflow may help you out" in err

    with open(mplibrc, "r") as ifile:
        assert ifile.read() == "Foo bar\n"

    def mock_raise_runtime2():
        raise RuntimeError("Unknown runtime error")

    monkeypatch.setattr(br, "dummy_func", mock_raise_runtime2)
    with mock.patch.dict('os.environ'):
        os.environ["DISPLAY"] = ":0"
        with pytest.raises(RuntimeError) as err:
            Pb.display_trees(pb_resources.get_one("o k"))
        assert "Unknown runtime error" in str(err)


def test_display_trees_error(pb_resources, monkeypatch):
    # Unknown program
    with pytest.raises(AttributeError) as err:
        Pb.display_trees("PhyloBuddy", program="foo")
    assert "Unknown program 'foo' selected for display. Please select between ['figtree', 'system']" in str(err)

    monkeypatch.setattr(shutil, "which", lambda *_: "")

    with mock.patch.dict('os.environ'):
        os.environ["DISPLAY"] = ":0"
        with pytest.raises(SystemError) as err:
            Pb.display_trees("PhyloBuddy", program="figtree")
        assert "Unable to find FigTree on your system. Please install it and ensure it is present in" \
               " your $PATH. http://tree.bio.ed.ac.uk/software/figtree/" in str(err)

    # noinspection PyUnresolvedReferences
    monkeypatch.setattr("builtins.input", lambda *_: "")
    monkeypatch.setattr(webbrowser, "open_new_tab", lambda *_: "")
    monkeypatch.setattr(os, "name", "posix")
    monkeypatch.setattr(Pb.sys, "platform", "blahh")

    with mock.patch.dict('os.environ'):
        if 'DISPLAY' in os.environ:
            del os.environ['DISPLAY']

        with pytest.raises(SystemError):
            Pb.display_trees(pb_resources.get_one("o k"))


# ###################### 'dis', '--distance' ###################### #
hashes = [('m k', 'aed02a04021dcbc99d41e34592ebeed0'), ('m n', 'cc23d0f2f6ae9b5a14425caef2d8ccb2'),
          ('m l', 'aed02a04021dcbc99d41e34592ebeed0')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_distance_wrf(key, next_hash, pb_resources, hf):
    tester = pb_resources.get_one(key)
    output = str(Pb.distance(tester, method='wrf'))
    assert hf.string2hash(output) == next_hash


hashes = [('m k', '49173bc33d89cc3d912a6af0fd51801d'), ('m n', '56574d4305bb5f094363a0fad351fd42'),
          ('m l', '49173bc33d89cc3d912a6af0fd51801d')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_distance_uwrf(key, next_hash, pb_resources, hf):
    tester = pb_resources.get_one(key)
    output = str(Pb.distance(tester, method='uwrf'))
    assert hf.string2hash(output) == next_hash, print(output)


hashes = [('m k', 'bae2c660250d42d6ba9bac7d311d6ffb'), ('m n', '0596ee4e2d4e76b27fe20b66a8fbea51'),
          ('m l', 'bae2c660250d42d6ba9bac7d311d6ffb')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_distance_ed(key, next_hash, pb_resources, hf):
    tester = pb_resources.get_one(key)
    output = str(Pb.distance(tester, method='ed'))
    assert hf.string2hash(output) == next_hash, print(output)


def test_distance_errors(pb_resources):
    with pytest.raises(AttributeError) as err:
        Pb.distance(pb_resources.get_one("m k"), method='foo')
    assert "foo is an invalid comparison method." in str(err)

    with pytest.raises(ValueError) as err:
        Pb.distance(pb_resources.get_one("o k"), method='ed')
    assert "Distance requires at least two trees." in str(err)


# ######################  'gt', '--generate_trees' ###################### #
class MockPopen(object):
    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs
        if " -v" in args[0]:
            self.output = ["This is RAxML version".encode("utf-8"), "".encode("utf-8")]
        else:
            self.output = None

    def communicate(self):
        return self.output

    @staticmethod
    def wait():
        return True


def mock_check_output(*args, **kwargs):
    print([args, kwargs])
    with open(join(RES_PATH, "mock_resources", "test_fasttree_inputs", "result.tre"), "r") as ifile:
        return ifile.read()


def test_raxml(alb_resources, hf, monkeypatch):
    mock_tmp_dir = br.TempDir()
    tmp_dir = br.TempDir()
    root, dirs, files = next(os.walk(join(RES_PATH, "mock_resources", "test_raxml_inputs")))
    for _file in files:
        shutil.copyfile(join(root, _file), join(mock_tmp_dir.path, _file))
    monkeypatch.setattr(Pb.shutil, "which", lambda *_: True)
    monkeypatch.setattr(Pb, "Popen", MockPopen)
    monkeypatch.setattr(br, "TempDir", lambda: mock_tmp_dir)

    # basic
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'raxml')
    assert hf.buddy2hash(tester) == "1cede6c576bb88125e2387d850f813ab", tester.write("temp.del")

    # quiet
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'raxml', quiet=True)
    assert hf.buddy2hash(tester) == "1cede6c576bb88125e2387d850f813ab", tester.write("temp.del")

    # params
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'raxml', "-w %s" % join("path", "to", "nowhere"), quiet=True)
    assert hf.buddy2hash(tester) == "1cede6c576bb88125e2387d850f813ab", tester.write("temp.del")

    # slight edges
    tester = alb_resources.get_one("o p n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'raxml-HPC', tmp_dir.path, quiet=True)
    assert hf.buddy2hash(tester) == "1cede6c576bb88125e2387d850f813ab", tester.write("temp.del")

    # bootstraps
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'raxml', "-b 1234 -N 4", quiet=True)
    assert hf.buddy2hash(tester) == "b8a3b6068aa06bd8a70fcc9bda0efad9", tester.write("temp.del")

    # keep
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    Pb.generate_tree(tester, 'raxml', keep_temp=join(tmp_dir.path, "keep_files"))
    root, dirs, files = next(os.walk(join(tmp_dir.path, "keep_files")))
    kept_output = ""
    for file in sorted(files):
        with open(join(root, file), "r") as ifile:
            kept_output += ifile.read()
    if os.name == "nt":
        assert hf.string2hash(kept_output) == "7f2fdfe55dbe805bd994f3f56c79bb1b"
    else:
        assert hf.string2hash(kept_output) == "393353fb47861460aecaefa69a6ec55c"

    # multi-run
    os.remove(join(mock_tmp_dir.path, "RAxML_bestTree.result"))
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'raxml', "-N 3", quiet=True)
    assert hf.buddy2hash(tester) == "e8ce19bba744f0188df2ebb4ffced4d8", tester.write("temp.del")

    # bipartitions
    shutil.copy(join(RES_PATH, "mock_resources", "test_raxml_inputs", "bipartitions", "RAxML_bipartitions.result"),
                join(mock_tmp_dir.path, "RAxML_bipartitions.result"))
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'raxml', "-f b -z {0}{1}RAxML_bootstrap.result "
                                               "-t {0}{1}RAxML_bestTree.result".format(tmp_dir.path, os.path.sep))
    assert hf.buddy2hash(tester) == "457533ada8e987fd0c50a41aabe1700b"


def test_phyml(alb_resources, hf, monkeypatch):
    mock_tmp_dir = br.TempDir()
    tmp_dir = br.TempDir()
    for root, dirs, files in os.walk(join(RES_PATH, "mock_resources", "test_phyml_inputs")):
        for _file in files:
            shutil.copyfile(join(root, _file), join(mock_tmp_dir.path, _file))
    monkeypatch.setattr(Pb.shutil, "which", lambda *_: True)
    monkeypatch.setattr(Pb, "Popen", MockPopen)
    monkeypatch.setattr(br, "TempDir", lambda: mock_tmp_dir)

    # basic
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'phyml')
    assert hf.buddy2hash(tester) == "9fbcfe1d565b9fd23e2c5fca86019f8e", tester.write("temp.del")

    # quiet
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'phyml', quiet=True)
    assert hf.buddy2hash(tester) == "9fbcfe1d565b9fd23e2c5fca86019f8e", tester.write("temp.del")

    # params
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'phyml', "--sequential", quiet=True)
    assert hf.buddy2hash(tester) == "9fbcfe1d565b9fd23e2c5fca86019f8e", tester.write("temp.del")

    # slight edges
    tester = alb_resources.get_one("o p n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'phyml-foo', tmp_dir.path, quiet=True)
    assert hf.buddy2hash(tester) == "9fbcfe1d565b9fd23e2c5fca86019f8e", tester.write("temp.del")

    # keep
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    Pb.generate_tree(tester, 'phyml', keep_temp=join(tmp_dir.path, "keep_files"))
    root, dirs, files = next(os.walk(join(tmp_dir.path, "keep_files")))
    kept_output = ""
    for _file in sorted(files):
        with open(join(root, _file), "r") as ifile:
            kept_output += ifile.read()
    if os.name == "nt":
        assert hf.string2hash(kept_output) == "a795da6869c5e3a34962a52ec35006ed"
    else:
        assert hf.string2hash(kept_output) == "5a3559c264cb4c4779f15a515aaf2286"


def test_fasttree(alb_resources, hf, monkeypatch):
    mock_tmp_dir = br.TempDir()
    tmp_dir = br.TempDir()
    for root, dirs, files in os.walk(join(RES_PATH, "mock_resources", "test_fasttree_inputs")):
        for _file in files:
            shutil.copyfile(join(root, _file), join(mock_tmp_dir.path, _file))
    monkeypatch.setattr(Pb.shutil, "which", lambda *_: True)
    monkeypatch.setattr(Pb, "Popen", MockPopen)
    monkeypatch.setattr(Pb, "check_output", mock_check_output)
    monkeypatch.setattr(br, "TempDir", lambda: mock_tmp_dir)

    # basic
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'fasttree')
    assert hf.buddy2hash(tester) == "4b2ab8c39f27b9871f9a370ca7d0c4b3", tester.write("temp.del")

    # quiet
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'fasttree', quiet=True)
    assert hf.buddy2hash(tester) == "4b2ab8c39f27b9871f9a370ca7d0c4b3", tester.write("temp.del")

    # params
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'fasttree', "-pseudo", quiet=True)
    assert hf.buddy2hash(tester) == "4b2ab8c39f27b9871f9a370ca7d0c4b3", tester.write("temp.del")

    # slight edges
    tester = alb_resources.get_one("m p s")
    tester.hash_map = HASH_MAP
    tester = Pb.generate_tree(tester, 'fasttree-foo', quiet=True)
    assert hf.buddy2hash(tester) == "7692ee7edb8df4f91d2bdca6c3767796", tester.write("temp.del")

    # keep
    tester = alb_resources.get_one("o d n")
    tester.hash_map = HASH_MAP
    Pb.generate_tree(tester, 'fasttree', keep_temp=join(tmp_dir.path, "keep_files"))
    root, dirs, files = next(os.walk(join(tmp_dir.path, "keep_files")))
    kept_output = ""
    for _file in sorted(files):
        with open(join(root, _file), "r") as ifile:
            kept_output += ifile.read()
    if os.name == "nt":
        assert hf.string2hash(kept_output) == "bcd034f0db63a7b41f4b3b6661200ef3"
    else:
        assert hf.string2hash(kept_output) == "10df99cd1696b002f841aa18b78477ca"


def test_generate_tree_edges(alb_resources, monkeypatch):
    mock_tmp_dir = br.TempDir()
    tmp_dir = br.TempDir()
    for root, dirs, files in os.walk(join(RES_PATH, "mock_resources", "test_fasttree_inputs")):
        for _file in files:
            shutil.copyfile(join(root, _file), join(mock_tmp_dir.path, _file))
    monkeypatch.setattr(Pb.shutil, "which", lambda *_: True)
    monkeypatch.setattr(Pb, "Popen", MockPopen)
    monkeypatch.setattr(br, "TempDir", lambda: mock_tmp_dir)

    with pytest.raises(RuntimeError) as err:
        Pb.generate_tree(alb_resources.get_one("o d n"), 'fasttree', "--nonsense")
    assert "fasttree threw an error. Scroll up for more info." in str(err)

    with pytest.raises(FileExistsError) as err:
        Pb.generate_tree(alb_resources.get_one("o d n"), 'raxml', keep_temp=tmp_dir.path)
    assert "Execution of raxml was halted to prevent files in" in str(err)

    with pytest.raises(FileNotFoundError) as err:
        Pb.generate_tree(alb_resources.get_one("o d n"), 'Foo')
    assert 'Foo failed to generate a tree' in str(err)

    monkeypatch.setattr(Pb.shutil, "which", lambda *_: None)
    with pytest.raises(ProcessLookupError) as err:
        Pb.generate_tree(alb_resources.get_one("o d n"), 'Foo')
    assert "ProcessLookupError: #### Could not find Foo on your system. ####" in str(err)

    monkeypatch.setattr(Pb.Popen, "communicate", lambda *_: ["".encode(), "".encode()])
    with pytest.raises(AttributeError) as err:
        Pb.generate_tree(alb_resources.get_one("o d n"), 'Foo')
    assert "Foo is not a recognized tree inference tool. Please check your spelling (case sensitive)" in str(err)


# ###################### 'hi', '--hash_ids' ###################### #
def test_hash_ids(pb_resources, hf):
    for phylobuddy in pb_resources.get_list("m o k n l"):
        orig_hash = hf.buddy2hash(phylobuddy)
        Pb.hash_ids(phylobuddy)
        assert hf.buddy2hash(phylobuddy) != orig_hash


def test_hash_ids_edges(monkeypatch, pb_resources, hf, pb_odd_resources):
    class MockRandom(object):
        def __init__(self):
            self.values = ["YCva8wcKxN", "r5iHz5XItJ", "b0dfUBWelI", "Tv0XxZziIQ", "hkwU50UQfq", "Jbn4JWgxXy",
                           "OeGgev9Mdu", "pQRpDbjomQ", "z9JvFAa1vZ", "vHwNvExi6A", "i6tluyP88x", "mXvIO8GYAg",
                           "clDGuGrku8", "atmvoy1uN6", "0R59tE6LgP", "eOUHJPbVaf", "2U0FarU2Lv", "4tITUmacIV",
                           "UND58iKTiH", "aqQWMJnjeu", "rOFRXh4ftJ", "HeAcnkJ8LV", "bdWtyB4OP4", "L0aFqF3oMc",
                           "4HN9dP8uqL", "TtdeGIxxY6", "s4S1RfO9sO", "mIScdjmMzq", "ohGB7fpJGq", "CtLTzCg6RB",
                           "QufRP2p8jN", "NA98oE9g9u", "6IbtWWrdDr", "ZqnxsFXwpT", "hZd46Bc15s", "iL936ARfV9",
                           "iOcjJIWv7Y", "yVkDeDhFpK", "d2amahGHzt", "H9VzjZZdbW", "kxFL0xMfFx", "glhRhmReup",
                           "kxFL0xMfFx"]
            self.current_value = self.values.pop()

        def choice(self, *_):
            if not self.current_value:
                self.current_value = self.values.pop()
            next_char = self.current_value[0]
            self.current_value = self.current_value[1:]
            return next_char

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
    test_hash = hf.buddy2hash(tester)
    tester = Pb.hash_ids(tester, hash_length=5, nodes=True)
    assert hf.buddy2hash(tester) != test_hash

    monkeypatch.setattr(Pb.random, "Random", MockRandom)
    tester = Pb.hash_ids(pb_resources.get_one("o n"))
    assert hf.buddy2hash(tester) == "48b1b2b0e1f7012ea1a964269300ac6b"


# ###################### 'ld', '--ladderize' ###################### #
def test_ladderize(pb_resources, hf):
    tester = pb_resources.get_one("m k")
    tester = Pb.ladderize(tester)
    assert hf.buddy2hash(tester) == "63ee71da75031d09f953932a1f0195b5"

    tester = Pb.ladderize(tester, ascending=False)
    assert hf.buddy2hash(tester) == "0dfa9fbb23428d2992b982776777428c"


# ###################### 'li', '--list_ids' ###################### #
hashes = [('m k', '514675543e958d5177f248708405224d'), ('m n', '229e5d7cd8bb2bfc300fd45ec18e8424'),
          ('m l', '514675543e958d5177f248708405224d')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_list_ids(key, next_hash, pb_resources, hf):
    tester = str(Pb.list_ids(pb_resources.get_one(key)))
    assert hf.string2hash(tester) == next_hash


# ###################### 'ptr', '--print_trees' ###################### #
hashes = [('m k', ['16d1fa2a370fc41160bf06532e6f0a04', '1b276812fec15fc5f3ec21e680473994']),
          ('m n', ['16d1fa2a370fc41160bf06532e6f0a04', '1b276812fec15fc5f3ec21e680473994']),
          ('m l', ['16d1fa2a370fc41160bf06532e6f0a04', '1b276812fec15fc5f3ec21e680473994']),
          ('o k', ['c684c870ef1e11ee311942232d8b2a3b', 'd78151efc2909afa92d6521dad122624']),
          ('o n', ['c684c870ef1e11ee311942232d8b2a3b', 'd78151efc2909afa92d6521dad122624']),
          ('o l', ['c684c870ef1e11ee311942232d8b2a3b', 'd78151efc2909afa92d6521dad122624'])]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_print_trees(key, next_hash, pb_resources, hf):
    tester = Pb.trees_to_ascii(pb_resources.get_one(key))
    tester = "\n".join([tree for indx, tree in tester.items()])
    if os.name == "nt":
        assert hf.string2hash(tester) == next_hash[1]
    else:
        assert hf.string2hash(tester) == next_hash[0]


# ###################### 'pr', '--prune_taxa' ###################### #
pt_hashes = [('m k', '99635c6dbf708f94cf4dfdca87113c44'), ('m n', 'fc03b4f100f038277edf6a9f48913dd0'),
             ('m l', '001db76033cba463a0f187266855e8d5')]


@pytest.mark.parametrize("key, next_hash", pt_hashes)
def test_prune_taxa(key, next_hash, pb_resources, hf):
    tester = pb_resources.get_one(key)
    Pb.prune_taxa(tester, 'fir')
    assert hf.buddy2hash(tester) == next_hash


def test_prune_taxa_edge(pb_resources, hf):
    # Test for when all leaves are deleted with prune
    fir_only = pb_resources.get_one("o n")
    Pb.rename(fir_only, "pen|ovi", "fir")
    mix = pb_resources.get_one("o n")
    mix.trees += fir_only.trees
    Pb.prune_taxa(mix, "fir")
    assert hf.buddy2hash(mix) == "6dd711f4330e7f088563045616085ba5"


# ######################  'ri', '--rename_ids' ###################### #
hashes = [('m k', '6843a620b725a3a0e0940d4352f2036f'), ('m n', '543d2fc90ca1f391312d6b8fe896c59c'),
          ('m l', '6ce146e635c20ad62e21a1ed6fddbd3a'), ('o k', '4dfed97b2a23b8957ee5141bf4681fe4'),
          ('o n', '77d00fdc512fa09bd1146037d25eafa0'), ('o l', '9b1014be1b38d27f6b7ef73d17003dae')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_rename_ids(key, next_hash, pb_resources, hf):
    tester = pb_resources.get_one(key)
    tester = Pb.rename(tester, 'Mle', 'Phylo')
    assert hf.buddy2hash(tester) == next_hash


def test_rename_nodes(pb_odd_resources, hf):
    tester = Pb.PhyloBuddy(pb_odd_resources["node_lables"])
    Pb.rename(tester, "Equus|Ruminantiamorpha|Canis", "Family")
    assert hf.buddy2hash(tester) == "46ba6ce0f3a1859b4c7326a6f5e69263"


# ###################### 'rt', '--root' ###################### #
# Note that there's a weird behaviour with dendropy here. If running pytest with -n, it returns slightly different trees
hashes = [('m k', ['8a7fdd9421e0752c9cd58a1e073186c7', '25ea14c2e89530a0fb48163c0ef2a102']),
          ('m n', ['c16b19aaa11678595f6ed4e7c6b77955', 'e3711259d579cbf0511a5ded66dfd437']),
          ('m l', ['773e71730ccf270ea3a7cd37a7c0990d', '8135cec021240619c27d61288885d8e1']),
          ('o k', ['eacf232776eea70b5de156328e10ecc7', '05a83105f54340839dca64a62a22026e']),
          ('o n', ['53caffda3fed5b9004b79effc6d29c36', 'f0e26202274a191c9939835b25c1fae4']),
          ('o l', ['3137d568fe07d88620c08480a15006d3', 'e97c246dc7ebf4d80363f836beff4a81'])]


@pytest.mark.parametrize("key, next_hashes", hashes)
def test_root_middle(key, next_hashes, pb_resources, hf):
    tester = pb_resources.get_one(key)
    tester = Pb.root(tester)
    assert hf.buddy2hash(tester) in next_hashes


hashes = [('m k', 'f32bdc34bfe127bb0453a80cf7b01302'), ('m n', 'a7003478d75ad76ef61fcdc643ccdab8'),
          ('m l', 'c490e3a937b6ee2073c74119984a896e'), ('o k', 'eacf232776eea70b5de156328e10ecc7'),
          ('o n', '53caffda3fed5b9004b79effc6d29c36'), ('o l', '3137d568fe07d88620c08480a15006d3')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_root_leaf(key, next_hash, pb_resources, hf):
    tester = pb_resources.get_one(key)
    tester = Pb.root(tester, "firSA25a")
    assert hf.buddy2hash(tester) == next_hash


hashes = [('m k', 'edcc2400de6b0a4fb05c0a5159215ecd'), ('m n', '09e4c41d22f43b847677eec8be899a72'),
          ('m l', '89d62bd49d89daa14d2986fe8b826221'), ('o k', 'b6be77f1d16776554c5a61598ddb6899'),
          ('o n', '034f6fcb778284bd1bb9db34525f108e'), ('o l', '21d6e13646df1eb75cecb7b10e913eb0')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_root_mrca1(key, next_hash, pb_resources, hf):
    tester = pb_resources.get_one(key)
    tester = Pb.root(tester, "ovi47[ab]", "penIT11b")
    assert hf.buddy2hash(tester) == next_hash


def test_root_mrca2(pb_resources, hf):
    tester = pb_resources.get_one("o k")
    starting_hash = hf.string2hash(re.sub("\[&U\]", "[&R]", str(tester)))
    tester = Pb.root(tester, "NoMatchHere!")
    assert hf.buddy2hash(tester) == starting_hash

    # Step through a few nodes to find MRCA
    tester = pb_resources.get_one("m k")
    tester = Pb.root(tester, "ovi47[ab]", "penIT12b")
    assert hf.buddy2hash(tester) == "292f739318987701aedac62115be38c8"


# ######################  'su', '--show_unique' ###################### #
def test_show_unique(pb_odd_resources, pb_resources, hf):
    tester = Pb.PhyloBuddy(pb_odd_resources['compare'])
    Pb.show_unique(tester)
    assert hf.buddy2hash(tester) == "396e27e3c7c5aa126ec07f31307a288e"

    with pytest.raises(AssertionError):  # If only a single tree is present
        tester = Pb.PhyloBuddy(pb_resources.get_one("m k"))
        Pb.show_unique(tester)


def test_show_unique_unrooted(pb_odd_resources, hf):
    tester = Pb.PhyloBuddy(pb_odd_resources['compare'])
    Pb.unroot(tester)
    Pb.show_unique(tester)
    assert hf.buddy2hash(tester) == "50f1d86072989ac61dcf95ee5fc19f7e"


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
def test_unroot(pb_odd_resources, hf):
    tester = Pb.PhyloBuddy(pb_odd_resources['figtree'])
    Pb.unroot(tester)
    assert hf.buddy2hash(tester) == "e24e85fdc2f877f9340f147d9fed5fef", print(tester)
