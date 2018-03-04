#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: test_pb_ui.py
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Collection of PyTest unit tests for the PhyloBuddy.py program
"""


import pytest
import os
import sys
import argparse
from copy import deepcopy
from unittest import mock
import shutil
import webbrowser
import pylab

import buddy_resources as br
import PhyloBuddy as Pb

VERSION = Pb.VERSION
WRITE_FILE = br.TempFile()


def fmt(prog):
    return br.CustomHelpFormatter(prog)


parser = argparse.ArgumentParser(prog="PhyloBuddy.py", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                 description='''\
\033[1mPhyloBuddy\033[m
Put a little bonsai into your phylogeny.

\033[1mUsage examples\033[m:
PhyloBuddy.py "/path/to/tree_file" -<cmd>
PhyloBuddy.py "/path/to/tree_file" -<cmd> | PhyloBuddy.py -<cmd>
PhyloBuddy.py "(A,(B,C));" -f "raw" -<cmd>
''')

br.flags(parser, ("trees", "Supply file path(s) or raw tree string, If piping trees into PhyloBuddy "
                           "this argument can be left blank."), br.pb_flags, br.pb_modifiers, VERSION)

# This is to allow py.test to work with its own flags
in_args = parser.parse_args([])


def mock_fileexistserror(*args, **kwargs):
    raise FileExistsError(args, kwargs)


def mock_filenotfounderror(*args, **kwargs):
    raise FileNotFoundError(args, kwargs)


def mock_keyboardinterrupt(*args, **kwargs):
    raise KeyboardInterrupt(args, kwargs)


def mock_guesserror(*args, **kwargs):
    raise br.GuessError("%s, %s" % (args, kwargs))


def mock_assertionerror(*args, **kwargs):
    raise AssertionError("%s, %s" % (args, kwargs))


def mock_systemexit(*args, **kwargs):
    sys.exit("%s, %s" % (args, kwargs))


# ###################### argparse_init() ###################### #
def test_argparse_init(capsys, monkeypatch, pb_odd_resources, hf):
    monkeypatch.setattr(sys, "argv", ['PhyloBuddy.py', pb_odd_resources["compare"]])
    temp_in_args, phylobuddy = Pb.argparse_init()
    assert hf.string2hash(str(phylobuddy)) == "d8e14a2bfc8e9c0ac3c524f5fb478c67"

    monkeypatch.setattr(sys, "argv", ['PhyloBuddy.py', pb_odd_resources["compare"], "-f", "foo"])
    with pytest.raises(SystemExit):
        Pb.argparse_init()

    out, err = capsys.readouterr()
    assert "Error: The format 'foo' passed in with the -f flag is not recognized." in err

    monkeypatch.setattr(sys, "argv", ['PhyloBuddy.py', pb_odd_resources["compare"], "-o", "foo"])
    with pytest.raises(SystemExit):
        Pb.argparse_init()
    out, err = capsys.readouterr()
    assert "Error: Output type foo is not recognized/supported" in err


# ###################### INTERNAL FUNCTIONS ###################### #
def test_print_trees_internal_ui(capsys, pb_resources):
    # Most of this function is covered repeatedly below, so just test the -t flag
    test_in_args = deepcopy(in_args)
    test_in_args.test = True
    test_in_args.unroot = True
    Pb.command_line_ui(test_in_args, pb_resources.get_one('m k'), skip_exit=True)
    out, err = capsys.readouterr()
    assert err == "*** Test passed ***\n"


def test_in_place_ui(capsys, pb_resources):
    # Some of this function is covered below, so just test an edge
    test_in_args = deepcopy(in_args)
    test_in_args.in_place = True
    test_in_args.trees = [pb_resources.get_one("m k")]
    test_in_args.screw_formats = "nexus"
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Warning: The -i flag was passed in, but the positional" in err


# ###################### 'ab', '--add_branch' ###################### #
def test_add_branch_ui(capsys, pb_odd_resources, hf):
    test_in_args = deepcopy(in_args)

    test_in_args.add_branch = [["Foo", "α8"]]
    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(pb_odd_resources['lengths']), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "b83c4cfb07bd5df9f4bc41a7f9ebe718"

    test_in_args.add_branch = [["((Foo:0.78,Bar:0.34):1.1,Baz:0.2);", "α2", "α5"]]
    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(pb_odd_resources['lengths']), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "828a289b16f3a818230f54a0d4aa8c5d"

    test_in_args.add_branch = [["Foo"]]
    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(pb_odd_resources['lengths']), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Add branch tool requires at least two arguments: new_branch and sister_taxa." in err


# ###################### 'cpt', '--collapse_polytomies' ###################### #
def test_collapse_polytomies_ui(capsys, pb_odd_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.collapse_polytomies = [[20]]
    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(pb_odd_resources['support']), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "1b0979265205b17ca7f34abbd02f6e26"

    test_in_args.collapse_polytomies = [[0.1, 'length']]
    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(pb_odd_resources['support']), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "252572f7b9566c62df24d57065412240"


# ###################### 'ct', '--consensus_tree' ###################### #
def test_consensus_tree_ui(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.consensus_tree = [False]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "acd3fb34cce867c37684244701f9f5bf"

    test_in_args.consensus_tree = [0.9]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "baf9b2def2c0fa2ff97d3a16d24b4738"

    test_in_args.consensus_tree = [1.5]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "acd3fb34cce867c37684244701f9f5bf"
    assert err == "Warning: The frequency value should be between 0 and 1. Defaulting to 0.5.\n\n"


# ###################### 'dt', '--display_trees' ###################### #
def test_display_trees_ui(monkeypatch, pb_resources, capsys):
    with mock.patch.dict('os.environ'):
        os.environ["DISPLAY"] = ":0"
        test_in_args = deepcopy(in_args)
        test_in_args.display_trees = [None]
        monkeypatch.setattr("builtins.input", lambda *_: "")
        monkeypatch.setattr(webbrowser, "open_new_tab", lambda *_: "")
        monkeypatch.setattr(pylab, "axis", lambda *_: True)
        monkeypatch.setattr(pylab, "savefig", lambda *_, **__: True)
        monkeypatch.setattr(Pb.Bio.Phylo, "draw", lambda *_, **__: True)
        Pb.command_line_ui(test_in_args, pb_resources.get_one("o k"), skip_exit=True)

        test_in_args.display_trees = ["fake_program"]
        Pb.command_line_ui(test_in_args, pb_resources.get_one("o k"), skip_exit=True)
        out, err = capsys.readouterr()
        assert "AttributeError: Unknown program 'fake_program' selected for display" in err


def test_display_trees_ui_no_display(capsys, monkeypatch, pb_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.display_trees = [None]
    monkeypatch.setattr("builtins.input", lambda *_: "")
    monkeypatch.setattr(webbrowser, "open_new_tab", lambda *_: "")
    monkeypatch.setattr(os, "name", "posix")

    # Non-graphical, not Mac
    monkeypatch.setattr(sys, "platform", "blahh")
    with mock.patch.dict('os.environ'):
        if 'DISPLAY' in os.environ:
            del os.environ['DISPLAY']
        Pb.command_line_ui(test_in_args, pb_resources.get_one("o k"), skip_exit=True)
        out, err = capsys.readouterr()
    assert "Error: Your system does not appear to be graphical" in err

    # Try XQuarts Mac
    def mock_systemerror1(*_, **__):
        raise SystemError("graphical fail")
    monkeypatch.setattr(sys, "platform", "darwin")
    monkeypatch.setattr(br, "dummy_func", mock_systemerror1)
    with mock.patch.dict('os.environ'):
        if 'DISPLAY' in os.environ:
            del os.environ['DISPLAY']
        Pb.command_line_ui(test_in_args, pb_resources.get_one("o k"), skip_exit=True)
        out, err = capsys.readouterr()
    assert "Try installing XQuartz (https://www.xquartz.org/)" in err, print(out, err)

    # Unknown system error
    def mock_systemerror2(*_, **__):
        raise SystemError("Foo")
    monkeypatch.setattr(br, "dummy_func", mock_systemerror2)
    with mock.patch.dict('os.environ'):
        os.environ["DISPLAY"] = ":0"
        with pytest.raises(SystemError) as err:
            Pb.command_line_ui(test_in_args, pb_resources.get_one("o k"), skip_exit=True)
    assert "Foo" in str(err)


# ###################### 'dis', '--distance' ###################### #
def test_distance_ui(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.distance = [False]  # Should default to wrf
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "66df36fa117e4c4660f04e2649c3fa6b"
    assert err == "Tree 1	Tree 2	Value\n"

    test_in_args.distance = ["wrf"]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "66df36fa117e4c4660f04e2649c3fa6b"
    assert err == "Tree 1	Tree 2	Value\n"

    test_in_args.distance = ["uwrf"]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "9dd162b4cc4fa28f402e6f31ef2fb349"
    assert err == "Tree 1	Tree 2	Value\n"

    test_in_args.distance = ["ed"]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "a49e54a6c14ab4dbbf73d3f7d1d6aa82"
    assert err == "Tree 1	Tree 2	Value\n"

    test_in_args.distance = ["foo"]
    with pytest.raises(AttributeError) as err:
        Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), pass_through=True)
    assert "foo is an invalid comparison method." in str(err)

    with pytest.raises(ValueError) as err:
        Pb.command_line_ui(test_in_args, pb_resources.get_one("o k"), pass_through=True)
    assert "Distance requires at least two trees." in str(err)


# ###################### 'gt', '--generate_tree' ###################### #
def test_generate_tree_ui1(capsys, pb_resources, alb_resources, monkeypatch):
    test_in_args = deepcopy(in_args)
    test_in_args.in_format, test_in_args.out_format = "nexus", "newick"
    test_in_args.trees = [alb_resources.get_one("o d n")]

    monkeypatch.setattr(Pb, "generate_tree", lambda *_, **__: pb_resources.get_one("o n"))
    monkeypatch.setattr(shutil, "which", lambda *_: True)

    test_in_args.generate_tree = [None]
    Pb.command_line_ui(test_in_args, [], skip_exit=True)
    out, err = capsys.readouterr()
    assert out == str(pb_resources.get_one("o k"))

    monkeypatch.setattr(Pb, "generate_tree", mock_fileexistserror)
    with pytest.raises(FileExistsError):
        Pb.command_line_ui(test_in_args, [], pass_through=True)

    monkeypatch.setattr(Pb, "generate_tree", mock_filenotfounderror)
    with pytest.raises(FileNotFoundError):
        Pb.command_line_ui(test_in_args, [], pass_through=True)

    monkeypatch.setattr(shutil, "which", lambda _: False)
    with pytest.raises(AttributeError) as err:
        Pb.command_line_ui(test_in_args, [], pass_through=True)
    assert "Unable to identify any supported phylogenetic inference software on your system." in str(err)


# ###################### 'hi', '--hash_ids' ###################### #
def test_hash_ids_ui(capsys, monkeypatch, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.hash_ids = [[1, "nodes"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("o n"), skip_exit=True)
    out, err = capsys.readouterr()

    assert hf.string2hash(out) != hf.buddy2hash(pb_resources.get_one("o n"))
    assert "Warning: The hash_length parameter was passed in with the value 1" in err

    test_in_args.hash_ids = [[-1, "nodes"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m n"), skip_exit=True)
    out, err = capsys.readouterr()

    assert hf.string2hash(out) != hf.buddy2hash(pb_resources.get_one("m n"))
    assert "Warning: The hash_length parameter was passed in with the value -1" in err

    def hash_ids(*args, **kwargs):
        if args or kwargs:
            pass
        raise ValueError("Foo bar")

    test_in_args.hash_ids = [[1, "nodes"]]
    monkeypatch.setattr(Pb, "hash_ids", hash_ids)
    with pytest.raises(ValueError):
        Pb.command_line_ui(test_in_args, pb_resources.get_one("o n"), pass_through=True)


# ###################### 'ld', '--ladderize' ###################### #
def test_ladderize_ui(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.ladderize = [None]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "63ee71da75031d09f953932a1f0195b5"

    test_in_args.ladderize = ['reverse']
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "0dfa9fbb23428d2992b982776777428c"

    test_in_args.ladderize = ['rev']
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "0dfa9fbb23428d2992b982776777428c"


# ###################### 'li', '--list_ids' ###################### #
def test_list_ids_ui(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.list_ids = [True]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "4f0857e0211ff0cd058d0cf7cbaf64d5"

    test_in_args.list_ids = [4]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "bf2fbfe1bd52e9b27ae21f5c06e7763a"

    test_in_args.list_ids = [4]
    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy("(, );"), skip_exit=True)
    out, err = capsys.readouterr()
    assert out == "#### tree_1 ####\nNone\n\n"


# ###################### '-nt', '--num_tips' ###################### #
def test_num_tips_ui(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.num_tips = True

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "f43920f4df66e76fbacae4af178eeebb"


# ###################### 'ptr', '--print_trees' ###################### #
def test_print_trees_ui(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.print_trees = True

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    if os.name == "nt":
        assert hf.string2hash(out) == "05d286a56bf98457c17830bb9c1766d0"
    else:
        assert hf.string2hash(out) == "fe340117cb8f573100c00fc897e6c8ce"


# ###################### 'pt', '--prune_taxa' ###################### #
def test_prune_taxa_ui(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.prune_taxa = [["fir"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "99635c6dbf708f94cf4dfdca87113c44"

    test_in_args.prune_taxa = [["fir", "ovi"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m n"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "2a385fa95024323fea412fd2b3c3e91f"


# ###################### 'ri', '--rename_ids' ###################### #
def test_rename_ids_ui(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.rename_ids = ['Mle', 'Phylo']

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "6843a620b725a3a0e0940d4352f2036f"


# ###################### 'rt', '--root' ###################### #
def test_root_ui_midpoint(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.root = [[]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) in ["8a7fdd9421e0752c9cd58a1e073186c7", "25ea14c2e89530a0fb48163c0ef2a102"]


def test_root_ui_leaf(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.root = [["firSA25a"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "f32bdc34bfe127bb0453a80cf7b01302"


def test_root_ui_mrca(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.root = [["firSA25a", "penSH31b"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "757489907bd5c82882d084ffcb22cfba"


# ###################### 'sf', '--screw_formats' ###################### #
def test_screw_formats_ui(capsys, pb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = "nexus"

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "543d2fc90ca1f391312d6b8fe896c59c"


def test_screw_formats_fail(capsys, pb_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = "foo"
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    foo_out, foo_err = capsys.readouterr()
    assert foo_err == "Error: unknown format 'foo'\n"


def test_screw_formats_inplace_ui(capsys, pb_odd_resources):
    temp_file = br.TempFile()
    with open(pb_odd_resources["compare"], "r") as ifile:
        temp_file.write(ifile.read())

    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = "nexus"
    test_in_args.in_place = True
    test_in_args.trees[0] = temp_file.path

    tester = Pb.PhyloBuddy(temp_file.path)
    Pb.command_line_ui(test_in_args, tester, skip_exit=True)
    out, err = capsys.readouterr()
    assert "File overwritten at:" in err
    check_file = os.path.isfile("%s.nex" % temp_file.path)
    assert check_file

    test_in_args.trees[0] = "%s.nex" % temp_file.path
    test_in_args.screw_formats = "newick"

    tester = Pb.PhyloBuddy("%s.nex" % temp_file.path)
    Pb.command_line_ui(test_in_args, tester, skip_exit=True)
    out, err = capsys.readouterr()
    assert "File overwritten at:" in err
    check_file = os.path.isfile("%s.nwk" % temp_file.path)
    assert check_file


# ###################### 'su', '--show_unique' ###################### #
def test_show_unique_ui(capsys, pb_resources, hf, pb_odd_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.show_unique = True

    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(pb_odd_resources["compare"]), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "396e27e3c7c5aa126ec07f31307a288e"

    with pytest.raises(AssertionError) as err:
        Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), pass_through=True)
    assert "PhyloBuddy object should have exactly 2 trees." in str(err)


# ###################### 'sp', '--split_polytomies' ###################### #
def test_split_polytomies_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.split_polytomies = True

    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(Pb.PhyloBuddy('(A,(B,C,D));')), skip_exit=True)
    out, err = capsys.readouterr()
    assert out in ['(A:1.0,(B:1.0,(C:1.0,D:1.0):1e-06):1.0):1.0;\n',
                   '(A:1.0,(B:1.0,(D:1.0,C:1.0):1e-06):1.0):1.0;\n',
                   '(A:1.0,(C:1.0,(B:1.0,D:1.0):1e-06):1.0):1.0;\n',
                   '(A:1.0,(C:1.0,(D:1.0,B:1.0):1e-06):1.0):1.0;\n',
                   '(A:1.0,(D:1.0,(B:1.0,C:1.0):1e-06):1.0):1.0;\n',
                   '(A:1.0,(D:1.0,(C:1.0,B:1.0):1e-06):1.0):1.0;\n']


# ###################### 'ur', '--unroot' ###################### #
def test_unroot_ui(capsys, pb_odd_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.unroot = True

    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(pb_odd_resources["figtree"]), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "e24e85fdc2f877f9340f147d9fed5fef", print(out)


# ###################### main() ###################### #
def test_main(monkeypatch, pb_resources):
    monkeypatch.setattr(sys, "argv", ["PhyloBuddy", pb_resources.get_one("o k", mode="paths"), "-li"])
    monkeypatch.setattr(Pb, "command_line_ui", lambda *_: True)
    assert Pb.main()

    monkeypatch.setattr(Pb, "command_line_ui", mock_keyboardinterrupt)
    assert not Pb.main()

    monkeypatch.setattr(Pb, "command_line_ui", mock_guesserror)
    assert not Pb.main()

    monkeypatch.setattr(Pb, "command_line_ui", mock_systemexit)
    assert not Pb.main()

    monkeypatch.setattr(Pb, "command_line_ui", mock_fileexistserror)
    monkeypatch.setattr(br, "send_traceback", lambda *_: True)
    assert not Pb.main()


# ######################  loose command line ui helpers ###################### #
def test_exit(monkeypatch, capsys, pb_resources):
    class MockUsage(object):
        @staticmethod
        def increment(*args):
            print(args)
            return True

        @staticmethod
        def save():
            return True

    monkeypatch.setattr(br, "Usage", MockUsage)
    test_in_args = deepcopy(in_args)
    test_in_args.list_ids = [True]

    with pytest.raises(SystemExit):
        Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"))
    out, err = capsys.readouterr()
    assert "('PhyloBuddy', '%s', 'list_ids', 2412)" % Pb.VERSION.short() in out


def test_error(capsys, pb_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.show_unique = [True]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("o k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "PhyloBuddy object should have exactly 2 trees." in err
