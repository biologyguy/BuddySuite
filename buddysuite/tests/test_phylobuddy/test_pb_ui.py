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
import ete3

from ... import buddy_resources as br
from ... import PhyloBuddy as Pb
from ... import AlignBuddy as Alb

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


# ###################### argparse_init() ###################### #
def test_argparse_init(capsys, pb_odd_resources, pb_helpers):
    sys.argv = ['PhyloBuddy.py', pb_odd_resources["compare"]]
    temp_in_args, phylobuddy = Pb.argparse_init()
    assert pb_helpers.string2hash(str(phylobuddy)) == "d8e14a2bfc8e9c0ac3c524f5fb478c67"

    sys.argv += ["-f", "foo"]
    with pytest.raises(SystemExit):
        Pb.argparse_init()

    out, err = capsys.readouterr()
    assert "Error: The format 'foo' passed in with the -f flag is not recognized." in err


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


# ###################### 'cpt', '--collapse_polytomies' ###################### #
def test_collapse_polytomies_ui(capsys, pb_odd_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.collapse_polytomies = [[20]]
    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(pb_odd_resources['support']), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "1b0979265205b17ca7f34abbd02f6e26"

    test_in_args.collapse_polytomies = [[0.1, 'length']]
    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(pb_odd_resources['support']), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "252572f7b9566c62df24d57065412240"


# ###################### 'ct', '--consensus_tree' ###################### #
def test_consensus_tree_ui(capsys, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.consensus_tree = [False]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "acd3fb34cce867c37684244701f9f5bf"

    test_in_args.consensus_tree = [0.9]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "baf9b2def2c0fa2ff97d3a16d24b4738"

    test_in_args.consensus_tree = [1.5]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "acd3fb34cce867c37684244701f9f5bf"
    assert err == "Warning: The frequency value should be between 0 and 1. Defaulting to 0.5.\n\n"


# ###################### 'dt', '--display_trees' ###################### #
def test_display_trees_ui(monkeypatch, pb_resources):
    if 'DISPLAY' in os.environ:
        test_in_args = deepcopy(in_args)
        test_in_args.display_trees = True
        show = mock.Mock(return_value=True)
        monkeypatch.setattr(ete3.TreeNode, "show", show)
        Pb.command_line_ui(test_in_args, pb_resources.get_one("o k"), skip_exit=True)


def test_display_trees_ui_no_display(capsys, monkeypatch, pb_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.display_trees = True
    show = mock.Mock(return_value=True)
    monkeypatch.setattr(ete3.TreeNode, "show", show)
    # noinspection PyUnresolvedReferences
    with mock.patch.dict('os.environ'):
        if 'DISPLAY' in os.environ:
            del os.environ['DISPLAY']
        Pb.command_line_ui(test_in_args, pb_resources.get_one("o k"), skip_exit=True)
        out, err = capsys.readouterr()

    assert "Error: Your system is non-graphical, so display_trees can not work." in err


# ###################### 'dis', '--distance' ###################### #
def test_distance_ui(capsys, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.distance = [False]  # Should default to wrf
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "66df36fa117e4c4660f04e2649c3fa6b"
    assert err == "Tree 1	Tree 2	Value\n"

    test_in_args.distance = ["wrf"]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "66df36fa117e4c4660f04e2649c3fa6b"
    assert err == "Tree 1	Tree 2	Value\n"

    test_in_args.distance = ["uwrf"]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "9dd162b4cc4fa28f402e6f31ef2fb349"
    assert err == "Tree 1	Tree 2	Value\n"

    test_in_args.distance = ["ed"]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "a49e54a6c14ab4dbbf73d3f7d1d6aa82"
    assert err == "Tree 1	Tree 2	Value\n"

    test_in_args.distance = ["foo"]
    with pytest.raises(AttributeError) as err:
        Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
        assert "foo is an invalid comparison method." in err

"""
# ###################### 'gt', '--generate_tree' ###################### #
@pytest.mark.generate_trees
def test_generate_tree_ui1(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.in_format, test_in_args.out_format = "nexus", "newick"
    test_in_args.trees = [resource("Mnemiopsis_cds.nex")]

    test_in_args.generate_tree = [["fasttree", "-seed 12345"]]
    Pb.command_line_ui(test_in_args, [], skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "d7f505182dd1a1744b45cc326096f70c"


@pytest.mark.generate_trees
def test_generate_tree_ui2():
    test_in_args = deepcopy(in_args)
    test_in_args.in_format, test_in_args.out_format = "nexus", "newick"
    test_in_args.trees = [resource("Mnemiopsis_cds.nex")]

    test_in_args.generate_tree = [["fasttree", "-s 12345"]]  # Command doesn't exist in fasttree
    with pytest.raises(SystemExit):
        Pb.command_line_ui(test_in_args, [])


@pytest.mark.generate_trees
def test_generate_tree_ui3():
    test_in_args = deepcopy(in_args)
    test_in_args.in_format, test_in_args.out_format = "nexus", "newick"
    test_in_args.trees = [resource("Mnemiopsis_cds.nex")]

    test_in_args.generate_tree = [["foo"]]
    with pytest.raises(SystemExit):
        Pb.command_line_ui(test_in_args, [])
"""


# ###################### 'hi', '--hash_ids' ###################### #
def test_hash_ids_ui(capsys, monkeypatch, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.hash_ids = [[1, "nodes"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("o n"), skip_exit=True)
    out, err = capsys.readouterr()

    assert pb_helpers.string2hash(out) != pb_helpers.phylo2hash(pb_resources.get_one("o n"))
    assert "Warning: The hash_length parameter was passed in with the value 1" in err

    test_in_args.hash_ids = [[-1, "nodes"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m n"), skip_exit=True)
    out, err = capsys.readouterr()

    assert pb_helpers.string2hash(out) != pb_helpers.phylo2hash(pb_resources.get_one("m n"))
    assert "Warning: The hash_length parameter was passed in with the value -1" in err

    def hash_ids(*args):
        if args:
            pass
        raise ValueError("Foo bar")

    test_in_args.hash_ids = [[1, "nodes"]]
    monkeypatch.setattr(Pb, "hash_ids", hash_ids)
    with pytest.raises(ValueError):
        Pb.command_line_ui(test_in_args, pb_resources.get_one("o n"))
    out, err = capsys.readouterr()
    assert not err


# ###################### 'li', '--list_ids' ###################### #
def test_list_ids_ui(capsys, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.list_ids = [True]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "4f0857e0211ff0cd058d0cf7cbaf64d5"

    test_in_args.list_ids = [4]
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "bf2fbfe1bd52e9b27ae21f5c06e7763a"

    test_in_args.list_ids = [4]
    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy("(, );"), skip_exit=True)
    out, err = capsys.readouterr()
    assert out == "#### tree_1 ####\nNone\n\n"


# ###################### '-nt', '--num_tips' ###################### #
def test_num_tips_ui(capsys, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.num_tips = True

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "f43920f4df66e76fbacae4af178eeebb"


# ###################### 'ptr', '--print_trees' ###################### #
def test_print_trees_ui(capsys, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.print_trees = True

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "fe340117cb8f573100c00fc897e6c8ce"


# ###################### 'pt', '--prune_taxa' ###################### #
def test_prune_taxa_ui(capsys, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.prune_taxa = [["fir"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "99635c6dbf708f94cf4dfdca87113c44"

    test_in_args.prune_taxa = [["fir", "ovi"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m n"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "2a385fa95024323fea412fd2b3c3e91f"


# ###################### 'ri', '--rename_ids' ###################### #
def test_rename_ids_ui(capsys, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.rename_ids = ['Mle', 'Phylo']

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "6843a620b725a3a0e0940d4352f2036f"


# ###################### 'rt', '--root' ###################### #
def test_root_ui_midpoint(capsys, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.root = [[]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "8a7fdd9421e0752c9cd58a1e073186c7"


def test_root_ui_leaf(capsys, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.root = [["firSA25a"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "f32bdc34bfe127bb0453a80cf7b01302"


def test_root_ui_mrca(capsys, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.root = [["firSA25a", "penSH31b"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "757489907bd5c82882d084ffcb22cfba"


# ###################### 'sf', '--screw_formats' ###################### #
def test_screw_formats_ui(capsys, pb_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = "nexus"

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "543d2fc90ca1f391312d6b8fe896c59c"


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
    assert "File over-written at:" in err
    check_file = os.path.isfile("%s.nex" % temp_file.path)
    assert check_file

    test_in_args.trees[0] = "%s.nex" % temp_file.path
    test_in_args.screw_formats = "newick"

    tester = Pb.PhyloBuddy("%s.nex" % temp_file.path)
    Pb.command_line_ui(test_in_args, tester, skip_exit=True)
    out, err = capsys.readouterr()
    assert "File over-written at:" in err
    check_file = os.path.isfile("%s.nwk" % temp_file.path)
    assert check_file


# ###################### 'su', '--show_unique' ###################### #
def test_show_unique_ui(capsys, pb_resources, pb_helpers, pb_odd_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.show_unique = True

    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(pb_odd_resources["compare"]), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "ea5b0d1fcd7f39cb556c0f5df96281cf"

    with pytest.raises(SystemExit):
        Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"))

    out, err = capsys.readouterr()
    assert err == "AssertionError: PhyloBuddy object should have exactly 2 trees.\n"


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
def test_unroot_ui(capsys, pb_odd_resources, pb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.unroot = True

    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(pb_odd_resources["figtree"]), skip_exit=True)
    out, err = capsys.readouterr()
    assert pb_helpers.string2hash(out) == "10e9024301b3178cdaed0b57ba33f615"
