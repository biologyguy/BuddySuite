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

sys.path.insert(0, os.path.abspath("../"))
try:
    from buddysuite import buddy_resources as br
    from buddysuite import PhyloBuddy as Pb
    from buddysuite import AlignBuddy as Alb
    from buddysuite import MyFuncs
except ImportError:
    import buddy_resources as br
    import PhyloBuddy as Pb
    import AlignBuddy as Alb
    import MyFuncs

VERSION = Pb.VERSION
WRITE_FILE = MyFuncs.TempFile()


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
