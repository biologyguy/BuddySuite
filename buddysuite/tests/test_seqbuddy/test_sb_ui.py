#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: tests.py
version: 1.1
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Collection of PyTest unit tests for the AlignBuddy.py program
"""

import pytest
import os
import sys
import argparse
from copy import deepcopy
from unittest import mock

sys.path.insert(0, os.path.abspath("../"))
try:
    from buddysuite import buddy_resources as br
    from buddysuite import SeqBuddy as Sb
    from buddysuite import AlignBuddy as Alb
    from buddysuite import MyFuncs
except ImportError:
    import MyFuncs
    import SeqBuddy as Sb
    import buddy_resources as br

TEMP_DIR = MyFuncs.TempDir()
VERSION = Sb.VERSION


def fmt(prog):
    return br.CustomHelpFormatter(prog)

parser = argparse.ArgumentParser(prog="SeqBuddy", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                 description='''\
\033[1mSeqBuddy\033[m
  See your sequence files. Be your sequence files.

033[1mUsage examples\033[m:
  SeqBuddy.py "/path/to/seq_file" -<cmd>
  SeqBuddy.py "/path/to/seq_file" -<cmd> | SeqBuddy.py -<cmd>
  SeqBuddy.py "ATGATGCTAGTC" -f "raw" -<cmd>
''')

br.flags(parser, ("sequence", "Supply file path(s) or raw sequence. If piping sequences "
                              "into SeqBuddy this argument can be left blank."),
         br.sb_flags, br.sb_modifiers, VERSION)

# This is to allow py.test to work with its own flags
in_args = parser.parse_args([])


# ##################### '-ano', '--annotate' ###################### ##
def test_annotate_ui(capsys, sb_resources, sb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.annotate = ["misc_feature", "1-100,200-250", "+"]
    Sb.command_line_ui(test_in_args, sb_resources.get_one("d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "d22c44cf1a53624b58a86b0fb98c33a6"

    test_in_args = deepcopy(in_args)
    test_in_args.annotate = ["misc_feature", "1-100,200-250", "foo=bar", "hello=world", "-", "α4"]
    Sb.command_line_ui(test_in_args, sb_resources.get_one("d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "8c56cb3d6950ea43ce25f0c402355834"

    test_in_args = deepcopy(in_args)
    test_in_args.annotate = ["unknown_feature_that_is_t0o_long", "1-100,200-250", "foo=bar", "hello=world", "-", "α4"]
    Sb.command_line_ui(test_in_args, sb_resources.get_one("d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Warning: The provided annotation type is not part of the GenBank format standard" in err
    assert "Warning: Feature type is longer than 16 characters" in err


# ######################  '-asl', '--ave_seq_length' ###################### #
def test_ave_seq_length_ui(capsys, sb_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.ave_seq_length = [False]
    Sb.command_line_ui(test_in_args, sb_resources.get_one("p f"), True)
    out, err = capsys.readouterr()
    assert out == '428.38\n'

    test_in_args.ave_seq_length = ['clean']
    Sb.command_line_ui(test_in_args, sb_resources.get_one("p f"), True)
    out, err = capsys.readouterr()
    assert out == '427.38\n'


# ######################  '-cc', '--count_codons' ###################### #
def test_count_codons_ui(capsys, sb_resources, sb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.count_codons = ["foo"]
    Sb.command_line_ui(test_in_args, sb_resources.get_one("d g"), True)
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "5661f0ce92bb6cfba4519a61e0a838ed"

    test_in_args.count_codons = ["conc"]
    Sb.command_line_ui(test_in_args, sb_resources.get_one("d g"), True)
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "3e76bd510de4a61efb17ffc186ef9e68"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, sb_resources.get_one("p g"))
    out, err = capsys.readouterr()
    assert "Nucleic acid sequence required, not protein" in err


# ######################  '-efs', '--extract_feature_sequences' ###################### #
def test_extact_feature_sequences_ui(capsys, sb_resources, sb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.extract_feature_sequences = [["CDS"]]
    Sb.command_line_ui(test_in_args, sb_resources.get_one("d g"), True)
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "7e8a80caf902575c5eb3fc6ba8563956"

    test_in_args.extract_feature_sequences = [["TMD"]]
    Sb.command_line_ui(test_in_args, sb_resources.get_one("d g"), True)
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "13944b21484d5ea22af4fe57cc8074df"

    test_in_args.extract_feature_sequences = [["TMD", "splice_a"]]
    Sb.command_line_ui(test_in_args, sb_resources.get_one("d g"), True)
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "78629d308a89b458fb02e71d5568c978"

    test_in_args.extract_feature_sequences = [["foo"]]
    Sb.command_line_ui(test_in_args, sb_resources.get_one("d g"), True)
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "3cdbd5c8790f12871f8e04e40e315c93"
