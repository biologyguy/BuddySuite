#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: test_alb_ui.py
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

from ... import AlignBuddy as Alb
from ... import SeqBuddy as Sb
from ... import buddy_resources as br

TEMP_DIR = br.TempDir()
VERSION = Sb.VERSION


def fmt(prog):
    return br.CustomHelpFormatter(prog)

parser = argparse.ArgumentParser(prog="alignBuddy", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                 description='''\
\033[1mAlignBuddy\033[m
  Sequence alignments with a splash of kava.

\033[1mUsage examples\033[m:
  AlignBuddy.py "/path/to/align_file" -<cmd>
  AlignBuddy.py "/path/to/align_file" -<cmd> | AlignBuddy.py -<cmd>
  AlignBuddy.py "/path/to/seq_file" -ga "mafft" -p "--auto --thread 8"
''')

br.flags(parser, ("alignments", "The file(s) you want to start working on"),
         br.alb_flags, br.alb_modifiers, VERSION)

# This is to allow py.test to work with its own flags
in_args = parser.parse_args([])


# ###################### argparse_init() ###################### #
def test_argparse_init(capsys, alb_resources, alb_helpers, alb_odd_resources):
    sys.argv = ['AlignBuddy.py', alb_resources.get_one("o p py", "paths"), "-con", "-o", "stockholm"]
    temp_in_args, alignbuddy = Alb.argparse_init()
    assert alb_helpers.align2hash(alignbuddy) == "5d9a03d9e1b4bf72d991257d3a696306"

    sys.argv = ['AlignBuddy.py', alb_resources.get_one("o p py", "paths"), "-con", "-o", "foo"]
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "Format type 'foo' is not recognized/supported" in err

    sys.argv = ['AlignBuddy.py', alb_odd_resources["dna"]["single"]["fasta"], "-con"]
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "GuessError: Could not determine format from _input file" in err

    sys.argv = ['AlignBuddy.py', alb_odd_resources["dna"]["single"]["fasta"], "-con", "-f", "phylip"]
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "ValueError: First line should have two integers" in err

    sys.argv = ['AlignBuddy.py', alb_odd_resources["dna"]["single"]["phylipss_recs"], "-con", "-f", "phylipss"]
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "PhylipError: Malformed Phylip --> 9 sequences expected, 8 found." in err

    sys.argv = ['AlignBuddy.py', alb_resources.get_one("o p py", "paths"), "-con", "-f", "foo"]
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "TypeError: Format type 'foo' is not recognized/supported" in err

    sys.argv = ['AlignBuddy.py', alb_resources.get_one("o p f", "paths"),
                "--quiet", "--generate_alignment", "mafft", "--reorder"]
    temp_in_args, alignbuddy = Alb.argparse_init()
    assert alignbuddy == []


# ##################### '-al', '--alignment_lengths' ###################### ##
def test_alignment_lengths_ui(capsys, alb_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.alignment_lengths = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p c"), skip_exit=True)
    out, err = capsys.readouterr()
    assert out == "481\n683\n"
    assert err == "# Alignment 1\n# Alignment 2\n"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("o p py"), skip_exit=True)
    out, err = capsys.readouterr()
    assert out == "681\n"
    assert err == ""


# ##################### '-bts', '--bootstrap' ###################### ##
def test_bootstrap_ui(capsys, alb_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.bootstrap = [False]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    tester = Alb.AlignBuddy(out)
    assert tester.lengths() == [481, 683]

    test_in_args.bootstrap = [3]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    tester = Alb.AlignBuddy(out)
    assert tester.lengths() == [481, 481, 481, 683, 683, 683]


# ##################### '-cs', '--clean_seqs' ###################### ##
def test_clean_seqs_ui(capsys, alb_resources, alb_odd_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.clean_seq = [[None]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p pr"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "73b5d11dd25dd100648870228ab10d3d"

    test_in_args.clean_seq = [['strict', 'X']]
    Alb.command_line_ui(test_in_args, Alb.AlignBuddy(alb_odd_resources['dna']['single']['ambiguous']), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "6755ea1408eddd0e5f267349c287d989"


# ##################### '-cta', '--concat_alignments' ###################### ##
def test_concat_alignments_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.concat_alignments = [[]]

    tester = Sb.SeqBuddy("%s/Cnidaria_pep.nexus" % alb_helpers.resource_path)
    Sb.pull_recs(tester, "Ccr|Cla|Hec")
    tester = Alb.AlignBuddy(str(tester))
    tester.alignments.append(tester.alignments[0])
    tester.set_format("genbank")
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "d21940f3dad2295dd647f632825d8541"

    test_in_args.concat_alignments = [["(.).(.)-Panx(.)"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "5ac908ebf7918a45664a31da480fda58"

    test_in_args.concat_alignments = [["...", "Panx.*"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "e754350b0397cf54f531421d1e85774f"

    test_in_args.concat_alignments = [[3, "Panx.*"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "e754350b0397cf54f531421d1e85774f"

    test_in_args.concat_alignments = [[-9, "Panx.*"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "9d2886afc640d35618754e05223032a2"

    test_in_args.concat_alignments = [[3, 3]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "4e4101f9b5a6d44d524a9783a8c4004b"

    test_in_args.concat_alignments = [[3, -3]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "5d9d9ac8fae604be74c436e5f0b5b6db"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("p o g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Please provide at least two alignments." in err

    test_in_args.concat_alignments = [["foo"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p c"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "No match found for record" in err


# ##################### '-con', '--consensus' ###################### ##
def test_consensus_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.consensus = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "7b0aa3cca159b276158cf98209be7dab"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "89130797253646e61b78ab7d91ad3fd9"


# ##################### '-dr', '--delete_records' ###################### ##
def test_delete_records_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_records = ["α[1-5]", "β[A-M]"]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "de5beddbc7f0a7f8e3dc2d5fd43b7b29"
    assert alb_helpers.string2hash(err) == "31bb4310333851964015e21562f602c2"

    test_in_args.delete_records = ["α[1-5]", "β[A-M]", 4]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(err) == "ce6c9b29c95ba853eb444de5c71aeca9"

    test_in_args.delete_records = ["foo"]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "No sequence identifiers match 'foo'\n" in err


# ##################### '-et', '--enforce_triplets' ###################### ##
def test_enforce_triplets_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.enforce_triplets = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "6ff2a8a7c58bb6ac0d98fe373981e220"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p c"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Nucleic acid sequence required, not protein." in err


# ##################### '-er', '--extract_regions' ###################### ##
def test_extract_regions_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.extract_regions = [10, 110]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "3929c5875a58e9a1e64425d4989e590a"

    test_in_args.extract_regions = [110, 10]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "3929c5875a58e9a1e64425d4989e590a"

    test_in_args.extract_regions = [-110, 10]
    with pytest.raises(ValueError) as err:
        Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), pass_through=True)
    assert "Please specify positive integer indices" in str(err)


# ##################### '-ga', '--generate_alignment' ###################### ##
@pytest.mark.generate_alignments
def test_generate_alignment_ui(capsys, sb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.generate_alignment = [[]]

    test_in_args.alignments = [sb_resources.get_one("d g", "paths")]
    test_in_args.out_format = "gb"
    Alb.command_line_ui(test_in_args, Alb.AlignBuddy, skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "00579d6d858b85a6662b4c29094bf205"

    test_in_args.generate_alignment = [["foo"]]
    with pytest.raises(AttributeError) as err:
        Alb.command_line_ui(test_in_args, Alb.AlignBuddy, pass_through=True)
    assert "foo is not a supported alignment tool" in str(err)


def test_generate_alignment_ui_patch_path(monkeypatch):
    test_in_args = deepcopy(in_args)
    test_in_args.generate_alignment = [[]]
    monkeypatch.setenv("PATH", [])
    with pytest.raises(AttributeError) as err:
        Alb.command_line_ui(test_in_args, Alb.AlignBuddy, pass_through=True)
    assert "Unable to identify any supported alignment tools on your system." in str(err)


# ######################  '-hsi', '--hash_ids' ###################### #
def test_hash_seq_ids_ui(capsys, alb_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.hash_ids = [None]
    tester = alb_resources.get_one("m p s")
    ids = [rec.id for rec in tester.records()]
    Alb.command_line_ui(test_in_args, tester, True)
    for indx, rec in enumerate(tester.records()):
        assert rec.id != ids[indx]
        assert ids[indx] == tester.hash_map[rec.id]

    test_in_args.hash_ids = [0]
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert "Warning: The hash_length parameter was passed in with the value 0. This is not a positive integer" in err

    tester.alignments *= 10
    test_in_args.hash_ids = [1]
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert "cover all sequences, so it has been increased to 2" in err


# ###############################  '-li', '--list_ids' ############################## #
def test_list_ids(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.list_ids = [False]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "f087f9c1413ba66c28fb0fccf7c974e6"

    test_in_args.list_ids = [3]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "4d85249a1f187d38d411a78ced65a98c"


# #################################### '-lc', '--lowercase' ################################### #
def test_lowercase_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.lowercase = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "00661f7afb419c6bb8c9ac654af7c976"


# ##################### '-mf2a', '--map_features2alignment' ###################### ##
def test_map_features2alignment_ui(capsys, alb_resources, sb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.mapfeat2align = [sb_resources.get_one("d g", "paths")]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d n"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "9fece109249f4d787c13e6fb2742843d"


# ##############################################  '-ns', '--num_seqs'  ############################################### #
def test_num_seqs_ui(capsys, alb_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.num_seqs = True
    for alignbuddy in alb_resources.get_list("m d c pr psr s"):
        Alb.command_line_ui(test_in_args, alignbuddy, skip_exit=True)
        out, err = capsys.readouterr()
        assert out == "# Alignment 1\n8\n\n# Alignment 2\n21\n" or out == "# Alignment 1\n13\n\n# Alignment 2\n21\n"

    for alignbuddy in alb_resources.get_list("m p c pr psr s"):
        Alb.command_line_ui(test_in_args, alignbuddy, skip_exit=True)
        out, err = capsys.readouterr()
        assert out == "# Alignment 1\n20\n\n# Alignment 2\n13\n" or out == "# Alignment 1\n13\n\n# Alignment 2\n21\n"

    for alignbuddy in alb_resources.get_list("o p c pr psr s"):
        Alb.command_line_ui(test_in_args, alignbuddy, skip_exit=True)
        out, err = capsys.readouterr()
        assert out == "13\n"


# ##############################################  '-oi', '--order_ids' ############################################### #
def test_order_ids_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.order_ids = [False]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "f36f73e973aa7a2dcd2fc86f239d5a23"

    test_in_args.order_ids = ['rev']
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "d4dcdc5059fd82c6b9cc44a66770b801"


# ##################### '-pr', '--pull_records' ###################### ##
def test_pull_records_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.pull_records = [["α[1-5]$", "β[A-M]"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "2de557d6fd3dc6cd1bf43a1995392a4c"
    assert err == ""

    test_in_args.pull_records = [["α[1-5]$", "ML218922a", "full"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "fb82ffec15ece60a74d9ac8db92d2999"


# ######################  '-ri', '--rename_ids' ###################### #
def test_rename_ids_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.rename_ids = [["[a-z](.)", "?\\1", 2]]
    tester = alb_resources.get_one("m p s")
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "888f2e3feb9e67f9bc008183082c822a"

    test_in_args.rename_ids = [["[a-z](.)"]]
    with pytest.raises(AttributeError) as err:
        Alb.command_line_ui(test_in_args, tester, pass_through=True)
    assert "rename_ids requires two or three argments:" in str(err)

    test_in_args.rename_ids = [["[a-z](.)", "?\\1", 2, "foo"]]
    with pytest.raises(AttributeError) as err:
        Alb.command_line_ui(test_in_args, tester, pass_through=True)
    assert "rename_ids requires two or three argments:" in str(err)

    test_in_args.rename_ids = [["[a-z](.)", "?\\1", "foo"]]
    with pytest.raises(ValueError) as err:
        Alb.command_line_ui(test_in_args, tester, pass_through=True)
    assert "Max replacements argument must be an integer" in str(err)

    test_in_args.rename_ids = [["[a-z](.)", "?\\1\\2", 2]]
    with pytest.raises(AttributeError) as err:
        Alb.command_line_ui(test_in_args, tester, pass_through=True)
    assert "There are more replacement" in str(err)


# ##################### '-r2d', '--reverse_transcribe' ###################### ##
def test_reverse_transcribe_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.reverse_transcribe = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o r n"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "f8c2b216fa65fef9c74c1d0c4abc2ada"

    with pytest.raises(TypeError) as err:
        Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), pass_through=True)
    assert "RNA sequence required, not IUPACAmbiguousDNA()." in str(err)


# ######################  '-sf', '--screw_formats' ###################### #
hashes = [("fasta", "cfa898d43918055b6a02041195874da9"), ("gb", "ceac7a2a57aa8e3f530f70e2765f9ab2"),
          ("nexus", "49bf9b3f56104e4f19048523d725f025"), ("phylip", "968ed9fa772e65750f201000d7da670f"),
          ("phylipr", "5064c1d6ae6192a829972b7ec0f129ed"), ("phylipss", "4bd927145de635c429b2917e0a1db176"),
          ("phylipsr", "b46b57ede57f12c3c3b906681882f81a"), ("stockholm", "5d9a03d9e1b4bf72d991257d3a696306"),
          ("clustal", "9d328711cf6f6750c33373a912efb521")]


@pytest.mark.parametrize("_format,next_hash", hashes)
def test_screw_formats_ui(_format, next_hash, capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = _format
    tester = alb_resources.get_one("o p py")
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == next_hash

hashes = [("clustal", "cf349d6061c602439b72b51368f694ed"), ("phylip", "2a77f5761d4f51b88cb86b079e564e3b"),
          ("phylipr", "1f172a3beef76e8e3d42698bb2c3c87d"), ("phylipss", "eb82cda31fcb2cf00e11d7e910fde695"),
          ("phylipsr", "368169cb86c6ddb7074ed89e2d42c4dd"), ("stockholm", "f221b9973aef4771169136a30bd030fa")]


@pytest.mark.parametrize("_format,next_hash", hashes)
def test_screw_formats_ui2(_format, next_hash, capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = _format
    tester = alb_resources.get_one("m p py")
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == next_hash


def test_screw_formats_ui3(capsys, alb_resources):
    test_in_args = deepcopy(in_args)
    tester = alb_resources.get_one("o p py")
    tester.write("%s/seq.phy" % TEMP_DIR.path)
    test_in_args.in_place = True
    test_in_args.alignments = ["%s/seq.phy" % TEMP_DIR.path]
    test_in_args.screw_formats = "genbank"
    Alb.command_line_ui(test_in_args, tester, True)
    assert os.path.isfile("%s/seq.gb" % TEMP_DIR.path)

    test_in_args.in_place = False
    test_in_args.screw_formats = "foo"
    with pytest.raises(TypeError) as err:
        Alb.command_line_ui(test_in_args, tester, pass_through=True)
    assert "Format type 'foo' is not recognized/supported" in str(err)

    test_in_args.screw_formats = "gb"
    tester = alb_resources.get_one("m p py")
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "gb format does not support multiple alignments in one file." in err


# ##################################  '-stf', '--split_to_files' ################################### #
hashes = [("clustal", "m p c"), ("phylip", "m p py"),
          ("phylip-relaxed", "m p pr"), ("phylipss", "m p pss"),
          ("phylipsr", "m p psr"), ("stockholm", "m p s")]


@pytest.mark.parametrize("_format,_code", hashes)
def test_split_alignment_ui(_format, _code, alb_resources):
    TEMP_DIR.subdir(_format)
    test_in_args = deepcopy(in_args)
    test_in_args.split_to_files = [["%s/%s" % (TEMP_DIR.path, _format), "Foo_"]]
    tester = alb_resources.get_one(_code)
    for num in range(7):
        tester.alignments += tester.alignments
    Alb.command_line_ui(test_in_args, tester, True)
    assert os.path.isfile("%s/%s/Foo_001.%s" % (TEMP_DIR.path, _format, br.format_to_extension[_format]))
    assert os.path.isfile("%s/%s/Foo_061.%s" % (TEMP_DIR.path, _format, br.format_to_extension[_format]))
    assert os.path.isfile("%s/%s/Foo_121.%s" % (TEMP_DIR.path, _format, br.format_to_extension[_format]))
    TEMP_DIR.del_subdir(_format)


def test_split_alignment_ui2(capsys, alb_resources):
    TEMP_DIR.subdir("split_alignment")
    test_in_args = deepcopy(in_args)
    test_in_args.split_to_files = [["%s/split_alignment" % TEMP_DIR.path, "Foo_", "Bar"]]
    tester = alb_resources.get_one("m p c")
    for num in range(4):
        tester.alignments += tester.alignments
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert os.path.isfile("%s/split_alignment/Foo_01.clus" % TEMP_DIR.path)
    assert os.path.isfile("%s/split_alignment/Foo_16.clus" % TEMP_DIR.path)

    assert "Warning: Only one prefix can be accepted, 2 where provided. Using the first.\n" in err
    assert "New file: " in err

    tester = alb_resources.get_one("o d c")
    test_in_args.quiet = False
    with pytest.raises(ValueError) as err:
        Alb.command_line_ui(test_in_args, tester, pass_through=True)
    assert "Only one alignment present, nothing written." in str(err)
    TEMP_DIR.del_subdir("split_alignment")


# ##################### '-d2r', '--transcribe' ###################### ##
def test_transcribe_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.transcribe = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d n"), skip_exit=True)
    out, err = capsys.readouterr()

    assert alb_helpers.string2hash(out) == "e531dc31f24192f90aa1f4b6195185b0"

    with pytest.raises(TypeError) as err:
        Alb.command_line_ui(test_in_args, alb_resources.get_one("o r n"), pass_through=True)
    assert "DNA sequence required, not IUPACAmbiguousRNA()." in str(err)


# ##################### '-tr', '--translate' ###################### ##
def test_translate_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.translate = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "a949edce98525924dbbc3ced03c18214"

    with pytest.raises(TypeError) as err:
        Alb.command_line_ui(test_in_args, alb_resources.get_one("o p n"), pass_through=True)
    assert "Nucleic acid sequence required, not protein." in str(err)


# ##################### '-trm', '--trimal' ###################### ##
def test_trimal_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.trimal = [False]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "362577b8b42f18c9a4fa557e785d17e1"

    test_in_args.trimal = ["gappyout"]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "362577b8b42f18c9a4fa557e785d17e1"

    test_in_args.trimal = [0.25]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d psr"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "5df948e4b2cb6c0d0740984445655135"

    test_in_args.trimal = ["foo"]
    with pytest.raises(NotImplementedError) as err:
        Alb.command_line_ui(test_in_args, alb_resources.get_one("o p n"), pass_through=True)
    assert "foo not an implemented trimal method" in str(err)


# #################################### '-uc', '--uppercase' ################################### #
def test_uppercase_ui(capsys, alb_resources, alb_helpers):
    test_in_args = deepcopy(in_args)
    test_in_args.uppercase = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert alb_helpers.string2hash(out) == "6f3f234d796520c521cb85c66a3e239a"
