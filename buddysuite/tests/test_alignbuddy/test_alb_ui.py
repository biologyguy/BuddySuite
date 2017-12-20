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
import re
from copy import deepcopy

import AlignBuddy as Alb
import SeqBuddy as Sb
import buddy_resources as br

TEMP_DIR = br.TempDir()
VERSION = Sb.VERSION


def mock_raisekeyboardinterrupt(*args, **kwargs):
    raise KeyboardInterrupt("Fake KeyboardInterrupt: %s, %s" % (args, kwargs))


def mock_raisesystemexit(*args, **kwargs):
    raise SystemExit("Fake SystemExit: %s, %s" % (args, kwargs))


def mock_raisesystemerror(*args, **kwargs):
    raise SystemError("Could not find: %s, %s" % (args, kwargs))


def mock_raiseattribute(*args, **kwargs):
    raise AttributeError("is not a supported alignment tool: %s, %s" % (args, kwargs))


def mock_raiseruntimeerror(*args, **kwargs):
    raise RuntimeError("Fake RuntimeError: %s, %s" % (args, kwargs))


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
def test_argparse_init(capsys, monkeypatch, alb_resources, hf, alb_odd_resources):
    monkeypatch.setattr(sys, 'argv', ['AlignBuddy.py', alb_resources.get_one("o p py", "paths"), "-con",
                                      "-o", "stockholm"])
    temp_in_args, alignbuddy = Alb.argparse_init()
    assert hf.buddy2hash(alignbuddy) == "5d9a03d9e1b4bf72d991257d3a696306"

    monkeypatch.setattr(sys, 'argv', ['AlignBuddy.py', alb_resources.get_one("o p py", "paths"), "-con", "-o", "foo"])
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "Format type 'foo' is not recognized/supported" in err

    monkeypatch.setattr(sys, 'argv', ['AlignBuddy.py', alb_odd_resources["dna"]["single"]["fasta"], "-con"])
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "GuessError: Could not determine format from _input file" in err

    monkeypatch.setattr(sys, 'argv', ['AlignBuddy.py', alb_odd_resources["dna"]["single"]["fasta"],
                                      "-con", "-f", "phylip"])
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "ValueError: First line should have two integers" in err

    monkeypatch.setattr(sys, 'argv', ['AlignBuddy.py', alb_odd_resources["dna"]["single"]["phylipss_recs"], "-con",
                                      "-f", "phylipss"])
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "PhylipError: Malformed Phylip --> 9 sequences expected, 8 found." in err

    monkeypatch.setattr(sys, 'argv', ['AlignBuddy.py', alb_resources.get_one("o p py", "paths"), "-con", "-f", "foo"])
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "TypeError: Format type 'foo' is not recognized/supported" in err

    monkeypatch.setattr(sys, 'argv', ['AlignBuddy.py', alb_resources.get_one("o p f", "paths"),
                                      "--quiet", "--generate_alignment", "mafft", "--reorder"])
    temp_in_args, alignbuddy = Alb.argparse_init()
    assert alignbuddy == []

    monkeypatch.setattr(sys, 'argv', ['AlignBuddy.py', alb_resources.get_one("o p f", "paths"),
                                      "--generate_alignment", "mafft", "--op", "5", "--quiet"])
    temp_in_args, alignbuddy = Alb.argparse_init()
    assert temp_in_args.generate_alignment == [['mafft', ' --op', '5']]

    monkeypatch.setattr(sys, 'argv', ['AlignBuddy.py', alb_resources.get_one("o p f", "paths"),
                                      "--generate_alignment", "mafft", "-q"])
    temp_in_args, alignbuddy = Alb.argparse_init()
    assert temp_in_args.generate_alignment == [['mafft']]


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
def test_clean_seqs_ui(capsys, alb_resources, alb_odd_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.clean_seq = [[None]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p pr"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "73b5d11dd25dd100648870228ab10d3d"

    test_in_args.clean_seq = [['strict', 'X']]
    Alb.command_line_ui(test_in_args, Alb.AlignBuddy(alb_odd_resources['dna']['single']['ambiguous']), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "6755ea1408eddd0e5f267349c287d989"


# ##################### '-cta', '--concat_alignments' ###################### ##
def test_concat_alignments_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.concat_alignments = [[]]

    tester = Sb.SeqBuddy("%s/Cnidaria_pep.nexus" % hf.resource_path)
    Sb.pull_recs(tester, "Ccr|Cla|Hec")
    tester = Alb.AlignBuddy(str(tester))
    tester.alignments.append(tester.alignments[0])
    tester.set_format("genbank")
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "86349e715f41e0bdd91bbd1dc0914769"

    test_in_args.concat_alignments = [["(.).(.)-Panx(.)"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "cd2b6594b22c431aea67fa45899f933a"

    test_in_args.concat_alignments = [["...", "Panx.*"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "e49b26f695c910a93f93d70563fd9dd9"

    test_in_args.concat_alignments = [[3, "Panx.*"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "e49b26f695c910a93f93d70563fd9dd9"

    test_in_args.concat_alignments = [[-9, "Panx.*"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "aaa9d9b717a5f79cfdf5d2666fb0f687"

    test_in_args.concat_alignments = [[3, 3]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "2f37a4e395162032bf43fab291c882f4"

    test_in_args.concat_alignments = [[3, -3]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "7fa8cd803df82414a5e1e190916456d8"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("p o g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Please provide at least two alignments." in err

    test_in_args.concat_alignments = [["foo"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p c"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "No match found for record" in err


# ##################### '-con', '--consensus' ###################### ##
def test_consensus_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)

    test_in_args.consensus = ['SiMpL']
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "9b69e5fb65ca1512de5a17472d105500"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "d1a8f7e629a020f5130373d7af65f9d9"

    test_in_args.consensus = ['WeiG']
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "2364532e0ec2465ea27f04acc5d0e61b"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "bf50c95916e9d62c95a460bbc517c053"

    test_in_args.consensus = ['foo']
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "No valid consensus mode" in err


# ######################################  '-dinv', '--delete_invariant_sites' ####################################### #
def test_delete_invariant_sites_ui(capsys, hf, alb_odd_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_invariant_sites = [[]]
    tester = Alb.AlignBuddy(alb_odd_resources['dna']['single']['ambiguous'])
    Alb.command_line_ui(test_in_args, tester, skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "27233a416437eabc72aa5d57cb695036"


# ##################### '-dr', '--delete_records' ###################### ##
def test_delete_records_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_records = ["α[1-5]", "β[A-M]"]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "de5beddbc7f0a7f8e3dc2d5fd43b7b29"
    assert hf.string2hash(err) == "31bb4310333851964015e21562f602c2"

    test_in_args.delete_records = ["α[1-5]", "β[A-M]", 4]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(err) == "ce6c9b29c95ba853eb444de5c71aeca9"

    test_in_args.delete_records = ["foo"]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "No sequence identifiers match 'foo'\n" in err

    tmp_file = br.TempFile()
    with open(tmp_file.path, "w", encoding="utf-8") as ofile:
        ofile.write('''\
α[1-5]
β[A-M]
''')
    test_in_args.delete_records = [tmp_file.path]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "de5beddbc7f0a7f8e3dc2d5fd43b7b29"


# ##################### '-et', '--enforce_triplets' ###################### ##
def test_enforce_triplets_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.enforce_triplets = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "d30529911c2ffdfb49152797225e3ff0"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p c"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Nucleic acid sequence required, not protein." in err


# ######################  '-efs', '--extract_feature_sequences' ###################### #
def test_extact_feature_sequences_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.extract_feature_sequences = [["CDS"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "842d9c911a33c0fd0484383eabefb0fe"

    test_in_args.extract_feature_sequences = [["TMD"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "cc7a1c6a22f721ec0668fc8ea6b23429"

    test_in_args.extract_feature_sequences = [["TMD1", "splice_a"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "497d536b1be9a90ef0ef75281d0c867f"

    test_in_args.extract_feature_sequences = [["TMD2:TMD3"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "07773f4fb1dc430c0c3ce6cd5a799439"

    test_in_args.extract_feature_sequences = [["foo"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "ac15492b38ca2ac4baa63e63a9b747f7"


# ##################### '-er', '--extract_regions' ###################### ##
def test_extract_regions_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.extract_regions = [["100:200", "250", ":10/50"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "7e202f9f6caaf83c558631a3577af20a"

    test_in_args.extract_regions = [["100:200", "250", ":10/foo"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), True)
    out, err = capsys.readouterr()
    assert "Unable to decode the positions string" in err


# ##################### '-fa', '--faux_align' ###################### ##
def test_faux_align_ui(capsys, alb_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.faux_align = [None]
    test_in_args.alignments = [alb_resources.get_one("o p g", "paths")]

    Alb.command_line_ui(test_in_args, Alb.AlignBuddy, skip_exit=True)
    out, err = capsys.readouterr()
    alignbuddy = Alb.AlignBuddy(out)
    assert len(alignbuddy.alignments[0][0]) == 625


# ##################### '-ga', '--generate_alignment' ###################### ##
@pytest.mark.generate_alignments
def test_generate_alignment_ui(capsys, monkeypatch, sb_resources, alb_resources, hf):
    monkeypatch.setattr(Alb, "generate_msa", lambda *_: alb_resources.get_one("o d g"))
    test_in_args = deepcopy(in_args)
    test_in_args.generate_alignment = [[]]

    test_in_args.alignments = [sb_resources.get_one("d g", "paths")]
    test_in_args.out_format = "gb"
    Alb.command_line_ui(test_in_args, Alb.AlignBuddy, skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "842d9c911a33c0fd0484383eabefb0fe"


@pytest.mark.generate_alignments
def test_generate_alignment_ui_patch_path(monkeypatch, capsys, sb_resources):
    class MockClass(object):
        @staticmethod
        def isatty():
            return True

    class TextIOWrapper(object):
        def __init__(self):
            self.buffer = MockClass()
            self.buffer.raw = MockClass()

    test_in_args = deepcopy(in_args)
    test_in_args.generate_alignment = [[]]
    monkeypatch.setenv("PATH", [])
    with pytest.raises(AttributeError) as err:
        Alb.command_line_ui(test_in_args, Alb.AlignBuddy, pass_through=True)
    assert "Unable to identify any supported alignment tools on your system." in str(err)

    monkeypatch.setattr(Alb, "which", lambda *_: True)
    test_in_args.alignments = [sb_resources.get_one("d g", "paths")]

    test_in_args.generate_alignment = [["mafft"]]
    monkeypatch.setattr(Alb, 'generate_msa', mock_raiseattribute)
    with pytest.raises(AttributeError):
        Alb.command_line_ui(test_in_args, Alb.AlignBuddy, pass_through=True)

    monkeypatch.setattr(Alb, 'generate_msa', mock_raisesystemerror)
    with pytest.raises(SystemError):
        Alb.command_line_ui(test_in_args, Alb.AlignBuddy, pass_through=True)

    monkeypatch.setattr(Alb, "TextIOWrapper", TextIOWrapper)
    text_wrapper = TextIOWrapper()
    test_in_args.alignments = [text_wrapper]
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, Alb.AlignBuddy)
    out, err = capsys.readouterr()
    assert "Warning: No input detected so AlignBuddy is aborting..." in err


# ######################  '-gh', '--generate_hmm' ###################### #
@br.skip_windows
def test_generate_hmm_ui(alb_resources, hf, capsys, monkeypatch):
    test_in_args = deepcopy(in_args)
    test_in_args.generate_hmm = [[]]

    tester = alb_resources.get_one("m p c")
    Alb.command_line_ui(test_in_args, tester, True)

    out, err = capsys.readouterr()
    hmm = re.findall("(HMM +A +C +D.+?//)", out, re.DOTALL)
    assert "COMPO   2.68250  3.88919  3.04853  2.78121  3.08118  3.13138  3.72607  2.65113  2.67024  2.34849  3.43272" \
           "  3.05512  3.56890  3.09868  3.03335  2.74953  2.90269  2.58958  4.30351  3.11199" in hmm[0]
    assert "COMPO   2.61975  3.93095  3.12640  2.80659  3.03969  2.94881  3.78599  2.73397  2.73613  2.36723  3.48106" \
           "  3.11755  3.38828  3.15135  3.06078  2.68581  2.82442  2.59321  4.24683  3.17591" in hmm[1]

    tester = alb_resources.get_one("m d c")
    Alb.command_line_ui(test_in_args, tester, True)

    out, err = capsys.readouterr()
    hmm = re.findall("(HMM +A +C +G.+?//)", out, re.DOTALL)
    assert """\
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   1.37149  1.47979  1.42806  1.27722
          1.38629  1.38629  1.38629  1.38629
          0.10249  4.35641  2.47000  1.46634  0.26236  0.00000        *
      1   0.06560  3.99042  3.73531  3.85677      1 A - - -
          1.38629  1.38629  1.38629  1.38629
          0.02802  4.28194  4.28194  1.46634  0.26236  2.15125  0.12368""" in hmm[0]

    assert """\
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   1.26618  1.65166  1.51488  1.18245
          0.93669  1.43354  1.84356  1.55419
          0.72237  1.59420  1.16691  3.55520  0.02899  0.00000        *
      1   0.60107  2.55837  1.28853  2.31598     85 a - - -
          1.38629  1.38629  1.38629  1.38629
          0.03300  4.12082  4.12082  1.46634  0.26236  3.38099  0.03461""" in hmm[1]

    test_in_args.generate_hmm = ["foo"]
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, tester)

    out, err = capsys.readouterr()
    assert "Could not find foo on your system. Please check your spelling or install HMMER3." in err

    def raise_syserror(*_, **__):
        raise SystemError("No output detected after running")

    monkeypatch.setattr(Alb, "generate_hmm", raise_syserror)
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, tester)

    out, err = capsys.readouterr()
    assert "No output detected after running" in err

    def raise_valueerror(*_, **__):
        raise ValueError("This is just a check for dealing with unexpected error")

    monkeypatch.setattr(Alb, "generate_hmm", raise_valueerror)

    with pytest.raises(ValueError):
        Alb.command_line_ui(test_in_args, tester)


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
def test_list_ids(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.list_ids = [False]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "f087f9c1413ba66c28fb0fccf7c974e6"

    test_in_args.list_ids = [3]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "4d85249a1f187d38d411a78ced65a98c"


# #################################### '-lc', '--lowercase' ################################### #
def test_lowercase_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.lowercase = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "00661f7afb419c6bb8c9ac654af7c976"


# ##################### '-mf2a', '--map_features2alignment' ###################### ##
def test_map_features2alignment_ui(capsys, alb_resources, sb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.mapfeat2align = [sb_resources.get_one("d g", "paths")]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d n"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "c1359470ad0916902c4e96facd088378"


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
def test_order_ids_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.order_ids = [False]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "f36f73e973aa7a2dcd2fc86f239d5a23"

    test_in_args.order_ids = ['rev']
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "d4dcdc5059fd82c6b9cc44a66770b801"


# ##################### '-pi', '--percent_id' ###################### ##
def test_percent_id_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.percent_id = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "8d8a52ebeacf68069773784162cf6d54"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "55553113f5ee206041f085488029d4b5"


# ###########################################  '-pfm', '--pos_freq_mat' ############################################ #
def test_position_frequency_matrix_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.pos_freq_mat = True

    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "a229988a3fdae4f089b91f5047ad816d", print(out)


# ##################### '-pr', '--pull_records' ###################### ##
def test_pull_records_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.pull_records = [["α[1-5]$", "β[A-M]"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "2de557d6fd3dc6cd1bf43a1995392a4c"
    assert err == ""

    test_in_args.pull_records = [["α[1-5]$", "ML218922a", "full"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "fb82ffec15ece60a74d9ac8db92d2999"


# ######################  '-ri', '--rename_ids' ###################### #
def test_rename_ids_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.rename_ids = [["[a-z](.)", "?\\1", 2]]
    tester = alb_resources.get_one("m p s")
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "888f2e3feb9e67f9bc008183082c822a"

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
def test_reverse_transcribe_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.reverse_transcribe = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o r n"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "f8c2b216fa65fef9c74c1d0c4abc2ada"

    with pytest.raises(TypeError) as err:
        Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), pass_through=True)
    assert "RNA sequence required, not IUPACAmbiguousDNA()." in str(err)


# ######################  '-sf', '--screw_formats' ###################### #
hashes = [("fasta", "cfa898d43918055b6a02041195874da9"), ("gb", "28c76bfc1ae74b2a55c3044287e074a8"),
          ("nexus", "49bf9b3f56104e4f19048523d725f025"), ("phylip", "968ed9fa772e65750f201000d7da670f"),
          ("phylipr", "5064c1d6ae6192a829972b7ec0f129ed"), ("phylipss", "4bd927145de635c429b2917e0a1db176"),
          ("phylipsr", "b46b57ede57f12c3c3b906681882f81a"), ("stockholm", "5d9a03d9e1b4bf72d991257d3a696306"),
          ("clustal", "9d328711cf6f6750c33373a912efb521")]


@pytest.mark.parametrize("_format,next_hash", hashes)
def test_screw_formats_ui(_format, next_hash, capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = _format
    tester = alb_resources.get_one("o p py")
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == next_hash


hashes = [("clustal", "cf349d6061c602439b72b51368f694ed"), ("phylip", "2a77f5761d4f51b88cb86b079e564e3b"),
          ("phylipr", "1f172a3beef76e8e3d42698bb2c3c87d"), ("phylipss", "eb82cda31fcb2cf00e11d7e910fde695"),
          ("phylipsr", "368169cb86c6ddb7074ed89e2d42c4dd"), ("stockholm", "f221b9973aef4771169136a30bd030fa")]


@pytest.mark.parametrize("_format,next_hash", hashes)
def test_screw_formats_ui2(_format, next_hash, capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = _format
    tester = alb_resources.get_one("m p py")
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == next_hash


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
def test_transcribe_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.transcribe = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d n"), skip_exit=True)
    out, err = capsys.readouterr()

    assert hf.string2hash(out) == "e531dc31f24192f90aa1f4b6195185b0"

    with pytest.raises(TypeError) as err:
        Alb.command_line_ui(test_in_args, alb_resources.get_one("o r n"), pass_through=True)
    assert "DNA sequence required, not IUPACAmbiguousRNA()." in str(err)


# ##################### '-tr', '--translate' ###################### ##
def test_translate_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.translate = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "542794541324d74ff636eaf4ee5e6b1a"

    with pytest.raises(TypeError) as err:
        Alb.command_line_ui(test_in_args, alb_resources.get_one("o p n"), pass_through=True)
    assert "Nucleic acid sequence required, not protein." in str(err)


# ##################### '-trm', '--trimal' ###################### ##
def test_trimal_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.trimal = [False]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "b5cb0f31ab3bed2cf93fb9c1eac52be7"

    test_in_args.trimal = ["gappyout"]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "b5cb0f31ab3bed2cf93fb9c1eac52be7"

    test_in_args.trimal = [0.25]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d psr"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "5df948e4b2cb6c0d0740984445655135"

    test_in_args.trimal = ["foo"]
    with pytest.raises(NotImplementedError) as err:
        Alb.command_line_ui(test_in_args, alb_resources.get_one("o p n"), pass_through=True)
    assert "foo not an implemented trimal method" in str(err)


# #################################### '-uc', '--uppercase' ################################### #
def test_uppercase_ui(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.uppercase = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "6f3f234d796520c521cb85c66a3e239a"


# ######################  main() ###################### #
def test_main(monkeypatch, capsys, alb_resources):
    test_in_args = deepcopy(in_args)
    test_in_args.enforce_triplets = True
    monkeypatch.setattr(Alb, "argparse_init", lambda: [in_args, alb_resources.get_one("o d f")])
    monkeypatch.setattr(Alb, "command_line_ui", lambda *_: True)
    assert Alb.main()

    monkeypatch.setattr(Alb, "command_line_ui", mock_raisekeyboardinterrupt)
    assert not Alb.main()
    out, err = capsys.readouterr()
    assert "Fake KeyboardInterrupt" in out

    monkeypatch.setattr(Alb, "command_line_ui", mock_raisesystemexit)
    assert not Alb.main()

    monkeypatch.setattr(Alb, "command_line_ui", mock_raiseruntimeerror)
    monkeypatch.setattr(br, "send_traceback", lambda *_: True)
    assert not Alb.main()


# ######################  loose command line ui helpers ###################### #
def test_test(capsys, alb_resources, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_records = ["α1"]
    test_in_args.test = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d f"), skip_exit=True)
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "d41d8cd98f00b204e9800998ecf8427e"


def test_inplace(capsys, alb_resources, hf):
    tmp_dir = br.TempDir()
    tester = alb_resources.get_one("o d f")
    tester.write("%s/align" % tmp_dir.path)

    test_in_args = deepcopy(in_args)
    test_in_args.transcribe = True
    test_in_args.in_place = True
    test_in_args.alignments = ["%s/align" % tmp_dir.path]

    Alb.command_line_ui(test_in_args, tester, skip_exit=True)
    out, err = capsys.readouterr()
    tester = Alb.AlignBuddy("%s/align" % tmp_dir.path)
    assert "File overwritten at:" in err
    assert hf.buddy2hash(tester) == "8f78e0c99e2d6d7d9b89b8d854e02bcd", tester.write("temp.del")

    test_in_args.alignments = ["I/do/not/exist"]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d f"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Warning: The -i flag was passed in, but the positional argument doesn't seem to be a file." in err
