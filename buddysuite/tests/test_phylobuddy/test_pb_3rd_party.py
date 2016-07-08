#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests 3rd party Phylogenetic tree estimation software that is wrapped by PhyloBuddy
"""

import pytest
import os
from unittest import mock
from subprocess import Popen, PIPE
import re

from ... import AlignBuddy as Alb
from ... import PhyloBuddy as Pb
from ... import buddy_resources as br


# ######################  'gt', '--generate_trees' ###################### #
# ######### RAxML ########### #
raxml_version = Popen("raxml -v", shell=True, stdout=PIPE).communicate()[0].decode()
raxml_version = re.search("([0-9]+\.[0-9]+\.[0-9]+)", raxml_version).group(0)
if raxml_version not in ["8.2.4", "8.2.8"]:
    raise ValueError("Untested RAxML version (%s). Please update the tests as necessary." % raxml_version)


def test_raxml_inputs_nuc(alb_resources, pb_helpers):
    tester = alb_resources.get_one("o d n")
    tester = Pb.generate_tree(tester, 'raxml')
    assert pb_helpers.phylo2hash(tester) in ['4e083fc6d8b4f4342674efd93d5e313c', '1cede6c576bb88125e2387d850f813ab']


def test_raxml_inputs_quiet(alb_resources, pb_helpers):
    tester = alb_resources.get_one("o d n")
    tester = Pb.generate_tree(tester, 'raxml', quiet=True)
    assert pb_helpers.phylo2hash(tester) in ['4e083fc6d8b4f4342674efd93d5e313c', '1cede6c576bb88125e2387d850f813ab']


def test_raxml_inputs_pep(alb_resources, pb_helpers):
    tester = alb_resources.get_one("o p py")
    tester = Pb.generate_tree(tester, 'raxml')
    assert pb_helpers.phylo2hash(tester) in ['41ccd8852002f3b98c032378403c38b0', '6443e4dddb9b64a783bd9e97b369e0d4']


def test_raxml_multi_param(alb_resources, pb_helpers):
    tester = alb_resources.get_one("o d n")
    tester = Pb.generate_tree(tester, 'raxml', '-m GTRCAT -p 112358 -K MK -N 2')
    assert pb_helpers.phylo2hash(tester) in ['1999ef930c0d6cadbe5844b7c6355029', '53ea2002e19d1d88c684f0ddc02de187']


def test_raxml_multiple_searches(alb_resources, pb_helpers):
    tester = alb_resources.get_one("o d py")
    tester = Pb.generate_tree(tester, 'raxml', '-N 2')
    assert pb_helpers.phylo2hash(tester) in ['4f722d2a473604c6fab0bc63366128dc', '59a3f6cd366d9918b36ac95c04572a2f']


# def test_raxml_multiple_trees(alb_resources, pb_helpers):
#    tester = alb_resources.get_one("m p n")
#    tester = Pb.generate_tree(tester, 'raxml')
#    assert pb_helpers.phylo2hash(tester) == ''


# ######### PhyML ######### #
phyml_version = Popen("phyml --version", shell=True, stdout=PIPE).communicate()[0].decode()
phyml_version = re.search("([0-9]+\.[0-9]+\.[0-9]+)|([0-9]+)", phyml_version).group(0)
if phyml_version not in ["20120412", "20160207", "3.2.20160701"]:
    raise ValueError("Untested PhyML version (%s). Please update the tests as necessary." % phyml_version)


def test_phyml_dna(alb_resources, pb_helpers):
    # Nucleotide
    tester = alb_resources.get_one("o d py")
    tester = Pb.generate_tree(tester, 'phyml', '-m GTR --r_seed 12345')
    assert pb_helpers.phylo2hash(tester) in ['b61e75e4706d35e92f2208d438f52771',
                                             'b0bdb3f5bf1fb2e44bef3c16f80c38f2',
                                             'b9d3f11e332c3589110322e939aa41cc'], print(tester)


def test_phyml_pep(alb_resources, pb_helpers):
    # Peptide
    tester = alb_resources.get_one("o p py")
    tester = Pb.generate_tree(tester, 'phyml', '-m Blosum62 --r_seed 12345')
    assert pb_helpers.phylo2hash(tester) in ['7caa5c641fa83085c2980efca875112a',
                                             '2bf0a204b2de7bc5132aa7073ecfb011',
                                             '981d16e8f02989a8642095016c88af90'], print(tester)


# ######### FastTree ######### #
fasttree_version = Popen("fasttree", shell=True, stderr=PIPE).communicate()[1].decode()
fasttree_version = re.search("([0-9]+\.[0-9]+\.[0-9]+)", fasttree_version).group(0)
if fasttree_version not in ["2.1.8", "2.1.9"]:
    raise ValueError("Untested FastTree version (%s). Please update the tests as necessary." % fasttree_version)


def test_fasttree_inputs(alb_resources, pb_helpers):
    temp_dir = br.TempDir()
    # Nucleotide
    alignbuddy = alb_resources.get_one("o d n")

    tester = Pb.generate_tree(Alb.make_copy(alignbuddy), 'FastTree', '-seed 12345')
    assert pb_helpers.phylo2hash(tester) in ['da8a67cae6f3f70668f7cf04060b7cd8', '732c5e9a978cebb1cfce6af6d64950c2']

    tester = Pb.generate_tree(alignbuddy, 'fasttree', '-seed 12345', quiet=True)
    assert pb_helpers.phylo2hash(tester) in ['da8a67cae6f3f70668f7cf04060b7cd8', '732c5e9a978cebb1cfce6af6d64950c2']

    alignbuddy = alb_resources.get_one("o p n")
    tester = Pb.generate_tree(alignbuddy, 'fasttree', '-seed 12345', keep_temp="%s/new_dir" % temp_dir.path)
    assert pb_helpers.phylo2hash(tester) in ['82d5a9d4f44fbedf29565686a7cdfcaa', '682210ef16beedee0e9f43c05edac112']


# def test_fasttree_multi_param(alb_resources, pb_helpers):
#    temp_file = br.TempFile()
#    tester = alb_resources.get_one("m d pr")
#    tester = Pb.generate_tree(tester, 'FastTree', '-seed 12345 -wag -fastest -log %s' % temp_file.path)
#    assert pb_helpers.phylo2hash(tester) == 'c6f02fe52d111a89120878d36b8fc506'


def test_generate_trees_edge_cases(alb_resources):
    temp_file = br.TempFile()
    tester = alb_resources.get_one("o d n")
    with pytest.raises(FileExistsError):
        Pb.generate_tree(tester, "raxml", keep_temp=temp_file.path)

    with pytest.raises(AttributeError):
        Pb.generate_tree(tester, "foo")

    with pytest.raises(FileNotFoundError):
        Pb.generate_tree(tester, "raxml", "-m BINCATLG")

    with pytest.raises(RuntimeError):
        Pb.generate_tree(tester, "fasttree", "-s 12345")

    # noinspection PyUnresolvedReferences
    with mock.patch.dict(os.environ, {"PATH": ""}):
        with pytest.raises(ProcessLookupError):
            Pb.generate_tree(tester, "raxml")
