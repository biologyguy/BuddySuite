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

import AlignBuddy as Alb
import PhyloBuddy as Pb
import buddy_resources as br


# ######################  'gt', '--generate_trees' ###################### #
# ######### RAxML ########### #
@br.skip_windows
def get_raxml_version():
    raxml_version = Popen("raxml -v", shell=True, stdout=PIPE).communicate()[0].decode()
    raxml_version = re.search("([0-9]+\.[0-9]+\.[0-9]+)", raxml_version).group(0)
    if raxml_version not in ["7.6.6", "8.2.4", "8.2.8", "8.2.9", "8.2.10"]:
        raise ValueError("Untested RAxML version (%s). Please update the tests as necessary." % raxml_version)


get_raxml_version()


@br.skip_windows
def test_raxml_inputs_nuc(alb_resources, hf):
    tester = alb_resources.get_one("o d n")
    tester = Pb.generate_tree(tester, 'raxml')
    assert hf.buddy2hash(tester) in ['7569f9f6c7f8079579bfb77291b99616', '4e083fc6d8b4f4342674efd93d5e313c',
                                     '706ba436f8657ef3aee7875217dd07c0', '1cede6c576bb88125e2387d850f813ab',
                                     '3280f0b404d62fbedd2e6090cfd6fedc', 'b954440b4baa338eb4a63b0a5a7f7942',
                                     '629fa27f6e636c137257abc6a8a40956', '878e1ca74539e9b969c539bd8743c27e']


@br.skip_windows
def test_raxml_inputs_quiet(alb_resources, hf):
    tester = alb_resources.get_one("o d n")
    tester = Pb.generate_tree(tester, 'raxml', quiet=True)
    assert hf.buddy2hash(tester) in ['7569f9f6c7f8079579bfb77291b99616', '4e083fc6d8b4f4342674efd93d5e313c',
                                     '706ba436f8657ef3aee7875217dd07c0', '1cede6c576bb88125e2387d850f813ab',
                                     '3280f0b404d62fbedd2e6090cfd6fedc', 'b954440b4baa338eb4a63b0a5a7f7942',
                                     '629fa27f6e636c137257abc6a8a40956', '878e1ca74539e9b969c539bd8743c27e']


@br.skip_windows
def test_raxml_inputs_pep(alb_resources, hf):
    tester = alb_resources.get_one("o p py")
    tester = Pb.generate_tree(tester, 'raxml')
    assert hf.buddy2hash(tester) in ['3e6ab2efc088d5547fa8244462f7cc03', '832f3b301a2a320affb9864b4c6a3979',
                                     '41ccd8852002f3b98c032378403c38b0', '6443e4dddb9b64a783bd9e97b369e0d4',
                                     'caff20f9bb5192a799ec18db0faa8569', '31218fd9e1b803df09fce22eee8da62e']


@br.skip_windows
def test_raxml_multi_param(alb_resources, hf):
    tester = alb_resources.get_one("o d n")
    tester = Pb.generate_tree(tester, 'raxml', '-m GTRCAT -p 112358 -K MK -N 3')
    assert hf.buddy2hash(tester) in ['08fa932a0cbb33d936ef4c8aef3c0095', '53ea2002e19d1d88c684f0ddc02de187',
                                     '02e5ea7b756b68b8622636ba7e28e15b', 'c24b5e9c1899246b9a459a61efe0aad5',
                                     '1999ef930c0d6cadbe5844b7c6355029', 'ead385a55d3cad48ecfdd1ebddb00b8b',
                                     'ae6e53a6e7c03786928b087a15e0435e', '135f244c3de9ba417fe685715fef6807',
                                     '886a9d841e68c2fd4c6052ad85c4eaf7']


@br.skip_windows
def test_raxml_multiple_searches(alb_resources, hf):
    tester = alb_resources.get_one("o d py")
    tester = Pb.generate_tree(tester, 'raxml', '-N 3')
    assert hf.buddy2hash(tester) in ['76356987f7e2368cdf13c42567cb7453', 'ee223e46a9f9753203553e6fd7473ec9',
                                     'b3a8359c62e9d29b952782df53a4782a', '807de55171690b2af1b724d86390cbc7',
                                     '46e0a6d38e3c44b93a04d2b5b7f4aca1', 'ecf0cabfab48bc2449ca3e43645b3d36',
                                     'ce60f6627ff665994f0c6836d7469632', 'c8f5fe6c50150e55b43d545e44c7c5ce']


# def test_raxml_multiple_trees(alb_resources, hf):
#    tester = alb_resources.get_one("m p n")
#    tester = Pb.generate_tree(tester, 'raxml')
#    assert hf.buddy2hash(tester) == ''


# ######### PhyML ######### #
@br.skip_windows
def get_phyml_version():
    phyml_version = Popen("phyml --version", shell=True, stdout=PIPE).communicate()[0].decode()
    phyml_version = re.search("([0-9]+\.[0-9]+\.[0-9]+)|([0-9]+)", phyml_version).group(0)
    if phyml_version not in ["20111216", "20120412", "20131022", "20160207",
                             "3.2.20160701", "3.2.20160531", "3.3.20170530"]:
        raise ValueError("Untested PhyML version (%s). Please update the tests as necessary." % phyml_version)


get_phyml_version()


@br.skip_windows
def test_phyml_dna(alb_resources, hf):
    # Nucleotide
    tester = alb_resources.get_one("o d py")
    tester = Pb.generate_tree(tester, 'phyml', '-m GTR --r_seed 12345', r_seed=12345)
    assert hf.buddy2hash(tester) in {'b61e75e4706d35e92f2208d438f52771': "",
                                     'b0bdb3f5bf1fb2e44bef3c16f80c38f2': "",
                                     'b9d3f11e332c3589110322e939aa41cc': "",
                                     '754c38fab99c01c68a68c0a59248d242': "",
                                     '3ca772c34cdcf0a22c09e1592aba9ebf': "",
                                     'd7ae1badd31d48487276495bad4522e5': "",
                                     'e84fb949f1a6ed0296cda4e5a8422423': "3.3.20170530"}, print(tester)


@br.skip_windows
def test_phyml_pep(alb_resources, hf):
    # Peptide
    tester = alb_resources.get_one("o p py")
    tester = Pb.generate_tree(tester, 'phyml', '-m Blosum62 --r_seed 12345', r_seed=12345)
    assert hf.buddy2hash(tester) in {'7caa5c641fa83085c2980efca875112a': "",
                                     '2bf0a204b2de7bc5132aa7073ecfb011': "",
                                     '981d16e8f02989a8642095016c88af90': "",
                                     'd8ee3631002b6603d08272c2b44fd21c': "",
                                     '03acc8e899955f7e838852d7d71049ad': "",
                                     'abe46f3bac533ad2f510bd4657aa9505': "",
                                     '06f5ec5db5e1a27c07e4e8b5f5850685': "3.3.20170530"}, print(tester)


# ######### FastTree ######### #
@br.skip_windows
def get_fasttree_version():
    fasttree_version = Popen("fasttree", shell=True, stderr=PIPE).communicate()[1].decode()
    fasttree_version = re.search("([0-9]+\.[0-9]+\.[0-9]+)", fasttree_version).group(0)
    if fasttree_version not in ["2.1.4", "2.1.8", "2.1.9", "2.1.10"]:
        raise ValueError("Untested FastTree version (%s). Please update the tests as necessary." % fasttree_version)


get_fasttree_version()


@br.skip_windows
def test_fasttree_inputs(alb_resources, hf):
    temp_dir = br.TempDir()
    # Nucleotide
    alignbuddy = alb_resources.get_one("o d n")

    tester = Pb.generate_tree(Alb.make_copy(alignbuddy), 'FastTree', '-seed 12345')
    assert hf.buddy2hash(tester) in ['d7f505182dd1a1744b45cc326096f70c', 'da8a67cae6f3f70668f7cf04060b7cd8',
                                     '732c5e9a978cebb1cfce6af6d64950c2']

    tester = Pb.generate_tree(alignbuddy, 'fasttree', '-seed 12345', quiet=True)
    assert hf.buddy2hash(tester) in ['d7f505182dd1a1744b45cc326096f70c', 'da8a67cae6f3f70668f7cf04060b7cd8',
                                     '732c5e9a978cebb1cfce6af6d64950c2']

    alignbuddy = alb_resources.get_one("o p n")
    tester = Pb.generate_tree(alignbuddy, 'fasttree', '-seed 12345', keep_temp="%s%snew_dir" % (temp_dir.path, os.sep))
    assert hf.buddy2hash(tester) in ['57eace9bdd2074297cbd2692c1f4cd38', '82d5a9d4f44fbedf29565686a7cdfcaa',
                                     '682210ef16beedee0e9f43c05edac112']


# def test_fasttree_multi_param(alb_resources, hf):
#    temp_file = br.TempFile()
#    tester = alb_resources.get_one("m d pr")
#    tester = Pb.generate_tree(tester, 'FastTree', '-seed 12345 -wag -fastest -log %s' % temp_file.path)
#    assert hf.buddy2hash(tester) == 'c6f02fe52d111a89120878d36b8fc506'


# ######### IQ-TREE ######### #
@br.skip_windows
def get_iqtree_version():
    iqtree_version = Popen("iqtree -h", shell=True, stderr=PIPE, stdout=PIPE).communicate()[0].decode()
    iqtree_version = re.search("version ([0-9]+\.[0-9]+\.[0-9]+)", iqtree_version).group(1)
    if iqtree_version not in ["1.5.5", "1.6.1"]:
        raise ValueError("Untested IQ-TREE version (%s). Please update the tests as necessary." % iqtree_version)


get_iqtree_version()


@br.skip_windows
def test_iqtree_inputs(alb_resources, hf, capsys):
    temp_dir = br.TempDir()
    # Nucleotide
    alignbuddy = alb_resources.get_one("o d py")

    tester = Pb.generate_tree(alignbuddy, 'iqtree', '-m JC -seed 12345 -nt 1', r_seed=12345)
    capsys.readouterr()
    assert hf.buddy2hash(tester) in {'08d1c7473503e52af5750caa79db1a0e': "1.5.5",
                                     '476efd414108a677c8ae0b8dba55c677': "1.6.1"}, print(tester)

    alignbuddy = alb_resources.get_one("o p py")
    tester = Pb.generate_tree(alignbuddy, 'iqtree', '-m LG -seed 12345 -nt 1',
                              keep_temp=os.path.join(temp_dir.path, "new_dir"), r_seed=12345)
    assert hf.buddy2hash(tester) in {'4b4c79ddcae0836db33bc118f0b17292': '1.5.5',
                                     '2abb71ad0bc5b61bcf1437e2485669b4': "1.6.1"}, print(tester)

    root, dirs, files = next(os.walk(os.path.join(temp_dir.path, 'new_dir')))
    assert sorted(files) == ['pb_input.aln', 'pb_input.aln.bionj', 'pb_input.aln.ckp.gz',
                             'pb_input.aln.iqtree', 'pb_input.aln.log', 'pb_input.aln.mldist', 'pb_input.aln.treefile']


@br.skip_windows
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
