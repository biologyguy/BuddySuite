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

try:
    from buddysuite import AlignBuddy as Alb
    from buddysuite import PhyloBuddy as Pb
    from buddysuite import buddy_resources as br
except ImportError:
    import AlignBuddy as Alb
    import PhyloBuddy as Pb
    import buddy_resources as br


# ######################  'gt', '--generate_trees' ###################### #
# ######### RAxML ########### #
raxml_version = Popen("raxml -v", shell=True, stdout=PIPE).communicate()[0].decode()
raxml_version = re.search("([0-9]+\.[0-9]+\.[0-9]+)", raxml_version).group(0)
if raxml_version not in ["8.2.4", "8.2.8"]:
    raise ValueError("Untested RAxML version (%s). Please update the tests as necessary." % raxml_version)


def test_raxml_inputs_nuc(alb_resources, pb_helpers):
    tester = alb_resources.get_one("o d n")
    tester = Pb.generate_tree(tester, 'raxml')
    assert pb_helpers.phylo_to_hash(tester) in ['4e083fc6d8b4f4342674efd93d5e313c', '1cede6c576bb88125e2387d850f813ab']


def test_raxml_inputs_quiet(alb_resources, pb_helpers):
    tester = alb_resources.get_one("o d n")
    tester = Pb.generate_tree(tester, 'raxml', quiet=True)
    assert pb_helpers.phylo_to_hash(tester) in ['4e083fc6d8b4f4342674efd93d5e313c', '1cede6c576bb88125e2387d850f813ab']


def test_raxml_inputs_pep(alb_resources, pb_helpers):
    tester = alb_resources.get_one("o p py")
    tester = Pb.generate_tree(tester, 'raxml')
    assert pb_helpers.phylo_to_hash(tester) in ['41ccd8852002f3b98c032378403c38b0', '6443e4dddb9b64a783bd9e97b369e0d4']


def test_raxml_multi_param(alb_resources, pb_helpers):
    tester = alb_resources.get_one("o d n")
    tester = Pb.generate_tree(tester, 'raxml', '-m GTRCAT -p 112358 -K MK -N 2')
    assert pb_helpers.phylo_to_hash(tester) in ['1999ef930c0d6cadbe5844b7c6355029', '53ea2002e19d1d88c684f0ddc02de187']


def test_raxml_multiple_searches(alb_resources, pb_helpers):
    tester = alb_resources.get_one("o d py")
    tester = Pb.generate_tree(tester, 'raxml', '-N 2')
    assert pb_helpers.phylo_to_hash(tester) in ['4f722d2a473604c6fab0bc63366128dc', '59a3f6cd366d9918b36ac95c04572a2f']


#def test_raxml_multiple_trees(alb_resources, pb_helpers):
#    tester = alb_resources.get_one("m p n")
#    tester = Pb.generate_tree(tester, 'raxml')
#    assert pb_helpers.phylo_to_hash(tester) == ''

"""  Need to sort this out for version 20160207
# ######### PhyML ######### #
phyml_version = Popen("phyml --version", shell=True, stdout=PIPE).communicate()[0].decode()
phyml_version = re.search("([0-9]+)", phyml_version).group(0)
if phyml_version not in ["20120412"]:
    raise ValueError("Untested PhyML version (%s). Please update the tests as necessary." % phyml_version)


def test_phyml_inputs(alb_resources, pb_helpers):
    # Nucleotide
    tester = alb_resources.get_one("o d n")
    tester = Pb.generate_tree(tester, 'phyml', '-m GTR --r_seed 12345')
    assert pb_helpers.phylo_to_hash(tester) == 'd3a4e7601998885f333ddd714ca764db'
    # Peptide
    tester = alb_resources.get_one("o p n")
    tester = Pb.generate_tree(tester, 'phyml', '-m Blosum62 --r_seed 12345')
    assert pb_helpers.phylo_to_hash(tester) == '52c7d028341b250bcc867d57a68c794c'


def test_phyml_multi_param(alb_resources, pb_helpers):
    tester = alb_resources.get_one("o d n")
    tester = Pb.generate_tree(tester, 'phyml', '-m GTR -o tl -b 2 --r_seed 12345')
    assert pb_helpers.phylo_to_hash(tester) == '5434f29509eab76dd52dd69d2c0e186f'
"""

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
    assert pb_helpers.phylo_to_hash(tester) in ['da8a67cae6f3f70668f7cf04060b7cd8', '732c5e9a978cebb1cfce6af6d64950c2']

    tester = Pb.generate_tree(alignbuddy, 'fasttree', '-seed 12345', quiet=True)
    assert pb_helpers.phylo_to_hash(tester) in ['da8a67cae6f3f70668f7cf04060b7cd8', '732c5e9a978cebb1cfce6af6d64950c2']

    alignbuddy = alb_resources.get_one("o p n")
    tester = Pb.generate_tree(alignbuddy, 'fasttree', '-seed 12345', keep_temp="%s/new_dir" % temp_dir.path)
    assert pb_helpers.phylo_to_hash(tester) in ['82d5a9d4f44fbedf29565686a7cdfcaa', '682210ef16beedee0e9f43c05edac112']


#def test_fasttree_multi_param(alb_resources, pb_helpers):
#    temp_file = br.TempFile()
#    tester = alb_resources.get_one("m d pr")
#    tester = Pb.generate_tree(tester, 'FastTree', '-seed 12345 -wag -fastest -log %s' % temp_file.path)
#    assert pb_helpers.phylo_to_hash(tester) == 'c6f02fe52d111a89120878d36b8fc506'


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
