#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: alignbuddy_tests.py
version: 1, alpha
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Collection of PyTest unit tests for the buddy_resources.py module
"""

import pytest
from subprocess import Popen, PIPE

import buddy_resources as Br


def test_versions():
    sb_ver = Popen("sb -v", stdout=PIPE, shell=True).communicate()
    assert sb_ver[0].decode() == str(Br.VERSIONS["SeqBuddy"])

    alb_ver = Popen("alb -v", stdout=PIPE, shell=True).communicate()
    assert alb_ver[0].decode() == str(Br.VERSIONS["AlignBuddy"])

    pb_ver = Popen("pb -v", stdout=PIPE, shell=True).communicate()
    assert pb_ver[0].decode() == str(Br.VERSIONS["PhyloBuddy"])

    db_ver = Popen("db -v", stdout=PIPE, shell=True).communicate()
    assert db_ver[0].decode() == str(Br.VERSIONS["DatabaseBuddy"])


# ######################  '_format_to_extension' ###################### #
def test_format_to_extension():
    ext_dict = {'fasta': 'fa', 'fa': 'fa', 'genbank': 'gb', 'gb': 'gb', 'newick': 'nwk', 'nwk': 'nwk', 'nexus': 'nex',
                'nex': 'nex', 'phylip': 'phy', 'phy': 'phy', 'phylip-relaxed': 'phyr', 'phyr': 'phyr',
                'phylipss': 'physs', 'physs': 'physs', 'phylipsr': 'physr', 'physr': 'physr', 'stockholm': 'stklm',
                'stklm': 'stklm', 'clustal': 'clus', 'clus': 'clus'}
    for i, j in ext_dict.items():
        assert j == Br.format_to_extension(i)


def testphylip_sequential_out():
    tester = alignments.get_one("m p py")
    output = Br.phylip_sequential_out(tester, relaxed=False)
    assert string2hash("%s\n" % output.rstrip()) == "eb82cda31fcb2cf00e11d7e910fde695"

    tester = alignments.get_one("m p pr")
    output = Br.phylip_sequential_out(tester, relaxed=True)
    assert string2hash("%s\n" % output.rstrip()) == "a16c6e617e5a88fef080eea54e54e8a8"

    records = tester.records()
    with pytest.raises(br.PhylipError) as e:
        records[0].id = str(records[1].id)
        Br.phylip_sequential_out(tester)
    assert "Malformed Phylip --> Repeat id" in str(e)

    with pytest.raises(br.PhylipError) as e:
        records[1].id = "foo"
        records[0].seq = Seq(str(records[0].seq[:-2]), records[0].seq.alphabet)
        Br.phylip_sequential_out(tester)
    assert "Malformed Phylip --> The length of record 'foo' is incorrect" in str(e)

    with pytest.raises(br.PhylipError) as e:
        records[0].seq = Seq("AT%s" % records[0].seq, records[0].seq.alphabet)
        Br.phylip_sequential_out(tester, relaxed=False)
    assert "Malformed Phylip --> Repeat id" in str(e) and "after strict truncation." in str(e)


def testphylip_sequential_read():
    tester = alignments.get_one("m p pss")
    aligns = Br.phylip_sequential_read(str(tester), relaxed=False)
    assert len(Br.AlignBuddy(aligns).records()) == 29

    with pytest.raises(br.PhylipError) as e:
        Br.phylip_sequential_read(str(tester))
    assert "Malformed Phylip --> 8 sequences expected, 4 found" in str(e)

    tester = alignments.get_one("m p psr")
    aligns = Br.phylip_sequential_read(str(tester))
    assert len(Br.AlignBuddy(aligns).records()) == 34

    with pytest.raises(br.PhylipError) as e:
        with open(resource("malformed_phylip_columns.physs"), "r") as ifile:
            Br.phylip_sequential_read(ifile.read(), relaxed=False)
    assert "Malformed Phylip --> Less sequence found than expected" in str(e)

    with pytest.raises(br.PhylipError) as e:
        with open(resource("malformed_phylip_records_more.physs"), "r") as ifile:
            Br.phylip_sequential_read(ifile.read(), relaxed=False)
    assert "Malformed Phylip --> 7 sequences expected, 8 found." in str(e)

    with pytest.raises(br.PhylipError) as e:
        with open(resource("malformed_phylip_columns_more.physs"), "r") as ifile:
            Br.phylip_sequential_read(ifile.read(), relaxed=False)
    assert "Malformed Phylip --> Sequence Mle-Panxα2 has 2046 columns, 2043 expected." in str(e)

    with pytest.raises(br.PhylipError) as e:
        with open(resource("malformed_phylip_repeat_id.physr"), "r") as ifile:
            Br.phylip_sequential_read(ifile.read(), relaxed=True)
    assert "Malformed Phylip --> Repeat ID Mle-Pxα1." in str(e)

    with pytest.raises(br.PhylipError) as e:
        with open(resource("malformed_phylip_repeat_id.physr"), "r") as ifile:
            Br.phylip_sequential_read(ifile.read(), relaxed=False)
    assert "Malformed Phylip --> Repeat id 'Mle-Pxα1  ' after strict truncation." in str(e)


    # ######################  'shift_features' ###################### #
def test_shift_features():
    tester = Sb._make_copy(sb_objects[1])
    features = tester.records[0].features
    assert string2hash(str(Sb._shift_features(features, 3, 1203))) == "cc002df59db8a04f47cb5c764f8a1e1f"

    tester = Sb._make_copy(sb_objects[1])
    features = tester.records[0].features
    assert string2hash(str(Sb._shift_features(features, -50, 1203))) == "86dcc19cd51ba3acdaf48e8c15bf7f1a"

    tester = Sb._make_copy(sb_objects[1])
    features = tester.records[0].features
    features[0].location.parts = []
    assert string2hash(str(Sb._shift_features(features, 3, 1203))) == "467ee934600573e659053b8117447986"

    tester = Sb._make_copy(sb_objects[1])
    features = tester.records[0].features
    features[0].location.parts = [FeatureLocation(3, 50, strand=+1)]
    assert string2hash(str(Sb._shift_features(features, 3, 1203))) == "94bf6774e97fcecaef858cbdf811def4"

    with pytest.raises(TypeError):
        tester = Sb._make_copy(sb_objects[1])
        features = tester.records[0].features
        features[0].location = [dict]
        Sb._shift_features(features, 3, 1203)