#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

import pytest
from hashlib import md5
import os
import re
from copy import deepcopy
from Bio.Alphabet import IUPAC

try:
    import workshop.AlignBuddy as Alb
except ImportError:
    import AlignBuddy as Alb
import MyFuncs

write_file = MyFuncs.TempFile()


def align_to_hash(_alignbuddy, mode='hash'):
    if _alignbuddy.out_format == "phylipi":
        write_file.write(Alb.phylipi(_alignbuddy, "relaxed"))
    elif _alignbuddy.out_format == "phylipis":
        write_file.write(Alb.phylipi(_alignbuddy, "strict"))
    else:
        _alignbuddy.write(write_file.path)

    align_string = "{0}\n".format(write_file.read().strip())

    if mode != "hash":
        return align_string

    _hash = md5(align_string.encode()).hexdigest()
    return _hash


root_dir = os.getcwd()


def resource(file_name):
    return "{0}/unit_test_resources/{1}".format(root_dir, file_name)


align_files = ["Mnemiopsis_cds.nex", "Mnemiopsis_cds.phy", "Mnemiopsis_cds.phyr", "Mnemiopsis_cds.stklm",
               "Mnemiopsis_pep.nex", "Mnemiopsis_pep.phy", "Mnemiopsis_pep.phyr", "Mnemiopsis_pep.stklm",

               "Alignments_pep.phy", "Alignments_pep.phyr", "Alignments_pep.stklm"]


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_file(align_file):
    assert type(Alb.AlignBuddy(resource(align_file))) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_handle(align_file):
    with open(resource(align_file), 'r') as ifile:
        assert type(Alb.AlignBuddy(ifile)) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_raw(align_file):
    with open(resource(align_file), 'r') as ifile:
        assert type(Alb.AlignBuddy(ifile.read())) == Alb.AlignBuddy


# Now that we know that all the files are being turned into SeqBuddy objects okay, make them all objects so it doesn't
# need to be done over and over for each subsequent test.
sb_objects = [Alb.AlignBuddy(resource(x)) for x in align_files]