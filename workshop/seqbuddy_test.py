#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

import pytest
from hashlib import md5
from Bio import SeqIO
import os
import workshop.SeqBuddy as Sb
import MyFuncs

write_file = MyFuncs.TempFile()

def seqs_to_hash(_seqbuddy):
    if _seqbuddy.out_format == "phylipi":
        write_file.write(Sb.phylipi(_seqbuddy, "relaxed"))
    elif _seqbuddy.out_format == "phylipis":
        write_file.write(Sb.phylipi(_seqbuddy, "strict"))
    else:
        _seqbuddy.write(write_file.path)

    seqs_string = "{0}\n".format(write_file.read().strip()).encode()
    _hash = md5(seqs_string).hexdigest()
    return _hash

root_dir = os.getcwd()

def resource(file_name):
    return "{0}/unit_test_resources/{1}".format(root_dir, file_name)

seq_files = ["Mnemiopsis_cds.fa", "Mnemiopsis_cds.gb", "Mnemiopsis_cds.nex", "Mnemiopsis_cds.phy",
             "Mnemiopsis_cds.phyr", "Mnemiopsis_cds.stklm",
             "Mnemiopsis_pep.fa", "Mnemiopsis_pep.gb", "Mnemiopsis_pep.nex", "Mnemiopsis_pep.phy",
             "Mnemiopsis_pep.phyr", "Mnemiopsis_pep.stklm"]

@pytest.mark.parametrize("seq_file", seq_files)
def test_instantiate_seqbuddy_from_file(seq_file):
    assert type(Sb.SeqBuddy(resource(seq_file))) == Sb.SeqBuddy

@pytest.mark.parametrize("seq_file", seq_files)
def test_instantiate_seqbuddy_from_handle(seq_file):
    with open(resource(seq_file), 'r') as ifile:
        assert type(Sb.SeqBuddy(ifile)) == Sb.SeqBuddy


@pytest.mark.parametrize("seq_file", seq_files)
def test_instantiate_seqbuddy_from_raw(seq_file):
    with open(resource(seq_file), 'r') as ifile:
        assert type(Sb.SeqBuddy(ifile.read())) == Sb.SeqBuddy

# Now that we know that all the files are being turned into SeqBuddy objects okay, make them all objects so it doesn't
# need to be done over and over for each subsequent test. Note: take care that objects are not modified in place
sb_objects = [Sb.SeqBuddy(resource(x)) for x in seq_files]
formats = ["fasta", "gb", "nexus", "phylip-relaxed", "stockholm",
           "fasta", "gb", "nexus", "phylip-relaxed", "stockholm"]

# fa gb nex phy phyr stklm

#'-ofa', '--order_features_alphabetically'
hashes = ["25073539df4a982b7f99c72dd280bb8f", "0ca1c910e2d566556e839c3ca29e162d", "cb1169c2dd357771a97a02ae2160935d",
          "d1524a20ef968d53a41957d696bfe7ad", "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
          "c10d136c93f41db280933d5b3468f187", "2cd6ab6a922a6695ba2ffdbfa81f740d", "8b6737fe33058121fd99d2deee2f9a76",
          "40f10dc94d85b32155af7446e6402dea", "b229db9c07ff3e4bc049cea73d3ebe2c", "f35cbc6e929c51481e4ec31e95671638"]

hashes = [(sb_objects[indx], value) for indx, value in enumerate(hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_order_features_alphabetically(seqbuddy, next_hash):
    tester = Sb.order_features_alphabetically(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

