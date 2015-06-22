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
        SeqIO.write(_seqbuddy.records, write_file.path, _seqbuddy.out_format)

    seqs_string = write_file.read().encode()
    _hash = md5(seqs_string).hexdigest()
    return _hash

root_dir = os.getcwd()

def resource(file_name):
    return "{0}/unit_test_resources/{1}".format(root_dir, file_name)

seq_files = ["Mnemiopsis_cds.fa", "Mnemiopsis_cds.gb", "Mnemiopsis_cds.nex", "Mnemiopsis_cds.phyr",
             "Mnemiopsis_cds.stklm",
             "Mnemiopsis_pep.fa", "Mnemiopsis_pep.gb", "Mnemiopsis_pep.nex", "Mnemiopsis_pep.phyr",
             "Mnemiopsis_pep.stklm"]

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
hashes = ["25073539df4a982b7f99c72dd280bb8f", "ffa7cb60cb98e50bc4741eed7c88e553", "cb1169c2dd357771a97a02ae2160935d",
          "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b"]

hashes = [(sb_objects[indx], value) for indx, value in enumerate(hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_order_features_alphabetically(seqbuddy, next_hash):
    tester = Sb.order_features_alphabetically(seqbuddy)
    if seqbuddy.out_format == "phylip-relaxed":
        tester.write("troubleshooting.txt")
    assert seqs_to_hash(tester) == next_hash


"""
@pytest.mark.parametrize("seqbuddy,hash", seq_files)
def test_order_features_alphabetically_reverse(fasta_cds, fasta_pep, gb_cds, gb_pep):
    tester = Sb.order_features_alphabetically(gb_cds)
    assert seqs_to_hash(tester, "gb") == "0ca1c910e2d566556e839c3ca29e162d"

    tester = sb.order_features_alphabetically(gb_cds, True)
    assert seqs_to_hash(tester, "gb") == "71bc8a8d19e5c7a6c54db4196f632cea"

    tester = sb.order_features_alphabetically(gb_pep)
    assert seqs_to_hash(tester, "gb") == "2cd6ab6a922a6695ba2ffdbfa81f740d"

    tester = sb.order_features_alphabetically(gb_pep, True)
    assert seqs_to_hash(tester, "gb") == "172f886830e8f8880c31417ec1821490"

    tester = sb.order_features_alphabetically(fasta_cds)
    assert seqs_to_hash(tester, "fasta") == "25073539df4a982b7f99c72dd280bb8f"
"""