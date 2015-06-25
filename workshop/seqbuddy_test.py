#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

import pytest
from hashlib import md5
from Bio import SeqIO
import os
try:
    import workshop.SeqBuddy as Sb
except ImportError:
    import SeqBuddy as Sb
import MyFuncs

write_file = MyFuncs.TempFile()

def seqs_to_hash(_seqbuddy, mode='hash'):
    if _seqbuddy.out_format == "phylipi":
        write_file.write(Sb.phylipi(_seqbuddy, "relaxed"))
    elif _seqbuddy.out_format == "phylipis":
        write_file.write(Sb.phylipi(_seqbuddy, "strict"))
    else:
        _seqbuddy.write(write_file.path)

    seqs_string = "{0}\n".format(write_file.read().strip())

    if mode != "hash":
        return seqs_string

    _hash = md5(seqs_string.encode()).hexdigest()
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
hashes = ["25073539df4a982b7f99c72dd280bb8f", "ffa7cb60cb98e50bc4741eed7c88e553", "cb1169c2dd357771a97a02ae2160935d",
          "d1524a20ef968d53a41957d696bfe7ad", "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
          "c10d136c93f41db280933d5b3468f187", "f6b3c090ab6ac147cf6d87881a8ce5dc", "8b6737fe33058121fd99d2deee2f9a76",
          "40f10dc94d85b32155af7446e6402dea", "b229db9c07ff3e4bc049cea73d3ebe2c", "f35cbc6e929c51481e4ec31e95671638"]

hashes = [(sb_objects[indx], value) for indx, value in enumerate(hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_order_features_alphabetically(seqbuddy, next_hash):
    tester = Sb.order_features_alphabetically(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

#'-mw', '--molecular_weight'
mw_files = ["mw_test_pep.fa", "mw_test_cds_a.fa", "mw_test_cds_u.fa", "mw_test_rna_a.fa", "mw_test_rna_u.fa"]
mw_formats = ["protein", "dna", "dna", "rna", "rna"]
mw_objects = [(Sb.SeqBuddy(resource(value), "fasta", "fasta", mw_formats[indx])) for indx, value in enumerate(mw_files)]
expected_mw = [[2505.75, None], [5022.19, 10044.28], [3168.0, 6335.9], [4973.0, None], [3405.0, None]]
expected_mw = [(mw_objects[indx], value) for indx, value in enumerate(expected_mw)]
@pytest.mark.parametrize("seqbuddy,next_mw", expected_mw)
def test_molecular_weight(seqbuddy, next_mw):
    tester = Sb.molecular_weight(seqbuddy)
    masses_ss = tester['masses_ss']
    masses_ds = tester['masses_ds']
    print(tester)
    print(seqbuddy)
    print(seqbuddy.records)
    print(masses_ss)
    assert masses_ss[0] == next_mw[0]
    if len(masses_ds) != 0:
        assert masses_ds[0] == next_mw[1]

#'cs', '--clean_seq'
cs_files = ["cs_test_pep.fa", "cs_test_cds_a.fa", "cs_test_cds_u.fa"]
cs_formats = ["protein","dna","dna"]
cs_objects = [(Sb.SeqBuddy(resource(value), "fasta", "fasta", cs_formats[indx])) for indx, value in enumerate(cs_files)]
cs_hashes = ['91a49ad673238e5454227f5d0f8c29ba', "11f39968ab6e6841960b7dd7a4ca124b", "d71a59aaa3716cdc4b3934a2c04d932a"]
cs_hashes = [(cs_objects[indx], value) for indx, value in enumerate(cs_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash",cs_hashes)
def test_clean_seq(seqbuddy, next_hash):
    tester = Sb.clean_seq(seqbuddy)
    print(tester.records[0].seq)
    assert seqs_to_hash(tester) == next_hash

if __name__ == '__main__':
    debug = Sb.order_features_alphabetically(sb_objects[1])
    print(seqs_to_hash(debug, "string"))