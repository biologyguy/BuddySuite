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
    assert masses_ss[0] == next_mw[0]
    if len(masses_ds) != 0:
        assert masses_ds[0] == next_mw[1]

#'cs', '--clean_seq'
cs_files = ["cs_test_pep.fa", "cs_test_cds_a.fa", "cs_test_cds_u.fa"]
cs_formats = ["protein","dna","dna"]
cs_objects = [(Sb.SeqBuddy(resource(value), "fasta", "fasta", cs_formats[indx])) for indx, value in enumerate(cs_files)]
cs_hashes = ['9289d387b1c8f990b44a9cb15e12443b', "8e161d5e4115bf483f5196adf7de88f0", "2e873cee6f807fe17cb0ff9437d698fb"]
cs_hashes = [(cs_objects[indx], value) for indx, value in enumerate(cs_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash",cs_hashes)
def test_clean_seq(seqbuddy, next_hash):
    tester = Sb.clean_seq(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

lc_files = ["lower_Mnemiopsis_cds.fa", "lower_Mnemiopsis_cds.gb", "lower_Mnemiopsis_cds.nex",
            "lower_Mnemiopsis_cds.phyr", "lower_Mnemiopsis_cds.stklm", "lower_Mnemiopsis_pep.fa",
            "lower_Mnemiopsis_pep.gb", "lower_Mnemiopsis_pep.nex", "lower_Mnemiopsis_pep.phyr",
            "lower_Mnemiopsis_pep.stklm"]
lc_objects = [Sb.SeqBuddy(resource(file)) for file in lc_files]
lc_hashes = ["25073539df4a982b7f99c72dd280bb8f", "d41d8cd98f00b204e9800998ecf8427e", "52e74a09c305d031fc5263d1751e265d",
             "3d17ebd1f6edd528a153ea48dc37ce7d", "b82538a4630810c004dc8a4c2d5165ce", "c10d136c93f41db280933d5b3468f187",
             "d41d8cd98f00b204e9800998ecf8427e", "8b6737fe33058121fd99d2deee2f9a76", "b229db9c07ff3e4bc049cea73d3ebe2c",
             "f35cbc6e929c51481e4ec31e95671638"]
lc_hashes = [(lc_objects[indx], value) for indx, value in enumerate(lc_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", lc_hashes)
def test_uppercase(seqbuddy, next_hash): # genbank should fail right now
    tester = Sb.uppercase(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

uc_files = ["upper_Mnemiopsis_cds.fa", "upper_Mnemiopsis_cds.gb", "upper_Mnemiopsis_cds.nex",
            "upper_Mnemiopsis_cds.phyr", "upper_Mnemiopsis_cds.stklm", "upper_Mnemiopsis_pep.fa",
            "upper_Mnemiopsis_pep.gb", "upper_Mnemiopsis_pep.nex", "upper_Mnemiopsis_pep.phyr",
            "upper_Mnemiopsis_pep.stklm"]
uc_objects = [Sb.SeqBuddy(resource(file)) for file in lc_files]
uc_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "4ccc2d108eb01614351bcbeb21932ceb", "cb1169c2dd357771a97a02ae2160935d",
             "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b", "14227e77440e75dd3fbec477f6fd8bdc",
             "0e575609db51ae25d4d41333c56f5661", "17ff1b919cac899c5f918ce8d71904f6", "6a3ee818e2711995c95372afe073490b",
             "c0dce60745515b31a27de1f919083fe9"]
uc_hashes = [(uc_objects[indx], value) for indx, value in enumerate(uc_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", uc_hashes)
def test_lowercase(seqbuddy, next_hash): # not sure why genbank is failing here
    tester = Sb.lowercase(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

if __name__ == '__main__':
    debug = Sb.order_features_alphabetically(sb_objects[1])
    print(seqs_to_hash(debug, "string"))
