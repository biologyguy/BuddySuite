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
    import workshop.SeqBuddy as Sb
except ImportError:
    import SeqBuddy as Sb
import MyFuncs

write_file = MyFuncs.TempFile()


def seqs_to_hash(_seqbuddy, mode='hash'):
    if _seqbuddy.out_format in ["gb", "genbank"]:
            for _rec in _seqbuddy.records:
                try:
                    if re.search("(\. )+", _rec.annotations['organism']):
                        _rec.annotations['organism'] = "."
                except KeyError:
                    pass

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


seq_files = ["Mnemiopsis_cds.fa", "Mnemiopsis_cds.gb", "Mnemiopsis_cds.nex",
             "Mnemiopsis_cds.phy", "Mnemiopsis_cds.phyr", "Mnemiopsis_cds.stklm",
             "Mnemiopsis_pep.fa", "Mnemiopsis_pep.gb", "Mnemiopsis_pep.nex",
             "Mnemiopsis_pep.phy", "Mnemiopsis_pep.phyr", "Mnemiopsis_pep.stklm",
             "ambiguous_dna.fa", "ambiguous_rna.fa"]


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

@pytest.mark.parametrize("seq_file", seq_files)
def test_instantiate_seqbuddy_from_seqbuddy(seq_file):
    input_buddy = Sb.SeqBuddy(resource(seq_file))
    tester = Sb.SeqBuddy(input_buddy)
    assert seqs_to_hash(input_buddy) == seqs_to_hash(tester)


def test_alpha_arg_dna():
    tester = Sb.SeqBuddy(resource(seq_files[0]), _alpha='dna')
    assert tester.alpha is IUPAC.ambiguous_dna

def test_alpha_arg_rna():
    tester = Sb.SeqBuddy(resource(seq_files[0]), _alpha='rna')
    assert tester.alpha is IUPAC.ambiguous_rna

def test_alpha_arg_prot():
    tester = Sb.SeqBuddy(resource(seq_files[6]), _alpha='prot')
    assert tester.alpha is IUPAC.protein

def test_alpha_arg_guess():
    tester = Sb.SeqBuddy(resource(seq_files[0]), _alpha='dgasldfkjhgaljhetlfdjkfg')
    assert tester.alpha is IUPAC.ambiguous_dna

def test_to_string():
    tester = deepcopy(sb_objects[0])
    assert seqs_to_hash(tester) == md5(str(tester).encode()).hexdigest()

def test_to_dict():
    # still figuring out how to do this one
    assert 0



# Now that we know that all the files are being turned into SeqBuddy objects okay, make them all objects so it doesn't
# need to be done over and over for each subsequent test.
sb_objects = [Sb.SeqBuddy(resource(x)) for x in seq_files]
# ######################  'sh', '--shuffle' ###################### #
@pytest.mark.parametrize("seqbuddy", [deepcopy(x) for x in sb_objects])
def test_shuffle(seqbuddy):
    tester = Sb.shuffle(deepcopy(seqbuddy))
    assert seqs_to_hash(seqbuddy) != seqs_to_hash(tester)
    assert seqs_to_hash(Sb.order_ids(tester)) == seqs_to_hash(tester)

# ######################  'rs', '--raw_seq' ###################### #
hashes = ["6f0ff2d43706380d92817e644e5b78a5", "5d00d481e586e287f32d2d29916374ca", "6f0ff2d43706380d92817e644e5b78a5",
          "cda59127d6598f44982a2d1875064bb1", "6f0ff2d43706380d92817e644e5b78a5", "6f0ff2d43706380d92817e644e5b78a5",
          "cdfe71aefecc62c5f5f2f45e9800922c", "4dd913ee3f73ba4bb5dc90d612d8447f", "cdfe71aefecc62c5f5f2f45e9800922c",
          "3f48f81ab579a389947641f36889901a", "cdfe71aefecc62c5f5f2f45e9800922c", "cdfe71aefecc62c5f5f2f45e9800922c"]
hashes = [(sb_objects[indx], value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_raw_seq(seqbuddy, next_hash):
    tester = Sb.raw_seq(seqbuddy)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash

# ######################  'uc', '--uppercase'  and 'lc', '--lowercase' ###################### #
uc_hashes = ["25073539df4a982b7f99c72dd280bb8f", "2e02a8e079267bd9add3c39f759b252c", "52e74a09c305d031fc5263d1751e265d",
             "7117732590f776836cbabdda05f9a982", "3d17ebd1f6edd528a153ea48dc37ce7d", "b82538a4630810c004dc8a4c2d5165ce",
             "c10d136c93f41db280933d5b3468f187", "7a8e25892dada7eb45e48852cbb6b63d", "8b6737fe33058121fd99d2deee2f9a76",
             "40f10dc94d85b32155af7446e6402dea", "b229db9c07ff3e4bc049cea73d3ebe2c", "f35cbc6e929c51481e4ec31e95671638"]

lc_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "2e02a8e079267bd9add3c39f759b252c", "cb1169c2dd357771a97a02ae2160935d",
             "d1524a20ef968d53a41957d696bfe7ad", "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
             "14227e77440e75dd3fbec477f6fd8bdc", "7a8e25892dada7eb45e48852cbb6b63d", "17ff1b919cac899c5f918ce8d71904f6",
             "c934f744c4dac95a7544f9a814c3c22a", "6a3ee818e2711995c95372afe073490b", "c0dce60745515b31a27de1f919083fe9"]

hashes = [(deepcopy(sb_objects[indx]), uc_hash, lc_hashes[indx]) for indx, uc_hash in enumerate(uc_hashes)]


@pytest.mark.parametrize("seqbuddy,uc_hash,lc_hash", hashes)
def test_cases(seqbuddy, uc_hash, lc_hash):  # NOTE: Biopython always writes genbank to spec in lower case
    tester = Sb.uppercase(seqbuddy)
    assert seqs_to_hash(tester) == uc_hash
    tester = Sb.lowercase(tester)
    assert seqs_to_hash(tester) == lc_hash

# ######################  '-ofa', '--order_features_alphabetically' ###################### #
fwd_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "21547b4b35e49fa37e5c5b858808befb",
              "cb1169c2dd357771a97a02ae2160935d", "d1524a20ef968d53a41957d696bfe7ad",
              "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "d0297078b4c480a49b6da5b719310d0e",
              "17ff1b919cac899c5f918ce8d71904f6", "c934f744c4dac95a7544f9a814c3c22a",
              "6a3ee818e2711995c95372afe073490b", "c0dce60745515b31a27de1f919083fe9"]

rev_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "3b718ec3cb794bcb658d900e517110cc",
              "cb1169c2dd357771a97a02ae2160935d", "d1524a20ef968d53a41957d696bfe7ad",
              "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "c6a788d8ea916964605ac2942c459c9b",
              "17ff1b919cac899c5f918ce8d71904f6", "c934f744c4dac95a7544f9a814c3c22a",
              "6a3ee818e2711995c95372afe073490b", "c0dce60745515b31a27de1f919083fe9"]
hashes = [(deepcopy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)  # modifies in place?
def test_order_features_alphabetically(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_features_alphabetically(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_features_alphabetically(seqbuddy, _reverse=True)
    assert seqs_to_hash(tester) == rev_hash

# ######################  'ofp', '--order_features_by_position' ###################### #
fwd_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "2e02a8e079267bd9add3c39f759b252c",
              "cb1169c2dd357771a97a02ae2160935d", "d1524a20ef968d53a41957d696bfe7ad",
              "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "7a8e25892dada7eb45e48852cbb6b63d",
              "17ff1b919cac899c5f918ce8d71904f6", "c934f744c4dac95a7544f9a814c3c22a",
              "6a3ee818e2711995c95372afe073490b", "c0dce60745515b31a27de1f919083fe9"]

rev_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "4345a14fe27570b3c837c30a8cb55ea9",
              "cb1169c2dd357771a97a02ae2160935d", "d1524a20ef968d53a41957d696bfe7ad",
              "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "9e7c2571db1386bba5983365ae235e1b",
              "17ff1b919cac899c5f918ce8d71904f6", "c934f744c4dac95a7544f9a814c3c22a",
              "6a3ee818e2711995c95372afe073490b", "c0dce60745515b31a27de1f919083fe9"]
hashes = [(deepcopy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)  # modifies in place?
def test_order_features_position(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_features_by_position(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_features_by_position(seqbuddy, _reverse=True)
    assert seqs_to_hash(tester) == rev_hash

# ######################  '-mw', '--molecular_weight' ###################### #
def test_molecular_weight():
    # Unambiguous DNA
    tester = Sb.molecular_weight(sb_objects[0])
    assert tester['masses_ds'][0] == 743477.1
    assert tester['masses_ss'][0] == 371242.6
    # Ambiguous DNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("ambiguous_dna.fa")))
    assert tester['masses_ds'][0] == 743477.08
    assert tester['masses_ss'][0] == 371202.59
    # Unambiguous RNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("Mnemiopsis_rna.fa")))
    assert tester['masses_ss'][0] == 387372.6
    # Ambiguous RNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("ambiguous_rna.fa")))
    assert tester['masses_ss'][0] == 387371.6
    # Protein
    tester = Sb.molecular_weight(sb_objects[6])
    assert tester['masses_ss'][0] == 45692.99

# ######################  'cs', '--clean_seq'  ###################### #
def test_clean_seq():
    # Protein
    tester = deepcopy(sb_objects[6])
    tester.alpha = IUPAC.protein
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert seqs_to_hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert seqs_to_hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"

    # DNA
    tester = deepcopy(sb_objects[12])
    tester.alpha = IUPAC.ambiguous_dna
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert seqs_to_hash(tester) == "71b28ad2730a9849f2ba0f70e9e51a9f"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert seqs_to_hash(tester) == "5fd0b78e37c81e0fa727db34a37cc743"

    # RNA
    tester = deepcopy(sb_objects[13])
    tester.alpha = IUPAC.ambiguous_rna
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert seqs_to_hash(tester) == "cdb1b963536d57efc7b7f87d2bf4ad22"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert seqs_to_hash(tester) == "ef61174330f2d9cf5da4f087d12ca201"

    # Alignment formats should raise an error because seq lengths change
    with pytest.raises(ValueError):
        tester = Sb.clean_seq(deepcopy(sb_objects[2]))
        tester.write("temp.del")
        tester = Sb.clean_seq(deepcopy(sb_objects[3]))
        tester.write("temp.del")
        tester = Sb.clean_seq(deepcopy(sb_objects[4]))
        tester.write("temp.del")
        tester = Sb.clean_seq(deepcopy(sb_objects[5]))
        tester.write("temp.del")

# ######################  'dm', '--delete_metadata' ###################### #
hashes = ["aa92396a9bb736ae6a669bdeaee36038", "544ab887248a398d6dd1aab513bae5b1", "cb1169c2dd357771a97a02ae2160935d",
          "d1524a20ef968d53a41957d696bfe7ad", "99d522e8f52e753b4202b1c162197459", "a50943ccd028b6f5fa658178fa8cf54d",
          "bac5dc724b1fee092efccd2845ff2513", "858e8475f7bc6e6a24681083a8635ef9", "17ff1b919cac899c5f918ce8d71904f6",
          "c934f744c4dac95a7544f9a814c3c22a", "6a3ee818e2711995c95372afe073490b", "e224c16f6c27267b5f104c827e78df33"]
hashes = [(deepcopy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_delete_metadata(seqbuddy, next_hash):
    tester = Sb.delete_metadata(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

# ######################  'tr', '--translate' ###################### #
hashes = ["3de7b7be2f2b92cf166b758625a1f316", "c841658e657b4b21b17e4613ac27ea0e", ]
# NOTE: the first 6 sb_objects are DNA.
hashes = [(deepcopy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_translate(seqbuddy, next_hash):
    tester = Sb.translate_cds(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_translate_pep_exception():
    with pytest.raises(TypeError):
        Sb.translate_cds(sb_objects[6])

# ######################  'sfr', '--select_frame' ###################### #
# Only fasta
hashes = ["b831e901d8b6b1ba52bad797bad92d14", "a518e331fb29e8be0fdd5f3f815f5abb", "2cbe39bea876030da6d6bd45e514ae0e"]
frame = [1, 2, 3]
hashes = [(deepcopy(sb_objects[0]), _hash, frame[indx]) for indx, _hash in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash,shift", hashes)
def test_select_frame(seqbuddy, next_hash, shift):
    tester = Sb.select_frame(seqbuddy, shift)
    assert seqs_to_hash(tester) == next_hash


def test_select_frame_pep_exception():
    with pytest.raises(TypeError):  # If protein is input
        Sb.select_frame(sb_objects[6], 2)

# ######################  'tr6', '--translate6frames' ###################### #
# Only fasta and genbank
hashes = ["d5d39ae9212397f491f70d6928047341", "42bb6caf86d2d8be8ab0defabc5af477"]
hashes = [(deepcopy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_translate6frames(seqbuddy, next_hash):
    tester = Sb.translate6frames(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_translate6frames_pep_exception():
    with pytest.raises(TypeError):
        Sb.translate6frames(sb_objects[6])

# ######################  'btr', '--back_translate' ###################### #
# Only fasta and genbank
hashes = ["1b14489a78bfe8255c777138877b9648", "b6bcb4e5104cb202db0ec4c9fc2eaed2",
          "859ecfb88095f51bfaee6a1d1abeb50f", "ba5c286b79a3514fba0b960ff81af25b",
          "952a91a4506afb57f27136aa1f2a8af9", "40c4a3e08c811b6bf3be8bedcb5d65a0"]
organisms = ['human', 'human', 'yeast', 'yeast', 'ecoli', 'ecoli']
hashes = [(deepcopy(sb_objects[sb_obj_indx]), organisms[indx], hashes[indx]) for indx, sb_obj_indx in
          enumerate([6, 7, 6, 7, 6, 7])]


@pytest.mark.parametrize("seqbuddy,_organism,next_hash", hashes)
def test_back_translate(seqbuddy, _organism, next_hash):
    seqbuddy.alpha = IUPAC.protein
    tester = Sb.back_translate(seqbuddy, 'OPTIMIZED', _organism)
    assert seqs_to_hash(tester) == next_hash


def test_back_translate_nucleotide_exception():
    with pytest.raises(TypeError):
        Sb.back_translate(sb_objects[1])

def test_back_translate_bad_mode():
    with pytest.raises(AttributeError):
        Sb.back_translate(deepcopy(sb_objects[6]), 'fgsdjkghjdalgsdf', 'human')

def test_back_translate_bad_organism():
    seqbuddy = deepcopy(sb_objects[6])
    seqbuddy.alpha = IUPAC.protein
    with pytest.raises(AttributeError):
        Sb.back_translate(seqbuddy, 'OPTIMIZED', 'fgsdjkghjdalgsdf')

# ######################  'd2r', '--transcribe' and 'r2d', '--back_transcribe' ###################### #
d2r_hashes = ["d2db9b02485e80323c487c1dd6f1425b", "9ef3a2311a80f05f21b289ff7f401fff",
              "f3bd73151645359af5db50d2bdb6a33d", "1371b536e41e3bca304794512122cf17",
              "866aeaca326891b9ebe5dc9d762cba2c", "45b511f34653e3b984e412182edee3ca"]
r2d_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "2e02a8e079267bd9add3c39f759b252c",
              "cb1169c2dd357771a97a02ae2160935d", "d1524a20ef968d53a41957d696bfe7ad",
              "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b"]

hashes = [(deepcopy(sb_objects[indx]), d2r_hash, r2d_hashes[indx]) for indx, d2r_hash in enumerate(d2r_hashes)]


@pytest.mark.parametrize("seqbuddy,d2r_hash,r2d_hash", hashes)
def test_transcribe(seqbuddy, d2r_hash, r2d_hash):
    tester = Sb.dna2rna(seqbuddy)
    assert seqs_to_hash(tester) == d2r_hash
    tester = Sb.rna2dna(tester)
    assert seqs_to_hash(tester) == r2d_hash


def test_transcribe_pep_exception():  # Asserts that a ValueError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.dna2rna(sb_objects[6])


def test_back_transcribe_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.rna2dna(sb_objects[6])


# ######################  'cmp', '--complement' ###################### #
hashes = ["e4a358ca57aca0bbd220dc6c04c88795", "3366fcc6ead8f1bba4a3650e21db4ec3",
          "365bf5d08657fc553315aa9a7f764286", "10ce87a53aeb5bd4f911380ebf8e7a85",
          "8e5995813da43c7c00e98d15ea466d1a", "5891348e8659290c2355fabd0f3ba4f4"]
hashes = [(deepcopy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_complement(seqbuddy, next_hash):
    tester = Sb.complement(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_complement_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.complement(sb_objects[6])

# ######################  'rc', '--reverse_complement' ###################### #
hashes = ["e77be24b8a7067ed54f06e0db893ce27", "47941614adfcc5bd107f71abef8b3e00", "f549c8dc076f6b3b4cf5a1bc47bf269d",
          "a62edd414978f91f7391a59fc1a72372", "08342be5632619fd1b1251b7ad2b2c84", "0d6b7deda824b4fc42b65cb87e1d4d14"]
hashes = [(deepcopy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_reverse_complement(seqbuddy, next_hash):
    tester = Sb.reverse_complement(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_reverse_complement_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.reverse_complement(sb_objects[6])

# ######################  'li', '--list_ids' ###################### #
# first test that 1 column works for all file types
hashes = ["1c4a395d8aa3496d990c611c3b6c4d0a", "1c4a395d8aa3496d990c611c3b6c4d0a", "1c4a395d8aa3496d990c611c3b6c4d0a",
          "78a9289ab2d508a13c76cf9f5a308cc5", "1c4a395d8aa3496d990c611c3b6c4d0a", "1c4a395d8aa3496d990c611c3b6c4d0a"]
hashes = [(sb_objects[indx], value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_list_ids_one_col(seqbuddy, next_hash):
    tester = Sb.list_ids(seqbuddy, 1)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash

# Now test different numbers of columns
hashes = ["6fcee2c407bc4f7f70e0ae2a7e101761", "1c4a395d8aa3496d990c611c3b6c4d0a", "6fcee2c407bc4f7f70e0ae2a7e101761",
          "bd177e4db7dd772c5c42199b0dff49a5", "6b595a436a38e353a03e36a9af4ba1f9", "c57028374ed3fc474009e890acfb041e"]
columns = [-2, 0, 2, 5, 10, 100]
hashes = [(sb_objects[0], value, columns[indx]) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash,cols", hashes)
def test_list_ids_multi_col(seqbuddy, next_hash, cols):
    tester = Sb.list_ids(seqbuddy, cols)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash

#  are num_seqs and ave_seq_length even worth testing? YES!!!! Just feed in numbers though, no need to calculate hashes

# ######################  'cts', '--concat_seqs' ###################### #
hashes = ["2e46edb78e60a832a473397ebec3d187", "7421c27be7b41aeedea73ff41869ac47",
          "494988ffae2ef3072c1619eca8a0ff3b", "710cad348c5560446daf2c916ff3b3e4",
          "494988ffae2ef3072c1619eca8a0ff3b", "494988ffae2ef3072c1619eca8a0ff3b",
          "46741638cdf7abdf53c55f79738ee620", "8d0bb4e5004fb6a1a0261c30415746b5",
          "2651271d7668081cde8012db4f9a6574", "36526b8e0360e259d8957fa2261cf45a",
          "2651271d7668081cde8012db4f9a6574", "2651271d7668081cde8012db4f9a6574"]
hashes = [(deepcopy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_concat_seqs(seqbuddy, next_hash):
    tester = Sb.concat_seqs(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

# ToDo: Test the _clean parameter

# ######################  'fd2p', '--map_features_dna2prot' ###################### #
# Map the genbank DNA file to all protein files, and the fasta DNA file to fasta protein
hashes = ["5216ef85afec36d5282578458a41169a", "a8f7c129cf57a746c20198bf0a6b9cf4", "0deeea532d6dcbc0486e9b74d0d6aca8",
          "d595fabb157d5c996357b6a7058af4e8", "bb06e94456f99efc2068f5a52f0e0462", "a287e0054df7f5df76e792e0e0ab6756"]
prot_indx = [6, 7, 8, 9, 10, 11]
hashes = [(deepcopy(sb_objects[1]), deepcopy(sb_objects[prot_indx[indx]]), value) for indx, value in enumerate(hashes)]
hashes.append((deepcopy(sb_objects[0]), deepcopy(sb_objects[6]), "854566b485af0f277294bbfb15f7dd0a"))


@pytest.mark.parametrize("_dna,_prot,next_hash", hashes)
def test_map_features_dna2prot(_dna, _prot, next_hash):
    _prot.alpha = IUPAC.protein
    _dna.alpha = IUPAC.ambiguous_dna
    tester = Sb.map_features_dna2prot(_dna, _prot)
    assert seqs_to_hash(tester) == next_hash


# ######################  'fp2d', '--map_features_prot2dna' ###################### #
# Map the genbank protein file to all dna files, and the fasta protein file to fasta DNA
hashes = ["3ebc92ca11505489cab2453d2ebdfcf2", "6b4dd3fc66cb7419acaf064b589f4dd1",
          "8d403750ef83d60e31de0dee79a8f5d1", "74c6c4b5531c41f55f7349ed6c6b2f43",
          "9133ab0becbec95ce7ed31e02dc17ef5", "3ebc92ca11505489cab2453d2ebdfcf2"]
dna_indx = [0, 1, 2, 3, 4, 5]
hashes = [(deepcopy(sb_objects[7]), deepcopy(sb_objects[dna_indx[indx]]), value) for indx, value in enumerate(hashes)]
hashes.append((deepcopy(sb_objects[6]), deepcopy(sb_objects[0]), "720f36544f9c11855ac2673e63282f89"))


@pytest.mark.parametrize("_prot,_dna,next_hash", hashes)
def test_map_features_prot2dna(_prot, _dna, next_hash):
    _prot.alpha = IUPAC.protein
    _dna.alpha = IUPAC.ambiguous_dna
    tester = Sb.map_features_prot2dna(_prot, _dna)
    assert seqs_to_hash(tester) == next_hash


# ######################  'ri', '--rename_ids' ###################### #
hashes = ["59bea136d93d30e3f11fd39d73a9adff", "78c73f97117bd937fd5cf52f4bd6c26e", "243024bfd2f686e6a6e0ef65aa963494",
          "83f10d1be7a5ba4d363eb406c1c84ac7", "973e3d7138b78db2bb3abda8a9323226", "4289f03afb6c9f8a8b0d8a75bb60a2ce"]
hashes = [(deepcopy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_rename_ids(seqbuddy, next_hash):
    tester = Sb.rename(seqbuddy, 'Panx', 'Test', 0)
    assert seqs_to_hash(tester) == next_hash

# ######################  'cf', '--combine_features' ###################### #
# Only one test at the moment. Maybe add more later.
dummy_feats = Sb.SeqBuddy(resource("Mnemiopsis_cds_dummy_features.gb"))


def test_combine_features():
    tester = Sb.combine_features(dummy_feats, sb_objects[1])
    assert seqs_to_hash(tester) == "11ae234528ba6871035c1a3641d71736"

# ######################  'oi', '--order_ids' ###################### #
fwd_hashes = ["ccc59629741421fb717b9b2403614c62", "2f9bc0dd9d79fd8160a621280be0b0aa",
              "60bbc6306cbb4eb903b1212718bb4592", "4188a065adb5b8e80acd3073afc1c7f9",
              "433520b63864022d82973f493dbf804b", "4078182a81382b815528fdd5c158fbec"]
rev_hashes = ["503a71fc2e8d143361cbe8f4611527fd", "dd269961d4d5301d1bf87e0093568851",
              "82fea6e3d3615ac75ec5022abce255da", "9d0910f3d303297283bace2718f60d61",
              "8af06c3523a1bf7cde4fc2b8c64a388c", "3b83a3c73a6cdded6635ffa10c4a16e1"]

hashes = [(deepcopy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)
def test_order_ids(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_ids(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_ids(seqbuddy, _reverse=True)
    assert seqs_to_hash(tester) == rev_hash

# ######################  'er', '--extract_range' ###################### #
hashes = ["201235ed91ad0ed9a7021136487fed94", "3e791c6a6683516aff9572c24f38f0b3", "4063ab66ced2fafb080ceba88965d2bb",
          "0c857970ebef51b4bbd9c7b3229d7088", "e0e256cebd6ead99ed3a2a20b7417ba1", "d724df01ae688bfac4c6dfdc90027440",
          "904a188282f19599a78a9d7af4169de6", "b8413624b9e684a14fc9f398a62e3965", "6a27222d8f60ee8496cbe0c41648a116",
          "9ecc1d83eff77c61284869b088c833a1", "9c85530cd3e3aa628b0e8297c0c9f977", "38d571c9681b4fa420e3d8b54c507f9c"]
hashes = [(deepcopy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_extract_range(seqbuddy, next_hash):
    tester = Sb.extract_range(seqbuddy, 50, 300)
    assert seqs_to_hash(tester) == next_hash


def test_extract_range_end_less_than_start():
    with pytest.raises(ValueError):
        Sb.extract_range(sb_objects[0], 500, 50)

# ######################  'ns', '--num_seqs' ###################### #
seq_counts = [(sb_objects[0], 13), (sb_objects[1], 13), (sb_objects[2], 13), (sb_objects[3], 8),
              (sb_objects[4], 13), (sb_objects[5], 13), (sb_objects[6], 13), (sb_objects[9], 8)]


@pytest.mark.parametrize("seqbuddy, num", seq_counts)
def test_num_seqs(seqbuddy, num):
    assert Sb.num_seqs(seqbuddy) == num


def test_empty_file():
    with pytest.raises(SystemExit):
        Sb.SeqBuddy(resource("blank.fa"))


# ######################  'asl', '--ave_seq_length' ###################### #
@pytest.mark.parametrize("seqbuddy", [Sb.SeqBuddy(resource(x)) for x in seq_files[0:6] if ".phy" not in x])
def test_ave_seq_length_dna(seqbuddy):
    print(seqbuddy.in_format)
    assert round(Sb.ave_seq_length(seqbuddy, _clean=True), 2) == 1285.15

@pytest.mark.parametrize("seqbuddy", [Sb.SeqBuddy(resource(x)) for x in seq_files[6:12] if ".phy" not in x])
def test_ave_seq_length_pep(seqbuddy):
    assert round(Sb.ave_seq_length(seqbuddy, _clean=True), 2) == 427.38


# ######################  'pr', '--pull_recs' ###################### #
pr_hashes = ["5b4154c2662b66d18776cdff5af89fc0", "e196fdc5765ba2c47f97807bafb6768c", "bc7dbc612bc8139eba58bf896b7eaf2f",
             "e33828908fa836f832ee915957823039", "e33828908fa836f832ee915957823039", "b006b40ff17ba739929448ae2f9133a6"]
pr_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(pr_hashes)]
@pytest.mark.parametrize("seqbuddy, next_hash", pr_hashes)
def test_pull_recs(seqbuddy, next_hash):
    tester = Sb.pull_recs(seqbuddy, 'α2')
    assert seqs_to_hash(tester) == next_hash


# ######################  'dr', '--delete_records' ###################### #
dr_hashes = ["54bdb42423b1d331acea18218101e5fc", "e2c03f1fa21fd27b2ff55f7f721a1a99", "6bc8a9409b1ef38e4f6f12121368883e",
             "a63320b6a679f97368329396e3b72bdd", "fe927a713e92b7de0c3a63996a4ce7c8", "4c97144c5337f8a40c4fd494e622bf0d"]
dr_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(dr_hashes)]
@pytest.mark.parametrize("seqbuddy, next_hash", dr_hashes)
def test_pull_recs(seqbuddy, next_hash):
    tester = Sb.delete_records(seqbuddy, 'α2')
    assert seqs_to_hash(tester) == next_hash


# ######################  'dr', '--delete_records' ###################### #
@pytest.mark.parametrize("seqbuddy", [Sb.SeqBuddy(resource(x)) for x in [seq_files[6], seq_files[7], seq_files[8],
                                                                         seq_files[10], seq_files[11]]])
def test_isoelectric_point(seqbuddy):
    output = Sb.isoelectric_point(Sb.clean_seq(seqbuddy))
    assert output[0][1] == 6.0117797852

def test_isoelectric_point_type_error():
    with pytest.raises(TypeError):
        Sb.isoelectric_point(Sb.SeqBuddy(resource("/Mnemiopsis_cds.fa")))


# ######################  'GuessError' ###################### #
def test_guesserror_raw_seq():
    with pytest.raises(Sb.GuessError):
        Sb.SeqBuddy("JSKHGLHGLSDKFLSDYUIGJVSBDVHJSDKGIUSUEWUIOIFUBCVVVBVNNJS{QF(*&#@$(*@#@*(*(%")

def test_guesserror_infile():
    with pytest.raises(Sb.GuessError):
        Sb.SeqBuddy(resource("gibberish.fa"))

def test_no__input():
    with pytest.raises(TypeError):
        Sb.SeqBuddy()

