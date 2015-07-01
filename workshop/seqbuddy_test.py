#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

import pytest
from hashlib import md5
from Bio import SeqIO
import os
import subprocess
import copy

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


seq_files = ["Mnemiopsis/Mnemiopsis_cds.fa", "Mnemiopsis/Mnemiopsis_cds.gb", "Mnemiopsis/Mnemiopsis_cds.nex",
             "Mnemiopsis/Mnemiopsis_cds.phy", "Mnemiopsis/Mnemiopsis_cds.phyr", "Mnemiopsis/Mnemiopsis_cds.stklm",
             "Mnemiopsis/Mnemiopsis_pep.fa", "Mnemiopsis/Mnemiopsis_pep.gb", "Mnemiopsis/Mnemiopsis_pep.nex",
             "Mnemiopsis/Mnemiopsis_pep.phy", "Mnemiopsis/Mnemiopsis_pep.phyr", "Mnemiopsis/Mnemiopsis_pep.stklm"]


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

# '-ofa', '--order_features_alphabetically'
hashes = ["25073539df4a982b7f99c72dd280bb8f", "ffa7cb60cb98e50bc4741eed7c88e553", "cb1169c2dd357771a97a02ae2160935d",
          "d1524a20ef968d53a41957d696bfe7ad", "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
          "c10d136c93f41db280933d5b3468f187", "f6b3c090ab6ac147cf6d87881a8ce5dc", "8b6737fe33058121fd99d2deee2f9a76",
          "40f10dc94d85b32155af7446e6402dea", "b229db9c07ff3e4bc049cea73d3ebe2c", "f35cbc6e929c51481e4ec31e95671638"]

hashes = [(sb_objects[indx], value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)  # might modify in place
def test_order_features_alphabetically(seqbuddy, next_hash):
    tester = Sb.order_features_alphabetically(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

# '-mw', '--molecular_weight'
mw_files = ["mw/mw_test_pep.fa", "mw/mw_test_cds_a.fa", "mw/mw_test_cds_u.fa", "mw/mw_test_rna_a.fa",
            "mw/mw_test_rna_u.fa"]
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

# 'cs', '--clean_seq'
cs_files = ["cs/cs_test_pep.fa", "cs/cs_test_cds_a.fa", "cs/cs_test_cds_u.fa"]
cs_formats = ["protein", "dna", "dna"]
cs_objects = [(Sb.SeqBuddy(resource(value), "fasta", "fasta", cs_formats[indx])) for indx, value in enumerate(cs_files)]
cs_hashes = ['9289d387b1c8f990b44a9cb15e12443b', "8e161d5e4115bf483f5196adf7de88f0", "2e873cee6f807fe17cb0ff9437d698fb"]
cs_hashes = [(cs_objects[indx], value) for indx, value in enumerate(cs_hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", cs_hashes)
def test_clean_seq(seqbuddy, next_hash):
    tester = Sb.clean_seq(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

# 'uc', '--uppercase'
lc_files = ["lower/lower_Mnemiopsis_cds.fa", "lower/lower_Mnemiopsis_cds.gb", "lower/lower_Mnemiopsis_cds.nex",
            "lower/lower_Mnemiopsis_cds.phyr", "lower/lower_Mnemiopsis_cds.stklm", "lower/lower_Mnemiopsis_pep.fa",
            "lower/lower_Mnemiopsis_pep.gb", "lower/lower_Mnemiopsis_pep.nex", "lower/lower_Mnemiopsis_pep.phyr",
            "lower/lower_Mnemiopsis_pep.stklm"]
lc_objects = [Sb.SeqBuddy(resource(file)) for file in lc_files]
lc_hashes = ["25073539df4a982b7f99c72dd280bb8f", "d41d8cd98f00b204e9800998ecf8427e", "52e74a09c305d031fc5263d1751e265d",
             "3d17ebd1f6edd528a153ea48dc37ce7d", "b82538a4630810c004dc8a4c2d5165ce", "c10d136c93f41db280933d5b3468f187",
             "d41d8cd98f00b204e9800998ecf8427e", "8b6737fe33058121fd99d2deee2f9a76", "b229db9c07ff3e4bc049cea73d3ebe2c",
             "f35cbc6e929c51481e4ec31e95671638"]
lc_hashes = [(lc_objects[indx], value) for indx, value in enumerate(lc_hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", lc_hashes)
def test_uppercase(seqbuddy, next_hash):  # genbank should fail right now
    tester = Sb.uppercase(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

# 'lc', '--lowercase'
uc_files = ["upper/upper_Mnemiopsis_cds.fa", "upper/upper_Mnemiopsis_cds.gb", "upper/upper_Mnemiopsis_cds.nex",
            "upper/upper_Mnemiopsis_cds.phyr", "upper/upper_Mnemiopsis_cds.stklm", "upper/upper_Mnemiopsis_pep.fa",
            "upper/upper_Mnemiopsis_pep.gb", "upper/upper_Mnemiopsis_pep.nex", "upper/upper_Mnemiopsis_pep.phyr",
            "upper/upper_Mnemiopsis_pep.stklm"]
uc_objects = [Sb.SeqBuddy(resource(file)) for file in uc_files]
uc_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "4ccc2d108eb01614351bcbeb21932ceb", "cb1169c2dd357771a97a02ae2160935d",
             "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b", "14227e77440e75dd3fbec477f6fd8bdc",
             "0e575609db51ae25d4d41333c56f5661", "17ff1b919cac899c5f918ce8d71904f6", "6a3ee818e2711995c95372afe073490b",
             "c0dce60745515b31a27de1f919083fe9"]
uc_hashes = [(uc_objects[indx], value) for indx, value in enumerate(uc_hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", uc_hashes)
def test_lowercase(seqbuddy, next_hash):
    tester = Sb.lowercase(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

# 'dm', '--delete_metadata'

dm_hashes = ["8b98dc863c2483ee2fa813df983f9941", "544ab887248a398d6dd1aab513bae5b1", "cb1169c2dd357771a97a02ae2160935d",
             "d1524a20ef968d53a41957d696bfe7ad", "99d522e8f52e753b4202b1c162197459", "a50943ccd028b6f5fa658178fa8cf54d",
             "31b3aa8d247175ae31d5f17b146daf46", "858e8475f7bc6e6a24681083a8635ef9", "8b6737fe33058121fd99d2deee2f9a76",
             "40f10dc94d85b32155af7446e6402dea", "b229db9c07ff3e4bc049cea73d3ebe2c", "45ee83084fbd886d5ab1875f8af4560e"]
dm_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(dm_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", dm_hashes)
def test_delete_metadata(seqbuddy, next_hash):
    tester = Sb.delete_metadata(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

# 'rs', '--raw_seq'

seq_path = root_dir+"/unit_test_resources/"
rs_hashes = ["5d00d481e586e287f32d2d29916374ca", "5d00d481e586e287f32d2d29916374ca", "5d00d481e586e287f32d2d29916374ca",
             "2602037afcfaa467b77db42a0f25a9c8", "5d00d481e586e287f32d2d29916374ca", "5d00d481e586e287f32d2d29916374ca",
             "4dd913ee3f73ba4bb5dc90d612d8447f", "4dd913ee3f73ba4bb5dc90d612d8447f", "4dd913ee3f73ba4bb5dc90d612d8447f",
             "215c09fec462c202989b416ebc47cccc", "4dd913ee3f73ba4bb5dc90d612d8447f", "4dd913ee3f73ba4bb5dc90d612d8447f"]
rs_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(rs_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", rs_hashes)
def test_raw_seq(seqbuddy, next_hash):
    tester = Sb.raw_seq(seqbuddy)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash

# 'tr', '--translate'

tr_hashes = ["fasta_placeholder", "genbank_placeholder"]  # need to replace hashes later
tr_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(tr_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", tr_hashes)
def test_translate(seqbuddy,next_hash):
    tester = Sb.translate_cds(seqbuddy)
    tester = Sb.order_features_alphabetically(tester)
    assert seqs_to_hash(tester) == next_hash

def test_translate_pep_exception():
    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis/Mnemiopsis_pep.fa"))
    with pytest.raises(ValueError):
        Sb.translate_cds(seqbuddy)

# 'sfr', '--select_frame'

sfr_hashes = ["25073539df4a982b7f99c72dd280bb8f", "91ade6dd5aa97dfb14826a44f0497e17",
              "ba029cbc1a52fccddb4304e9b6a4c3db"]
sfr_hashes = [(value, indx+1) for indx, value in enumerate(sfr_hashes)]
@pytest.mark.parametrize("next_hash, shift", sfr_hashes)  # only tests fasta, shouldn't matter
def test_select_frame(next_hash, shift):
    sfr_buddy = Sb.SeqBuddy(resource("Mnemiopsis/Mnemiopsis_cds.fa"))
    tester = Sb.select_frame(sfr_buddy, shift)
    assert seqs_to_hash(tester) == next_hash

def test_select_frame_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein into -sfr
    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis/Mnemiopsis_pep.fa"))
    with pytest.raises(ValueError):
        Sb.select_frame(seqbuddy, 3)

# 'tr6', '--translate6frames' (inconsistent with underscore naming convention)

tr6_hashes = ["fasta_placeholder", "genbank_placeholder"]  # need to replace hashes later
tr6_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(tr6_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", tr6_hashes)
def test_translate6frames(seqbuddy,next_hash):
    tester = Sb.translate6frames(seqbuddy)
    tester = Sb.order_features_alphabetically(tester)
    assert seqs_to_hash(tester) == next_hash

def test_translate6frames_pep_exception():
    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis/Mnemiopsis_pep.fa"))
    with pytest.raises(ValueError):
        Sb.translate6frames(seqbuddy)

# 'btr', '--back_translate'

btr_hashes = ["fasta_human_placeholder", "fasta_yeast_placeholder", "fasta_ecoli_placeholder",
              "genbank_human_placeholder", "genbank_yeast_placeholder", "genbank_ecoli_placeholder"]
btr_hashes = [(Sb.SeqBuddy(resource(seq_files[indx+6])), value) for indx, value in enumerate(btr_hashes)]
btr_organisms = ['human', 'yeast', 'ecoli']
btr_objects = []
for indx, md5_hash in enumerate(btr_hashes):
    for organism in btr_organisms:
        obj = Sb.SeqBuddy(resource(seq_files[indx+6]))
        btr_objects.append((obj, organism, md5_hash))

@pytest.mark.parametrize("seqbuddy,_organism,next_hash", btr_objects)
def test_back_translate(seqbuddy, _organism, next_hash):
    tester = Sb.back_translate(seqbuddy, 'OPTIMIZED', _organism)
    tester = Sb.order_features_alphabetically(tester)
    assert seqs_to_hash(tester) == next_hash

def test_back_translate_nucleotide_exception():
    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis/Mnemiopsis_cds.fa"))
    with pytest.raises(ValueError):
        Sb.back_translate(seqbuddy)

# 'd2r', '--transcribe'

d2r_objects = [Sb.SeqBuddy(resource(x)) for x in seq_files[:6]]
d2r_hashes = ["013ebe2bc7d83c44f58344b865e1f55b", "7464605c739e23d34ce08d3ef51e6a0a",
              "f3bd73151645359af5db50d2bdb6a33d", "1371b536e41e3bca304794512122cf17",
              "866aeaca326891b9ebe5dc9d762cba2c", "45b511f34653e3b984e412182edee3ca"]
d2r_hashes = [(d2r_objects[indx],value) for indx, value in enumerate(d2r_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", d2r_hashes)  # might modify in place
def test_transcribe(seqbuddy, next_hash):
    tester = Sb.dna2rna(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

def test_transcribe_pep_exception():  # Asserts that a ValueError will be thrown if user inputs protein into -sfr
    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis/Mnemiopsis_pep.fa"))
    with pytest.raises(ValueError):
        Sb.dna2rna(seqbuddy)

# 'r2d', '--back_transcribe'

r2d_hashes = ["25073539df4a982b7f99c72dd280bb8f", "57bd873f08cd69f00333eac94c120000",
              "cb1169c2dd357771a97a02ae2160935d", "d1524a20ef968d53a41957d696bfe7ad",
              "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b"]
r2d_objects = [Sb.SeqBuddy(resource(x)) for x in seq_files[:6]]
r2d_hashes = [(Sb.dna2rna(r2d_objects[indx]), value) for indx, value in enumerate(r2d_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", r2d_hashes)  # once again fails for genbank???
def test_back_transcribe(seqbuddy, next_hash):
    tester = Sb.rna2dna(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

#TODO make errors consistent
def test_back_transcribe_pep_exception():  # Asserts that a ValueError will be thrown if user inputs protein into -sfr
    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis/Mnemiopsis_pep.fa"))
    with pytest.raises(ValueError):
        Sb.rna2dna(seqbuddy)

# 'cmp', '--complement'

cmp_objects = [Sb.SeqBuddy(resource(x)) for x in seq_files[:6]]
cmp_hashes = ["cd4a98936eef4ebb05f58a8c614a0f7c", "366e0a2b28623c51591047768e8ddb08",
              "365bf5d08657fc553315aa9a7f764286", "10ce87a53aeb5bd4f911380ebf8e7a85",
              "8e5995813da43c7c00e98d15ea466d1a", "5891348e8659290c2355fabd0f3ba4f4"]
cmp_hashes = [(cmp_objects[indx], value) for indx, value in enumerate(cmp_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", cmp_hashes)
def test_complement(seqbuddy, next_hash):
    tester = Sb.complement(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

def test_complement_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein into -sfr
    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis/Mnemiopsis_pep.fa"))
    with pytest.raises(ValueError):
        Sb.complement(seqbuddy)

# 'rc', '--reverse_complement'

rc_objects = [Sb.SeqBuddy(resource(x)) for x in seq_files[:6]]
rc_hashes = ["cb3ad86bdaed7dd0bcfbca0a46cdfbf9", "d3c70b16443606ef4f51296f67420231", "f549c8dc076f6b3b4cf5a1bc47bf269d",
             "a62edd414978f91f7391a59fc1a72372", "08342be5632619fd1b1251b7ad2b2c84", "0d6b7deda824b4fc42b65cb87e1d4d14"]
rc_hashes = [(rc_objects[indx], value) for indx, value in enumerate(rc_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", rc_hashes)
def test_reverse_complement(seqbuddy, next_hash):
    tester = Sb.reverse_complement(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

def test_reverse_complement_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein into -sfr
    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis/Mnemiopsis_pep.fa"))
    with pytest.raises(ValueError):
        Sb.reverse_complement(seqbuddy)

li_hashes = ["1c4a395d8aa3496d990c611c3b6c4d0a", "1c4a395d8aa3496d990c611c3b6c4d0a", "1c4a395d8aa3496d990c611c3b6c4d0a",
             "78a9289ab2d508a13c76cf9f5a308cc5", "1c4a395d8aa3496d990c611c3b6c4d0a", "1c4a395d8aa3496d990c611c3b6c4d0a"]
li_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(li_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", li_hashes)
def test_list_ids(seqbuddy, next_hash):
    tester = Sb.list_ids(seqbuddy, 1)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash

#  are num_seqs and ave_seq_length even worth testing?

# 'cts', '--concat_seqs'

cts_hashes = ["2e46edb78e60a832a473397ebec3d187", "7421c27be7b41aeedea73ff41869ac47",
              "494988ffae2ef3072c1619eca8a0ff3b", "710cad348c5560446daf2c916ff3b3e4",
              "494988ffae2ef3072c1619eca8a0ff3b", "494988ffae2ef3072c1619eca8a0ff3b",
              "46741638cdf7abdf53c55f79738ee620", "8d0bb4e5004fb6a1a0261c30415746b5",
              "2651271d7668081cde8012db4f9a6574", "36526b8e0360e259d8957fa2261cf45a",
              "2651271d7668081cde8012db4f9a6574", "2651271d7668081cde8012db4f9a6574"]
cts_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(cts_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", cts_hashes)
def test_concat_seqs(seqbuddy, next_hash):
    tester = Sb.concat_seqs(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

# TODO fd2p and fp2d

# 'ri', '--rename_ids'

ri_hashes = ["0672973b580581f15cf2ce467b89144e", "3847ad5210a85b8db59e256261552ee7", "243024bfd2f686e6a6e0ef65aa963494",
             "83f10d1be7a5ba4d363eb406c1c84ac7", "973e3d7138b78db2bb3abda8a9323226", "4289f03afb6c9f8a8b0d8a75bb60a2ce"]
ri_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(ri_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", ri_hashes)
def test_rename_ids(seqbuddy, next_hash):  # Probably don't need to test ALL the files
    tester = Sb.rename(seqbuddy, 'Panx', 'Test', 0)
    assert seqs_to_hash(tester) == next_hash

#TODO combine_features

# How do you test the 'shuffle' method?

# 'oi', '--order_ids'

oi_files = [seq_files[0],seq_files[6]]
oi_hashes = ["d103a7a41a5644614d59e0c44469d1b2", "945b2b43a423a5371ad7e90adda6e703"]
oi_hashes = [(Sb.SeqBuddy(resource(oi_files[indx])), value) for indx, value in enumerate(oi_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", oi_hashes)
def test_order_ids(seqbuddy, next_hash):
    tester = Sb.order_ids(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

oi_rev_hashes = ["3658c4b79cd7e8dfd6afe1a9cddc2dfa", "8c4e450b72410c37683f83e528c9a610"]
oi_rev_hashes = [(Sb.SeqBuddy(resource(oi_files[indx])), value) for indx, value in enumerate(oi_rev_hashes)]
@pytest.mark.parametrize("seqbuddy,next_hash", oi_rev_hashes)
def test_order_ids_rev(seqbuddy, next_hash):
    tester = Sb.order_ids(seqbuddy, _reverse=True)
    assert seqs_to_hash(tester) == next_hash

# 'ofp', '--order_features_by_position'

ofp_hashes = ["25073539df4a982b7f99c72dd280bb8f", "c10d136c93f41db280933d5b3468f187"]
ofp_hashes = [(Sb.SeqBuddy(resource(oi_files[indx])), value) for indx, value in enumerate(ofp_hashes)]
@pytest.mark.parametrize("seqbuddy, next_hash", ofp_hashes)
def test_order_features_by_position(seqbuddy, next_hash):
    tester = Sb.order_features_by_position(seqbuddy)
    assert seqs_to_hash(tester) == next_hash

ofp_rev_hashes = ["25073539df4a982b7f99c72dd280bb8f", "c10d136c93f41db280933d5b3468f187"]
ofp_rev_hashes = [(Sb.SeqBuddy(resource(oi_files[indx])), value) for indx, value in enumerate(ofp_rev_hashes)]
@pytest.mark.parametrize("seqbuddy, next_hash", ofp_rev_hashes)
def test_order_features_by_position_rev(seqbuddy, next_hash):
    tester = Sb.order_features_by_position(seqbuddy, _reverse=True)
    assert seqs_to_hash(tester) == next_hash

# 'sf', '--screw_formats'

# last three formats rebuilt each time to work correctly in parallel
fasta_files = ["Mnemiopsis/Mnemiopsis_cds.fa", "Mnemiopsis/Mnemiopsis_pep.fa"]
fasta_files = [Sb.SeqBuddy(resource(file)) for file in fasta_files]
gb_files = ["Mnemiopsis/Mnemiopsis_cds.gb", "Mnemiopsis/Mnemiopsis_pep.gb"]
gb_files = [Sb.SeqBuddy(resource(file)) for file in gb_files]
nex_files = ["Mnemiopsis/Mnemiopsis_cds.nex", "Mnemiopsis/Mnemiopsis_pep.nex"]
# nex_files = [Sb.SeqBuddy(resource(file)) for file in nex_files]
phyr_files = ["Mnemiopsis/Mnemiopsis_cds.phyr", "Mnemiopsis/Mnemiopsis_pep.phyr"]
# phyr_files = [Sb.SeqBuddy(resource(file)) for file in phyr_files]
stklm_files = ["Mnemiopsis/Mnemiopsis_cds.stklm", "Mnemiopsis/Mnemiopsis_pep.stklm"]
# stklm_files = [Sb.SeqBuddy(resource(file)) for file in stklm_files]
@pytest.mark.parametrize("indx", [0,1])
def test_screw_formats_fa_gb(indx):
    fasta = fasta_files[indx]
    gb = gb_files[indx]
    assert seqs_to_hash(fasta) == seqs_to_hash(Sb.screw_formats(gb, "fasta"))

other_files = [nex_files, phyr_files, stklm_files]
screw_files = []
for l1 in other_files:
    for l2 in other_files:
        if l1 is not l2:
            for indx in range(2):
                screw_files.append((indx, l1, l2))

@pytest.mark.parametrize("indx,l1,l2", screw_files) # fails when stklm is the out_format
def test_screw_formats_other(indx, l1, l2):
    sb1 = Sb.SeqBuddy(resource(l1[indx]))
    sb2 = Sb.SeqBuddy(resource(l2[indx]))
    sb2 = Sb.screw_formats(sb2, sb1.out_format)
    assert seqs_to_hash(sb1) == seqs_to_hash(sb2)


if __name__ == '__main__':
    debug = Sb.order_features_alphabetically(sb_objects[1])
    print(seqs_to_hash(debug, "string"))
