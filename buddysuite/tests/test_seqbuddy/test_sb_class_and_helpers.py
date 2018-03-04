#!/usr/bin/env python3
# coding=utf-8
""" tests basic functionality of SeqBuddy class """
import pytest
from Bio.Alphabet import IUPAC
from collections import OrderedDict
import os
import buddy_resources as br
import SeqBuddy as Sb


def test_instantiate_seqbuddy_from_file(sb_resources):
    for _path in sb_resources.get_list("", mode="paths"):
        assert type(Sb.SeqBuddy(_path)) == Sb.SeqBuddy


def test_instantiate_seqbuddy_from_handle(sb_resources):
    for _path in sb_resources.get_list("", mode="paths"):
        with open(_path, 'r', encoding="utf-8") as ifile:
            assert type(Sb.SeqBuddy(ifile)) == Sb.SeqBuddy


def test_instantiate_seqbuddy_from_raw(sb_resources):
    for _path in sb_resources.get_list("", mode="paths"):
        with open(_path, 'r') as ifile:
            assert type(Sb.SeqBuddy(ifile.read(), in_format="raw")) == Sb.SeqBuddy

        with open(_path, 'r') as ifile:
            assert type(Sb.SeqBuddy(ifile, in_format="raw")) == Sb.SeqBuddy


def test_instantiate_seqbuddy_from_seqbuddy(sb_resources, hf):
    for _path in sb_resources.get_list("", mode="paths"):
        input_buddy = Sb.SeqBuddy(_path)
        tester = Sb.SeqBuddy(input_buddy)
        assert hf.buddy2hash(input_buddy) == hf.buddy2hash(tester)


def test_alpha_arg_dna(sb_resources):
    tester = Sb.SeqBuddy(sb_resources.get_one("d f", mode="paths"), alpha='dna')
    assert tester.alpha is IUPAC.ambiguous_dna


def test_alpha_arg_rna(sb_resources):
    tester = Sb.SeqBuddy(sb_resources.get_one("r f", mode="paths"), alpha='rna')
    assert tester.alpha is IUPAC.ambiguous_rna


def test_alpha_arg_prot(sb_resources):
    tester = Sb.SeqBuddy(sb_resources.get_one("p f", mode="paths"), alpha='prot')
    assert tester.alpha is IUPAC.protein


def test_alpha_arg_guess(sb_resources):
    tester = Sb.SeqBuddy(sb_resources.get_one("d f", mode="paths"), alpha='foo')
    assert tester.alpha is IUPAC.ambiguous_dna


def test_seqlist_error():
    with pytest.raises(TypeError):
        Sb.SeqBuddy([str, dict])


# ##################### SeqBuddy methods ###################### ##
def test_to_dict(sb_resources, hf):
    tester = str(sb_resources.get_one("o d f").to_dict())
    assert hf.string2hash(tester) == '2311d1712d41c5ec9c23ad107c8a06c3'

    with pytest.raises(RuntimeError):
        tester = Sb.SeqBuddy(">duplicate_id\nATGCTCGTA\n>duplicate_id\nATGCTCGTCGATGCT\n")
        tester.to_dict()


def test_to_string(sb_resources, hf, capsys):
    tester = sb_resources.get_one("d f")
    assert hf.string2hash(str(tester)) == "b831e901d8b6b1ba52bad797bad92d14"

    tester = sb_resources.get_one("d g")
    assert hf.string2hash(str(tester)) == "2e02a8e079267bd9add3c39f759b252c"

    tester.out_format = "raw"
    assert hf.string2hash(str(tester)) == "5d00d481e586e287f32d2d29916374ca"

    tester = sb_resources.get_one("d n")
    Sb.pull_recs(tester, "α[2-9]")

    tester.out_format = "phylip"
    assert hf.string2hash(str(tester)) == "6a4d62e1ee130b324cce48323c6d1d41"

    tester.out_format = "phylip-relaxed"
    assert hf.string2hash(str(tester)) == "4c2c5900a57aad343cfdb8b35a8f8442"

    tester.out_format = "phylipss"
    assert hf.string2hash(str(tester)) == "089cfb52076e63570597a74b2b000660"

    tester.out_format = "phylipsr"
    assert hf.string2hash(str(tester)) == "58a74f5e08afa0335ccfed0bdd94d3f2"

    tester.records = []
    assert str(tester) == "Error: No sequences in object.\n"

    with pytest.raises(ValueError):
        str(Sb.SeqBuddy(">foo\nATGATGATGTAGT\n>bar\nATGATGATGTAGT\n>foo\nATGATGATGTAGT\n", out_format="phylip"))
    out, err = capsys.readouterr()
    assert "Warning: Phylip format returned a 'repeat name' error, probably due to truncation" in err

    tester = Sb.SeqBuddy(">fooooooooooooooobar\nATGATGATGTAGT\n>bar\nATGATGATGTAGT\n"
                         ">fooooooooooooooobaz\nATGATGATGTAGT\n", out_format="phylip")
    assert hf.string2hash(str(tester)) == "d8b82d5eea15918aac180e5d1095d5ca"
    out, err = capsys.readouterr()
    assert "Attempting phylip-relaxed." in err

    tester = Sb.SeqBuddy(">fooooooooooooooobaaaaaaaaaaaaaaaaaar\nATGATGATGTAGT\n>bar\nATGATGATGTAGT\n",
                         out_format="gb")
    assert hf.string2hash(str(tester)) == "e10a8872a05242a32f2f29c309d150f9", print(str(tester))
    out, err = capsys.readouterr()
    assert "Warning: Genbank format returned an 'ID too long' error. Format changed to EMBL." in err


def test_write(sb_resources, hf):
    temp_dir = br.TempDir()
    tester = sb_resources.get_one("d g")
    tester.write("%s/sequences.gb" % temp_dir.path)
    with open("%s/sequences.gb" % temp_dir.path, encoding="utf-8") as ifile:
        assert hf.string2hash(ifile.read()) == "2e02a8e079267bd9add3c39f759b252c"

    tester.write("%s/sequences.fa" % temp_dir.path, out_format="fasta")
    with open("%s/sequences.fa" % temp_dir.path, encoding="utf-8") as ifile:
        data = ifile.read()
        assert hf.string2hash(data) == "6a9b3b554aa9ddb90ea62967bd26d5b7"


def test_print_hashmap(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    Sb.hash_ids(tester)
    test_hashes = ["FEhFs96uVr", "5dOVoJsEaC", "muOhKHqlRK", "99id32X9JY", "hflijfeJXB", "0m9x7xeSqC", "qwgaHU3fms",
                   "uD7zXF2uEp", "btvnHXOJbc", "GiHvUV1n55", "dJm5uViNsC", "to4ctKvNG7", "VN579cevl3"]
    orig_ids = [rec_id for _hash, rec_id in tester.hash_map.items()]
    tester.hash_map = OrderedDict(zip(test_hashes, orig_ids))
    assert hf.string2hash(tester.print_hashmap()) == "ab38a224d002a5b227265b8211c9f7bc"


def test_reverse_hashmap(sb_resources):
    tester = sb_resources.get_one("d f")
    tester_str = str(tester)
    Sb.hash_ids(tester)
    assert tester_str != str(tester)
    tester.reverse_hashmap()
    assert tester_str == str(tester)


# ################################################# HELPER FUNCTIONS ################################################# #
# ToDo: Missing tests for --> _add_buddy_data, FeatureReMapper
# ######################  '_check_for_blast_bin' ###################### #
def test_check_blast_bin(monkeypatch, capsys):
    monkeypatch.setattr(Sb, "which", lambda *_: True)
    assert Sb._check_for_blast_bin("blastp")

    monkeypatch.setattr(Sb, "which", lambda *_: False)
    assert not Sb._check_for_blast_bin("blastn")
    out, err = capsys.readouterr()
    assert "blastn binary not found. Please install BLAST+ executables.\n" in err


# ######################  '_feature_rc' ###################### #
def test_feature_rc(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    seq1 = tester.records[0]
    hashed = ["51eb205f066387e146fb7a3ce5b122c1", "9ab8296fb3443198674d90abe3311ba6",
              "10018d1b15c7f76a6333ac3bf96d2d07", "273463b9eace12d2eeadbf272692d73e",
              "c452d66d13120cd6eb5f041b7c37dd27", "1811b0695dba1fc3fe431a6ee00ef359"]
    for feature in zip(seq1.features, hashed):
        feature[0].sub_features = []
        assert hf.string2hash(str(Sb._feature_rc(feature[0], 1203))) == feature[1]

    with pytest.raises(TypeError):
        feature = seq1.features[0]
        feature.location = {}
        Sb._feature_rc(feature, 1203)


# ######################  'guess_alphabet' ###################### #
def test_guess_alphabet(sb_resources):
    tester = sb_resources.get_one("d f")
    assert Sb.guess_alphabet(tester) == IUPAC.ambiguous_dna

    tester = sb_resources.get_one("p f")
    assert Sb.guess_alphabet(tester) == IUPAC.protein

    tester = sb_resources.get_one("r f")
    assert Sb.guess_alphabet(tester) == IUPAC.ambiguous_rna

    tester = Sb.SeqBuddy(">Seq1", in_format="fasta")
    assert not Sb.guess_alphabet(tester)


# ######################  'GuessError' ###################### #
def test_guesserror_raw_seq():
    with pytest.raises(br.GuessError):
        Sb.SeqBuddy("JSKHGLHGLSDKFLSDYUIGJVSBDVHJSDKGIUSUEWUIOIFUBCVVVBVNNJS{QF(*&#@$(*@#@*(*(%")
    try:
        Sb.SeqBuddy("JSKHGLHGLSDKFLSDYUIGJVSBDVHJSDKGIUSUEWUIOIFUBCVVVBVNNJS{QF(*&#@$(*@#@*(*(%")
    except br.GuessError as e:
        assert "File not found, or could not determine format from raw input" in str(e)


def test_guesserror_infile(sb_odd_resources):
    with pytest.raises(br.GuessError):
        Sb.SeqBuddy(sb_odd_resources["gibberish"])


def test_guesserror_in_handle(sb_odd_resources):
    with pytest.raises(br.GuessError):
        with open(sb_odd_resources["gibberish"], "r", encoding="utf-8") as ifile:
            Sb.SeqBuddy(ifile)


def test_no__input():
    with pytest.raises(TypeError):
        # noinspection PyArgumentList
        Sb.SeqBuddy()


# ######################  'make_copy' ###################### #
def test_make_copy(sb_resources, hf):
    tester = Sb.SeqBuddy(sb_resources.get_one("d f", mode="paths"))
    tester_copy = Sb.make_copy(tester)
    assert hf.buddy2hash(tester) == hf.buddy2hash(tester_copy)
