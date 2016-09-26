#!/usr/bin/env python3
# coding=utf-8
""" tests basic functionality of SeqBuddy class """
import pytest
from Bio.Alphabet import IUPAC
from collections import OrderedDict
from unittest import mock
import os
from shutil import which
import subprocess

from ... import buddy_resources as br
from ...SeqBuddy import SeqBuddy, hash_ids, pull_recs, make_copy,\
    _guess_alphabet, _guess_format, _feature_rc, _check_for_blast_bin, Popen


def test_instantiate_seqbuddy_from_file(sb_resources):
    for _path in sb_resources.get_list("", mode="paths"):
        assert type(SeqBuddy(_path)) == SeqBuddy


def test_instantiate_seqbuddy_from_handle(sb_resources):
    for _path in sb_resources.get_list("", mode="paths"):
        with open(_path, 'r') as ifile:
            assert type(SeqBuddy(ifile)) == SeqBuddy


def test_instantiate_seqbuddy_from_raw(sb_resources):
    for _path in sb_resources.get_list("", mode="paths"):
        with open(_path, 'r') as ifile:
            assert type(SeqBuddy(ifile.read(), in_format="raw")) == SeqBuddy

        with open(_path, 'r') as ifile:
            assert type(SeqBuddy(ifile, in_format="raw")) == SeqBuddy


def test_instantiate_seqbuddy_from_seqbuddy(sb_resources, hf):
    for _path in sb_resources.get_list("", mode="paths"):
        input_buddy = SeqBuddy(_path)
        tester = SeqBuddy(input_buddy)
        assert hf.buddy2hash(input_buddy) == hf.buddy2hash(tester)


def test_alpha_arg_dna(sb_resources):
    tester = SeqBuddy(sb_resources.get_one("d f", mode="paths"), alpha='dna')
    assert tester.alpha is IUPAC.ambiguous_dna


def test_alpha_arg_rna(sb_resources):
    tester = SeqBuddy(sb_resources.get_one("r f", mode="paths"), alpha='rna')
    assert tester.alpha is IUPAC.ambiguous_rna


def test_alpha_arg_prot(sb_resources):
    tester = SeqBuddy(sb_resources.get_one("p f", mode="paths"), alpha='prot')
    assert tester.alpha is IUPAC.protein


def test_alpha_arg_guess(sb_resources):
    tester = SeqBuddy(sb_resources.get_one("d f", mode="paths"), alpha='foo')
    assert tester.alpha is IUPAC.ambiguous_dna


def test_seqlist_error():
    with pytest.raises(TypeError):
        SeqBuddy([str, dict])


# ##################### SeqBuddy methods ###################### ##
def test_to_dict(sb_resources, hf):
    tester = str(sb_resources.get_one("o d f").to_dict())
    assert hf.string2hash(tester) == '2311d1712d41c5ec9c23ad107c8a06c3'

    with pytest.raises(RuntimeError):
        tester = SeqBuddy(">duplicate_id\nATGCTCGTA\n>duplicate_id\nATGCTCGTCGATGCT\n")
        tester.to_dict()


def test_to_string(sb_resources, hf, capsys):
    tester = sb_resources.get_one("d f")
    assert hf.string2hash(str(tester)) == "b831e901d8b6b1ba52bad797bad92d14"

    tester = sb_resources.get_one("d g")
    assert hf.string2hash(str(tester)) == "2e02a8e079267bd9add3c39f759b252c"

    tester.out_format = "raw"
    assert hf.string2hash(str(tester)) == "5d00d481e586e287f32d2d29916374ca"

    tester = sb_resources.get_one("d n")
    pull_recs(tester, "α[2-9]")

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
        str(SeqBuddy(">foo\nATGATGATGTAGT\n>bar\nATGATGATGTAGT\n>foo\nATGATGATGTAGT\n", out_format="phylip"))
    out, err = capsys.readouterr()
    assert "Warning: Phylip format returned a 'repeat name' error, probably due to truncation" in err

    tester = SeqBuddy(">fooooooooooooooobar\nATGATGATGTAGT\n>bar\nATGATGATGTAGT\n>fooooooooooooooobaz\nATGATGATGTAGT\n",
                      out_format="phylip")
    assert hf.string2hash(str(tester)) == "d8b82d5eea15918aac180e5d1095d5ca"
    out, err = capsys.readouterr()
    assert "Attempting phylip-relaxed." in err

    tester = SeqBuddy(">fooooooooooooooobaaaaaaaaaaaaaaaaaar\nATGATGATGTAGT\n>bar\nATGATGATGTAGT\n",
                      out_format="gb")
    assert hf.string2hash(str(tester)) == "473777dd047b46ef727a4b680247e374"
    out, err = capsys.readouterr()
    assert "Warning: Genbank format returned an 'ID too long' error. Format changed to EMBL." in err


def test_write(sb_resources, hf):
    temp_dir = br.TempDir()
    tester = sb_resources.get_one("d g")
    tester.write("%s/sequences.gb" % temp_dir.path)
    with open("%s/sequences.gb" % temp_dir.path) as ifile:
        assert hf.string2hash(ifile.read()) == "2e02a8e079267bd9add3c39f759b252c"

    tester.write("%s/sequences.fa" % temp_dir.path, out_format="fasta")
    with open("%s/sequences.fa" % temp_dir.path) as ifile:
        assert hf.string2hash(ifile.read()) == "25073539df4a982b7f99c72dd280bb8f"


def test_print_hashmap(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    hash_ids(tester)
    test_hashes = ["FEhFs96uVr", "5dOVoJsEaC", "muOhKHqlRK", "99id32X9JY", "hflijfeJXB", "0m9x7xeSqC", "qwgaHU3fms",
                   "uD7zXF2uEp", "btvnHXOJbc", "GiHvUV1n55", "dJm5uViNsC", "to4ctKvNG7", "VN579cevl3"]
    orig_ids = [rec_id for _hash, rec_id in tester.hash_map.items()]
    tester.hash_map = OrderedDict(zip(test_hashes, orig_ids))
    assert hf.string2hash(tester.print_hashmap()) == "ab38a224d002a5b227265b8211c9f7bc"


def test_reverse_hashmap(sb_resources):
    tester = sb_resources.get_one("d f")
    tester_str = str(tester)
    hash_ids(tester)
    assert tester_str != str(tester)
    tester.reverse_hashmap()
    assert tester_str == str(tester)


# ################################################# HELPER FUNCTIONS ################################################# #
# ToDo: Missing tests for --> _add_buddy_data, FeatureReMapper
# ######################  '_check_for_blast_bin' ###################### #
def test_check_blast_bin_success():
    for _bin in ["blastn", "blastp", "blastdbcmd", "makeblastdb"]:
        assert _check_for_blast_bin(_bin)


def test_check_blast_bin_error(capsys):
    with mock.patch.dict(os.environ, {"PATH": ""}):
        assert not _check_for_blast_bin("blastp")
        out, err = capsys.readouterr()
        assert "blastp binary not found." in err
        assert "Please install BLAST+ executables." in err


def test_check_blast_bin_download(monkeypatch, capsys):
    def patch_ask(*args):
        print(args)
        _blastp = temp_dir.subfile("blastp")
        os.chmod(_blastp, 0o777)
        return True

    monkeypatch.setattr(subprocess, "Popen", lambda *_: True)
    temp_dir = br.TempDir()
    monkeypatch.setitem(os.environ, "PATH", temp_dir.path)
    assert not which("blastp")

    assert not _check_for_blast_bin("blastp")
    out, err = capsys.readouterr()
    assert "Please install BLAST+ executables." in err

    blastp = temp_dir.subfile("blastp")
    os.chmod(blastp, 0o777)
    assert which("blastp")
    assert _check_for_blast_bin("blastp")

    temp_dir.del_subfile("blastp")
    assert not which("blastp")

    conda = temp_dir.subfile("conda")
    os.chmod(conda, 0o777)
    monkeypatch.setattr(br, "ask", lambda _: False)
    assert not _check_for_blast_bin("blastp")
    _out, _err = capsys.readouterr()
    assert "blastp binary not found." in _err
    assert "Abort..." in _err

    monkeypatch.setattr(br, "ask", lambda _: True)
    assert not _check_for_blast_bin("blastp")
    out, err = capsys.readouterr()
    assert "Failed to install BLAST+ executables with Conda." in err

    monkeypatch.setattr(br, "ask", patch_ask)
    assert _check_for_blast_bin("blastp")


# ######################  '_feature_rc' ###################### #
def test_feature_rc(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    seq1 = tester.records[0]
    hashed = ["51eb205f066387e146fb7a3ce5b122c1", "9ab8296fb3443198674d90abe3311ba6",
              "10018d1b15c7f76a6333ac3bf96d2d07", "273463b9eace12d2eeadbf272692d73e",
              "c452d66d13120cd6eb5f041b7c37dd27", "1811b0695dba1fc3fe431a6ee00ef359"]
    for feature in zip(seq1.features, hashed):
        assert hf.string2hash(str(_feature_rc(feature[0], 1203))) == feature[1]

    with pytest.raises(TypeError):
        feature = seq1.features[0]
        feature.location = {}
        _feature_rc(feature, 1203)


# ######################  'guess_alphabet' ###################### #
def test_guess_alphabet(sb_resources):
    tester = sb_resources.get_one("d f")
    assert _guess_alphabet(tester) == IUPAC.ambiguous_dna

    tester = sb_resources.get_one("p f")
    assert _guess_alphabet(tester) == IUPAC.protein

    tester = sb_resources.get_one("r f")
    assert _guess_alphabet(tester) == IUPAC.ambiguous_rna

    tester = SeqBuddy(">Seq1", in_format="fasta")
    assert not _guess_alphabet(tester)


# ######################  'guess_format' ###################### #
def test_guess_stockholm(hf):
    assert _guess_format("%s%sMnemiopsis_cds.stklm" % (hf.resource_path, os.path.sep)) == "stockholm"

    with open("%s%sMnemiopsis_cds.stklm" % (hf.resource_path, os.path.sep), "r") as ifile:
        assert _guess_format(ifile) == "stockholm"

    seqbuddy = SeqBuddy("%s%sMnemiopsis_cds.stklm" % (hf.resource_path, os.path.sep))
    assert _guess_format(seqbuddy) == "stockholm"


def test_guess_fasta(hf):
    assert _guess_format("%s%sMnemiopsis_cds.fa" % (hf.resource_path, os.path.sep)) == "fasta"

    with open("%s%sMnemiopsis_cds.fa" % (hf.resource_path, os.path.sep), "r") as ifile:
        assert _guess_format(ifile) == "fasta"

    seqbuddy = SeqBuddy("%s%sMnemiopsis_cds.fa" % (hf.resource_path, os.path.sep))
    assert _guess_format(seqbuddy) == "fasta"


def test_guess_gb(hf):
    assert _guess_format("%s%sMnemiopsis_cds.gb" % (hf.resource_path, os.path.sep)) == "gb"

    with open("%s%sMnemiopsis_cds.gb" % (hf.resource_path, os.path.sep), "r") as ifile:
        assert _guess_format(ifile) == "gb"

    seqbuddy = SeqBuddy("%s%sMnemiopsis_cds.gb" % (hf.resource_path, os.path.sep))
    assert _guess_format(seqbuddy) == "gb"

'''
def test_guess_phylipss(hf):
    print("%s%sMnemiopsis_cds.physs" % (hf.resource_path, os.path.sep))
    assert _guess_format("%s%sMnemiopsis_cds.physs" % (hf.resource_path, os.path.sep)) == "phylipss"

    with open("%s%sMnemiopsis_cds.physs" % (hf.resource_path, os.path.sep), "r") as ifile:
        assert _guess_format(ifile) == "phylipss"

    seqbuddy = SeqBuddy("%s%sMnemiopsis_cds.physs" % (hf.resource_path, os.path.sep))
    assert _guess_format(seqbuddy) == "phylipss"
'''

def test_guess_format(sb_resources, sb_odd_resources):
    assert _guess_format(["foo", "bar"]) == "gb"
    assert _guess_format(sb_resources.get_one("d f")) == "fasta"
    assert _guess_format(sb_resources.get_one("d f", mode="paths")) == "fasta"
    assert _guess_format(sb_odd_resources["blank"]) == "empty file"
    with pytest.raises(br.GuessError):
        _guess_format("foo")

    temp_file = br.TempFile()
    temp_file.write('''\
<?xml version="1.0" encoding="ISO-8859-1"?>
<nex:nexml
    version="0.9"
    xsi:schemaLocation="http://www.nexml.org/2009 ../xsd/nexml.xsd"
    xmlns="http://www.nexml.org/2009"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xmlns:xml="http://www.w3.org/XML/1998/namespace"
    xmlns:nex="http://www.nexml.org/2009"
    xmlns:xsd="http://www.w3.org/2001/XMLSchema#"
''')

    assert not _guess_format(temp_file.path)


# ######################  'GuessError' ###################### #
def test_guesserror_raw_seq():
    with pytest.raises(br.GuessError):
        SeqBuddy("JSKHGLHGLSDKFLSDYUIGJVSBDVHJSDKGIUSUEWUIOIFUBCVVVBVNNJS{QF(*&#@$(*@#@*(*(%")
    try:
        SeqBuddy("JSKHGLHGLSDKFLSDYUIGJVSBDVHJSDKGIUSUEWUIOIFUBCVVVBVNNJS{QF(*&#@$(*@#@*(*(%")
    except br.GuessError as e:
        assert "File not found, or could not determine format from raw input" in str(e)


def test_guesserror_infile(sb_odd_resources):
    with pytest.raises(br.GuessError):
        SeqBuddy(sb_odd_resources["gibberish"])


def test_guesserror_in_handle(sb_odd_resources):
    with pytest.raises(br.GuessError):
        with open(sb_odd_resources["gibberish"], "r") as ifile:
            SeqBuddy(ifile)


def test_no__input():
    with pytest.raises(TypeError):
        # noinspection PyArgumentList
        SeqBuddy()


# ######################  'make_copy' ###################### #
def test_make_copy(sb_resources, hf):
    tester = SeqBuddy(sb_resources.get_one("d f", mode="paths"))
    tester_copy = make_copy(tester)
    assert hf.buddy2hash(tester) == hf.buddy2hash(tester_copy)
