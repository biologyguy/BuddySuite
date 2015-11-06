#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: seqbuddy_tests.py
version: 1, alpha
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Collection of PyTest unit tests for the SeqBuddy.py program
"""

import pytest
from hashlib import md5
import os
import re
import sys
import argparse
import io
from copy import deepcopy

from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC
from unittest import mock

sys.path.insert(0, "./")
import buddy_resources as br
import SeqBuddy as Sb
import MyFuncs

VERSION = Sb.VERSION
WRITE_FILE = MyFuncs.TempFile()
TEMP_DIR = MyFuncs.TempDir()
BACKUP_PATH = os.environ["PATH"]


def fmt(prog):
    return br.CustomHelpFormatter(prog)

parser = argparse.ArgumentParser(prog="SeqBuddy.py", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                 description='''\
\033[1mSeqBuddy\033[m
  See your sequence files. Be your sequence files.

\033[1mUsage examples\033[m:
  SeqBuddy.py "/path/to/seq_file" -<cmd>
  SeqBuddy.py "/path/to/seq_file" -<cmd> | SeqBuddy.py -<cmd>
  SeqBuddy.py "ATGATGCTAGTC" -f "raw" -<cmd>
''')

br.flags(parser, ("sequence", "Supply file path(s) or raw sequence. If piping sequences "
                              "into SeqBuddy this argument can be left blank."),
         br.sb_flags, br.sb_modifiers, VERSION)

# This is to allow py.test to work with the -x flag
parser.add_argument("-x", nargs="?")
parser.add_argument("-m", nargs="?")
parser.add_argument("-n", nargs="?")
parser.add_argument("--cov", nargs="?")
parser.add_argument("--cov-report", nargs="?")
in_args = parser.parse_args()


def seqs_to_hash(_seqbuddy, mode='hash'):
    if _seqbuddy.out_format in ["gb", "genbank"]:
        for _rec in _seqbuddy.records:
            try:
                if re.search("(\. )+", _rec.annotations['organism']):
                    _rec.annotations['organism'] = "."
            except KeyError:
                pass

    if _seqbuddy.out_format == "phylipsr":
        WRITE_FILE.write(br.phylip_sequential_out(_seqbuddy, relaxed=True, _type="seqbuddy"))
    elif _seqbuddy.out_format == "phylipss":
        WRITE_FILE.write(br.phylip_sequential_out(_seqbuddy, relaxed=False, _type="seqbuddy"))
    else:
        _seqbuddy.write(WRITE_FILE.path)

    seqs_string = "{0}\n".format(WRITE_FILE.read().rstrip())

    if mode != "hash":
        return seqs_string

    _hash = md5(seqs_string.encode()).hexdigest()
    return _hash

root_dir = os.getcwd()


def string2hash(_input):
    return md5(_input.encode()).hexdigest()


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
        assert type(Sb.SeqBuddy(ifile.read(), in_format="raw")) == Sb.SeqBuddy

    with open(resource(seq_file), 'r') as ifile:
        assert type(Sb.SeqBuddy(ifile, in_format="raw")) == Sb.SeqBuddy


@pytest.mark.parametrize("seq_file", seq_files)
def test_instantiate_seqbuddy_from_seqbuddy(seq_file):
    input_buddy = Sb.SeqBuddy(resource(seq_file))
    tester = Sb.SeqBuddy(input_buddy)
    assert seqs_to_hash(input_buddy) == seqs_to_hash(tester)


def test_alpha_arg_dna():
    tester = Sb.SeqBuddy(resource(seq_files[0]), alpha='dna')
    assert tester.alpha is IUPAC.ambiguous_dna


def test_alpha_arg_rna():
    tester = Sb.SeqBuddy(resource(seq_files[0]), alpha='rna')
    assert tester.alpha is IUPAC.ambiguous_rna


def test_alpha_arg_prot():
    tester = Sb.SeqBuddy(resource(seq_files[6]), alpha='prot')
    assert tester.alpha is IUPAC.protein


def test_alpha_arg_guess():
    tester = Sb.SeqBuddy(resource(seq_files[0]), alpha='dgasldfkjhgaljhetlfdjkfg')
    assert tester.alpha is IUPAC.ambiguous_dna


def test_seqlist_error():
    with pytest.raises(TypeError):
        Sb.SeqBuddy([str, dict])


# ##################### SeqBuddy methods ###################### ##
def test_to_dict():
    tester = str(Sb.SeqBuddy(resource("Mnemiopsis_cds.fa")).to_dict())
    assert string2hash(tester) == '2311d1712d41c5ec9c23ad107c8a06c3'

    with pytest.raises(RuntimeError):
        tester = Sb.SeqBuddy(">duplicate_id\nATGCTCGTA\n>duplicate_id\nATGCTCGTCGATGCT\n")
        tester.to_dict()


def test_to_string(capsys):
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    assert string2hash(str(tester)) == "b831e901d8b6b1ba52bad797bad92d14"
    tester.print()
    out, err = capsys.readouterr()
    assert string2hash(out) == "b831e901d8b6b1ba52bad797bad92d14"

    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.gb"))
    assert string2hash(str(tester)) == "2e02a8e079267bd9add3c39f759b252c"

    tester.out_format = "raw"
    assert string2hash(str(tester)) == "5d00d481e586e287f32d2d29916374ca"

    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.nex"))
    Sb.pull_recs(tester, "α[2-9]")

    tester.out_format = "phylip"
    assert string2hash(str(tester)) == "6a4d62e1ee130b324cce48323c6d1d41"

    tester.out_format = "phylip-relaxed"
    assert string2hash(str(tester)) == "4c2c5900a57aad343cfdb8b35a8f8442"

    tester.out_format = "phylipss"
    assert string2hash(str(tester)) == "089cfb52076e63570597a74b2b000660"

    tester.out_format = "phylipsr"
    assert string2hash(str(tester)) == "58a74f5e08afa0335ccfed0bdd94d3f2"

    tester.records = []
    assert str(tester) == "Error: No sequences in object.\n"


# Now that we know that all the files are being turned into SeqBuddy objects okay, make them all objects so it doesn't
# need to be done over and over for each subsequent test.
sb_objects = [Sb.SeqBuddy(resource(x)) for x in seq_files]


# ################################################# HELPER FUNCTIONS ################################################# #
# ######################  'GuessError' ###################### #
def test_guesserror_raw_seq():
    with pytest.raises(br.GuessError):
        Sb.SeqBuddy("JSKHGLHGLSDKFLSDYUIGJVSBDVHJSDKGIUSUEWUIOIFUBCVVVBVNNJS{QF(*&#@$(*@#@*(*(%")
    try:
        Sb.SeqBuddy("JSKHGLHGLSDKFLSDYUIGJVSBDVHJSDKGIUSUEWUIOIFUBCVVVBVNNJS{QF(*&#@$(*@#@*(*(%")
    except br.GuessError as e:
        assert "File not found, or could not determine format from raw input" in str(e)


def test_guesserror_infile():
    with pytest.raises(br.GuessError):
        Sb.SeqBuddy(resource("gibberish.fa"))


def test_guesserror_in_handle():
    with pytest.raises(br.GuessError):
        with open(resource("gibberish.fa"), "r") as ifile:
            Sb.SeqBuddy(ifile)


def test_no__input():
    with pytest.raises(TypeError):
        Sb.SeqBuddy()


# ######################  'make_copy' ###################### #
def test_make_copy():
    assert seqs_to_hash(Sb._make_copy(sb_objects[0])) == seqs_to_hash(sb_objects[0])


# ######################  '_check_for_blast_bin' ###################### #
@pytest.mark.internet
@pytest.mark.slow
def test_check_blast_bin(capsys):
    for _bin in ["blastn", "blastp", "blastdbcmd", "makeblastdb"]:
        assert Sb._check_for_blast_bin(_bin)

    with mock.patch.dict(os.environ, {"PATH": ""}):
        with mock.patch('MyFuncs.ask', return_value=False):
            assert not Sb._check_for_blast_bin("blastp")

        with mock.patch('MyFuncs.ask', return_value=True):
            with mock.patch("SeqBuddy._download_blast_binaries", return_value=False):
                assert not Sb._check_for_blast_bin("foo")
                out, err = capsys.readouterr()
                assert "Failed to download foo" in err

            assert not Sb._check_for_blast_bin("blastp")
            out, err = capsys.readouterr()
            assert "blastp downloaded" in err
            assert os.path.isfile("./blastp")
            os.remove("./blastp")


# ######################  '_download_blast_binaries' ###################### #
@pytest.mark.internet
@pytest.mark.slow
def test_dl_blast_bins():
    sys_memory = int(sys.maxsize)
    for platform in ["darwin", "linux32", "linux64", "win"]:
        sys.maxsize = 2000000000 if platform == "linux32" else sys_memory
        platform = "linux" if "linux" in platform else platform
        Sb._download_blast_binaries(ignore_pre_install=True, system=platform)
        if platform != "win":
            assert os.path.isfile('./blastdbcmd')
            assert os.path.isfile('./blastn')
            assert os.path.isfile('./blastp')
            os.remove("./blastdbcmd")
            os.remove("./blastn")
            os.remove("./blastp")
        else:
            assert os.path.isfile('./blastdbcmd.exe')
            assert os.path.isfile('./blastn.exe')
            assert os.path.isfile('./blastp.exe')
            os.remove("./blastdbcmd.exe")
            os.remove("./blastn.exe")
            os.remove("./blastp.exe")


# ######################  '_feature_rc' ###################### #
def test_feature_rc():
    tester = Sb._make_copy(sb_objects[1])
    seq1 = tester.records[0]
    hashed = ["dc71a33b64a766da8653c19f22fc4caa", "9ab8296fb3443198674d90abe3311ba6",
              "10018d1b15c7f76a6333ac3bf96d2d07", "273463b9eace12d2eeadbf272692d73e",
              "c452d66d13120cd6eb5f041b7c37dd27", "1811b0695dba1fc3fe431a6ee00ef359"]
    for feature in zip(seq1.features, hashed):
        assert string2hash(str(Sb._feature_rc(feature[0], 1203))) == feature[1]

    with pytest.raises(TypeError):
        feature = seq1.features[0]
        feature.location = {}
        Sb._feature_rc(feature, 1203)


# ######################  'guess_alphabet' ###################### #
def test_guess_alphabet():
    tester = Sb._make_copy(sb_objects[0])
    assert Sb._guess_alphabet(tester) == IUPAC.ambiguous_dna

    tester = Sb._make_copy(sb_objects[6])
    assert Sb._guess_alphabet(tester) == IUPAC.protein

    tester = Sb.SeqBuddy(resource("Mnemiopsis_rna.fa"))
    assert Sb._guess_alphabet(tester) == IUPAC.ambiguous_rna

    tester = Sb.SeqBuddy(">Seq1", in_format="fasta")
    assert not Sb._guess_alphabet(tester)


# ######################  'guess_format' ###################### #
def test_guess_format():
    assert Sb._guess_format(["foo", "bar"]) == "gb"
    assert Sb._guess_format(sb_objects[0]) == "fasta"
    assert Sb._guess_format(resource("Mnemiopsis_cds.fa")) == "fasta"
    assert Sb._guess_format(resource("blank.fa")) == "empty file"
    with pytest.raises(br.GuessError):
        Sb._guess_format("foo")

    temp_file = MyFuncs.TempFile()
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

    assert not Sb._guess_format(temp_file.path)


# ######################  '_stdout and _stderr' ###################### #
def test_stdout(capsys):
    Sb._stdout("Hello std_out", quiet=False)
    out, err = capsys.readouterr()
    assert out == "Hello std_out"

    Sb._stdout("Hello std_out", quiet=True)
    out, err = capsys.readouterr()
    assert out == ""


def test_stderr(capsys):
    Sb._stderr("Hello std_err", quiet=False)
    out, err = capsys.readouterr()
    assert err == "Hello std_err"

    Sb._stderr("Hello std_err", quiet=True)
    out, err = capsys.readouterr()
    assert err == ""


# ################################################ MAIN API FUNCTIONS ################################################ #
# ##################### '-ano', '--annotate' ###################### ##
def test_annotate_pattern():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), pattern='α4')
    assert seqs_to_hash(tester) == '527b2730d23d2ae4d708dbc55f560604'


def test_annotate_no_pattern():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100))
    assert seqs_to_hash(tester) == 'e8afc1fd9b0d5833972acfd2b719e42a'


def test_annotate_compoundlocation():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', [(1, 100), (200, 250)])
    assert seqs_to_hash(tester) == '4a294c1ea5783fe04e909d767e7dfce2'


def test_annotate_list_str():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', ['1-100', '(200-250)'])
    assert seqs_to_hash(tester) == '4a294c1ea5783fe04e909d767e7dfce2'


def test_annotate_str():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', '1-100, (200-250)')
    assert seqs_to_hash(tester) == '4a294c1ea5783fe04e909d767e7dfce2'


def test_annotate_fl_obj():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', FeatureLocation(start=0, end=100))
    assert seqs_to_hash(tester) == 'e8afc1fd9b0d5833972acfd2b719e42a'


def test_annotate_cl_obj():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', CompoundLocation([FeatureLocation(start=0, end=100),
                                                                   FeatureLocation(start=199, end=250)],
                                                                  operator='order'))
    assert seqs_to_hash(tester) == '4a294c1ea5783fe04e909d767e7dfce2'


def test_annotate_typerror():
    with pytest.raises(TypeError):
        tester = Sb._make_copy(sb_objects[1])
        Sb.annotate(tester, 'misc_feature', 5)


def test_annotate_pos_strand():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand='+')
    assert seqs_to_hash(tester) == 'e8afc1fd9b0d5833972acfd2b719e42a'


def test_annotate_neg_strand():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand='-')
    assert seqs_to_hash(tester) == '9de19b80180cdc65ca8a5b8e9c0e61e9'


def test_annotate_no_strand():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand=0)
    assert seqs_to_hash(tester) == 'e8afc1fd9b0d5833972acfd2b719e42a'


def test_annotate_qualifier_dict():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers={'foo': 'bar', 'hello': 'world'})
    assert seqs_to_hash(tester) == 'a948567058575491dbb3ed6ebee1f9b6'


def test_annotate_qualifier_list():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers=['foo=bar', 'hello=Blahh, =world'])
    assert seqs_to_hash(tester) == '5f8d9b4874754e50877126994141cc92'


def test_annotate_qualifier_error():
    tester = Sb._make_copy(sb_objects[1])
    with pytest.raises(TypeError):
        Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers=tuple)


def test_annotate_out_of_range():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', [(-10, 100), (200, 10000)])
    assert seqs_to_hash(tester) == '57f57b5e46c547542ea30e932ec4a42b'

    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', [(1, 10000)])
    assert seqs_to_hash(tester) == 'faa0e801376b3daba753890c18f0f029'

    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', [(1, 10000), (20000, 30000)])
    assert seqs_to_hash(tester) == 'faa0e801376b3daba753890c18f0f029'

    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', FeatureLocation(start=-10, end=100))
    assert seqs_to_hash(tester) == 'e8afc1fd9b0d5833972acfd2b719e42a'


def test_annotate_protein():
    tester = Sb._make_copy(sb_objects[7])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100))
    assert seqs_to_hash(tester) == '34f9601d72d5c37a8fd301ff90088cd8'


def test_annotate_unrec_strand(capsys):
    tester = Sb._make_copy(sb_objects[1])
    Sb.annotate(tester, 'misc_feature', (1, 100), strand='foo')
    out, err = capsys.readouterr()
    assert err == "Warning: strand input not recognized. Value set to None."


# ######################  '-asl', '--ave_seq_length' ###################### #
@pytest.mark.parametrize("seqbuddy", [Sb._make_copy(x) for x in sb_objects[0:6] if len(x.records) == 13])
def test_ave_seq_length_dna(seqbuddy):
    assert round(Sb.ave_seq_length(seqbuddy, clean=True), 2) == 1285.15


@pytest.mark.parametrize("seqbuddy", [Sb._make_copy(x) for x in sb_objects[6:12] if len(x.records) == 13])
def test_ave_seq_length_pep(seqbuddy):
    assert round(Sb.ave_seq_length(seqbuddy, clean=True), 2) == 427.38


# ######################  '-btr', '--back_translate' ###################### #
# Only fasta and genbank
hashes = ["1b14489a78bfe8255c777138877b9648", "b6bcb4e5104cb202db0ec4c9fc2eaed2",
          "859ecfb88095f51bfaee6a1d1abeb50f", "ba5c286b79a3514fba0b960ff81af25b",
          "952a91a4506afb57f27136aa1f2a8af9", "40c4a3e08c811b6bf3be8bedcb5d65a0",
          "3a3ee57f8dcde25c99a655494b218928", "bc9f1ec6ec92c30b5053cd9bb6bb6f53"]
organisms = ['human', 'human', 'yeast', 'yeast', 'ecoli', 'ecoli', 'mouse', 'mouse']
hashes = [(Sb._make_copy(sb_objects[sb_obj_indx]), organisms[indx], hashes[indx]) for indx, sb_obj_indx in
          enumerate([6, 7, 6, 7, 6, 7, 6, 7])]


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
        Sb.back_translate(Sb._make_copy(sb_objects[6]), 'fgsdjkghjdalgsdf', 'human')


def test_back_translate_bad_organism():
    seqbuddy = Sb._make_copy(sb_objects[6])
    with pytest.raises(AttributeError):
        Sb.back_translate(seqbuddy, 'OPTIMIZED', 'fgsdjkghjdalgsdf')


# ######################  '-bl2s', '--bl2seq' ###################### #
def test_bl2seq():
    seqbuddy = Sb._make_copy(sb_objects[0])
    result = Sb.bl2seq(seqbuddy)
    assert string2hash(str(result)) == '87dbd3baeb59285ad25e6473c87bb5bb'

    seqbuddy = Sb._make_copy(sb_objects[6])
    result = Sb.bl2seq(seqbuddy)
    assert string2hash(str(result)) == '248d4c53d7947c4c8dfd7c415bfbfbf2'


def test_bl2_no_binary():
    os.environ["PATH"] = ""
    with mock.patch('builtins.input', return_value="n"):
        with pytest.raises(RuntimeError):
            seqbuddy = Sb._make_copy(sb_objects[0])
            Sb.bl2seq(seqbuddy)

        with pytest.raises(RuntimeError):
            seqbuddy = Sb._make_copy(sb_objects[6])
            Sb.bl2seq(seqbuddy)
    os.environ["PATH"] = BACKUP_PATH


# ######################  '-bl', '--blast' ###################### #
def test_blastn():
    tester = Sb.pull_recs(Sb._make_copy(sb_objects[0]), '8', True)
    tester = Sb.blast(tester, blast_db=resource("blast/Mnemiopsis_cds.n"))
    assert seqs_to_hash(tester) == "95c417b6c2846d1b7a1a07f50c62ff8a"

    with pytest.raises(RuntimeError) as e:
        tester = Sb._make_copy(sb_objects[0])
        Sb.blast(tester, blast_db=resource("Mnemiopsis_cds.nhr"))
    assert "The .nhr file of your blast database was not found" in str(e.value)

    with mock.patch("SeqBuddy._check_for_blast_bin", return_value=False):
        with pytest.raises(SystemError) as e:
            Sb.blast(tester, blast_db=resource("blast/Mnemiopsis_cds.n"))
        assert 'blastn not found in system path' in str(e.value)

    tester = Sb.SeqBuddy(">Seq1\nATGCGCGCTACGCTAGCTAGCTAGCTCGCATGCAT")
    tester = Sb.blast(tester, blast_db=resource("blast/Mnemiopsis_cds.n"))
    assert len(tester.records) == 0


def test_blastp():
    seqbuddy = Sb.pull_recs(Sb.SeqBuddy(resource(seq_files[6])), '8', True)
    tester = Sb.blast(seqbuddy, blast_db=resource("blast/Mnemiopsis_pep.p"))
    assert seqs_to_hash(tester) == "4237c79672c1cf1d4a9bdb160a53a4b9"

    with pytest.raises(RuntimeError) as e:
        tester = Sb._make_copy(sb_objects[6])
        Sb.blast(tester, blast_db=resource("Mnemiopsis_cds.phr"))
    assert "The .phr file of your blast database was not found" in str(e.value)

    with mock.patch("SeqBuddy._check_for_blast_bin", return_value=False):
        with pytest.raises(SystemError) as e:
            Sb.blast(tester, blast_db=resource("blast/Mnemiopsis_pep.n"))
        assert 'blastp not found in system path' in str(e.value)


# ######################  '-cs', '--clean_seq'  ###################### #
def test_clean_seq():
    # Protein
    tester = Sb._make_copy(sb_objects[6])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert seqs_to_hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert seqs_to_hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"

    # DNA
    tester = Sb._make_copy(sb_objects[12])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert seqs_to_hash(tester) == "71b28ad2730a9849f2ba0f70e9e51a9f"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert seqs_to_hash(tester) == "1912fadb5ec52a38ec707c58085b86ad"
    tester = Sb._make_copy(sb_objects[12])
    tester = Sb.clean_seq(tester, ambiguous=False, rep_char="X")
    assert seqs_to_hash(tester) == "4c10ba4474d7484652cb633f03db1be1"

    # RNA
    tester = Sb._make_copy(sb_objects[13])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert seqs_to_hash(tester) == "cdb1b963536d57efc7b7f87d2bf4ad22"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert seqs_to_hash(tester) == "e66b785c649ad5086bcefd22e9ef9b41"
    tester = Sb._make_copy(sb_objects[13])
    tester = Sb.clean_seq(tester, ambiguous=False, rep_char="X")
    assert seqs_to_hash(tester) == "8ae19ab51b04076112d2f649353a4a79"

    # Alignment formats are converted to fasta to prevent errors with sequence lengths
    for tester in sb_objects[2:5]:
        tester = Sb.clean_seq(Sb._make_copy(tester))
        seqs_to_hash(tester) == "aa92396a9bb736ae6a669bdeaee36038"

# ######################  '-cmp', '--complement' ###################### #
hashes = ["e4a358ca57aca0bbd220dc6c04c88795", "3366fcc6ead8f1bba4a3650e21db4ec3",
          "365bf5d08657fc553315aa9a7f764286", "520036b49dd7c70b9dbf4ce4d2c0e1d8",
          "dc1f7a3769a1e0b007969db1ab405e89", "5891348e8659290c2355fabd0f3ba4f4"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_complement(seqbuddy, next_hash):
    tester = Sb.complement(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_complement_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.complement(sb_objects[6])


# ######################  '-cts', '--concat_seqs' ###################### #
# ToDo: Test the _clean parameter
hashes = ["2e46edb78e60a832a473397ebec3d187", "7421c27be7b41aeedea73ff41869ac47",
          "494988ffae2ef3072c1619eca8a0ff3b", "710cad348c5560446daf2c916ff3b3e4",
          "494988ffae2ef3072c1619eca8a0ff3b", "494988ffae2ef3072c1619eca8a0ff3b",
          "46741638cdf7abdf53c55f79738ee620", "8d0bb4e5004fb6a1a0261c30415746b5",
          "2651271d7668081cde8012db4f9a6574", "7846b2d080f09b60efc6ee43cd6d8502",
          "5d1f8db03d6be30a7d77b00a0fba0b43", "2651271d7668081cde8012db4f9a6574"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_concat_seqs(seqbuddy, next_hash):
    tester = Sb.concat_seqs(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_concat_seqs_clean():
    tester = Sb.concat_seqs(Sb._make_copy(sb_objects[2]), clean=True)
    assert seqs_to_hash(tester) == "2e46edb78e60a832a473397ebec3d187"


# ##################### '-cc', 'count_codons' ###################### ##
def test_count_codons_dna():
    tester = Sb._make_copy(sb_objects[0])
    tester.records[0].seq = tester.records[0].seq[:-4]
    counter = Sb.count_codons(tester)
    assert seqs_to_hash(counter[0]) == '4dbd9a5c68d85bb200c75b309fdaeeca'
    assert string2hash(str(counter[1])) == '8d47313f3db02aee48575ff8ff4741b4'


def test_count_codons_rna():
    tester = Sb.count_codons(Sb.dna2rna(Sb._make_copy(sb_objects[0])))[1]
    assert md5(str(tester).encode()).hexdigest() == 'b91daa8905533b5885d2067d9d6ffe36'


def test_count_codons_dna_badchar():
    tester = Sb.count_codons(Sb.insert_sequence(Sb._make_copy(sb_objects[0]), 'PPP', -1))[1]
    assert md5(str(tester).encode()).hexdigest() == '9aba116675fe0e9eaaf43e5c6e0ba99d'


def test_count_codons_pep_exception():
    tester = Sb._make_copy(sb_objects[6])
    with pytest.raises(TypeError):
        Sb.count_codons(tester)


# ######################  '-cr', '--count_residues' ###################### #
def test_count_residues():
    # Unambiguous DNA
    tester = Sb.SeqBuddy(">seq1\nACGCGAAGCGAACGCGCAGACGACGCGACGACGACGACGCA", in_format="fasta")
    tester = Sb.count_residues(tester).to_dict()
    assert tester["seq1"].buddy_data["res_count"]['A'] == [13, 0.3170731707317073]

    tester = Sb._make_copy(sb_objects[0])
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert res_count['G'] == [265, 0.21703521703521703]
    assert "% Ambiguous" not in res_count and "U" not in res_count

    # Unambiguous RNA
    Sb.dna2rna(tester)
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert res_count['U'] == [356, 0.2915642915642916]
    assert "% Ambiguous" not in res_count and "T" not in res_count

    # Ambiguous DNA
    tester = Sb._make_copy(sb_objects[12])
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert "U" not in res_count
    assert res_count['Y'] == [1, 0.000819000819000819]
    assert res_count['% Ambiguous'] == 0.98

    # Ambiguous RNA
    tester = Sb._make_copy(sb_objects[13])
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert "T" not in res_count
    assert res_count['U'] == [353, 0.2891072891072891]
    assert res_count["Y"] == [1, 0.000819000819000819]
    assert res_count['% Ambiguous'] == 0.98

    # Protein
    tester = Sb._make_copy(sb_objects[6])
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert res_count['P'] == [17, 0.04176904176904177]
    assert res_count["G"] == [23, 0.056511056511056514]
    assert "% Ambiguous" not in res_count
    res_count = tester.to_dict()['Mle-Panxα8'].buddy_data["res_count"]
    assert res_count["% Ambiguous"] == 1.2
    assert res_count["% Positive"] == 12.23
    assert res_count["% Negative"] == 12.71
    assert res_count["% Uncharged"] == 73.62
    assert res_count["% Hyrdophilic"] == 36.93
    assert res_count["% Hyrdophobic"] == 55.4


# ######################  '-dgn' '--degenerate_sequence'################### #
def test_degenerate_sequence_without_arguments():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Sb.degenerate_sequence(tester)
    assert seqs_to_hash(tester) == '0638bc6546eebd9d50f771367d6d7855'

dgn_hashes = ['0638bc6546eebd9d50f771367d6d7855', '72373f8356051e2c6b67642451379054',
              '9172ad5947c0961b54dc5adbd03d4249', 'b45ac94ee6a98e495e115bfeb5bd9bcd',
              '76c45b4de8f7527b4139446b4551712b', 'baa5b48938cc5cae953c9083a5b21b12',
              '0ca67c4740fefbc7a20d806715c3ca12', 'd43ad8f328ff1d30eb1fb7bcd667a345',
              'd9d0f5cd8f0c25a0042527cc1cea802e', '4b9790f3f4eeeae1a9667b62b93bc961',
              '7ec4365c3571813d63cee4b70ba5dcf5']

hashes = [(Sb._make_copy(sb_objects[0]), dgn_hash, indx + 1) for indx, dgn_hash in enumerate(dgn_hashes)]


@pytest.mark.parametrize("seqbuddy, dgn_hash, tables", hashes)
def test_degenerate_sequence_with_different_codon_tables(seqbuddy, dgn_hash, tables):
    tester = Sb.degenerate_sequence(seqbuddy, table=tables)
    assert seqs_to_hash(tester) == dgn_hash


def test_degenerate_sequence_edges():
    tester = Sb._make_copy(sb_objects[6])

    # Bad table reference
    with pytest.raises(KeyError) as e:
        Sb.degenerate_sequence(tester, 100)
    assert "Could not locate codon dictionary" in str(e.value)

    # Protein input
    with pytest.raises(TypeError) as e:
        Sb.degenerate_sequence(tester, 1)
    assert "Nucleic acid sequence required, not protein." in str(e.value)

    # RNA input
    tester = Sb.degenerate_sequence(Sb._make_copy(sb_objects[13]))
    assert seqs_to_hash(tester) == "ff52e05971aeafd24c73a3b543901e4b"


# ######################  '-df', '--delete_features' ###################### #
def test_delete_features():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.delete_features(tester, 'donor')
    assert seqs_to_hash(tester) == 'f84df6a77063c7def13babfaa0555bbf'


# ######################  '-dl', '--delete_large' ###################### #
def test_delete_large():
    tester = Sb._make_copy(sb_objects[0])
    tester = Sb.delete_large(tester, 1285)
    assert seqs_to_hash(tester) == '25859dc69d46651a1e04a70c07741b35'


# ######################  '-dm', '--delete_metadata' ###################### #
hashes = ["aa92396a9bb736ae6a669bdeaee36038", "544ab887248a398d6dd1aab513bae5b1", "cb1169c2dd357771a97a02ae2160935d",
          "503e23720beea201f8fadf5dabda75e4", "52c23bd793c9761b7c0f897d3d757c12", "a50943ccd028b6f5fa658178fa8cf54d",
          "bac5dc724b1fee092efccd2845ff2513", "858e8475f7bc6e6a24681083a8635ef9", "17ff1b919cac899c5f918ce8d71904f6",
          "968ed9fa772e65750f201000d7da670f", "ce423d5b99d5917fbef6f3b47df40513", "e224c16f6c27267b5f104c827e78df33"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_delete_metadata(seqbuddy, next_hash):
    tester = Sb.delete_metadata(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


# ######################  '-dr', '--delete_records' ###################### #
dr_hashes = ["54bdb42423b1d331acea18218101e5fc", "e2c03f1fa21fd27b2ff55f7f721a1a99", "6bc8a9409b1ef38e4f6f12121368883e",
             "bda7be10061b0dcaeb66bebe3d736fee", "6e2fce2239e2669b23f290049f87fbc4", "4c97144c5337f8a40c4fd494e622bf0d"]
dr_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(dr_hashes)]


@pytest.mark.parametrize("seqbuddy, next_hash", dr_hashes)
def test_delete_records(seqbuddy, next_hash):
    tester = Sb.delete_records(seqbuddy, 'α2')
    assert seqs_to_hash(tester) == next_hash


def test_delete_records2():
    tester = Sb.delete_records(Sb._make_copy(sb_objects[0]), ['α1', 'α2'])
    assert seqs_to_hash(tester) == "eca4f181dae3d7998464ff71e277128f"

    with pytest.raises(ValueError) as e:
        Sb.delete_records(Sb._make_copy(sb_objects[0]), dict)
    assert "'patterns' must be a list or a string." in str(e.value)


# #####################  '-drp', '--delete_repeats' ###################### ##
def test_delete_repeats():
    tester = Sb.SeqBuddy(resource("Duplicate_seqs.fa"))
    tester = Sb.delete_repeats(tester)
    tester = Sb.find_repeats(tester)
    assert len(tester.repeat_ids) == 0
    assert len(tester.repeat_seqs) == 0


# ######################  '-ds', '--delete_small' ###################### #
def test_delete_small():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Sb.delete_small(tester, 1285)
    assert seqs_to_hash(tester) == '196adf08d4993c51050289e5167dacdf'


# ######################  '-d2r', '--transcribe' and 'r2d', '--back_transcribe' ###################### #
d2r_hashes = ["d2db9b02485e80323c487c1dd6f1425b", "9ef3a2311a80f05f21b289ff7f401fff",
              "f3bd73151645359af5db50d2bdb6a33d", "e55bd18b6d82a7fc3150338173e57e6a",
              "a083e03b4e0242fa3c23afa80424d670", "45b511f34653e3b984e412182edee3ca"]
r2d_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "2e02a8e079267bd9add3c39f759b252c",
              "cb1169c2dd357771a97a02ae2160935d", "503e23720beea201f8fadf5dabda75e4",
              "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b"]

hashes = [(Sb._make_copy(sb_objects[indx]), d2r_hash, r2d_hashes[indx]) for indx, d2r_hash in enumerate(d2r_hashes)]


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


# ######################  '-er', '--extract_range' ###################### #
hashes = ["201235ed91ad0ed9a7021136487fed94", "3e791c6a6683516aff9572c24f38f0b3", "4063ab66ced2fafb080ceba88965d2bb",
          "33e6347792aead3c454bac0e05a292c6", "9a5c491aa293c6cedd48c4c249d55aff", "d724df01ae688bfac4c6dfdc90027440",
          "904a188282f19599a78a9d7af4169de6", "b8413624b9e684a14fc9f398a62e3965", "6a27222d8f60ee8496cbe0c41648a116",
          "c9a1dd913190f95bba5eca6a89685c75", "6f579144a43dace285356ce6eb326d3b", "38d571c9681b4fa420e3d8b54c507f9c"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_extract_range(seqbuddy, next_hash):
    tester = Sb.extract_range(seqbuddy, 50, 300)
    assert seqs_to_hash(tester) == next_hash


def test_extract_range_end_less_than_start():
    with pytest.raises(ValueError):
        Sb.extract_range(sb_objects[0], 500, 50)


# #####################  '-fcpg', '--find_CpG' ###################### ##
def test_find_cpg():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.gb"))
    tester = Sb.find_cpg(tester)
    assert seqs_to_hash(tester) == "9499f524da0c35a60502031e94864928"


# #####################  '-fp', '--find_pattern' ###################### ##
def test_find_pattern():
    tester = Sb.find_pattern(Sb._make_copy(sb_objects[1]), "ATGGT")
    assert seqs_to_hash(tester) == "a2e110d3ba14ea5d1b236a8abef7341a"
    tester = Sb.find_pattern(Sb._make_copy(sb_objects[1]), "ATg{2}T")
    assert seqs_to_hash(tester) == "c8a3e2e8b97f47c5cc41f04a44243e34"
    tester = Sb.find_pattern(Sb._make_copy(sb_objects[1]), "ATg{2}T", "tga.{1,6}tg")
    assert seqs_to_hash(tester) == "0e11b2c0e9451fbfcbe39e3b5be2cf60"


# #####################  '-frp', '--find_repeats' ###################### ##
def test_find_repeats():
    tester = Sb.SeqBuddy(resource("Duplicate_seqs.fa"))
    tester.unique_seqs, tester.repeat_ids, tester.repeat_seqs = {}, {}, {}
    Sb.find_repeats(tester)
    assert 'Seq1' in tester.unique_seqs
    assert 'Seq12' in tester.repeat_ids

    for key in tester.repeat_seqs:
        assert 'Seq12' in tester.repeat_seqs[key] or 'Seq10A' in tester.repeat_seqs[key]


# ######################  '-frs', '--find_restriction_sites' ###################### #
def test_restriction_sites(capsys):
    # No arguments passed in = commercial REs and any number of cut sites
    tester = Sb.find_restriction_sites(Sb._make_copy(sb_objects[1]))
    assert seqs_to_hash(tester) == 'fce4c8bee040d5ea6fa4bf9985f7310f'
    assert md5(str(tester.restriction_sites).encode()).hexdigest() == "741ca6ca9204a067dce7398f15c6e350"

    # All enzymes
    tester = Sb.find_restriction_sites(Sb._make_copy(sb_objects[1]), enzyme_group=["all"])
    assert seqs_to_hash(tester) == '254e854862142c54a53079ae224b4180'
    assert md5(str(tester.restriction_sites).encode()).hexdigest() == "0a5012b1a3d83f719cdd8cecf48258f4"

    # Specify a few REs and limit the number of cuts
    tester = Sb.find_restriction_sites(Sb._make_copy(sb_objects[1]), min_cuts=2, max_cuts=4,
                                       enzyme_group=["EcoRI", "KspI", "TasI", "Bme1390I", "FooBR"])
    out, err = capsys.readouterr()
    assert seqs_to_hash(tester) == 'c42b3bf0367557383000b897432fed2d'
    assert md5(str(tester.restriction_sites).encode()).hexdigest() == "0d2e5fdba6fed434495481397a91e56a"
    assert "Warning: FooBR not a known enzyme" in err

    with pytest.raises(TypeError) as e:
        Sb.find_restriction_sites(sb_objects[7])
    assert str(e.value) == "Unable to identify restriction sites in protein sequences."

    with pytest.raises(ValueError) as e:
        Sb.find_restriction_sites(tester, min_cuts=4, max_cuts=2)
    assert str(e.value) == "min_cuts parameter has been set higher than max_cuts."

    # 2-cutters and non-cutters
    Sb.find_restriction_sites(tester, enzyme_group=["AjuI", "AlwFI"])
    out, err = capsys.readouterr()
    assert "Warning: Double-cutters not supported." in err
    assert "Warning: No-cutters not supported." in err


# ######################  '-hsi', '--hash_sequence_ids' ###################### #
def test_hash_seq_ids():
    tester = Sb.SeqBuddy(Sb._make_copy(sb_objects[0]))
    Sb.hash_sequence_ids(tester)
    assert len(tester.records[0].id) == 10

    tester = Sb.hash_sequence_ids(tester, 25)
    assert len(tester.records[0].id) == 25
    assert len(tester.hash_map) == 13


def test_hash_seq_ids_errors():
    tester = Sb.SeqBuddy(Sb._make_copy(sb_objects[0]))
    with pytest.raises(TypeError) as e:
        Sb.hash_sequence_ids(tester, "foo")
    assert str(e.value) == "Hash length argument must be an integer, not <class 'str'>"

    with pytest.raises(ValueError) as e:
        Sb.hash_sequence_ids(tester, 0)
    assert str(e.value) == "Hash length must be greater than 0"

    tester.records *= 10
    with pytest.raises(ValueError) as e:
        Sb.hash_sequence_ids(tester, 1)
    assert "Insufficient number of hashes available to cover all sequences." in str(e.value)


# ##################### '-is', 'insert_seq' ###################### ##
def test_insert_seqs_start():
    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA')) == 'f65fee08b892af5ef93caa1bf3cb3980'

    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', -9000)) == 'f65fee08b892af5ef93caa1bf3cb3980'

    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', -1)) == '792397e2e32e95b56ddc15b8b2310ec0'

    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', 9000)) == '792397e2e32e95b56ddc15b8b2310ec0'

    tester = Sb._make_copy(sb_objects[0])
    Sb.insert_sequence(tester, 'AACAGGTCGAGCA', 100, ["α[23]", "α5"])
    assert seqs_to_hash(tester) == 'edcd7934eb026ac3ea4b603ac85ca79f'

    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', -25)) == '29cab1e72ba95572c3aec469270071e9'


# ######################  '-ip', '--isoelectric_point' ###################### #
@pytest.mark.parametrize("seqbuddy", [Sb._make_copy(x) for x in [sb_objects[6], sb_objects[7],
                                                                 sb_objects[8], sb_objects[10], sb_objects[11]]])
def test_isoelectric_point(seqbuddy):
    tester = Sb.isoelectric_point(Sb.clean_seq(seqbuddy))
    tester = tester.to_dict()
    assert tester["Mle-Panxα12"].features[-1].qualifiers["value"] == 6.0117797852
    if seqbuddy.out_format == "gb":
        assert seqs_to_hash(seqbuddy) == "8bc299e31f436d192bf8cf8b7af671a8"

    with pytest.raises(TypeError):
        Sb.isoelectric_point(Sb._make_copy(sb_objects[0]))


# ######################  '-lc', '--lowercase' and 'uc', '--uppercase'  ###################### #
uc_hashes = ["25073539df4a982b7f99c72dd280bb8f", "2e02a8e079267bd9add3c39f759b252c", "52e74a09c305d031fc5263d1751e265d",
             "cfe6cb9c80aebd353cf0378a6d284239", "6e5542f41d17ff33afb530b4d07408a3", "b82538a4630810c004dc8a4c2d5165ce",
             "c10d136c93f41db280933d5b3468f187", "7a8e25892dada7eb45e48852cbb6b63d", "8b6737fe33058121fd99d2deee2f9a76",
             "968ed9fa772e65750f201000d7da670f", "ce423d5b99d5917fbef6f3b47df40513", "f35cbc6e929c51481e4ec31e95671638"]

lc_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "2e02a8e079267bd9add3c39f759b252c", "cb1169c2dd357771a97a02ae2160935d",
             "503e23720beea201f8fadf5dabda75e4", "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
             "14227e77440e75dd3fbec477f6fd8bdc", "7a8e25892dada7eb45e48852cbb6b63d", "17ff1b919cac899c5f918ce8d71904f6",
             "aacda2f5d4077f23926400f74afa2f46", "e3dc2e0347f40fffec45d053f4f34c96", "c0dce60745515b31a27de1f919083fe9"]

hashes = [(Sb._make_copy(sb_objects[indx]), uc_hash, lc_hashes[indx]) for indx, uc_hash in enumerate(uc_hashes)]


@pytest.mark.parametrize("seqbuddy,uc_hash,lc_hash", hashes)
def test_cases(seqbuddy, uc_hash, lc_hash):  # NOTE: Biopython always writes genbank to spec in lower case
    tester = Sb.uppercase(seqbuddy)
    assert seqs_to_hash(tester) == uc_hash
    tester = Sb.lowercase(tester)
    assert seqs_to_hash(tester) == lc_hash


# ######################  '-mui', '--make_ids_unique' ###################### #
def test_make_ids_unique():
    tester = Sb.SeqBuddy(resource("Duplicate_seqs.fa"))
    Sb.make_ids_unique(tester)
    assert seqs_to_hash(tester) == "3f84819b81185d49a7357947ca7f0707"


# ######################  '-fn2p', '--map_features_nucl2prot' ###################### #
# Map the genbank DNA file to all protein files, and the fasta DNA file to fasta protein
hashes = ["5216ef85afec36d5282578458a41169a", "a8f7c129cf57a746c20198bf0a6b9cf4", "bb0c9da494b5418fb87862dab2a66cfa",
          "3c0e3ec45abd774813a274fda1b4a5f2", "a7f6c4bb410f17cfc3e8966ccbe3e065", "1b8c44f4ace877b568c0915033980bed"]
prot_indx = [6, 7, 8, 9, 10, 11]
hashes = [(Sb._make_copy(sb_objects[1]), Sb._make_copy(sb_objects[prot_indx[indx]]), value)
          for indx, value in enumerate(hashes)]
hashes.append((Sb._make_copy(sb_objects[0]), Sb._make_copy(sb_objects[6]), "854566b485af0f277294bbfb15f7dd0a"))


@pytest.mark.parametrize("_dna,_prot,next_hash", hashes)
def test_map_features_nucl2prot(_dna, _prot, next_hash):
    tester = Sb.map_features_nucl2prot(_dna, _prot)
    tester.out_format = "gb"
    assert seqs_to_hash(tester) == next_hash


def test_map_features_nucl2prot_2(capsys):
    tester = Sb._make_copy(sb_objects[1])
    Sb.pull_record_ends(tester, 1300)
    tester = Sb.annotate(tester, "foo", [(20, 40), (50, 60)], pattern="α9")
    mapped = Sb.map_features_nucl2prot(Sb._make_copy(tester), Sb._make_copy(sb_objects[6]))
    mapped.out_format = "gb"
    assert seqs_to_hash(mapped) == "807025489aadf98d851501f49d463e4a"
    out, err = capsys.readouterr()
    assert string2hash(err) == "6c840e4acaaf4328672ca164f854000a"

    prot_tester = Sb._make_copy(sb_objects[7])
    prot_tester.records = sorted(prot_tester.records, key=lambda x: x.id)
    dna_tester = Sb._make_copy(sb_objects[1])
    dna_tester.records = sorted(dna_tester.records, key=lambda x: x.id)
    Sb.rename(dna_tester, "α4", "A4")

    Sb.map_features_nucl2prot(Sb._make_copy(dna_tester), prot_tester, mode="list")
    assert seqs_to_hash(prot_tester) == "9fdb606ea65d6c050540a94137ae6e0d"

    Sb.map_features_nucl2prot(Sb._make_copy(dna_tester), prot_tester, mode="key")
    assert seqs_to_hash(prot_tester) == "9fdb606ea65d6c050540a94137ae6e0d"
    out, err = capsys.readouterr()
    assert string2hash(err) == "6fd2b5f2a7a3995d3f49c4919c3358b0"

    mock_obj = type('test', (object,), {})()
    mock_obj.start = 2
    mock_obj.end = 6
    dna_tester.records[0].features[0].location = mock_obj
    with pytest.raises(TypeError) as e:
        Sb.map_features_nucl2prot(Sb._make_copy(dna_tester), Sb._make_copy(sb_objects[6]))
    assert "FeatureLocation or CompoundLocation object required." in str(e.value)

    tester.records[0].features[0].location = tester.records[1].features[0].location
    with pytest.raises(ValueError) as e:
        Sb.map_features_nucl2prot(Sb._make_copy(tester), Sb._make_copy(sb_objects[6]), mode="foo")
    assert "'mode' must be either 'key' or 'position'" in str(e.value)

    with pytest.raises(ValueError) as e:
        Sb.pull_recs(tester, "α[1-8]")
        Sb.map_features_nucl2prot(Sb._make_copy(tester), Sb._make_copy(sb_objects[6]), mode="list")
    assert "The two input files do not contain the same number of sequences" in str(e.value)

# ######################  '-fp2n', '--map_features_prot2nucl' ###################### #
hashes = ["3ebc92ca11505489cab2453d2ebdfcf2", "feceaf5e17935afb100b4b6030e27fee",
          "bfd36942768cf65c473b3aaebb83e4fa", "9ba4af4e5dd0bf4a445d173604b92996",
          "c178763aa9596e341bbbc088f1f791c9", "84cc7ecb54603c5032737e5263a52bd3"]

hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("_dna,next_hash", hashes)
def test_map_features_prot2nucl(_dna, next_hash):
    tester = Sb.map_features_prot2nucl(Sb._make_copy(sb_objects[7]), _dna)
    tester.out_format = "gb"
    assert seqs_to_hash(tester) == next_hash


def test_map_features_prot2nucl_2(capsys):
    prot_tester = Sb._make_copy(sb_objects[7])
    Sb.pull_record_ends(prot_tester, 450)
    prot_tester = Sb.annotate(prot_tester, "foo", [(20, 40), (50, 60)], pattern="α9")
    mapped = Sb.map_features_prot2nucl(Sb._make_copy(prot_tester), Sb._make_copy(sb_objects[0]))
    mapped.out_format = "gb"
    assert seqs_to_hash(mapped) == "552b31b8a7068d12d6c55c7e5d293c54"
    out, err = capsys.readouterr()
    assert err == "Warning: size mismatch between aa and nucl seqs for Mle-Panxα7A --> 450, 1875\n"

    prot_tester = Sb._make_copy(sb_objects[7])
    prot_tester.records = sorted(prot_tester.records, key=lambda x: x.id)
    dna_tester = Sb._make_copy(sb_objects[1])
    dna_tester.records = sorted(dna_tester.records, key=lambda x: x.id)
    Sb.rename(prot_tester, "α4", "A4")
    Sb.map_features_prot2nucl(Sb._make_copy(prot_tester), dna_tester, mode="list")
    assert seqs_to_hash(dna_tester) == "c32f43cc5205867c0eb1d3873e27319b"

    Sb.map_features_prot2nucl(Sb._make_copy(prot_tester), dna_tester, mode="key")
    assert seqs_to_hash(dna_tester) == "c32f43cc5205867c0eb1d3873e27319b"
    out, err = capsys.readouterr()
    assert string2hash(err) == "c1f51587b07d2cc6156b8aac07384834"

    mock_obj = type('test', (object,), {})()
    mock_obj.start = 2
    mock_obj.end = 6
    prot_tester.records[0].features[0].location = mock_obj
    with pytest.raises(TypeError) as e:
        Sb.map_features_prot2nucl(Sb._make_copy(prot_tester), Sb._make_copy(sb_objects[0]))
    assert "FeatureLocation or CompoundLocation object required." in str(e.value)

    prot_tester.records[0].features[0].location = prot_tester.records[1].features[0].location
    with pytest.raises(ValueError) as e:
        Sb.map_features_prot2nucl(Sb._make_copy(prot_tester), Sb._make_copy(sb_objects[6]), mode="foo")
    assert "'mode' must be either 'key' or 'position'" in str(e.value)

    with pytest.raises(ValueError) as e:
        Sb.pull_recs(prot_tester, "α[1-8]")
        Sb.map_features_prot2nucl(Sb._make_copy(prot_tester), Sb._make_copy(sb_objects[6]), mode="list")
    assert "The two input files do not contain the same number of sequences" in str(e.value)


# #####################  '-mg', '--merge' ###################### ##
def test_merge():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds_dummy_features.gb"))
    tester = Sb.merge(tester, Sb._make_copy(sb_objects[1]))
    assert seqs_to_hash(tester) == "5fb43a869bb9f766625d1336081766fc"

    tester.records[0].seq = tester.records[1].seq
    with pytest.raises(RuntimeError) as e:
        Sb.merge(tester, Sb._make_copy(sb_objects[1]))
    assert "Record mismatch: ID" in str(e.value)


# ######################  '-mw', '--molecular_weight' ###################### #
def test_molecular_weight():
    # Unambiguous DNA
    tester = Sb.molecular_weight(Sb._make_copy(sb_objects[1]))
    assert tester.molecular_weights['masses_ds'][0] == 743477.1
    assert tester.molecular_weights['masses_ss'][0] == 371242.6
    assert seqs_to_hash(tester) == "e080cffef0ec6c5e8eada6f57bbc35f9"
    # Ambiguous DNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("ambiguous_dna.fa")))
    assert tester.molecular_weights['masses_ds'][0] == 743477.08
    assert tester.molecular_weights['masses_ss'][0] == 371202.59
    # Unambiguous RNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("Mnemiopsis_rna.fa")))
    assert tester.molecular_weights['masses_ss'][0] == 387372.6
    # Ambiguous RNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("ambiguous_rna.fa")))
    assert tester.molecular_weights['masses_ss'][0] == 387371.6
    # Protein
    tester = Sb.molecular_weight(Sb._make_copy(sb_objects[7]))
    assert tester.molecular_weights['masses_ss'][0] == 45692.99
    assert seqs_to_hash(tester) == "fb1a66b7eb576c0584fc7988c45b6a18"


# ######################  '-ns', '--num_seqs' ###################### #
seq_counts = [(sb_objects[0], 13), (sb_objects[1], 13), (sb_objects[2], 13), (sb_objects[3], 8),
              (sb_objects[4], 13), (sb_objects[5], 13), (sb_objects[6], 13), (sb_objects[9], 8)]


@pytest.mark.parametrize("seqbuddy, num", seq_counts)
def test_num_seqs(seqbuddy, num):
    assert Sb.num_seqs(seqbuddy) == num


def test_empty_file():
    tester = Sb.SeqBuddy(resource("blank.fa"))
    assert type(tester) == Sb.SeqBuddy
    assert len(tester.records) == 0

# ######################  '-ofa', '--order_features_alphabetically' ###################### #
fwd_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "21547b4b35e49fa37e5c5b858808befb",
              "cb1169c2dd357771a97a02ae2160935d", "503e23720beea201f8fadf5dabda75e4",
              "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "d0297078b4c480a49b6da5b719310d0e",
              "17ff1b919cac899c5f918ce8d71904f6", "968ed9fa772e65750f201000d7da670f",
              "ce423d5b99d5917fbef6f3b47df40513", "c0dce60745515b31a27de1f919083fe9"]

rev_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "3b718ec3cb794bcb658d900e517110cc",
              "cb1169c2dd357771a97a02ae2160935d", "503e23720beea201f8fadf5dabda75e4",
              "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "c6a788d8ea916964605ac2942c459c9b",
              "17ff1b919cac899c5f918ce8d71904f6", "968ed9fa772e65750f201000d7da670f",
              "ce423d5b99d5917fbef6f3b47df40513", "c0dce60745515b31a27de1f919083fe9"]
hashes = [(Sb._make_copy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)
def test_order_features_alphabetically(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_features_alphabetically(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_features_alphabetically(seqbuddy, reverse=True)
    assert seqs_to_hash(tester) == rev_hash


# ######################  '-ofp', '--order_features_by_position' ###################### #
fwd_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "2e02a8e079267bd9add3c39f759b252c",
              "cb1169c2dd357771a97a02ae2160935d", "503e23720beea201f8fadf5dabda75e4",
              "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "7a8e25892dada7eb45e48852cbb6b63d",
              "17ff1b919cac899c5f918ce8d71904f6", "968ed9fa772e65750f201000d7da670f",
              "ce423d5b99d5917fbef6f3b47df40513", "c0dce60745515b31a27de1f919083fe9"]

rev_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "4345a14fe27570b3c837c30a8cb55ea9",
              "cb1169c2dd357771a97a02ae2160935d", "503e23720beea201f8fadf5dabda75e4",
              "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "9e7c2571db1386bba5983365ae235e1b",
              "17ff1b919cac899c5f918ce8d71904f6", "968ed9fa772e65750f201000d7da670f",
              "ce423d5b99d5917fbef6f3b47df40513", "c0dce60745515b31a27de1f919083fe9"]
hashes = [(Sb._make_copy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)
def test_order_features_position(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_features_by_position(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_features_by_position(seqbuddy, reverse=True)
    assert seqs_to_hash(tester) == rev_hash


# ######################  '-oi', '--order_ids' ###################### #
fwd_hashes = ["ccc59629741421fb717b9b2403614c62", "2f9bc0dd9d79fd8160a621280be0b0aa",
              "60bbc6306cbb4eb903b1212718bb4592", "3c49bdc1b0fe4e1d6bfc148eb0293e21",
              "54b25170fbda0845b1537a674d37cbd4", "4078182a81382b815528fdd5c158fbec"]
rev_hashes = ["503a71fc2e8d143361cbe8f4611527fd", "dd269961d4d5301d1bf87e0093568851",
              "82fea6e3d3615ac75ec5022abce255da", "d6e79a5faeaff396aa7eab0b460c3eb9",
              "269e151ec39dd5b0643fea5f801150a6", "3b83a3c73a6cdded6635ffa10c4a16e1"]

hashes = [(Sb._make_copy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)
def test_order_ids(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_ids(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_ids(seqbuddy, reverse=True)
    assert seqs_to_hash(tester) == rev_hash


# ######################  '-oir', '--order_ids_randomly' ###################### #
@pytest.mark.parametrize("seqbuddy", [Sb._make_copy(x) for x in sb_objects])
def test_order_ids_randomly(seqbuddy):
    tester = Sb.order_ids_randomly(Sb._make_copy(seqbuddy))
    assert seqs_to_hash(seqbuddy) != seqs_to_hash(tester)
    assert seqs_to_hash(Sb.order_ids(tester)) == seqs_to_hash(tester)


def test_order_ids_randomly2():
    tester = Sb._make_copy(sb_objects[0])
    for _ in range(15):  # This will fail to repeat the while loop only ~5% of the time
        Sb.pull_recs(tester, "α[789]")
        assert seqs_to_hash(tester) != seqs_to_hash(Sb.order_ids_randomly(tester))

    Sb.pull_recs(tester, "α[89]")
    assert seqs_to_hash(tester) != seqs_to_hash(Sb.order_ids_randomly(tester))

    Sb.pull_recs(tester, "α[9]")
    assert seqs_to_hash(tester) == seqs_to_hash(Sb.order_ids_randomly(tester))

    tester = Sb.SeqBuddy(tester.records * 3)
    assert seqs_to_hash(tester) == seqs_to_hash(Sb.order_ids_randomly(tester))


# #####################  '-prr', '--pull_random_recs' ###################### ##
@pytest.mark.parametrize("seqbuddy", sb_objects)
def test_pull_random_recs(seqbuddy):
    tester = Sb.pull_random_recs(Sb._make_copy(seqbuddy))
    orig_seqs = tester.to_dict()
    assert len(tester.records) == 1
    assert tester.records[0].id in orig_seqs


# #####################  '-pre', '--pull_record_ends' ###################### ##
def test_pull_record_ends():
    tester = Sb.pull_record_ends(Sb._make_copy(sb_objects[1]), 10)
    assert seqs_to_hash(tester) == 'd46867e4ca7a9f474c45473fc3495413'

    tester = Sb.pull_record_ends(Sb._make_copy(sb_objects[1]), 2000)
    assert seqs_to_hash(tester) == '908744b00d9f3392a64b4b18f0db9fee'

    tester = Sb.pull_record_ends(Sb._make_copy(sb_objects[1]), -10)
    assert seqs_to_hash(tester) == 'd7970570d65872993df8a3e1d80f9ff5'

    tester = Sb.pull_record_ends(Sb._make_copy(sb_objects[1]), -2000)
    assert seqs_to_hash(tester) == '908744b00d9f3392a64b4b18f0db9fee'

    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    with pytest.raises(ValueError):
        Sb.pull_record_ends(Sb._make_copy(seqbuddy), 'foo')


# ######################  '-pr', '--pull_records' ###################### #
pr_hashes = ["5b4154c2662b66d18776cdff5af89fc0", "e196fdc5765ba2c47f97807bafb6768c", "bc7dbc612bc8139eba58bf896b7eaf2f",
             "7bb4aac2bf50381ef1d27d82b7dd5a53", "72328eddc751fd79406bb911dafa57a2", "b006b40ff17ba739929448ae2f9133a6"]
pr_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(pr_hashes)]


@pytest.mark.parametrize("seqbuddy, next_hash", pr_hashes)
def test_pull_recs(seqbuddy, next_hash):
    tester = Sb.pull_recs(seqbuddy, 'α2')
    assert seqs_to_hash(tester) == next_hash


# #####################  '-prg', '--purge' ###################### ##
def test_purge():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_pep.fa"))
    Sb.purge(tester, 200)
    assert seqs_to_hash(tester) == 'b21b2e2f0ca1fcd7b25efbbe9c08858c'


# ######################  '-ri', '--rename_ids' ###################### #
hashes = ["8b4a9e3d3bb58cf8530ee18b9df67ff1", "144cc5ed20678a818bce908c475ae450", "243024bfd2f686e6a6e0ef65aa963494",
          "98bb9b57f97555d863054ddb526055b4", "2443c47a712f19099e94fc015dc980a9", "65196fd4f2a4e339e1545f6ed2a6acc3"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_rename_ids(seqbuddy, next_hash):
    tester = Sb.rename(seqbuddy, 'Panx', 'Test', 0)
    assert seqs_to_hash(tester) == next_hash


def test_rename_ids2():
    tester = Sb._make_copy(sb_objects[0])
    Sb.rename(tester, "Panxα([1-6])", "testing\\1")
    assert seqs_to_hash(tester) == "f0d7ec055b3b2d9ec11a86634b32e9ef"

    tester = Sb._make_copy(sb_objects[0])
    Sb.rename(tester, "Panxα([1-6])", "testing\\1", -1)
    assert seqs_to_hash(tester) == "f0d7ec055b3b2d9ec11a86634b32e9ef"

    tester = Sb._make_copy(sb_objects[0])
    Sb.rename(tester, "[a-z]", "?", 2)
    assert seqs_to_hash(tester) == "08aaee2e7b9997b512c7d1b2fe748d40"

    tester = Sb._make_copy(sb_objects[0])
    Sb.rename(tester, "[a-z]", "?", -2)
    assert seqs_to_hash(tester) == "9b3946afde20c991099463d099be22e0"

    tester = Sb._make_copy(sb_objects[0])
    Sb.rename(tester, "[A-Z]", "?", -20)
    assert seqs_to_hash(tester) == "451993d7e816881e2700697263b1d8fa"

    with pytest.raises(AttributeError) as e:
        Sb.rename(tester, "[a-z]", "\\1?", -2)
    assert "There are more replacement match" in str(e)


# ##################### '-rs', 'replace_subseq' ###################### ##
def test_replace_subsequence():
    tester = Sb._make_copy(sb_objects[0])
    Sb.replace_subsequence(tester, "atg(.{5}).{3}", "FOO\\1BAR")
    assert seqs_to_hash(tester) == "f12707c2b0ef866f0039bac96abb29e0"


# ######################  '-rc', '--reverse_complement' ###################### #
hashes = ["e77be24b8a7067ed54f06e0db893ce27", "47941614adfcc5bd107f71abef8b3e00", "f549c8dc076f6b3b4cf5a1bc47bf269d",
          "0dd20827c011a0a7a0e78881b38ae06a", "0b954f4a263bf38ddeac61ab54f77dc2", "0d6b7deda824b4fc42b65cb87e1d4d14"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_reverse_complement(seqbuddy, next_hash):
    tester = Sb.reverse_complement(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_reverse_complement_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError) as e:
        Sb.reverse_complement(sb_objects[6])
    assert str(e.value) == "SeqBuddy object is protein. Nucleic acid sequences required."

# ######################  '-sfr', '--select_frame' ###################### #
hashes = [(0, 1, "b831e901d8b6b1ba52bad797bad92d14"), (0, 2, "2de033b2bf2327f2795fe425db0bd78f"),
          (0, 3, "1c29898d4964e0d1b03207d7e67e1958"), (1, 1, "908744b00d9f3392a64b4b18f0db9fee"),
          (1, 2, "08fe54a87249f5fb9ba22ff6d0053787"), (1, 3, "cfe2d405487d69dceb2a11dd44ceec59"),
          (2, 1, "cb1169c2dd357771a97a02ae2160935d"), (2, 2, "87d784f197b55f812d2fc82774da43d1"),
          (2, 3, "5d6d2f337ecdc6f9a85e981c975f3e08")]
hashes = [(Sb._make_copy(sb_objects[sb_indx]), _hash, frame) for sb_indx, frame, _hash in hashes]


@pytest.mark.parametrize("seqbuddy,next_hash,shift", hashes)
def test_select_frame(seqbuddy, next_hash, shift):
    tester = Sb.select_frame(seqbuddy, shift)
    assert seqs_to_hash(tester) == next_hash


def test_select_frame_edges():
    tester = Sb.select_frame(Sb._make_copy(sb_objects[0]), 2)
    temp_file = MyFuncs.TempFile()
    tester.write(temp_file.path)
    tester = Sb.select_frame(Sb.SeqBuddy(temp_file.path), 1)
    assert seqs_to_hash(tester) == "b831e901d8b6b1ba52bad797bad92d14"

    tester = Sb.select_frame(Sb._make_copy(sb_objects[1]), 2)
    tester = Sb.select_frame(tester, 1)
    assert seqs_to_hash(tester) == "908744b00d9f3392a64b4b18f0db9fee"

    with pytest.raises(TypeError) as e:  # If protein is input
        Sb.select_frame(sb_objects[6], 2)
    assert "Select frame requires nucleic acid, not protein." in str(e.value)


# ##################### '-ss', 'shuffle_seqs' ###################### ##
def test_shuffle_seqs():
    for seqbuddy in sb_objects:
        tester1 = Sb.shuffle_seqs(Sb._make_copy(seqbuddy))
        tester2 = Sb._make_copy(seqbuddy)
        assert seqs_to_hash(tester1) != seqs_to_hash(tester2)

        for indx, record in enumerate(tester1.records):
            assert sorted(record.seq) == sorted(tester2.records[indx].seq)


# #####################  make_groups' ###################### ##
def test_make_groups():
    tester = Sb.SeqBuddy(resource("Cnidaria_pep.nexus"))
    sb_list = Sb.make_groups(tester)
    assert len(sb_list) == 20
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier == seqbuddy.records[0].id

    sb_list = Sb.make_groups(tester, split_patterns=["u", "h"])
    assert len(sb_list) == 4
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Unknown", "Hv", "C", "Pp"]

    sb_list = Sb.make_groups(tester, num_chars=1)
    assert len(sb_list) == 5
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["A", "H", "C", "P", "N"]

    sb_list = Sb.make_groups(tester, split_patterns=["l"], num_chars=3)
    assert len(sb_list) == 3
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Unknown", "Ae", "C"]

    sb_list = Sb.make_groups(tester, regex="([ACH]).*([βγ])")
    assert len(sb_list) == 5
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Unknown", "Aβ", "Cβ", "Hβ", "Cγ"]

    sb_list = Sb.make_groups(tester, regex="Ate")
    assert len(sb_list) == 2
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Ate", "Unknown"]

    sb_list = Sb.make_groups(tester, regex="Panx.(G*)")
    assert len(sb_list) == 2
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["G", "Unknown"]


# ######################  '-tr6', '--translate6frames' ###################### #
# Only fasta and genbank
hashes = ["95cf24202007399e6ccd6e6f33ae012e", "0b5daa810e1589c3973e1436c40baf08"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_translate6frames(seqbuddy, next_hash):
    tester = Sb.translate6frames(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_translate6frames_pep_exception():
    with pytest.raises(TypeError):
        Sb.translate6frames(Sb._make_copy(sb_objects[6]))


# ######################  '-tr', '--translate' ###################### #
hashes = ["06893e14839dc0448e6f522c1b8f8957", "e8840e22096e933ce10dbd91036f3fa5", "f3339e0193c10427f017dd8f6bd81d7e"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]
hashes.append((Sb._make_copy(sb_objects[13]), "648ccc7c3400882be5bf6e8d9781f74e"))


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_translate(seqbuddy, next_hash):
    tester = Sb.translate_cds(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_translate_edges_and_exceptions(capsys):
    with pytest.raises(TypeError):
        Sb.translate_cds(Sb._make_copy(sb_objects[6]))

    tester = Sb.SeqBuddy("ATTCGTTAACGCTAGCGTCG", in_format="raw")
    tester = Sb.translate_cds(tester)
    assert str(tester) == ">raw_input\nIR*R*R\n"
    out, err = capsys.readouterr()
    assert err == "Warning: size mismatch between aa and nucl seqs for raw_input --> 20, 6\n"

    tester = Sb.select_frame(Sb._make_copy(sb_objects[1]), 3)
    tester = Sb.translate_cds(tester)
    assert seqs_to_hash(tester) == "68ca15f5ac737e4a4ca65a67ad2dc897"
    out, err = capsys.readouterr()
    assert string2hash(err) == "9e2a0b4b03f54c209d3a9111792762df"


# ################################################# COMMAND LINE UI ################################################## #
# ##################### '-ano', '--annotate' ###################### ##
def test_annotate_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.annotate = ["misc_feature", "1-100,200-250", "+"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "4a294c1ea5783fe04e909d767e7dfce2"

    test_in_args = deepcopy(in_args)
    test_in_args.annotate = ["misc_feature", "1-100,200-250", "foo=bar", "hello=world", "-", "α4"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "42ca93bcf3103124287fffbde49f2b26"

    test_in_args = deepcopy(in_args)
    test_in_args.annotate = ["unknown_feature_that_is_t0o_long", "1-100,200-250", "foo=bar", "hello=world", "-", "α4"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Warning: The provided annotation type is not part of the GenBank format standard" in err
    assert "Warning: Feature type is longer than 16 characters" in err


# ######################  '-asl', '--ave_seq_length' ###################### #
def test_ave_seq_length_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.ave_seq_length = [False]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert out == '428.38\n'

    test_in_args.ave_seq_length = ['clean']
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert out == '427.38\n'


# ######################  '-r2d', '--back_transcribe' ###################### #
def test_back_transcribe_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.back_transcribe = True
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("Mnemiopsis_rna.fa")), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b831e901d8b6b1ba52bad797bad92d14"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]))

    out, err = capsys.readouterr()
    assert err == "TypeError: RNA sequence required, not IUPACAmbiguousDNA().\n"


# ######################  '-btr', '--back_translate' ###################### #
def test_back_translate_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.back_translate = [False]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert "Panxα4" in out

    test_in_args = deepcopy(in_args)
    test_in_args.back_translate = [["human", "o"]]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[7]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b6bcb4e5104cb202db0ec4c9fc2eaed2"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]))

    out, err = capsys.readouterr()
    assert err == "TypeError: The input sequence needs to be protein, not nucleotide\n"


# ######################  '-bl2s', '--bl2seq' ###################### #
def test_bl2s_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.bl2seq = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "339377aee781fb9d01456f04553e3923"

    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("Duplicate_seqs.fa")), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "d24495bd87371cd0720084b5d723a4fc"
    assert err == "Warning: There are records with duplicate ids which will be renamed.\n"

    os.environ["PATH"] = ""
    with mock.patch('builtins.input', return_value="n"):
        with pytest.raises(SystemExit):
            Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]))
    out, err = capsys.readouterr()
    assert "not present in $PATH or working directory" in err
    os.environ["PATH"] = BACKUP_PATH


# ######################  '-bl', '--blast' ###################### #
def test_blast_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.blast = resource("blast/Mnemiopsis_cds")
    tester = Sb._make_copy(sb_objects[0])
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "a56ec76a64b25b7ca8587c7aa8554412"

    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("blank.fa")), True)
    out, err = capsys.readouterr()
    assert out == "No significant matches found\n"

    with pytest.raises(SystemExit):
        test_in_args.blast = resource("./Mnemiopsis_cds")
        Sb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "RuntimeError:" in err


# ######################  '-cs', '--clean_seq' ###################### #
def test_clean_seq_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.clean_seq = [[None]]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "dc53f3be7a7c24425dddeea26ea0ebb5"

    test_in_args.clean_seq = [["strict"]]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("ambiguous_dna.fa")), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "1912fadb5ec52a38ec707c58085b86ad"

    test_in_args.clean_seq = [["strict", "X"]]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("ambiguous_dna.fa")), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "4c10ba4474d7484652cb633f03db1be1"


# ######################  '-cmp', '--complement' ###################### #
def test_complement_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.complement = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "e4a358ca57aca0bbd220dc6c04c88795"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]))

    out, err = capsys.readouterr()
    assert "Nucleic acid sequence required, not protein." in err


# ######################  'cts', '--concat_seqs' ###################### #
def test_concat_seqs_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.concat_seqs = [True]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "7421c27be7b41aeedea73ff41869ac47"

    test_in_args.concat_seqs = ["clean"]
    test_in_args.out_format = "embl"
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[2]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "15b9c79dea034cef74e3a622bd357705"


# ######################  '-cc', '--count_codons' ###################### #
def test_count_codons_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.count_codons = ["foo"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "5661f0ce92bb6cfba4519a61e0a838ed"

    test_in_args.count_codons = ["conc"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "3e76bd510de4a61efb17ffc186ef9e68"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[7]))
    out, err = capsys.readouterr()
    assert "Nucleic acid sequence required, not protein" in err


# ######################  '-cr', '--count_residues' ###################### #
def test_count_residues_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.count_residues = ["foo"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "c9bc54835fb54232d4d346d04344bf8b"

    test_in_args.count_residues = ["conc"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "c7062408f939f4b310f2f97c3e94eb37"


# ######################  '-dgn' '--degenerate_sequence'################### #
def test_degenerate_sequence_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.degenerate_sequence = [False]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == '0638bc6546eebd9d50f771367d6d7855'        

    test_in_args.degenerate_sequence = [2, 7]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == '72373f8356051e2c6b67642451379054'

    test_in_args.degenerate_sequence = [100]
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]))
    out, err = capsys.readouterr()
    assert "Could not locate codon dictionary" in err

    test_in_args.degenerate_sequence = [1]
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[7]))
    out, err = capsys.readouterr()
    assert "Nucleic acid sequence required, not protein" in err


# ######################  '-df', '--delete_features' ###################### #
def test_delete_features_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_features = ["donor"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "f84df6a77063c7def13babfaa0555bbf"


# ######################  '-dlg', '--delete_large' ###################### #
def test_delete_large_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_large = 1285
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "25859dc69d46651a1e04a70c07741b35"


# ######################  'dm', '--delete_metadata' ###################### #
def test_delete_metadata_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_metadata = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "544ab887248a398d6dd1aab513bae5b1"


# ######################  '-dr', '--delete_records' ###################### #
def test_delete_records_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_records = ['α1']
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "54e810265a6ecf7a3d140fc806597f93"
    assert string2hash(err) == "4f420c9128e515dc24031b5075c034e3"

    test_in_args.delete_records = ['α1', 'α2']
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "eca4f181dae3d7998464ff71e277128f"
    assert string2hash(err) == "b3983fac3c2cf15f83650a34a17151da"

    test_in_args.delete_records = ['α1', 'α2', "3"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "eca4f181dae3d7998464ff71e277128f"
    assert string2hash(err) == "7e0929af515502484feb4b1b2c35eaba"

    test_in_args.delete_records = ['foo']
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b831e901d8b6b1ba52bad797bad92d14"
    assert string2hash(err) == "553348fa37d9c67f4ce0c8c53b578481"


# ######################  '-drp', '--delete_repeats' ###################### #
def test_delete_repeats_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_repeats = [None]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("Duplicate_seqs.fa")), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "df8e5f139ff41e1d81a082b83c208e12"
    assert string2hash(err) == "3c27f0df0e892a1c66ed8fef047162ae"

    test_in_args.delete_repeats = [[2, "all"]]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b831e901d8b6b1ba52bad797bad92d14"
    assert err == "No duplicate records found\n"

    test_in_args.quiet = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert not err


# ######################  '-ds', '--delete_small' ###################### #
def test_delete_small_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_small = 1285
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "196adf08d4993c51050289e5167dacdf"


# ######################  '-er', '--extract_region' ###################### #
def test_extract_region_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.extract_region = [50, 300]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "3e791c6a6683516aff9572c24f38f0b3"

    with pytest.raises(SystemExit):
        test_in_args.extract_region = [300, 50]
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]))
    out, err = capsys.readouterr()
    assert "The value given for end of range is smaller than for" in err


# ######################  '-fcpg', '--find_cpg' ###################### #
def test_find_cpg_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.find_CpG = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "9499f524da0c35a60502031e94864928"
    assert string2hash(err) == "1fc07e5884a2a4c1344865f385b1dc79"

    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(">seq1\nATGCCTAGCTAGCT", in_format="fasta"), True)
    out, err = capsys.readouterr()
    assert err == "# No Islands identified\n\n"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[7]))
    out, err = capsys.readouterr()
    assert "DNA sequence required, not protein or RNA" in err


# ######################  '-fp', '--find_pattern' ###################### #
def test_find_pattern_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.find_pattern = ["ATg{2}T", "tga.{1,6}tg"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()

    assert string2hash(out) == "0e11b2c0e9451fbfcbe39e3b5be2cf60"
    assert string2hash(err) == "f2bb95f89e7b9e198f18a049afbe4a93"


# ######################  '-frp', '--find_repeats' ###################### #
def test_find_repeats_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.find_repeats = [True]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert "#### No records with duplicate IDs ####" in out and "#### No records with duplicate sequences ####" in out

    tester = Sb.SeqBuddy(resource("Duplicate_seqs.fa"))
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "58a57c8151c3591fbac2b94353038a55"

    test_in_args.find_repeats = [2]
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b34b99828596a5a46c6ab244c6ccc6f6"


# ######################  '-frs', '--find_restriction_sites' ###################### #
def test_find_restriction_sites_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.find_restriction_sites = [["MaeI", "BseRI", "BccI", "MboII", 3, 4, 2, 5, "alpha"]]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b06ef2b0a4814fc43a0688f05825486a"
    assert string2hash(err) == "a240a6db9dfc1f2257faa80bc4b1445b"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[7]))
    out, err = capsys.readouterr()
    assert "Unable to identify restriction sites in protein sequences." in err


# ######################  '-gbp', '--group_by_prefix' ###################### #
def test_group_by_prefix_ui(capsys):
    tester = Sb.SeqBuddy(resource("Cnidaria_pep.nexus"))
    test_in_args = deepcopy(in_args)
    test_in_args.group_by_prefix = [[TEMP_DIR.path]]
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert out == ""
    assert "New file: " in err
    assert "Ate.nex" in err
    for prefix in ["Ate", "Hvu", "Che", "Ael", "Cla", "Hec", "Pph", "Nbi", "Ccr"]:
        assert os.path.isfile("%s/%s.nex" % (TEMP_DIR.path, prefix))
        os.unlink("%s/%s.nex" % (TEMP_DIR.path, prefix))

    test_in_args.group_by_prefix = [[TEMP_DIR.path, "u", "h"]]
    Sb.command_line_ui(test_in_args, tester, True)
    for prefix in ["Unknown", "Hv", "C", "Pp"]:
        assert os.path.isfile("%s/%s.nex" % (TEMP_DIR.path, prefix))
        os.unlink("%s/%s.nex" % (TEMP_DIR.path, prefix))

    test_in_args.group_by_prefix = [[TEMP_DIR.path, 1]]
    Sb.command_line_ui(test_in_args, tester, True)
    for prefix in ["A", "H", "C", "P", "N"]:
        assert os.path.isfile("%s/%s.nex" % (TEMP_DIR.path, prefix))
        os.unlink("%s/%s.nex" % (TEMP_DIR.path, prefix))

    test_in_args.group_by_prefix = [[TEMP_DIR.path, "l", 3]]
    Sb.command_line_ui(test_in_args, tester, True)
    for prefix in ["Unknown", "Ae", "C"]:
        assert os.path.isfile("%s/%s.nex" % (TEMP_DIR.path, prefix))
        os.unlink("%s/%s.nex" % (TEMP_DIR.path, prefix))


# ######################  '-gbr', '--group_by_regex' ###################### #
def test_group_by_regex_ui(capsys):
    tester = Sb.SeqBuddy(resource("Cnidaria_pep.nexus"))
    test_in_args = deepcopy(in_args)
    test_in_args.group_by_regex = [[TEMP_DIR.path]]
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert err == "ValueError: You must provide at least one regular expression.\n"

    test_in_args.group_by_regex = [[TEMP_DIR.path, "Ate"]]
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert "New file: " in err
    assert "Ate.nex" in err
    for prefix in ["Unknown", "Ate"]:
        assert os.path.isfile("%s/%s.nex" % (TEMP_DIR.path, prefix))
        os.unlink("%s/%s.nex" % (TEMP_DIR.path, prefix))


# ######################  '-gf', '--guess_format' ###################### #
def test_guess_alpha_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.guess_alphabet = True
    paths = ["Mnemiopsis_%s.fa" % x for x in ["cds", "pep", "rna"]]
    paths += ["gibberish.fa", "figtree.nexus"]
    test_in_args.sequence = [resource(x) for x in paths]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "3a2edec76860e60f4f5b8b16b6d32b82"

    text_io = io.open(resource("Mnemiopsis_cds.embl"), "r")
    test_in_args.sequence = [text_io]
    tester = Sb.SeqBuddy(text_io)
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert out == "PIPE\t-->\tdna\n"

    temp_file = MyFuncs.TempFile()
    temp_file.write(">seq1\n123456789")
    text_io = io.open(temp_file.path, "r")
    test_in_args.sequence = [text_io]
    tester = Sb.SeqBuddy(text_io)
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert out == "PIPE\t-->\tUndetermined\n"


# ######################  '-gf', '--guess_format' ###################### #
def test_guess_format_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.guess_format = True
    paths = ["Mnemiopsis_cds.%s" % x for x in ["embl", "fa", "gb", "nex", "phy", "phyr", "physs", "physr",
                                               "seqxml", "stklm", "clus"]]
    paths += ["gibberish.fa", "figtree.nexus"]
    test_in_args.sequence = [resource(x) for x in paths]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "c0b299566f102a6644b49744bc75ddcf"

    text_io = io.open(resource("Mnemiopsis_cds.embl"), "r")
    test_in_args.sequence = [text_io]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy, True)
    out, err = capsys.readouterr()
    assert out == "PIPE\t-->\tembl\n"


# ######################  '-hsi', '--hash_seq_ids' ###################### #
def test_hash_seq_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.hash_seq_ids = [None]
    tester = Sb._make_copy(sb_objects[0])
    ids = [rec.id for rec in tester.records]
    Sb.command_line_ui(test_in_args, tester, True)
    for indx, rec in enumerate(tester.records):
        assert rec.id != ids[indx]
        assert ids[indx] == tester.hash_map[rec.id]

    test_in_args.hash_seq_ids = [0]
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert "Warning: The hash_length parameter was passed in with the value 0. This is not a positive integer" in err

    tester.records *= 10
    test_in_args.hash_seq_ids = [1]
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert "cover all sequences, so it has been increased to 2" in err


# ######################  '-is', '--insert_seq' ###################### #
def test_insert_seqs_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.insert_seq = [["DYKDDDDK"]]
    tester = Sb._make_copy(sb_objects[6])
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "The insert_seq tool requires at least two arguments (sequence and position)" in err

    test_in_args.insert_seq = [[4, "DYKDDDDK"]]
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "The first argment must be your insert sequence, not location." in err

    test_in_args.insert_seq = [["DYKDDDDK", "Foo"]]
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "The second argment must be location, not insert sequence or regex." in err

    test_in_args.insert_seq = [["DYKDDDDK", "10", "α[23]", "α6"]]
    tester = Sb._make_copy(sb_objects[6])
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "345836c75922e5e2a7367c7f7748b591"


# ######################  '-ip', '--isoelectric_point' ###################### #
def test_isoelectric_point_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.isoelectric_point = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "2b0de14da980b6b2c155c34f41da814e"
    assert string2hash(err) == "402411565abc86649581bf7ab65535b8"

    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "d1ba12963ee508bc64b64f63464bfb4a"
    assert err == "ID\tpI\n"


# ######################  '-li', '--list_ids' ###################### #
def test_list_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.list_ids = [3]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "53d5d7afd8f15a1a0957f5d5a29cbdc4"


# ######################  '-lf', '--list_features' ###################### #
def test_list_features_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.list_features = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b99acb13c76f86bcd4e8dc15b97fa11d"

    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "4e37613d1916aa7653d3fec37fc9e368"


# ######################  '-lc', '--lowercase' and 'uc', '--uppercase'  ###################### #
def test_lower_and_upper_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.uppercase = True
    tester = Sb._make_copy(sb_objects[0])
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "25073539df4a982b7f99c72dd280bb8f"

    test_in_args.uppercase = False
    test_in_args.lowercase = True
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b831e901d8b6b1ba52bad797bad92d14"


# ######################  '-mui', '--make_ids_unique' ###################### #
def test_make_ids_unique_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.make_ids_unique = True
    tester = Sb.SeqBuddy(resource("Duplicate_seqs.fa"))
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "3f84819b81185d49a7357947ca7f0707"


# ######################  '-fn2p', '--map_features_nucl2prot' ###################### #
def test_map_features_nucl2prot_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.map_features_nucl2prot = True
    test_in_args.sequence = [resource("Mnemiopsis_cds.gb"), resource("Mnemiopsis_pep.fa")]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "5216ef85afec36d5282578458a41169a"

    test_in_args.sequence = [resource("Mnemiopsis_pep.fa"), resource("Mnemiopsis_cds.gb")]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "5216ef85afec36d5282578458a41169a"

    test_in_args.sequence = [resource("Mnemiopsis_pep.fa"), resource("Mnemiopsis_cds.gb")]
    test_in_args.out_format = "embl"
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "4f86356e79fa4beb79961ce37b5aa19a"

    with pytest.raises(SystemExit):
        test_in_args.sequence = [resource("Mnemiopsis_cds.gb"), resource("Duplicate_seqs.fa")]
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "There are repeat IDs in self.records" in err

    with pytest.raises(SystemExit):
        test_in_args.sequence = [resource("Mnemiopsis_cds.gb")]
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "You must provide one DNA file and one protein file" in err

    with pytest.raises(SystemExit):
        test_in_args.sequence = [resource("Mnemiopsis_cds.gb"), resource("Mnemiopsis_cds.fa")]
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "You must provide one DNA file and one protein file" in err


# ######################  '-fp2n', '--map_features_prot2nucl' ###################### #
def test_map_features_prot2nucl_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.map_features_prot2nucl = True
    test_in_args.sequence = [resource("Mnemiopsis_cds.fa"), resource("Mnemiopsis_pep.gb")]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "3ebc92ca11505489cab2453d2ebdfcf2"

    test_in_args.sequence = [resource("Mnemiopsis_pep.gb"), resource("Mnemiopsis_cds.fa")]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "3ebc92ca11505489cab2453d2ebdfcf2"

    test_in_args.sequence = [resource("Mnemiopsis_pep.gb"), resource("Mnemiopsis_cds.fa")]
    test_in_args.out_format = "embl"
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "bbbfc9ebc83d3abe3bb3160a38d208e3"

    with pytest.raises(SystemExit):
        temp_file = MyFuncs.TempFile()
        duplicate_seqs = Sb.SeqBuddy(resource("Duplicate_seqs.fa"))
        Sb.back_translate(duplicate_seqs)
        duplicate_seqs.write(temp_file.path)
        test_in_args.sequence = [resource("Mnemiopsis_pep.gb"), temp_file.path]
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "There are repeat IDs in self.records" in err

    with pytest.raises(SystemExit):
        test_in_args.sequence = [resource("Mnemiopsis_pep.gb")]
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "You must provide one DNA file and one protein file" in err

    with pytest.raises(SystemExit):
        test_in_args.sequence = [resource("Mnemiopsis_pep.gb"), resource("Mnemiopsis_pep.fa")]
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "You must provide one DNA file and one protein file" in err


# ######################  '-mg', '--merge' ###################### #
def test_merge_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.merge = True
    test_in_args.sequence = [resource("Mnemiopsis_cds_dummy_features.gb"), resource("Mnemiopsis_cds.gb")]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "5fb43a869bb9f766625d1336081766fc"

    with pytest.raises(SystemExit):
        test_in_args.sequence = [resource("Mnemiopsis_pep.gb"), resource("Mnemiopsis_cds.gb")]
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "RuntimeError" in err


# ######################  '-mw', '--molecular_weight' ###################### #
def test_molecular_weight_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.molecular_weight = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "7a456f37b9d7a780dfe81e453f2e9525"
    assert err == "ID\tssDNA\tdsDNA\n"

    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "d1014a98fe227f8847ed7478bbdfc857"
    assert err == "ID\tProtein\n"

    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[13]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "55ff25f26504c5360557c2dfeb041036"
    assert err == "ID\tssRNA\n"


# ######################  '-ns', '--num_seqs' ###################### #
def test_num_seqs_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.num_seqs = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert out == '13\n'


# ######################  '-ofa', '--order_features_alphabetically' ###################### #
def test_order_features_alphabetically_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.order_features_alphabetically = [True]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == '21547b4b35e49fa37e5c5b858808befb'

    test_in_args.order_features_alphabetically = ["rev"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == '3b718ec3cb794bcb658d900e517110cc'


# ######################  '-ofp', '--order_features_by_position' ###################### #
def test_order_features_by_position_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.order_features_by_position = [True]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == '2e02a8e079267bd9add3c39f759b252c'

    test_in_args.order_features_by_position = ["rev"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == '4345a14fe27570b3c837c30a8cb55ea9'


# ######################  '-oi', '--order_ids' ###################### #
def test_order_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.order_ids = [True]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == '2f9bc0dd9d79fd8160a621280be0b0aa'

    test_in_args.order_ids = ["rev"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == 'dd269961d4d5301d1bf87e0093568851'


# ######################  '-oir', '--order_ids_randomly' ###################### #
def test_order_ids_randomly_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.order_ids_randomly = [True]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) != seqs_to_hash(sb_objects[0])

    tester = Sb.order_ids(Sb.SeqBuddy(out))
    assert seqs_to_hash(tester) == seqs_to_hash(Sb.order_ids(Sb._make_copy(sb_objects[0])))


# ######################  '-prr', '--pull_random_recs' ###################### #
def test_pull_random_recs_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.pull_random_record = [True]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    tester = Sb.SeqBuddy(out)
    assert len(tester.records) == 1
    assert tester.records[0].id in sb_objects[0].to_dict()

    test_in_args.pull_random_record = [20]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    tester = Sb.SeqBuddy(out)
    assert len(tester.records) == 13
    assert sorted([rec.id for rec in tester.records]) == sorted([rec.id for rec in sb_objects[0].records])


# ######################  '-pr', '--pull_record_ends' ###################### #
def test_pull_record_ends_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.pull_record_ends = 10
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "754d6868030d1122b35386118612db72"

    test_in_args.pull_record_ends = -10
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "9cfc91c3fdc5cd9daabce0ef9bac2db7"


# ######################  '-pr', '--pull_records' ###################### #
def test_pull_records_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.pull_records = ["α1"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "db52337c628fd8d8d70a5581355c51a5"

    test_in_args.pull_records = ["α1", "α2"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "cd8d7284f039233e090c16e8aa6b5035"


# ######################  '-prg', '--purge' ###################### #
def test_purge_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.purge = 200
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b21b2e2f0ca1fcd7b25efbbe9c08858c"
    assert string2hash(err) == "fbfde496ae179f83e3d096da15d90920"


# ######################  '-ri', '--rename_ids' ###################### #
def test_rename_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.rename_ids = [["[a-z](.)", "?\\1", 2]]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "f12c44334b507117439928c529eb2944"

    test_in_args.rename_ids = [["[a-z](.)"]]
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "rename_ids requires two or three argments:" in err

    test_in_args.rename_ids = [["[a-z](.)", "?\\1", 2, "foo"]]
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "rename_ids requires two or three argments:" in err

    test_in_args.rename_ids = [["[a-z](.)", "?\\1", "foo"]]
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "Max replacements argument must be an integer" in err

    test_in_args.rename_ids = [["[a-z](.)", "?\\1\\2", 2]]
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]))
    out, err = capsys.readouterr()
    assert "There are more replacement" in err


# ######################  '-rs', '--replace_subseq' ###################### #
def test_replace_subseq_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.replace_subseq = [["atg(.{5}).{3}", "FOO\\1BAR"]]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "4e1b13745d256331ccb46dd275627edb"


# ######################  '-rc', '--reverse_complement' ###################### #
def test_reverse_complement_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.reverse_complement = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "47941614adfcc5bd107f71abef8b3e00"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[7]))
    out, err = capsys.readouterr()
    assert "SeqBuddy object is protein. Nucleic acid sequences required." in err

    tester = Sb.SeqBuddy(resource("mixed_alpha.fa"))
    tester.alpha = IUPAC.ambiguous_dna
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "efbcc71f3b2820ea05bf32038012b883"

    tester.records[0].seq.alphabet = IUPAC.protein
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert err == "TypeError: Record 'Mle-Panxα12' is protein. Nucleic acid sequences required.\n"

# ######################  '-sf', '--screw_formats' ###################### #
hashes = [("fasta", "09f92be10f39c7ce3f5671ef2534ac17"), ("gb", "26718f0a656116bfd0a7f6c03d270ecf"),
          ("nexus", "2822cc00c2183a0d01e3b79388d344b3"), ("phylip", "6a4d62e1ee130b324cce48323c6d1d41"),
          ("phylip-relaxed", "4c2c5900a57aad343cfdb8b35a8f8442"), ("phylipss", "089cfb52076e63570597a74b2b000660"),
          ("phylipsr", "58a74f5e08afa0335ccfed0bdd94d3f2"), ("stockholm", "8c0f5e2aea7334a0f2774b0366d6da0b"),
          ("raw", "f0ce73f4d05a5fb3d222fb0277ff61d2")]


@pytest.mark.parametrize("_format,next_hash", hashes)
def test_screw_formats_ui(_format, next_hash, capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = _format
    tester = Sb.pull_recs(sb_objects[2], "α[2-9]")
    Sb.command_line_ui(test_in_args, Sb._make_copy(tester), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == next_hash


def test_screw_formats_ui2(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = "foo"
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "Error: unknown format" in err

    sb_objects[0].write("%s/seq.fa" % TEMP_DIR.path)
    test_in_args.sequence = ["%s/seq.fa" % TEMP_DIR.path]
    test_in_args.screw_formats = "genbank"
    test_in_args.in_place = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    assert os.path.isfile("%s/seq.gb" % TEMP_DIR.path)


# ######################  '-sfr', '--select_frame' ###################### #
def test_select_frame_ui(capsys):
    test_in_args = deepcopy(in_args)
    tester = Sb._make_copy(sb_objects[1])
    test_in_args.select_frame = 1
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "908744b00d9f3392a64b4b18f0db9fee"

    test_in_args.select_frame = 2
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "08fe54a87249f5fb9ba22ff6d0053787"

    test_in_args.select_frame = 3
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "cfe2d405487d69dceb2a11dd44ceec59"

    test_in_args.select_frame = 1
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "908744b00d9f3392a64b4b18f0db9fee"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]))

    out, err = capsys.readouterr()
    assert "Select frame requires nucleic acid, not protein" in err


# ######################  '-ss', '--shuffle_seqs' ###################### #
def test_shuffle_seqs_ui(capsys):
    test_in_args = deepcopy(in_args)
    tester = Sb._make_copy(sb_objects[0])
    test_in_args.shuffle_seqs = True
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) != "b831e901d8b6b1ba52bad797bad92d14"


# ######################  '-d2r', '--transcribe' ###################### #
def test_transcribe_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.transcribe = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "d2db9b02485e80323c487c1dd6f1425b"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]))

    out, err = capsys.readouterr()
    assert err == "TypeError: DNA sequence required, not IUPACProtein().\n"


# ######################  '-tr', '--translate' ###################### #
def test_translate_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.translate = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "06893e14839dc0448e6f522c1b8f8957"

    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("ambiguous_rna.fa")), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "648ccc7c3400882be5bf6e8d9781f74e"

    tester = Sb.SeqBuddy(resource("mixed_alpha.fa"))
    tester.alpha = IUPAC.ambiguous_dna
    tester.records[0].seq.alphabet = IUPAC.protein
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert err == 'TypeError: Record Mle-Panxα12 is protein.\n'

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]))
    out, err = capsys.readouterr()
    assert "Nucleic acid sequence required, not protein." in err


# ######################  '-tr6', '--translate6frames' ###################### #
def test_translate6frames_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.translate6frames = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "95cf24202007399e6ccd6e6f33ae012e"

    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("ambiguous_rna.fa")), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "6d3fcad0ea417014bc825cedd354fd26"
    assert err == ""

    tester = Sb.SeqBuddy(resource("mixed_alpha.fa"))
    tester.alpha = IUPAC.ambiguous_dna
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b54ec5c49bd88126de337e1eb3d2ad23"

    tester.records[0].seq.alphabet = IUPAC.protein
    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert err == "TypeError: Record 'Mle-Panxα12' is protein. Nucleic acid sequences required.\n"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]))
    out, err = capsys.readouterr()
    assert "TypeError: You need to supply DNA or RNA sequences to translate" in err
