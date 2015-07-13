#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

import pytest
from hashlib import md5
import os
from io import StringIO
import re
from copy import deepcopy
from Bio.Alphabet import IUPAC
import SeqBuddy as SB

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


def test_guess_format():
    assert Alb.guess_format(["dummy", "list"]) == "stockholm"
    assert Alb.guess_format(alb_objects[0]) == "nexus"

    with open(resource("Alignments_pep.stklm"), "r") as ifile:
        assert Alb.guess_format(ifile) == "stockholm"
        ifile.seek(0)
        string_io = StringIO(ifile.read())
    assert Alb.guess_format(string_io) == "stockholm"
    with pytest.raises(Alb.GuessError):
        Alb.guess_format(("Dummy-Tuple"))


align_files = ["Mnemiopsis_cds.nex", "Mnemiopsis_cds.phy", "Mnemiopsis_cds.phyr", "Mnemiopsis_cds.stklm",
               "Mnemiopsis_pep.nex", "Mnemiopsis_pep.phy", "Mnemiopsis_pep.phyr", "Mnemiopsis_pep.stklm",
               "Alignments_pep.phy", "Alignments_pep.phyr", "Alignments_pep.stklm"]

file_types = ["nexus", "phylip", "phylip-relaxed", "stockholm",
              "nexus", "phylip", "phylip-relaxed", "stockholm",
              "phylip", "phylip-relaxed", "stockholm"]

input_tuples = [(next_file, file_types[indx]) for indx, next_file in enumerate(align_files)]


@pytest.mark.parametrize("align_file,file_type", input_tuples)
def test_instantiate_alignbuddy_from_file(align_file, file_type):
    assert type(Alb.AlignBuddy(resource(align_file), _in_format=file_type)) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_file_guess(align_file):
    assert type(Alb.AlignBuddy(resource(align_file))) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_handle(align_file):
    with open(resource(align_file), 'r') as ifile:
        assert type(Alb.AlignBuddy(ifile)) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_raw(align_file):
    with open(resource(align_file), 'r') as ifile:
        assert type(Alb.AlignBuddy(ifile.read())) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_alignbuddy(align_file):
    tester = Alb.AlignBuddy(resource(align_file))
    assert type(Alb.AlignBuddy(tester)) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_list(align_file):
    tester = Alb.AlignBuddy(resource(align_file))
    assert type(Alb.AlignBuddy(tester.alignments)) == Alb.AlignBuddy

    with pytest.raises(TypeError):  # When non-MultipleSeqAlignment objects are in the .alignments list
        tester.alignments.append("Dummy string object")
        Alb.AlignBuddy(tester.alignments)


def test_empty_file(capsys):
    with open(resource("blank.fa"), "r") as ifile:
        with pytest.raises(SystemExit):
            Alb.AlignBuddy(ifile)


def test_guess_error():
    # File path
    with pytest.raises(Alb.GuessError):
        Alb.AlignBuddy(resource("unrecognizable.txt"))

    with open(resource("unrecognizable.txt"), 'r') as ifile:
        # Raw
        with pytest.raises(Alb.GuessError):
            Alb.AlignBuddy(ifile.read())

        # Handle
        with pytest.raises(Alb.GuessError):
            ifile.seek(0)
            Alb.AlignBuddy(ifile)

    # GuessError output
    try:
        Alb.AlignBuddy(resource("unrecognizable.txt"))
    except Alb.GuessError as e:
        assert str(e) == "Could not determine format from _input file '/Users/bondsr/Documents/BuddySuite/workshop/" \
                         "unit_test_resources/unrecognizable.txt'.\nTry explicitly setting with -f flag."


def test_stderr(capsys):
    Alb._stderr("Hello std_err", quiet=False)
    out, err = capsys.readouterr()
    assert err == "Hello std_err"

    Alb._stderr("Hello std_err", quiet=True)
    out, err = capsys.readouterr()
    assert err == ""


def test_stdout(capsys):
    Alb._stdout("Hello std_out", quiet=False)
    out, err = capsys.readouterr()
    assert out == "Hello std_out"

    Alb._stdout("Hello std_out", quiet=True)
    out, err = capsys.readouterr()
    assert out == ""

# Now that we know that all the files are being turned into AlignBuddy objects okay, make them all objects so it doesn't
# need to be done over and over for each subsequent test.
alb_objects = [Alb.AlignBuddy(resource(x)) for x in align_files]


# AlignBuddy print() and __str__() methods
hashes = ["cb1169c2dd357771a97a02ae2160935d", "f59e28493949f78637691caeb617ab50",
          "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
          "17ff1b919cac899c5f918ce8d71904f6", "5af1cf061f003d3351417458c0d23811",
          "f3e98897f1bbb3d3d6c5507afa9f814e", "c0dce60745515b31a27de1f919083fe9",
          "90578980479ad235338dbb767444b05b", "9c6773e7d24000f8b72dd9d25620cff1",
          "3fd5805f61777f7f329767c5f0fb7467"]
hashes = [(alb_objects[indx], value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_print(alignbuddy, next_hash, capsys):
    alignbuddy.print()
    out, err = capsys.readouterr()
    out = "{0}\n".format(out.rstrip())
    tester = md5(out.encode()).hexdigest()
    assert tester == next_hash


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_str(alignbuddy, next_hash):
    tester = str(alignbuddy)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_write(alignbuddy, next_hash):
    alignbuddy.write("/tmp/alignbuddywritetest")
    with open("/tmp/alignbuddywritetest", "r") as ifile:
        out = "{0}\n".format(ifile.read().rstrip())
    tester = md5(out.encode()).hexdigest()
    assert tester == next_hash

def test_get_seq_recs():
    tester = str(Alb._get_seq_recs(alb_objects[8]))
    tester = md5(tester.encode()).hexdigest()
    assert tester == "6168f8b57d0ff78d70fd22ee09d713b5"


def test_phylipi():
    tester = Alb.phylipi(alb_objects[0], _format="relaxed")
    tester = "{0}\n".format(tester.rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "c5fb6a5ce437afa1a4004e4f8780ad68"

    tester = Alb.phylipi(alb_objects[8], _format="relaxed").rstrip()
    tester = "{0}\n".format(tester.rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "af97ddb03817ff050d3dfb42472c91e0"

    tester = Alb.phylipi(alb_objects[0], _format="strict").rstrip()
    tester = "{0}\n".format(tester.rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "270f1bac51b2e29c0e163d261795c5fe"

    tester = Alb.phylipi(alb_objects[8], _format="strict").rstrip()
    tester = "{0}\n".format(tester.rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "af97ddb03817ff050d3dfb42472c91e0"


def test_guess_alphabet():
    assert str(type(Alb.guess_alphabet(alb_objects[0]))) == "<class 'Bio.Alphabet.IUPAC.IUPACAmbiguousDNA'>"
    assert str(type(Alb.guess_alphabet(alb_objects[4]))) == "<class 'Bio.Alphabet.IUPAC.IUPACProtein'>"
    tester = Alb.AlignBuddy(resource("Mnemiopsis_rna.nex"))
    assert str(type(Alb.guess_alphabet(tester))) == "<class 'Bio.Alphabet.IUPAC.IUPACAmbiguousRNA'>"
    assert not Alb.guess_alphabet(Alb.AlignBuddy("", _in_format="fasta"))


def test_list_ids():
    tester = Alb.list_ids(alb_objects[0])
    tester = md5(tester.encode()).hexdigest()
    assert tester == "1c4a395d8aa3496d990c611c3b6c4d0a"

    tester = Alb.list_ids(alb_objects[0], _columns=4)
    tester = md5(tester.encode()).hexdigest()
    assert tester == "26fa56b979d009015612647c87a47a51"

    tester = Alb.list_ids(alb_objects[8])
    tester = md5(tester.encode()).hexdigest()
    assert tester == "7f7cc5d09164cb2f5deb915193b06639"

    tester = Alb.list_ids(alb_objects[8], _columns=4)
    tester = md5(tester.encode()).hexdigest()
    assert tester == "74bd0e70fd325d59c0399c4f8a0ea7c9"
