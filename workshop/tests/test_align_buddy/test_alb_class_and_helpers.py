#!/usr/bin/env python3
# coding=utf-8
""" tests basic functionality of AlignBuddy class """
import pytest
from AlignBuddy import AlignBuddy
from buddy_resources import GuessError, PhylipError


def test_instantiate_alignbuddy_from_file(alb_resources):
    for key, _path in alb_resources.get(mode="paths").items():
        in_format = alb_resources.parse_code(key, strict=True)
        in_format = alb_resources.single_letter_codes[in_format["format"][0]]
        assert type(AlignBuddy(_path, in_format=in_format)) == AlignBuddy


def test_instantiate_alignbuddy_from_file_guess(alb_resources):
    for _path in alb_resources.get_list(mode="paths"):
        assert type(AlignBuddy(_path)) == AlignBuddy


def test_instantiate_alignbuddy_from_handle(alb_resources):
    for _path in alb_resources.get_list(mode="paths"):
        with open(_path, 'r') as ifile:
            assert type(AlignBuddy(ifile)) == AlignBuddy


def test_instantiate_alignbuddy_from_raw(alb_resources):
    for _path in alb_resources.get_list(mode="paths"):
        with open(_path, 'r') as ifile:
            assert type(AlignBuddy(ifile.read())) == AlignBuddy


def test_instantiate_alignbuddy_from_alignbuddy(alb_resources):
    for alignbuddy in alb_resources.get_list():
        assert type(AlignBuddy(alignbuddy)) == AlignBuddy


def test_instantiate_alignbuddy_from_list(alb_resources):
    for alignbuddy in alb_resources.get_list():
        assert type(AlignBuddy(alignbuddy.alignments)) == AlignBuddy

    with pytest.raises(TypeError):  # When non-MultipleSeqAlignment objects are in the .alignments list
        alignbuddy = alb_resources.get_one("p s m")
        alignbuddy.alignments.append("Dummy string object")
        AlignBuddy(alignbuddy.alignments)


def test_instantiation_alignbuddy_errors(alignment_bad_resources):
    with pytest.raises(GuessError) as e:
        AlignBuddy(alignment_bad_resources["dna"]["single"]["fasta"])
    assert "Could not determine format from _input file" in str(e)

    tester = open(alignment_bad_resources["dna"]["single"]["fasta"], "r")
    with pytest.raises(GuessError) as e:
        AlignBuddy(tester.read())
    assert "Could not determine format from raw" in str(e)

    tester.seek(0)
    with pytest.raises(GuessError) as e:
        AlignBuddy(tester)
    assert "Could not determine format from input file-like object" in str(e)


def test_empty_file(alignment_bad_resources):
    with open(alignment_bad_resources["blank"], "r") as ifile:
        with pytest.raises(GuessError) as e:
            AlignBuddy(ifile)
        assert "Empty file" in str(e)


def test_throws_errors_on_invalid_files(alignment_bad_resources):
    """ expect AlignBuddy to raise errors on invalid filesr """
    with pytest.raises(GuessError):
        AlignBuddy(alignment_bad_resources['dna']['single']['fasta'])
