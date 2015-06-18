#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

import pytest
from hashlib import md5
import os
import workshop.SeqBuddy as Sb

def seqs_to_string(_seqbuddy, _format):
    _output = ""
    for rec in _seqbuddy.records:
        _output += rec.format(_format)
    return _output


def seqs_to_hash(_seqbuddy, _format):
    seqs_string = seqs_to_string(_seqbuddy, _format)
    seqs_string = seqs_string.encode()
    _hash = md5(seqs_string).hexdigest()
    return _hash

root_dir = os.getcwd()

def resource(file_name):
    return "{0}/unit_test_resources/{1}".format(root_dir, file_name)

seq_files = ["Mnemiopsis_cds.fa", "Mnemiopsis_cds.gb", "Mnemiopsis_cds.nex", "Mnemiopsis_cds.phy",
             "Mnemiopsis_cds.phyr", "Mnemiopsis_cds.stklm",
             "Mnemiopsis_pep.fa", "Mnemiopsis_pep.gb", "Mnemiopsis_pep.nex", "Mnemiopsis_pep.phy",
             "Mnemiopsis_pep.phyr", "Mnemiopsis_pep.stklm"]


def test_instantiate_seqbuddy_from_file():
    for next_file in seq_files:
        assert Sb.SeqBuddy(resource(next_file))


def test_instantiate_seqbuddy_from_handle():
    for next_file in seq_files:
        with open(resource(next_file), 'r') as ifile:
            assert Sb.SeqBuddy(ifile)


def test_instantiate_seqbuddy_from_raw():
    for next_file in seq_files:
        with open(resource(next_file), 'r') as ifile:
            assert Sb.SeqBuddy(ifile.read())


"""
@pytest.fixture(scope="session")
def fasta_cds(request):
    return sb.SeqBuddy("unit_test_resources/Mnemiopsis_cds.fa")


@pytest.fixture(scope="session")
def fasta_pep(request):
    return sb.SeqBuddy("unit_test_resources/Mnemiopsis_pep.fa")


@pytest.fixture(scope="session")
def gb_cds(request):
    return sb.SeqBuddy("unit_test_resources/Mnemiopsis_cds.gb")


@pytest.fixture(scope="session")
def gb_pep(request):
    return sb.SeqBuddy("unit_test_resources/Mnemiopsis_pep.gb")


def test_order_features_alphabetically(fasta_cds, fasta_pep, gb_cds, gb_pep):
    tester = sb.order_features_alphabetically(gb_cds)
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