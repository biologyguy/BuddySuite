#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: Unit BioPython 16.6 + required.
import unittest
import warnings
from hashlib import md5

import SeqBuddy


def seqs_to_string(_seqbuddy, _format):
    _output = ""
    for rec in _seqbuddy.records:
        _output += rec.format(_format)
    return _output

fasta_cds = SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_cds.fa")
fasta_pep = SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_pep.fa")
gb_cds = SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_cds.gb")
gb_pep = SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_pep.gb")


def seqs_to_hash(_seqbuddy, _format):
    seqs_string = seqs_to_string(_seqbuddy, _format)
    seqs_string = seqs_string.encode()
    _hash = md5(seqs_string).hexdigest()
    return _hash


class SeqbuddyTest(unittest.TestCase):
    # Ensure that SeqBuddy objects are correctly loaded from file path, file handle, or raw sequence
    def test_instantiate_seqbuddy(self):
        warnings.simplefilter("ignore", ResourceWarning)
        SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_cds.fa")
        SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_pep.fa")
        SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_cds.gb")
        SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_pep.gb")
        # Need to add lots more here...

    def test_order_features_alphabetically(self):
        tester = SeqBuddy.order_features_alphabetically(gb_cds)
        self.assertEqual(seqs_to_hash(tester, "gb"), "0ca1c910e2d566556e839c3ca29e162d")

        tester = SeqBuddy.order_features_alphabetically(gb_cds, True)
        self.assertEqual(seqs_to_hash(tester, "gb"), "71bc8a8d19e5c7a6c54db4196f632cea")

        tester = SeqBuddy.order_features_alphabetically(gb_pep)
        self.assertEqual(seqs_to_hash(tester, "gb"), "2cd6ab6a922a6695ba2ffdbfa81f740d")

        tester = SeqBuddy.order_features_alphabetically(gb_pep, True)
        self.assertEqual(seqs_to_hash(tester, "gb"), "172f886830e8f8880c31417ec1821490")

        tester = SeqBuddy.order_features_alphabetically(fasta_cds)
        self.assertEqual(seqs_to_hash(tester, "fasta"), "25073539df4a982b7f99c72dd280bb8f")

if __name__ == '__main__':


    unittest.main()