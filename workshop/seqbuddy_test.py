#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest
import warnings

import SeqBuddy


class SeqbuddyTest(unittest.TestCase):
    fasta_cds = SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_cds.fa")
    fasta_pep = SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_pep.fa")
    gb_cds = SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_cds.gb")
    gb_pep = SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_pep.gb")

    # Ensure that SeqBuddy objects are correctly loaded from file path, file handle, or raw sequence
    def test_instantiate_seqbuddy(self):
        warnings.simplefilter("ignore", ResourceWarning)
        SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_cds.fa")
        SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_pep.fa")
        SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_cds.gb")
        SeqBuddy.SeqBuddy("unit_test_resources/Mnemiopsis_pep.gb")
        # Need to add lots more here...


if __name__ == '__main__':
    unittest.main()