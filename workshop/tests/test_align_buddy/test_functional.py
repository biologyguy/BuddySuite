# coding=utf-8
""" tests basic functionality of AlignBuddy class """
import pytest
from AlignBuddy import AlignBuddy
from buddy_resources import (GuessError,
                             PhylipError)


class TestAlignBuddyInit:
    """ Instantiate AlignBuddy objects in various ways """

    def test_can_create_align_buddy(self, alignment_valid_resources):
        """ Create an AlignBuddy obj from all valid file types

            Valid file types for dna, rna, and pep include:
                'clustal'
                'fasta'
                'gb'
                'nexus'
                'phylip'
                'phylipr'
                'phylipss'
                'phylipsr'
                'stockholm'
        """

        for molecule in alignment_valid_resources:  # dna, rna, pep
            for quantity in alignment_valid_resources[molecule]:  # single, multi
                # rna resource only has a single alignment file so skip multi
                if molecule == 'rna' and quantity == 'multi':
                    continue

                for file_type in alignment_valid_resources[molecule][quantity]:  # clustal, fasta, etc
                    assert bool(AlignBuddy(alignment_valid_resources[molecule][quantity][file_type]))

    def test_throws_errors_on_invalid_files(self, alignment_bad_resources):
        """ expect AlignBuddy to raise errors on invalid filesr """
        with pytest.raises(GuessError):
            AlignBuddy(alignment_bad_resources['dna']['single']['fasta'])
