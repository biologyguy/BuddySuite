# coding=utf-8
""" tests basic functionality of AlignBuddy class """
import pytest
from AlignBuddy import AlignBuddy
from buddy_resources import GuessError, PhylipError


def test_can_create_align_buddy(alb_resources):
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

    for molecule in alb_resources.resource_list:  # dna, rna, pep
        for quantity in alb_resources.resource_list[molecule]:  # single, multi
            # rna resource only has a single alignment file so skip multi
            if molecule == 'rna' and quantity == 'multi':
                continue

            for file_type in alb_resources.resource_list[molecule][quantity]:  # clustal, fasta, etc
                assert bool(AlignBuddy(alb_resources.resource_list[molecule][quantity][file_type]))


def test_throws_errors_on_invalid_files(alignment_bad_resources):
    """ expect AlignBuddy to raise errors on invalid filesr """
    with pytest.raises(GuessError):
        AlignBuddy(alignment_bad_resources['dna']['single']['fasta'])
