# -*- coding: utf-8 -*-

"""
verify that fixtures are working as expected. I.e., test the tests before testing ;)
"""
import os
import pytest
from buddysuite import SeqBuddy, AlignBuddy, PhyloBuddy, DatabaseBuddy


# ToDo: Full coverage of the conftest.py file
class TestAlignmentResourceFixture:
    def test_alignment_valid_resources_has_values(self, alb_resources):
        """ checks that the alignment resource fixture has values """
        assert alb_resources
        assert isinstance(alb_resources.resources, dict)

        # check subject type
        assert alb_resources.resources['dna']
        assert alb_resources.resources['rna']
        assert alb_resources.resources['pep']

        assert alb_resources.alb_objs['dna']
        assert alb_resources.alb_objs['rna']
        assert alb_resources.alb_objs['pep']

        # check quantity types
        for molecule in alb_resources.resources:
            for quantity in ['multi', 'single']:
                # rna resource only has a single alignment file so skip multi
                if molecule == 'rna' and quantity == 'multi':
                    continue

                assert bool(alb_resources.resources[molecule][quantity])

    def test_alignment_valid_resources_files_exist(self, alb_resources):
        """ verifies that the alb_resources fixture points to real files """
        for molecule in alb_resources.resources:  # dna, rna, pep
            for quantity in ['multi', 'single']:
                # rna resource only has a single alignment file so skip multi
                if molecule == 'rna' and quantity == 'multi':
                    continue

                for _, path in alb_resources.resources[molecule][quantity].items():
                    assert os.path.isfile(path)

    def test_alignment_bad_resources_file_exists(self, alb_odd_resources):
        """ ensure that our bad test files exist """
        assert os.path.isfile(alb_odd_resources['dna']['single']['fasta'])


def test_sb_resources(sb_resources, capsys):
    assert sb_resources
    for molecule in [('dna', 11), ('rna', 2), ('pep', 8)]:
        assert molecule[0] in sb_resources.resources
        assert len(sb_resources.resources[molecule[0]]) == molecule[1]
        assert molecule[0] in sb_resources.sb_objs

    assert sb_resources.resources['dna']["clustal"] == "%s/Mnemiopsis_cds.clus" % sb_resources.res_path
    assert type(sb_resources.sb_objs['dna']["clustal"]) == SeqBuddy.SeqBuddy

    for key in ["molecule", "format"]:
        assert key in sb_resources.code_dict
    for key in ["p", "d", "r"]:
        assert key in sb_resources.code_dict["molecule"]
    for key in ["c", "f", "g", "n", "py", "pr", "pss", "psr", "s", "e", "x"]:
        assert key in sb_resources.code_dict["format"]

    for key in ["c", "f", "g", "n", "py", "pr", "pss", "psr", "s", "e", "x"]:
        assert key in sb_resources.single_letter_codes

    codes = sb_resources.parse_code()
    assert sorted(codes["molecule"]) == ['d', 'p', 'r']
    assert sorted(codes["format"]) == ['c', 'e', 'f', 'g', 'n', 'pr', 'psr', 'pss', 'py', 's', 'x']

    key = sb_resources.get_key("e d")
    assert key["molecule"] == "dna"
    assert key["format"] == "embl"

    with pytest.raises(AttributeError) as err:
        sb_resources.get_key("d e f")
    assert "Only explicit two-component codes are accepted" in str(err)

    with pytest.raises(AttributeError) as err:
        sb_resources.get_key("d z")
    assert "Malformed letter code, 'z' not recognized" in str(err)
    capsys.readouterr()
    with pytest.raises(AttributeError) as err:
        sb_resources.get_key("d p")
    assert "Malformed letter code, trying to append multiple values to the molecule component of the key" in str(err)


def test_sb_odd_resources(sb_odd_resources):
    for key in ["blank", "figtree", "unrecognizable", "gibberish", "phylipss_cols", "duplicate",
                "ambiguous_dna", "ambiguous_rna", "blastn", "blastp", "dummy_feats", "cnidaria_pep", "mixed"]:
        assert key in sb_odd_resources
