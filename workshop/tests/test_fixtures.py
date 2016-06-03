# coding=utf-8
"""
verify that fixtures are working as expected. I.e., test the tests before testing ;)
"""
import os


class TestAlignmentResourceFixture:
    def test_alignment_valid_resources_has_values(self, alb_resources):
        """ checks that the alignment resource fixture has values """
        assert alb_resources
        assert isinstance(alb_resources.resource_list, dict)

        # check subject type
        assert alb_resources.resource_list['dna']
        assert alb_resources.resource_list['rna']
        assert alb_resources.resource_list['pep']

        # check quantity types
        for molecule in alb_resources.resource_list:
            for quantity in ['multi', 'single']:
                # rna resource only has a single alignment file so skip multi
                if molecule == 'rna' and quantity == 'multi':
                    continue

                assert bool(alb_resources.resource_list[molecule][quantity])

    def test_alignment_valid_resources_files_exist(self, alb_resources):
        """ verifies that the alb_resources fixture points to real files """
        for molecule in alb_resources.resource_list:  # dna, rna, pep
            for quantity in ['multi', 'single']:
                # rna resource only has a single alignment file so skip multi
                if molecule == 'rna' and quantity == 'multi':
                    continue

                for _, path in alb_resources.resource_list[molecule][quantity].items():
                    assert os.path.isfile(path)

    def test_alignment_bad_resources_file_exists(self, alignment_bad_resources):
        """ ensure that our bad test files exist """
        assert os.path.isfile(alignment_bad_resources['dna']['single']['fasta'])