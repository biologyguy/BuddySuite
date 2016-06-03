# coding=utf-8
""" verify that fixtures are working as expected """
import os


class TestAlignmentResourceFixture:
    def test_alignment_valid_resources_has_values(self, alignment_valid_resources):
        """ checks that the alignment resource fixture has values """
        assert alignment_valid_resources
        assert isinstance(alignment_valid_resources, dict)

        # check subject type
        assert alignment_valid_resources['dna']
        assert alignment_valid_resources['rna']
        assert alignment_valid_resources['pep']

        # check quantity types
        for molecule in alignment_valid_resources:
            for quantity in ['multi', 'single']:
                # rna resource only has a single alignment file so skip multi
                if molecule == 'rna' and quantity == 'multi':
                    continue

                assert bool(alignment_valid_resources[molecule][quantity])

    def test_alignment_valid_resources_files_exist(self, alignment_valid_resources):
        """ verifies that the alignment_valid_resources fixture points to real files """
        for molecule in alignment_valid_resources:  # dna, rna, pep
            for quantity in ['multi', 'single']:
                # rna resource only has a single alignment file so skip multi
                if molecule == 'rna' and quantity == 'multi':
                    continue

                for _, path in alignment_valid_resources[molecule][quantity].items():
                    assert os.path.isfile(path)

    def test_alignment_bad_resources_file_exists(self, alignment_bad_resources):
        """ ensure that our bad test files exist """
        assert os.path.isfile(alignment_bad_resources['dna']['single']['fasta'])