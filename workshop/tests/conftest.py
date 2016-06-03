# coding=utf-8
""" Fixtures for py.test  """
import os
import pytest
from collections import OrderedDict

resource_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'unit_test_resources')


@pytest.fixture()
def alignment_valid_resources():
    """ A dict of alignment file by molecule, quantity, and file_type

    alignment_valid_resources[<molecule_type>][<file_type>][<quantity type>]

        <molecule_type>:
            'dna', 'rna', or 'pep'

        <quantity>:
            'multi'
            'single'

        <file_type>:
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
    resource_list = {
        'dna': {
            'single': OrderedDict([(file_format, name.format(path=resource_path)) for file_format, name in [
                ("clustal", "{path}/Mnemiopsis_cds.clus"),
                ("fasta", "{path}/Mnemiopsis_cds_aln.fa"),
                ("gb", "{path}/Mnemiopsis_cds_aln.gb"),
                ("nexus", "{path}/Mnemiopsis_cds.nex"),
                ("phylip", "{path}/Mnemiopsis_cds.phy"),
                ("phylipr", "{path}/Mnemiopsis_cds.phyr"),
                ("phylipss", "{path}/Mnemiopsis_cds.physs"),
                ("phylipsr", "{path}/Mnemiopsis_cds.physr"),
                ("stockholm", "{path}/Mnemiopsis_cds.stklm")]]),

            'multi': OrderedDict([(file_format, name.format(path=resource_path)) for file_format, name in [
                ("clustal", "{path}/Alignments_cds.clus"),
                ("phylip", "{path}/Alignments_cds.phy"),
                ("phylipr", "{path}/Alignments_cds.phyr"),
                ("phylipss", "{path}/Alignments_cds.physs"),
                ("phylipsr", "{path}/Alignments_cds.physr"),
                ("stockholm", "{path}/Alignments_cds.stklm")]])
        },

        'rna': {
            'single': OrderedDict([("nexus", "{path}/Mnemiopsis_rna.nex".format(path=resource_path))]),
        },
        'pep': {
            'single': OrderedDict([(file_format, name.format(path=resource_path)) for file_format, name in [
                ("gb", "{path}/Mnemiopsis_pep_aln.gb"),
                ("nexus", "{path}/Mnemiopsis_pep.nex"),
                ("phylip", "{path}/Mnemiopsis_pep.phy"),
                ("phylipr", "{path}/Mnemiopsis_pep.phyr"),
                ("phylipss", "{path}/Mnemiopsis_pep.physs"),
                ("phylipsr", "{path}/Mnemiopsis_pep.physr"),
                ("stockholm", "{path}/Mnemiopsis_pep.stklm")]]),

            'multi': OrderedDict([(file_format, name.format(path=resource_path)) for file_format, name in [
                ("clustal", "{path}/Alignments_pep.clus"),
                ("phylip", "{path}/Alignments_pep.phy"),
                ("phylipr", "{path}/Alignments_pep.phyr"),
                ("phylipss", "{path}/Alignments_pep.physs"),
                ("phylipsr", "{path}/Alignments_pep.physr"),
                ("stockholm", "{path}/Alignments_pep.stklm")]])
        }
    }
    return resource_list


@pytest.fixture()
def alignment_bad_resources():
    """ A dict of invalid file resources """
    resource_list = {
        'dna': {
            'single': OrderedDict([(file_format, name.format(path=resource_path)) for file_format, name in [
                ('fasta', '{path}/gibberish.fa')]])
        },
    }
    return resource_list
