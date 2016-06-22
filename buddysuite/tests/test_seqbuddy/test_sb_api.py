#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" tests basic functionality of AlignBuddy class """
import pytest
from Bio.SeqFeature import FeatureLocation, CompoundLocation

try:
    from buddysuite import SeqBuddy as Sb
    from buddysuite import MyFuncs
except ImportError:
    import SeqBuddy as Sb
    import MyFuncs


TEMPDIR = MyFuncs.TempDir()


# ##################### '-ano', '--annotate' ###################### ##
def test_annotate_pattern(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), pattern='α4')
    assert sb_helpers.seqs_to_hash(tester) == '4f9b3f82ce3c8d6c953f60a5a0e9574e'


def test_annotate_no_pattern(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100))
    assert sb_helpers.seqs_to_hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_compoundlocation(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', [(1, 100), (200, 250)])
    assert sb_helpers.seqs_to_hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_list_str(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', ['1-100', '(200-250)'])
    assert sb_helpers.seqs_to_hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_str(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', '1-100, (200-250)')
    assert sb_helpers.seqs_to_hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_fl_obj(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', FeatureLocation(start=0, end=100))
    assert sb_helpers.seqs_to_hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_cl_obj(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', CompoundLocation([FeatureLocation(start=0, end=100),
                                                                   FeatureLocation(start=199, end=250)],
                                                                  operator='order'))
    assert sb_helpers.seqs_to_hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_typerror(sb_resources):
    with pytest.raises(TypeError):
        tester = Sb.make_copy(sb_resources.get_one("d g"))
        Sb.annotate(tester, 'misc_feature', 5)


def test_annotate_pos_strand(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand='+')
    assert sb_helpers.seqs_to_hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_neg_strand(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand='-')
    assert sb_helpers.seqs_to_hash(tester) == '08524707f09d7eb273775f791d92964c'


def test_annotate_no_strand(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand=0)
    assert sb_helpers.seqs_to_hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_qualifier_dict(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers={'foo': 'bar', 'hello': 'world'})
    assert sb_helpers.seqs_to_hash(tester) == '34e9dfb9cfe62f0a4657c977eda45688'


def test_annotate_qualifier_list(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers=['foo=bar', 'hello=world'])
    assert sb_helpers.seqs_to_hash(tester) == '34e9dfb9cfe62f0a4657c977eda45688'


def test_annotate_qualifier_error(sb_resources):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    with pytest.raises(TypeError):
        Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers=tuple)


def test_annotate_out_of_range(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', [(-10, 100), (200, 10000)])
    assert sb_helpers.seqs_to_hash(tester) == 'a8f90863b2bbeaa519e8230a187532ca'

    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', [(1, 10000)])
    assert sb_helpers.seqs_to_hash(tester) == 'f5c90e3458fbca9b9565dac7877cc248'

    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', [(1, 10000), (20000, 30000)])
    assert sb_helpers.seqs_to_hash(tester) == 'f5c90e3458fbca9b9565dac7877cc248'

    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', FeatureLocation(start=-10, end=100))
    assert sb_helpers.seqs_to_hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_protein(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("p g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100))
    assert sb_helpers.seqs_to_hash(tester) == '93a248cacdaa1a58697c16827fe8709d'


def test_annotate_unrec_strand(capsys, sb_resources):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    Sb.annotate(tester, 'misc_feature', (1, 100), strand='foo')
    out, err = capsys.readouterr()
    assert err == "Warning: strand input not recognized. Value set to None."
    
    
# ######################  '-asl', '--ave_seq_length' ###################### #
def test_ave_seq_length_dna(sb_resources):
    for seqbuddy in sb_resources.get_list("d f g n pr s"):
        assert round(Sb.ave_seq_length(seqbuddy, clean=True), 2) == 1285.15


def test_ave_seq_length_pep(sb_resources):
    for seqbuddy in sb_resources.get_list("p f g n pr s"):
        assert round(Sb.ave_seq_length(seqbuddy, clean=True), 2) == 427.38


# ##################### '-cc', 'count_codons' ###################### ##
def test_count_codons_dna(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d f"))
    tester.records[0].seq = tester.records[0].seq[:-4]
    counter = Sb.count_codons(tester)
    assert sb_helpers.seqs_to_hash(counter[0]) == '4dbd9a5c68d85bb200c75b309fdaeeca'
    assert sb_helpers.string2hash(str(counter[1])) == '8d47313f3db02aee48575ff8ff4741b4'


def test_count_codons_rna(sb_resources, sb_helpers):
    tester = Sb.count_codons(sb_resources.get_one("r f"))[1]
    assert sb_helpers.string2hash(str(tester)) == 'b91daa8905533b5885d2067d9d6ffe36'


def test_count_codons_dna_badchar(sb_resources, sb_helpers):
    tester = Sb.count_codons(Sb.insert_sequence(sb_resources.get_one("d f"), 'PPP', -1))[1]
    assert sb_helpers.string2hash(str(tester)) == '9aba116675fe0e9eaaf43e5c6e0ba99d'


def test_count_codons_pep_exception(sb_resources):
    tester = Sb.make_copy(sb_resources.get_one("p f"))
    with pytest.raises(TypeError):
        Sb.count_codons(tester)


# ######################  '-efs', '--extract_feature_sequences' ###################### #
def test_extract_feature_sequences(sb_resources, sb_helpers):
    tester = sb_resources.get_one("d g")
    tester = Sb.extract_feature_sequences(tester, "CDS")
    assert sb_helpers.seqs_to_hash(tester) == "7e8a80caf902575c5eb3fc6ba8563956"

    tester = sb_resources.get_one("d g")
    tester = Sb.extract_feature_sequences(tester, ["TMD"])
    assert sb_helpers.seqs_to_hash(tester) == "13944b21484d5ea22af4fe57cc8074df"

    tester = sb_resources.get_one("d g")
    tester = Sb.extract_feature_sequences(tester, ["TMD", "splice_a"])
    assert sb_helpers.seqs_to_hash(tester) == "78629d308a89b458fb02e71d5568c978"

    tester = sb_resources.get_one("d g")
    tester = Sb.extract_feature_sequences(tester, "foo")
    assert sb_helpers.seqs_to_hash(tester) == "3cdbd5c8790f12871f8e04e40e315c93"

    tester = sb_resources.get_one("d g")
    tester = Sb.extract_feature_sequences(tester, [])
    assert sb_helpers.seqs_to_hash(tester) == "3cdbd5c8790f12871f8e04e40e315c93"
