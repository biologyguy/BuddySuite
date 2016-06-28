#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" tests basic functionality of AlignBuddy class """
import pytest
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from unittest import mock
import os

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

# ######################  '-btr', '--back_translate' ###################### #
# Only fasta and genbank
hashes = [('p f', 'human', '1b14489a78bfe8255c777138877b9648'), ('p g', 'human', 'b6bcb4e5104cb202db0ec4c9fc2eaed2'),
          ('p f', 'yeast', '859ecfb88095f51bfaee6a1d1abeb50f'), ('p g', 'yeast', 'ba5c286b79a3514fba0b960ff81af25b'),
          ('p f', 'ecoli', '952a91a4506afb57f27136aa1f2a8af9'), ('p g', 'ecoli', '40c4a3e08c811b6bf3be8bedcb5d65a0'),
          ('p f', 'mouse', '3a3ee57f8dcde25c99a655494b218928'), ('p g', 'mouse', 'bc9f1ec6ec92c30b5053cd9bb6bb6f53')]


@pytest.mark.parametrize("key,organism,next_hash", hashes)
def test_back_translate(key, organism, next_hash, sb_resources, sb_helpers):
    tester = sb_resources.get_one(key)
    # tester.alpha = IUPAC.protein
    tester = Sb.back_translate(tester, 'OPTIMIZED', organism)
    assert sb_helpers.seqs_to_hash(tester) == next_hash


def test_back_translate_nucleotide_exception(sb_resources):
    with pytest.raises(TypeError):
        Sb.back_translate(sb_resources.get_one("d g"))


def test_back_translate_bad_mode(sb_resources):
    with pytest.raises(AttributeError):
        Sb.back_translate(Sb.make_copy(sb_resources.get_one("p f")), 'fgsdjkghjdalgsdf', 'human')


def test_back_translate_bad_organism(sb_resources):
    seqbuddy = Sb.make_copy(sb_resources.get_one("p f"))
    with pytest.raises(AttributeError):
        Sb.back_translate(seqbuddy, 'OPTIMIZED', 'fgsdjkghjdalgsdf')


# ######################  '-bl2s', '--bl2seq' ###################### #
# ToDo: Mock output data
def test_bl2seq(sb_resources, sb_helpers):
    pass


def test_bl2_no_binary(sb_resources):
    # noinspection PyUnresolvedReferences
    with mock.patch.dict(os.environ, {"PATH": ""}):
        with mock.patch('builtins.input', return_value="n"):
            with pytest.raises(RuntimeError):
                Sb.bl2seq(sb_resources.get_one("d f"))

            with pytest.raises(RuntimeError):
                Sb.bl2seq(sb_resources.get_one("p f"))


# ######################  '-bl', '--blast' ###################### #
def test_blastn(sb_resources, sb_odd_resources, sb_helpers):
    pass


def test_blastp(sb_resources, sb_odd_resources, sb_helpers):
    pass


# ######################  '-cs', '--clean_seq'  ###################### #
def test_clean_seq_prot(sb_resources, sb_helpers):
    # Protein
    tester = sb_resources.get_one("p f")
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert sb_helpers.seqs_to_hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert sb_helpers.seqs_to_hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"


def test_clean_seq_dna(sb_odd_resources, sb_helpers):
    # DNA
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_dna"])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert sb_helpers.seqs_to_hash(tester) == "71b28ad2730a9849f2ba0f70e9e51a9f"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert sb_helpers.seqs_to_hash(tester) == "1912fadb5ec52a38ec707c58085b86ad"
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_dna"])
    tester = Sb.clean_seq(tester, ambiguous=False, rep_char="X")
    assert sb_helpers.seqs_to_hash(tester) == "4c10ba4474d7484652cb633f03db1be1"


def test_clean_seq_rna(sb_odd_resources, sb_helpers):
    # RNA
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_rna"])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert sb_helpers.seqs_to_hash(tester) == "cdb1b963536d57efc7b7f87d2bf4ad22"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert sb_helpers.seqs_to_hash(tester) == "e66b785c649ad5086bcefd22e9ef9b41"
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_rna"])
    tester = Sb.clean_seq(tester, ambiguous=False, rep_char="X")
    assert sb_helpers.seqs_to_hash(tester) == "8ae19ab51b04076112d2f649353a4a79"


def test_clean_seq_align(sb_resources, sb_helpers):
    # Alignment formats are converted to fasta to prevent errors with sequence lengths
    for tester in sb_resources.get_list("d n py pr s"):
        tester = Sb.clean_seq(tester)
        sb_helpers.seqs_to_hash(tester) == "aa92396a9bb736ae6a669bdeaee36038"

# ######################  '-cmp', '--complement' ###################### #
hashes = [('d f', 'e4a358ca57aca0bbd220dc6c04c88795'), ('d g', '3366fcc6ead8f1bba4a3650e21db4ec3'),
          ('d n', '365bf5d08657fc553315aa9a7f764286'), ('d py', '520036b49dd7c70b9dbf4ce4d2c0e1d8'),
          ('d pr', 'dc1f7a3769a1e0b007969db1ab405e89'), ('d s', '5891348e8659290c2355fabd0f3ba4f4')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_complement(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.complement(sb_resources.get_one(key))
    assert sb_helpers.seqs_to_hash(tester) == next_hash


def test_complement_pep_exception(sb_resources):  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.complement(sb_resources.get_one("p f"))


# ######################  '-cts', '--concat_seqs' ###################### #
# ToDo: Test the _clean parameter
hashes = [('d f', '2e46edb78e60a832a473397ebec3d187'), ('d g', '7421c27be7b41aeedea73ff41869ac47'),
          ('d n', '494988ffae2ef3072c1619eca8a0ff3b'), ('d py', '710cad348c5560446daf2c916ff3b3e4'),
          ('d pr', '494988ffae2ef3072c1619eca8a0ff3b'), ('d s', '494988ffae2ef3072c1619eca8a0ff3b'),
          ('p f', '46741638cdf7abdf53c55f79738ee620'), ('p g', '8d0bb4e5004fb6a1a0261c30415746b5'),
          ('p n', '2651271d7668081cde8012db4f9a6574'), ('p py', '7846b2d080f09b60efc6ee43cd6d8502'),
          ('p pr', '5d1f8db03d6be30a7d77b00a0fba0b43'), ('p s', '2651271d7668081cde8012db4f9a6574')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_concat_seqs(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.concat_seqs(sb_resources.get_one(key))
    assert sb_helpers.seqs_to_hash(tester) == next_hash


def test_concat_seqs_clean(sb_resources, sb_helpers):
    tester = Sb.concat_seqs(sb_resources.get_one("d n"), clean=True)
    assert sb_helpers.seqs_to_hash(tester) == "2e46edb78e60a832a473397ebec3d187"


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
