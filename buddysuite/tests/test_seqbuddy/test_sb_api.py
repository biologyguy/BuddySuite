#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" tests basic functionality of AlignBuddy class """
import pytest
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.Seq import Seq
from unittest import mock
import os
import urllib.request
import suds.client
import shutil
import time
from collections import OrderedDict

from ... import SeqBuddy as Sb
from ... import buddy_resources as br


TEMPDIR = br.TempDir()


# ##################### '-ano', '--annotate' ###################### ##
def test_annotate_pattern(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), pattern='α4')
    assert sb_helpers.seqs2hash(tester) == '4f9b3f82ce3c8d6c953f60a5a0e9574e'


def test_annotate_no_pattern(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100))
    assert sb_helpers.seqs2hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_compoundlocation(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', [(1, 100), (200, 250)])
    assert sb_helpers.seqs2hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_list_str(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', ['1-100', '(200-250)'])
    assert sb_helpers.seqs2hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_str(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', '1-100, (200-250)')
    assert sb_helpers.seqs2hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_fl_obj(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', FeatureLocation(start=0, end=100))
    assert sb_helpers.seqs2hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_cl_obj(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', CompoundLocation([FeatureLocation(start=0, end=100),
                                                                   FeatureLocation(start=199, end=250)],
                                                                  operator='order'))
    assert sb_helpers.seqs2hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_typerror(sb_resources):
    with pytest.raises(TypeError):
        tester = Sb.make_copy(sb_resources.get_one("d g"))
        Sb.annotate(tester, 'misc_feature', 5)


def test_annotate_pos_strand(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand='+')
    assert sb_helpers.seqs2hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_neg_strand(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand='-')
    assert sb_helpers.seqs2hash(tester) == '08524707f09d7eb273775f791d92964c'


def test_annotate_no_strand(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand=0)
    assert sb_helpers.seqs2hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_qualifier_dict(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers={'foo': 'bar', 'hello': 'world'})
    assert sb_helpers.seqs2hash(tester) == '34e9dfb9cfe62f0a4657c977eda45688'


def test_annotate_qualifier_list(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers=['foo=bar', 'hello=world'])
    assert sb_helpers.seqs2hash(tester) == '34e9dfb9cfe62f0a4657c977eda45688'


def test_annotate_qualifier_error(sb_resources):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    with pytest.raises(TypeError):
        Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers=tuple)


def test_annotate_out_of_range(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', [(-10, 100), (200, 10000)])
    assert sb_helpers.seqs2hash(tester) == 'a8f90863b2bbeaa519e8230a187532ca'

    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', [(1, 10000)])
    assert sb_helpers.seqs2hash(tester) == 'f5c90e3458fbca9b9565dac7877cc248'

    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', [(1, 10000), (20000, 30000)])
    assert sb_helpers.seqs2hash(tester) == 'f5c90e3458fbca9b9565dac7877cc248'

    tester = Sb.make_copy(sb_resources.get_one("d g"))
    tester = Sb.annotate(tester, 'misc_feature', FeatureLocation(start=-10, end=100))
    assert sb_helpers.seqs2hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_protein(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("p g"))
    tester = Sb.annotate(tester, 'misc_feature', (1, 100))
    assert sb_helpers.seqs2hash(tester) == '93a248cacdaa1a58697c16827fe8709d'


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
    assert sb_helpers.seqs2hash(tester) == next_hash


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
def test_bl2seq():
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
# ToDo: mock output data
def test_blastn():
    pass


def test_blastp():
    pass


# ######################  '-cs', '--clean_seq'  ###################### #
def test_clean_seq_prot(sb_resources, sb_helpers):
    # Protein
    tester = sb_resources.get_one("p f")
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert sb_helpers.seqs2hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert sb_helpers.seqs2hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"


def test_clean_seq_dna(sb_odd_resources, sb_helpers):
    # DNA
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_dna"])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert sb_helpers.seqs2hash(tester) == "71b28ad2730a9849f2ba0f70e9e51a9f"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert sb_helpers.seqs2hash(tester) == "1912fadb5ec52a38ec707c58085b86ad"
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_dna"])
    tester = Sb.clean_seq(tester, ambiguous=False, rep_char="X")
    assert sb_helpers.seqs2hash(tester) == "4c10ba4474d7484652cb633f03db1be1"


def test_clean_seq_rna(sb_odd_resources, sb_helpers):
    # RNA
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_rna"])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert sb_helpers.seqs2hash(tester) == "cdb1b963536d57efc7b7f87d2bf4ad22"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert sb_helpers.seqs2hash(tester) == "e66b785c649ad5086bcefd22e9ef9b41"
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_rna"])
    tester = Sb.clean_seq(tester, ambiguous=False, rep_char="X")
    assert sb_helpers.seqs2hash(tester) == "8ae19ab51b04076112d2f649353a4a79"


def test_clean_seq_align(sb_resources, sb_helpers):
    # Alignment formats are converted to fasta to prevent errors with sequence lengths
    for tester in sb_resources.get_list("d n py pr s"):
        tester = Sb.clean_seq(tester)
        sb_helpers.seqs2hash(tester) == "aa92396a9bb736ae6a669bdeaee36038"

# ######################  '-cmp', '--complement' ###################### #
hashes = [('d f', 'e4a358ca57aca0bbd220dc6c04c88795'), ('d g', '3366fcc6ead8f1bba4a3650e21db4ec3'),
          ('d n', '365bf5d08657fc553315aa9a7f764286'), ('d py', '520036b49dd7c70b9dbf4ce4d2c0e1d8'),
          ('d pr', 'dc1f7a3769a1e0b007969db1ab405e89'), ('d s', '5891348e8659290c2355fabd0f3ba4f4')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_complement(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.complement(sb_resources.get_one(key))
    assert sb_helpers.seqs2hash(tester) == next_hash


def test_complement_pep_exception(sb_resources):  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.complement(sb_resources.get_one("p f"))


# ######################  '-cts', '--concat_seqs' ###################### #
hashes = [('d f', '2e46edb78e60a832a473397ebec3d187'), ('d g', '7421c27be7b41aeedea73ff41869ac47'),
          ('d n', '494988ffae2ef3072c1619eca8a0ff3b'), ('d py', '710cad348c5560446daf2c916ff3b3e4'),
          ('d pr', '494988ffae2ef3072c1619eca8a0ff3b'), ('d s', '494988ffae2ef3072c1619eca8a0ff3b'),
          ('p f', '46741638cdf7abdf53c55f79738ee620'), ('p g', '8d0bb4e5004fb6a1a0261c30415746b5'),
          ('p n', '2651271d7668081cde8012db4f9a6574'), ('p py', '7846b2d080f09b60efc6ee43cd6d8502'),
          ('p pr', '5d1f8db03d6be30a7d77b00a0fba0b43'), ('p s', '2651271d7668081cde8012db4f9a6574')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_concat_seqs(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.concat_seqs(sb_resources.get_one(key))
    assert sb_helpers.seqs2hash(tester) == next_hash


hashes = [('d g', '7421c27be7b41aeedea73ff41869ac47'), ('d n', '2e46edb78e60a832a473397ebec3d187'),
          ('p g', '9cb5443d90e64693ad1bd74f29169ac5'), ('p py', 'e902a4faf44ebd28a43ca8103df7b828')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_concat_seqs_clean(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.concat_seqs(sb_resources.get_one(key), clean=True)
    tester.out_format = "gb"
    assert sb_helpers.seqs2hash(tester) == next_hash


# ##################### '-cc', 'count_codons' ###################### ##
def test_count_codons_dna(sb_resources, sb_helpers):
    tester = Sb.make_copy(sb_resources.get_one("d f"))
    tester.records[0].seq = tester.records[0].seq[:-4]
    counter = Sb.count_codons(tester)
    assert sb_helpers.seqs2hash(counter[0]) == '4dbd9a5c68d85bb200c75b309fdaeeca'
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


# ######################  '-cr', '--count_residues' ###################### #
def test_count_residues_unambig_dna(sb_resources):
    tester = Sb.SeqBuddy(">seq1\nACGCGAAGCGAACGCGCAGACGACGCGACGACGACGACGCA", in_format="fasta")
    tester = Sb.count_residues(tester).to_dict()
    assert tester["seq1"].buddy_data["res_count"]['A'] == [13, 0.3170731707317073]

    tester = sb_resources.get_one("d f")
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert res_count['G'] == [265, 0.21703521703521703]
    assert "% Ambiguous" not in res_count and "U" not in res_count


def test_count_residues_unambig_rna(sb_resources):
    tester = sb_resources.get_one("r f")
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert res_count['U'] == [356, 0.2915642915642916]
    assert "% Ambiguous" not in res_count and "T" not in res_count


def test_count_residues_ambig_dna(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources['ambiguous_dna'])
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert "U" not in res_count
    assert res_count['Y'] == [1, 0.000819000819000819]
    assert res_count['% Ambiguous'] == 0.98


def test_count_residues_ambig_rna(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources['ambiguous_rna'])
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert "T" not in res_count
    assert res_count['U'] == [353, 0.2891072891072891]
    assert res_count["Y"] == [1, 0.000819000819000819]
    assert res_count['% Ambiguous'] == 0.98


def test_count_residues_prot(sb_resources):
    tester = sb_resources.get_one("f p")
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert res_count['P'] == [17, 0.04176904176904177]
    assert res_count["G"] == [23, 0.056511056511056514]
    assert "% Ambiguous" not in res_count
    res_count = tester.to_dict()['Mle-Panxα8'].buddy_data["res_count"]
    assert res_count["% Ambiguous"] == 1.2
    assert res_count["% Positive"] == 12.23
    assert res_count["% Negative"] == 12.71
    assert res_count["% Uncharged"] == 73.62
    assert res_count["% Hyrdophilic"] == 36.93
    assert res_count["% Hyrdophobic"] == 55.4


# ######################  '-dgn' '--degenerate_sequence'################### #
def test_degenerate_sequence_without_arguments(sb_resources, sb_helpers):
    tester = sb_resources.get_one("f d")
    tester = Sb.degenerate_sequence(tester)
    assert sb_helpers.seqs2hash(tester) == '0638bc6546eebd9d50f771367d6d7855'


def test_degenerate_sequence_with_different_codon_tables(sb_resources, sb_helpers):
    for indx, next_hash in enumerate(['0638bc6546eebd9d50f771367d6d7855', '72373f8356051e2c6b67642451379054',
                                      '9172ad5947c0961b54dc5adbd03d4249', 'b45ac94ee6a98e495e115bfeb5bd9bcd',
                                      '76c45b4de8f7527b4139446b4551712b', 'baa5b48938cc5cae953c9083a5b21b12',
                                      '0ca67c4740fefbc7a20d806715c3ca12', 'd43ad8f328ff1d30eb1fb7bcd667a345',
                                      'd9d0f5cd8f0c25a0042527cc1cea802e', '4b9790f3f4eeeae1a9667b62b93bc961',
                                      '7ec4365c3571813d63cee4b70ba5dcf5']):
        tester = Sb.degenerate_sequence(sb_resources.get_one("d f"), table=(indx + 1))
        assert sb_helpers.seqs2hash(tester) == next_hash


def test_degenerate_sequence_edges(sb_resources, sb_odd_resources, sb_helpers):
    tester = sb_resources.get_one("p f")

    # Bad table reference
    with pytest.raises(KeyError) as e:
        Sb.degenerate_sequence(tester, 100)
    assert "Could not locate codon dictionary" in str(e.value)

    # Protein input
    with pytest.raises(TypeError) as e:
        Sb.degenerate_sequence(tester, 1)
    assert "Nucleic acid sequence required, not protein." in str(e.value)

    # RNA input
    tester = Sb.degenerate_sequence(Sb.SeqBuddy(sb_odd_resources["ambiguous_rna"]))
    assert sb_helpers.seqs2hash(tester) == "ff52e05971aeafd24c73a3b543901e4b"


# ######################  '-df', '--delete_features' ###################### #
def test_delete_features(sb_resources, sb_helpers):
    tester = sb_resources.get_one("d g")
    tester = Sb.delete_features(tester, 'donor')
    assert sb_helpers.seqs2hash(tester) == 'f84df6a77063c7def13babfaa0555bbf'


# ######################  '-dl', '--delete_large' ###################### #
def test_delete_large(sb_resources, sb_helpers):
    tester = sb_resources.get_one("d f")
    tester = Sb.delete_large(tester, 1285)
    assert sb_helpers.seqs2hash(tester) == '25859dc69d46651a1e04a70c07741b35'


# ######################  '-dm', '--delete_metadata' ###################### #
hashes = [('d f', 'aa92396a9bb736ae6a669bdeaee36038'), ('d g', '544ab887248a398d6dd1aab513bae5b1'),
          ('d n', 'cb1169c2dd357771a97a02ae2160935d'), ('d py', '503e23720beea201f8fadf5dabda75e4'),
          ('d pr', '52c23bd793c9761b7c0f897d3d757c12'), ('d s', 'a50943ccd028b6f5fa658178fa8cf54d'),
          ('p f', 'bac5dc724b1fee092efccd2845ff2513'), ('p g', '858e8475f7bc6e6a24681083a8635ef9'),
          ('p n', '17ff1b919cac899c5f918ce8d71904f6'), ('p py', '968ed9fa772e65750f201000d7da670f'),
          ('p pr', 'ce423d5b99d5917fbef6f3b47df40513'), ('p s', 'e224c16f6c27267b5f104c827e78df33')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_delete_metadata(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.delete_metadata(sb_resources.get_one(key))
    assert sb_helpers.seqs2hash(tester) == next_hash


# ######################  '-dr', '--delete_records' ###################### #
hashes = [('d f', '54bdb42423b1d331acea18218101e5fc'), ('d g', 'e2c03f1fa21fd27b2ff55f7f721a1a99'),
          ('d n', '6bc8a9409b1ef38e4f6f12121368883e'), ('d py', 'bda7be10061b0dcaeb66bebe3d736fee'),
          ('d pr', '6e2fce2239e2669b23f290049f87fbc4'), ('d s', '4c97144c5337f8a40c4fd494e622bf0d')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_delete_records(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.delete_records(sb_resources.get_one(key), 'α2')
    assert sb_helpers.seqs2hash(tester) == next_hash


def test_delete_records2(sb_resources, sb_helpers):
    tester = Sb.delete_records(sb_resources.get_one("d f"), ['α1', 'α2'])
    assert sb_helpers.seqs2hash(tester) == "eca4f181dae3d7998464ff71e277128f"

    tester = Sb.delete_records(sb_resources.get_one("d f"), 'α1|α2')
    assert sb_helpers.seqs2hash(tester) == "eca4f181dae3d7998464ff71e277128f"

    with pytest.raises(ValueError) as e:
        Sb.delete_records(sb_resources.get_one("d f"), dict)
    assert "'patterns' must be a list or a string." in str(e.value)


# #####################  '-drp', '--delete_repeats' ###################### ##
def test_delete_repeats(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources["duplicate"])
    tester = Sb.delete_repeats(tester)
    tester = Sb.find_repeats(tester)
    assert len(tester.repeat_ids) == 0
    assert len(tester.repeat_seqs) == 0


# ######################  '-ds', '--delete_small' ###################### #
def test_delete_small(sb_resources, sb_helpers):
    tester = sb_resources.get_one("d f")
    tester = Sb.delete_small(tester, 1285)
    assert sb_helpers.seqs2hash(tester) == '196adf08d4993c51050289e5167dacdf'


# ######################  '-d2r', '--transcribe' and 'r2d', '--back_transcribe' ###################### #
hashes = [('d f', 'd2db9b02485e80323c487c1dd6f1425b', 'b831e901d8b6b1ba52bad797bad92d14'),
          ('d g', '9ef3a2311a80f05f21b289ff7f401fff', '2e02a8e079267bd9add3c39f759b252c'),
          ('d n', 'f3bd73151645359af5db50d2bdb6a33d', 'cb1169c2dd357771a97a02ae2160935d'),
          ('d py', 'e55bd18b6d82a7fc3150338173e57e6a', '503e23720beea201f8fadf5dabda75e4'),
          ('d pr', 'a083e03b4e0242fa3c23afa80424d670', '52c23bd793c9761b7c0f897d3d757c12'),
          ('d s', '45b511f34653e3b984e412182edee3ca', '228e36a30e8433e4ee2cd78c3290fa6b')]


@pytest.mark.parametrize("key,d2r_hash,r2d_hash", hashes)
def test_transcribe(key, d2r_hash, r2d_hash, sb_resources, sb_helpers):
    tester = Sb.dna2rna(sb_resources.get_one(key))
    assert sb_helpers.seqs2hash(tester) == d2r_hash
    tester = Sb.rna2dna(tester)
    assert sb_helpers.seqs2hash(tester) == r2d_hash


def test_transcribe_pep_exception(sb_resources):  # Asserts that a ValueError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.dna2rna(sb_resources.get_one("p f"))


def test_back_transcribe_pep_exception(sb_resources):  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.rna2dna(sb_resources.get_one("p f"))


# ######################  '-efs', '--extract_feature_sequences' ###################### #
def test_extract_feature_sequences(sb_resources, sb_helpers):
    tester = sb_resources.get_one("d g")
    tester = Sb.extract_feature_sequences(tester, "CDS")
    assert sb_helpers.seqs2hash(tester) == "7e8a80caf902575c5eb3fc6ba8563956"

    tester = sb_resources.get_one("d g")
    tester = Sb.extract_feature_sequences(tester, ["TMD"])
    assert sb_helpers.seqs2hash(tester) == "13944b21484d5ea22af4fe57cc8074df"

    tester = sb_resources.get_one("d g")
    tester = Sb.extract_feature_sequences(tester, ["TMD", "splice_a"])
    assert sb_helpers.seqs2hash(tester) == "78629d308a89b458fb02e71d5568c978"

    tester = sb_resources.get_one("d g")
    tester = Sb.extract_feature_sequences(tester, "foo")
    assert sb_helpers.seqs2hash(tester) == "3cdbd5c8790f12871f8e04e40e315c93"

    tester = sb_resources.get_one("d g")
    tester = Sb.extract_feature_sequences(tester, [])
    assert sb_helpers.seqs2hash(tester) == "3cdbd5c8790f12871f8e04e40e315c93"


# ######################  '-er', '--extract_regions' ###################### #
hashes = [('d f', '8c2fac57aedf6b0dab3d0f5bcf88e99f'), ('d g', '25ad9670e8a6bac7962ab46fd79251e5'),
          ('d n', '4063ab66ced2fafb080ceba88965d2bb'), ('d py', '33e6347792aead3c454bac0e05a292c6'),
          ('d pr', '9a5c491aa293c6cedd48c4c249d55aff'), ('d s', 'cd8d857feba9b6e459b8a9d56f11b7f5'),
          ('p f', '2586d1e3fc283e6f5876251c1c57efce'), ('p g', 'a9c22659967916dcdae499c06d1aaafb'),
          ('p n', '6a27222d8f60ee8496cbe0c41648a116'), ('p py', 'c9a1dd913190f95bba5eca6a89685c75'),
          ('p pr', '6f579144a43dace285356ce6eb326d3b'), ('p s', '727099e0abb89482760eeb20f7edd0cd')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_extract_regions_multiformat(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.extract_regions(sb_resources.get_one(key), "50:300")
    assert sb_helpers.seqs2hash(tester) == next_hash

    tester = Sb.extract_regions(sb_resources.get_one(key), "300:50")
    assert sb_helpers.seqs2hash(tester) == next_hash


def test_extract_regions_singlets(sb_resources, sb_helpers):
    tester = Sb.extract_regions(sb_resources.get_one("p g"), "0")
    assert sb_helpers.seqs2hash(tester) == "73e7f4d6eafba9e3fa57cf21be46ca62"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "1")
    assert sb_helpers.seqs2hash(tester) == "73e7f4d6eafba9e3fa57cf21be46ca62"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "-10000000")
    assert sb_helpers.seqs2hash(tester) == "73e7f4d6eafba9e3fa57cf21be46ca62"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), ",1/")
    assert sb_helpers.seqs2hash(tester) == "73e7f4d6eafba9e3fa57cf21be46ca62"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "1000000")
    assert sb_helpers.seqs2hash(tester) == "998b272614ceb172ea57137b66d2669d"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "2,5,9,-5")
    assert sb_helpers.seqs2hash(tester) == "8540d5c2f6d7ef050886f1192da1396f"


def test_extract_regions_ranges(sb_resources, sb_helpers):
    tester = Sb.extract_regions(sb_resources.get_one("p g"), "0:10")
    assert sb_helpers.seqs2hash(tester) == "083702438b2577c414eba3c812443249"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "1:10")
    assert sb_helpers.seqs2hash(tester) == "083702438b2577c414eba3c812443249"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "10:1")
    assert sb_helpers.seqs2hash(tester) == "083702438b2577c414eba3c812443249"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), ":10")
    assert sb_helpers.seqs2hash(tester) == "083702438b2577c414eba3c812443249"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "-10:")
    assert sb_helpers.seqs2hash(tester) == "e3f2b1995d47a08429767c559f41691c"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "40:75,89:100,432:-45")
    assert sb_helpers.seqs2hash(tester) == "6059e115128c8f211f7e41c2745b5d34"


def test_extract_regions_mth_of_nth(sb_resources, sb_helpers):
    tester = Sb.extract_regions(sb_resources.get_one("p g"), "1/50")
    assert sb_helpers.seqs2hash(tester) == "869e2e07dfeae4e174aa9bb4c135ff25"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "-1/50")
    assert sb_helpers.seqs2hash(tester) == "0a893feff83f6dc2f1f105c01e162409"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "1/-50")
    assert sb_helpers.seqs2hash(tester) == "4e1f1c8c6caacac7e65dbf152807baa1"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "50/1")
    assert sb_helpers.seqs2hash(tester) == "87e7701ff6ab2da7b6cd46e5cc48e0a3"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "50/25")
    assert sb_helpers.seqs2hash(tester) == "20239fbef24e40b903a75f941b43b9c8"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "1:5/50")
    assert sb_helpers.seqs2hash(tester) == "cdb538506a2db396b8d67d2864f8109f"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "-5:/50")
    assert sb_helpers.seqs2hash(tester) == "368dac090406e6e0fcf6d025f5c0e069"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), ":5/50")
    assert sb_helpers.seqs2hash(tester) == "cdb538506a2db396b8d67d2864f8109f"

    tester = Sb.extract_regions(sb_resources.get_one("p g"), "1:10,1/50,-1")
    assert sb_helpers.seqs2hash(tester) == "9204dd879fa59cf35a253fb0ff82758c"


def test_extract_regions_edges(sb_resources):
    with pytest.raises(ValueError) as err:
        Sb.extract_regions(sb_resources.get_one("p g"), "foo")
    assert "Unable to decode the positions string 'foo'" in str(err)


# #####################  '-fcpg', '--find_CpG' ###################### ##
def test_find_cpg(sb_resources, sb_helpers):
    tester = sb_resources.get_one("d g")
    tester = Sb.find_cpg(tester)
    assert sb_helpers.seqs2hash(tester) == "9499f524da0c35a60502031e94864928"


# #####################  '-orf', '--find_orf' ###################### ##
def test_find_orf_basic(sb_helpers):
    tester = Sb.SeqBuddy("ATGAAATTTCCCGGGTAG", in_format='raw', out_format='gb')
    tester = Sb.find_orfs(tester)
    assert sb_helpers.seqs2hash(tester) == "100c4f29e758b33b8803f46dba1a75ef"
    assert sb_helpers.features2hash(tester) == "e47d6ea2a9572729024b55258411a579"


def test_find_orf_standard(sb_resources, sb_helpers):
    tester = Sb.find_orfs(sb_resources.get_one("d g"))
    assert sb_helpers.seqs2hash(tester) == "4474e4d21c02eb171e0cce128fb6bf5e"
    assert sb_helpers.features2hash(tester) == "f8abc7372fe7a8bbfdfce715f6bf6f0e"
    tester.out_format = "fasta"
    assert sb_helpers.seqs2hash(tester) == "25073539df4a982b7f99c72dd280bb8f"
    assert sb_helpers.features2hash(tester) == "f8abc7372fe7a8bbfdfce715f6bf6f0e"
    tester = Sb.find_orfs(sb_resources.get_one("d g"), include_feature=False)
    assert sb_helpers.seqs2hash(tester) == "908744b00d9f3392a64b4b18f0db9fee"
    assert sb_helpers.features2hash(tester) == "969ab7e60bba0d958ed17f6980b19059"
    tester = Sb.find_orfs(sb_resources.get_one("r f"))
    assert sb_helpers.seqs2hash(tester) == "013ebe2bc7d83c44f58344b865e1f55b"
    assert sb_helpers.features2hash(tester) == "a80f082f0649bb037c34d0f268fec901"
    tester.out_format = "gb"
    assert sb_helpers.seqs2hash(tester) == "76287ca18781224027e9ed8327a97558"
    assert sb_helpers.features2hash(tester) == "a80f082f0649bb037c34d0f268fec901"
    tester = Sb.find_orfs(sb_resources.get_one("r f"), include_feature=False)
    assert sb_helpers.features2hash(tester) == "aa0c0c27a9b8cde585da2c42a9280f3c"
    tester = sb_resources.get_one("p g")
    with pytest.raises(TypeError):
        Sb.find_orfs(tester)


# #####################  '-fp', '--find_pattern' ###################### ##
def test_find_pattern(sb_resources, sb_helpers):
    tester = Sb.find_pattern(sb_resources.get_one("d g"), "ATGGT")
    assert sb_helpers.seqs2hash(tester) == "ca129f98c6c719d50f0cf43eaf6dc90a"
    tester.out_format = "fasta"
    assert sb_helpers.seqs2hash(tester) == "6f23f80b52ffb736bbecc9f4c72d8fab"
    tester = Sb.find_pattern(sb_resources.get_one("d g"), "ATg{2}T")
    assert sb_helpers.seqs2hash(tester) == "9ec8561c264bff6f7166855d60457df1"
    tester = Sb.find_pattern(sb_resources.get_one("d g"), "ATg{2}T", "tga.{1,6}tg")
    assert sb_helpers.seqs2hash(tester) == "ec43ce98c9ae577614403933b2c5f37a"
    tester = Sb.find_pattern(sb_resources.get_one("d g"), "ATg{2}T", "tga.{1,6}tg", include_feature=False)
    assert sb_helpers.seqs2hash(tester) == "2e02a8e079267bd9add3c39f759b252c"
    tester = Sb.find_pattern(sb_resources.get_one("p g"), "[bz]{2}x{50,100}[bz]{2}", ambig=True)
    assert sb_helpers.seqs2hash(tester) == "339ff26803a2d12267d873458d40bca2"
    tester = Sb.find_pattern(sb_resources.get_one("d g"), "ATGGN{6}", ambig=True)
    assert sb_helpers.seqs2hash(tester) == "ac9adb42fbfa9cf22f033e9a02130985"
    tester = Sb.find_pattern(sb_resources.get_one("r f"), "AUGGN{6}", ambig=True)
    assert sb_helpers.seqs2hash(tester) == "b7abcb4334232e38dfbac9f46234501a"


# #####################  '-frp', '--find_repeats' ###################### ##
def test_find_repeats(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources["duplicate"])
    tester.unique_seqs, tester.repeat_ids, tester.repeat_seqs = {}, {}, {}
    Sb.find_repeats(tester)
    assert 'Seq1' in tester.unique_seqs
    assert 'Seq12' in tester.repeat_ids

    for key in tester.repeat_seqs:
        assert 'Seq12' in tester.repeat_seqs[key] or 'Seq10A' in tester.repeat_seqs[key]


# ######################  '-frs', '--find_restriction_sites' ###################### #
def test_restriction_sites_no_args(sb_resources, sb_helpers):
    # No arguments passed in = commercial REs and any number of cut sites
    tester = Sb.find_restriction_sites(sb_resources.get_one("d g"))
    assert sb_helpers.seqs2hash(tester) == 'a48fc20dc07b6bf03b0cef32ed27c5d2'
    assert sb_helpers.string2hash(str(tester.restriction_sites)) == "646d1026fc5b245ad7130dab3f027489"


def test_restriction_sites_all_emzymes(sb_resources, sb_helpers):
    # All enzymes
    tester = Sb.find_restriction_sites(sb_resources.get_one("d g"), enzyme_group=["all"])
    assert sb_helpers.seqs2hash(tester) == '52a175a264a2c101dca7dcbc9e8d01f0'
    assert sb_helpers.string2hash(str(tester.restriction_sites)) == "3c2d3fbabfeae48c9ec4cbfa7f67afb1"


def test_restriction_sites_limit_cuts(capsys, sb_resources, sb_helpers):
    # Specify a few REs and limit the number of cuts
    tester = Sb.find_restriction_sites(sb_resources.get_one("d g"), min_cuts=2, max_cuts=4,
                                       enzyme_group=["EcoRI", "KspI", "TasI", "Bme1390I", "FooBR"])
    out, err = capsys.readouterr()
    assert sb_helpers.seqs2hash(tester) == 'c42b3bf0367557383000b897432fed2d'
    assert sb_helpers.string2hash(str(tester.restriction_sites)) == "0d2e5fdba6fed434495481397a91e56a"
    assert "Warning: FooBR not a known enzyme" in err

    with pytest.raises(TypeError) as e:
        Sb.find_restriction_sites(sb_resources.get_one("p g"))
    assert str(e.value) == "Unable to identify restriction sites in protein sequences."

    with pytest.raises(ValueError) as e:
        Sb.find_restriction_sites(tester, min_cuts=4, max_cuts=2)
    assert str(e.value) == "min_cuts parameter has been set higher than max_cuts."

    # 2-cutters and non-cutters
    Sb.find_restriction_sites(tester, enzyme_group=["AjuI", "AlwFI"])
    out, err = capsys.readouterr()
    assert "Warning: Double-cutters not supported." in err
    assert "Warning: No-cutters not supported." in err


# ######################  '-hsi', '--hash_sequence_ids' ###################### #
def test_hash_seq_ids(sb_resources):
    tester = sb_resources.get_one("d f")
    Sb.hash_ids(tester)
    assert len(tester.records[0].id) == 10

    tester = sb_resources.get_one("d f")
    tester = Sb.hash_ids(tester, 25)
    assert len(tester.records[0].id) == 25
    assert len(tester.hash_map) == 13

    tester_copy = Sb.make_copy(tester)
    tester_copy.reverse_hashmap()
    tester_copy = Sb.hash_ids(tester_copy, 25)
    assert tester_copy.records[0].id == tester.records[0].id

    tester_copy.reverse_hashmap()
    tester_copy.records[0].id = tester_copy.records[0].id[:-1]
    tester_copy = Sb.hash_ids(tester_copy, 25)
    assert tester_copy.records[0].id != tester.records[0].id


def test_hash_seq_ids_errors(sb_resources):
    tester = sb_resources.get_one("d f")
    with pytest.raises(TypeError) as e:
        Sb.hash_ids(tester, "foo")
    assert str(e.value) == "Hash length argument must be an integer, not <class 'str'>"

    with pytest.raises(ValueError) as e:
        Sb.hash_ids(tester, 0)
    assert str(e.value) == "Hash length must be greater than 0"

    tester.records *= 10
    with pytest.raises(ValueError) as e:
        Sb.hash_ids(tester, 1)
    assert "Insufficient number of hashes available to cover all sequences." in str(e.value)


# ##################### '-is', 'insert_seq' ###################### ##
def test_insert_seqs_start(sb_resources, sb_helpers):
    insert = 'AACAGGTCGAGCA'
    tester = sb_resources.get_one("d f")
    assert sb_helpers.seqs2hash(Sb.insert_sequence(tester, insert)) == 'f65fee08b892af5ef93caa1bf3cb3980'

    tester = sb_resources.get_one("d f")
    assert sb_helpers.seqs2hash(Sb.insert_sequence(tester, insert, -9000)) == 'f65fee08b892af5ef93caa1bf3cb3980'

    tester = sb_resources.get_one("d f")
    assert sb_helpers.seqs2hash(Sb.insert_sequence(tester, insert, -1)) == '792397e2e32e95b56ddc15b8b2310ec0'

    tester = sb_resources.get_one("d f")
    assert sb_helpers.seqs2hash(Sb.insert_sequence(tester, insert, 9000)) == '792397e2e32e95b56ddc15b8b2310ec0'

    tester = sb_resources.get_one("d f")
    Sb.insert_sequence(tester, insert, 100, ["α[23]", "α5"])
    assert sb_helpers.seqs2hash(tester) == 'edcd7934eb026ac3ea4b603ac85ca79f'

    tester = sb_resources.get_one("d f")
    assert sb_helpers.seqs2hash(Sb.insert_sequence(tester, insert, -25)) == '29cab1e72ba95572c3aec469270071e9'


# ######################  '-ip', '--isoelectric_point' ###################### #
def test_isoelectric_point(sb_resources, sb_helpers):
    for tester in sb_resources.get_list("p f g n pr s"):
        tester = Sb.isoelectric_point(Sb.clean_seq(tester))
        assert tester.to_dict()["Mle-Panxα12"].features[-1].qualifiers["value"] == 6.0117797852
        if tester.out_format == "gb":
            assert sb_helpers.seqs2hash(tester) == "8bc299e31f436d192bf8cf8b7af671a8"

        with pytest.raises(TypeError):
            Sb.isoelectric_point(sb_resources.get_one("d f"))


# ######################  '-lc', '--lowercase' and 'uc', '--uppercase'  ###################### #
hashes = [('d f', '25073539df4a982b7f99c72dd280bb8f', 'b831e901d8b6b1ba52bad797bad92d14'),
          ('d g', '2e02a8e079267bd9add3c39f759b252c', '2e02a8e079267bd9add3c39f759b252c'),
          ('d n', '52e74a09c305d031fc5263d1751e265d', 'cb1169c2dd357771a97a02ae2160935d'),
          ('d py', 'cfe6cb9c80aebd353cf0378a6d284239', '503e23720beea201f8fadf5dabda75e4'),
          ('d pr', '6e5542f41d17ff33afb530b4d07408a3', '52c23bd793c9761b7c0f897d3d757c12'),
          ('d s', 'b82538a4630810c004dc8a4c2d5165ce', '228e36a30e8433e4ee2cd78c3290fa6b'),
          ('p f', 'c10d136c93f41db280933d5b3468f187', '14227e77440e75dd3fbec477f6fd8bdc'),
          ('p g', '7a8e25892dada7eb45e48852cbb6b63d', '7a8e25892dada7eb45e48852cbb6b63d'),
          ('p n', '8b6737fe33058121fd99d2deee2f9a76', '17ff1b919cac899c5f918ce8d71904f6'),
          ('p py', '968ed9fa772e65750f201000d7da670f', 'aacda2f5d4077f23926400f74afa2f46'),
          ('p pr', 'ce423d5b99d5917fbef6f3b47df40513', 'e3dc2e0347f40fffec45d053f4f34c96'),
          ('p s', 'f35cbc6e929c51481e4ec31e95671638', 'c0dce60745515b31a27de1f919083fe9')]


@pytest.mark.parametrize("key,uc_hash,lc_hash", hashes)
def test_cases(key, uc_hash, lc_hash, sb_resources, sb_helpers):  # NOTE: Biopython always writes genbank to lower case
    tester = Sb.uppercase(sb_resources.get_one(key))
    assert sb_helpers.seqs2hash(tester) == uc_hash
    tester = Sb.lowercase(tester)
    assert sb_helpers.seqs2hash(tester) == lc_hash


# ######################  '-mui', '--make_ids_unique' ###################### #
def test_make_ids_unique(sb_odd_resources, sb_helpers):
    tester = Sb.SeqBuddy(sb_odd_resources["duplicate"])
    Sb.make_ids_unique(tester)
    assert sb_helpers.seqs2hash(tester) == "363c7ed14be59bcacede092b8f334a52"

    tester = Sb.SeqBuddy(sb_odd_resources["duplicate"])
    Sb.make_ids_unique(tester, sep="-", padding=4)
    assert sb_helpers.seqs2hash(tester) == "0054df3003ba16287159147f3b85dc7b"

# ######################  '-fn2p', '--map_features_nucl2prot' ###################### #
# Map the genbank DNA file to all protein files, and the fasta DNA file to fasta protein
hashes = [('p f', '5216ef85afec36d5282578458a41169a'), ('p g', 'a8f7c129cf57a746c20198bf0a6b9cf4'),
          ('p n', 'bb0c9da494b5418fb87862dab2a66cfa'), ('p py', '3c0e3ec45abd774813a274fda1b4a5f2'),
          ('p pr', 'a7f6c4bb410f17cfc3e8966ccbe3e065'), ('p s', '1b8c44f4ace877b568c0915033980bed')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_map_features_nucl2prot(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.map_features_nucl2prot(sb_resources.get_one("d g"), sb_resources.get_one(key))
    tester.out_format = "gb"
    assert sb_helpers.seqs2hash(tester) == next_hash

    if key == "p f":
        tester = Sb.map_features_nucl2prot(sb_resources.get_one("d f"), sb_resources.get_one(key))
        tester.out_format = "gb"
        assert sb_helpers.seqs2hash(tester) == "854566b485af0f277294bbfb15f7dd0a"


def test_map_features_nucl2prot_2(capsys, sb_resources, sb_helpers):
    tester = sb_resources.get_one("d g")
    Sb.pull_record_ends(tester, 1300)
    tester = Sb.annotate(tester, "foo", [(20, 40), (50, 60)], pattern="α9")
    mapped = Sb.map_features_nucl2prot(Sb.make_copy(tester), sb_resources.get_one("p f"))
    mapped.out_format = "gb"
    assert sb_helpers.seqs2hash(mapped) == "807025489aadf98d851501f49d463e4a"
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(err) == "6c840e4acaaf4328672ca164f854000a"

    prot_tester = sb_resources.get_one("p g")
    prot_tester.records = sorted(prot_tester.records, key=lambda x: x.id)
    dna_tester = sb_resources.get_one("d g")
    dna_tester.records = sorted(dna_tester.records, key=lambda x: x.id)
    Sb.rename(dna_tester, "α4", "A4")

    Sb.map_features_nucl2prot(Sb.make_copy(dna_tester), prot_tester, mode="list")
    assert sb_helpers.seqs2hash(prot_tester) == "9fdb606ea65d6c050540a94137ae6e0d"

    Sb.map_features_nucl2prot(Sb.make_copy(dna_tester), prot_tester, mode="key")
    assert sb_helpers.seqs2hash(prot_tester) == "9fdb606ea65d6c050540a94137ae6e0d"
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(err) == "6fd2b5f2a7a3995d3f49c4919c3358b0"

    mock_obj = type('test', (object,), {})()
    mock_obj.start = 2
    mock_obj.end = 6
    dna_tester.records[0].features[0].location = mock_obj
    with pytest.raises(TypeError) as e:
        Sb.map_features_nucl2prot(Sb.make_copy(dna_tester), sb_resources.get_one("p f"))
    assert "FeatureLocation or CompoundLocation object required." in str(e.value)

    tester.records[0].features[0].location = tester.records[1].features[0].location
    with pytest.raises(ValueError) as e:
        Sb.map_features_nucl2prot(Sb.make_copy(tester), sb_resources.get_one("p f"), mode="foo")
    assert "'mode' must be either 'key' or 'position'" in str(e.value)

    with pytest.raises(ValueError) as e:
        Sb.pull_recs(tester, "α[1-8]")
        Sb.map_features_nucl2prot(Sb.make_copy(tester), sb_resources.get_one("p f"), mode="list")
    assert "The two input files do not contain the same number of sequences" in str(e.value)

# ######################  '-fp2n', '--map_features_prot2nucl' ###################### #
hashes = [('d f', '3ebc92ca11505489cab2453d2ebdfcf2'), ('d g', 'feceaf5e17935afb100b4b6030e27fee'),
          ('d n', 'bfd36942768cf65c473b3aaebb83e4fa'), ('d py', '9ba4af4e5dd0bf4a445d173604b92996'),
          ('d pr', 'c178763aa9596e341bbbc088f1f791c9'), ('d s', '84cc7ecb54603c5032737e5263a52bd3')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_map_features_prot2nucl(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.map_features_prot2nucl(sb_resources.get_one("p g"), sb_resources.get_one(key))
    tester.out_format = "gb"
    assert sb_helpers.seqs2hash(tester) == next_hash


def test_map_features_prot2nucl_2(capsys, sb_resources, sb_helpers):
    prot_tester = sb_resources.get_one("p g")
    Sb.pull_record_ends(prot_tester, 450)
    prot_tester = Sb.annotate(prot_tester, "foo", [(20, 40), (50, 60)], pattern="α9")
    mapped = Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), sb_resources.get_one("d f"))
    mapped.out_format = "gb"
    assert sb_helpers.seqs2hash(mapped) == "552b31b8a7068d12d6c55c7e5d293c54"
    out, err = capsys.readouterr()
    assert err == "Warning: size mismatch between aa and nucl seqs for Mle-Panxα7A --> 450, 1875\n"

    prot_tester = sb_resources.get_one("p g")
    prot_tester.records = sorted(prot_tester.records, key=lambda x: x.id)
    dna_tester = sb_resources.get_one("d g")
    dna_tester.records = sorted(dna_tester.records, key=lambda x: x.id)
    Sb.rename(prot_tester, "α4", "A4")
    Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), dna_tester, mode="list")
    assert sb_helpers.seqs2hash(dna_tester) == "c32f43cc5205867c0eb1d3873e27319b"

    Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), dna_tester, mode="key")
    assert sb_helpers.seqs2hash(dna_tester) == "c32f43cc5205867c0eb1d3873e27319b"
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(err) == "c1f51587b07d2cc6156b8aac07384834"

    mock_obj = type('test', (object,), {})()
    mock_obj.start = 2
    mock_obj.end = 6
    prot_tester.records[0].features[0].location = mock_obj
    with pytest.raises(TypeError) as e:
        Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), sb_resources.get_one("d f"))
    assert "FeatureLocation or CompoundLocation object required." in str(e.value)

    prot_tester.records[0].features[0].location = prot_tester.records[1].features[0].location
    with pytest.raises(ValueError) as e:
        Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), sb_resources.get_one("p f"), mode="foo")
    assert "'mode' must be either 'key' or 'position'" in str(e.value)

    with pytest.raises(ValueError) as e:
        Sb.pull_recs(prot_tester, "α[1-8]")
        Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), sb_resources.get_one("p f"), mode="list")
    assert "The two input files do not contain the same number of sequences" in str(e.value)


# #####################  '-mg', '--merge' ###################### ##
def test_merge(sb_resources, sb_odd_resources, sb_helpers):
    tester = Sb.SeqBuddy(sb_odd_resources['dummy_feats'])
    tester = Sb.merge(tester, sb_resources.get_one("d g"))
    assert sb_helpers.seqs2hash(tester) == "bae5aeb130b3d5319378a122a6f61df5"

    tester.records[0].seq = tester.records[1].seq
    with pytest.raises(RuntimeError) as e:
        Sb.merge(tester, sb_resources.get_one("d g"))
    assert "Sequence mismatch for record 'Mle-Panxα1'" in str(e.value)


# ######################  '-mw', '--molecular_weight' ###################### #
def test_molecular_weight(sb_resources, sb_odd_resources, sb_helpers):
    # Unambiguous DNA
    tester = Sb.molecular_weight(sb_resources.get_one("d g"))
    assert tester.molecular_weights['masses_ds'][0] == 743477.1
    assert tester.molecular_weights['masses_ss'][0] == 371242.6
    assert sb_helpers.seqs2hash(tester) == "e080cffef0ec6c5e8eada6f57bbc35f9"
    # Ambiguous DNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(sb_odd_resources['ambiguous_dna']))
    assert tester.molecular_weights['masses_ds'][0] == 743477.08
    assert tester.molecular_weights['masses_ss'][0] == 371202.59
    # Unambiguous RNA
    tester = Sb.molecular_weight(sb_resources.get_one("r f"))
    assert tester.molecular_weights['masses_ss'][0] == 387372.6
    # Ambiguous RNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(sb_odd_resources['ambiguous_rna']))
    assert tester.molecular_weights['masses_ss'][0] == 387371.6
    # Protein
    tester = Sb.molecular_weight(sb_resources.get_one("p g"))
    assert tester.molecular_weights['masses_ss'][0] == 45692.99
    assert sb_helpers.seqs2hash(tester) == "fb1a66b7eb576c0584fc7988c45b6a18"


# ######################  '-ns', '--num_seqs' ###################### #
def test_num_seqs(sb_resources):
    for tester in sb_resources.get_list("d p f g pr n s"):
        assert Sb.num_seqs(tester) == 13
    for tester in sb_resources.get_list("d p py ps"):
        assert Sb.num_seqs(tester) == 8


def test_empty_file(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources["blank"])
    assert type(tester) == Sb.SeqBuddy
    assert len(tester.records) == 0

# ######################  '-ofa', '--order_features_alphabetically' ###################### #
hashes = [('d f', 'b831e901d8b6b1ba52bad797bad92d14', 'b831e901d8b6b1ba52bad797bad92d14'),
          ('d g', '21547b4b35e49fa37e5c5b858808befb', '3b718ec3cb794bcb658d900e517110cc'),
          ('d n', 'cb1169c2dd357771a97a02ae2160935d', 'cb1169c2dd357771a97a02ae2160935d'),
          ('d py', '503e23720beea201f8fadf5dabda75e4', '503e23720beea201f8fadf5dabda75e4'),
          ('d pr', '52c23bd793c9761b7c0f897d3d757c12', '52c23bd793c9761b7c0f897d3d757c12'),
          ('d s', '228e36a30e8433e4ee2cd78c3290fa6b', '228e36a30e8433e4ee2cd78c3290fa6b'),
          ('p f', '14227e77440e75dd3fbec477f6fd8bdc', '14227e77440e75dd3fbec477f6fd8bdc'),
          ('p g', 'd0297078b4c480a49b6da5b719310d0e', 'c6a788d8ea916964605ac2942c459c9b'),
          ('p n', '17ff1b919cac899c5f918ce8d71904f6', '17ff1b919cac899c5f918ce8d71904f6'),
          ('p py', '968ed9fa772e65750f201000d7da670f', '968ed9fa772e65750f201000d7da670f'),
          ('p pr', 'ce423d5b99d5917fbef6f3b47df40513', 'ce423d5b99d5917fbef6f3b47df40513'),
          ('p s', 'c0dce60745515b31a27de1f919083fe9', 'c0dce60745515b31a27de1f919083fe9')]


@pytest.mark.parametrize("key,fwd_hash,rev_hash", hashes)
def test_order_features_alphabetically(key, fwd_hash, rev_hash, sb_resources, sb_helpers):
    tester = Sb.order_features_alphabetically(sb_resources.get_one(key))
    assert sb_helpers.seqs2hash(tester) == fwd_hash
    tester = Sb.order_features_alphabetically(sb_resources.get_one(key), reverse=True)
    assert sb_helpers.seqs2hash(tester) == rev_hash


# ######################  '-ofp', '--order_features_by_position' ###################### #
hashes = [('d f', 'b831e901d8b6b1ba52bad797bad92d14', 'b831e901d8b6b1ba52bad797bad92d14'),
          ('d g', '2e02a8e079267bd9add3c39f759b252c', '4345a14fe27570b3c837c30a8cb55ea9'),
          ('d n', 'cb1169c2dd357771a97a02ae2160935d', 'cb1169c2dd357771a97a02ae2160935d'),
          ('d py', '503e23720beea201f8fadf5dabda75e4', '503e23720beea201f8fadf5dabda75e4'),
          ('d pr', '52c23bd793c9761b7c0f897d3d757c12', '52c23bd793c9761b7c0f897d3d757c12'),
          ('d s', '228e36a30e8433e4ee2cd78c3290fa6b', '228e36a30e8433e4ee2cd78c3290fa6b'),
          ('p f', '14227e77440e75dd3fbec477f6fd8bdc', '14227e77440e75dd3fbec477f6fd8bdc'),
          ('p g', '7a8e25892dada7eb45e48852cbb6b63d', '9e7c2571db1386bba5983365ae235e1b'),
          ('p n', '17ff1b919cac899c5f918ce8d71904f6', '17ff1b919cac899c5f918ce8d71904f6'),
          ('p py', '968ed9fa772e65750f201000d7da670f', '968ed9fa772e65750f201000d7da670f'),
          ('p pr', 'ce423d5b99d5917fbef6f3b47df40513', 'ce423d5b99d5917fbef6f3b47df40513'),
          ('p s', 'c0dce60745515b31a27de1f919083fe9', 'c0dce60745515b31a27de1f919083fe9')]


@pytest.mark.parametrize("key,fwd_hash,rev_hash", hashes)
def test_order_features_position(key, fwd_hash, rev_hash, sb_resources, sb_helpers):
    tester = Sb.order_features_by_position(sb_resources.get_one(key))
    assert sb_helpers.seqs2hash(tester) == fwd_hash
    tester = Sb.order_features_by_position(sb_resources.get_one(key), reverse=True)
    assert sb_helpers.seqs2hash(tester) == rev_hash


# ######################  '-oi', '--order_ids' ###################### #
hashes = [('d f', 'e9090efdd362d527a115049dfced42cd', 'fcb016fff87d26822fa518d62e355c65'),
          ('d g', 'c0d656543aa5d20a266cffa790c035ce', '2507c667a304fdc003bc68255e094d7b'),
          ('d n', '132757da01b3caf174d024efdb2c3acd', '286bac7a213997924203622c3357457c'),
          ('d py', '3c49bdc1b0fe4e1d6bfc148eb0293e21', 'd6e79a5faeaff396aa7eab0b460c3eb9'),
          ('d pr', '684a99151e8cbf1e5eb96e44b875ba08', 'a5fcbd95e837b4807408124c2396db6e'),
          ('d s', 'c06985b566277a29e598dea0dc41baef', '20b43a724670a1151b1f0418478046ef')]


@pytest.mark.parametrize("key,fwd_hash,rev_hash", hashes)
def test_order_ids1(key, fwd_hash, rev_hash, sb_resources, sb_helpers):
    tester = Sb.order_ids(sb_resources.get_one(key))
    assert sb_helpers.seqs2hash(tester) == fwd_hash
    tester = Sb.order_ids(sb_resources.get_one(key), reverse=True)
    assert sb_helpers.seqs2hash(tester) == rev_hash


def test_order_ids2(sb_resources, sb_helpers):
    seqbuddy = sb_resources.get_one("p n")
    Sb.rename(seqbuddy, "Mle-Panxα4", "Mle004-Panxα4")
    Sb.rename(seqbuddy, "Mle-Panxα5", "Mle05-Panxα5")
    Sb.rename(seqbuddy, "Mle-Panxα9", "aMle-PanxαBlahh")
    Sb.order_ids(seqbuddy)
    assert sb_helpers.seqs2hash(seqbuddy) == "5c1316e18205432b044101e720646cd5"


# ######################  '-oir', '--order_ids_randomly' ###################### #
def test_order_ids_randomly(sb_resources, sb_helpers):
    for seqbuddy in sb_resources.get_list(""):
        tester = Sb.order_ids_randomly(Sb.make_copy(seqbuddy))
        assert sb_helpers.seqs2hash(seqbuddy) != sb_helpers.seqs2hash(tester)
        assert sb_helpers.seqs2hash(Sb.order_ids(tester)) == sb_helpers.seqs2hash(tester)


def test_order_ids_randomly2(sb_resources, sb_helpers):
    tester = sb_resources.get_one("d f")
    for _ in range(15):  # This will fail to repeat the while loop only ~5% of the time
        Sb.pull_recs(tester, "α[789]")
        assert sb_helpers.seqs2hash(tester) != sb_helpers.seqs2hash(Sb.order_ids_randomly(tester))

    Sb.pull_recs(tester, "α[89]")
    assert sb_helpers.seqs2hash(tester) != sb_helpers.seqs2hash(Sb.order_ids_randomly(tester))

    Sb.pull_recs(tester, "α[9]")
    assert sb_helpers.seqs2hash(tester) == sb_helpers.seqs2hash(Sb.order_ids_randomly(tester))

    tester = Sb.SeqBuddy(tester.records * 3)
    assert sb_helpers.seqs2hash(tester) == sb_helpers.seqs2hash(Sb.order_ids_randomly(tester))


# #####################  '-psc', '--prosite_scan' ###################### ##
def test_prosite_scan_init(sb_resources):
    seqbuddy = sb_resources.get_one("d f")
    ps_scan = Sb.PrositeScan(seqbuddy)
    assert hash(ps_scan.seqbuddy) == hash(seqbuddy)
    assert ps_scan.common_match
    assert not ps_scan.quiet
    assert ps_scan.base_url == 'http://www.ebi.ac.uk/Tools/services/rest/ps_scan'
    assert ps_scan.check_interval == 10
    assert len(ps_scan.http_headers) == 1
    assert 'User-Agent' in ps_scan.http_headers
    for key in ['data_dir', 'diagnostics', 'email', 'user_hash']:
        assert key in ps_scan.user_deets


def test_prosite_scan_rest_request(sb_resources, monkeypatch):
    def mock_urlopen(req, request_data=None):
        tmp_file = br.TempFile(byte_mode=True)
        file_text = "Hello world\n%s\n%s" % (req.full_url, request_data)
        tmp_file.write(file_text.encode())
        return tmp_file.get_handle("r")

    monkeypatch.setattr(urllib.request, "urlopen", mock_urlopen)
    seqbuddy = sb_resources.get_one("d f")
    ps_scan = Sb.PrositeScan(seqbuddy)
    assert ps_scan._rest_request("http://www.foo.bar", {"test1": 1}) == "Hello world\nhttp://www.foo.bar\n{'test1': 1}"
    assert ps_scan._rest_request("http://www.foo.bar") == "Hello world\nhttp://www.foo.bar\nNone"


def test_prosite_scan_mc_run_prosite(sb_resources, sb_helpers, monkeypatch):
    def status():
        for next_status in ["RUNNING", "PENDING", "FINISHED"]:
            yield next_status
    status_obj = status()

    def mock_rest_request(self, url, *args):
        if "run" in url:
            file_text = "job_1"
        elif "status" in url:
            file_text = next(status_obj)
        elif "result" in url:
            file_text = """\
>EMBOSS_001 : PS00001 ASN_GLYCOSYLATION N-glycosylation site.
    329 - 332  NNTA
>EMBOSS_001 : PS00004 CAMP_PHOSPHO_SITE cAMP- and cGMP-dependent protein kinase phosphorylation site.
    111 - 114  RRgS
    137 - 140  KKmT
>EMBOSS_001 : PS00005 PKC_PHOSPHO_SITE Protein kinase C phosphorylation site.
        4 - 6  SeK
      56 - 58  TvR
>EMBOSS_001 : PS00006 CK2_PHOSPHO_SITE Casein kinase II phosphorylation site.
    197 - 200  SigD
    241 - 244  SgiE
>EMBOSS_001 : PS00008 MYRISTYL N-myristoylation site.
      62 - 67  GSviSC
      74 - 79  GStfAE
>EMBOSS_001 : PS51013 PANNEXIN Pannexin family profile.
     28 - 353  WGITIDDGWDQLNRSFMFGLLVVMGTTVTVRQYTGSVISCDGFKKFGS---TFAEDYCWT L=0
 QGQYTVLEGYDQP-------NQNIPCPVPRPPSRRGSTLNTMSQTQGFLHNPV--ESDQE
 LKKMTDKAA------TWLFYKFDLYMSEQSLLASLTNKHG-------LGLSVVFVKILYA
 AVSFGCFLLTADMFSiGDFKTYGSEWINKLKLeDNLATEEKDKLFPKMVACEV-KRWGAS
 GIEEEQGMCVLAPNVINQYLFLILWFCLVFVMFCNIVSIFASLIKLLFTYG-----SYRR
 LLSTAFLRDDSAIKHMYFNVGSSGRLILHVLANNTAPRVFEDILLTLAPKLIQRKLR
"""
        else:
            raise RuntimeError("This shouldn't ever happen", self, args)

        return file_text

    monkeypatch.setattr(Sb.PrositeScan, "_rest_request", mock_rest_request)
    monkeypatch.setattr(Sb.time, "sleep", lambda _: True)
    out_file = br.TempFile()
    seqbuddy = sb_resources.get_one("d f")
    Sb.pull_recs(seqbuddy, "Mle-Panxα10B")
    ps_scan = Sb.PrositeScan(seqbuddy)
    ps_scan._mc_run_prosite(seqbuddy.records[0], [out_file.path, Sb.Lock()])
    assert sb_helpers.string2hash(out_file.read()) == "7ced43edaee481ac149d6ece152c4621"


def test_prosite_scan_run(sb_resources, sb_helpers, monkeypatch):
    def mock_mc_run_prosite(self, _rec, args):
        print(self)
        out_file_path, lock = args
        temp_seq = Sb.SeqBuddy([_rec], out_format="gb")
        Sb.annotate(temp_seq, "Foo", "1-100")
        with lock:
            with open(out_file_path, "a") as out_file:
                out_file.write("%s\n" % str(temp_seq))

    monkeypatch.setattr(Sb.PrositeScan, "_mc_run_prosite", mock_mc_run_prosite)
    seqbuddy = sb_resources.get_one("d g")
    Sb.delete_features(seqbuddy, "splice")
    ps_scan = Sb.PrositeScan(seqbuddy)
    seqbuddy = ps_scan.run()
    assert sb_helpers.seqs2hash(seqbuddy) == "bc477b683784a24524b72422e04ff949"

    seqbuddy = sb_resources.get_one("p g")
    Sb.delete_features(seqbuddy, "splice")
    ps_scan = Sb.PrositeScan(seqbuddy)
    seqbuddy = ps_scan.run()
    seqbuddy.write("temp.del")
    assert sb_helpers.seqs2hash(seqbuddy) == "e8cd292ada589ddde4747bd9f9ebfb17"


# #####################  '-prr', '--pull_random_recs' ###################### ##
def test_pull_random_recs(sb_resources):
    for seqbuddy in sb_resources.get_list(""):
        tester = Sb.pull_random_recs(Sb.make_copy(seqbuddy))
        orig_seqs = tester.to_dict()
        assert len(tester.records) == 1
        assert tester.records[0].id in orig_seqs


# #####################  '-pre', '--pull_record_ends' ###################### ##
def test_pull_record_ends(sb_resources, sb_helpers):
    tester = Sb.pull_record_ends(sb_resources.get_one("d g"), 10)
    assert sb_helpers.seqs2hash(tester) == 'd46867e4ca7a9f474c45473fc3495413'

    tester = Sb.pull_record_ends(sb_resources.get_one("d g"), 2000)
    assert sb_helpers.seqs2hash(tester) == '908744b00d9f3392a64b4b18f0db9fee'

    tester = Sb.pull_record_ends(sb_resources.get_one("d g"), -10)
    assert sb_helpers.seqs2hash(tester) == 'd7970570d65872993df8a3e1d80f9ff5'

    tester = Sb.pull_record_ends(sb_resources.get_one("d g"), -2000)
    assert sb_helpers.seqs2hash(tester) == '908744b00d9f3392a64b4b18f0db9fee'

    with pytest.raises(ValueError):
        Sb.pull_record_ends(sb_resources.get_one("d f"), 'foo')

# ######################  '-pr', '--pull_records' ###################### #
hashes = [('d f', '5b4154c2662b66d18776cdff5af89fc0'), ('d g', 'e196fdc5765ba2c47f97807bafb6768c'),
          ('d n', 'bc7dbc612bc8139eba58bf896b7eaf2f'), ('d py', '7bb4aac2bf50381ef1d27d82b7dd5a53'),
          ('d pr', '72328eddc751fd79406bb911dafa57a2'), ('d s', 'b006b40ff17ba739929448ae2f9133a6')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_pull_recs(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.pull_recs(sb_resources.get_one(key), 'α2')
    assert sb_helpers.seqs2hash(tester) == next_hash

# ######################  '-pr', '--pull_records_with_feature' ###################### #
hashes = [('p g', '83d15851d489e89761c8faa31e5263f2'), ('d g', '36757409966ede91ab19deb56045d584')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_pull_records_with_feature(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.pull_recs_with_feature(sb_resources.get_one(key), 'splice_acceptor')
    assert sb_helpers.seqs2hash(tester) == next_hash


# #####################  '-prg', '--purge' ###################### ##
def test_purge(sb_resources, sb_helpers):
    tester = sb_resources.get_one("p f")
    Sb.purge(tester, 200)
    assert sb_helpers.seqs2hash(tester) == 'b21b2e2f0ca1fcd7b25efbbe9c08858c'


# ######################  '-ri', '--rename_ids' ###################### #
hashes = [('d f', '8b4a9e3d3bb58cf8530ee18b9df67ff1'), ('d g', '78c73f97117bd937fd5cf52f4bd6c26e'),
          ('d n', '243024bfd2f686e6a6e0ef65aa963494'), ('d py', '98bb9b57f97555d863054ddb526055b4'),
          ('d pr', '2443c47a712f19099e94fc015dc980a9'), ('d s', '65196fd4f2a4e339e1545f6ed2a6acc3')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_rename_ids(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.rename(sb_resources.get_one(key), 'Panx', 'Test', 0)
    assert sb_helpers.seqs2hash(tester) == next_hash


def test_rename_ids2(sb_resources, sb_helpers):
    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "Panxα([1-6])", "testing\\1")
    assert sb_helpers.seqs2hash(tester) == "f0d7ec055b3b2d9ec11a86634b32e9ef"

    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "Panxα([1-6])", "testing\\1", -1)
    assert sb_helpers.seqs2hash(tester) == "f0d7ec055b3b2d9ec11a86634b32e9ef"

    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "[a-z]", "?", 2)
    assert sb_helpers.seqs2hash(tester) == "08aaee2e7b9997b512c7d1b2fe748d40"

    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "[a-z]", "?", -2)
    assert sb_helpers.seqs2hash(tester) == "9b3946afde20c991099463d099be22e0"

    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "[A-Z]", "?", -20)
    assert sb_helpers.seqs2hash(tester) == "451993d7e816881e2700697263b1d8fa"

    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "[a-z]", "?", 2, store_old_id=True)
    assert sb_helpers.seqs2hash(tester) == "959fe04c0366c9a143052a02f090707e"

    with pytest.raises(AttributeError) as e:
        Sb.rename(tester, "[a-z]", "\\1?", -2)
    assert "There are more replacement match" in str(e)


# ##################### '-rs', 'replace_subseq' ###################### ##
def test_replace_subsequence(sb_resources, sb_helpers):
    tester = sb_resources.get_one("d f")
    Sb.replace_subsequence(tester, "atg(.{5}).{3}", "FOO\\1BAR")
    assert sb_helpers.seqs2hash(tester) == "f12707c2b0ef866f0039bac96abb29e0"


# ######################  '-rc', '--reverse_complement' ###################### #
hashes = [('d f', 'e77be24b8a7067ed54f06e0db893ce27'), ('d g', '47941614adfcc5bd107f71abef8b3e00'),
          ('d n', 'f549c8dc076f6b3b4cf5a1bc47bf269d'), ('d py', '0dd20827c011a0a7a0e78881b38ae06a'),
          ('d pr', '0b954f4a263bf38ddeac61ab54f77dc2'), ('d s', '0d6b7deda824b4fc42b65cb87e1d4d14')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_reverse_complement(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.reverse_complement(sb_resources.get_one(key))
    assert sb_helpers.seqs2hash(tester) == next_hash


def test_reverse_complement_pep_exception(sb_resources):  # Asserts a TypeError will be thrown if user inputs protein
    tester = sb_resources.get_one('p f')
    with pytest.raises(TypeError) as e:
        Sb.reverse_complement(tester)
    assert str(e.value) == "SeqBuddy object is protein. Nucleic acid sequences required."

    prot_record = tester.records[0]
    tester = sb_resources.get_one('d f')
    tester.records.append(prot_record)
    with pytest.raises(TypeError) as e:
        Sb.reverse_complement(tester)
    assert str(e.value) == "Record 'Mle-Panxα12' is protein. Nucleic acid sequences required."

# ######################  '-sfr', '--select_frame' ###################### #
hashes = [('d f', 1, "b831e901d8b6b1ba52bad797bad92d14"), ('d f', 2, "2de033b2bf2327f2795fe425db0bd78f"),
          ('d f', 3, "1c29898d4964e0d1b03207d7e67e1958"), ('d g', 1, "908744b00d9f3392a64b4b18f0db9fee"),
          ('d g', 2, "08fe54a87249f5fb9ba22ff6d0053787"), ('d g', 3, "cfe2d405487d69dceb2a11dd44ceec59"),
          ('d n', 1, "cb1169c2dd357771a97a02ae2160935d"), ('d n', 2, "87d784f197b55f812d2fc82774da43d1"),
          ('d n', 3, "5d6d2f337ecdc6f9a85e981c975f3e08")]


@pytest.mark.parametrize("key,shift,next_hash", hashes)
def test_select_frame(key, shift, next_hash, sb_resources, sb_helpers):
    tester = Sb.select_frame(sb_resources.get_one(key), shift)
    assert sb_helpers.seqs2hash(tester) == next_hash


def test_select_frame_edges(sb_resources, sb_helpers):
    tester = Sb.select_frame(sb_resources.get_one("d f"), 2)
    temp_file = br.TempFile()
    tester.write(temp_file.path)
    tester = Sb.select_frame(Sb.SeqBuddy(temp_file.path), 1)
    assert sb_helpers.seqs2hash(tester) == "b831e901d8b6b1ba52bad797bad92d14"

    tester = Sb.select_frame(sb_resources.get_one("d g"), 2)
    tester = Sb.select_frame(tester, 1)
    assert sb_helpers.seqs2hash(tester) == "908744b00d9f3392a64b4b18f0db9fee"

    with pytest.raises(TypeError) as e:  # If protein is input
        Sb.select_frame(sb_resources.get_one("p f"), 2)
    assert "Select frame requires nucleic acid, not protein." in str(e.value)


# ##################### '-ss', 'shuffle_seqs' ###################### ##
def test_shuffle_seqs(sb_resources, sb_helpers):
    for seqbuddy in sb_resources.get_list(""):
        tester1 = Sb.make_copy(seqbuddy)
        tester2 = Sb.shuffle_seqs(seqbuddy)
        assert sb_helpers.seqs2hash(tester1) != sb_helpers.seqs2hash(tester2)

        for indx, record in enumerate(tester1.records):
            assert sorted(record.seq) == sorted(tester2.records[indx].seq)


# #####################  make_groups' ###################### ##
def test_make_groups(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources["cnidaria_pep"])
    sb_list = Sb.make_groups(tester)
    assert len(sb_list) == 20
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier == seqbuddy.records[0].id

    sb_list = Sb.make_groups(tester, split_patterns=["u", "h"])
    assert len(sb_list) == 4
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Unknown", "Hv", "C", "Pp"]

    sb_list = Sb.make_groups(tester, num_chars=1)
    assert len(sb_list) == 5
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["A", "H", "C", "P", "N"]

    sb_list = Sb.make_groups(tester, split_patterns=["l"], num_chars=3)
    assert len(sb_list) == 3
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Unknown", "Ae", "C"]

    sb_list = Sb.make_groups(tester, regex="([ACH]).*([βγ])")
    assert len(sb_list) == 5
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Unknown", "Aβ", "Cβ", "Hβ", "Cγ"]

    sb_list = Sb.make_groups(tester, regex="Ate")
    assert len(sb_list) == 2
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Ate", "Unknown"]

    sb_list = Sb.make_groups(tester, regex="Panx.(G*)")
    assert len(sb_list) == 2
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["G", "Unknown"]


# ######################  '-tr6', '--translate6frames' ###################### #
def test_translate6frames(sb_resources, sb_helpers):
    tester = Sb.translate6frames(sb_resources.get_one("d f"))
    assert sb_helpers.seqs2hash(tester) == '95cf24202007399e6ccd6e6f33ae012e'

    tester = Sb.translate6frames(sb_resources.get_one("d g"))
    assert sb_helpers.seqs2hash(tester) == '0b5daa810e1589c3973e1436c40baf08'


def test_translate6frames_pep_exception(sb_resources):
    with pytest.raises(TypeError):
        Sb.translate6frames(sb_resources.get_one("p f"))

# ######################  '-tr', '--translate' ###################### #
hashes = [('d f', '06893e14839dc0448e6f522c1b8f8957'), ('d g', 'e8840e22096e933ce10dbd91036f3fa5'),
          ('d n', 'f3339e0193c10427f017dd8f6bd81d7e'), ('r f', '06893e14839dc0448e6f522c1b8f8957'),
          ('r n', 'f3339e0193c10427f017dd8f6bd81d7e')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_translate(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.translate_cds(sb_resources.get_one(key))
    assert sb_helpers.seqs2hash(tester) == next_hash


def test_translate_ambig(sb_odd_resources, sb_helpers):
    tester = Sb.SeqBuddy(sb_odd_resources['ambiguous_dna'])
    tester = Sb.translate_cds(tester)
    assert sb_helpers.seqs2hash(tester) == '648ccc7c3400882be5bf6e8d9781f74e'

    tester = Sb.SeqBuddy(sb_odd_resources['ambiguous_rna'])
    tester = Sb.translate_cds(tester)
    assert sb_helpers.seqs2hash(tester) == '648ccc7c3400882be5bf6e8d9781f74e'


def test_translate_edges_and_exceptions(capsys, sb_resources, sb_helpers):
    with pytest.raises(TypeError):
        Sb.translate_cds(sb_resources.get_one("p f"))

    tester = Sb.SeqBuddy("ATTCGTTAACGCTAGCGTCG", in_format="raw")
    tester = Sb.translate_cds(tester)
    assert str(tester) == ">raw_input\nIR*R*R\n"
    out, err = capsys.readouterr()
    assert err == "Warning: size mismatch between aa and nucl seqs for raw_input --> 20, 6\n"

    tester = Sb.select_frame(sb_resources.get_one("d g"), 3)
    tester = Sb.translate_cds(tester)
    assert sb_helpers.seqs2hash(tester) == "68ca15f5ac737e4a4ca65a67ad2dc897"
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(err) == "9e2a0b4b03f54c209d3a9111792762df"


# ######################  '-tmd', '--transmembrane_domains' ###################### #
def test_transmembrane_domains_pep(sb_resources, sb_helpers, monkeypatch, capsys):
    def mock_runtimeerror(*args, **kwargs):
        raise RuntimeError("%s %s" % (args, kwargs))

    def mock_contenttooshorterror(*args, **kwargs):
        raise urllib.error.ContentTooShortError("%s %s" % (args, kwargs), "content")

    def mock_httperror(*args, **kwargs):
        raise urllib.error.HTTPError(url="http://fake.come", code=503, msg=args, hdrs=kwargs, fp="Bar")

    class MockSudsClient(object):
        def __init__(self):
            self.service = MockServiceSelector()

    class MockServiceSelector(object):
        def __init__(self):
            self.job_id_generator = self.job_id_gen()
            self.current_job_id = next(self.job_id_generator)
            self.result_url = "www.something.com"
            self.numseq_str = "doesn't matter"
            self.errinfo = "Did I fail?"
            self.warninfo = "Also doesn't matter"
            self.job_check = self.status()
            self.job_statuses = ["Queued", "Running", "Finished"]

        def status(self):
            while True:
                for next_status in self.job_statuses:
                    yield next_status

        @staticmethod
        def job_id_gen():
            for job_id in ["rst_MFhyxO", "rst_lE27A5"]:
                yield job_id

        def submitjob(self, *args):
            return [[self.current_job_id, self.result_url, self.numseq_str, self.errinfo, self.warninfo]]

        def checkjob(self, jobid):
            #raise ValueError
            return [[next(self.job_check), self.result_url, self.errinfo]]

    def mock_urlretrieve(result_url, filename, reporthook):
        job_id = filename.split("/")[-1].split(".")[0]
        if os.path.isfile("%s/%s.hashmap" % (work_dir.path, job_id)):
            os.remove("%s/%s.hashmap" % (work_dir.path, job_id))
        shutil.copy("%s/topcons/%s.zip" % (sb_resources.res_path, job_id), "%s/" % work_dir.path)
        reporthook(2, 10, 100)
        return

    def mock_hash_ids(seqbuddy):
        hashmap = OrderedDict()
        with open("%s/topcons/%s.hashmap" % (sb_resources.res_path, suds_client.service.current_job_id), "r") as ifile:
            for line in ifile:
                if line:
                    line = line.strip().split("\t")
                    hashmap[line[0]] = line[1]
        seqbuddy.hash_map = hashmap
        for rec, hash_id in zip(seqbuddy.records, list(hashmap.items())):
            rec.id = hash_id[0]
            rec.name = hash_id[0]
        return seqbuddy

    work_dir = br.TempDir()
    keep_dir = br.TempDir()
    suds_client = MockSudsClient()

    monkeypatch.setattr(time, "sleep", lambda _: True)
    monkeypatch.setattr(suds.client, "Client", lambda *_, **__: suds_client)
    monkeypatch.setattr(urllib.request, "urlretrieve", mock_urlretrieve)
    monkeypatch.setattr(br, "TempDir", lambda: work_dir)
    monkeypatch.setattr(Sb, "hash_ids", mock_hash_ids)

    tester = sb_resources.get_one("d g")
    Sb.pull_recs(tester, "α[56]")
    Sb.delete_features(tester, "splice|TMD")
    tester = Sb.transmembrane_domains(tester)
    assert sb_helpers.seqs2hash(tester) == "443462d4a7d7ed3121378fca55491d5c"

    suds_client.service.current_job_id = next(suds_client.service.job_id_generator)
    tester = sb_resources.get_one("p g")
    Sb.pull_recs(tester, "α[56]")
    Sb.delete_features(tester, "splice|TMD")
    tester = Sb.transmembrane_domains(tester)
    assert sb_helpers.seqs2hash(tester) == "eb31602e292e5a056b956f13dbb0d590"

    tester = sb_resources.get_one("p g")
    Sb.pull_recs(tester, "α[56]")
    Sb.delete_features(tester, "splice|TMD")
    tester = Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    assert sb_helpers.seqs2hash(tester) == "eb31602e292e5a056b956f13dbb0d590"

    tester = sb_resources.get_one("p g")
    Sb.pull_recs(tester, "α[56]")
    Sb.delete_features(tester, "splice|TMD")
    tester = Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"], keep_temp=keep_dir.path)
    _root, dirs, files = next(br.walklevel(keep_dir.path))
    assert dirs == ['rst_lE27A5', 'rst_MFhyxO', 'topcons']
    assert files == ['seqs.tmp', work_dir.path.split("/")[-1]]

    with pytest.raises(FileNotFoundError) as err:
        Sb.transmembrane_domains(tester, job_ids=["rst_BLAHHH!!"])
    assert "SeqBuddy does not have the necessary hash-map to process job id 'rst_BLAHHH!!'." in str(err)

    suds_client.service.job_statuses = ["Failed"]
    with pytest.raises(ConnectionError) as err:
        Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    assert "Job failed..." in str(err)

    suds_client.service.job_statuses = ["None"]
    with pytest.raises(ConnectionError) as err:
        Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    assert "The job seems to have been lost by the server." in str(err)

    suds_client.service.job_statuses = ["Finished"]
    monkeypatch.setattr(urllib.request, "urlretrieve", mock_runtimeerror)
    Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    out, err = capsys.readouterr()
    assert "Error: Failed to download TOPCONS job rst_lE27A5 after 5 attempts." in err

    monkeypatch.setattr(urllib.request, "urlretrieve", mock_contenttooshorterror)
    Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    out, err = capsys.readouterr()
    assert "Error: Failed to download TOPCONS job rst_lE27A5 after 5 attempts." in err

    monkeypatch.setattr(urllib.request, "urlretrieve", mock_httperror)
    Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    out, err = capsys.readouterr()
    assert "Error: Failed to download TOPCONS job rst_lE27A5 after 5 attempts." in err

