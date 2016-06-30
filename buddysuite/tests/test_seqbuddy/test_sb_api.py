#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" tests basic functionality of AlignBuddy class """
import pytest
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from unittest import mock
import os

try:
    from buddysuite import SeqBuddy as Sb
    from buddysuite import buddy_resources as br
except ImportError:
    import SeqBuddy as Sb
    import buddy_resources as br


TEMPDIR = br.TempDir()


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
# ToDo: mock output data
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


hashes = [('d g', '7421c27be7b41aeedea73ff41869ac47'), ('d n', '2e46edb78e60a832a473397ebec3d187'),
          ('p g', '9cb5443d90e64693ad1bd74f29169ac5'), ('p py', 'e902a4faf44ebd28a43ca8103df7b828')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_concat_seqs_clean(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.concat_seqs(sb_resources.get_one(key), clean=True)
    tester.out_format = "gb"
    assert sb_helpers.seqs_to_hash(tester) == next_hash


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
    assert sb_helpers.seqs_to_hash(tester) == '0638bc6546eebd9d50f771367d6d7855'


def test_degenerate_sequence_with_different_codon_tables(sb_resources, sb_helpers):
    for indx, next_hash in enumerate(['0638bc6546eebd9d50f771367d6d7855', '72373f8356051e2c6b67642451379054',
                                      '9172ad5947c0961b54dc5adbd03d4249', 'b45ac94ee6a98e495e115bfeb5bd9bcd',
                                      '76c45b4de8f7527b4139446b4551712b', 'baa5b48938cc5cae953c9083a5b21b12',
                                      '0ca67c4740fefbc7a20d806715c3ca12', 'd43ad8f328ff1d30eb1fb7bcd667a345',
                                      'd9d0f5cd8f0c25a0042527cc1cea802e', '4b9790f3f4eeeae1a9667b62b93bc961',
                                      '7ec4365c3571813d63cee4b70ba5dcf5']):
        tester = Sb.degenerate_sequence(sb_resources.get_one("d f"), table=(indx + 1))
        assert sb_helpers.seqs_to_hash(tester) == next_hash


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
    assert sb_helpers.seqs_to_hash(tester) == "ff52e05971aeafd24c73a3b543901e4b"


# ######################  '-df', '--delete_features' ###################### #
def test_delete_features(sb_resources, sb_helpers):
    tester = sb_resources.get_one("d g")
    tester = Sb.delete_features(tester, 'donor')
    assert sb_helpers.seqs_to_hash(tester) == 'f84df6a77063c7def13babfaa0555bbf'


# ######################  '-dl', '--delete_large' ###################### #
def test_delete_large(sb_resources, sb_helpers):
    tester = sb_resources.get_one("d f")
    tester = Sb.delete_large(tester, 1285)
    assert sb_helpers.seqs_to_hash(tester) == '25859dc69d46651a1e04a70c07741b35'


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
    assert sb_helpers.seqs_to_hash(tester) == next_hash


# ######################  '-dr', '--delete_records' ###################### #
hashes = [('d f', '54bdb42423b1d331acea18218101e5fc'), ('d g', 'e2c03f1fa21fd27b2ff55f7f721a1a99'),
          ('d n', '6bc8a9409b1ef38e4f6f12121368883e'), ('d py', 'bda7be10061b0dcaeb66bebe3d736fee'),
          ('d pr', '6e2fce2239e2669b23f290049f87fbc4'), ('d s', '4c97144c5337f8a40c4fd494e622bf0d')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_delete_records(key, next_hash, sb_resources, sb_helpers):
    tester = Sb.delete_records(sb_resources.get_one(key), 'α2')
    assert sb_helpers.seqs_to_hash(tester) == next_hash


def test_delete_records2(sb_resources, sb_helpers):
    tester = Sb.delete_records(sb_resources.get_one("d f"), ['α1', 'α2'])
    assert sb_helpers.seqs_to_hash(tester) == "eca4f181dae3d7998464ff71e277128f"

    tester = Sb.delete_records(sb_resources.get_one("d f"), 'α1|α2')
    assert sb_helpers.seqs_to_hash(tester) == "eca4f181dae3d7998464ff71e277128f"

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
    assert sb_helpers.seqs_to_hash(tester) == '196adf08d4993c51050289e5167dacdf'


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
    assert sb_helpers.seqs_to_hash(tester) == d2r_hash
    tester = Sb.rna2dna(tester)
    assert sb_helpers.seqs_to_hash(tester) == r2d_hash


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


"""
# ######################  '-er', '--extract_regions' ###################### #
hashes = ["8c2fac57aedf6b0dab3d0f5bcf88e99f", "25ad9670e8a6bac7962ab46fd79251e5", "4063ab66ced2fafb080ceba88965d2bb",
          "33e6347792aead3c454bac0e05a292c6", "9a5c491aa293c6cedd48c4c249d55aff", "cd8d857feba9b6e459b8a9d56f11b7f5",
          "2586d1e3fc283e6f5876251c1c57efce", "a9c22659967916dcdae499c06d1aaafb", "6a27222d8f60ee8496cbe0c41648a116",
          "c9a1dd913190f95bba5eca6a89685c75", "6f579144a43dace285356ce6eb326d3b", "727099e0abb89482760eeb20f7edd0cd"]
hashes = [(Sb.make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_extract_regions(seqbuddy, next_hash):
    tester = Sb.extract_regions(Sb.make_copy(seqbuddy), "50:300")
    assert seqs_to_hash(tester) == next_hash

    tester = Sb.extract_regions(Sb.make_copy(seqbuddy), "300:50")
    assert seqs_to_hash(tester) == next_hash


# #####################  '-fcpg', '--find_CpG' ###################### ##
def test_find_cpg():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.gb"))
    tester = Sb.find_cpg(tester)
    assert seqs_to_hash(tester) == "9499f524da0c35a60502031e94864928"


# #####################  '-fp', '--find_pattern' ###################### ##
def test_find_pattern():
    tester = Sb.find_pattern(sb_resources.get_one("d g"), "ATGGT")
    assert seqs_to_hash(tester) == "ca129f98c6c719d50f0cf43eaf6dc90a"
    tester.out_format = "fasta"
    assert seqs_to_hash(tester) == "6f23f80b52ffb736bbecc9f4c72d8fab"
    tester = Sb.find_pattern(sb_resources.get_one("d g"), "ATg{2}T")
    assert seqs_to_hash(tester) == "9ec8561c264bff6f7166855d60457df1"
    tester = Sb.find_pattern(sb_resources.get_one("d g"), "ATg{2}T", "tga.{1,6}tg")
    assert seqs_to_hash(tester) == "ec43ce98c9ae577614403933b2c5f37a"
    tester = Sb.find_pattern(sb_resources.get_one("d g"), "ATg{2}T", "tga.{1,6}tg", include_feature=False)
    assert seqs_to_hash(tester) == "2e02a8e079267bd9add3c39f759b252c"
    tester = Sb.find_pattern(sb_resources.get_one("p g"), "[bz]{2}x{50,100}[bz]{2}", ambig=True)
    assert seqs_to_hash(tester) == "339ff26803a2d12267d873458d40bca2"
    tester = Sb.find_pattern(sb_resources.get_one("d g"), "ATGGN{6}", ambig=True)
    assert seqs_to_hash(tester) == "ac9adb42fbfa9cf22f033e9a02130985"
    tester = Sb.find_pattern(sb_resources.get_one("r f"), "AUGGN{6}", ambig=True)
    assert seqs_to_hash(tester) == "b7abcb4334232e38dfbac9f46234501a"


# #####################  '-frp', '--find_repeats' ###################### ##
def test_find_repeats():
    tester = Sb.SeqBuddy(resource("Duplicate_seqs.fa"))
    tester.unique_seqs, tester.repeat_ids, tester.repeat_seqs = {}, {}, {}
    Sb.find_repeats(tester)
    assert 'Seq1' in tester.unique_seqs
    assert 'Seq12' in tester.repeat_ids

    for key in tester.repeat_seqs:
        assert 'Seq12' in tester.repeat_seqs[key] or 'Seq10A' in tester.repeat_seqs[key]


# ######################  '-frs', '--find_restriction_sites' ###################### #
def test_restriction_sites(capsys):
    # No arguments passed in = commercial REs and any number of cut sites
    tester = Sb.find_restriction_sites(Sb.make_copy(sb_objects[1]))
    assert seqs_to_hash(tester) == 'fce4c8bee040d5ea6fa4bf9985f7310f'
    assert md5(str(tester.restriction_sites).encode()).hexdigest() == "741ca6ca9204a067dce7398f15c6e350"

    # All enzymes
    tester = Sb.find_restriction_sites(Sb.make_copy(sb_objects[1]), enzyme_group=["all"])
    assert seqs_to_hash(tester) == '254e854862142c54a53079ae224b4180'
    assert md5(str(tester.restriction_sites).encode()).hexdigest() == "0a5012b1a3d83f719cdd8cecf48258f4"

    # Specify a few REs and limit the number of cuts
    tester = Sb.find_restriction_sites(Sb.make_copy(sb_objects[1]), min_cuts=2, max_cuts=4,
                                       enzyme_group=["EcoRI", "KspI", "TasI", "Bme1390I", "FooBR"])
    out, err = capsys.readouterr()
    assert seqs_to_hash(tester) == 'c42b3bf0367557383000b897432fed2d'
    assert md5(str(tester.restriction_sites).encode()).hexdigest() == "0d2e5fdba6fed434495481397a91e56a"
    assert "Warning: FooBR not a known enzyme" in err

    with pytest.raises(TypeError) as e:
        Sb.find_restriction_sites(sb_objects[7])
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
def test_hash_seq_ids():
    tester = Sb.SeqBuddy(Sb.make_copy(sb_objects[0]))
    Sb.hash_ids(tester)
    assert len(tester.records[0].id) == 10

    tester = Sb.SeqBuddy(Sb.make_copy(sb_objects[0]))
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


def test_hash_seq_ids_errors():
    tester = Sb.SeqBuddy(Sb.make_copy(sb_objects[0]))
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
def test_insert_seqs_start():
    tester = Sb.make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA')) == 'f65fee08b892af5ef93caa1bf3cb3980'

    tester = Sb.make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', -9000)) == 'f65fee08b892af5ef93caa1bf3cb3980'

    tester = Sb.make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', -1)) == '792397e2e32e95b56ddc15b8b2310ec0'

    tester = Sb.make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', 9000)) == '792397e2e32e95b56ddc15b8b2310ec0'

    tester = Sb.make_copy(sb_objects[0])
    Sb.insert_sequence(tester, 'AACAGGTCGAGCA', 100, ["α[23]", "α5"])
    assert seqs_to_hash(tester) == 'edcd7934eb026ac3ea4b603ac85ca79f'

    tester = Sb.make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', -25)) == '29cab1e72ba95572c3aec469270071e9'


# ######################  '-ip', '--isoelectric_point' ###################### #
@pytest.mark.parametrize("seqbuddy", [Sb.make_copy(x) for x in [sb_objects[6], sb_objects[7],
                                                                sb_objects[8], sb_objects[10], sb_objects[11]]])
def test_isoelectric_point(seqbuddy):
    tester = Sb.isoelectric_point(Sb.clean_seq(seqbuddy))
    tester = tester.to_dict()
    assert tester["Mle-Panxα12"].features[-1].qualifiers["value"] == 6.0117797852
    if seqbuddy.out_format == "gb":
        assert seqs_to_hash(seqbuddy) == "8bc299e31f436d192bf8cf8b7af671a8"

    with pytest.raises(TypeError):
        Sb.isoelectric_point(Sb.make_copy(sb_objects[0]))


# ######################  '-lc', '--lowercase' and 'uc', '--uppercase'  ###################### #
uc_hashes = ["25073539df4a982b7f99c72dd280bb8f", "2e02a8e079267bd9add3c39f759b252c", "52e74a09c305d031fc5263d1751e265d",
             "cfe6cb9c80aebd353cf0378a6d284239", "6e5542f41d17ff33afb530b4d07408a3", "b82538a4630810c004dc8a4c2d5165ce",
             "c10d136c93f41db280933d5b3468f187", "7a8e25892dada7eb45e48852cbb6b63d", "8b6737fe33058121fd99d2deee2f9a76",
             "968ed9fa772e65750f201000d7da670f", "ce423d5b99d5917fbef6f3b47df40513", "f35cbc6e929c51481e4ec31e95671638"]

lc_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "2e02a8e079267bd9add3c39f759b252c", "cb1169c2dd357771a97a02ae2160935d",
             "503e23720beea201f8fadf5dabda75e4", "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
             "14227e77440e75dd3fbec477f6fd8bdc", "7a8e25892dada7eb45e48852cbb6b63d", "17ff1b919cac899c5f918ce8d71904f6",
             "aacda2f5d4077f23926400f74afa2f46", "e3dc2e0347f40fffec45d053f4f34c96", "c0dce60745515b31a27de1f919083fe9"]

hashes = [(Sb.make_copy(sb_objects[indx]), uc_hash, lc_hashes[indx]) for indx, uc_hash in enumerate(uc_hashes)]


@pytest.mark.parametrize("seqbuddy,uc_hash,lc_hash", hashes)
def test_cases(seqbuddy, uc_hash, lc_hash):  # NOTE: Biopython always writes genbank to spec in lower case
    tester = Sb.uppercase(seqbuddy)
    assert seqs_to_hash(tester) == uc_hash
    tester = Sb.lowercase(tester)
    assert seqs_to_hash(tester) == lc_hash


# ######################  '-mui', '--make_ids_unique' ###################### #
def test_make_ids_unique():
    tester = Sb.SeqBuddy(resource("Duplicate_seqs.fa"))
    Sb.make_ids_unique(tester)
    assert seqs_to_hash(tester) == "363c7ed14be59bcacede092b8f334a52"

    tester = Sb.SeqBuddy(resource("Duplicate_seqs.fa"))
    Sb.make_ids_unique(tester, sep="-", padding=4)
    assert seqs_to_hash(tester) == "0054df3003ba16287159147f3b85dc7b"


# ######################  '-fn2p', '--map_features_nucl2prot' ###################### #
# Map the genbank DNA file to all protein files, and the fasta DNA file to fasta protein
hashes = ["5216ef85afec36d5282578458a41169a", "a8f7c129cf57a746c20198bf0a6b9cf4", "bb0c9da494b5418fb87862dab2a66cfa",
          "3c0e3ec45abd774813a274fda1b4a5f2", "a7f6c4bb410f17cfc3e8966ccbe3e065", "1b8c44f4ace877b568c0915033980bed"]
prot_indx = [6, 7, 8, 9, 10, 11]
hashes = [(Sb.make_copy(sb_objects[1]), Sb.make_copy(sb_objects[prot_indx[indx]]), value)
          for indx, value in enumerate(hashes)]
hashes.append((Sb.make_copy(sb_objects[0]), Sb.make_copy(sb_objects[6]), "854566b485af0f277294bbfb15f7dd0a"))


@pytest.mark.parametrize("_dna,_prot,next_hash", hashes)
def test_map_features_nucl2prot(_dna, _prot, next_hash):
    tester = Sb.map_features_nucl2prot(_dna, _prot)
    tester.out_format = "gb"
    assert seqs_to_hash(tester) == next_hash


def test_map_features_nucl2prot_2(capsys):
    tester = Sb.make_copy(sb_objects[1])
    Sb.pull_record_ends(tester, 1300)
    tester = Sb.annotate(tester, "foo", [(20, 40), (50, 60)], pattern="α9")
    mapped = Sb.map_features_nucl2prot(Sb.make_copy(tester), Sb.make_copy(sb_objects[6]))
    mapped.out_format = "gb"
    assert seqs_to_hash(mapped) == "807025489aadf98d851501f49d463e4a"
    out, err = capsys.readouterr()
    assert string2hash(err) == "6c840e4acaaf4328672ca164f854000a"

    prot_tester = Sb.make_copy(sb_objects[7])
    prot_tester.records = sorted(prot_tester.records, key=lambda x: x.id)
    dna_tester = Sb.make_copy(sb_objects[1])
    dna_tester.records = sorted(dna_tester.records, key=lambda x: x.id)
    Sb.rename(dna_tester, "α4", "A4")

    Sb.map_features_nucl2prot(Sb.make_copy(dna_tester), prot_tester, mode="list")
    assert seqs_to_hash(prot_tester) == "9fdb606ea65d6c050540a94137ae6e0d"

    Sb.map_features_nucl2prot(Sb.make_copy(dna_tester), prot_tester, mode="key")
    assert seqs_to_hash(prot_tester) == "9fdb606ea65d6c050540a94137ae6e0d"
    out, err = capsys.readouterr()
    assert string2hash(err) == "6fd2b5f2a7a3995d3f49c4919c3358b0"

    mock_obj = type('test', (object,), {})()
    mock_obj.start = 2
    mock_obj.end = 6
    dna_tester.records[0].features[0].location = mock_obj
    with pytest.raises(TypeError) as e:
        Sb.map_features_nucl2prot(Sb.make_copy(dna_tester), Sb.make_copy(sb_objects[6]))
    assert "FeatureLocation or CompoundLocation object required." in str(e.value)

    tester.records[0].features[0].location = tester.records[1].features[0].location
    with pytest.raises(ValueError) as e:
        Sb.map_features_nucl2prot(Sb.make_copy(tester), Sb.make_copy(sb_objects[6]), mode="foo")
    assert "'mode' must be either 'key' or 'position'" in str(e.value)

    with pytest.raises(ValueError) as e:
        Sb.pull_recs(tester, "α[1-8]")
        Sb.map_features_nucl2prot(Sb.make_copy(tester), Sb.make_copy(sb_objects[6]), mode="list")
    assert "The two input files do not contain the same number of sequences" in str(e.value)

# ######################  '-fp2n', '--map_features_prot2nucl' ###################### #
hashes = ["3ebc92ca11505489cab2453d2ebdfcf2", "feceaf5e17935afb100b4b6030e27fee",
          "bfd36942768cf65c473b3aaebb83e4fa", "9ba4af4e5dd0bf4a445d173604b92996",
          "c178763aa9596e341bbbc088f1f791c9", "84cc7ecb54603c5032737e5263a52bd3"]

hashes = [(Sb.make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("_dna,next_hash", hashes)
def test_map_features_prot2nucl(_dna, next_hash):
    tester = Sb.map_features_prot2nucl(Sb.make_copy(sb_objects[7]), _dna)
    tester.out_format = "gb"
    assert seqs_to_hash(tester) == next_hash


def test_map_features_prot2nucl_2(capsys):
    prot_tester = Sb.make_copy(sb_objects[7])
    Sb.pull_record_ends(prot_tester, 450)
    prot_tester = Sb.annotate(prot_tester, "foo", [(20, 40), (50, 60)], pattern="α9")
    mapped = Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), Sb.make_copy(sb_objects[0]))
    mapped.out_format = "gb"
    assert seqs_to_hash(mapped) == "552b31b8a7068d12d6c55c7e5d293c54"
    out, err = capsys.readouterr()
    assert err == "Warning: size mismatch between aa and nucl seqs for Mle-Panxα7A --> 450, 1875\n"

    prot_tester = Sb.make_copy(sb_objects[7])
    prot_tester.records = sorted(prot_tester.records, key=lambda x: x.id)
    dna_tester = Sb.make_copy(sb_objects[1])
    dna_tester.records = sorted(dna_tester.records, key=lambda x: x.id)
    Sb.rename(prot_tester, "α4", "A4")
    Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), dna_tester, mode="list")
    assert seqs_to_hash(dna_tester) == "c32f43cc5205867c0eb1d3873e27319b"

    Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), dna_tester, mode="key")
    assert seqs_to_hash(dna_tester) == "c32f43cc5205867c0eb1d3873e27319b"
    out, err = capsys.readouterr()
    assert string2hash(err) == "c1f51587b07d2cc6156b8aac07384834"

    mock_obj = type('test', (object,), {})()
    mock_obj.start = 2
    mock_obj.end = 6
    prot_tester.records[0].features[0].location = mock_obj
    with pytest.raises(TypeError) as e:
        Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), Sb.make_copy(sb_objects[0]))
    assert "FeatureLocation or CompoundLocation object required." in str(e.value)

    prot_tester.records[0].features[0].location = prot_tester.records[1].features[0].location
    with pytest.raises(ValueError) as e:
        Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), Sb.make_copy(sb_objects[6]), mode="foo")
    assert "'mode' must be either 'key' or 'position'" in str(e.value)

    with pytest.raises(ValueError) as e:
        Sb.pull_recs(prot_tester, "α[1-8]")
        Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), Sb.make_copy(sb_objects[6]), mode="list")
    assert "The two input files do not contain the same number of sequences" in str(e.value)


# #####################  '-mg', '--merge' ###################### ##
def test_merge():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds_dummy_features.gb"))
    tester = Sb.merge(tester, Sb.make_copy(sb_objects[1]))
    assert seqs_to_hash(tester) == "bae5aeb130b3d5319378a122a6f61df5"

    tester.records[0].seq = tester.records[1].seq
    with pytest.raises(RuntimeError) as e:
        Sb.merge(tester, Sb.make_copy(sb_objects[1]))
    assert "Sequence mismatch for record 'Mle-Panxα1'" in str(e.value)


# ######################  '-mw', '--molecular_weight' ###################### #
def test_molecular_weight():
    # Unambiguous DNA
    tester = Sb.molecular_weight(Sb.make_copy(sb_objects[1]))
    assert tester.molecular_weights['masses_ds'][0] == 743477.1
    assert tester.molecular_weights['masses_ss'][0] == 371242.6
    assert seqs_to_hash(tester) == "e080cffef0ec6c5e8eada6f57bbc35f9"
    # Ambiguous DNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("ambiguous_dna.fa")))
    assert tester.molecular_weights['masses_ds'][0] == 743477.08
    assert tester.molecular_weights['masses_ss'][0] == 371202.59
    # Unambiguous RNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("Mnemiopsis_rna.fa")))
    assert tester.molecular_weights['masses_ss'][0] == 387372.6
    # Ambiguous RNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("ambiguous_rna.fa")))
    assert tester.molecular_weights['masses_ss'][0] == 387371.6
    # Protein
    tester = Sb.molecular_weight(Sb.make_copy(sb_objects[7]))
    assert tester.molecular_weights['masses_ss'][0] == 45692.99
    assert seqs_to_hash(tester) == "fb1a66b7eb576c0584fc7988c45b6a18"


# ######################  '-ns', '--num_seqs' ###################### #
seq_counts = [(sb_objects[0], 13), (sb_objects[1], 13), (sb_objects[2], 13), (sb_objects[3], 8),
              (sb_objects[4], 13), (sb_objects[5], 13), (sb_objects[6], 13), (sb_objects[9], 8)]


@pytest.mark.parametrize("seqbuddy, num", seq_counts)
def test_num_seqs(seqbuddy, num):
    assert Sb.num_seqs(seqbuddy) == num


def test_empty_file():
    tester = Sb.SeqBuddy(resource("blank.fa"))
    assert type(tester) == Sb.SeqBuddy
    assert len(tester.records) == 0

# ######################  '-ofa', '--order_features_alphabetically' ###################### #
fwd_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "21547b4b35e49fa37e5c5b858808befb",
              "cb1169c2dd357771a97a02ae2160935d", "503e23720beea201f8fadf5dabda75e4",
              "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "d0297078b4c480a49b6da5b719310d0e",
              "17ff1b919cac899c5f918ce8d71904f6", "968ed9fa772e65750f201000d7da670f",
              "ce423d5b99d5917fbef6f3b47df40513", "c0dce60745515b31a27de1f919083fe9"]

rev_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "3b718ec3cb794bcb658d900e517110cc",
              "cb1169c2dd357771a97a02ae2160935d", "503e23720beea201f8fadf5dabda75e4",
              "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "c6a788d8ea916964605ac2942c459c9b",
              "17ff1b919cac899c5f918ce8d71904f6", "968ed9fa772e65750f201000d7da670f",
              "ce423d5b99d5917fbef6f3b47df40513", "c0dce60745515b31a27de1f919083fe9"]
hashes = [(Sb.make_copy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)
def test_order_features_alphabetically(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_features_alphabetically(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_features_alphabetically(seqbuddy, reverse=True)
    assert seqs_to_hash(tester) == rev_hash


# ######################  '-ofp', '--order_features_by_position' ###################### #
fwd_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "2e02a8e079267bd9add3c39f759b252c",
              "cb1169c2dd357771a97a02ae2160935d", "503e23720beea201f8fadf5dabda75e4",
              "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "7a8e25892dada7eb45e48852cbb6b63d",
              "17ff1b919cac899c5f918ce8d71904f6", "968ed9fa772e65750f201000d7da670f",
              "ce423d5b99d5917fbef6f3b47df40513", "c0dce60745515b31a27de1f919083fe9"]

rev_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "4345a14fe27570b3c837c30a8cb55ea9",
              "cb1169c2dd357771a97a02ae2160935d", "503e23720beea201f8fadf5dabda75e4",
              "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "9e7c2571db1386bba5983365ae235e1b",
              "17ff1b919cac899c5f918ce8d71904f6", "968ed9fa772e65750f201000d7da670f",
              "ce423d5b99d5917fbef6f3b47df40513", "c0dce60745515b31a27de1f919083fe9"]
hashes = [(Sb.make_copy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)
def test_order_features_position(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_features_by_position(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_features_by_position(seqbuddy, reverse=True)
    assert seqs_to_hash(tester) == rev_hash


# ######################  '-oi', '--order_ids' ###################### #
fwd_hashes = ["e9090efdd362d527a115049dfced42cd", "c0d656543aa5d20a266cffa790c035ce",
              "132757da01b3caf174d024efdb2c3acd", "3c49bdc1b0fe4e1d6bfc148eb0293e21",
              "684a99151e8cbf1e5eb96e44b875ba08", "c06985b566277a29e598dea0dc41baef"]
rev_hashes = ["fcb016fff87d26822fa518d62e355c65", "2507c667a304fdc003bc68255e094d7b",
              "286bac7a213997924203622c3357457c", "d6e79a5faeaff396aa7eab0b460c3eb9",
              "a5fcbd95e837b4807408124c2396db6e", "20b43a724670a1151b1f0418478046ef"]

hashes = [(Sb.make_copy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)
def test_order_ids1(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_ids(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_ids(seqbuddy, reverse=True)
    assert seqs_to_hash(tester) == rev_hash


def test_order_ids2():
    seqbuddy = sb_resources.get_one("p n")
    Sb.rename(seqbuddy, "Mle-Panxα4", "Mle004-Panxα4")
    Sb.rename(seqbuddy, "Mle-Panxα5", "Mle05-Panxα5")
    Sb.rename(seqbuddy, "Mle-Panxα9", "aMle-PanxαBlahh")
    Sb.order_ids(seqbuddy)
    assert seqs_to_hash(seqbuddy) == "5c1316e18205432b044101e720646cd5"


# ######################  '-oir', '--order_ids_randomly' ###################### #
@pytest.mark.parametrize("seqbuddy", [Sb.make_copy(x) for x in sb_objects])
def test_order_ids_randomly(seqbuddy):
    tester = Sb.order_ids_randomly(Sb.make_copy(seqbuddy))
    assert seqs_to_hash(seqbuddy) != seqs_to_hash(tester)
    assert seqs_to_hash(Sb.order_ids(tester)) == seqs_to_hash(tester)


def test_order_ids_randomly2():
    tester = Sb.make_copy(sb_objects[0])
    for _ in range(15):  # This will fail to repeat the while loop only ~5% of the time
        Sb.pull_recs(tester, "α[789]")
        assert seqs_to_hash(tester) != seqs_to_hash(Sb.order_ids_randomly(tester))

    Sb.pull_recs(tester, "α[89]")
    assert seqs_to_hash(tester) != seqs_to_hash(Sb.order_ids_randomly(tester))

    Sb.pull_recs(tester, "α[9]")
    assert seqs_to_hash(tester) == seqs_to_hash(Sb.order_ids_randomly(tester))

    tester = Sb.SeqBuddy(tester.records * 3)
    assert seqs_to_hash(tester) == seqs_to_hash(Sb.order_ids_randomly(tester))


# #####################  '-prr', '--pull_random_recs' ###################### ##
@pytest.mark.parametrize("seqbuddy", sb_objects)
def test_pull_random_recs(seqbuddy):
    tester = Sb.pull_random_recs(Sb.make_copy(seqbuddy))
    orig_seqs = tester.to_dict()
    assert len(tester.records) == 1
    assert tester.records[0].id in orig_seqs


# #####################  '-pre', '--pull_record_ends' ###################### ##
def test_pull_record_ends():
    tester = Sb.pull_record_ends(Sb.make_copy(sb_objects[1]), 10)
    assert seqs_to_hash(tester) == 'd46867e4ca7a9f474c45473fc3495413'

    tester = Sb.pull_record_ends(Sb.make_copy(sb_objects[1]), 2000)
    assert seqs_to_hash(tester) == '908744b00d9f3392a64b4b18f0db9fee'

    tester = Sb.pull_record_ends(Sb.make_copy(sb_objects[1]), -10)
    assert seqs_to_hash(tester) == 'd7970570d65872993df8a3e1d80f9ff5'

    tester = Sb.pull_record_ends(Sb.make_copy(sb_objects[1]), -2000)
    assert seqs_to_hash(tester) == '908744b00d9f3392a64b4b18f0db9fee'

    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    with pytest.raises(ValueError):
        Sb.pull_record_ends(Sb.make_copy(seqbuddy), 'foo')


# ######################  '-pr', '--pull_records' ###################### #
pr_hashes = ["5b4154c2662b66d18776cdff5af89fc0", "e196fdc5765ba2c47f97807bafb6768c", "bc7dbc612bc8139eba58bf896b7eaf2f",
             "7bb4aac2bf50381ef1d27d82b7dd5a53", "72328eddc751fd79406bb911dafa57a2", "b006b40ff17ba739929448ae2f9133a6"]
pr_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(pr_hashes)]


@pytest.mark.parametrize("seqbuddy, next_hash", pr_hashes)
def test_pull_recs(seqbuddy, next_hash):
    tester = Sb.pull_recs(seqbuddy, 'α2')
    assert seqs_to_hash(tester) == next_hash


# #####################  '-prg', '--purge' ###################### ##
def test_purge():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_pep.fa"))
    Sb.purge(tester, 200)
    assert seqs_to_hash(tester) == 'b21b2e2f0ca1fcd7b25efbbe9c08858c'


# ######################  '-ri', '--rename_ids' ###################### #
hashes = ["8b4a9e3d3bb58cf8530ee18b9df67ff1", "78c73f97117bd937fd5cf52f4bd6c26e", "243024bfd2f686e6a6e0ef65aa963494",
          "98bb9b57f97555d863054ddb526055b4", "2443c47a712f19099e94fc015dc980a9", "65196fd4f2a4e339e1545f6ed2a6acc3"]
hashes = [(Sb.make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_rename_ids(seqbuddy, next_hash):
    tester = Sb.rename(seqbuddy, 'Panx', 'Test', 0)
    assert seqs_to_hash(tester) == next_hash


def test_rename_ids2():
    tester = Sb.make_copy(sb_objects[0])
    Sb.rename(tester, "Panxα([1-6])", "testing\\1")
    assert seqs_to_hash(tester) == "f0d7ec055b3b2d9ec11a86634b32e9ef"

    tester = Sb.make_copy(sb_objects[0])
    Sb.rename(tester, "Panxα([1-6])", "testing\\1", -1)
    assert seqs_to_hash(tester) == "f0d7ec055b3b2d9ec11a86634b32e9ef"

    tester = Sb.make_copy(sb_objects[0])
    Sb.rename(tester, "[a-z]", "?", 2)
    assert seqs_to_hash(tester) == "08aaee2e7b9997b512c7d1b2fe748d40"

    tester = Sb.make_copy(sb_objects[0])
    Sb.rename(tester, "[a-z]", "?", -2)
    assert seqs_to_hash(tester) == "9b3946afde20c991099463d099be22e0"

    tester = Sb.make_copy(sb_objects[0])
    Sb.rename(tester, "[A-Z]", "?", -20)
    assert seqs_to_hash(tester) == "451993d7e816881e2700697263b1d8fa"

    tester = Sb.make_copy(sb_objects[0])
    Sb.rename(tester, "[a-z]", "?", 2, store_old_id=True)
    assert seqs_to_hash(tester) == "959fe04c0366c9a143052a02f090707e"

    with pytest.raises(AttributeError) as e:
        Sb.rename(tester, "[a-z]", "\\1?", -2)
    assert "There are more replacement match" in str(e)


# ##################### '-rs', 'replace_subseq' ###################### ##
def test_replace_subsequence():
    tester = Sb.make_copy(sb_objects[0])
    Sb.replace_subsequence(tester, "atg(.{5}).{3}", "FOO\\1BAR")
    assert seqs_to_hash(tester) == "f12707c2b0ef866f0039bac96abb29e0"


# ######################  '-rc', '--reverse_complement' ###################### #
hashes = ["e77be24b8a7067ed54f06e0db893ce27", "47941614adfcc5bd107f71abef8b3e00", "f549c8dc076f6b3b4cf5a1bc47bf269d",
          "0dd20827c011a0a7a0e78881b38ae06a", "0b954f4a263bf38ddeac61ab54f77dc2", "0d6b7deda824b4fc42b65cb87e1d4d14"]
hashes = [(Sb.make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_reverse_complement(seqbuddy, next_hash):
    tester = Sb.reverse_complement(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_reverse_complement_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError) as e:
        Sb.reverse_complement(sb_objects[6])
    assert str(e.value) == "SeqBuddy object is protein. Nucleic acid sequences required."

# ######################  '-sfr', '--select_frame' ###################### #
hashes = [(0, 1, "b831e901d8b6b1ba52bad797bad92d14"), (0, 2, "2de033b2bf2327f2795fe425db0bd78f"),
          (0, 3, "1c29898d4964e0d1b03207d7e67e1958"), (1, 1, "908744b00d9f3392a64b4b18f0db9fee"),
          (1, 2, "08fe54a87249f5fb9ba22ff6d0053787"), (1, 3, "cfe2d405487d69dceb2a11dd44ceec59"),
          (2, 1, "cb1169c2dd357771a97a02ae2160935d"), (2, 2, "87d784f197b55f812d2fc82774da43d1"),
          (2, 3, "5d6d2f337ecdc6f9a85e981c975f3e08")]
hashes = [(Sb.make_copy(sb_objects[sb_indx]), _hash, frame) for sb_indx, frame, _hash in hashes]


@pytest.mark.parametrize("seqbuddy,next_hash,shift", hashes)
def test_select_frame(seqbuddy, next_hash, shift):
    tester = Sb.select_frame(seqbuddy, shift)
    assert seqs_to_hash(tester) == next_hash


def test_select_frame_edges():
    tester = Sb.select_frame(Sb.make_copy(sb_objects[0]), 2)
    temp_file = br.TempFile()
    tester.write(temp_file.path)
    tester = Sb.select_frame(Sb.SeqBuddy(temp_file.path), 1)
    assert seqs_to_hash(tester) == "b831e901d8b6b1ba52bad797bad92d14"

    tester = Sb.select_frame(Sb.make_copy(sb_objects[1]), 2)
    tester = Sb.select_frame(tester, 1)
    assert seqs_to_hash(tester) == "908744b00d9f3392a64b4b18f0db9fee"

    with pytest.raises(TypeError) as e:  # If protein is input
        Sb.select_frame(sb_objects[6], 2)
    assert "Select frame requires nucleic acid, not protein." in str(e.value)


# ##################### '-ss', 'shuffle_seqs' ###################### ##
def test_shuffle_seqs():
    for seqbuddy in sb_objects:
        tester1 = Sb.shuffle_seqs(Sb.make_copy(seqbuddy))
        tester2 = Sb.make_copy(seqbuddy)
        assert seqs_to_hash(tester1) != seqs_to_hash(tester2)

        for indx, record in enumerate(tester1.records):
            assert sorted(record.seq) == sorted(tester2.records[indx].seq)


# #####################  make_groups' ###################### ##
def test_make_groups():
    tester = Sb.SeqBuddy(resource("Cnidaria_pep.nexus"))
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
# Only fasta and genbank
hashes = ["95cf24202007399e6ccd6e6f33ae012e", "0b5daa810e1589c3973e1436c40baf08"]
hashes = [(Sb.make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_translate6frames(seqbuddy, next_hash):
    tester = Sb.translate6frames(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_translate6frames_pep_exception():
    with pytest.raises(TypeError):
        Sb.translate6frames(Sb.make_copy(sb_objects[6]))


# ######################  '-tr', '--translate' ###################### #
hashes = ["06893e14839dc0448e6f522c1b8f8957", "e8840e22096e933ce10dbd91036f3fa5", "f3339e0193c10427f017dd8f6bd81d7e"]
hashes = [(Sb.make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]
hashes.append((Sb.make_copy(sb_objects[13]), "648ccc7c3400882be5bf6e8d9781f74e"))


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_translate(seqbuddy, next_hash):
    tester = Sb.translate_cds(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_translate_edges_and_exceptions(capsys):
    with pytest.raises(TypeError):
        Sb.translate_cds(Sb.make_copy(sb_objects[6]))

    tester = Sb.SeqBuddy("ATTCGTTAACGCTAGCGTCG", in_format="raw")
    tester = Sb.translate_cds(tester)
    assert str(tester) == ">raw_input\nIR*R*R\n"
    out, err = capsys.readouterr()
    assert err == "Warning: size mismatch between aa and nucl seqs for raw_input --> 20, 6\n"

    tester = Sb.select_frame(Sb.make_copy(sb_objects[1]), 3)
    tester = Sb.translate_cds(tester)
    assert seqs_to_hash(tester) == "68ca15f5ac737e4a4ca65a67ad2dc897"
    out, err = capsys.readouterr()
    assert string2hash(err) == "9e2a0b4b03f54c209d3a9111792762df"


# ######################  '-tmd', '--transmembrane_domains' ###################### #
@pytest.mark.internet
@pytest.mark.slow
def test_transmembrane_domains_pep():
    tester = sb_resources.get_one("p f")
    Sb.pull_recs(tester, "Panxα[234]")
    tester = Sb.transmembrane_domains(tester, quiet=True)
    tester.out_format = "gb"
    assert seqs_to_hash(tester) == "7285d3c6d60ccb656e39d6f134d1df8b"


@pytest.mark.internet
@pytest.mark.slow
def test_transmembrane_domains_cds():
    TEMP_DIR.subdir("topcons")
    tester = sb_resources.get_one("d f")
    Sb.pull_recs(tester, "Panxα[234]")
    tester = Sb.transmembrane_domains(tester, quiet=True, keep_temp="%s/topcons" % TEMP_DIR.path)
    tester.out_format = "gb"
    assert seqs_to_hash(tester) == "e5c9bd89810a39090fc3326e51e1ac6a"
    _root, dirs, files = next(br.walklevel("%s/topcons" % TEMP_DIR.path))
    _root, dirs, files = next(br.walklevel("%s/topcons/%s" % (TEMP_DIR.path, dirs[0])))
    assert files
"""