#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" tests basic functionality of AlignBuddy class """
import pytest
import os
from os.path import join
import shutil
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.AlignIO import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
import AlignBuddy as Alb
import SeqBuddy as Sb
import buddy_resources as br
from .. import __init__
from collections import OrderedDict

TEMPDIR = br.TempDir()
RES_PATH = __init__.RESOURCE_PATH


def mock_hash_ids(sb, *args):
    # This will only work for Mnemiopsis full files
    sb.hash_map = OrderedDict([('KWpHHnwG', 'Mle-Panxα9'), ('ui9bSgqc', 'Mle-Panxα7A'), ('xDsAi9jU', 'Mle-Panxα1'),
                               ('RraJkRhu', 'Mle-Panxα3'), ('PQKaRX5O', 'Mle-Panxα12'), ('hEZa4U1t', 'Mle-Panxα11'),
                               ('oh8QNn5a', 'Mle-Panxα4'), ('ejCQN1iZ', 'Mle-Panxα8'), ('9PYuN1Vz', 'Mle-Panxα6'),
                               ('yUf5nZfJ', 'Mle-Panxα10B'), ('xwOHuotu', 'Mle-Panxα5'), ('eDI1PEGu', 'Mle-Panxα2'),
                               ('sfdhEBIk', 'Mle-Panxα10A')])
    return sb


# ##########################################  '-al', '--alignment_lengths' ########################################### #
def test_alignment_lengths(alb_resources):
    lengths = Alb.alignment_lengths(alb_resources.get_one("m p c"))
    assert lengths[0] == 481
    assert lengths[1] == 683

    lengths = Alb.alignment_lengths(alb_resources.get_one("m d s"))
    assert lengths[0] == 2043
    assert lengths[1] == 1440


# ##############################################  '-bts', '--bootstrap' ############################################## #
def test_bootstrap(alb_resources, hf):
    # Test an amino acid file
    tester = Alb.bootstrap(alb_resources.get_one("m p py"), r_seed=12345)
    assert hf.buddy2hash(tester) == "ce9009b430ad0fe7f151477985e6f922"

    tester = Alb.bootstrap(alb_resources.get_one("m p py"), 3, r_seed=12345)
    assert hf.buddy2hash(tester) == "e6d5f30c3a7f53ec7899618a045b017d"


# ##############################################  '-cs', '--clean_seqs' ############################################## #
def test_clean_seqs(alb_resources, hf):
    # Test an amino acid file
    tester = Alb.clean_seq(alb_resources.get_one("m p py"))
    assert hf.buddy2hash(tester) == "07a861a1c80753e7f89f092602271072"

    tester = Alb.clean_seq(Alb.AlignBuddy("%sambiguous_dna_alignment.fa" % hf.resource_path),
                           ambiguous=False, rep_char="X")
    assert hf.buddy2hash(tester) == "6755ea1408eddd0e5f267349c287d989"


# ###########################################  '-cta', '--concat_alignments' ######################################### #
def test_concat_alignments(alb_resources, hf):
    with pytest.raises(AttributeError) as e:
        Alb.concat_alignments(alb_resources.get_one("p o g"), '.*')
    assert "Please provide at least two alignments." in str(e)

    tester = alb_resources.get_one("o p g")
    tester.alignments.append(alb_resources.get_one("o p g").alignments[0])

    with pytest.raises(ValueError) as e:
        Alb.concat_alignments(tester, 'foo')
    assert "No match found for record" in str(e)

    with pytest.raises(ValueError) as e:
        Alb.concat_alignments(tester, 'Panx')
    assert "Replicate matches" in str(e)

    tester = Sb.SeqBuddy("%sCnidaria_pep.nexus" % hf.resource_path)
    Sb.pull_recs(tester, "Ccr|Cla|Hec")
    tester = Alb.AlignBuddy(str(tester))
    tester.alignments.append(tester.alignments[0])

    tester.set_format("gb")
    tester2 = Alb.concat_alignments(Alb.make_copy(tester), suppress_position=True)
    assert hf.buddy2hash(tester2) == 'a4bf16d88352935848dac29f4afdc871', print(tester2)

    tester.set_format("nexus")
    tester2 = Alb.concat_alignments(Alb.make_copy(tester))
    assert hf.buddy2hash(tester2) == '32a507107b7dcd044ea7760c8812441c'

    tester.set_format("gb")
    tester2 = Alb.concat_alignments(Alb.make_copy(tester), "(.).(.)-Panx(.)")
    assert hf.buddy2hash(tester2) == 'cd2b6594b22c431aea67fa45899f933a'

    tester.set_format("gb")
    tester2 = Alb.concat_alignments(Alb.make_copy(tester), "(.).(.)-Panx(.)")
    assert hf.buddy2hash(tester2) == 'cd2b6594b22c431aea67fa45899f933a'

    tester.set_format("gb")
    tester2 = Alb.concat_alignments(Alb.make_copy(tester), "...", "Panx.*")
    assert hf.buddy2hash(tester2) == 'e49b26f695c910a93f93d70563fd9dd9'

    tester.set_format("gb")
    tester2 = Alb.concat_alignments(Alb.make_copy(tester), "...", "(P)an(x)(.)")
    assert hf.buddy2hash(tester2) == '3abfdf2217050ac2170c0de27352a8c6'

    shorten = Alb.delete_records(Alb.make_copy(tester), "Ccr")
    tester.alignments[1] = shorten.alignments[1]
    assert hf.buddy2hash(Alb.concat_alignments(Alb.make_copy(tester))) == '685f24ee1fc88860dd9465035040c91e'


# ###########################################  '-con', '--consensus' ############################################ #
hashes = [('o d g', '36ed26ded1f4e146ad8018a2dc31f0b2'), ('o d n', 'b7a49ccc640b088f2fe7de67bad30f06'),
          ('o p py', 'c32119f34633d3956b3a1d3ac578869c'), ('m p s', 'd1a8f7e629a020f5130373d7af65f9d9')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_simple_consensus(alb_resources, hf, key, next_hash):
    tester = Alb.consensus_sequence(alb_resources.get_one(key), "simple")
    assert hf.buddy2hash(tester) == next_hash


hashes = [('o d g', 'c85a1726c3a6f8c3a402475d4f8fce27'), ('o d n', '15bb6c1959a6a5fee1a385921c39a4f5'),
          ('o p py', '714bc7c67332fa96ef12cfeada88834c'), ('m p s', 'bf50c95916e9d62c95a460bbc517c053')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_weighted_consensus(alb_resources, hf, key, next_hash):
    tester = Alb.consensus_sequence(alb_resources.get_one(key))
    assert hf.buddy2hash(tester) == next_hash


def test_consensus_error(alb_resources):
    tester = alb_resources.get_one('o d f')
    with pytest.raises(ValueError) as err:
        Alb.consensus_sequence(tester, mode="foo")
    assert "No valid consensus" in str(err)


# ######################################  '-dinv', '--delete_invariant_sites' ####################################### #
def test_delete_invariant_sites(alb_resources, hf, alb_odd_resources):
    tester = Alb.AlignBuddy(alb_odd_resources['dna']['single']['ambiguous'])
    tester = Alb.delete_invariant_sites(tester)
    assert hf.buddy2hash(tester) == "27233a416437eabc72aa5d57cb695036"

    tester = alb_resources.get_one("o p py")
    tester = Alb.delete_invariant_sites(tester, consider_ambiguous=False)
    assert hf.buddy2hash(tester) == "f0b16bb8133bfc9e29ad43bdfc4ad2ee"

    tester = Alb.delete_invariant_sites(tester)
    assert hf.buddy2hash(tester) == "c13031016c1f7382e808bd4e68d8f406"

    tester.alignments.append([])  # Catch empty alignment
    tester = Alb.delete_invariant_sites(tester)
    assert hf.buddy2hash(tester) == "c13031016c1f7382e808bd4e68d8f406"


# ###########################################  '-dr', '--delete_records' ############################################ #
hashes = [('o d g', 'c22d5cbef500d8baed8cead1d5fe9628'), ('o d n', '355a98dad5cf382797eb907e83940978'),
          ('o d py', 'fe9a2776558f3fe9a1732c777c4bc9ac'), ('o d s', '35dc92c4f4697fb508eb1feca43d9d75'),
          ('o r n', '96e6964115200d46c7cb4eb975718304'), ('o p g', 'ad384f639efda97e300dbe6a70dd9f7c'),
          ('o p n', '1cfaa4109c5db8fbfeaedabdc57af655'), ('o p py', '1d0e7b4d8e89b42b0ef7cc8c40ed1a93'),
          ('o p s', '1578d98739d2aa6196463957c7b408fa'), ('m d py', 'db4ed247b40707e8e1f0622bb420733b'),
          ('m d s', 'de5beddbc7f0a7f8e3dc2d5fd43b7b29'), ('m p py', '31f91f7dc548e4b075bfb0fdd7d5c82c'),
          ('m p s', '043e35023b355ed560166db9130cfe30')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_delete_records(alb_resources, hf, key, next_hash):
    tester = Alb.delete_records(alb_resources.get_one(key), "α[1-5]|β[A-M]")
    assert hf.buddy2hash(tester) == next_hash


# ######################  'd2r', '--transcribe' and 'r2d', '--reverse_transcribe' ###################### #
hashes = [('o d g', '81a17d48b752c5a760e98cf8e665d086', '842d9c911a33c0fd0484383eabefb0fe'),
          ('o d n', 'e531dc31f24192f90aa1f4b6195185b0', 'cb1169c2dd357771a97a02ae2160935d'),
          ('o d py', 'e55bd18b6d82a7fc3150338173e57e6a', '503e23720beea201f8fadf5dabda75e4'),
          ('o d s', '45b511f34653e3b984e412182edee3ca', '228e36a30e8433e4ee2cd78c3290fa6b'),
          ('m d py', '16cb082f5cd9f103292ccea0c4d65a06', '42679a32ebd93b628303865f68b0293d'),
          ('m d s', 'd81dae9714a553bddbf38084f7a8e00e', 'ae352b908be94738d6d9cd54770e5b5d')]


@pytest.mark.parametrize("key,d2r_hash,r2d_hash", hashes)
def test_transcribe(alb_resources, hf, key, d2r_hash, r2d_hash):
    tester = Alb.dna2rna(alb_resources.get_one(key))
    assert hf.buddy2hash(tester) == d2r_hash
    tester = Alb.rna2dna(tester)
    assert hf.buddy2hash(tester) == r2d_hash


def test_transcribe_exceptions(alb_resources):
    with pytest.raises(TypeError) as e:
        Alb.dna2rna(alb_resources.get_one("o p s"))
    assert "TypeError: DNA sequence required, not IUPACProtein()." in str(e)

    with pytest.raises(TypeError) as e:
        Alb.dna2rna(alb_resources.get_one("o r n"))
    assert "TypeError: DNA sequence required, not IUPACAmbiguousRNA()." in str(e)


def test_reverse_transcribe_exceptions(alb_resources):  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError) as e:
        Alb.rna2dna(alb_resources.get_one("o p s"))
    assert "TypeError: RNA sequence required, not IUPACProtein()." in str(e)

    with pytest.raises(TypeError) as e:
        Alb.rna2dna(alb_resources.get_one("o d s"))
    assert "TypeError: RNA sequence required, not IUPACAmbiguousDNA()." in str(e)


# ###########################################  '-et', '--enforce_triplets' ########################################### #
hashes = [('o d g', 'd30529911c2ffdfb49152797225e3ff0'), ('o d n', 'c907d29434fe2b45db60f1a9b70f110d'),
          ('o d py', 'b6cf61c86588023b58257c9008c862b5'), ('o r n', '0ed7383ab2897f8350c2791739f0b0a4'),
          ('m d py', '669ffc4fa602fb101c559cb576bddee1')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_enforce_triplets(key, next_hash, alb_resources, hf):
    tester = Alb.enforce_triplets(alb_resources.get_one(key))
    assert hf.buddy2hash(tester) == next_hash


def test_enforce_triplets_error(alb_resources):
    with pytest.raises(TypeError) as e:
        Alb.enforce_triplets(alb_resources.get_one("m p c"))
    assert "Nucleic acid sequence required, not protein." in str(e)

    with pytest.raises(TypeError) as e:
        tester = Alb.enforce_triplets(alb_resources.get_one("m d pr"))
        tester.alignments[0][0].seq = Seq("MLDILSKFKGVTPFKGITIDDGWDQLNRSFMFVLLVVMGTTVTVRQYTGSVISCDGFKKFGSTFAEDYCWTQGLY",
                                          alphabet=IUPAC.protein)
        Alb.enforce_triplets(tester)
    assert "Record 'Mle-Panxα9' is protein. Nucleic acid sequence required." in str(e)


# ######################  '-efs', '--extract_feature_sequences' ###################### #
def test_extract_feature_sequences(alb_resources, hf):
    tester = alb_resources.get_one("o d g")
    tester = Alb.extract_feature_sequences(tester, "CDS")
    assert hf.buddy2hash(tester) == "842d9c911a33c0fd0484383eabefb0fe"

    tester = alb_resources.get_one("o d g")
    tester = Alb.extract_feature_sequences(tester, ["TMD"])
    assert hf.buddy2hash(tester) == "cc7a1c6a22f721ec0668fc8ea6b23429"

    tester = alb_resources.get_one("o d g")
    tester = Alb.extract_feature_sequences(tester, ["TMD1", "splice_a"])
    assert hf.buddy2hash(tester) == "497d536b1be9a90ef0ef75281d0c867f"

    tester = alb_resources.get_one("o d g")
    tester = Alb.extract_feature_sequences(tester, ["TMD2:TMD3"])
    assert hf.buddy2hash(tester) == "07773f4fb1dc430c0c3ce6cd5a799439"

    tester = alb_resources.get_one("o d g")
    tester = Alb.extract_feature_sequences(tester, ["TMD3:TMD2"])
    assert hf.buddy2hash(tester) == "07773f4fb1dc430c0c3ce6cd5a799439"

    tester = alb_resources.get_one("o d g")
    tester = Alb.extract_feature_sequences(tester, ["TMD2:foo"])
    assert hf.buddy2hash(tester) == "ac15492b38ca2ac4baa63e63a9b747f7"

    tester = alb_resources.get_one("o d g")
    tester = Alb.extract_feature_sequences(tester, "foo")
    assert hf.buddy2hash(tester) == "ac15492b38ca2ac4baa63e63a9b747f7"

    tester = alb_resources.get_one("o d g")
    tester = Alb.extract_feature_sequences(tester, [])
    assert hf.buddy2hash(tester) == "ac15492b38ca2ac4baa63e63a9b747f7"


# ###########################################  'er', '--extract_regions' ############################################ #
hashes = [('o d g', 'd06c46d7458f8b9a90aba14b83cdb329'), ('o d n', '10ca718b74f3b137c083a766cb737f31'),
          ('o d py', 'd738a9ab3ab200a7e013177e1042e86c'), ('o p g', '19b5ae01128233d622c4a638f00d612e'),
          ('o p n', '5f400edc6f0990c0cd6eb52ae7687e39'), ('o p py', '69c9ad73ae02525150d4682f9dd68093'),
          ('m d py', 'd06ba679c8a686c8f077bb460a4193b0'), ('m p py', '8151eeda36b9a170512709829d70230b')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_extract_regions(key, next_hash, alb_resources, hf):
    tester = Alb.extract_regions(alb_resources.get_one(key), "0:50")
    assert hf.buddy2hash(tester) == next_hash, print(tester)


def test_extract_regions_singlets(alb_resources, hf):
    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "0")
    assert hf.buddy2hash(tester) == "4e6b813255000b181357dca293c6e9c4"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "1")
    assert hf.buddy2hash(tester) == "4e6b813255000b181357dca293c6e9c4"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "-10000000")
    assert hf.buddy2hash(tester) == "4e6b813255000b181357dca293c6e9c4"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), ",1/")
    assert hf.buddy2hash(tester) == "4e6b813255000b181357dca293c6e9c4"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "1000000")
    assert hf.buddy2hash(tester) == "da0849f4fbdf833223d010427f25d5e6"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "2,5,9,-5")
    assert hf.buddy2hash(tester) == "f26c7778136905caadbe5f5f368c63f8"


def test_extract_regions_ranges(alb_resources, hf):
    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "0:10")
    assert hf.buddy2hash(tester) == "65214ccec8e0be7459b8cac7686949ed"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "1:10")
    assert hf.buddy2hash(tester) == "65214ccec8e0be7459b8cac7686949ed"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "10:1")
    assert hf.buddy2hash(tester) == "65214ccec8e0be7459b8cac7686949ed"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), ":10")
    assert hf.buddy2hash(tester) == "65214ccec8e0be7459b8cac7686949ed"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "-10:")
    assert hf.buddy2hash(tester) == "e1092c191c8a1721d399f0e1a03af015"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "40:75,89:100,432:-45")
    assert hf.buddy2hash(tester) == "2fe7d2b15d28a3778b52e30171cea6bf"


def test_extract_regions_mth_of_nth(alb_resources, hf):
    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "1/50")
    assert hf.buddy2hash(tester) == "b4f86edc6869a50b9ac273abaa22b894"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "-1/50")
    assert hf.buddy2hash(tester) == "9c98929737aa6673eba9d6ad1486d50b"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "1/-500")
    assert hf.buddy2hash(tester) == "76cc5b13425726b6818de1ec8a09d4e8"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "50/1")
    assert hf.buddy2hash(tester) == "46388b175b31b81f47199ae6327768af"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "50/25")
    assert hf.buddy2hash(tester) == "4f42e487d6f8917451f7e838ff168663"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "1:5/50")
    assert hf.buddy2hash(tester) == "2d0dca7b78db5b5a51f987ca60aa2fd0"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "-5:/50")
    assert hf.buddy2hash(tester) == "518557f3e66d0e45747deab55cffcace"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), ":5/50")
    assert hf.buddy2hash(tester) == "2d0dca7b78db5b5a51f987ca60aa2fd0"

    tester = Alb.extract_regions(alb_resources.get_one("o p g"), "1:10,1/50,-1")
    assert hf.buddy2hash(tester) == "6858e0cca31d254d7336bde3c60f5622"


# ###########################################  'fa', '--faux_alignment' ############################################ #
hashes = [('d g', '3a9aa82548f3690485ddcbe1fff28476'), ('d n', 'e0e56affb50efe8b2305a83f348064c1'),
          ('d py', '7fe42a710e91979bd40a93207548bbda'), ('p g', '65fcda412cd548bb9c37f00429cb6fb7'),
          ('p n', '97b76e936b623ec31b5bc96626af5a18'), ('p py', '15c05d7ad6366919fa23e9e9228d1f07')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_faux_alignment(key, next_hash, sb_resources, hf):
    tester = Alb.faux_alignment(sb_resources.get_one(key), r_seed=12345)
    assert hf.buddy2hash(tester) == next_hash


# ###########################################  'ga', '--generate_alignment' ########################################## #
class MockPopen(object):
    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs
        if "my_mucsle" in args[0]:
            self.output = ["Robert C. Edgar".encode("utf-8"), "".encode("utf-8")]
        elif "clustalo" in args[0] and "-h" in args[0]:
            stdout = "Clustal Omega - 1.2.3 (AndreaGiacomo)"
            self.output = [stdout.encode("utf-8"), ''.encode("utf-8")]
        elif "clustalo" in args[0]:
            self.output = [None, None]
        elif "clustalw2" in args[0] and "-help" in args[0]:
            stdout = "CLUSTAL 2.1 Multiple Sequence Alignments"
            self.output = [stdout.encode("utf-8"), ''.encode("utf-8")]
        elif "clustalw2" in args[0]:
            _file = join(RES_PATH, "mock_resources", "test_clustalw2", "stdout.txt")
            with open(_file, "r", encoding="utf-8") as ifile:
                stdout = ifile.read()
            self.output = [stdout.encode("utf-8"), ''.encode("utf-8")]
        elif "mafft" in args[0] and "--help" in args[0]:
            stdout = "MAFFT v7.310"
            self.output = [stdout.encode("utf-8"), ''.encode("utf-8")]
        elif "mafft" in args[0]:
            _file = join(RES_PATH, "mock_resources", "test_mafft", "result")
            with open(_file, "r", encoding="utf-8") as ifile:
                stdout = ifile.read()
            self.output = [stdout.encode("utf-8"), ''.encode("utf-8")]
        elif "muscle" in args[0] and "-version" in args[0]:
            stdout = "MUSCLE v3.8.31 by Robert C. Edgar"
            self.output = [stdout.encode("utf-8"), ''.encode("utf-8")]
        elif "muscle" in args[0]:
            _file = join(RES_PATH, "mock_resources", "test_muscle", "result")
            with open(_file, "r", encoding="utf-8") as ifile:
                stdout = ifile.read()
            self.output = [stdout.encode("utf-8"), ''.encode("utf-8")]
        elif "pagan" in args[0] and "-v" in args[0]:
            stdout = "This is PAGAN v.0.61."
            self.output = [stdout.encode("utf-8"), ''.encode("utf-8")]
        elif "pagan" in args[0]:
            _file = join(RES_PATH, "mock_resources", "test_pagan", "stdout.txt")
            with open(_file, "r", encoding="utf-8") as ifile:
                stdout = ifile.read()
            self.output = [stdout, ""]
        elif "prank" in args[0] and "-help" in args[0]:
            stdout = "prank v.140603. Minimal usage: 'prank sequence_file'"
            self.output = [stdout.encode("utf-8"), ''.encode("utf-8")]
        elif "prank" in args[0]:
            self.output = [None, None]
        elif "-h" in args[0] or "-v" in args[0]:
            self.output = ["Nothing".encode("utf-8"), "here".encode("utf-8")]

    def communicate(self):
        return self.output


def test_clustalomega(sb_resources, hf, monkeypatch):
    mock_tmp_dir = br.TempDir()
    tmp_dir = br.TempDir()
    shutil.copy(join(RES_PATH, "mock_resources", "test_clustalo", "result"), join(mock_tmp_dir.path, "result"))
    monkeypatch.setattr(Alb, "which", lambda *_: True)
    monkeypatch.setattr(Alb, "Popen", MockPopen)
    monkeypatch.setattr(br, "TempDir", lambda: mock_tmp_dir)

    # basic
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalomega')
    assert hf.buddy2hash(tester) == "f5afdc7c76ab822bdc95230329766aba", tester.write("temp.del")

    # quiet
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalomega', quiet=True)
    assert hf.buddy2hash(tester) == "f5afdc7c76ab822bdc95230329766aba", tester.write("temp.del")

    # params
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalomega', "--outfmt=nexus", quiet=True)
    assert hf.buddy2hash(tester) == "23d7c9fa33454ed551a5896e532cf552", tester.write("temp.del")

    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalomega', "--outfmt=foobar", quiet=True)
    assert hf.buddy2hash(tester) == "f5afdc7c76ab822bdc95230329766aba", tester.write("temp.del")

    # keep
    monkeypatch.setattr(Sb, "hash_ids", mock_hash_ids)
    tester = sb_resources.get_one("d f")
    Alb.generate_msa(tester, 'clustalomega', keep_temp=join(tmp_dir.path, "keep_files"))
    root, dirs, files = next(os.walk(join(tmp_dir.path, "keep_files")))
    assert "result" in files
    assert "tmp.fa" in files


def test_clustalw2(sb_resources, hf, monkeypatch):
    mock_tmp_dir = br.TempDir()
    tmp_dir = br.TempDir()
    shutil.copy(join(RES_PATH, "mock_resources", "test_clustalw2", "result"), join(mock_tmp_dir.path, "result"))
    monkeypatch.setattr(Alb, "which", lambda *_: True)
    monkeypatch.setattr(Alb, "Popen", MockPopen)
    monkeypatch.setattr(br, "TempDir", lambda: mock_tmp_dir)

    # basic
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalw2')
    assert hf.buddy2hash(tester) == "955440b5139c8e6d7d3843b7acab8446", tester.write("temp.del")

    # quiet
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalw2', quiet=True)
    assert hf.buddy2hash(tester) == "955440b5139c8e6d7d3843b7acab8446", tester.write("temp.del")

    # params
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalw2', "-output=nexus", quiet=True)
    assert hf.buddy2hash(tester) == "f4a61a8c2d08a1d84a736231a4035e2e", tester.write("temp.del")

    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalw2', "-output=foobar", quiet=True)
    assert hf.buddy2hash(tester) == "955440b5139c8e6d7d3843b7acab8446", tester.write("temp.del")

    # keep
    monkeypatch.setattr(Sb, "hash_ids", mock_hash_ids)
    tester = sb_resources.get_one("d f")
    Alb.generate_msa(tester, 'clustalw2', keep_temp=join(tmp_dir.path, "keep_files"))
    root, dirs, files = next(os.walk(join(tmp_dir.path, "keep_files")))
    assert "result" in files
    assert "tmp.fa" in files


def test_pagan(sb_resources, hf, monkeypatch):
    mock_tmp_dir = br.TempDir()
    tmp_dir = br.TempDir()
    shutil.copy(join(RES_PATH, "mock_resources", "test_pagan", "result.fas"), join(mock_tmp_dir.path, "result.fas"))
    open("warnings", "w").close()
    assert os.path.isfile("warnings")
    monkeypatch.setattr(Alb, "which", lambda *_: True)
    monkeypatch.setattr(Alb, "Popen", MockPopen)
    monkeypatch.setattr(br, "TempDir", lambda: mock_tmp_dir)

    # basic
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan')
    assert hf.buddy2hash(tester) == "da1c6bb365e2da8cb4e7fad32d7dafdb"
    assert not os.path.isfile("warnings")

    # quiet
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', quiet=True)
    assert hf.buddy2hash(tester) == "da1c6bb365e2da8cb4e7fad32d7dafdb"

    # params
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', "-f nexus", quiet=True)
    assert hf.buddy2hash(tester) == "f93607e234441a2577fa7d8a387ef7ec", tester.write("temp.del")

    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', "-f foobar", quiet=True)
    assert hf.buddy2hash(tester) == "da1c6bb365e2da8cb4e7fad32d7dafdb", tester.write("temp.del")

    # keep
    monkeypatch.setattr(Sb, "hash_ids", mock_hash_ids)
    tester = sb_resources.get_one("d f")
    Alb.generate_msa(tester, 'pagan', keep_temp=join(tmp_dir.path, "keep_files"))
    root, dirs, files = next(os.walk(join(tmp_dir.path, "keep_files")))
    assert "result.fas" in files
    assert "tmp.fa" in files


def test_prank(sb_resources, hf, monkeypatch):
    mock_tmp_dir = br.TempDir()
    tmp_dir = br.TempDir()
    shutil.copy(join(RES_PATH, "mock_resources", "test_prank", "result.best.fas"),
                join(mock_tmp_dir.path, "result.best.fas"))
    monkeypatch.setattr(Alb, "which", lambda *_: True)
    monkeypatch.setattr(Alb, "Popen", MockPopen)
    monkeypatch.setattr(br, "TempDir", lambda: mock_tmp_dir)

    # basic
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'prank')
    assert hf.buddy2hash(tester) == "eff3e6728b5126e285a422863567294f", tester.write("temp.del")

    # quiet
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'prank', quiet=True)
    assert hf.buddy2hash(tester) == "eff3e6728b5126e285a422863567294f", tester.write("temp.del")

    # params
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'prank', "-f=nexus", quiet=True)
    assert hf.buddy2hash(tester) == "4dcaa948e109487ee12512b6ac02183c", tester.write("temp.del")

    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'prank', "-f=foobar", quiet=True)
    assert hf.buddy2hash(tester) == "eff3e6728b5126e285a422863567294f", tester.write("temp.del")

    # keep
    monkeypatch.setattr(Sb, "hash_ids", mock_hash_ids)
    tester = sb_resources.get_one("d f")
    Alb.generate_msa(tester, 'prank', keep_temp=join(tmp_dir.path, "keep_files"))
    root, dirs, files = next(os.walk(join(tmp_dir.path, "keep_files")))
    assert "result.best.fas" in files
    assert "tmp.fa" in files


def test_muscle(sb_resources, hf, monkeypatch):
    tmp_dir = br.TempDir()
    monkeypatch.setattr(Alb, "which", lambda *_: True)
    monkeypatch.setattr(Alb, "Popen", MockPopen)

    # basic
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'muscle')
    assert hf.buddy2hash(tester) == "5ec18f3e0c9f5cf96944a1abb130232f", tester.write("temp.del")

    # quiet
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'muscle', quiet=True)
    assert hf.buddy2hash(tester) == "5ec18f3e0c9f5cf96944a1abb130232f", tester.write("temp.del")

    # keep
    monkeypatch.setattr(Sb, "hash_ids", mock_hash_ids)
    tester = sb_resources.get_one("d f")
    Alb.generate_msa(tester, 'muscle', keep_temp=join(tmp_dir.path, "keep_files"))
    root, dirs, files = next(os.walk(join(tmp_dir.path, "keep_files")))
    assert files == ["tmp.fa"]
    with open(join(root, "tmp.fa"), "r", encoding="utf-8") as ifile:
        kept_file = ifile.read()
    assert hf.string2hash(kept_file) == "b831e901d8b6b1ba52bad797bad92d14", print(kept_file)


def test_mafft(sb_resources, hf, monkeypatch):
    tmp_dir = br.TempDir()
    tmp_dir.subdir("keep_files")
    monkeypatch.setattr(Alb, "which", lambda *_: True)
    monkeypatch.setattr(Alb, "Popen", MockPopen)

    # basic
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'mafft')
    assert hf.buddy2hash(tester) == "f94e0fd591dad83bd94201f0af038904", tester.write("temp.del")

    # quiet
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'mafft', quiet=True)
    assert hf.buddy2hash(tester) == "f94e0fd591dad83bd94201f0af038904", tester.write("temp.del")

    # keep
    monkeypatch.setattr(Sb, "hash_ids", mock_hash_ids)
    monkeypatch.setattr(br, "ask", lambda *_: True)
    tester = sb_resources.get_one("d f")
    Alb.generate_msa(tester, 'mafft', keep_temp=join(tmp_dir.path, "keep_files"))
    root, dirs, files = next(os.walk(join(tmp_dir.path, "keep_files")))
    assert files == ["tmp.fa"]
    with open(join(root, "tmp.fa"), "r", encoding="utf-8") as ifile:
        kept_file = ifile.read()
    assert hf.string2hash(kept_file) == "b831e901d8b6b1ba52bad797bad92d14", print(kept_file)


def test_alignment_edges(monkeypatch, sb_resources):
    mock_tmp_dir = br.TempDir()
    tmp_dir = br.TempDir()
    shutil.copy(join(RES_PATH, "mock_resources", "test_muscle", "result"), join(mock_tmp_dir.path, "result"))
    monkeypatch.setattr(Alb, "which", lambda *_: True)
    monkeypatch.setattr(br, "Popen", MockPopen)
    monkeypatch.setattr(Alb, "Popen", MockPopen)
    monkeypatch.setattr(br, "TempDir", lambda: mock_tmp_dir)

    # Weird binary given, but muscle found
    tester = sb_resources.get_one("d f")
    with pytest.raises(br.GuessError) as err:
        tester = Alb.generate_msa(tester, "my_mucsle")
    assert "Could not determine format from raw input" in str(err)

    with pytest.raises(AttributeError) as err:
        Alb.generate_msa(tester, "foo")
    assert "foo is not a recognized alignment tool. Please check your spelling (case sensitive)" in str(err)

    monkeypatch.setattr(br, "ask", lambda *_: False)
    tmp_dir.subdir("keep_files")
    with pytest.raises(SystemExit):
        Alb.generate_msa(tester, "mafft", keep_temp=join(tmp_dir.path, "keep_files"))


# ######################  '-hi', '--hash_ids' ###################### #
def test_hash_seq_ids(alb_resources):
    tester = alb_resources.get_one("o p g")
    Alb.hash_ids(tester)
    assert len(tester.records()[0].id) == 10

    tester = alb_resources.get_one("m d pr")
    Alb.hash_ids(tester, 25)
    assert len(tester.records()[0].id) == 25
    assert len(tester.hash_map) == 34


def test_hash_seq_ids_errors(alb_resources):
    tester = alb_resources.get_one("o d f")

    with pytest.raises(TypeError) as e:
        Alb.hash_ids(tester, "foo")
    assert str(e.value) == "Hash length argument must be an integer, not <class 'str'>"

    with pytest.raises(ValueError) as e:
        Alb.hash_ids(tester, 0)
    assert str(e.value) == "Hash length must be greater than 0"

    tester.alignments *= 10
    with pytest.raises(ValueError) as e:
        Alb.hash_ids(tester, 1)
    assert "Insufficient number of hashes available to cover all sequences." in str(e.value)


# #################################### 'lc', '--lowercase' and 'uc', '--uppercase' ################################### #
hashes = [('o d g', '842d9c911a33c0fd0484383eabefb0fe', '842d9c911a33c0fd0484383eabefb0fe'),
          ('o d n', '52e74a09c305d031fc5263d1751e265d', 'cb1169c2dd357771a97a02ae2160935d'),
          ('o d py', 'cfe6cb9c80aebd353cf0378a6d284239', '503e23720beea201f8fadf5dabda75e4'),
          ('o d s', 'b82538a4630810c004dc8a4c2d5165ce', '228e36a30e8433e4ee2cd78c3290fa6b'),
          ('o p g', '46388b175b31b81f47199ae6327768af', '46388b175b31b81f47199ae6327768af'),
          ('o p n', '8b6737fe33058121fd99d2deee2f9a76', '17ff1b919cac899c5f918ce8d71904f6'),
          ('o p py', '968ed9fa772e65750f201000d7da670f', 'aacda2f5d4077f23926400f74afa2f46'),
          ('o p s', 'f35cbc6e929c51481e4ec31e95671638', 'c0dce60745515b31a27de1f919083fe9'),
          ('m d py', '6259e675def07bd4977f4ab1f5ffc26d', '0974ac9aefb2fb540957f15c4869c242'),
          ('m d s', 'f3f7b66ef034d3807683a2d5a0c44cad', 'a217b9f6000f9eeff98faeb9fd09efe4'),
          ('m p py', '2a77f5761d4f51b88cb86b079e564e3b', 'd13551548c9c1e966d0519755a8fb4eb'),
          ('m p s', '6f3f234d796520c521cb85c66a3e239a', '00661f7afb419c6bb8c9ac654af7c976')]


@pytest.mark.parametrize("key,uc_hash,lc_hash", hashes)
def test_cases(key, uc_hash, lc_hash, alb_resources, hf):
    tester = Alb.uppercase(alb_resources.get_one(key))
    assert hf.buddy2hash(tester) == uc_hash, tester.write(join("error_files", uc_hash))
    tester = Alb.lowercase(tester)
    assert hf.buddy2hash(tester) == lc_hash, tester.write(join("error_files", lc_hash))


# ##################### '-mf2a', '--map_features2alignment' ###################### ##
hashes = [('o p n', '46f48b266b424b5e43d11d589812f1bf'), ('o p pr', 'c377cbacc34121e9460ec048f4c72e1b'),
          ('o p psr', 'c377cbacc34121e9460ec048f4c72e1b'), ('o p s', '125179d2965d57d4beac5cd91c5fd218'),
          ('o d n', 'c1359470ad0916902c4e96facd088378'), ('o d pr', 'a61222fad78fd4ddee6108cda8d4bdc7'),
          ('o d psr', 'a61222fad78fd4ddee6108cda8d4bdc7'), ('o d s', 'de9ed2e4c71c8882ee527e1d0b405016')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_map_features2alignment(key, next_hash, alb_resources, hf):
    alignbuddy = alb_resources.get_one(key)
    if alignbuddy.alpha == IUPAC.protein:
        seqbuddy = Sb.SeqBuddy(join(RES_PATH, "Mnemiopsis_pep.gb"))
    else:
        seqbuddy = Sb.SeqBuddy(join(RES_PATH, "Mnemiopsis_cds.gb"))
    tester = Alb.map_features2alignment(seqbuddy, alignbuddy)
    tester.set_format("genbank")
    assert hf.buddy2hash(tester) == next_hash, tester.write(join("error_files", next_hash))


# ###########################################  '-pi', '--percent_id' ############################################ #
def test_percent_id(alb_resources):
    alignbuddy = alb_resources.get_one("o p s")
    alignbuddy = Alb.percent_id(alignbuddy)

    assert alignbuddy.alignments[0].percent_ids
    assert alignbuddy.alignments[0].percent_ids["Mle-Panxα12"]["Mle-Panxα8"] == 0.47877358490566035

    alignbuddy = alb_resources.get_one("o d s")
    Alb.percent_id(alignbuddy)
    assert alignbuddy.alignments[0].percent_ids
    assert alignbuddy.alignments[0].percent_ids["Mle-Panxα9"]["Mle-Panxα11"] == 0.5563725490196079


# ###########################################  '-pfm', '--pos_freq_mat' ############################################ #
def test_position_frequency_matrix(alb_resources):
    for alignbuddy in alb_resources.get_list("o p d r n"):
        Alb.position_frequency_matrix(alignbuddy)
        for align, length in zip(alignbuddy.alignments, alignbuddy.lengths()):
            assert len(align.pfm) == length
            if alignbuddy.alpha == IUPAC.ambiguous_dna:
                assert align.pfm[-1]["-"] == 0.769231
            if alignbuddy.alpha == IUPAC.ambiguous_rna:
                assert align.pfm[-1]["-"] == 0.769231
            if alignbuddy.alpha == IUPAC.protein:
                assert align.pfm[-1]["-"] == 0.692308


# ###########################################  '-oi', '--order_ids' ############################################ #
hashes = [('o d g', '8f1846922f3c4d955c42964ba0c24649', '982e66fa5eeba8de5c570a770042ec10'),
          ('o d n', '132757da01b3caf174d024efdb2c3acd', '286bac7a213997924203622c3357457c'),
          ('o d py', '3c49bdc1b0fe4e1d6bfc148eb0293e21', 'd6e79a5faeaff396aa7eab0b460c3eb9'),
          ('o p g', '80cc4fc343879ce5eb9c0ec0f430d965', '3a00afabf18266b813a99ebc34140d18'),
          ('o p n', '197ba12e799ab2a1dadfe1b254381e00', 'f32fabc627615d25d8cd57553e7281af'),
          ('o p py', 'ffae954adc0d362354e43c1b70d9be29', 'f4c0924087fdb624823d02e909d94e95'),
          ('m d py', 'a44938e26e4b35967ed8e17a0eaebe4c', '9d6b6087d07f7d1fd701591ab7cb576d'),
          ('m p py', '5bdda310b29b18057e056f3c982446b2', '439f57b891dd2a72724b10c124f96378')]


@pytest.mark.parametrize("key,fwd_hash,rev_hash", hashes)
def test_order_ids1(key, fwd_hash, rev_hash, alb_resources, hf):
    alignbuddy = alb_resources.get_one(key)
    Alb.order_ids(alignbuddy)
    assert hf.buddy2hash(alignbuddy) == fwd_hash, alignbuddy.write(join("error_files", fwd_hash))

    Alb.order_ids(alignbuddy, reverse=True)
    assert hf.buddy2hash(alignbuddy) == rev_hash, alignbuddy.write(join("error_files", rev_hash))


def test_order_ids2(alb_resources, hf):
    alignbuddy = alb_resources.get_one("o p n")
    Alb.rename(alignbuddy, "Mle-Panxα4", "Mle004-Panxα4")
    Alb.rename(alignbuddy, "Mle-Panxα5", "Mle05-Panxα5")
    Alb.rename(alignbuddy, "Mle-Panxα9", "aMle-PanxαBlahh")
    Alb.order_ids(alignbuddy)
    assert hf.buddy2hash(alignbuddy) == "5c1316e18205432b044101e720646cd5"


# ##################### '-pr', '--pull_records' ###################### ##
hashes = [('o d g', '8488d218201ef84c8fe458576d38e3be'), ('o d n', 'd82e66c57548bcf8cba202b13b070ead'),
          ('o d py', 'd141752c38a892ccca800c637f609608'), ('o p g', '11de462227fb680510227b822f8c375e'),
          ('o p n', '027bbc7e34522f9521f83ee7d03793a1'), ('o p py', '2cd74d7ede4d1fb6e18363567426437e'),
          ('m d py', '7c77c6f3245c21842f4be585714ec6ce'), ('m p py', 'f34fa4c34cfe5c1e6b228949557c9483')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_pull_records(key, next_hash, alb_resources, hf):
    alignbuddy = alb_resources.get_one(key)
    Alb.pull_records(alignbuddy, "α[1-5]$|β[A-M]")
    assert hf.buddy2hash(alignbuddy) == next_hash, alignbuddy.write(join("error_files", next_hash))


# ###########################################  '-ri', '--rename_ids' ############################################ #
hashes = [('o d g', '98f69c2d39c9a4ca0cb5f7da026095cd'), ('o d n', '243024bfd2f686e6a6e0ef65aa963494'),
          ('o d py', '98bb9b57f97555d863054ddb526055b4'), ('o p g', '7a72ab9a2ef49a97ee60862aab1c88c3'),
          ('o p n', '3598e85169ed3bcdbcb676bb2eb6cef0'), ('o p py', 'd49eb4de01d727b9e3ad648d6a04a3c9'),
          ('m d py', 'ddfffd9b999637abf7f5926f017de987'), ('m p py', '0a68229bd13439040f045cd8c72d7cc9')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_rename_ids(key, next_hash, alb_resources, hf):
    alignbuddy = alb_resources.get_one(key)
    Alb.rename(alignbuddy, 'Panx', 'Test', 0)
    assert hf.buddy2hash(alignbuddy) == next_hash, alignbuddy.write(join("error_files", next_hash))


# ###########################################  'tr', '--translate' ############################################ #
hashes = [('o d f', 'b7fe22a87fb78ce747d80e1d73e39c35'), ('o d g', '542794541324d74ff636eaf4ee5e6b1a'),
          ('o d n', 'a2586af672ad71f16bbd54f359b323ff'), ('o d py', 'd0d4dd408e559215b2780f4f0ae0c418'),
          ('o d pr', 'f77705e32cd753267916539ee0936e1f'), ('o d pss', 'ede672b15221ec60981287ca1e286c52'),
          ('o d psr', '623fe1634752e812f482cfa7b7ea20ee'), ('o d s', '4ff563c39229d30aa3eda193cb290344'),
          ('o d c', '150179326629fffadb7aef7796bd1cec'), ('m d py', '7fd28236f491c38ba261dfde20919595'),
          ('m d pr', '0de676236eda864172b73b6abe4d7a05'), ('m d pss', 'c7feff60c16b2b187e03db5d160a4748'),
          ('m d psr', 'b5945f0317fe9ce8fc03ac7f4c0d5932'), ('m d s', 'ee5a41b6f8b32645359beafc72efe825'),
          ('m d c', '09251e1f4fc0e07a5bba4c64c22bac9b')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_translate1(key, next_hash, alb_resources, hf):
    alignbuddy = alb_resources.get_one(key)
    Alb.translate_cds(alignbuddy)
    assert hf.buddy2hash(alignbuddy) == next_hash, alignbuddy.write(join("error_files", next_hash))


def test_translate2(alb_resources):
    # Protein input
    with pytest.raises(TypeError) as e:
        Alb.translate_cds(alb_resources.get_one('o p s'))
    assert "Nucleic acid sequence required, not protein." in str(e)

    tester = alb_resources.get_one('o d s')
    tester.records()[0].seq.alphabet = IUPAC.protein
    with pytest.raises(TypeError) as e:
        Alb.translate_cds(tester)
    assert "Record 'Mle-Panxα9' is protein." in str(e)


# ###########################################  'tm', '--trimal' ############################################ #
hashes = [('o d psr', '5df948e4b2cb6c0d0740984445655135', '384563eb411713e90cb2fea0c799bf0d'),
          ('m d psr', '0e93f0a8c77da8ec974eeca311ca6636', 'b15f333416e9dd44834f468d5cd4ca8d'),
          ('o p psr', 'b87f927511aade73bc795e024af8975e', 'e0f5ce9201249daf4bb3b4f70a7b5ce8'),
          ('m p psr', 'f0f2115e29f6dfcb75036d90b06edab4', 'f443fbe1831fe368a11edc51e25fa330')]


@pytest.mark.parametrize("key,hash3,hash07", hashes)
def test_trimal(key, hash3, hash07, alb_resources, hf):
    alignbuddy = alb_resources.get_one(key)
    tester1, tester2 = Alb.make_copy(alignbuddy), Alb.make_copy(alignbuddy)
    Alb.trimal(tester1, 3)
    assert hf.buddy2hash(tester1) == hash3, alignbuddy.write(join("error_files", hash3))

    tester1, tester2 = Alb.make_copy(alignbuddy), Alb.make_copy(alignbuddy)
    Alb.trimal(tester1, 0.7)
    assert hf.buddy2hash(tester1) == hash07, alignbuddy.write(join("error_files", hash07))


def test_trimal2(alb_resources, hf):
    tester = Alb.trimal(alb_resources.get_one("o p n"), 'all')
    assert hf.buddy2hash(tester) == "8faaf09741ddb3137653cb77ee66974a"
    tester = alb_resources.get_one("o p n")
    tester.alignments[0]._records = tester.alignments[0]._records[:5]
    Alb.trimal(tester, 'clean')
    assert hf.buddy2hash(tester) == "93a2aa21e6baf5ca70eb2de52ae8dbea"
    tester = alb_resources.get_one("o p n")
    tester_dir = TEMPDIR.subdir()
    tester.write(join(tester_dir, "trimal"))
    assert hf.buddy2hash(Alb.trimal(tester, 'gappyout')) == "2877ecfb201fc35211a4625f34c7afdd"
    """ Probably not a good idea to be calling binaries like this...
    real_trimal = Popen("trimal -in %s%strimal -gappyout" % (tester_dir, os.path.sep),
                        stdout=PIPE, shell=True).communicate()
    real_trimal = real_trimal[0].decode()
    with open("%s%strimal" % (tester_dir, os.path.sep), "w") as ofile:
        ofile.write(real_trimal)
    tester = Alb.AlignBuddy("%s%strimal" % (tester_dir, os.path.sep))
    assert hf.buddy2hash(tester) == "2877ecfb201fc35211a4625f34c7afdd"
    """
    records = [SeqRecord(Seq("A--G-")), SeqRecord(Seq("--T--")), SeqRecord(Seq("--TG-")), SeqRecord(Seq("A---C"))]
    tester = Alb.AlignBuddy([MultipleSeqAlignment(records)])
    Alb.trimal(tester, "gappyout")
    assert "".join([str(rec.seq) for rec in tester.records()]) == ""
