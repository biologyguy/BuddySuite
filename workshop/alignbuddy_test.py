#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

import pytest
from hashlib import md5
import os
from io import StringIO
import re
from copy import deepcopy
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import sys

sys.path.insert(0, "./")
import buddy_resources as br
from MyFuncs import TempFile
import AlignBuddy as Alb
import SeqBuddy as Sb


def align_to_hash(_alignbuddy, mode='hash'):
    if mode != "hash":
        return str(_alignbuddy)
    _hash = md5("{0}\n".format(str(_alignbuddy).rstrip()).encode()).hexdigest()
    return _hash

root_dir = os.getcwd()


def resource(file_name):
    return "{0}/unit_test_resources/{1}".format(root_dir, file_name)


def test_guess_format():
    assert Alb._guess_format(["dummy", "list"]) == "stockholm"
    assert Alb._guess_format(alb_objects[0]) == "nexus"

    with open(resource("Alignments_pep.stklm"), "r") as ifile:
        assert Alb._guess_format(ifile) == "stockholm"
        ifile.seek(0)
        string_io = StringIO(ifile.read())
    assert Alb._guess_format(string_io) == "stockholm"
    with pytest.raises(Alb.GuessError):
        Alb._guess_format({"Dummy dict": "Type not recognized by guess_format()"})


align_files = ["Mnemiopsis_cds.nex", "Mnemiopsis_cds.phy", "Mnemiopsis_cds.phyr", "Mnemiopsis_cds.stklm",
               "Mnemiopsis_pep.nex", "Mnemiopsis_pep.phy", "Mnemiopsis_pep.phyr", "Mnemiopsis_pep.stklm",
               "Alignments_pep.phy", "Alignments_pep.phyr", "Alignments_pep.stklm",
               "Alignments_cds.phyr", "Alignments_cds.stklm"]

file_types = ["nexus", "phylip", "phylip-relaxed", "stockholm",
              "nexus", "phylip", "phylip-relaxed", "stockholm",
              "phylip", "phylip-relaxed", "stockholm",
              "phylip-relaxed", "stockholm"]

nucl_indices = [0, 1, 2, 3, 11, 12]

input_tuples = [(next_file, file_types[indx]) for indx, next_file in enumerate(align_files)]


@pytest.mark.parametrize("align_file,file_type", input_tuples)
def test_instantiate_alignbuddy_from_file(align_file, file_type):
    assert type(Alb.AlignBuddy(resource(align_file), _in_format=file_type)) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_file_guess(align_file):
    assert type(Alb.AlignBuddy(resource(align_file))) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_handle(align_file):
    with open(resource(align_file), 'r') as ifile:
        assert type(Alb.AlignBuddy(ifile)) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_raw(align_file):
    with open(resource(align_file), 'r') as ifile:
        assert type(Alb.AlignBuddy(ifile.read())) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_alignbuddy(align_file):
    tester = Alb.AlignBuddy(resource(align_file))
    assert type(Alb.AlignBuddy(tester)) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", align_files)
def test_instantiate_alignbuddy_from_list(align_file):
    tester = Alb.AlignBuddy(resource(align_file))
    assert type(Alb.AlignBuddy(tester.alignments)) == Alb.AlignBuddy

    with pytest.raises(TypeError):  # When non-MultipleSeqAlignment objects are in the .alignments list
        tester.alignments.append("Dummy string object")
        Alb.AlignBuddy(tester.alignments)


def test_empty_file():
    with open(resource("blank.fa"), "r") as ifile:
        with pytest.raises(SystemExit):
            Alb.AlignBuddy(ifile)


# Now that we know that all the files are being turned into AlignBuddy objects okay, make them all objects so it doesn't
# need to be done over and over for each subsequent test.
alb_objects = [Alb.AlignBuddy(resource(x)) for x in align_files]


# ################################################# HELPER FUNCTIONS ################################################# #
def test_get_seq_recs():
    tester = str(Alb._get_seq_recs(alb_objects[8]))
    tester = md5(tester.encode()).hexdigest()
    assert tester == "6168f8b57d0ff78d70fd22ee09d713b5"


def test_guess_alphabet():
    assert str(type(Alb._guess_alphabet(alb_objects[0]))) == "<class 'Bio.Alphabet.IUPAC.IUPACAmbiguousDNA'>"
    assert str(type(Alb._guess_alphabet(alb_objects[4]))) == "<class 'Bio.Alphabet.IUPAC.IUPACProtein'>"
    tester = Alb.AlignBuddy(resource("Mnemiopsis_rna.nex"))
    assert str(type(Alb._guess_alphabet(tester))) == "<class 'Bio.Alphabet.IUPAC.IUPACAmbiguousRNA'>"
    assert not Alb._guess_alphabet(Alb.AlignBuddy("", _in_format="fasta"))


def test_guess_error():
    # File path
    with pytest.raises(Alb.GuessError):
        Alb.AlignBuddy(resource("unrecognizable.txt"))

    with open(resource("unrecognizable.txt"), 'r') as ifile:
        # Raw
        with pytest.raises(Alb.GuessError):
            Alb.AlignBuddy(ifile.read())

        # Handle
        with pytest.raises(Alb.GuessError):
            ifile.seek(0)
            Alb.AlignBuddy(ifile)

    # GuessError output
    try:
        Alb.AlignBuddy(resource("unrecognizable.txt"))
    except Alb.GuessError as e:
        assert "Could not determine format from _input file" in str(e.value) and \
               "\nTry explicitly setting with -f flag." in str(e.value)


def test_phylip():
    tester = Alb.AlignBuddy(resource('Mnemiopsis_cds_hashed_ids.nex'), _out_format='phylip-strict')
    tester = "{0}\n".format(str(tester).rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "16b3397d6315786e8ad8b66e0d9c798f"

    tester = Alb.AlignBuddy(resource('Mnemiopsis_cds.nex'), _out_format='phylip-relaxed')
    tester = "{0}\n".format(str(tester).rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "52c23bd793c9761b7c0f897d3d757c12"

    tester = Alb.AlignBuddy(resource('Mnemiopsis_cds_hashed_ids.nex'), _out_format='phylip-sequential-strict')
    tester = "{0}\n".format(str(tester).rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "cc9278c0dbe6b315162545115978dda3"

    tester = Alb.AlignBuddy(resource('Mnemiopsis_cds.nex'), _out_format='phylip-sequential')
    tester = "{0}\n".format(str(tester).rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "90267e3a80f29bc966dc580e26c1fc0a"

    tester = Alb.AlignBuddy(resource('Mnemiopsis_cds.nex'), _out_format='phylip-sequential-strict')
    with pytest.raises(ValueError):
        str(tester)

    tester = Alb.AlignBuddy(resource('Mnemiopsis_cds.nex'), _out_format='phylip-interleaved-strict')
    with pytest.raises(ValueError):
        str(tester)


# AlignBuddy print() and __str__() methods
hashes = ["cb1169c2dd357771a97a02ae2160935d", "f59e28493949f78637691caeb617ab50",
          "52c23bd793c9761b7c0f897d3d757c12", "228e36a30e8433e4ee2cd78c3290fa6b",
          "17ff1b919cac899c5f918ce8d71904f6", "5af1cf061f003d3351417458c0d23811",
          "f3e98897f1bbb3d3d6c5507afa9f814e", "c0dce60745515b31a27de1f919083fe9",
          "90578980479ad235338dbb767444b05b", "9c6773e7d24000f8b72dd9d25620cff1",
          "3fd5805f61777f7f329767c5f0fb7467"]
hashes = [(alb_objects[indx], value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_print(alignbuddy, next_hash, capsys):
    alignbuddy.print()
    out, err = capsys.readouterr()
    out = "{0}\n".format(out.rstrip())
    tester = md5(out.encode()).hexdigest()
    assert tester == next_hash


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_str(alignbuddy, next_hash):
    tester = str(alignbuddy)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash


def test_stderr(capsys):
    Alb._stderr("Hello std_err", quiet=False)
    out, err = capsys.readouterr()
    assert err == "Hello std_err"

    Alb._stderr("Hello std_err", quiet=True)
    out, err = capsys.readouterr()
    assert err == ""


def test_stdout(capsys):
    Alb._stdout("Hello std_out", quiet=False)
    out, err = capsys.readouterr()
    assert out == "Hello std_out"

    Alb._stdout("Hello std_out", quiet=True)
    out, err = capsys.readouterr()
    assert out == ""


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_write1(alignbuddy, next_hash):
    temp_file = TempFile()
    alignbuddy.write(temp_file.path)
    out = "{0}\n".format(temp_file.read().rstrip())
    tester = md5(out.encode()).hexdigest()
    assert tester == next_hash


def test_write2():  # Unloopable components
    tester = Alb._make_copies(alb_objects[8])
    tester.out_format = "fasta"
    with pytest.raises(ValueError):
        str(tester)

    tester.alignments = []
    assert str(tester) == "AlignBuddy object contains no alignments.\n"

    tester = Alb._make_copies(alb_objects[2])
    tester.out_format = "phylipi"
    assert md5(str(tester).encode()).hexdigest() == "52c23bd793c9761b7c0f897d3d757c12"

    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds_hashed_ids.nex"))
    tester.out_format = "phylip-strict"
    assert md5(str(tester).encode()).hexdigest() == "16b3397d6315786e8ad8b66e0d9c798f"


# ################################################ MAIN API FUNCTIONS ################################################ #
# ############################################  'al', '--alignment_lengths' ############################################## #
def test_alignment_lengths():
    tester_100 = Alb.alignment_lengths(Alb.AlignBuddy(resource('concat_alignment_file.phyr')))
    assert tester_100[0] == 100
    assert tester_100[1] == 100

    tester_variable = Alb.alignment_lengths(Alb.AlignBuddy(resource('Alignments_cds.phyr')))
    assert tester_variable[0] == 2043
    assert tester_variable[1] == 1440


# ##############################################  'cs', '--clean_seqs' ############################################### #
def test_clean_seqs():
    # Test an amino acid file
    tester = Alb._make_copies(alb_objects[8])
    Alb.clean_seq(tester)

    with pytest.raises(ValueError):
        str(tester)

    tester.out_format = "fasta"
    with pytest.raises(ValueError):
        str(tester)

    tester.out_format = "genbank"
    assert align_to_hash(tester) == "b613060f43f66248cda0088c06de0949"

    # Test nucleotide files, but skip the errors
    tester = Alb._make_copies(alb_objects[11])
    Alb.clean_seq(tester, skip_list="RYWSMKHBVDNXrywsmkhbvdnx")
    tester.out_format = "genbank"
    assert align_to_hash(tester) == "81189bf7962cc2664815dc7fca8cd95d"

    Alb.clean_seq(tester)
    assert align_to_hash(tester) == "2f2cbf227d45aa49d15971ce214d3191"

# ###########################################  'ca', '--codon_alignment' ############################################ #
hashes = ["c907d29434fe2b45db60f1a9b70f110d", "f150b94234e93f354839d7c2ca8dae24", "6a7a5416f2ce1b3161c8b5b8b4b9e901",
          "54d412fbca5baa60e4c31305d35dd79a", "3ddc2109b15655ef0eed3908713510de", "f728ab606602ed67357f78194e500664"]

hashes = [(Alb._make_copies(alb_objects[nucl_indices[indx]]), next_hash) for indx, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_codon_alignment1(alignbuddy, next_hash):
    tester = Alb.codon_alignment(alignbuddy)
    assert align_to_hash(tester) == next_hash

    if next_hash == "f728ab606602ed67357f78194e500664":
        with pytest.raises(TypeError):
            Alb.codon_alignment(alb_objects[8])


def test_codon_alignment2():
    with pytest.raises(TypeError):
        tester = Alb._make_copies(alb_objects[3])
        tester.alignments[0][0].seq = Seq("MLDILSKFKGVTPFKGITIDDGWDQLNRSFMFVLLVVMGTTVTVRQYTGSVISCDGFKKFGSTFAEDYCWTQGLY",
                                          alphabet=IUPAC.protein)
        Alb.translate_cds(tester)


# ###########################################  'cta', '--concat_alignments' ######################################### #
def test_concat_alignments_identical():
    tester = Alb.AlignBuddy(resource("duplicate_alignment.nex"))
    tester.alignments.append(tester.alignments[0])
    Alb.concat_alignments(tester, '-')
    assert align_to_hash(tester) == '1b656db96d33973fe6b1368afc974148'


def test_concat_alignments_duplicate_taxa():
    tester = Alb.AlignBuddy(resource("concat_alignment_file.phyr"))
    Alb.concat_alignments(tester, '-')
    assert align_to_hash(tester) == 'd961ac6293597c91f495936fc0ea85f9'


# ###########################################  'dr', '--delete_rows' ############################################ #
hashes = ['23d3e0b42b6d63bcb810a5abb0a06ad7', '3458c1f790ff90f8403476af021e97e4', '3c09226535f46299ffd1132c9dd336d8',
          '7455b360c4f1155b1aa0ba7e1219485f', 'b3bf1be13a0b75374905f89aba0302c9', '819c954cc80aaf622734a61c09d301f4',
          'f2492c76baa277ba885e5232fa611fea', '67d81be82ffbcffe2bf0a6f78e361cc9', 'b34588123a2b0b56afdf5d719c61f677']
hashes = [(Alb._make_copies(alb_objects[x]), value) for x, value in enumerate(hashes)]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_delete_rows(alignbuddy, next_hash):
    Alb.delete_rows(alignbuddy, 'Mle-Panxα[567]')
    assert align_to_hash(alignbuddy) == next_hash

# ###########################################  'd2r', '--transcribe' ############################################ #
d2r_hashes = ['e531dc31f24192f90aa1f4b6195185b0', 'b34e4d1dcf0a3a36d36f2be630934d29',
              'a083e03b4e0242fa3c23afa80424d670', '45b511f34653e3b984e412182edee3ca']
d2r_hashes = [(Alb._make_copies(alb_objects[x]), value) for x, value in enumerate(d2r_hashes)]


@pytest.mark.parametrize("alignbuddy, next_hash", d2r_hashes)
def test_transcribe(alignbuddy, next_hash):
    assert align_to_hash(Alb.dna2rna(alignbuddy)) == next_hash


def test_transcribe_pep_exception():
    with pytest.raises(ValueError):
        Alb.dna2rna(deepcopy(alb_objects[4]))


# ###########################################  'er', '--extract_range' ############################################ #
er_hashes = ['10ca718b74f3b137c083a766cb737f31', '637582d17b5701896d80f0380fa00c12', '1bf71754ad8b1e157a9541e27ddd72e6',
             '8873381a66add6c680f54e379ef98c95', '5f400edc6f0990c0cd6eb52ae7687e39', '240a56a273b5049901177284a9240ac3',
             'f0817e10ba740992409193f1bc6600b2', '29a1d24a36d0f26b17aab7faa5b9ad9b', '06fde20b8bcb8bbe6ce5217324096911',
             '146d0550bb24f4f19f31c294c48c7e49', '64a848cc73d909ca6424ce049e2700cd', '9d9bbdf6f20274f26be802fa2361d459',
             '78ceb2e63d8c9d8800b5fa9335a87a30']
er_hashes = [(Alb._make_copies(alb_objects[x]), value) for x, value in enumerate(er_hashes)]


@pytest.mark.parametrize("alignbuddy, next_hash", er_hashes)
def test_extract_range(alignbuddy, next_hash):
    tester = Alb.extract_range(alignbuddy, 0, 50)
    assert align_to_hash(tester) == next_hash


# ###########################################  'ga', '--generate_alignment' ########################################## #
# This is tested for PAGAN version 0.61
@pytest.mark.generate_alignments
def test_pagan_inputs():
    # FASTA
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'pagan')
    assert align_to_hash(tester) == 'da1c6bb365e2da8cb4e7fad32d7dafdb'


@pytest.mark.generate_alignments
def test_pagan_outputs():
    # NEXUS
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'pagan', '-f nexus')
    assert align_to_hash(tester) == 'f93607e234441a2577fa7d8a387ef7ec'
    # PHYLIPI
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'pagan', '-f phylipi')
    assert align_to_hash(tester) == '09dd492fde598670d7cfee61d4e2eab8'
    # PHYLIPS
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'pagan', '-f phylips')
    assert align_to_hash(tester) == 'f079eddd44ffbe038e1418ab03ff7e64'


@pytest.mark.generate_alignments
def test_pagan_multi_param():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'pagan', '-f nexus --translate')
    assert align_to_hash(tester) == 'dd140ec4eb895ce75d574498a58aa28a'


# PRANK is not deterministic, so just test that something reasonable is returned
@pytest.mark.generate_alignments
def test_prank_inputs():
    # FASTA
    tester = Sb.pull_recs(Sb.SeqBuddy(resource("Mnemiopsis_cds.fa")), 'α1')
    tester = Alb.generate_msa(tester, 'prank', '-once')
    assert tester.out_format == 'fasta'


@pytest.mark.generate_alignments
def test_prank_outputs():
    # NEXUS
    tester = Sb.pull_recs(Sb.SeqBuddy(resource("Mnemiopsis_cds.fa")), 'α1')
    tester = Alb.generate_msa(tester, 'prank', '-f=nexus -once')
    assert tester.out_format == 'nexus'
    # PHYLIPI
    tester = Sb.pull_recs(Sb.SeqBuddy(resource("Mnemiopsis_cds.fa")), 'α1')
    tester = Alb.generate_msa(tester, 'prank', '-f=phylipi -once')
    assert tester.out_format == 'phylip-relaxed'
    # PHYLIPS
    tester = Sb.pull_recs(Sb.SeqBuddy(resource("Mnemiopsis_cds.fa")), 'α1')
    tester = Alb.generate_msa(tester, 'prank', '-f=phylips -once')
    assert tester.out_format == 'phylip-sequential'


@pytest.mark.generate_alignments
def test_muscle_inputs():
    # FASTA
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'muscle')
    assert align_to_hash(tester) == '60ada1630165a40be9d5700cc228b1e1'


@pytest.mark.generate_alignments
def test_muscle_outputs():
    # FASTA
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'muscle', '-clw')
    assert align_to_hash(tester) == 'ff8d81f75dfd6249ba1e91e5bbc8bdce'


@pytest.mark.generate_alignments
def test_muscle_multi_param():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'muscle', '-clw -diags')
    assert align_to_hash(tester) == 'ff8d81f75dfd6249ba1e91e5bbc8bdce'


@pytest.mark.generate_alignments
def test_clustalw2_inputs():
    # FASTA
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalw2')
    assert align_to_hash(tester) == 'd744b9cadf592a6d4e8d5eefef90e7c7'


@pytest.mark.generate_alignments
def test_clustalw2_outputs():
    # NEXUS
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalw2', '-output=nexus')
    assert align_to_hash(tester) == 'f4a61a8c2d08a1d84a736231a4035e2e'
    # PHYLIP
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalw2', '-output=phylip')
    assert align_to_hash(tester) == 'a9490f124039c6a2a6193d27d3d01205'
    # FASTA
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalw2', '-output=fasta')
    assert align_to_hash(tester) == '955440b5139c8e6d7d3843b7acab8446'


@pytest.mark.generate_alignments
def test_clustalw2_multi_param():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalw2', '-output=phylip -noweights')
    assert align_to_hash(tester) == 'ae9126eb8c482a82d4060d175803c478'


@pytest.mark.generate_alignments
def test_clustalomega_inputs():
    # FASTA
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalomega')
    assert align_to_hash(tester) == 'c041c78d3d3a62a027490a139ad435e4'
    # PHYLIP
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.phy"))
    tester = Alb.generate_msa(tester, 'clustalomega')
    assert align_to_hash(tester) == '734e93bac16fd2fe49a3340086bde048'
    # STOCKHOLM
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.stklm"))
    tester = Alb.generate_msa(tester, 'clustalomega')
    assert align_to_hash(tester) == '5c7a21e173f8bf54a26ed9d49764bf80'


@pytest.mark.generate_alignments
def test_clustalomega_outputs():
    # CLUSTAL
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=clustal')
    assert align_to_hash(tester) == 'ce25de1a84cc7bfbcd946c88b65cf3e8'
    # PHYLIP
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=phylip')
    assert align_to_hash(tester) == '692c6af848bd90966f15908903894dbd'
    # STOCKHOLM
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=stockholm')
    assert align_to_hash(tester) == '47cc879b68719b8de0eb031d2f0e9fcc'


@pytest.mark.generate_alignments
def test_clustalomega_multi_param():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=clustal --iter=1')
    assert align_to_hash(tester) == '294d8c0260eb81d2039ce8be7289dfcc'


@pytest.mark.generate_alignments
def test_mafft_inputs():
    # FASTA
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'mafft')
    assert align_to_hash(tester) == '8dda0524aaffb326aff09143a1df8a45'


@pytest.mark.generate_alignments
def test_mafft_outputs():
    # CLUSTAL
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'mafft', '--clustalout')
    assert align_to_hash(tester) == '2b8bf89e7459fe9d0b1f29628df6307e'


@pytest.mark.generate_alignments
def test_mafft_multi_param():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'mafft', '--clustalout --noscore')
    assert align_to_hash(tester) == '2b8bf89e7459fe9d0b1f29628df6307e'


# ###############################################  'li', '--list_ids' ################################################ #
def test_list_ids():
    tester = Alb.list_ids(alb_objects[0])
    tester = md5(str(tester).encode()).hexdigest()
    assert tester == "bdb865d00bef4eb3becaaae8f7bb94cd"

    tester = Alb.list_ids(alb_objects[8])
    tester = md5(str(tester).encode()).hexdigest()
    assert tester == "0f398f3ef5b73804d03b550d66b5462c"


# #################################### 'lc', '--lowercase' and 'uc', '--uppercase' ################################### #
lc_hashes = ["cb1169c2dd357771a97a02ae2160935d", "f59e28493949f78637691caeb617ab50", "52c23bd793c9761b7c0f897d3d757c12",
             "228e36a30e8433e4ee2cd78c3290fa6b", "17ff1b919cac899c5f918ce8d71904f6", "5af1cf061f003d3351417458c0d23811",
             "f3e98897f1bbb3d3d6c5507afa9f814e", "c0dce60745515b31a27de1f919083fe9", "e4d6766b7544557b9ddbdcbf0cde0c16",
             "12716bad78b2f7a40882df3ce183735b", "00661f7afb419c6bb8c9ac654af7c976"]

uc_hashes = ["52e74a09c305d031fc5263d1751e265d", "34eefeafabfc55811a5c9fe958b61490", "6e5542f41d17ff33afb530b4d07408a3",
             "b82538a4630810c004dc8a4c2d5165ce", "8b6737fe33058121fd99d2deee2f9a76", "b804b3c0077f342f8a5e8c36b8af627f",
             "747ca137dc659f302a07b0c39e989e54", "f35cbc6e929c51481e4ec31e95671638", "73e3da29aa78f4abb4bc6392b81cd279",
             "46e049a1be235d17f8379c293e1e393f", "6f3f234d796520c521cb85c66a3e239a"]

hashes = [(Alb._make_copies(alb_objects[indx]), uc_hash, lc_hashes[indx]) for indx, uc_hash in enumerate(uc_hashes)]


@pytest.mark.parametrize("alignbuddy,uc_hash,lc_hash", hashes)
def test_cases(alignbuddy, uc_hash, lc_hash):
    tester = Alb.uppercase(alignbuddy)
    assert align_to_hash(tester) == uc_hash
    tester = Alb.lowercase(tester)
    assert align_to_hash(tester) == lc_hash


# ###############################################  'ns', '--num_seqs' ################################################ #
def test_num_seqs():
    for i in [0, 2, 3, 4, 6, 7]:
        assert Alb.num_seqs(alb_objects[i]) == [13]
    for i in [1, 5]:
        assert Alb.num_seqs(alb_objects[i]) == [8]
    assert Alb.num_seqs(alb_objects[8]) == [20, 8]
    for i in [9, 10]:
        assert Alb.num_seqs(alb_objects[i]) == [20, 13]
    for i in [11, 12]:
        assert Alb.num_seqs(alb_objects[i]) == [8, 21]


# ###########################################  'oi', '--order_ids' ############################################ #
@pytest.mark.parametrize("alignbuddy", deepcopy(alb_objects))
def test_order_ids(alignbuddy):
    expected = Alb.list_ids(alignbuddy)
    for _id, align in enumerate(expected):
        expected[_id] = sorted(align)
    actual = Alb.list_ids(Alb.order_ids(alignbuddy))
    assert expected == actual

# ###########################################  'pr', '--pull_rows' ############################################ #
pr_hashes = ['2c0a60cd3f534d46662ed61272481898', '0f8e6552d9ac5a2bf7b3bd76fa54c9ca', '9ba71d3045d1929e34ae49d84816292e',
             'f41422a10e3b5a7e53c6ba5cc9f28875', '9a1692f1762e67ca0364de22b124dfee', '80fe4fa125e3d944f914fd0c8b923076',
             '52ee82ee31e0d0f9c84a577b01580f25', '068ad86e5e017bd88678b8f5b24512c1']
pr_hashes = [(Alb._make_copies(alb_objects[x]), value) for x, value in enumerate(pr_hashes)]


@pytest.mark.parametrize("alignbuddy,next_hash", pr_hashes)
def test_pull_rows(alignbuddy, next_hash):
    Alb.pull_rows(alignbuddy, 'Mle-Panxα[567]')
    assert align_to_hash(alignbuddy) == next_hash


# ###########################################  'ri', '--rename_ids' ############################################ #
def test_rename_ids_several():
    tester = Alb.AlignBuddy(resource("concat_alignment_file.phyr"))
    Alb.rename(tester, 'Mle', 'Xle')
    assert align_to_hash(tester) == '61d03559088c5bdd0fdebd7a8a2061fd'


def test_rename_ids_all_same():
    tester = Alb._make_copies(alb_objects[0])
    Alb.rename(tester, 'Mle', 'Xle')
    assert align_to_hash(tester) == '5a0c20a41fea9054f5476e6fad7c81f6'


# ###########################################  'r2d', '--back_transcribe' ############################################ #
@pytest.mark.parametrize("alignbuddy", deepcopy(alb_objects[0:4]))
def test_back_transcribe(alignbuddy):
    print(alignbuddy)
    original = align_to_hash(alignbuddy)
    final = align_to_hash(Alb.rna2dna(Alb.dna2rna(alignbuddy)))
    assert original == final


def test_back_transcribe_pep_exception():
    with pytest.raises(ValueError):
        Alb.rna2dna(deepcopy(alb_objects[4]))


# ###########################################  'stf', '--split_alignbuddy' ########################################### #
def test_split_alignment():
    tester = Alb.AlignBuddy(resource("concat_alignment_file.phyr"))
    output = Alb.split_alignbuddy(tester)
    for buddy in output:
        assert buddy.alignments[0] in tester.alignments

# ###########################################  'tr', '--translate' ############################################ #
hashes = ["fa915eafb9eb0bfa0ed8563f0fdf0ef9", "5064c1d6ae6192a829972b7ec0f129ed", "ce423d5b99d5917fbef6f3b47df40513",
          "2340addad40e714268d2523cdb17a78c", "6c66f5f63c5fb98f5855fb1c847486ad", "d9527fe1dfd2ea639df267bb8ee836f7"]
hashes = [(Alb._make_copies(alb_objects[nucl_indices[indx]]), next_hash) for indx, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_translate1(alignbuddy, next_hash, capsys):
    tester = Alb.translate_cds(alignbuddy)
    assert align_to_hash(tester) == next_hash

    if next_hash == "6c66f5f63c5fb98f5855fb1c847486ad":
        out, err = capsys.readouterr()
        assert err == "Warning: First codon 'GGT' is not a start codon in Ael_PanxβQ___CDS\n"


def test_translate2():
    # Protein input
    with pytest.raises(TypeError):
        Alb.translate_cds(alb_objects[4])

    # Non-standard length
    tester = Alb._make_copies(alb_objects[0])
    tester.alignments[0][0].seq = Seq(re.sub("-$", "t", str(tester.alignments[0][0].seq)),
                                      alphabet=tester.alignments[0][0].seq.alphabet)
    Alb.translate_cds(tester)
    assert align_to_hash(tester) == "fa915eafb9eb0bfa0ed8563f0fdf0ef9"

    # Final codon not stop and stop codon not at end of seq
    tester = Alb._make_copies(alb_objects[0])
    tester.alignments[0][0].seq = Seq(re.sub("---$", "ttt", str(tester.alignments[0][0].seq)),
                                      alphabet=tester.alignments[0][0].seq.alphabet)
    Alb.translate_cds(tester)
    assert align_to_hash(tester) == "2ca5b0bac98226ee6a53e17503f12197"

    # Non-standard codon
    tester = Alb.AlignBuddy(resource("ambiguous_dna_alignment.fa"))
    Alb.translate_cds(tester)
    assert align_to_hash(tester) == "ab8fb45a38a6e5d553a29f3613bbc1a1"


# ###########################################  'tm', '--trimal' ############################################ #
hashes = ['063e7f52f7c7f19da3a624730cd725a5', '0e2c65d9fc8d4b31cc352b3894837dc1', '4c08805a44e9a0cef08abc53c80f6b4c',
          '9478e5441b88680470f3a4a4db537467', "b1bd83debd405431f290e2e2306a206e", '049b7f9170505ea0799004d964ef33fb',
          '97ea3e55425894d3fe4b817deab003c3', '1f2a937d2c3b020121000822e62c4b96', '4690f58abd6d7f4f3c0d610ea25461c8',
          'b8e65b0a00f55286d57aa76ccfcf04ab', '49a17364aa4bd086c7c432de7baabd07', '82aaad11c9496d8f7e959dd5ce06df4d',
          '4e709de80531c358197f6e1f626c9c58']
hashes = [(Alb._make_copies(alb_objects[x]), value) for x, value in enumerate(hashes)]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_trimal(alignbuddy, next_hash):
    Alb.trimal(alignbuddy, .7)
    assert align_to_hash(alignbuddy) == next_hash


def test_trimal2():
    tester = Alb._make_copies(alb_objects[8])
    assert align_to_hash(Alb.trimal(tester, 'clean')) == "a94edb2636a7c9cf177993549371b3e6"
    assert align_to_hash(Alb.trimal(tester, 'gappyout')) == "2ac19a0eecb9901211c1c98a7a203cc2"
    assert align_to_hash(Alb.trimal(tester, 'all')) == "caebb7ace4940cea1b87667e5e113acb"
    with pytest.raises(ValueError):
        Alb.trimal(tester, "Foo")
