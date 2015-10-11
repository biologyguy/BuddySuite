#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

import pytest
from hashlib import md5
import os
import re
import sys
import argparse
import io
from copy import deepcopy

from Bio.Alphabet import IUPAC
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC
from unittest import mock

sys.path.insert(0, "./")
import buddy_resources as br
import SeqBuddy as Sb
import MyFuncs

VERSION = Sb.VERSION
WRITE_FILE = MyFuncs.TempFile()
BACKUP_PATH = os.environ["PATH"]


def fmt(prog):
    return br.CustomHelpFormatter(prog)

parser = argparse.ArgumentParser(prog="SeqBuddy.py", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                 description='''\
\033[1mSeqBuddy\033[m
  See your sequence files. Be your sequence files.

\033[1mUsage examples\033[m:
  SeqBuddy.py "/path/to/seq_file" -<cmd>
  SeqBuddy.py "/path/to/seq_file" -<cmd> | SeqBuddy.py -<cmd>
  SeqBuddy.py "ATGATGCTAGTC" -f "raw" -<cmd>
''')

br.flags(parser, ("sequence", "Supply file path(s) or raw sequence. If piping sequences "
                              "into SeqBuddy this argument can be left blank."),
         br.sb_flags, br.sb_modifiers, VERSION)

in_args = parser.parse_args()


def seqs_to_hash(_seqbuddy, mode='hash'):
    if _seqbuddy.out_format in ["gb", "genbank"]:
        for _rec in _seqbuddy.records:
            try:
                if re.search("(\. )+", _rec.annotations['organism']):
                    _rec.annotations['organism'] = "."
            except KeyError:
                pass

    if _seqbuddy.out_format == "phylipi":
        WRITE_FILE.write(Sb._phylipi(_seqbuddy, "relaxed"))
    elif _seqbuddy.out_format == "phylipis":
        WRITE_FILE.write(Sb._phylipi(_seqbuddy, "strict"))
    else:
        _seqbuddy.write(WRITE_FILE.path)

    seqs_string = "{0}\n".format(WRITE_FILE.read().strip())

    if mode != "hash":
        return seqs_string

    _hash = md5(seqs_string.encode()).hexdigest()
    return _hash

root_dir = os.getcwd()


def string2hash(_input):
    return md5(_input.encode()).hexdigest()


def resource(file_name):
    return "{0}/unit_test_resources/{1}".format(root_dir, file_name)


seq_files = ["Mnemiopsis_cds.fa", "Mnemiopsis_cds.gb", "Mnemiopsis_cds.nex",
             "Mnemiopsis_cds.phy", "Mnemiopsis_cds.phyr", "Mnemiopsis_cds.stklm",
             "Mnemiopsis_pep.fa", "Mnemiopsis_pep.gb", "Mnemiopsis_pep.nex",
             "Mnemiopsis_pep.phy", "Mnemiopsis_pep.phyr", "Mnemiopsis_pep.stklm",
             "ambiguous_dna.fa", "ambiguous_rna.fa"]


@pytest.mark.parametrize("seq_file", seq_files)
def test_instantiate_seqbuddy_from_file(seq_file):
    assert type(Sb.SeqBuddy(resource(seq_file))) == Sb.SeqBuddy


@pytest.mark.parametrize("seq_file", seq_files)
def test_instantiate_seqbuddy_from_handle(seq_file):
    with open(resource(seq_file), 'r') as ifile:
        assert type(Sb.SeqBuddy(ifile)) == Sb.SeqBuddy


@pytest.mark.parametrize("seq_file", seq_files)
def test_instantiate_seqbuddy_from_raw(seq_file):
    with open(resource(seq_file), 'r') as ifile:
        assert type(Sb.SeqBuddy(ifile.read(), in_format="raw")) == Sb.SeqBuddy

    with open(resource(seq_file), 'r') as ifile:
        assert type(Sb.SeqBuddy(ifile, in_format="raw")) == Sb.SeqBuddy


@pytest.mark.parametrize("seq_file", seq_files)
def test_instantiate_seqbuddy_from_seqbuddy(seq_file):
    input_buddy = Sb.SeqBuddy(resource(seq_file))
    tester = Sb.SeqBuddy(input_buddy)
    assert seqs_to_hash(input_buddy) == seqs_to_hash(tester)


def test_alpha_arg_dna():
    tester = Sb.SeqBuddy(resource(seq_files[0]), alpha='dna')
    assert tester.alpha is IUPAC.ambiguous_dna


def test_alpha_arg_rna():
    tester = Sb.SeqBuddy(resource(seq_files[0]), alpha='rna')
    assert tester.alpha is IUPAC.ambiguous_rna


def test_alpha_arg_prot():
    tester = Sb.SeqBuddy(resource(seq_files[6]), alpha='prot')
    assert tester.alpha is IUPAC.protein


def test_alpha_arg_guess():
    tester = Sb.SeqBuddy(resource(seq_files[0]), alpha='dgasldfkjhgaljhetlfdjkfg')
    assert tester.alpha is IUPAC.ambiguous_dna


def test_seqlist_error():
    with pytest.raises(TypeError):
        Sb.SeqBuddy([str, dict])


# ##################### SeqBuddy methods ###################### ##
def test_to_dict():
    tester = str(Sb.SeqBuddy(resource("Mnemiopsis_cds.fa")).to_dict())
    tester = ''.join(sorted(tester))
    assert string2hash(tester) == '06f50839f94e8f917311b682837461fd'

    with pytest.raises(RuntimeError):
        tester = Sb.SeqBuddy(">duplicate_id\nATGCTCGTA\n>duplicate_id\nATGCTCGTCGATGCT\n")
        tester.to_dict()


def test_to_string(capsys):
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    assert string2hash(str(tester)) == "b831e901d8b6b1ba52bad797bad92d14"
    tester.print()
    out, err = capsys.readouterr()
    assert string2hash(out) == "b831e901d8b6b1ba52bad797bad92d14"

    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.gb"))
    assert string2hash(str(tester)) == "2e02a8e079267bd9add3c39f759b252c"

    tester.out_format = "phylipi"
    assert string2hash(str(tester)) == "da16cb90b53e9f29ddfe3bf0e3faee22"

    tester.out_format = "phylipis"
    assert string2hash(str(tester)) == "0b00c1521259b0e1389e58cb5d414048"

    tester.records = []
    assert str(tester) == "Error: No sequences in object.\n"


# Now that we know that all the files are being turned into SeqBuddy objects okay, make them all objects so it doesn't
# need to be done over and over for each subsequent test.
sb_objects = [Sb.SeqBuddy(resource(x)) for x in seq_files]


# ################################################# HELPER FUNCTIONS ################################################# #
# ######################  'GuessError' ###################### #
def test_guesserror_raw_seq():
    with pytest.raises(Sb.GuessError):
        Sb.SeqBuddy("JSKHGLHGLSDKFLSDYUIGJVSBDVHJSDKGIUSUEWUIOIFUBCVVVBVNNJS{QF(*&#@$(*@#@*(*(%")
    try:
        Sb.SeqBuddy("JSKHGLHGLSDKFLSDYUIGJVSBDVHJSDKGIUSUEWUIOIFUBCVVVBVNNJS{QF(*&#@$(*@#@*(*(%")
    except Sb.GuessError as e:
        assert str(e) == "File not found, or could not determine format from raw input"


def test_guesserror_infile():
    with pytest.raises(Sb.GuessError):
        Sb.SeqBuddy(resource("gibberish.fa"))


def test_guesserror_in_handle():
    with pytest.raises(Sb.GuessError):
        with open(resource("gibberish.fa"), "r") as ifile:
            Sb.SeqBuddy(ifile)


def test_no__input():
    with pytest.raises(TypeError):
        Sb.SeqBuddy()


# ######################  '_make_copy' ###################### #
def test_make_copy():
    assert seqs_to_hash(Sb._make_copy(sb_objects[0])) == seqs_to_hash(sb_objects[0])


# ######################  '_check_for_blast_bin' ###################### #
@pytest.mark.internet
@pytest.mark.slow
def test_check_blast_bin(capsys):
    for _bin in ["blastn", "blastp", "blastdbcmd", "makeblastdb"]:
        assert Sb._check_for_blast_bin(_bin)

    with mock.patch.dict(os.environ, {"PATH": ""}):
        with mock.patch('MyFuncs.ask', return_value=False):
            assert not Sb._check_for_blast_bin("blastp")

        with mock.patch('MyFuncs.ask', return_value=True):
            with mock.patch("SeqBuddy._download_blast_binaries", return_value=False):
                assert not Sb._check_for_blast_bin("foo")
                out, err = capsys.readouterr()
                assert "Failed to download foo" in err

            assert not Sb._check_for_blast_bin("blastp")
            out, err = capsys.readouterr()
            assert "blastp downloaded" in err
            assert os.path.isfile("./blastp")
            os.remove("./blastp")


# ######################  '_download_blast_binaries' ###################### #
@pytest.mark.internet
@pytest.mark.slow
def test_dl_blast_bins():
    sys_memory = int(sys.maxsize)
    for platform in ["darwin", "linux32", "linux64", "win"]:
        sys.maxsize = 2000000000 if platform == "linux32" else sys_memory
        platform = "linux" if "linux" in platform else platform
        Sb._download_blast_binaries(ignore_pre_install=True, system=platform)
        if platform != "win":
            assert os.path.isfile('./blastdbcmd')
            assert os.path.isfile('./blastn')
            assert os.path.isfile('./blastp')
            os.remove("./blastdbcmd")
            os.remove("./blastn")
            os.remove("./blastp")
        else:
            assert os.path.isfile('./blastdbcmd.exe')
            assert os.path.isfile('./blastn.exe')
            assert os.path.isfile('./blastp.exe')
            os.remove("./blastdbcmd.exe")
            os.remove("./blastn.exe")
            os.remove("./blastp.exe")


# ######################  '_feature_rc' ###################### #
def test_feature_rc():
    tester = Sb._make_copy(sb_objects[1])
    seq1 = tester.records[0]
    hashed = ["dc71a33b64a766da8653c19f22fc4caa", "9ab8296fb3443198674d90abe3311ba6",
              "10018d1b15c7f76a6333ac3bf96d2d07", "273463b9eace12d2eeadbf272692d73e",
              "c452d66d13120cd6eb5f041b7c37dd27", "1811b0695dba1fc3fe431a6ee00ef359"]
    for feature in zip(seq1.features, hashed):
        assert string2hash(str(Sb._feature_rc(feature[0], 1203))) == feature[1]

    with pytest.raises(TypeError):
        feature = seq1.features[0]
        feature.location = {}
        Sb._feature_rc(feature, 1203)


# ######################  '_format_to_extension' ###################### #
def test_format_to_extension():
    ext_dict = {'fasta': 'fa', 'fa': 'fa', 'genbank': 'gb', 'gb': 'gb', 'nexus': 'nex', 'nex': 'nex', 'phylip': 'phy',
                'phy': 'phy', 'phylip-relaxed': 'phyr', 'phyr': 'phyr', 'stockholm': 'stklm', 'stklm': 'stklm'}
    for i, j in ext_dict.items():
        assert j == Sb._format_to_extension(i)


# ######################  '_guess_alphabet' ###################### #
def test_guess_alphabet():
    tester = Sb._make_copy(sb_objects[0])
    assert Sb._guess_alphabet(tester) == IUPAC.ambiguous_dna

    tester = Sb._make_copy(sb_objects[6])
    assert Sb._guess_alphabet(tester) == IUPAC.protein

    tester = Sb.SeqBuddy(resource("Mnemiopsis_rna.fa"))
    assert Sb._guess_alphabet(tester) == IUPAC.ambiguous_rna

    tester = Sb.SeqBuddy(">Seq1", in_format="fasta")
    assert not Sb._guess_alphabet(tester)


# ######################  '_guess_format' ###################### #
def test_guess_format():
    assert Sb._guess_format(["foo", "bar"]) == "gb"
    assert Sb._guess_format(sb_objects[0]) == "fasta"
    assert Sb._guess_format(resource("Mnemiopsis_cds.fa")) == "fasta"
    assert Sb._guess_format(resource("blank.fa")) == "empty file"
    with pytest.raises(Sb.GuessError):
        Sb._guess_format("foo")

    temp_file = MyFuncs.TempFile()
    temp_file.write('''\
<?xml version="1.0" encoding="ISO-8859-1"?>
<nex:nexml
    version="0.9"
    xsi:schemaLocation="http://www.nexml.org/2009 ../xsd/nexml.xsd"
    xmlns="http://www.nexml.org/2009"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xmlns:xml="http://www.w3.org/XML/1998/namespace"
    xmlns:nex="http://www.nexml.org/2009"
    xmlns:xsd="http://www.w3.org/2001/XMLSchema#"
''')

    assert not Sb._guess_format(temp_file.path)


# ######################  'phylipi' ###################### #
def test_phylipi():
    tester = Sb._phylipi(Sb.SeqBuddy(resource("Mnemiopsis_cds.nex")), _format="relaxed")
    tester = "{0}\n".format(tester.rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "c5fb6a5ce437afa1a4004e4f8780ad68"

    tester = Sb._phylipi(Sb.SeqBuddy(resource("Alignments_pep.phy")), _format="relaxed")
    tester = "{0}\n".format(tester.rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "5f70d8a339b922f27c308d48280d715f"

    tester = Sb._phylipi(Sb.SeqBuddy(resource("Mnemiopsis_cds.nex")), _format="strict")
    tester = "{0}\n".format(tester.rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "270f1bac51b2e29c0e163d261795c5fe"

    tester = Sb._phylipi(Sb.SeqBuddy(resource("Alignments_pep.phy")), _format="strict")
    tester = "{0}\n".format(tester.rstrip())
    tester = md5(tester.encode()).hexdigest()
    assert tester == "5f70d8a339b922f27c308d48280d715f"


# ######################  '_shift_features' ###################### #
def test_shift_features():
    tester = Sb._make_copy(sb_objects[1])
    features = tester.records[0].features
    assert string2hash(str(Sb._shift_features(features, 3, 1203))) == "cc002df59db8a04f47cb5c764f8a1e1f"

    tester = Sb._make_copy(sb_objects[1])
    features = tester.records[0].features
    assert string2hash(str(Sb._shift_features(features, -50, 1203))) == "86dcc19cd51ba3acdaf48e8c15bf7f1a"

    tester = Sb._make_copy(sb_objects[1])
    features = tester.records[0].features
    features[0].location.parts = []
    assert string2hash(str(Sb._shift_features(features, 3, 1203))) == "467ee934600573e659053b8117447986"

    tester = Sb._make_copy(sb_objects[1])
    features = tester.records[0].features
    features[0].location.parts = [FeatureLocation(3, 50, strand=+1)]
    assert string2hash(str(Sb._shift_features(features, 3, 1203))) == "94bf6774e97fcecaef858cbdf811def4"

    with pytest.raises(TypeError):
        tester = Sb._make_copy(sb_objects[1])
        features = tester.records[0].features
        features[0].location = [dict]
        Sb._shift_features(features, 3, 1203)


# ######################  '_stdout and _stderr' ###################### #
def test_stdout(capsys):
    Sb._stdout("Hello std_out", quiet=False)
    out, err = capsys.readouterr()
    assert out == "Hello std_out"

    Sb._stdout("Hello std_out", quiet=True)
    out, err = capsys.readouterr()
    assert out == ""


def test_stderr(capsys):
    Sb._stderr("Hello std_err", quiet=False)
    out, err = capsys.readouterr()
    assert err == "Hello std_err"

    Sb._stderr("Hello std_err", quiet=True)
    out, err = capsys.readouterr()
    assert err == ""


# ################################################ MAIN API FUNCTIONS ################################################ #
# ##################### 'ano', '--annotate' ###################### ##
def test_annotate_pattern():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), pattern='α4')
    assert seqs_to_hash(tester) == '527b2730d23d2ae4d708dbc55f560604'


def test_annotate_no_pattern():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100))
    assert seqs_to_hash(tester) == 'e8afc1fd9b0d5833972acfd2b719e42a'


def test_annotate_compoundlocation():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', [(1, 100), (200, 250)])
    assert seqs_to_hash(tester) == '4a294c1ea5783fe04e909d767e7dfce2'


def test_annotate_list_str():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', ['1-100', '(200-250)'])
    assert seqs_to_hash(tester) == '4a294c1ea5783fe04e909d767e7dfce2'


def test_annotate_str():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', '1-100, (200-250)')
    assert seqs_to_hash(tester) == '4a294c1ea5783fe04e909d767e7dfce2'


def test_annotate_fl_obj():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', FeatureLocation(start=0, end=100))
    assert seqs_to_hash(tester) == 'e8afc1fd9b0d5833972acfd2b719e42a'


def test_annotate_cl_obj():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', CompoundLocation([FeatureLocation(start=0, end=100),
                                                                   FeatureLocation(start=199, end=250)],
                                                                  operator='order'))
    assert seqs_to_hash(tester) == '4a294c1ea5783fe04e909d767e7dfce2'


def test_annotate_typerror():
    with pytest.raises(TypeError):
        tester = Sb._make_copy(sb_objects[1])
        Sb.annotate(tester, 'misc_feature', 5)


def test_annotate_pos_strand():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand='+')
    assert seqs_to_hash(tester) == 'e8afc1fd9b0d5833972acfd2b719e42a'


def test_annotate_neg_strand():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand='-')
    assert seqs_to_hash(tester) == '9de19b80180cdc65ca8a5b8e9c0e61e9'


def test_annotate_no_strand():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand=0)
    assert seqs_to_hash(tester) == 'e8afc1fd9b0d5833972acfd2b719e42a'


def test_annotate_qualifier_dict():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers={'foo': 'bar', 'hello': 'world'})
    assert seqs_to_hash(tester) == 'a948567058575491dbb3ed6ebee1f9b6'


def test_annotate_qualifier_list():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers=['foo=bar', 'hello=Blahh, =world'])
    assert seqs_to_hash(tester) == '5f8d9b4874754e50877126994141cc92'


def test_annotate_qualifier_error():
    tester = Sb._make_copy(sb_objects[1])
    with pytest.raises(TypeError):
        Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers=tuple)


def test_annotate_out_of_range():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', [(-10, 100), (200, 10000)])
    assert seqs_to_hash(tester) == '57f57b5e46c547542ea30e932ec4a42b'

    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', [(1, 10000)])
    assert seqs_to_hash(tester) == 'faa0e801376b3daba753890c18f0f029'

    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', [(1, 10000), (20000, 30000)])
    assert seqs_to_hash(tester) == 'faa0e801376b3daba753890c18f0f029'

    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.annotate(tester, 'misc_feature', FeatureLocation(start=-10, end=100))
    assert seqs_to_hash(tester) == 'e8afc1fd9b0d5833972acfd2b719e42a'


def test_annotate_protein():
    tester = Sb._make_copy(sb_objects[7])
    tester = Sb.annotate(tester, 'misc_feature', (1, 100))
    assert seqs_to_hash(tester) == '34f9601d72d5c37a8fd301ff90088cd8'


def test_annotate_unrec_strand(capsys):
    tester = Sb._make_copy(sb_objects[1])
    Sb.annotate(tester, 'misc_feature', (1, 100), strand='foo')
    out, err = capsys.readouterr()
    assert err == "Warning: strand input not recognized. Value set to None."


# ######################  'asl', '--ave_seq_length' ###################### #
@pytest.mark.parametrize("seqbuddy", [Sb._make_copy(x) for x in sb_objects[0:6] if len(x.records) == 13])
def test_ave_seq_length_dna(seqbuddy):
    assert round(Sb.ave_seq_length(seqbuddy, clean=True), 2) == 1285.15


@pytest.mark.parametrize("seqbuddy", [Sb._make_copy(x) for x in sb_objects[6:12] if len(x.records) == 13])
def test_ave_seq_length_pep(seqbuddy):
    assert round(Sb.ave_seq_length(seqbuddy, clean=True), 2) == 427.38


# ######################  'btr', '--back_translate' ###################### #
# Only fasta and genbank
hashes = ["1b14489a78bfe8255c777138877b9648", "b6bcb4e5104cb202db0ec4c9fc2eaed2",
          "859ecfb88095f51bfaee6a1d1abeb50f", "ba5c286b79a3514fba0b960ff81af25b",
          "952a91a4506afb57f27136aa1f2a8af9", "40c4a3e08c811b6bf3be8bedcb5d65a0",
          "3a3ee57f8dcde25c99a655494b218928", "bc9f1ec6ec92c30b5053cd9bb6bb6f53"]
organisms = ['human', 'human', 'yeast', 'yeast', 'ecoli', 'ecoli', 'mouse', 'mouse']
hashes = [(Sb._make_copy(sb_objects[sb_obj_indx]), organisms[indx], hashes[indx]) for indx, sb_obj_indx in
          enumerate([6, 7, 6, 7, 6, 7, 6, 7])]


@pytest.mark.parametrize("seqbuddy,_organism,next_hash", hashes)
def test_back_translate(seqbuddy, _organism, next_hash):
    seqbuddy.alpha = IUPAC.protein
    tester = Sb.back_translate(seqbuddy, 'OPTIMIZED', _organism)
    assert seqs_to_hash(tester) == next_hash


def test_back_translate_nucleotide_exception():
    with pytest.raises(TypeError):
        Sb.back_translate(sb_objects[1])


def test_back_translate_bad_mode():
    with pytest.raises(AttributeError):
        Sb.back_translate(Sb._make_copy(sb_objects[6]), 'fgsdjkghjdalgsdf', 'human')


def test_back_translate_bad_organism():
    seqbuddy = Sb._make_copy(sb_objects[6])
    with pytest.raises(AttributeError):
        Sb.back_translate(seqbuddy, 'OPTIMIZED', 'fgsdjkghjdalgsdf')


# ######################  'bl2s', '--bl2seq' ###################### #
def test_bl2seq_cds():
    seqbuddy = Sb._make_copy(sb_objects[0])
    result = Sb.bl2seq(seqbuddy)[1]
    assert md5(result.encode()).hexdigest() == '339377aee781fb9d01456f04553e3923'


def test_bl2seq_pep():
    seqbuddy = Sb._make_copy(sb_objects[6])
    result = Sb.bl2seq(seqbuddy)[1]
    assert md5(result.encode()).hexdigest() == '4c722c4db8bd5c066dc76ebb94583a37'


def test_bl2_no_binary():
    os.environ["PATH"] = ""
    with mock.patch('builtins.input', return_value="n"):
        with pytest.raises(RuntimeError):
            seqbuddy = Sb._make_copy(sb_objects[0])
            Sb.bl2seq(seqbuddy)

        with pytest.raises(RuntimeError):
            seqbuddy = Sb._make_copy(sb_objects[6])
            Sb.bl2seq(seqbuddy)
    os.environ["PATH"] = BACKUP_PATH


# ######################  'bl', '--blast' ###################### #
def test_blastn():
    tester = Sb.pull_recs(Sb._make_copy(sb_objects[0]), '8')
    tester = Sb.blast(tester, blast_db=resource("blast/Mnemiopsis_cds.n"))
    assert seqs_to_hash(tester) == "95c417b6c2846d1b7a1a07f50c62ff8a"

    with pytest.raises(RuntimeError) as e:
        tester = Sb._make_copy(sb_objects[0])
        Sb.blast(tester, blast_db=resource("Mnemiopsis_cds.nhr"))
    assert "The .nhr file of your blast database was not found" in str(e.value)

    with mock.patch("SeqBuddy._check_for_blast_bin", return_value=False):
        with pytest.raises(SystemError) as e:
            Sb.blast(tester, blast_db=resource("blast/Mnemiopsis_cds.n"))
        assert 'blastn not found in system path' in str(e.value)

    tester = Sb.SeqBuddy(">Seq1\nATGCGCGCTACGCTAGCTAGCTAGCTCGCATGCAT")
    tester = Sb.blast(tester, blast_db=resource("blast/Mnemiopsis_cds.n"))
    assert len(tester.records) == 0


def test_blastp():
    seqbuddy = Sb.pull_recs(Sb.SeqBuddy(resource(seq_files[6])), '8')
    tester = Sb.blast(seqbuddy, blast_db=resource("blast/Mnemiopsis_pep.p"))
    assert seqs_to_hash(tester) == "4237c79672c1cf1d4a9bdb160a53a4b9"

    with pytest.raises(RuntimeError) as e:
        tester = Sb._make_copy(sb_objects[6])
        Sb.blast(tester, blast_db=resource("Mnemiopsis_cds.phr"))
    assert "The .phr file of your blast database was not found" in str(e.value)

    with mock.patch("SeqBuddy._check_for_blast_bin", return_value=False):
        with pytest.raises(SystemError) as e:
            Sb.blast(tester, blast_db=resource("blast/Mnemiopsis_pep.n"))
        assert 'blastp not found in system path' in str(e.value)


# ######################  'cs', '--clean_seq'  ###################### #
def test_clean_seq():
    # Protein
    tester = Sb._make_copy(sb_objects[6])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert seqs_to_hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert seqs_to_hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"

    # DNA
    tester = Sb._make_copy(sb_objects[12])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert seqs_to_hash(tester) == "71b28ad2730a9849f2ba0f70e9e51a9f"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert seqs_to_hash(tester) == "5fd0b78e37c81e0fa727db34a37cc743"

    # RNA
    tester = Sb._make_copy(sb_objects[13])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert seqs_to_hash(tester) == "cdb1b963536d57efc7b7f87d2bf4ad22"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert seqs_to_hash(tester) == "ef61174330f2d9cf5da4f087d12ca201"

    # Alignment formats should raise an error because seq lengths change
    with pytest.raises(ValueError):
        tester = Sb.clean_seq(Sb._make_copy(sb_objects[2]))
        tester.write("/dev/null")
        tester = Sb.clean_seq(Sb._make_copy(sb_objects[3]))
        tester.write("/dev/null")
        tester = Sb.clean_seq(Sb._make_copy(sb_objects[4]))
        tester.write("/dev/null")
        tester = Sb.clean_seq(Sb._make_copy(sb_objects[5]))
        tester.write("/dev/null")


# ######################  'cmp', '--complement' ###################### #
hashes = ["e4a358ca57aca0bbd220dc6c04c88795", "3366fcc6ead8f1bba4a3650e21db4ec3",
          "365bf5d08657fc553315aa9a7f764286", "10ce87a53aeb5bd4f911380ebf8e7a85",
          "8e5995813da43c7c00e98d15ea466d1a", "5891348e8659290c2355fabd0f3ba4f4"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_complement(seqbuddy, next_hash):
    tester = Sb.complement(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_complement_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.complement(sb_objects[6])


# ######################  'cts', '--concat_seqs' ###################### #
# ToDo: Test the _clean parameter
hashes = ["2e46edb78e60a832a473397ebec3d187", "7421c27be7b41aeedea73ff41869ac47",
          "494988ffae2ef3072c1619eca8a0ff3b", "710cad348c5560446daf2c916ff3b3e4",
          "494988ffae2ef3072c1619eca8a0ff3b", "494988ffae2ef3072c1619eca8a0ff3b",
          "46741638cdf7abdf53c55f79738ee620", "8d0bb4e5004fb6a1a0261c30415746b5",
          "2651271d7668081cde8012db4f9a6574", "36526b8e0360e259d8957fa2261cf45a",
          "2651271d7668081cde8012db4f9a6574", "2651271d7668081cde8012db4f9a6574"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_concat_seqs(seqbuddy, next_hash):
    tester = Sb.concat_seqs(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_concat_seqs_clean():
    tester = Sb.concat_seqs(Sb._make_copy(sb_objects[2]), clean=True)
    assert seqs_to_hash(tester) == "2e46edb78e60a832a473397ebec3d187"


# ##################### 'cc', 'count_codons' ###################### ##
def test_count_codons_dna():
    tester = Sb._make_copy(sb_objects[0])
    tester.records[0].seq = tester.records[0].seq[:-4]
    counter = Sb.count_codons(tester)
    assert seqs_to_hash(counter[0]) == '4dbd9a5c68d85bb200c75b309fdaeeca'
    assert string2hash(str(counter[1])) == '8d47313f3db02aee48575ff8ff4741b4'


def test_count_codons_rna():
    tester = Sb.count_codons(Sb.dna2rna(Sb._make_copy(sb_objects[0])))[1]
    assert md5(str(tester).encode()).hexdigest() == 'b91daa8905533b5885d2067d9d6ffe36'


def test_count_codons_dna_badchar():
    tester = Sb.count_codons(Sb.insert_sequence(Sb._make_copy(sb_objects[0]), 'PPP', 'end'))[1]
    assert md5(str(tester).encode()).hexdigest() == '9aba116675fe0e9eaaf43e5c6e0ba99d'


def test_count_codons_pep_exception():
    tester = Sb._make_copy(sb_objects[6])
    with pytest.raises(TypeError):
        Sb.count_codons(tester)


# ######################  'cr', '--count_residues' ###################### #
def test_count_residues():
    # Unambiguous DNA
    tester = Sb.SeqBuddy(">seq1\nACGCGAAGCGAACGCGCAGACGACGCGACGACGACGACGCA", in_format="fasta")
    tester = Sb.count_residues(tester).to_dict()
    assert tester["seq1"].buddy_data["res_count"]['A'] == [13, 0.3170731707317073]

    tester = Sb._make_copy(sb_objects[0])
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert res_count['G'] == [265, 0.21703521703521703]
    assert "% Ambiguous" not in res_count and "U" not in res_count

    # Unambiguous RNA
    Sb.dna2rna(tester)
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert res_count['U'] == [356, 0.2915642915642916]
    assert "% Ambiguous" not in res_count and "T" not in res_count

    # Ambiguous DNA
    tester = Sb._make_copy(sb_objects[12])
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert "U" not in res_count
    assert res_count['Y'] == [1, 0.000819000819000819]
    assert res_count['% Ambiguous'] == 0.98

    # Ambiguous RNA
    tester = Sb._make_copy(sb_objects[13])
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert "T" not in res_count
    assert res_count['U'] == [353, 0.2891072891072891]
    assert res_count["Y"] == [1, 0.000819000819000819]
    assert res_count['% Ambiguous'] == 0.98

    # Protein
    tester = Sb._make_copy(sb_objects[6])
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


# ######################  'df', '--delete_features' ###################### #
def test_delete_features():
    tester = Sb._make_copy(sb_objects[1])
    tester = Sb.delete_features(tester, 'donor')
    assert seqs_to_hash(tester) == 'f84df6a77063c7def13babfaa0555bbf'


# ######################  'dl', '--delete_large' ###################### #
def test_delete_large():
    tester = Sb._make_copy(sb_objects[0])
    tester = Sb.delete_large(tester, 1285)
    assert seqs_to_hash(tester) == '25859dc69d46651a1e04a70c07741b35'


# ######################  'dm', '--delete_metadata' ###################### #
hashes = ["aa92396a9bb736ae6a669bdeaee36038", "544ab887248a398d6dd1aab513bae5b1", "cb1169c2dd357771a97a02ae2160935d",
          "d1524a20ef968d53a41957d696bfe7ad", "99d522e8f52e753b4202b1c162197459", "a50943ccd028b6f5fa658178fa8cf54d",
          "bac5dc724b1fee092efccd2845ff2513", "858e8475f7bc6e6a24681083a8635ef9", "17ff1b919cac899c5f918ce8d71904f6",
          "c934f744c4dac95a7544f9a814c3c22a", "6a3ee818e2711995c95372afe073490b", "e224c16f6c27267b5f104c827e78df33"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_delete_metadata(seqbuddy, next_hash):
    tester = Sb.delete_metadata(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


# ######################  'dr', '--delete_records' ###################### #
dr_hashes = ["54bdb42423b1d331acea18218101e5fc", "e2c03f1fa21fd27b2ff55f7f721a1a99", "6bc8a9409b1ef38e4f6f12121368883e",
             "a63320b6a679f97368329396e3b72bdd", "fe927a713e92b7de0c3a63996a4ce7c8", "4c97144c5337f8a40c4fd494e622bf0d"]
dr_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(dr_hashes)]


@pytest.mark.parametrize("seqbuddy, next_hash", dr_hashes)
def test_delete_records(seqbuddy, next_hash):
    tester = Sb.delete_records(seqbuddy, 'α2')
    assert seqs_to_hash(tester) == next_hash


def test_delete_records2():
    tester = Sb.delete_records(Sb._make_copy(sb_objects[0]), ['α1', 'α2'])
    assert seqs_to_hash(tester) == "eca4f181dae3d7998464ff71e277128f"

    with pytest.raises(ValueError) as e:
        Sb.delete_records(Sb._make_copy(sb_objects[0]), dict)
    assert "'patterns' must be a list or a string." in str(e.value)


# #####################  'drp', '--delete_repeats' ###################### ##
def test_delete_repeats_ids():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_dup_id.fa"))
    tester = Sb.delete_repeats(tester)
    tester = Sb.find_repeats(tester)
    assert len(tester[1]) == 0
    assert len(tester[2]) == 0

def test_delete_repeats_seqs():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_dup_seq.fa"))
    tester = Sb.delete_repeats(tester)
    tester = Sb.find_repeats(tester)
    assert len(tester[1]) == 0
    assert len(tester[2]) == 0


# ######################  'ds', '--delete_small' ###################### #
def test_delete_small():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Sb.delete_small(tester, 1285)
    assert seqs_to_hash(tester) == '196adf08d4993c51050289e5167dacdf'


# ######################  'd2r', '--transcribe' and 'r2d', '--back_transcribe' ###################### #
d2r_hashes = ["d2db9b02485e80323c487c1dd6f1425b", "9ef3a2311a80f05f21b289ff7f401fff",
              "f3bd73151645359af5db50d2bdb6a33d", "1371b536e41e3bca304794512122cf17",
              "866aeaca326891b9ebe5dc9d762cba2c", "45b511f34653e3b984e412182edee3ca"]
r2d_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "2e02a8e079267bd9add3c39f759b252c",
              "cb1169c2dd357771a97a02ae2160935d", "d1524a20ef968d53a41957d696bfe7ad",
              "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b"]

hashes = [(Sb._make_copy(sb_objects[indx]), d2r_hash, r2d_hashes[indx]) for indx, d2r_hash in enumerate(d2r_hashes)]


@pytest.mark.parametrize("seqbuddy,d2r_hash,r2d_hash", hashes)
def test_transcribe(seqbuddy, d2r_hash, r2d_hash):
    tester = Sb.dna2rna(seqbuddy)
    assert seqs_to_hash(tester) == d2r_hash
    tester = Sb.rna2dna(tester)
    assert seqs_to_hash(tester) == r2d_hash


def test_transcribe_pep_exception():  # Asserts that a ValueError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.dna2rna(sb_objects[6])


def test_back_transcribe_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.rna2dna(sb_objects[6])


# ######################  'er', '--extract_range' ###################### #
hashes = ["201235ed91ad0ed9a7021136487fed94", "3e791c6a6683516aff9572c24f38f0b3", "4063ab66ced2fafb080ceba88965d2bb",
          "0c857970ebef51b4bbd9c7b3229d7088", "e0e256cebd6ead99ed3a2a20b7417ba1", "d724df01ae688bfac4c6dfdc90027440",
          "904a188282f19599a78a9d7af4169de6", "b8413624b9e684a14fc9f398a62e3965", "6a27222d8f60ee8496cbe0c41648a116",
          "9ecc1d83eff77c61284869b088c833a1", "9c85530cd3e3aa628b0e8297c0c9f977", "38d571c9681b4fa420e3d8b54c507f9c"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_extract_range(seqbuddy, next_hash):
    tester = Sb.extract_range(seqbuddy, 50, 300)
    assert seqs_to_hash(tester) == next_hash


def test_extract_range_end_less_than_start():
    with pytest.raises(ValueError):
        Sb.extract_range(sb_objects[0], 500, 50)


# #####################  'fcpg', '--find_CpG' ###################### ##
def test_find_cpg():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.gb"))
    tester = Sb.find_cpg(tester)
    assert seqs_to_hash(tester) == "ce6aff066c03651401db627951862154"


# #####################  'fp', '--find_pattern' ###################### ##
def test_find_pattern():
    tester = Sb.find_pattern(Sb._make_copy(sb_objects[1]), "ATGGT")
    assert seqs_to_hash(tester) == "a2e110d3ba14ea5d1b236a8abef7341a"
    tester =  Sb.find_pattern(Sb._make_copy(sb_objects[1]), "ATg{2}T")
    assert seqs_to_hash(tester) == "c8a3e2e8b97f47c5cc41f04a44243e34"
    tester = Sb.find_pattern(Sb._make_copy(sb_objects[1]), "ATg{2}T", "tga.{1,6}tg")
    assert seqs_to_hash(tester) == "0e11b2c0e9451fbfcbe39e3b5be2cf60"


# #####################  'frp', '--find_repeats' ###################### ##
def test_find_repeats_ids():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_dup_id.fa"))
    assert 'Mle-Panxα12' in Sb.find_repeats(tester)[1]


def test_find_repeats_seqs():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_dup_seq.fa"))
    result = Sb.find_repeats(tester)[2]
    for key in result:
        assert 'Mle-Panxα1' in result[key]
        assert 'Mle-Dupα' in result[key]


def test_find_repeats_none():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_pep.fa"))
    tester = Sb.find_repeats(tester)
    assert len(tester[1]) == 0
    assert len(tester[2]) == 0


# ######################  'frs', '--find_restriction_sites' ###################### #
def test_restriction_sites():
    # No arguments passed in = commercial REs and any number of cut sites
    tester = Sb.find_restriction_sites(Sb._make_copy(sb_objects[1]))
    assert seqs_to_hash(tester[0]) == 'fce4c8bee040d5ea6fa4bf9985f7310f'
    assert md5(str(tester[1]).encode()).hexdigest() == "741ca6ca9204a067dce7398f15c6e350"

    # Specify a few REs and limit the number of cuts
    tester = Sb.find_restriction_sites(Sb._make_copy(sb_objects[1]), enzyme_group=["EcoRI", "KspI", "TasI", "Bme1390I"],
                                       min_cuts=2, max_cuts=4)
    assert seqs_to_hash(tester[0]) == 'c42b3bf0367557383000b897432fed2d'
    assert md5(str(tester[1]).encode()).hexdigest() == "0d2e5fdba6fed434495481397a91e56a"

    with pytest.raises(TypeError):
        Sb.find_restriction_sites(sb_objects[7])


# ######################  'hsi', '--hash_sequence_ids' ###################### #
def test_hash_seq_ids():
    tester = Sb.SeqBuddy(Sb._make_copy(sb_objects[0]))
    Sb.hash_sequence_ids(tester)
    assert len(tester.records[0].id) == 10

    tester = Sb.hash_sequence_ids(tester, 25)
    assert len(tester.records[0].id) == 25
    assert len(tester.hash_map) == 13


def test_hash_seq_ids_errors():
    tester = Sb.SeqBuddy(Sb._make_copy(sb_objects[0]))
    with pytest.raises(TypeError) as e:
        Sb.hash_sequence_ids(tester, "foo")
    assert str(e.value) == "Hash length argument must be an integer, not <class 'str'>"

    with pytest.raises(ValueError) as e:
        Sb.hash_sequence_ids(tester, 0)
    assert str(e.value) == "Hash length must be greater than 0"

    tester.records *= 10
    with pytest.raises(ValueError) as e:
        Sb.hash_sequence_ids(tester, 1)
    assert "Insufficient number of hashes available to cover all sequences." in str(e.value)


# ##################### 'is', 'insert_seqs' ###################### ##
def test_insert_seqs_start():
    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', 'start')) == 'f65fee08b892af5ef93caa1bf3cb3980'


def test_insert_seqs_end():
    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', 'end')) == '792397e2e32e95b56ddc15b8b2310ec0'


def test_insert_seqs_index():
    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', 100)) == 'da2b2e0efb5807a51e925076857b189d'


def test_insert_seqs_endminus():
    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', -25)) == '0f4115d81cc5fa2cc381f17bada0f0ce'


def test_insert_seqs_startplus():
    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', 25)) == 'cb37efd3069227476306f9129efd4d05'


def test_insert_seqs_endminus_extreme():
    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', -9000)) == 'f65fee08b892af5ef93caa1bf3cb3980'


def test_insert_seqs_startplus_extreme():
    tester = Sb._make_copy(sb_objects[0])
    assert seqs_to_hash(Sb.insert_sequence(tester, 'AACAGGTCGAGCA', 9000)) == '792397e2e32e95b56ddc15b8b2310ec0'


# ######################  'ip', '--isoelectric_point' ###################### #
@pytest.mark.parametrize("seqbuddy", [Sb._make_copy(x) for x in [sb_objects[6], sb_objects[7], sb_objects[8],
                                                                 sb_objects[10], sb_objects[11]]])
def test_isoelectric_point(seqbuddy):
    output = Sb.isoelectric_point(Sb.clean_seq(seqbuddy))
    assert output[1]["Mle-Panxα12"] == 6.0117797852
    if seqbuddy.out_format == "gb":
        assert seqs_to_hash(seqbuddy) == "d3d22b310411419ad9383a83e0ab5893"


def test_isoelectric_point_type_error():
    with pytest.raises(TypeError):
        Sb.isoelectric_point(Sb.SeqBuddy(resource("/Mnemiopsis_cds.fa")))


# ##################### 'lf', 'list_features' ###################### ##
def test_list_features():
    tester = Sb._make_copy(sb_objects[1])
    output = Sb.list_features(tester)
    for record in tester.records:
        assert record.id in output
        if output[record.id] is not None:
            assert output[record.id] == record.features


# ######################  'li', '--list_ids' ###################### #
# first test that 1 column works for all file types
hashes = ["1c4a395d8aa3496d990c611c3b6c4d0a", "1c4a395d8aa3496d990c611c3b6c4d0a", "1c4a395d8aa3496d990c611c3b6c4d0a",
          "78a9289ab2d508a13c76cf9f5a308cc5", "1c4a395d8aa3496d990c611c3b6c4d0a", "1c4a395d8aa3496d990c611c3b6c4d0a"]
hashes = [(sb_objects[indx], value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_list_ids_one_col(seqbuddy, next_hash):
    tester = Sb.list_ids(seqbuddy, 1)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash

# Now test different numbers of columns
hashes = ["6fcee2c407bc4f7f70e0ae2a7e101761", "1c4a395d8aa3496d990c611c3b6c4d0a", "6fcee2c407bc4f7f70e0ae2a7e101761",
          "bd177e4db7dd772c5c42199b0dff49a5", "6b595a436a38e353a03e36a9af4ba1f9", "c57028374ed3fc474009e890acfb041e"]
columns = [-2, 0, 2, 5, 10, 100]
hashes = [(sb_objects[0], value, columns[indx]) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash,cols", hashes)
def test_list_ids_multi_col(seqbuddy, next_hash, cols):
    tester = Sb.list_ids(seqbuddy, cols)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash


# ######################  'uc', '--uppercase'  and 'lc', '--lowercase' ###################### #
uc_hashes = ["25073539df4a982b7f99c72dd280bb8f", "2e02a8e079267bd9add3c39f759b252c", "52e74a09c305d031fc5263d1751e265d",
             "7117732590f776836cbabdda05f9a982", "3d17ebd1f6edd528a153ea48dc37ce7d", "b82538a4630810c004dc8a4c2d5165ce",
             "c10d136c93f41db280933d5b3468f187", "7a8e25892dada7eb45e48852cbb6b63d", "8b6737fe33058121fd99d2deee2f9a76",
             "40f10dc94d85b32155af7446e6402dea", "b229db9c07ff3e4bc049cea73d3ebe2c", "f35cbc6e929c51481e4ec31e95671638"]

lc_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "2e02a8e079267bd9add3c39f759b252c", "cb1169c2dd357771a97a02ae2160935d",
             "d1524a20ef968d53a41957d696bfe7ad", "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
             "14227e77440e75dd3fbec477f6fd8bdc", "7a8e25892dada7eb45e48852cbb6b63d", "17ff1b919cac899c5f918ce8d71904f6",
             "c934f744c4dac95a7544f9a814c3c22a", "6a3ee818e2711995c95372afe073490b", "c0dce60745515b31a27de1f919083fe9"]

hashes = [(Sb._make_copy(sb_objects[indx]), uc_hash, lc_hashes[indx]) for indx, uc_hash in enumerate(uc_hashes)]


@pytest.mark.parametrize("seqbuddy,uc_hash,lc_hash", hashes)
def test_cases(seqbuddy, uc_hash, lc_hash):  # NOTE: Biopython always writes genbank to spec in lower case
    tester = Sb.uppercase(seqbuddy)
    assert seqs_to_hash(tester) == uc_hash
    tester = Sb.lowercase(tester)
    assert seqs_to_hash(tester) == lc_hash


# ######################  'fd2p', '--map_features_dna2prot' ###################### #
# Map the genbank DNA file to all protein files, and the fasta DNA file to fasta protein
hashes = ["5216ef85afec36d5282578458a41169a", "a8f7c129cf57a746c20198bf0a6b9cf4", "0deeea532d6dcbc0486e9b74d0d6aca8",
          "d595fabb157d5c996357b6a7058af4e8", "bb06e94456f99efc2068f5a52f0e0462", "a287e0054df7f5df76e792e0e0ab6756"]
prot_indx = [6, 7, 8, 9, 10, 11]
hashes = [(Sb._make_copy(sb_objects[1]), Sb._make_copy(sb_objects[prot_indx[indx]]), value)
          for indx, value in enumerate(hashes)]
hashes.append((Sb._make_copy(sb_objects[0]), Sb._make_copy(sb_objects[6]), "854566b485af0f277294bbfb15f7dd0a"))


@pytest.mark.parametrize("_dna,_prot,next_hash", hashes)
def test_map_features_dna2prot(_dna, _prot, next_hash):
    _prot.alpha = IUPAC.protein
    _dna.alpha = IUPAC.ambiguous_dna
    tester = Sb.map_features_dna2prot(_dna, _prot)
    assert seqs_to_hash(tester) == next_hash


# ######################  'fp2d', '--map_features_prot2dna' ###################### #
# Map the genbank protein file to all dna files, and the fasta protein file to fasta DNA
hashes = ["3ebc92ca11505489cab2453d2ebdfcf2", "6b4dd3fc66cb7419acaf064b589f4dd1",
          "8d403750ef83d60e31de0dee79a8f5d1", "74c6c4b5531c41f55f7349ed6c6b2f43",
          "9133ab0becbec95ce7ed31e02dc17ef5", "3ebc92ca11505489cab2453d2ebdfcf2"]
dna_indx = [0, 1, 2, 3, 4, 5]
hashes = [(Sb._make_copy(sb_objects[7]), Sb._make_copy(sb_objects[dna_indx[indx]]), value)
          for indx, value in enumerate(hashes)]
hashes.append((Sb._make_copy(sb_objects[6]), Sb._make_copy(sb_objects[0]), "720f36544f9c11855ac2673e63282f89"))


@pytest.mark.parametrize("_prot,_dna,next_hash", hashes)
def test_map_features_prot2dna(_prot, _dna, next_hash):
    tester = Sb.map_features_prot2dna(_prot, _dna)
    assert seqs_to_hash(tester) == next_hash


# #####################  'mg', '--merge' ###################### ##
def test_merge():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds_dummy_features.gb"))
    tester = Sb.merge(tester, Sb._make_copy(sb_objects[1]))
    assert seqs_to_hash(tester) == "5fb43a869bb9f766625d1336081766fc"

    tester.records[0].seq = tester.records[1].seq
    with pytest.raises(RuntimeError) as e:
        Sb.merge(tester, Sb._make_copy(sb_objects[1]))
    assert "Record mismatch: ID" in str(e.value)


# ######################  '-mw', '--molecular_weight' ###################### #
def test_molecular_weight():
    # Unambiguous DNA
    tester = Sb.molecular_weight(Sb._make_copy(sb_objects[1]))
    assert tester[1]['masses_ds'][0] == 743477.1
    assert tester[1]['masses_ss'][0] == 371242.6
    assert seqs_to_hash(tester[0]) == "e080cffef0ec6c5e8eada6f57bbc35f9"
    # Ambiguous DNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("ambiguous_dna.fa")))[1]
    assert tester['masses_ds'][0] == 743477.08
    assert tester['masses_ss'][0] == 371202.59
    # Unambiguous RNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("Mnemiopsis_rna.fa")))[1]
    assert tester['masses_ss'][0] == 387372.6
    # Ambiguous RNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(resource("ambiguous_rna.fa")))[1]
    assert tester['masses_ss'][0] == 387371.6
    # Protein
    tester = Sb.molecular_weight(Sb._make_copy(sb_objects[7]))
    assert tester[1]['masses_ss'][0] == 45692.99
    assert seqs_to_hash(tester[0]) == "fb1a66b7eb576c0584fc7988c45b6a18"


# ######################  'ns', '--num_seqs' ###################### #
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
              "cb1169c2dd357771a97a02ae2160935d", "d1524a20ef968d53a41957d696bfe7ad",
              "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "d0297078b4c480a49b6da5b719310d0e",
              "17ff1b919cac899c5f918ce8d71904f6", "c934f744c4dac95a7544f9a814c3c22a",
              "6a3ee818e2711995c95372afe073490b", "c0dce60745515b31a27de1f919083fe9"]

rev_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "3b718ec3cb794bcb658d900e517110cc",
              "cb1169c2dd357771a97a02ae2160935d", "d1524a20ef968d53a41957d696bfe7ad",
              "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "c6a788d8ea916964605ac2942c459c9b",
              "17ff1b919cac899c5f918ce8d71904f6", "c934f744c4dac95a7544f9a814c3c22a",
              "6a3ee818e2711995c95372afe073490b", "c0dce60745515b31a27de1f919083fe9"]
hashes = [(Sb._make_copy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)  # modifies in place?
def test_order_features_alphabetically(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_features_alphabetically(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_features_alphabetically(seqbuddy, reverse=True)
    assert seqs_to_hash(tester) == rev_hash


# ######################  'ofp', '--order_features_by_position' ###################### #
fwd_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "2e02a8e079267bd9add3c39f759b252c",
              "cb1169c2dd357771a97a02ae2160935d", "d1524a20ef968d53a41957d696bfe7ad",
              "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "7a8e25892dada7eb45e48852cbb6b63d",
              "17ff1b919cac899c5f918ce8d71904f6", "c934f744c4dac95a7544f9a814c3c22a",
              "6a3ee818e2711995c95372afe073490b", "c0dce60745515b31a27de1f919083fe9"]

rev_hashes = ["b831e901d8b6b1ba52bad797bad92d14", "4345a14fe27570b3c837c30a8cb55ea9",
              "cb1169c2dd357771a97a02ae2160935d", "d1524a20ef968d53a41957d696bfe7ad",
              "99d522e8f52e753b4202b1c162197459", "228e36a30e8433e4ee2cd78c3290fa6b",
              "14227e77440e75dd3fbec477f6fd8bdc", "9e7c2571db1386bba5983365ae235e1b",
              "17ff1b919cac899c5f918ce8d71904f6", "c934f744c4dac95a7544f9a814c3c22a",
              "6a3ee818e2711995c95372afe073490b", "c0dce60745515b31a27de1f919083fe9"]
hashes = [(Sb._make_copy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)  # modifies in place?
def test_order_features_position(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_features_by_position(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_features_by_position(seqbuddy, reverse=True)
    assert seqs_to_hash(tester) == rev_hash


# ######################  'oi', '--order_ids' ###################### #
fwd_hashes = ["ccc59629741421fb717b9b2403614c62", "2f9bc0dd9d79fd8160a621280be0b0aa",
              "60bbc6306cbb4eb903b1212718bb4592", "4188a065adb5b8e80acd3073afc1c7f9",
              "433520b63864022d82973f493dbf804b", "4078182a81382b815528fdd5c158fbec"]
rev_hashes = ["503a71fc2e8d143361cbe8f4611527fd", "dd269961d4d5301d1bf87e0093568851",
              "82fea6e3d3615ac75ec5022abce255da", "9d0910f3d303297283bace2718f60d61",
              "8af06c3523a1bf7cde4fc2b8c64a388c", "3b83a3c73a6cdded6635ffa10c4a16e1"]

hashes = [(Sb._make_copy(sb_objects[indx]), fwd_hash, rev_hashes[indx]) for indx, fwd_hash in enumerate(fwd_hashes)]


@pytest.mark.parametrize("seqbuddy,fwd_hash,rev_hash", hashes)
def test_order_ids(seqbuddy, fwd_hash, rev_hash):
    tester = Sb.order_ids(seqbuddy)
    assert seqs_to_hash(tester) == fwd_hash
    tester = Sb.order_ids(seqbuddy, reverse=True)
    assert seqs_to_hash(tester) == rev_hash


# #####################  'prr', '--pull_random_recs' ###################### ##
@pytest.mark.parametrize("seqbuddy", sb_objects)
def test_pull_random_recs(seqbuddy):
    tester = Sb.pull_random_recs(Sb._make_copy(seqbuddy))
    assert len(tester.records) == 1
    assert tester.records[0].id in Sb.list_ids(seqbuddy)


# #####################  'pre', '--pull_record_ends' ###################### ##
def test_pull_record_ends_front():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Sb.pull_record_ends(tester, 10, 'front')
    assert seqs_to_hash(tester) == '754d6868030d1122b35386118612db72'


def test_pull_record_ends_back():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Sb.pull_record_ends(tester, 10, 'rear')
    assert seqs_to_hash(tester) == '9cfc91c3fdc5cd9daabce0ef9bac2db7'


def test_pull_record_ends_zero():
    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Sb.pull_record_ends(Sb._make_copy(seqbuddy), 0, 'rear')
    assert seqs_to_hash(tester) == seqs_to_hash(seqbuddy)


def test_pull_record_ends_neg():
    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    with pytest.raises(ValueError):
        Sb.pull_record_ends(Sb._make_copy(seqbuddy), -1, 'rear')


def test_pull_record_ends_wrong_end():
    seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    with pytest.raises(AttributeError):
        Sb.pull_record_ends(Sb._make_copy(seqbuddy), 100, 'fghhgj')


# ######################  'pr', '--pull_recs' ###################### #
pr_hashes = ["5b4154c2662b66d18776cdff5af89fc0", "e196fdc5765ba2c47f97807bafb6768c", "bc7dbc612bc8139eba58bf896b7eaf2f",
             "e33828908fa836f832ee915957823039", "e33828908fa836f832ee915957823039", "b006b40ff17ba739929448ae2f9133a6"]
pr_hashes = [(Sb.SeqBuddy(resource(seq_files[indx])), value) for indx, value in enumerate(pr_hashes)]


@pytest.mark.parametrize("seqbuddy, next_hash", pr_hashes)
def test_pull_recs(seqbuddy, next_hash):
    tester = Sb.pull_recs(seqbuddy, 'α2')
    assert seqs_to_hash(tester) == next_hash


# #####################  'prg', '--purge' ###################### ##
def test_purge():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_pep.fa"))
    tester = Sb.purge(tester, 200)
    assert seqs_to_hash(tester[0]) == 'b21b2e2f0ca1fcd7b25efbbe9c08858c'


# ######################  'rs', '--raw_seq' ###################### #
hashes = ["6f0ff2d43706380d92817e644e5b78a5", "5d00d481e586e287f32d2d29916374ca", "6f0ff2d43706380d92817e644e5b78a5",
          "cda59127d6598f44982a2d1875064bb1", "6f0ff2d43706380d92817e644e5b78a5", "6f0ff2d43706380d92817e644e5b78a5",
          "cdfe71aefecc62c5f5f2f45e9800922c", "4dd913ee3f73ba4bb5dc90d612d8447f", "cdfe71aefecc62c5f5f2f45e9800922c",
          "3f48f81ab579a389947641f36889901a", "cdfe71aefecc62c5f5f2f45e9800922c", "cdfe71aefecc62c5f5f2f45e9800922c"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_raw_seq(seqbuddy, next_hash):
    tester = Sb.raw_seq(seqbuddy)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash


# ######################  'ri', '--rename_ids' ###################### #
hashes = ["59bea136d93d30e3f11fd39d73a9adff", "78c73f97117bd937fd5cf52f4bd6c26e", "243024bfd2f686e6a6e0ef65aa963494",
          "83f10d1be7a5ba4d363eb406c1c84ac7", "973e3d7138b78db2bb3abda8a9323226", "4289f03afb6c9f8a8b0d8a75bb60a2ce"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_rename_ids(seqbuddy, next_hash):
    tester = Sb.rename(seqbuddy, 'Panx', 'Test', 0)
    assert seqs_to_hash(tester) == next_hash


# ######################  'rc', '--reverse_complement' ###################### #
hashes = ["e77be24b8a7067ed54f06e0db893ce27", "47941614adfcc5bd107f71abef8b3e00", "f549c8dc076f6b3b4cf5a1bc47bf269d",
          "a62edd414978f91f7391a59fc1a72372", "08342be5632619fd1b1251b7ad2b2c84", "0d6b7deda824b4fc42b65cb87e1d4d14"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_reverse_complement(seqbuddy, next_hash):
    tester = Sb.reverse_complement(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_reverse_complement_pep_exception():  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.reverse_complement(sb_objects[6])


# ######################  'sfr', '--select_frame' ###################### #
# Only fasta
hashes = ["b831e901d8b6b1ba52bad797bad92d14", "a518e331fb29e8be0fdd5f3f815f5abb", "2cbe39bea876030da6d6bd45e514ae0e"]
frame = [1, 2, 3]
hashes = [(Sb._make_copy(sb_objects[0]), _hash, frame[indx]) for indx, _hash in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash,shift", hashes)
def test_select_frame(seqbuddy, next_hash, shift):
    tester = Sb.select_frame(seqbuddy, shift)
    assert seqs_to_hash(tester) == next_hash


def test_select_frame_pep_exception():
    with pytest.raises(TypeError):  # If protein is input
        Sb.select_frame(sb_objects[6], 2)


# ##################### 'ss', 'shuffle_seqs' ###################### ##
def test_shuffle_seqs():
    tester1 = Sb._make_copy(sb_objects[0])
    tester2 = Sb._make_copy(tester1)
    Sb.shuffle_seqs(tester2)
    assert seqs_to_hash(tester1) != seqs_to_hash(tester2)
    for indx, record in enumerate(tester1.records):
        assert sorted(record.seq) == sorted(tester2.records[indx].seq)


# ######################  'sh', '--shuffle' ###################### #
@pytest.mark.parametrize("seqbuddy", [Sb._make_copy(x) for x in sb_objects])
def test_shuffle(seqbuddy):
    tester = Sb.order_ids_randomly(Sb._make_copy(seqbuddy))
    for i in range(3):  # Sometimes shuffle doesn't actually shuffle, so repeat a few times if necessary
        if seqs_to_hash(seqbuddy) != seqs_to_hash(tester):
            break
        tester = Sb.order_ids_randomly(Sb._make_copy(seqbuddy))

    assert seqs_to_hash(seqbuddy) != seqs_to_hash(tester)
    assert seqs_to_hash(Sb.order_ids(tester)) == seqs_to_hash(tester)


# #####################  'sbt', '--split_by_taxa' ###################### ##
def test_split_by_taxa():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Sb.split_by_taxa(tester, 'Mle-')
    output = md5(str(sorted(tester)).encode()).hexdigest()
    assert output == '438cc47dd268453f38342ca276f77a73'


# #####################  'sf', '--split_file' ###################### ##
def test_split_file():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    output = Sb.split_file(tester)
    for buddy in output:
        assert buddy.records[0] in tester.records


# ######################  'tr6', '--translate6frames' ###################### #
# Only fasta and genbank
hashes = ["d5d39ae9212397f491f70d6928047341", "42bb6caf86d2d8be8ab0defabc5af477"]
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_translate6frames(seqbuddy, next_hash):
    tester = Sb.translate6frames(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_translate6frames_pep_exception():
    with pytest.raises(TypeError):
        Sb.translate6frames(Sb._make_copy(sb_objects[6]))


# ######################  'tr', '--translate' ###################### #
hashes = ["3de7b7be2f2b92cf166b758625a1f316", "c841658e657b4b21b17e4613ac27ea0e", ]
# NOTE: the first 6 sb_objects are DNA.
hashes = [(Sb._make_copy(sb_objects[indx]), value) for indx, value in enumerate(hashes)]


@pytest.mark.parametrize("seqbuddy,next_hash", hashes)
def test_translate(seqbuddy, next_hash):
    tester = Sb.translate_cds(seqbuddy)
    assert seqs_to_hash(tester) == next_hash


def test_translate_pep_exception():
    with pytest.raises(TypeError):
        Sb.translate_cds(Sb._make_copy(sb_objects[6]))


# ################################################# COMMAND LINE UI ################################################## #
# ##################### 'ano', '--annotate' ###################### ##
def test_annotate_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.annotate = ["misc_feature", "1-100,200-250", "+"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "4a294c1ea5783fe04e909d767e7dfce2"

    test_in_args = deepcopy(in_args)
    test_in_args.annotate = ["misc_feature", "1-100,200-250", "foo=bar", "hello=world", "-", "α4"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "42ca93bcf3103124287fffbde49f2b26"

    test_in_args = deepcopy(in_args)
    test_in_args.annotate = ["unknown_feature_that_is_t0o_long", "1-100,200-250", "foo=bar", "hello=world", "-", "α4"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Warning: The provided annotation type is not part of the GenBank format standard" in err
    assert "Warning: Feature type is longer than 16 characters" in err


# ######################  'asl', '--ave_seq_length' ###################### #
def test_ave_seq_length_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.ave_seq_length = [False]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert out == '428.38\n'

    test_in_args.ave_seq_length = ['clean']
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert out == '427.38\n'


# ######################  'r2d', '--back_transcribe' ###################### #
def test_back_transcribe_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.back_transcribe = True
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("Mnemiopsis_rna.fa")), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b831e901d8b6b1ba52bad797bad92d14"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]))

    out, err = capsys.readouterr()
    assert "You need to provide an RNA sequence." in err


# ######################  'btr', '--back_translate' ###################### #
def test_back_translate_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.back_translate = [False]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert "Panxα4" in out

    test_in_args = deepcopy(in_args)
    test_in_args.back_translate = [["human", "o"]]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[7]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b6bcb4e5104cb202db0ec4c9fc2eaed2"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]))

    out, err = capsys.readouterr()
    assert err == "TypeError: The input sequence needs to be protein, not nucleotide\n"


# ######################  'bl2s', '--bl2seq' ###################### #
def test_bl2s_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.bl2seq = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "339377aee781fb9d01456f04553e3923"

    os.environ["PATH"] = ""
    with mock.patch('builtins.input', return_value="n"):
        with pytest.raises(SystemExit):
            Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]))
    out, err = capsys.readouterr()
    assert "not present in $PATH or working directory" in err
    os.environ["PATH"] = BACKUP_PATH


# ######################  'bl', '--blast' ###################### #
def test_blast_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.blast = resource("blast/Mnemiopsis_cds")
    tester = Sb._make_copy(sb_objects[0])
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "a56ec76a64b25b7ca8587c7aa8554412"

    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("blank.fa")), True)
    out, err = capsys.readouterr()
    assert out == "No significant matches found\n"

    with pytest.raises(SystemExit):
        test_in_args.blast = resource("./Mnemiopsis_cds")
        Sb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "RuntimeError:" in err


# ######################  'cs', '--clean_seq' ###################### #
def test_clean_seq_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.clean_seq = [True]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "dc53f3be7a7c24425dddeea26ea0ebb5"

    test_in_args.clean_seq = ["strict"]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("ambiguous_dna.fa")), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "5fd0b78e37c81e0fa727db34a37cc743"


# ######################  'cmp', '--complement' ###################### #
def test_complement_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.complement = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "e4a358ca57aca0bbd220dc6c04c88795"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]))

    out, err = capsys.readouterr()
    assert "Nucleic acid sequence required, not protein." in err


# ######################  'cts', '--concat_seqs' ###################### #
def test_concat_seqs_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.concat_seqs = [True]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "7421c27be7b41aeedea73ff41869ac47"

    test_in_args.concat_seqs = ["clean"]
    test_in_args.out_format = "embl"
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[2]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "15b9c79dea034cef74e3a622bd357705"


# ######################  'cc', '--count_codons' ###################### #
def test_count_codons_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.count_codons = ["foo"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "5661f0ce92bb6cfba4519a61e0a838ed"

    test_in_args.count_codons = ["conc"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "3e76bd510de4a61efb17ffc186ef9e68"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[7]))
    out, err = capsys.readouterr()
    assert "Nucleic acid sequence required, not protein" in err


# ######################  'cr', '--count_residues' ###################### #
def test_count_residues_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.count_residues = ["foo"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "c9bc54835fb54232d4d346d04344bf8b"

    test_in_args.count_residues = ["conc"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "c7062408f939f4b310f2f97c3e94eb37"


# ######################  'df', '--delete_features' ###################### #
def test_delete_features_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_features = ["donor"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "f84df6a77063c7def13babfaa0555bbf"


# ######################  'dlg', '--delete_large' ###################### #
def test_delete_large_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_large = 1285
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "25859dc69d46651a1e04a70c07741b35"


# ######################  'dm', '--delete_metadata' ###################### #
def test_delete_metadata_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_metadata = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "544ab887248a398d6dd1aab513bae5b1"


# ######################  'dr', '--delete_records' ###################### #
def test_delete_records_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_records = ['α1']
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "54e810265a6ecf7a3d140fc806597f93"
    assert string2hash(err) == "4f420c9128e515dc24031b5075c034e3"

    test_in_args.delete_records = ['α1', 'α2']
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "eca4f181dae3d7998464ff71e277128f"
    assert string2hash(err) == "b3983fac3c2cf15f83650a34a17151da"

    test_in_args.delete_records = ['α1', 'α2', "3"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "eca4f181dae3d7998464ff71e277128f"
    assert string2hash(err) == "7e0929af515502484feb4b1b2c35eaba"

    test_in_args.delete_records = ['foo']
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b831e901d8b6b1ba52bad797bad92d14"
    assert string2hash(err) == "553348fa37d9c67f4ce0c8c53b578481"


# ######################  'drp', '--delete_repeats' ###################### #
def test_delete_repeats_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_repeats = [None]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("Mnemiopsis_dup_id.fa")), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "03d257a8242469173fd3a314f8fecf92"
    assert string2hash(err) == "f05dba6e411e011bd684a13c72468a53"

    test_in_args.delete_repeats = [None]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(resource("Mnemiopsis_dup_seq.fa")), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "5f101b52693ab391ed85f9e4375a16f5"
    assert string2hash(err) == "db2b2fd88be453216d47325cf5c6aa35"

    test_in_args.delete_repeats = [[2, "all"]]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "b831e901d8b6b1ba52bad797bad92d14"
    assert err == "No duplicate records found\n"

    test_in_args.quiet = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert not err


# ######################  'ds', '--delete_small' ###################### #
def test_delete_small_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_small = 1285
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "196adf08d4993c51050289e5167dacdf"


# ######################  'er', '--extract_region' ###################### #
def test_extract_region_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.extract_region = [50, 300]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "3e791c6a6683516aff9572c24f38f0b3"

    with pytest.raises(SystemExit):
        test_in_args.extract_region = [300, 50]
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]))
    out, err = capsys.readouterr()
    assert "The value given for end of range is smaller than for" in err


# ######################  'fcpg', '--find_cpg' ###################### #
def test_find_cpg_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.find_CpG = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "ce6aff066c03651401db627951862154"
    assert string2hash(err) == "1fc07e5884a2a4c1344865f385b1dc79"

    Sb.command_line_ui(test_in_args, Sb.SeqBuddy(">seq1\nATGCCTAGCTAGCT", in_format="fasta"), True)
    out, err = capsys.readouterr()
    assert err == "# No Islands identified\n\n"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[7]))
    out, err = capsys.readouterr()
    assert "DNA sequence required, not protein or RNA" in err


# ######################  'fp', '--find_pattern' ###################### #
def test_find_pattern_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.find_pattern = ["ATg{2}T", "tga.{1,6}tg"]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[1]), True)
    out, err = capsys.readouterr()

    assert string2hash(out) == "0e11b2c0e9451fbfcbe39e3b5be2cf60"
    assert string2hash(err) == "f2bb95f89e7b9e198f18a049afbe4a93"


# ######################  'gf', '--guess_format' ###################### #
def test_guess_format_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.guess_format = True
    paths = ["Mnemiopsis_cds.%s" % x for x in ["embl", "fa", "gb", "nex", "phy", "phyr", "seqxml", "stklm"]]
    paths.append("gibberish.fa")
    test_in_args.sequence = [resource(x) for x in paths]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "1fbee52b607dea0c8c767d9b3bb7e8b1"

    text_io = io.open(resource("Mnemiopsis_cds.embl"), "r")
    test_in_args.sequence = [text_io]
    tester = Sb.SeqBuddy(text_io)
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert out == "embl\n"

    tester.in_format = None
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert out == "Unknown\n"


# ######################  'hsi', '--hash_seq_ids' ###################### #
def test_hash_seq_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.hash_seq_ids = [None]
    tester = Sb._make_copy(sb_objects[0])
    ids = [rec.id for rec in tester.records]
    Sb.command_line_ui(test_in_args, tester, True)
    for indx, rec in enumerate(tester.records):
        assert rec.id != ids[indx]
        assert ids[indx] == tester.hash_map[rec.id]

    test_in_args.hash_seq_ids = [0]
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert "Warning: The hash_length parameter was passed in with the value 0. This is not a positive integer" in err

    tester.records *= 10
    test_in_args.hash_seq_ids = [1]
    Sb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert "cover all sequences, so it has been increased to 2" in err


# ######################  'mg', '--merge' ###################### #
def test_merge_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.merge = True
    test_in_args.sequence = [resource("Mnemiopsis_cds_dummy_features.gb"), resource("Mnemiopsis_cds.gb")]
    Sb.command_line_ui(test_in_args, Sb.SeqBuddy, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "5fb43a869bb9f766625d1336081766fc"

    with pytest.raises(SystemExit):
        test_in_args.sequence = [resource("Mnemiopsis_pep.gb"), resource("Mnemiopsis_cds.gb")]
        Sb.command_line_ui(test_in_args, Sb.SeqBuddy)
    out, err = capsys.readouterr()
    assert "RuntimeError" in err


# ######################  'd2r', '--transcribe' ###################### #
def test_transcribe_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.transcribe = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "d2db9b02485e80323c487c1dd6f1425b"

    with pytest.raises(SystemExit):
        Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[6]))

    out, err = capsys.readouterr()
    assert "You need to provide a DNA sequence." in err


# ######################  '-ofa', '--order_features_alphabetically' ###################### #
def test_order_features_alphabetically_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.order_features_alphabetically = [True]
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    tester = Sb.SeqBuddy(out)
    assert seqs_to_hash(tester) == 'b831e901d8b6b1ba52bad797bad92d14'


# ######################  'ns', '--num_seqs' ###################### #
def test_num_seqs_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.num_seqs = True
    Sb.command_line_ui(test_in_args, Sb._make_copy(sb_objects[0]), True)
    out, err = capsys.readouterr()
    assert out == '13\n'
