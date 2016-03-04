#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: alignbuddy_tests.py
version: 1.0
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Collection of PyTest unit tests for the AlignBuddy.py program
"""

import pytest
from hashlib import md5
import os
import sys
import argparse
import io
from copy import deepcopy
from collections import OrderedDict
from subprocess import Popen, PIPE
from unittest import mock

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

sys.path.insert(0, "./")
import MyFuncs
import AlignBuddy as Alb
import SeqBuddy as Sb
import buddy_resources as br

TEMP_DIR = MyFuncs.TempDir()
VERSION = Sb.VERSION


def fmt(prog):
    return br.CustomHelpFormatter(prog)

parser = argparse.ArgumentParser(prog="alignBuddy", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                 description='''\
\033[1mAlignBuddy\033[m
  Sequence alignments with a splash of kava.

\033[1mUsage examples\033[m:
  AlignBuddy.py "/path/to/align_file" -<cmd>
  AlignBuddy.py "/path/to/align_file" -<cmd> | AlignBuddy.py -<cmd>
  AlignBuddy.py "/path/to/seq_file" -ga "mafft" -p "--auto --thread 8"
''')

br.flags(parser, ("alignments", "The file(s) you want to start working on"),
         br.alb_flags, br.alb_modifiers, VERSION)

# This is to allow py.test to work with the -x flag
parser.add_argument("-x", nargs="?")
parser.add_argument("-m", nargs="?")
parser.add_argument("-n", nargs="?")
parser.add_argument("--cov", nargs="?")
parser.add_argument("--cov-report", nargs="?")
in_args = parser.parse_args()


def align_to_hash(_alignbuddy, mode='hash'):
    if mode != "hash":
        return str(_alignbuddy)
    _hash = md5("{0}".format(str(_alignbuddy)).encode()).hexdigest()
    return _hash

root_dir = os.getcwd()


def string2hash(_input):
    return md5(_input.encode()).hexdigest()


def resource(file_name):
    return "{0}/unit_test_resources/{1}".format(root_dir, file_name)


# ################### Alignment resources ################### #
class Resources(object):
    def __init__(self):
        one_dna = OrderedDict([("clustal", "Mnemiopsis_cds.clus"),
                               ("fasta", "Mnemiopsis_cds_aln.fa"),
                               ("gb", "Mnemiopsis_cds_aln.gb"),
                               ("nexus", "Mnemiopsis_cds.nex"),
                               ("phylip", "Mnemiopsis_cds.phy"),
                               ("phylipr", "Mnemiopsis_cds.phyr"),
                               ("phylipss", "Mnemiopsis_cds.physs"),
                               ("phylipsr", "Mnemiopsis_cds.physr"),
                               ("stockholm", "Mnemiopsis_cds.stklm")])
        one_rna = OrderedDict([("nexus", "Mnemiopsis_rna.nex")])
        one_pep = OrderedDict([("gb", "Mnemiopsis_pep_aln.gb"),
                               ("nexus", "Mnemiopsis_pep.nex"),
                               ("phylip", "Mnemiopsis_pep.phy"),
                               ("phylipr", "Mnemiopsis_pep.phyr"),
                               ("phylipss", "Mnemiopsis_pep.physs"),
                               ("phylipsr", "Mnemiopsis_pep.physr"),
                               ("stockholm", "Mnemiopsis_pep.stklm")])

        one_aligns = OrderedDict([("dna", one_dna),
                                  ("rna", one_rna),
                                  ("pep", one_pep)])

        multi_dna = OrderedDict([("clustal", "Alignments_cds.clus"),
                                 ("phylip", "Alignments_cds.phy"),
                                 ("phylipr", "Alignments_cds.phyr"),
                                 ("phylipss", "Alignments_cds.physs"),
                                 ("phylipsr", "Alignments_cds.physr"),
                                 ("stockholm", "Alignments_cds.stklm")])
        multi_pep = OrderedDict([("clustal", "Alignments_pep.clus"),
                                 ("phylip", "Alignments_pep.phy"),
                                 ("phylipr", "Alignments_pep.phyr"),
                                 ("phylipss", "Alignments_pep.physs"),
                                 ("phylipsr", "Alignments_pep.physr"),
                                 ("stockholm", "Alignments_pep.stklm")])

        multi_aligns = OrderedDict([("dna", multi_dna),
                                    ("pep", multi_pep)])

        self.resources = OrderedDict([("one", one_aligns), ("multi", multi_aligns)])

        self.alb_objs = OrderedDict()
        self.res_paths = OrderedDict()
        for num in self.resources:
            self.alb_objs.setdefault(num, OrderedDict())
            self.res_paths.setdefault(num, OrderedDict())
            for _type in self.resources[num]:
                self.alb_objs[num][_type] = OrderedDict([(key, Alb.AlignBuddy(resource(path)))
                                                         for key, path in self.resources[num][_type].items()])
                self.res_paths[num][_type] = OrderedDict([(key, resource(path))
                                                         for key, path in self.resources[num][_type].items()])

        self.code_dict = OrderedDict([("num_aligns", OrderedDict([("o", "one"), ("m", "multi")])),
                                      ("type", OrderedDict([("p", "pep"), ("d", "dna"), ("r", "rna")])),
                                      ("format", OrderedDict([("c", "clustal"), ("f", "fasta"), ("g", "gb"),
                                                              ("n", "nexus"), ("py", "phylip"), ("pr", "phylipr"),
                                                              ("pss", "phylipss"), ("psr", "phylipsr"),
                                                              ("s", "stockholm")]))])

    def _parse_code(self, code=""):
        results = OrderedDict([("num_aligns", []), ("type", []), ("format", [])])
        code = code.split()
        for i in code:
            for j in results:
                if i in self.code_dict[j]:
                    results[j].append(i)

        # Fill up a field with all possibilities if nothing is given
        results["num_aligns"] = [key for key in self.code_dict["num_aligns"]] \
            if not results["num_aligns"] else results["num_aligns"]
        results["type"] = [key for key in self.code_dict["type"]] if not results["type"] else results["type"]
        results["format"] = [key for key in self.code_dict["format"]] if not results["format"] else results["format"]
        return results

    def get(self, code="", mode="objs"):
        """
        Returns copies of AlignBuddy objects, the
        :param code:
        :param mode: {"objs", "paths"}
        :return: OrderedDict {key: resource}
        """
        files = self._parse_code(code)
        output = OrderedDict()
        key = ["", "", ""]
        for num_aligns in files["num_aligns"]:
            key[0] = num_aligns
            n = self.code_dict["num_aligns"][num_aligns]
            for _type in files["type"]:
                key[1] = _type
                t = self.code_dict["type"][_type]
                for _format in files["format"]:
                    key[2] = _format
                    f = self.code_dict["format"][_format]
                    try:
                        assert not " ".join(key) in output
                        if mode == "objs":
                            output[" ".join(key)] = Alb.make_copy(self.alb_objs[n][t][f])
                        elif mode == "paths":
                            output[" ".join(key)] = self.res_paths[n][t][f]
                        else:
                            raise ValueError("The 'mode' parameter only accepts 'objs' or 'paths' as input.")
                    except KeyError:
                        pass
        return output

    def get_list(self, code="", mode="objs"):
        return [value for key, value in self.get(code=code, mode=mode).items()]

    def get_one(self, code, mode="objs"):
        output = self.get_list(code, mode)
        return None if not output or len(output) > 1 else output[0]

    def deets(self, code):
        code = code.split()
        return {"num_aligns": self.code_dict["num_aligns"][code[0]],
                "type": self.code_dict["type"][code[1]],
                "format": br.parse_format(self.code_dict["format"][code[2]])}


alb_resources = Resources()


@pytest.mark.parametrize("key, align_file", alb_resources.get(mode="paths").items())
def test_instantiate_alignbuddy_from_file(key, align_file):
    key = key.split()
    assert type(Alb.AlignBuddy(align_file, in_format=alb_resources.code_dict["format"][key[2]])) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", alb_resources.get_list(mode="paths"))
def test_instantiate_alignbuddy_from_file_guess(align_file):
    assert type(Alb.AlignBuddy(align_file)) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", alb_resources.get_list(mode="paths"))
def test_instantiate_alignbuddy_from_handle(align_file):
    with open(align_file, 'r') as ifile:
        assert type(Alb.AlignBuddy(ifile)) == Alb.AlignBuddy


@pytest.mark.parametrize("align_file", alb_resources.get_list(mode="paths"))
def test_instantiate_alignbuddy_from_raw(align_file):
    with open(align_file, 'r') as ifile:
        assert type(Alb.AlignBuddy(ifile.read())) == Alb.AlignBuddy


@pytest.mark.parametrize("alignbuddy", alb_resources.get_list(mode="objs"))
def test_instantiate_alignbuddy_from_alignbuddy(alignbuddy):
    assert type(Alb.AlignBuddy(alignbuddy)) == Alb.AlignBuddy


@pytest.mark.parametrize("alignbuddy", alb_resources.get_list(mode="objs"))
def test_instantiate_alignbuddy_from_list(alignbuddy):
    assert type(Alb.AlignBuddy(alignbuddy.alignments)) == Alb.AlignBuddy

    with pytest.raises(TypeError):  # When non-MultipleSeqAlignment objects are in the .alignments list
        alignbuddy.alignments.append("Dummy string object")
        Alb.AlignBuddy(alignbuddy.alignments)


def test_instantiation_alignbuddy_errors():
    with pytest.raises(br.GuessError) as e:
        Alb.AlignBuddy(resource("gibberish.fa"))
    assert "Could not determine format from _input file" in str(e)

    tester = open(resource("gibberish.fa"), "r")
    with pytest.raises(br.GuessError) as e:
        Alb.AlignBuddy(tester.read())
    assert "Could not determine format from raw" in str(e)

    tester.seek(0)
    with pytest.raises(br.GuessError) as e:
        Alb.AlignBuddy(tester)
    assert "Could not determine format from input file-like object" in str(e)


def test_empty_file():
    with open(resource("blank.fa"), "r") as ifile:
        with pytest.raises(br.GuessError) as e:
            Alb.AlignBuddy(ifile)
        assert "Empty file" in str(e)


# ##################### AlignBuddy methods ###################### ##
def test_set_format():
    tester = alb_resources.get_list("o d g")[0]
    tester.set_format("fasta")
    assert tester._out_format == "fasta"


def test_records():
    tester = alb_resources.get_list("m p py")[0]
    assert len(tester.records()) == 29


def test_records_iter():
    tester = alb_resources.get_list("m p py")[0]
    counter = 0
    for rec in tester.records_iter():
        assert type(rec) == SeqRecord
        counter += 1
    assert counter == 29


hashes = {'o p g': 'bf8485cbd30ff8986c2f50b677da4332', 'o p n': '17ff1b919cac899c5f918ce8d71904f6',
          'o p py': '968ed9fa772e65750f201000d7da670f', 'o p pr': 'ce423d5b99d5917fbef6f3b47df40513',
          'o p pss': "4bd927145de635c429b2917e0a1db176", 'o p psr': "8ff80c7f0b8fc7f237060f94603c17be",
          'o p s': 'c0dce60745515b31a27de1f919083fe9',

          'o d c': '3c937c9fec251a42f0994caabb64420c', 'o d f': '98a3a08389284461ea9379c217e99770',
          'o d g': '2a42c56df314609d042bdbfa742871a3', 'o d n': 'cb1169c2dd357771a97a02ae2160935d',
          'o d py': '503e23720beea201f8fadf5dabda75e4', 'o d pr': '52c23bd793c9761b7c0f897d3d757c12',
          'o d pss': '4c0c1c0c63298786e6fb3db1385af4d5', 'o d psr': 'c5fb6a5ce437afa1a4004e4f8780ad68',
          'o d s': '228e36a30e8433e4ee2cd78c3290fa6b',

          'o r n': 'f3bd73151645359af5db50d2bdb6a33d',

          'm p c': 'f0e20a55f679ee492bb0b3be444b46f9', 'm p s': '3fd5805f61777f7f329767c5f0fb7467',
          'm p py': '2a77f5761d4f51b88cb86b079e564e3b', 'm p pr': '3fef9a05058a5259ebd517d1500388d4',
          'm p pss': 'eb82cda31fcb2cf00e11d7e910fde695', 'm p psr': 'a16c6e617e5a88fef080eea54e54e8a8',

          'm d c': '058ef1525cfc1364f26dd5a5bd6b97fb', 'm d s': 'ae352b908be94738d6d9cd54770e5b5d',
          'm d py': '42679a32ebd93b628303865f68b0293d', 'm d pr': '22c0f0c8f014a34be8edd394bf477a2d',
          'm d pss': 'c789860da8f0b59e0adc7bde6342b4b0', 'm d psr': '28b2861275e0a488042cff35393ac36d'}

albs = [(hashes[key], alignbuddy) for key, alignbuddy in alb_resources.get("").items()]


@pytest.mark.parametrize("next_hash,alignbuddy", albs)
def test_print(next_hash, alignbuddy, capsys):
    alignbuddy.print()
    out, err = capsys.readouterr()
    tester = string2hash(out[:-1])
    assert tester == next_hash


@pytest.mark.parametrize("next_hash,alignbuddy", albs)
def test_str(next_hash, alignbuddy):
    tester = str(alignbuddy)
    tester = string2hash(tester)
    assert tester == next_hash


@pytest.mark.parametrize("next_hash,alignbuddy", albs)
def test_write1(next_hash, alignbuddy):
    temp_file = MyFuncs.TempFile()
    alignbuddy.write(temp_file.path)
    out = temp_file.read()
    tester = string2hash(out)
    assert tester == next_hash


def test_write2():  # Unloopable components
    tester = alb_resources.get_one("m p py")
    tester.set_format("fasta")
    with pytest.raises(ValueError):
        str(tester)

    tester.alignments = []
    assert str(tester) == "AlignBuddy object contains no alignments.\n"

    tester = alb_resources.get_one("o d pr")
    tester.set_format("phylipi")
    assert align_to_hash(tester) == "52c23bd793c9761b7c0f897d3d757c12"

    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds_hashed_ids.nex"))
    tester.set_format("phylip-strict")
    assert align_to_hash(tester) == "16b3397d6315786e8ad8b66e0d9c798f"


# ################################################# HELPER FUNCTIONS ################################################# #
def test_guess_error():
    # File path
    with pytest.raises(br.GuessError):
        Alb.AlignBuddy(resource("unrecognizable.txt"))

    with open(resource("unrecognizable.txt"), 'r') as ifile:
        # Raw
        with pytest.raises(br.GuessError) as e:
            Alb.AlignBuddy(ifile.read())
        assert "Could not determine format from raw input" in str(e)

        # Handle
        with pytest.raises(br.GuessError) as e:
            ifile.seek(0)
            Alb.AlignBuddy(ifile)
        assert "Could not determine format from input file-like object" in str(e)

    # GuessError output
    try:
        Alb.AlignBuddy(resource("unrecognizable.txt"))
    except br.GuessError as e:
        assert "Could not determine format from _input file" in str(e) and \
               "\nTry explicitly setting with -f flag." in str(e)


def test_guess_alphabet():
    for alb in alb_resources.get_list("d"):
        assert Alb.guess_alphabet(alb) == IUPAC.ambiguous_dna
    for alb in alb_resources.get_list("p"):
        assert Alb.guess_alphabet(alb) == IUPAC.protein
    for alb in alb_resources.get_list("r"):
        assert Alb.guess_alphabet(alb) == IUPAC.ambiguous_rna

    assert not Alb.guess_alphabet(Alb.AlignBuddy("", in_format="fasta"))


def test_guess_format():
    assert Alb.guess_format(["dummy", "list"]) == "stockholm"

    for key, obj in alb_resources.get().items():
        assert Alb.guess_format(obj) == alb_resources.deets(key)["format"]

    for key, path in alb_resources.get(mode="paths").items():
        assert Alb.guess_format(path) == alb_resources.deets(key)["format"]
        with open(path, "r") as ifile:
            assert Alb.guess_format(ifile) == alb_resources.deets(key)["format"]
            ifile.seek(0)
            string_io = io.StringIO(ifile.read())
        assert Alb.guess_format(string_io) == alb_resources.deets(key)["format"]

    Alb.guess_format(resource("blank.fa")) == "empty file"
    assert not Alb.guess_format(resource("malformed_phylip_records.physs"))
    assert not Alb.guess_format(resource("malformed_phylip_columns.physs"))

    with pytest.raises(br.GuessError) as e:
        Alb.guess_format({"Dummy dict": "Type not recognized by guess_format()"})
    assert "Unsupported _input argument in guess_format()" in str(e)


def test_parse_format():
    for _format in ["phylip", "phylipis", "phylip-strict", "phylip-interleaved-strict"]:
        assert br.parse_format(_format) == "phylip"

    for _format in ["phylipi", "phylip-relaxed", "phylip-interleaved", "phylipr"]:
        assert br.parse_format(_format) == "phylip-relaxed"

    for _format in ["phylips", "phylipsr", "phylip-sequential", "phylip-sequential-relaxed"]:
        assert br.parse_format(_format) == "phylipsr"

    for _format in ["phylipss", "phylip-sequential-strict"]:
        assert br.parse_format(_format) == "phylipss"

    with pytest.raises(TypeError) as e:
        br.parse_format("foo")
    assert "Format type 'foo' is not recognized/supported" in str(e)


def test_make_copy():
    for alb in alb_resources.get_list():
        tester = Alb.make_copy(alb)
        align_to_hash(tester) == align_to_hash(alb)


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


# ################################################ MAIN API FUNCTIONS ################################################ #
# ##########################################  '-al', '--alignment_lengths' ########################################### #
def test_alignment_lengths():
    lengths = Alb.alignment_lengths(alb_resources.get_one("m p c"))
    assert lengths[0] == 481
    assert lengths[1] == 683

    lengths = Alb.alignment_lengths(alb_resources.get_one("m d s"))
    assert lengths[0] == 2043
    assert lengths[1] == 1440


# ##############################################  '-cs', '--clean_seqs' ############################################## #
def test_clean_seqs():
    # Test an amino acid file
    tester = Alb.clean_seq(alb_resources.get_one("m p py"))
    assert align_to_hash(tester) == "07a861a1c80753e7f89f092602271072"

    tester = Alb.clean_seq(Alb.AlignBuddy(resource("ambiguous_dna_alignment.fa")), ambiguous=False, rep_char="X")
    assert align_to_hash(tester) == "6755ea1408eddd0e5f267349c287d989"


# ###########################################  '-cta', '--concat_alignments' ######################################### #
def test_concat_alignments():
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

    tester = Sb.SeqBuddy(resource("Cnidaria_pep.nexus"))
    Sb.pull_recs(tester, "Ccr|Cla|Hec")
    tester = Alb.AlignBuddy(str(tester))
    tester.alignments.append(tester.alignments[0])
    assert align_to_hash(Alb.concat_alignments(Alb.make_copy(tester))) == '32a507107b7dcd044ea7760c8812441c'

    tester.set_format("gb")
    assert align_to_hash(Alb.concat_alignments(Alb.make_copy(tester),
                                               "(.).(.)-Panx(.)")) == '5ac908ebf7918a45664a31da480fda58'

    tester.set_format("gb")
    assert align_to_hash(Alb.concat_alignments(Alb.make_copy(tester),
                                               "(.).(.)-Panx(.)")) == '5ac908ebf7918a45664a31da480fda58'

    tester.set_format("gb")
    assert align_to_hash(Alb.concat_alignments(Alb.make_copy(tester),
                                               "...", "Panx.*")) == 'e754350b0397cf54f531421d1e85774f'

    tester.set_format("gb")
    assert align_to_hash(Alb.concat_alignments(Alb.make_copy(tester),
                                               "...", "(P)an(x)(.)")) == '5c6653aec09489cadcbed68fbd2f7465'

    shorten = Alb.delete_records(Alb.make_copy(tester), "Ccr")
    tester.alignments[1] = shorten.alignments[1]
    assert align_to_hash(Alb.concat_alignments(Alb.make_copy(tester))) == 'f3ed9139ab6f97042a244d3f791228b6'

# ###########################################  '-con', '--consensus' ############################################ #
hashes = {'o d g': '888a13e13666afb4d3d851ca9150b442', 'o d n': '560d4fc4be7af5d09eb57a9c78dcbccf',
          'o d py': '01f1181187ffdba4fb08f4011a962642', 'o d s': '51b5cf4bb7d591c9d04c7f6b6bd70692',
          'o r n': '1123b95374085b5bcd079880b7762801', 'o p g': '2c038a306713800301b6b4cdbcf61659',
          'o p n': '756a3334c70f9272e2d9cb74dba9ad52', 'o p py': 'aaf1d5aff561c1769dd267ada2fea8b0',
          'o p s': 'b6f72510eeef6be0752ae86d72a44283', 'm d py': '0ae422fa0fafbe0f2edab9a042fb7834',
          'm d s': '7b0aa3cca159b276158cf98209be7dab', 'm p py': '460033d892db36d4750bafc6998d42d0',
          'm p s': '89130797253646e61b78ab7d91ad3fd9'}

hashes = [(alignbuddy, hashes[key]) for key, alignbuddy in alb_resources.get("o m d r p g n py s").items()]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_consensus(alignbuddy, next_hash):
    tester = Alb.consensus_sequence(alignbuddy)
    assert align_to_hash(tester) == next_hash

# ###########################################  '-dr', '--delete_records' ############################################ #
hashes = {'o d g': 'b418ba198da2b4a268a962db32cc2a31', 'o d n': '355a98dad5cf382797eb907e83940978',
          'o d py': 'fe9a2776558f3fe9a1732c777c4bc9ac', 'o d s': '35dc92c4f4697fb508eb1feca43d9d75',
          'o r n': '96e6964115200d46c7cb4eb975718304', 'o p g': '50e09d37a92af595f6fe881d4e57bfc5',
          'o p n': '1cfaa4109c5db8fbfeaedabdc57af655', 'o p py': '1d0e7b4d8e89b42b0ef7cc8c40ed1a93',
          'o p s': '1578d98739d2aa6196463957c7b408fa', 'm d py': 'db4ed247b40707e8e1f0622bb420733b',
          'm d s': 'de5beddbc7f0a7f8e3dc2d5fd43b7b29', 'm p py': '31f91f7dc548e4b075bfb0fdd7d5c82c',
          'm p s': '043e35023b355ed560166db9130cfe30'}

hashes = [(alignbuddy, hashes[key]) for key, alignbuddy in alb_resources.get("o m d r p g n py s").items()]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_delete_records(alignbuddy, next_hash):
    tester = Alb.delete_records(alignbuddy, "α[1-5]|β[A-M]")
    assert align_to_hash(tester) == next_hash

# ######################  'd2r', '--transcribe' and 'r2d', '--reverse_transcribe' ###################### #
d2r_hashes = {'o d g': '4bf291d91d4b27923ef07c660b011c72', 'o d n': 'e531dc31f24192f90aa1f4b6195185b0',
              'o d py': 'e55bd18b6d82a7fc3150338173e57e6a', 'o d s': '45b511f34653e3b984e412182edee3ca',
              'm d py': '16cb082f5cd9f103292ccea0c4d65a06', 'm d s': 'd81dae9714a553bddbf38084f7a8e00e'}

r2d_hashes = {'o d g': '2a42c56df314609d042bdbfa742871a3', 'o d n': 'cb1169c2dd357771a97a02ae2160935d',
              'o d py': '503e23720beea201f8fadf5dabda75e4', 'o d s': '228e36a30e8433e4ee2cd78c3290fa6b',
              'm d py': '42679a32ebd93b628303865f68b0293d', 'm d s': 'ae352b908be94738d6d9cd54770e5b5d'}

hashes = [(alignbuddy, d2r_hashes[key], r2d_hashes[key]) for
          key, alignbuddy in alb_resources.get("o m d g py s").items()]


@pytest.mark.parametrize("alignbuddy,d2r_hash,r2d_hash", hashes)
def test_transcribe(alignbuddy, d2r_hash, r2d_hash):
    tester = Alb.dna2rna(alignbuddy)
    assert align_to_hash(tester) == d2r_hash
    tester = Alb.rna2dna(tester)
    assert align_to_hash(tester) == r2d_hash


def test_transcribe_exceptions():
    with pytest.raises(TypeError) as e:
        Alb.dna2rna(alb_resources.get_one("o p s"))
    assert "TypeError: DNA sequence required, not IUPACProtein()." in str(e)

    with pytest.raises(TypeError) as e:
        Alb.dna2rna(alb_resources.get_one("o r n"))
    assert "TypeError: DNA sequence required, not IUPACAmbiguousRNA()." in str(e)


def test_reverse_transcribe_exceptions():  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError) as e:
        Alb.rna2dna(alb_resources.get_one("o p s"))
    assert "TypeError: RNA sequence required, not IUPACProtein()." in str(e)

    with pytest.raises(TypeError) as e:
        Alb.rna2dna(alb_resources.get_one("o d s"))
    assert "TypeError: RNA sequence required, not IUPACAmbiguousDNA()." in str(e)

# ###########################################  '-et', '--enforce_triplets' ########################################### #
hashes = {'o d g': '6ff2a8a7c58bb6ac0d98fe373981e220', 'o d n': 'c907d29434fe2b45db60f1a9b70f110d',
          'o d py': 'b6cf61c86588023b58257c9008c862b5', 'o r n': '0ed7383ab2897f8350c2791739f0b0a4',
          "m d py": "669ffc4fa602fb101c559cb576bddee1"}
hashes = [(alignbuddy, hashes[key]) for key, alignbuddy in alb_resources.get("m o d r g n py").items()]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_enforce_triplets(alignbuddy, next_hash):
    tester = Alb.enforce_triplets(alignbuddy)
    assert align_to_hash(tester) == next_hash


def test_enforce_triplets_error():
    with pytest.raises(TypeError) as e:
        Alb.enforce_triplets(alb_resources.get_one("m p c"))
    assert "Nucleic acid sequence required, not protein." in str(e)

    with pytest.raises(TypeError) as e:
        tester = Alb.enforce_triplets(alb_resources.get_one("m d pr"))
        tester.alignments[0][0].seq = Seq("MLDILSKFKGVTPFKGITIDDGWDQLNRSFMFVLLVVMGTTVTVRQYTGSVISCDGFKKFGSTFAEDYCWTQGLY",
                                          alphabet=IUPAC.protein)
        Alb.enforce_triplets(tester)
    assert "Record 'Mle-Panxα9' is protein. Nucleic acid sequence required." in str(e)

# ###########################################  'er', '--extract_range' ############################################ #
hashes = {'o d g': '4cef071777cfa87c45302f01b661b2c9', 'o d n': '10ca718b74f3b137c083a766cb737f31',
          'o d py': 'd738a9ab3ab200a7e013177e1042e86c', 'o p g': '61646ba8e9e0adae2f57f078af9d0ad3',
          'o p n': '5f400edc6f0990c0cd6eb52ae7687e39', 'o p py': '69c9ad73ae02525150d4682f9dd68093',
          "m d py": "d06ba679c8a686c8f077bb460a4193b0", "m p py": "8151eeda36b9a170512709829d70230b"}
hashes = [(alignbuddy, hashes[key]) for key, alignbuddy in alb_resources.get("m o d p g n py").items()]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
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
    assert align_to_hash(tester) == '249c88cb64d41c47388514c65bf8fff1'


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
    assert tester._out_format == 'fasta'


@pytest.mark.generate_alignments
def test_prank_outputs():
    # NEXUS
    tester = Sb.pull_recs(Sb.SeqBuddy(resource("Mnemiopsis_cds.fa")), 'α1')
    tester = Alb.generate_msa(tester, 'prank', '-f=nexus -once')
    assert tester._out_format == 'nexus'
    # PHYLIPI
    tester = Sb.pull_recs(Sb.SeqBuddy(resource("Mnemiopsis_cds.fa")), 'α1')
    tester = Alb.generate_msa(tester, 'prank', '-f=phylipi -once')
    assert tester._out_format == 'phylip-relaxed'
    # PHYLIPS
    tester = Sb.pull_recs(Sb.SeqBuddy(resource("Mnemiopsis_cds.fa")), 'α1')
    tester = Alb.generate_msa(tester, 'prank', '-f=phylips -once')
    assert tester._out_format == 'phylipsr'


@pytest.mark.generate_alignments
def test_muscle_inputs():
    # FASTA
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'muscle')
    assert align_to_hash(tester) == '91542667cef761ccaf39d8cb4e877944'


@pytest.mark.generate_alignments
def test_muscle_outputs():
    # FASTA
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'muscle', '-clw')
    assert align_to_hash(tester) == '91542667cef761ccaf39d8cb4e877944'


@pytest.mark.generate_alignments
def test_muscle_multi_param():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'muscle', '-clw -diags')
    assert align_to_hash(tester) == '91542667cef761ccaf39d8cb4e877944'


@pytest.mark.generate_alignments
def test_clustalw2_inputs():
    # FASTA
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalw2')
    assert align_to_hash(tester) == '955440b5139c8e6d7d3843b7acab8446'


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
    assert align_to_hash(tester) == '337380bfd03af4ce33f1246beba1082f'
    # PHYLIP
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.phy"))
    tester = Alb.generate_msa(tester, 'clustalomega')
    assert align_to_hash(tester) == '8299780bf9485b89a2f3462ead666142'
    # STOCKHOLM
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.stklm"))
    tester = Alb.generate_msa(tester, 'clustalomega')
    assert align_to_hash(tester) == '0f0a438dba121b3dbbd844bc213c11c7'


@pytest.mark.generate_alignments
def test_clustalomega_outputs():
    # CLUSTAL
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=clustal')
    assert align_to_hash(tester) == '970f6e4389f77a30563763937a3d32bc'
    # PHYLIP
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=phylip')
    assert align_to_hash(tester) == '692c6af848bd90966f15908903894dbd'
    # STOCKHOLM
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=stockholm')
    assert align_to_hash(tester) == 'ea79bce347b5590dc0f1d38b7710a17e'


@pytest.mark.generate_alignments
def test_clustalomega_multi_param():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=clustal --iter=1')
    assert align_to_hash(tester) == '25480f7a9340ff643bb7eeb326e8f981'


@pytest.mark.generate_alignments
def test_mafft_inputs():
    # FASTA
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'mafft')
    assert align_to_hash(tester) == 'd6046c77e2bdb5683188e5de653affe5'


@pytest.mark.generate_alignments
def test_mafft_outputs():
    # CLUSTAL
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'mafft', '--clustalout')
    assert align_to_hash(tester) == 'd6046c77e2bdb5683188e5de653affe5'


@pytest.mark.generate_alignments
def test_mafft_multi_param():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Alb.generate_msa(tester, 'mafft', '--clustalout --noscore')
    assert align_to_hash(tester) == 'd6046c77e2bdb5683188e5de653affe5'


@pytest.mark.generate_alignments
def test_generate_alignment_keep_temp(monkeypatch):
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    TEMP_DIR.subdir("ga_temp_files")

    def ask_false(*ask_args):
        if ask_args:
            pass
        return False

    def ask_true(*ask_args):
        if ask_args:
            pass
        return True

    monkeypatch.setattr("MyFuncs.ask", ask_false)
    with pytest.raises(SystemExit):
        Alb.generate_msa(tester, "pagan", keep_temp="%s/ga_temp_files" % TEMP_DIR.path)

    monkeypatch.setattr("MyFuncs.ask", ask_true)
    Alb.generate_msa(tester, "pagan", keep_temp="%s/ga_temp_files" % TEMP_DIR.path)
    assert os.path.isfile("%s/ga_temp_files/result.fas" % TEMP_DIR.path)
    assert os.path.isfile("%s/ga_temp_files/result.tre" % TEMP_DIR.path)
    assert os.path.isfile("%s/ga_temp_files/tmp.fa" % TEMP_DIR.path)


@pytest.mark.generate_alignments
def test_generate_alignments_genbank():
    tester = Sb.SeqBuddy(resource("Mnemiopsis_pep.gb"))
    tester = Alb.generate_msa(tester, "mafft")
    assert align_to_hash(tester) == "759d0e208f799fe46d35c15804ea3b5a"


@pytest.mark.generate_alignments
def test_generate_alignments_edges1(capsys):
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))

    with pytest.raises(AttributeError) as e:
        Alb.generate_msa(tester, "foo")
    assert "foo is not a supported alignment tool." in str(e)

    # noinspection PyUnresolvedReferences
    with mock.patch.dict('os.environ'):
        del os.environ['PATH']
        with pytest.raises(SystemExit):
            Alb.generate_msa(tester, "mafft")
        out, err = capsys.readouterr()
        assert "#### Could not find mafft in $PATH. ####\n" in err


args = [("prank", "-f=phylipi"), ("clustalomega", "--outfmt=foo"), ("clustalw2", "-output=foo"),
        ("prank", "-f=nexus"), ("prank", "-f=foo"), ("pagan", "-f foo"), ("pagan", "-f nexus"), ("pagan", "-f phylipi")]


@pytest.mark.generate_alignments
@pytest.mark.parametrize("tool,params", args)
def test_generate_alignments_edges2(tool, params):
    tester = Sb.SeqBuddy(resource("Mnemiopsis_cds.fa"))
    tester = Sb.pull_recs(tester, "α[2345]")
    Alb.generate_msa(tester, tool, params, quiet=True)


# ######################  '-hi', '--hash_ids' ###################### #
def test_hash_seq_ids():
    tester = alb_resources.get_one("o p g")
    Alb.hash_ids(tester)
    assert len(tester.records()[0].id) == 10

    tester = alb_resources.get_one("m d pr")
    Alb.hash_ids(tester, 25)
    assert len(tester.records()[0].id) == 25
    assert len(tester.hash_map) == 34


def test_hash_seq_ids_errors():
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
uc_hashes = {'o d g': '2a42c56df314609d042bdbfa742871a3', 'o d n': '52e74a09c305d031fc5263d1751e265d',
             'o d py': 'cfe6cb9c80aebd353cf0378a6d284239', 'o d s': 'b82538a4630810c004dc8a4c2d5165ce',
             'o p g': 'bf8485cbd30ff8986c2f50b677da4332', 'o p n': '8b6737fe33058121fd99d2deee2f9a76',
             'o p py': '968ed9fa772e65750f201000d7da670f', 'o p s': 'f35cbc6e929c51481e4ec31e95671638',
             'm d py': '6259e675def07bd4977f4ab1f5ffc26d', 'm d s': 'f3f7b66ef034d3807683a2d5a0c44cad',
             'm p py': '2a77f5761d4f51b88cb86b079e564e3b', 'm p s': '6f3f234d796520c521cb85c66a3e239a'}

lc_hashes = {'o d g': '2a42c56df314609d042bdbfa742871a3', 'o d n': 'cb1169c2dd357771a97a02ae2160935d',
             'o d py': '503e23720beea201f8fadf5dabda75e4', 'o d s': '228e36a30e8433e4ee2cd78c3290fa6b',
             'o p g': 'bf8485cbd30ff8986c2f50b677da4332', 'o p n': '17ff1b919cac899c5f918ce8d71904f6',
             'o p py': 'aacda2f5d4077f23926400f74afa2f46', 'o p s': 'c0dce60745515b31a27de1f919083fe9',
             'm d py': '0974ac9aefb2fb540957f15c4869c242', 'm d s': 'a217b9f6000f9eeff98faeb9fd09efe4',
             'm p py': 'd13551548c9c1e966d0519755a8fb4eb', 'm p s': '00661f7afb419c6bb8c9ac654af7c976'}

hashes = [(alignbuddy, uc_hashes[key], lc_hashes[key],)
          for key, alignbuddy in alb_resources.get("o m d p g py s").items()]


@pytest.mark.parametrize("alignbuddy,uc_hash,lc_hash", hashes)
def test_cases(alignbuddy, uc_hash, lc_hash):
    tester = Alb.uppercase(alignbuddy)
    assert align_to_hash(tester) == uc_hash
    tester = Alb.lowercase(tester)
    assert align_to_hash(tester) == lc_hash


# ##################### '-mf2a', '--map_features2alignment' ###################### ##
hashes = {"o p n": "79078260e8725a0d7ccbed9400c78eae", "o p pr": "02b977e5b086125255b792788014708a",
          "o p psr": "02b977e5b086125255b792788014708a", "o p s": "372daf72435e2f1a06531b5c030995c6",
          "o d n": "9fece109249f4d787c13e6fb2742843d", "o d pr": "899c6ab534f7af07e744eb173e94bd50",
          "o d psr": "899c6ab534f7af07e744eb173e94bd50", "o d s": "79cf44688165842eba1bb45b3543d458"}
hashes = [(alignbuddy, hashes[key]) for key, alignbuddy in alb_resources.get("o p d n pr psr s").items()]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_map_features2alignment(alignbuddy, next_hash):
    if alignbuddy.alpha == IUPAC.protein:
        seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis_pep.gb"))
    else:
        seqbuddy = Sb.SeqBuddy(resource("Mnemiopsis_cds.gb"))
    tester = Alb.map_features2alignment(seqbuddy, alignbuddy)
    tester.set_format("genbank")
    assert align_to_hash(tester) == next_hash


# ###########################################  '-oi', '--order_ids' ############################################ #
fwd_hashes = {'o d g': 'acf68c9196faa0960007abc40ba87244', 'o d n': '60bbc6306cbb4eb903b1212718bb4592',
              'o d py': '3c49bdc1b0fe4e1d6bfc148eb0293e21', 'o p g': '7fb69602ffac95a5eecd106876640fcc',
              'o p n': '9a790b9525ca8b1ac3cae3b98ca24b30', 'o p py': 'ffae954adc0d362354e43c1b70d9be29',
              'm d py': 'a44938e26e4b35967ed8e17a0eaebe4c', 'm p py': '5bdda310b29b18057e056f3c982446b2'}

rev_hashes = {'o d g': 'a593a2cd979f52c356c61e10ca9a1317', 'o d n': '82fea6e3d3615ac75ec5022abce255da',
              'o d py': 'd6e79a5faeaff396aa7eab0b460c3eb9', 'o p g': '39af830e6d3605ea1dd04979a4a33f54',
              'o p n': '85b3562b0eb0246d7dab56a4bcc6e2ae', 'o p py': 'f4c0924087fdb624823d02e909d94e95',
              'm d py': '9d6b6087d07f7d1fd701591ab7cb576d', 'm p py': '439f57b891dd2a72724b10c124f96378'}
hashes = [(alignbuddy, fwd_hashes[key], rev_hashes[key])
          for key, alignbuddy in alb_resources.get("m o d p g n py").items()]


@pytest.mark.parametrize("alignbuddy,fwd_hash,rev_hash", hashes)
def test_order_ids(alignbuddy, fwd_hash, rev_hash):
    Alb.order_ids(alignbuddy)
    assert align_to_hash(alignbuddy) == fwd_hash

    Alb.order_ids(alignbuddy, reverse=True)
    assert align_to_hash(alignbuddy) == rev_hash


# ##################### '-pr', '--pull_records' ###################### ##
hashes = {'o d g': '7d1091e16adc09e658563867e7c6bc35', 'o d n': 'd82e66c57548bcf8cba202b13b070ead',
          'o d py': 'd141752c38a892ccca800c637f609608', 'o p g': '862eb17cff2f6cd20b9cf7835879c90b',
          'o p n': '027bbc7e34522f9521f83ee7d03793a1', 'o p py': '2cd74d7ede4d1fb6e18363567426437e',
          'm d py': '7c77c6f3245c21842f4be585714ec6ce', 'm p py': 'f34fa4c34cfe5c1e6b228949557c9483'}

hashes = [(alignbuddy, hashes[key]) for key, alignbuddy in alb_resources.get("m o d p g n py").items()]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_pull_records(alignbuddy, next_hash):
    Alb.pull_records(alignbuddy, "α[1-5]$|β[A-M]")
    assert align_to_hash(alignbuddy) == next_hash


# ###########################################  '-ri', '--rename_ids' ############################################ #
hashes = {'o d g': 'c35db8b8353ef2fb468b0981bd960a38', 'o d n': '243024bfd2f686e6a6e0ef65aa963494',
          'o d py': '98bb9b57f97555d863054ddb526055b4', 'o p g': '2e4e3c365cd011821bcdc6275a3559af',
          'o p n': '3598e85169ed3bcdbcb676bb2eb6cef0', 'o p py': 'd49eb4de01d727b9e3ad648d6a04a3c9',
          'm d py': 'ddfffd9b999637abf7f5926f017de987', 'm p py': '0a68229bd13439040f045cd8c72d7cc9'}

hashes = [(alignbuddy, hashes[key]) for key, alignbuddy in alb_resources.get("m o d p g n py").items()]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_rename_ids(alignbuddy, next_hash):
    Alb.rename(alignbuddy, 'Panx', 'Test', 0)
    assert align_to_hash(alignbuddy) == next_hash


# ###########################################  'tr', '--translate' ############################################ #
hashes = {'o d f': 'b7fe22a87fb78ce747d80e1d73e39c35', 'o d g': 'a949edce98525924dbbc3ced03c18214',
          'o d n': 'a2586af672ad71f16bbd54f359b323ff', 'o d py': 'd0d4dd408e559215b2780f4f0ae0c418',
          'o d pr': 'f77705e32cd753267916539ee0936e1f', 'o d pss': 'ede672b15221ec60981287ca1e286c52',
          'o d psr': '623fe1634752e812f482cfa7b7ea20ee', 'o d s': '4ff563c39229d30aa3eda193cb290344',
          'o d c': '150179326629fffadb7aef7796bd1cec', 'm d py': '7fd28236f491c38ba261dfde20919595',
          'm d pr': '0de676236eda864172b73b6abe4d7a05', 'm d pss': 'c7feff60c16b2b187e03db5d160a4748',
          'm d psr': 'b5945f0317fe9ce8fc03ac7f4c0d5932', 'm d s': 'ee5a41b6f8b32645359beafc72efe825',
          'm d c': '09251e1f4fc0e07a5bba4c64c22bac9b'}

hashes = [(alignbuddy, hashes[key]) for key, alignbuddy in alb_resources.get("o m d c f g n py pr psr pss s").items()]


@pytest.mark.parametrize("alignbuddy,next_hash", hashes)
def test_translate1(alignbuddy, next_hash):
    Alb.translate_cds(alignbuddy)
    assert align_to_hash(alignbuddy) == next_hash


def test_translate2():
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
hashes = {'o d psr': {3: '5df948e4b2cb6c0d0740984445655135', 0.7: '384563eb411713e90cb2fea0c799bf0d'},
          'm d psr': {3: '0e93f0a8c77da8ec974eeca311ca6636', 0.7: 'b15f333416e9dd44834f468d5cd4ca8d'},
          'o p psr': {3: 'b87f927511aade73bc795e024af8975e', 0.7: 'e0f5ce9201249daf4bb3b4f70a7b5ce8'},
          'm p psr': {3: 'f0f2115e29f6dfcb75036d90b06edab4', 0.7: 'f443fbe1831fe368a11edc51e25fa330'}}

hashes = [(alignbuddy, hashes[key]) for key, alignbuddy in alb_resources.get("o m d p psr").items()]


@pytest.mark.parametrize("alignbuddy,hash_dict", hashes)
def test_trimal(alignbuddy, hash_dict):
    tester1, tester2 = Alb.make_copy(alignbuddy), Alb.make_copy(alignbuddy)
    Alb.trimal(tester1, 3)
    assert align_to_hash(tester1) == hash_dict[3]

    tester1, tester2 = Alb.make_copy(alignbuddy), Alb.make_copy(alignbuddy)
    Alb.trimal(tester1, 0.7)
    assert align_to_hash(tester1) == hash_dict[0.7]


def test_trimal2():
    assert align_to_hash(Alb.trimal(alb_resources.get_one("o p n"), 'all')) == "8faaf09741ddb3137653cb77ee66974a"
    tester = alb_resources.get_one("o p n")
    tester.alignments[0]._records = tester.alignments[0]._records[:5]
    Alb.trimal(tester, 'clean')
    assert align_to_hash(tester) == "93a2aa21e6baf5ca70eb2de52ae8dbea"
    tester = alb_resources.get_one("o p n")
    tester_dir = TEMP_DIR.subdir()
    tester.write("%s/trimal" % tester_dir)
    assert align_to_hash(Alb.trimal(tester, 'gappyout')) == "2877ecfb201fc35211a4625f34c7afdd"
    real_trimal = Popen("trimal -in %s/trimal -gappyout" % tester_dir, stdout=PIPE, shell=True).communicate()
    real_trimal = real_trimal[0].decode()
    with open("%s/trimal" % tester_dir, "w") as ofile:
        ofile.write(real_trimal)
    tester = Alb.AlignBuddy("%s/trimal" % tester_dir)
    assert align_to_hash(tester) == "2877ecfb201fc35211a4625f34c7afdd"

    records = [SeqRecord(Seq("A--G-")), SeqRecord(Seq("--T--")), SeqRecord(Seq("--TG-")), SeqRecord(Seq("A---C"))]
    tester = Alb.AlignBuddy([MultipleSeqAlignment(records)])
    Alb.trimal(tester, "gappyout")
    assert "".join([str(rec.seq) for rec in tester.records()]) == ""


# ################################################# COMMAND LINE UI ################################################## #
# ###################### argparse_init() ###################### #
def test_argparse_init(capsys):
    sys.argv = ['AlignBuddy.py', resource("Mnemiopsis_pep.phy"), "-con", "-o", "stockholm"]
    temp_in_args, alignbuddy = Alb.argparse_init()
    assert align_to_hash(alignbuddy) == "5d9a03d9e1b4bf72d991257d3a696306"

    sys.argv = ['AlignBuddy.py', resource("Mnemiopsis_pep.phy"), "-con", "-o", "foo"]
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "Format type 'foo' is not recognized/supported" in err

    sys.argv = ['AlignBuddy.py', resource("gibberish.fa"), "-con"]
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "GuessError: Could not determine format from _input file" in err

    sys.argv = ['AlignBuddy.py', resource("gibberish.fa"), "-con", "-f", "phylip"]
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "ValueError: First line should have two integers" in err

    sys.argv = ['AlignBuddy.py', resource("malformed_phylip_records.physs"), "-con", "-f", "phylipss"]
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "PhylipError: Malformed Phylip --> 9 sequences expected, 8 found." in err

    sys.argv = ['AlignBuddy.py', resource("Mnemiopsis_pep.phy"), "-con", "-f", "foo"]
    with pytest.raises(SystemExit):
        Alb.argparse_init()
    out, err = capsys.readouterr()
    assert "TypeError: Format type 'foo' is not recognized/supported" in err

    sys.argv = ['AlignBuddy.py', resource("Mnemiopsis_pep.fa"), "--quiet", "--generate_alignment", "mafft", "--reorder"]
    temp_in_args, alignbuddy = Alb.argparse_init()
    assert alignbuddy == []


# ##################### '-al', '--alignment_lengths' ###################### ##
def test_alignment_lengths_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.alignment_lengths = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p c"), skip_exit=True)
    out, err = capsys.readouterr()
    assert out == "481\n683\n"
    assert err == "# Alignment 1\n# Alignment 2\n"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("o p py"), skip_exit=True)
    out, err = capsys.readouterr()
    assert out == "681\n"
    assert err == ""


# ##################### '-cs', '--clean_seqs' ###################### ##
def test_clean_seqs_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.clean_seq = [[None]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p pr"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "73b5d11dd25dd100648870228ab10d3d"

    test_in_args.clean_seq = [['strict', 'X']]
    Alb.command_line_ui(test_in_args, Alb.AlignBuddy(resource("ambiguous_dna_alignment.fa")), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "6755ea1408eddd0e5f267349c287d989"


# ##################### '-cta', '--concat_alignments' ###################### ##
def test_concat_alignments_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.concat_alignments = [[]]

    tester = Sb.SeqBuddy(resource("Cnidaria_pep.nexus"))
    Sb.pull_recs(tester, "Ccr|Cla|Hec")
    tester = Alb.AlignBuddy(str(tester))
    tester.alignments.append(tester.alignments[0])
    tester.set_format("genbank")
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "d21940f3dad2295dd647f632825d8541"

    test_in_args.concat_alignments = [["(.).(.)-Panx(.)"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "5ac908ebf7918a45664a31da480fda58"

    test_in_args.concat_alignments = [["...", "Panx.*"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "e754350b0397cf54f531421d1e85774f"

    test_in_args.concat_alignments = [[3, "Panx.*"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "e754350b0397cf54f531421d1e85774f"

    test_in_args.concat_alignments = [[-9, "Panx.*"]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "9d2886afc640d35618754e05223032a2"

    test_in_args.concat_alignments = [[3, 3]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "4e4101f9b5a6d44d524a9783a8c4004b"

    test_in_args.concat_alignments = [[3, -3]]
    Alb.command_line_ui(test_in_args, Alb.make_copy(tester), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "5d9d9ac8fae604be74c436e5f0b5b6db"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("p o g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Please provide at least two alignments." in err

    test_in_args.concat_alignments = [["foo"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p c"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "No match found for record" in err


# ##################### '-con', '--consensus' ###################### ##
def test_consensus_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.consensus = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "7b0aa3cca159b276158cf98209be7dab"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "89130797253646e61b78ab7d91ad3fd9"


# ##################### '-dr', '--delete_records' ###################### ##
def test_delete_records_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.delete_records = [["α[1-5]", "β[A-M]"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "de5beddbc7f0a7f8e3dc2d5fd43b7b29"
    assert string2hash(err) == "31bb4310333851964015e21562f602c2"

    test_in_args.delete_records = [["α[1-5]", "β[A-M]", 4]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(err) == "ce6c9b29c95ba853eb444de5c71aeca9"

    test_in_args.delete_records = [["foo"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "No sequence identifiers match 'foo'\n" in err


# ##################### '-et', '--enforce_triplets' ###################### ##
def test_enforce_triplets_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.enforce_triplets = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "6ff2a8a7c58bb6ac0d98fe373981e220"

    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p c"), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Nucleic acid sequence required, not protein." in err


# ##################### '-er', '--extract_range' ###################### ##
def test_extract_range_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.extract_range = [10, 110]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "3929c5875a58e9a1e64425d4989e590a"

    test_in_args.extract_range = [110, 10]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "3929c5875a58e9a1e64425d4989e590a"

    test_in_args.extract_range = [-110, 10]
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"))
    out, err = capsys.readouterr()
    assert err == "ValueError: Please specify positive integer indices\n"


# ##################### '-ga', '--generate_alignment' ###################### ##
@pytest.mark.generate_alignments
def test_generate_alignment_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.generate_alignment = [[]]

    test_in_args.alignments = [resource("Mnemiopsis_cds.gb")]
    test_in_args.out_format = "gb"
    Alb.command_line_ui(test_in_args, Alb.AlignBuddy, skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "9979f319dfec89a52267bae86fd3e4ad"

    # noinspection PyUnresolvedReferences
    with mock.patch.dict('os.environ'):
        del os.environ['PATH']
        with pytest.raises(SystemExit):
            Alb.command_line_ui(test_in_args, Alb.AlignBuddy)
        out, err = capsys.readouterr()
        assert "Unable to identify any supported alignment tools on your system." in err

    test_in_args.generate_alignment = [["foo"]]
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, Alb.AlignBuddy)
    out, err = capsys.readouterr()
    assert "foo is not a supported alignment tool" in err


# ######################  '-hsi', '--hash_ids' ###################### #
def test_hash_seq_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.hash_ids = [None]
    tester = alb_resources.get_one("m p s")
    ids = [rec.id for rec in tester.records()]
    Alb.command_line_ui(test_in_args, tester, True)
    for indx, rec in enumerate(tester.records()):
        assert rec.id != ids[indx]
        assert ids[indx] == tester.hash_map[rec.id]

    test_in_args.hash_ids = [0]
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert "Warning: The hash_length parameter was passed in with the value 0. This is not a positive integer" in err

    tester.alignments *= 10
    test_in_args.hash_ids = [1]
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert "cover all sequences, so it has been increased to 2" in err


# ###############################  '-li', '--list_ids' ############################## #
def test_list_ids(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.list_ids = [False]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "f087f9c1413ba66c28fb0fccf7c974e6"

    test_in_args.list_ids = [3]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "4d85249a1f187d38d411a78ced65a98c"


# #################################### '-lc', '--lowercase' ################################### #
def test_lowercase_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.lowercase = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "00661f7afb419c6bb8c9ac654af7c976"


# ##################### '-mf2a', '--map_features2alignment' ###################### ##
def test_map_features2alignment_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.mapfeat2align = [resource("Mnemiopsis_cds.gb")]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d n"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "9fece109249f4d787c13e6fb2742843d"


# ##############################################  '-ns', '--num_seqs'  ############################################### #
def test_num_seqs_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.num_seqs = True
    for alignbuddy in alb_resources.get_list("m d c pr psr s"):
        Alb.command_line_ui(test_in_args, alignbuddy, skip_exit=True)
        out, err = capsys.readouterr()
        assert out == "# Alignment 1\n8\n\n# Alignment 2\n21\n" or out == "# Alignment 1\n13\n\n# Alignment 2\n21\n"

    for alignbuddy in alb_resources.get_list("m p c pr psr s"):
        Alb.command_line_ui(test_in_args, alignbuddy, skip_exit=True)
        out, err = capsys.readouterr()
        assert out == "# Alignment 1\n20\n\n# Alignment 2\n13\n" or out == "# Alignment 1\n13\n\n# Alignment 2\n21\n"

    for alignbuddy in alb_resources.get_list("o p c pr psr s"):
        Alb.command_line_ui(test_in_args, alignbuddy, skip_exit=True)
        out, err = capsys.readouterr()
        assert out == "13\n"


# ##############################################  '-oi', '--order_ids' ############################################### #
def test_order_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.order_ids = [False]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "0bce6aeaab76feda1aea7f5e79608c72"

    test_in_args.order_ids = ['rev']
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "f87f78dc9d7a0b76854f15c52130e3a7"


# ##################### '-pr', '--pull_records' ###################### ##
def test_pull_records_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.pull_records = [["α[1-5]$", "β[A-M]"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "2de557d6fd3dc6cd1bf43a1995392a4c"
    assert err == ""

    test_in_args.pull_records = [["α[1-5]$", "ML218922a", "full"]]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "fb82ffec15ece60a74d9ac8db92d2999"


# ######################  '-ri', '--rename_ids' ###################### #
def test_rename_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.rename_ids = [["[a-z](.)", "?\\1", 2]]
    tester = alb_resources.get_one("m p s")
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "888f2e3feb9e67f9bc008183082c822a"

    test_in_args.rename_ids = [["[a-z](.)"]]
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "rename_ids requires two or three argments:" in err

    test_in_args.rename_ids = [["[a-z](.)", "?\\1", 2, "foo"]]
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "rename_ids requires two or three argments:" in err

    test_in_args.rename_ids = [["[a-z](.)", "?\\1", "foo"]]
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "Max replacements argument must be an integer" in err

    test_in_args.rename_ids = [["[a-z](.)", "?\\1\\2", 2]]
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "There are more replacement" in err


# ##################### '-r2d', '--reverse_transcribe' ###################### ##
def test_reverse_transcribe_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.reverse_transcribe = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o r n"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "f8c2b216fa65fef9c74c1d0c4abc2ada"

    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, alb_resources.get_one("m d s"))
    out, err = capsys.readouterr()
    assert err == "TypeError: RNA sequence required, not IUPACAmbiguousDNA().\n"


# ######################  '-sf', '--screw_formats' ###################### #
hashes = [("fasta", "cfa898d43918055b6a02041195874da9"), ("gb", "ceac7a2a57aa8e3f530f70e2765f9ab2"),
          ("nexus", "49bf9b3f56104e4f19048523d725f025"), ("phylip", "968ed9fa772e65750f201000d7da670f"),
          ("phylipr", "5064c1d6ae6192a829972b7ec0f129ed"), ("phylipss", "4bd927145de635c429b2917e0a1db176"),
          ("phylipsr", "b46b57ede57f12c3c3b906681882f81a"), ("stockholm", "5d9a03d9e1b4bf72d991257d3a696306"),
          ("clustal", "9d328711cf6f6750c33373a912efb521")]


@pytest.mark.parametrize("_format,next_hash", hashes)
def test_screw_formats_ui(_format, next_hash, capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = _format
    tester = alb_resources.get_one("o p py")
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == next_hash

hashes = [("clustal", "cf349d6061c602439b72b51368f694ed"), ("phylip", "2a77f5761d4f51b88cb86b079e564e3b"),
          ("phylipr", "1f172a3beef76e8e3d42698bb2c3c87d"), ("phylipss", "eb82cda31fcb2cf00e11d7e910fde695"),
          ("phylipsr", "368169cb86c6ddb7074ed89e2d42c4dd"), ("stockholm", "f221b9973aef4771169136a30bd030fa")]


@pytest.mark.parametrize("_format,next_hash", hashes)
def test_screw_formats_ui2(_format, next_hash, capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = _format
    tester = alb_resources.get_one("m p py")
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert string2hash(out) == next_hash


def test_screw_formats_ui3(capsys):
    test_in_args = deepcopy(in_args)
    tester = alb_resources.get_one("o p py")
    tester.write("%s/seq.phy" % TEMP_DIR.path)
    test_in_args.in_place = True
    test_in_args.alignments = ["%s/seq.phy" % TEMP_DIR.path]
    test_in_args.screw_formats = "genbank"
    Alb.command_line_ui(test_in_args, tester, True)
    assert os.path.isfile("%s/seq.gb" % TEMP_DIR.path)

    test_in_args.in_place = False
    test_in_args.screw_formats = "foo"
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "TypeError: Format type 'foo' is not recognized/supported\n" in err

    test_in_args.screw_formats = "gb"
    tester = alb_resources.get_one("m p py")
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "ValueError: gb format does not support multiple alignments in one file." in err


# ##################################  '-stf', '--split_to_files' ################################### #
hashes = [("clustal", "m p c"), ("phylip", "m p py"),
          ("phylip-relaxed", "m p pr"), ("phylipss", "m p pss"),
          ("phylipsr", "m p psr"), ("stockholm", "m p s")]


@pytest.mark.parametrize("_format,_code", hashes)
def test_split_alignment_ui(_format, _code):
    TEMP_DIR.subdir(_format)
    test_in_args = deepcopy(in_args)
    test_in_args.split_to_files = [["%s/%s" % (TEMP_DIR.path, _format), "Foo_"]]
    tester = alb_resources.get_one(_code)
    for num in range(7):
        tester.alignments += tester.alignments
    Alb.command_line_ui(test_in_args, tester, True)
    assert os.path.isfile("%s/%s/Foo_001.%s" % (TEMP_DIR.path, _format, br.format_to_extension[_format]))
    assert os.path.isfile("%s/%s/Foo_061.%s" % (TEMP_DIR.path, _format, br.format_to_extension[_format]))
    assert os.path.isfile("%s/%s/Foo_121.%s" % (TEMP_DIR.path, _format, br.format_to_extension[_format]))
    TEMP_DIR.del_subdir(_format)


def test_split_alignment_ui2(capsys):
    TEMP_DIR.subdir("split_alignment")
    test_in_args = deepcopy(in_args)
    test_in_args.split_to_files = [["%s/split_alignment" % TEMP_DIR.path, "Foo_", "Bar"]]
    tester = alb_resources.get_one("m p c")
    for num in range(4):
        tester.alignments += tester.alignments
    Alb.command_line_ui(test_in_args, tester, True)
    out, err = capsys.readouterr()
    assert os.path.isfile("%s/split_alignment/Foo_01.clus" % TEMP_DIR.path)
    assert os.path.isfile("%s/split_alignment/Foo_16.clus" % TEMP_DIR.path)

    assert "Warning: Only one prefix can be accepted, 2 where provided. Using the first.\n" in err
    assert "New file: " in err

    tester = alb_resources.get_one("o d c")
    test_in_args.quiet = False
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, tester)
    out, err = capsys.readouterr()
    assert "Only one alignment present, nothing written." in err
    TEMP_DIR.del_subdir("split_alignment")


# ##################### '-d2r', '--transcribe' ###################### ##
def test_transcribe_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.transcribe = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d n"), skip_exit=True)
    out, err = capsys.readouterr()

    assert string2hash(out) == "e531dc31f24192f90aa1f4b6195185b0"

    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, alb_resources.get_one("o r n"))
    out, err = capsys.readouterr()
    assert err == "TypeError: DNA sequence required, not IUPACAmbiguousRNA().\n"


# ##################### '-tr', '--translate' ###################### ##
def test_translate_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.translate = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "a949edce98525924dbbc3ced03c18214"

    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, alb_resources.get_one("o p n"))
    out, err = capsys.readouterr()
    assert "Nucleic acid sequence required, not protein." in err


# ##################### '-trm', '--trimal' ###################### ##
def test_trimal_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.trimal = [False]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "14d48179ced1973adfe06cf85f04cf27"

    test_in_args.trimal = ["gappyout"]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d g"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "14d48179ced1973adfe06cf85f04cf27"

    test_in_args.trimal = [0.25]
    Alb.command_line_ui(test_in_args, alb_resources.get_one("o d psr"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "5df948e4b2cb6c0d0740984445655135"

    test_in_args.trimal = ["foo"]
    with pytest.raises(SystemExit):
        Alb.command_line_ui(test_in_args, alb_resources.get_one("o p n"))
    out, err = capsys.readouterr()
    assert "NotImplementedError: foo not an implemented trimal method\n" == err


# #################################### '-uc', '--uppercase' ################################### #
def test_uppercase_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.uppercase = True
    Alb.command_line_ui(test_in_args, alb_resources.get_one("m p s"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "6f3f234d796520c521cb85c66a3e239a"
