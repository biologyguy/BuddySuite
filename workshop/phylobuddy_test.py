#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: phylobuddy_tests.py
version: 1.1
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Collection of PyTest unit tests for the PhyloBuddy.py program
"""

import pytest
from hashlib import md5
import os
import sys
import argparse
from copy import deepcopy
from collections import OrderedDict
from unittest import mock
import ete3

sys.path.insert(0, "./")
import buddy_resources as br
import PhyloBuddy as Pb
import AlignBuddy as Alb
import MyFuncs

VERSION = Pb.VERSION
WRITE_FILE = MyFuncs.TempFile()


def fmt(prog):
    return br.CustomHelpFormatter(prog)

parser = argparse.ArgumentParser(prog="PhyloBuddy.py", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                 description='''\
\033[1mPhyloBuddy\033[m
Put a little bonsai into your phylogeny.

\033[1mUsage examples\033[m:
PhyloBuddy.py "/path/to/tree_file" -<cmd>
PhyloBuddy.py "/path/to/tree_file" -<cmd> | PhyloBuddy.py -<cmd>
PhyloBuddy.py "(A,(B,C));" -f "raw" -<cmd>
''')

br.flags(parser, ("trees", "Supply file path(s) or raw tree string, If piping trees into PhyloBuddy "
                           "this argument can be left blank."), br.pb_flags, br.pb_modifiers, VERSION)

# This is to allow py.test to work with the -x flag
parser.add_argument("-x", nargs="?")
parser.add_argument("-m", nargs="?")
parser.add_argument("-n", nargs="?")
parser.add_argument("--cov", nargs="?")
parser.add_argument("--cov-report", nargs="?")
in_args = parser.parse_args([])


def phylo_to_hash(_phylobuddy, mode='hash'):
    if mode != "hash":
        return str(_phylobuddy)
    _hash = md5("{0}\n".format(str(_phylobuddy).rstrip()).encode('utf-8')).hexdigest()
    return _hash

root_dir = os.getcwd()


def string2hash(_input):
    return md5(_input.encode('utf-8')).hexdigest()


def resource(file_name):
    return "{0}/unit_test_resources/{1}".format(root_dir, file_name)


class Resources(object):
    def __init__(self):
        one_tree = OrderedDict([("newick", "single_tree.newick"),
                                ("nexus", "single_tree.nex"),
                                ("nexml", "single_tree.xml")])
        multi_tree = OrderedDict([("newick", "multi_tree.newick"),
                                  ("nexus", "multi_tree.nex"),
                                  ("nexml", "multi_tree.xml")])

        self.resources = OrderedDict([("one", one_tree),
                                      ("multi", multi_tree)])
        self.pb_objs = OrderedDict()
        self.res_paths = OrderedDict()
        for num in self.resources:
            self.pb_objs[num] = OrderedDict([(key, Pb.PhyloBuddy(resource(path)))
                                             for key, path in self.resources[num].items()])

            self.res_paths[num] = OrderedDict([(key, resource(path)) for key, path in self.resources[num].items()])

        self.code_dict = OrderedDict([("num_trees", OrderedDict([("o", "one"), ("m", "multi")])),
                                      ("format", OrderedDict([("k", "newick"), ("n", "nexus"), ("l", "nexml")]))])

    def _parse_code(self, code=""):
        results = OrderedDict([("num_trees", []), ("format", [])])
        code = code.split()
        for i in code:
            for j in results:
                if i in self.code_dict[j]:
                    results[j].append(i)

        # Fill up a field with all possibilities if nothing is given
        results["num_trees"] = [key for key in self.code_dict["num_trees"]] \
            if not results["num_trees"] else results["num_trees"]
        results["format"] = [key for key in self.code_dict["format"]] if not results["format"] else results["format"]
        return results

    def get(self, code="", mode="objs"):
        """
        Returns copies of PhyloBuddy objects, the
        :param code:
        :param mode: {"objs", "paths"}
        :return: OrderedDict {key: resource}
        """
        files = self._parse_code(code)
        output = OrderedDict()
        key = ["", ""]
        for num_aligns in files["num_trees"]:
            key[0] = num_aligns
            n = self.code_dict["num_trees"][num_aligns]
            for _format in files["format"]:
                key[1] = _format
                f = self.code_dict["format"][_format]
                try:
                    assert not " ".join(key) in output
                    if mode == "objs":
                        output[" ".join(key)] = Pb.make_copy(self.pb_objs[n][f])
                    elif mode == "paths":
                        output[" ".join(key)] = self.res_paths[n][f]
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
        return {"num_trees": self.code_dict["num_trees"][code[0]],
                "format": br.parse_format(self.code_dict["format"][code[1]])}

pb_resources = Resources()

# Deprecated! Remove once all tests are transitioned over to the Resource syntax
phylo_files = ['multi_tree.newick', 'multi_tree.nex', 'multi_tree.xml', 'single_tree.newick', 'single_tree.nex',
               'single_tree.xml']

file_types = ['newick', 'nexus', 'nexml', 'newick', 'nexus', 'nexml']


@pytest.mark.parametrize("key,pb_path", pb_resources.get("o m k n l", "paths").items())
def test_instantiate_phylobuddy_from_file(key, pb_path):
    _format = pb_resources.deets(key)["format"]
    assert type(Pb.PhyloBuddy(pb_path, _in_format=_format)) == Pb.PhyloBuddy


@pytest.mark.parametrize("phylo_file", phylo_files)
def test_instantiate_phylobuddy_from_file_guess(phylo_file):
    assert type(Pb.PhyloBuddy(resource(phylo_file))) == Pb.PhyloBuddy


@pytest.mark.parametrize("phylo_file", phylo_files)
def test_instantiate_phylobuddy_from_handle(phylo_file):
    with open(resource(phylo_file), 'r') as ifile:
        assert type(Pb.PhyloBuddy(ifile)) == Pb.PhyloBuddy


@pytest.mark.parametrize("phylo_file", phylo_files)
def test_instantiate_phylobuddy_from_raw(phylo_file):
    with open(resource(phylo_file), 'r') as ifile:
        assert type(Pb.PhyloBuddy(ifile.read())) == Pb.PhyloBuddy


@pytest.mark.parametrize("phylo_file", phylo_files)
def test_instantiate_phylobuddy_from_phylobuddy(phylo_file):
    tester = Pb.PhyloBuddy(resource(phylo_file))
    assert type(Pb.PhyloBuddy(tester)) == Pb.PhyloBuddy


@pytest.mark.parametrize("phylo_file", phylo_files)
def test_instantiate_phylobuddy_from_list(phylo_file):
    tester = Pb.PhyloBuddy(resource(phylo_file))
    assert type(Pb.PhyloBuddy(tester.trees)) == Pb.PhyloBuddy


def test_empty_file():
    with open(resource("blank.fa"), "r") as ifile:
        with pytest.raises(SystemExit):
            Pb.PhyloBuddy(ifile)


def test_guess_error():
    # File path
    with pytest.raises(br.GuessError):
        Pb.PhyloBuddy(resource("unrecognizable.txt"))

    with open(resource("unrecognizable.txt"), 'r') as ifile:
        # Raw
        with pytest.raises(br.GuessError):
            Pb.PhyloBuddy(ifile.read())

        # Handle
        with pytest.raises(br.GuessError):
            ifile.seek(0)
            Pb.PhyloBuddy(ifile)

    # GuessError output
    test_error = br.GuessError("This is a test")
    assert str(test_error) == "This is a test"

    try:
        Pb.PhyloBuddy(resource("unrecognizable.txt"))
    except br.GuessError as e:
        assert "Could not automatically determine the format of" in str(e.value) and \
               "\nTry explicitly setting it with the -f flag." in str(e.value)


def test_stderr(capsys):
    Pb._stderr("Hello std_err", quiet=False)
    out, err = capsys.readouterr()
    assert err == "Hello std_err"

    Pb._stderr("Hello std_err", quiet=True)
    out, err = capsys.readouterr()
    assert err == ""


def test_stdout(capsys):
    Pb._stdout("Hello std_out", quiet=False)
    out, err = capsys.readouterr()
    assert out == "Hello std_out"

    Pb._stdout("Hello std_out", quiet=True)
    out, err = capsys.readouterr()
    assert out == ""


def test_phylobuddy_edges():
    # If the input list isn't a list of PhyloBuddy objects
    with pytest.raises(TypeError):
        Pb.PhyloBuddy(["Foo", "Bar"])

    # Catch figtree metadata
    assert type(Pb.PhyloBuddy(resource("figtree.nexus"))) == Pb.PhyloBuddy
    with open(resource("figtree.nexus"), 'r') as ifile:
        tester = Pb.PhyloBuddy(ifile)
        assert type(tester) == Pb.PhyloBuddy

    # Unsupported output format
    tester.out_format = "foo"
    with pytest.raises(TypeError):
        str(tester)

    # No trees in PhyloBuddy object
    tester.trees = []
    assert str(tester) == "Error: No trees in object.\n"

pb_objects = [Pb.PhyloBuddy(resource(x)) for x in phylo_files]

# ################################################# HELPER FUNCTIONS ################################################# #
hashes = ['6843a620b725a3a0e0940d4352f2036f', '543d2fc90ca1f391312d6b8fe896c59c', '6ce146e635c20ad62e21a1ed6fddbd3a',
          '4dfed97b2a23b8957ee5141bf4681fe4', '77d00fdc512fa09bd1146037d25eafa0', '9b1014be1b38d27f6b7ef73d17003dae']

hashes = [(Pb.make_copy(pb_objects[x]), hashes[x]) for x in range(len(pb_objects))]


@pytest.mark.parametrize("phylobuddy,next_hash", hashes)
def test_print(phylobuddy, next_hash, capsys):
    phylobuddy.print()
    out, err = capsys.readouterr()
    out = "{0}\n".format(out.rstrip())
    tester = md5(out.encode()).hexdigest()
    assert tester == next_hash


@pytest.mark.parametrize("phylobuddy,next_hash", hashes)
def test_str(phylobuddy, next_hash):
    tester = str(phylobuddy)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash


@pytest.mark.parametrize("phylobuddy,next_hash", hashes)
def test_write1(phylobuddy, next_hash):
    temp_file = MyFuncs.TempFile()
    phylobuddy.write(temp_file.path)
    out = "{0}\n".format(temp_file.read().rstrip())
    tester = md5(out.encode()).hexdigest()
    assert tester == next_hash


def test_convert_to_ete():
    tester = Pb.make_copy(pb_objects[0])
    tester.trees[0].seed_node.annotations.add_new("pb_color", '#ff0000')
    ete_tree = Pb._convert_to_ete(tester.trees[0])
    assert ete_tree.pb_color == '#ff0000'


def test_guess_format():
    with pytest.raises(br.GuessError):
        Pb._guess_format(dict)


# ################################################ MAIN API FUNCTIONS ################################################ #
# ###################### 'ct', '--consensus_tree' ###################### #
hashes = ["a565a2b86bafbc06cee02197c4a8428a", "64d4b94950652cd754c22aff9ae0dc32",
          "3319ac115720a67f9941968e26b2762a", "59956aecc03e49a4d3cce8505b80b088"]
hashes = [(Pb.make_copy(pb_objects[x]), hashes[x]) for x in range(4)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_consensus_tree(phylobuddy, next_hash):
    tester = Pb.consensus_tree(phylobuddy)
    assert phylo_to_hash(tester) == next_hash


hashes = ["df1cbcf1ba7e01fc11930af1a4dd61ce", "fa8a14163e6e55baa5b648a9b60a44d0",
          "043f8088bde0644dd2c06076802493f5", "59956aecc03e49a4d3cce8505b80b088"]
hashes = [(Pb.make_copy(pb_objects[x]), hashes[x]) for x in range(4)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_consensus_tree_95(phylobuddy, next_hash):
    tester = Pb.consensus_tree(phylobuddy, frequency=0.95)
    assert phylo_to_hash(tester) == next_hash


# ###################### 'dt', '--display_trees' ###################### #
@pytest.mark.display
def test_display_trees(monkeypatch):
    show = mock.Mock(return_value=True)
    monkeypatch.setattr(ete3.TreeNode, "show", show)
    assert Pb.display_trees(pb_resources.get_one("o k"))


def test_display_trees_error():
    # noinspection PyUnresolvedReferences
    with mock.patch.dict('os.environ'):
        if 'DISPLAY' in os.environ:
            del os.environ['DISPLAY']
        with pytest.raises(SystemError):
            Pb.display_trees(Pb.make_copy(pb_objects[0]))


# ###################### 'dis', '--distance' ###################### #
hashes = ['e39a8aadaec77680ad0d9004bab824ea', '3c49c6a7f06244c0b5d45812f6791519', '7df609f2e6ee613d3bf3c3d2aae26ad4']
hashes = [(Pb.make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_distance_wrf(phylobuddy, next_hash):
    tester = str(Pb.distance(phylobuddy, method='wrf'))
    assert md5(tester.encode()).hexdigest() == next_hash

hashes = ['c15d06fc5344da3149e19b134ca31c62', '6d087b86aa9f5bc5013113972173fe0f', '7ef096e3c32dbf898d4b1a035d5c9ad4']
hashes = [(Pb.make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_distance_uwrf(phylobuddy, next_hash):
    tester = str(Pb.distance(phylobuddy, method='uwrf'))
    assert md5(tester.encode()).hexdigest() == next_hash

hashes = ['68942718c8baf4e4bdf5dd2992fbbf9d', '3dba6b10fdd04505b4e4482d926b67d3', '8d0b3a035015d62916b525f371684bf8']
hashes = [(Pb.make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_distance_ed(phylobuddy, next_hash):
    tester = str(Pb.distance(phylobuddy, method='ed'))
    assert md5(tester.encode()).hexdigest() == next_hash


def test_distance_unknown_method():
    with pytest.raises(AttributeError):
        Pb.distance(Pb.make_copy(pb_objects[0]), method='foo')


# ######################  'gt', '--generate_trees' ###################### #
# Hashes are for RAxML version 8.2.3
@pytest.mark.generate_trees
def test_raxml_inputs():
    # Nucleotide
    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))
    assert phylo_to_hash(Pb.generate_tree(tester, 'raxml')) == '706ba436f8657ef3aee7875217dd07c0'

    # Peptide
    tester = Alb.AlignBuddy(resource("Mnemiopsis_pep.nex"))
    assert phylo_to_hash(Pb.generate_tree(Alb.make_copy(tester), 'raxml')) == 'fc35569091eeba49ac4dcec7fc6890bf'

    # Quiet
    assert phylo_to_hash(Pb.generate_tree(tester, 'raxml', quiet=True)) == 'fc35569091eeba49ac4dcec7fc6890bf'


@pytest.mark.generate_trees
def test_raxml_multi_param():
    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))
    tester = Pb.generate_tree(tester, 'raxml', '-m GTRCAT -p 112358 -K MK -N 2')
    assert phylo_to_hash(tester) == '2bce58a9c6756fa68fd828a307850d7d'


# PhyML version 20120412
@pytest.mark.generate_trees
def test_phyml_inputs():
    # Nucleotide
    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))
    tester = Pb.generate_tree(tester, 'phyml', '-m GTR --r_seed 12345')
    assert phylo_to_hash(tester) == 'd3a4e7601998885f333ddd714ca764db'
    # Peptide
    tester = Alb.AlignBuddy(resource("Mnemiopsis_pep.nex"))
    tester = Pb.generate_tree(tester, 'phyml', '-m Blosum62 --r_seed 12345')
    assert phylo_to_hash(tester) == '52c7d028341b250bcc867d57a68c794c'


@pytest.mark.generate_trees
def test_phyml_multi_param():
    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))
    tester = Pb.generate_tree(tester, 'phyml', '-m GTR -o tl -b 2 --r_seed 12345')
    assert phylo_to_hash(tester) == '5434f29509eab76dd52dd69d2c0e186f'


@pytest.mark.generate_trees
def test_fasttree_inputs():
    temp_dir = MyFuncs.TempDir()
    # Nucleotide
    alignbuddy = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))

    tester = Pb.generate_tree(Alb.make_copy(alignbuddy), 'FastTree', '-seed 12345')
    assert phylo_to_hash(tester) == 'd7f505182dd1a1744b45cc326096f70c'

    tester = Pb.generate_tree(alignbuddy, 'fasttree', '-seed 12345', quiet=True)
    assert phylo_to_hash(tester) == 'd7f505182dd1a1744b45cc326096f70c'

    alignbuddy = Alb.AlignBuddy(resource("Mnemiopsis_pep.nex"))
    tester = Pb.generate_tree(alignbuddy, 'fasttree', '-seed 12345', keep_temp="%s/new_dir" % temp_dir.path)
    assert phylo_to_hash(tester) == '57eace9bdd2074297cbd2692c1f4cd38'


@pytest.mark.generate_trees
def test_fasttree_multi_param():
    temp_file = MyFuncs.TempFile()
    tester = Alb.AlignBuddy(resource("Alignments_cds.phyr"))
    tester = Pb.generate_tree(tester, 'FastTree', '-seed 12345 -wag -fastest -log %s' % temp_file.path)
    assert phylo_to_hash(tester) == '0877f4e8f46c3f77390dbf962d24ff71'


def test_generate_trees_edge_cases():
    temp_file = MyFuncs.TempFile()
    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))
    with pytest.raises(FileExistsError):
        Pb.generate_tree(tester, "raxml", keep_temp=temp_file.path)

    with pytest.raises(AttributeError):
        Pb.generate_tree(tester, "foo")

    with pytest.raises(FileNotFoundError):
        Pb.generate_tree(tester, "raxml", "-m BINCATLG")

    with pytest.raises(RuntimeError):
        Pb.generate_tree(tester, "fasttree", "-s 12345")

    # noinspection PyUnresolvedReferences
    with mock.patch.dict(os.environ, {"PATH": ""}):
        with pytest.raises(ProcessLookupError):
            Pb.generate_tree(tester, "raxml")


# ###################### 'hi', '--hash_ids' ###################### #
@pytest.mark.parametrize("phylobuddy", pb_resources.get_list("m o k n l"))
def test_hash_ids(phylobuddy):
    orig_hash = phylo_to_hash(phylobuddy)
    Pb.hash_ids(phylobuddy)
    assert phylo_to_hash(phylobuddy) != orig_hash


def test_hash_ids_edges():
    with pytest.raises(TypeError) as e:
        Pb.hash_ids(Pb.PhyloBuddy, hash_length="foo")
    assert "Hash length argument must be an integer, not <class 'str'>" in str(e)

    with pytest.raises(ValueError) as e:
        Pb.hash_ids(Pb.PhyloBuddy, hash_length=0)
    assert "Hash length must be greater than 0" in str(e)

    with pytest.raises(ValueError) as e:
        Pb.hash_ids(pb_resources.get_one("m n"), hash_length=1)
    assert "Insufficient number of hashes available to cover all sequences." in str(e)

    tester = Pb.PhyloBuddy(resource("tree_with_node_lables.nwk"))
    test_hash = phylo_to_hash(tester)
    tester = Pb.hash_ids(tester, hash_length=5, nodes=True)
    assert phylo_to_hash(tester) != test_hash

# ###################### 'li', '--list_ids' ###################### #
li_hashes = ['514675543e958d5177f248708405224d', '229e5d7cd8bb2bfc300fd45ec18e8424', '514675543e958d5177f248708405224d']
li_hashes = [(Pb.make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(li_hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", li_hashes)
def test_list_ids(phylobuddy, next_hash):
    tester = str(Pb.list_ids(phylobuddy))
    assert md5(tester.encode()).hexdigest() == next_hash

# ###################### 'ptr', '--print_trees' ###################### #
hashes = ['1f62f35b8fe64c5eaaf394c272febf4f', 'd5b502ae7e83f0c6ae2211b7d8281b90', '1f62f35b8fe64c5eaaf394c272febf4f',
          '8a5a1f3e409471184da4138ab45fc54e', '370b1c17c1d96ebb866277747f276add', '8a5a1f3e409471184da4138ab45fc54e']
hashes = [(Pb.make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_print_trees(phylobuddy, next_hash):
    tester = Pb.trees_to_ascii(phylobuddy)
    assert phylo_to_hash(str(tester)) == next_hash


# ###################### 'pr', '--prune_taxa' ###################### #
pt_hashes = ['99635c6dbf708f94cf4dfdca87113c44', 'fc03b4f100f038277edf6a9f48913dd0', '001db76033cba463a0f187266855e8d5']
pt_hashes = [(Pb.make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(pt_hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", pt_hashes)
def test_prune_taxa(phylobuddy, next_hash):
    Pb.prune_taxa(phylobuddy, 'fir')
    assert phylo_to_hash(phylobuddy) == next_hash


# ######################  'ri', '--rename_ids' ###################### #
hashes = ['6843a620b725a3a0e0940d4352f2036f', '543d2fc90ca1f391312d6b8fe896c59c', '6ce146e635c20ad62e21a1ed6fddbd3a',
          '4dfed97b2a23b8957ee5141bf4681fe4', '77d00fdc512fa09bd1146037d25eafa0', '9b1014be1b38d27f6b7ef73d17003dae']
hashes = [(Pb.make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_rename_ids(phylobuddy, next_hash):
    tester = Pb.rename(phylobuddy, 'Mle', 'Phylo')
    assert phylo_to_hash(tester) == next_hash


def test_rename_nodes():
    tester = Pb.PhyloBuddy(resource("tree_with_node_lables.nwk"))
    Pb.rename(tester, "Equus|Ruminantiamorpha|Canis", "Family")
    assert phylo_to_hash(tester) == "46ba6ce0f3a1859b4c7326a6f5e69263"

# ###################### 'rt', '--root' ###################### #
hashes = ['25ea14c2e89530a0fb48163c0ef2a102', 'e3711259d579cbf0511a5ded66dfd437', '8135cec021240619c27d61288885d8e1',
          '05a83105f54340839dca64a62a22026e', 'f0e26202274a191c9939835b25c1fae4', 'e97c246dc7ebf4d80363f836beff4a81']
hashes = [(Pb.make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_root_middle(phylobuddy, next_hash):
    tester = Pb.root(phylobuddy)
    assert phylo_to_hash(tester) == next_hash

hashes = ['f32bdc34bfe127bb0453a80cf7b01302', 'a7003478d75ad76ef61fcdc643ccdab8', 'c490e3a937b6ee2073c74119984a896e',
          'eacf232776eea70b5de156328e10ecc7', '53caffda3fed5b9004b79effc6d29c36', '3137d568fe07d88620c08480a15006d3']
hashes = [(Pb.make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_root_leaf(phylobuddy, next_hash):
    tester = Pb.root(phylobuddy, "firSA25a")
    assert phylo_to_hash(tester) == next_hash

hashes = ['edcc2400de6b0a4fb05c0a5159215ecd', '09e4c41d22f43b847677eec8be899a72', '89d62bd49d89daa14d2986fe8b826221',
          'b6be77f1d16776554c5a61598ddb6899', '034f6fcb778284bd1bb9db34525f108e', '21d6e13646df1eb75cecb7b10e913eb0']
hashes = [(Pb.make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_root_mrca(phylobuddy, next_hash):
    tester = Pb.root(phylobuddy, "ovi47[ab]", "penIT11b")
    assert phylo_to_hash(tester) == next_hash


# ######################  'su', '--show_unique' ###################### #
def test_show_unique():
    tester = Pb.PhyloBuddy(resource("compare_trees.newick"))
    Pb.show_unique(tester)
    assert phylo_to_hash(tester) == "ea5b0d1fcd7f39cb556c0f5df96281cf"

    with pytest.raises(AssertionError):  # If only a single tree is present
        tester = Pb.PhyloBuddy(Pb.make_copy(pb_objects[0]))
        Pb.show_unique(tester)


def test_show_unique_unrooted():
    tester = Pb.PhyloBuddy(resource("compare_trees.newick"))
    Pb.unroot(tester)
    Pb.show_unique(tester)
    assert phylo_to_hash(tester) == "2bba16e2c77102ba150adecc352407a9"


# ###################### 'sp', '--split_polytomies' ###################### #
def test_split_polytomies():
    tester = Pb.PhyloBuddy('(A,(B,C,D));')
    Pb.split_polytomies(tester)
    assert str(tester) in ['(A:1.0,(B:1.0,(C:1.0,D:1.0):1e-06):1.0):1.0;\n',
                           '(A:1.0,(B:1.0,(D:1.0,C:1.0):1e-06):1.0):1.0;\n',
                           '(A:1.0,(C:1.0,(B:1.0,D:1.0):1e-06):1.0):1.0;\n',
                           '(A:1.0,(C:1.0,(D:1.0,B:1.0):1e-06):1.0):1.0;\n',
                           '(A:1.0,(D:1.0,(B:1.0,C:1.0):1e-06):1.0):1.0;\n',
                           '(A:1.0,(D:1.0,(C:1.0,B:1.0):1e-06):1.0):1.0;\n']


# ###################### 'ur', '--unroot' ###################### #
def test_unroot():
    tester = Pb.PhyloBuddy(resource("figtree.nexus"))
    Pb.unroot(tester)
    assert phylo_to_hash(tester) == "10e9024301b3178cdaed0b57ba33f615"


# ################################################# COMMAND LINE UI ################################################## #
# ###################### argparse_init() ###################### #
def test_argparse_init(capsys):
    sys.argv = ['PhyloBuddy.py', resource("compare_trees.newick")]
    temp_in_args, phylobuddy = Pb.argparse_init()
    assert string2hash(str(phylobuddy)) == "d8e14a2bfc8e9c0ac3c524f5fb478c67"

    sys.argv += ["-f", "foo"]
    with pytest.raises(SystemExit):
        Pb.argparse_init()

    out, err = capsys.readouterr()
    assert "Error: The format 'foo' passed in with the -f flag is not recognized." in err


# ###################### INTERNAL FUNCTIONS ###################### #
def test_print_trees_internal_ui(capsys):
    # Most of this function is covered repeatedly below, so just test the -t flag
    test_in_args = deepcopy(in_args)
    test_in_args.test = True
    test_in_args.unroot = True
    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert err == "*** Test passed ***\n"


def test_in_place_ui(capsys):
    # Some of this function is covered below, so just test an edge
    test_in_args = deepcopy(in_args)
    test_in_args.in_place = True
    test_in_args.trees = [Pb.make_copy(pb_objects[0])]
    test_in_args.screw_formats = "nexus"
    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Warning: The -i flag was passed in, but the positional" in err


# ###################### 'ct', '--consensus_tree' ###################### #
def test_consensus_tree_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.consensus_tree = [False]
    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "f20cbd5aae5971cce8efbda15e4e0b7e"

    test_in_args.consensus_tree = [0.9]
    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "447862ed1ed6e98f2fb535ecce70218b"

    test_in_args.consensus_tree = [1.5]
    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "f20cbd5aae5971cce8efbda15e4e0b7e"


# ###################### 'dt', '--display_trees' ###################### #
@pytest.mark.display
def test_display_trees_ui(capsys, monkeypatch):
    test_in_args = deepcopy(in_args)
    test_in_args.display_trees = True
    show = mock.Mock(return_value=True)
    monkeypatch.setattr(ete3.TreeNode, "show", show)
    Pb.command_line_ui(test_in_args, pb_resources.get_one("o k"), skip_exit=True)

def test_display_trees_ui_no_display(capsys, monkeypatch):
    test_in_args = deepcopy(in_args)
    test_in_args.display_trees = True
    show = mock.Mock(return_value=True)
    monkeypatch.setattr(ete3.TreeNode, "show", show)
    # noinspection PyUnresolvedReferences
    with mock.patch.dict('os.environ'):
        if 'DISPLAY' in os.environ:
            del os.environ['DISPLAY']
        Pb.command_line_ui(test_in_args, pb_resources.get_one("o k"), skip_exit=True)
        out, err = capsys.readouterr()

    assert "Error: Your system is non-graphical, so display_trees can not work." in err


# ###################### 'dis', '--distance' ###################### #
def test_distance_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.distance = [False]
    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "7a9f9902ef30ca3d7ac97b63cfdd0b2e"

    test_in_args.distance = ["uwrf"]
    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "6e520f6911aabdf91dfd075d1545bc1e"

    test_in_args.distance = ["ed"]
    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "f0acdd9c751bce77fe14355cc77b69cc"


# ###################### 'gt', '--generate_tree' ###################### #
@pytest.mark.generate_trees
def test_generate_tree_ui1(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.in_format, test_in_args.out_format = "nexus", "newick"
    test_in_args.trees = [resource("Mnemiopsis_cds.nex")]

    test_in_args.generate_tree = [["fasttree", "-seed 12345"]]
    Pb.command_line_ui(test_in_args, [], skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "d7f505182dd1a1744b45cc326096f70c"


@pytest.mark.generate_trees
def test_generate_tree_ui2():
    test_in_args = deepcopy(in_args)
    test_in_args.in_format, test_in_args.out_format = "nexus", "newick"
    test_in_args.trees = [resource("Mnemiopsis_cds.nex")]

    test_in_args.generate_tree = [["fasttree", "-s 12345"]]  # Command doesn't exist in fasttree
    with pytest.raises(SystemExit):
        Pb.command_line_ui(test_in_args, [])


@pytest.mark.generate_trees
def test_generate_tree_ui3():
    test_in_args = deepcopy(in_args)
    test_in_args.in_format, test_in_args.out_format = "nexus", "newick"
    test_in_args.trees = [resource("Mnemiopsis_cds.nex")]

    test_in_args.generate_tree = [["foo"]]
    with pytest.raises(SystemExit):
        Pb.command_line_ui(test_in_args, [])


# ###################### 'hi', '--hash_ids' ###################### #
def test_hash_ids_ui(capsys, monkeypatch):
    test_in_args = deepcopy(in_args)
    test_in_args.hash_ids = [[1, "nodes"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("o n"), skip_exit=True)
    out, err = capsys.readouterr()

    assert string2hash(out) != phylo_to_hash(pb_resources.get_one("o n"))
    assert "Warning: The hash_length parameter was passed in with the value 1" in err

    test_in_args.hash_ids = [[-1, "nodes"]]

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m n"), skip_exit=True)
    out, err = capsys.readouterr()

    assert string2hash(out) != phylo_to_hash(pb_resources.get_one("m n"))
    assert "Warning: The hash_length parameter was passed in with the value -1" in err

    def hash_ids(*args):
        if args:
            pass
        raise ValueError("Foo bar")

    test_in_args.hash_ids = [[1, "nodes"]]
    monkeypatch.setattr(Pb, "hash_ids", hash_ids)
    with pytest.raises(ValueError):
        Pb.command_line_ui(test_in_args, pb_resources.get_one("o n"))
    out, err = capsys.readouterr()
    assert not err


# ###################### 'li', '--list_ids' ###################### #
def test_list_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.list_ids = [True]

    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "4f0857e0211ff0cd058d0cf7cbaf64d5"

    test_in_args.list_ids = [4]
    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "bf2fbfe1bd52e9b27ae21f5c06e7763a"

    test_in_args.list_ids = [4]
    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy("(, );"), skip_exit=True)
    out, err = capsys.readouterr()
    assert out == "#### tree_1 ####\nNone\n\n"


# ###################### '-nt', '--num_tips' ###################### #
def test_num_tips_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.num_tips = True

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "f43920f4df66e76fbacae4af178eeebb"


# ###################### 'ptr', '--print_trees' ###################### #
def test_print_trees_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.print_trees = True

    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "fe340117cb8f573100c00fc897e6c8ce"


# ###################### 'pt', '--prune_taxa' ###################### #
def test_prune_taxa_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.prune_taxa = [["fir"]]

    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "99635c6dbf708f94cf4dfdca87113c44"

    test_in_args.prune_taxa = [["fir", "ovi"]]

    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[1]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "2a385fa95024323fea412fd2b3c3e91f"


# ###################### 'ri', '--rename_ids' ###################### #
def test_rename_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.rename_ids = ['Mle', 'Phylo']

    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "6843a620b725a3a0e0940d4352f2036f"


# ###################### 'rt', '--root' ###################### #
def test_root_ui_midpoint(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.root = [[]]

    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "25ea14c2e89530a0fb48163c0ef2a102"


def test_root_ui_leaf(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.root = [["firSA25a"]]

    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "f32bdc34bfe127bb0453a80cf7b01302"


def test_root_ui_mrca(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.root = [["firSA25a", "penSH31b"]]

    Pb.command_line_ui(test_in_args, Pb.make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "757489907bd5c82882d084ffcb22cfba"


# ###################### 'sf', '--screw_formats' ###################### #
def test_screw_formats_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = "nexus"

    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "543d2fc90ca1f391312d6b8fe896c59c"


def test_screw_formats_fail(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = "foo"
    Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"), skip_exit=True)
    foo_out, foo_err = capsys.readouterr()
    assert foo_err == "Error: unknown format 'foo'\n"


def test_screw_formats_inplace_ui(capsys):
    temp_file = MyFuncs.TempFile()
    with open(resource("compare_trees.newick"), "r") as ifile:
        temp_file.write(ifile.read())

    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = "nexus"
    test_in_args.in_place = True
    test_in_args.trees[0] = temp_file.path

    tester = Pb.PhyloBuddy(temp_file.path)
    Pb.command_line_ui(test_in_args, tester, skip_exit=True)
    out, err = capsys.readouterr()
    assert "File over-written at:" in err
    check_file = os.path.isfile("%s.nex" % temp_file.path)
    assert check_file

    test_in_args.trees[0] = "%s.nex" % temp_file.path
    test_in_args.screw_formats = "newick"

    tester = Pb.PhyloBuddy("%s.nex" % temp_file.path)
    Pb.command_line_ui(test_in_args, tester, skip_exit=True)
    out, err = capsys.readouterr()
    assert "File over-written at:" in err
    check_file = os.path.isfile("%s.nwk" % temp_file.path)
    assert check_file


# ###################### 'su', '--show_unique' ###################### #
def test_show_unique_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.show_unique = True

    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(resource("compare_trees.newick")), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "ea5b0d1fcd7f39cb556c0f5df96281cf"

    with pytest.raises(SystemExit):
        Pb.command_line_ui(test_in_args, pb_resources.get_one("m k"))

    out, err = capsys.readouterr()
    assert err == "AssertionError: PhyloBuddy object should have exactly 2 trees.\n"


# ###################### 'sp', '--split_polytomies' ###################### #
def test_split_polytomies_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.split_polytomies = True

    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(Pb.PhyloBuddy('(A,(B,C,D));')), skip_exit=True)
    out, err = capsys.readouterr()
    assert out in ['(A:1.0,(B:1.0,(C:1.0,D:1.0):1e-06):1.0):1.0;\n',
                   '(A:1.0,(B:1.0,(D:1.0,C:1.0):1e-06):1.0):1.0;\n',
                   '(A:1.0,(C:1.0,(B:1.0,D:1.0):1e-06):1.0):1.0;\n',
                   '(A:1.0,(C:1.0,(D:1.0,B:1.0):1e-06):1.0):1.0;\n',
                   '(A:1.0,(D:1.0,(B:1.0,C:1.0):1e-06):1.0):1.0;\n',
                   '(A:1.0,(D:1.0,(C:1.0,B:1.0):1e-06):1.0):1.0;\n']


# ###################### 'ur', '--unroot' ###################### #
def test_unroot_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.unroot = True

    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(resource("figtree.nexus")), skip_exit=True)
    out, err = capsys.readouterr()
    assert string2hash(out) == "10e9024301b3178cdaed0b57ba33f615"
