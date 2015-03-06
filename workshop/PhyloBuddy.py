#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Dec 6 2014 

"""
DESCRIPTION OF PROGRAM
"""

from Bio import Phylo
import sys
import os
import argparse
from newick_utils import *


# ##################################################### WISH LIST #################################################### #
def unroot(_trees):
    x = 1


def root(_trees, position="guess"):
    for tree in _trees.trees:
        outgroup = [{'name': taxon_name} for taxon_name in ('E', 'F', 'G')]
        if tree.is_monophyletic(outgroup):
            tree.root_with_outgroup(*outgroup)
        else:
            raise ValueError("outgroup is paraphyletic")
    return _trees

def screw_formats(_trees, _format):
    x = 1

# Compare two trees, and add colour to the nodes that differ.

# Implement sum_bootstrap(), but generalize to any value.

# Prune taxa

# Regex taxa names

# See http://cegg.unige.ch/system/files/nwutils_tutorial.pdf for ideas
# Re-implement many or all of Phyultility commands: https://code.google.com/p/phyutility/
# ################################################ INTERNAL FUNCTIONS ################################################ #
def _print_trees(_trees):  # TODO: Remove the calls to in_args
    if len(_trees.trees) == 0:
        sys.stderr.write("Nothing returned.\n")
        return False

    _output = ""
    for _tree in _trees.trees:
        try:
            _output += _tree.format(_trees.out_format) + "\n"
        except ValueError as e:
            sys.stderr.write("Error: %s\n" % e)

    if in_args.in_place and in_place_allowed:
        sys.exit("In-place writing not implemented yet")
        if not os.path.exists(in_args.sequence[0]):
            sys.stderr.write("Warning: The -i flag was passed in, but the positional argument doesn't seem to be a "
                             "file. Nothing was written.\n")
            sys.stdout.write("%s\n" % _output.strip())
        else:
            with open(os.path.abspath(in_args.sequence[0]), "w") as ofile:
                ofile.write(_output)
            sys.stderr.write("File over-written at:\n%s\n" % os.path.abspath(in_args.sequence[0]))
    else:
        sys.stdout.write("%s\n" % _output.strip())


# ################################################# HELPER FUNCTIONS ################################################# #
class TreePreparer():
    def __init__(self, _input, _in_format=None, _out_format=None):
        self.in_format = "newick"
        self.out_format = "newick"

        if str(type(_input)) == "<class '__main__.TreePreparer'>":
            _trees = _input.trees

        elif isinstance(_input, list):
            # make sure that the list is actually Phylo records (just test a few...)
            for _tree in _input[:3]:
                if str(type(_tree)).split(".")[-1] != "Tree'>":
                    sys.exit("Error: TreeList is not populated with Tree objects.")
            _trees = _input

        elif str(type(_input)) == "<class '_io.TextIOWrapper'>":
            _trees = list(Phylo.parse(_input, self.in_format))

        elif os.path.isfile(_input):
            _input = open(_input, "r")
            _trees = list(Phylo.parse(_input, self.in_format))
            _input.close()
        else:
            # _trees = [SeqRecord(Seq(_input))]
            sys.exit("Raw tree not implemented in TreePreparer yet.")

        self.trees = _trees
# #################################################################################################################### #


def split_polytomies(_trees):
    _output = []
    for _tree in _trees.trees:
        print(Tree.parse_newick_input(_tree))
    return _trees

# ################################################# COMMAND LINE UI ################################################## #
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="phylobuddy", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("trees", help="Supply a file path or raw tree string", nargs="+")

    parser.add_argument("-spts", "--split_polys", action="store_true",
                        help="Create a binary tree by splitting polytomies randomly.")

    parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')
    parser.add_argument('-f', '--in_format', help="If PhyloBuddy can't guess the file format, just specify it directly.",
                        action='store')

    in_args = parser.parse_args()

    in_place_allowed = False

    trees = []
    tree_set = ""
    for tree_set in in_args.trees:
        tree_set = TreePreparer(tree_set, in_args.in_format)
        trees += tree_set.trees

    trees = TreePreparer(trees)

    if in_args.split_polys:
        # _print_trees(split_polytomies(trees))
        split_polytomies(trees)