#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Dec 6 2014 

"""
DESCRIPTION OF PROGRAM
"""

import sys
import os
import argparse
from io import StringIO
from random import sample

# Third party package imports
import Bio.Phylo
import Bio.Phylo.Newick
import Bio.Phylo.NeXML
import Bio.Phylo.PhyloXML
from newick_utils import *


# ##################################################### WISH LIST #################################################### #
def unroot(_trees):
    pass


def screw_formats(_phylobuddy, _format):
    pass

# Compare two trees, and add colour to the nodes that differ.

# Implement sum_bootstrap(), but generalize to any value.

# Prune taxa

# Regex taxa names

# List all leaf names

# 'Clean' a tree, as implemented in phyutility

# See http://cegg.unige.ch/system/files/nwutils_tutorial.pdf for ideas
# Re-implement many or all of Phyultility commands: https://code.google.com/p/phyutility/


# ################################################# HELPER FUNCTIONS ################################################# #
def _stderr(message, quiet=False):
    if not quiet:
        sys.stderr.write(message)
    return


class GuessError(Exception):
    """Raised when input format cannot be guessed"""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


# #################################################################################################################### #
class PhyloBuddy:
    def __init__(self, _input, _in_format=None, _out_format=None):
        # ####  IN AND OUT FORMATS  #### #
        # Holders for input type. Used for some error handling below
        in_handle = None
        raw_seq = None
        in_file = None
        tree_classes = [Bio.Phylo.Newick.Tree, Bio.Phylo.NeXML.Tree, Bio.Phylo.PhyloXML.Phylogeny]

        # Handles
        if str(type(_input)) == "<class '_io.TextIOWrapper'>":
            if not _input.seekable():  # Deal with input streams (e.g., stdout pipes)
                temp = StringIO(_input.read())
                _input = temp
            _input.seek(0)
            in_handle = _input.read()
            _input.seek(0)

        # Raw sequences
        if type(_input) == str and not os.path.isfile(_input):
            raw_seq = _input
            temp = StringIO(_input)
            _input = temp
            _input.seek(0)

        # File paths
        try:
            if os.path.isfile(_input):
                in_file = _input

        except TypeError:  # This happens when testing something other than a string.
            pass

        if not _in_format:
            self.in_format = guess_format(_input)
            self.out_format = str(self.in_format) if not _out_format else _out_format

        else:
            self.in_format = _in_format

        if not self.in_format:
            if in_file:
                raise GuessError("Could not determine format from _input file '{0}'.\n"
                                 "Try explicitly setting with -f flag.".format(in_file))
            elif raw_seq:
                raise GuessError("Could not determine format from raw input\n{0} ..."
                                 "Try explicitly setting with -f flag.".format(raw_seq)[:50])
            elif in_handle:
                raise GuessError("Could not determine format from input file-like object\n{0} ..."
                                 "Try explicitly setting with -f flag.".format(in_handle)[:50])
            else:
                raise GuessError("Unable to determine format or input type. Please check how PhyloBuddy is being called.")

        self.out_format = self.in_format if not _out_format else _out_format

        # ####  RECORDS  #### #
        if str(type(_input)) == "<class '__main__.PhyloBuddy'>":
            _trees = _input.trees

        elif isinstance(_input, list):
            # make sure that the list is actually Bio.Phylo records (just test a few...)
            _sample = _input if len(_input) < 5 else sample(_input, 5)
            for _tree in _sample:
                if type(_tree) not in tree_classes:
                    raise TypeError("Tree list is not populated with Phylo objects.")
            _trees = _input

        elif str(type(_input)) == "<class '_io.TextIOWrapper'>" or isinstance(_input, StringIO):
            _trees = list(Bio.Phylo.parse(_input, self.in_format))

        elif os.path.isfile(_input):
            with open(_input, "r") as _input:
                _trees = list(Bio.Phylo.parse(_input, self.in_format))
        else:
            raise GuessError("Not sure what type this is...")

        self.trees = _trees


def guess_format(_input):
    # If input is just a list, there is no BioPython in-format. Default to Newick.
    if isinstance(_input, list):
        return "newick"

    # Pull value directly from object if appropriate
    if type(_input) == PhyloBuddy:
        return _input.in_format

    # If input is a handle or path, try to read the file in each format, and assume success if not error and # trees > 0
    if os.path.isfile(str(_input)):
        _input = open(_input, "r")

    if str(type(_input)) == "<class '_io.TextIOWrapper'>" or isinstance(_input, StringIO):
        # Die if file is empty
        if _input.read() == "":
            sys.exit("Input file is empty.")
        _input.seek(0)

        possible_formats = ["newick", "nexus", "nexml", "phyloxml"]  # Look into cdao in future
        for _format in possible_formats:
            _input.seek(0)
            _seqs = Bio.Phylo.parse(_input, _format)
            if next(_seqs):
                _input.seek(0)
                return _format
            else:
                continue

        return None  # Unable to determine format from file handle

    else:
        raise GuessError("Unsupported _input argument in guess_format(). %s" % _input)
# #################################################### TOOL KIT ###################################################### #


def split_polytomies(_trees):
    _output = []
    for _tree in _trees.trees:
        # print(Tree.parse_newick_input(_tree))
        print(_tree)
    return _trees


def root(_phylobuddy, position="guess"):
    for tree in _phylobuddy.trees:
        outgroup = [{'name': taxon_name} for taxon_name in ('E', 'F', 'G')]
        if tree.is_monophyletic(outgroup):
            tree.root_with_outgroup(*outgroup)
        else:
            raise ValueError("outgroup is paraphyletic")
    return _phylobuddy

# ################################################# COMMAND LINE UI ################################################## #
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="phylobuddy", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("trees", help="Supply a file path or raw tree string", nargs="+")

    parser.add_argument("-spts", "--split_polys", action="store_true",
                        help="Create a binary tree by splitting polytomies randomly.")

    parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')
    parser.add_argument('-q', '--quiet', help="Suppress stderr messages", action='store_true')
    parser.add_argument('-t', '--test', action='store_true',
                        help="Run the function and return any stderr/stdout other than sequences.")
    parser.add_argument('-o', '--out_format', help="If you want a specific format output", action='store')
    parser.add_argument('-f', '--in_format', help="If PhyloBuddy can't guess the file format, just specify it directly.",
                        action='store')

    in_args = parser.parse_args()

    phylobuddy = []
    tree_set = ""

    for tree_set in in_args.trees:
        tree_set = PhyloBuddy(tree_set, in_args.in_format)
        phylobuddy += tree_set.trees

    phylobuddy = PhyloBuddy(phylobuddy)
    phylobuddy.out_format = in_args.out_format if in_args.out_format else tree_set.out_format

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

        if in_args.in_place:
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

# ############################################## COMMAND LINE LOGIC ############################################## #

    if in_args.split_polys:
        # _print_trees(split_polytomies(trees))
        split_polytomies(phylobuddy)
