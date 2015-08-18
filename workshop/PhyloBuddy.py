#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Dec 6 2014 

"""
DESCRIPTION OF PROGRAM
"""

import sys
import os
import argparse
import random
import re
from io import StringIO, TextIOWrapper
from collections import OrderedDict
from random import sample
from copy import deepcopy
from tempfile import TemporaryDirectory

# Third party package imports
import Bio.Phylo
from Bio.Phylo import PhyloXML, NeXML, Newick
import ete3
from dendropy.datamodel.treemodel import Tree
from dendropy.datamodel.treecollectionmodel import TreeList


# from newick_utils import *


# ##################################################### WISH LIST #################################################### #
def unroot(_trees):
    pass


def screw_formats(_phylobuddy, _format):
    pass

# Compare two trees, and add colour to the nodes that differ. [ ]

# Implement sum_bootstrap(), but generalize to any value.

# Prune taxa [X]

# Regex taxa names

# List all leaf names [ ]

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


def _stdout(message, quiet=False):
    if not quiet:
        sys.stdout.write(message)
    return


def _format_to_extension(_format):
    format_to_extension = {'fasta': 'fa', 'fa': 'fa', 'genbank': 'gb', 'gb': 'gb', 'nexus': 'nex',
                           'nex': 'nex', 'phylip': 'phy', 'phy': 'phy', 'phylip-relaxed': 'phyr', 'phyr': 'phyr',
                           'stockholm': 'stklm', 'stklm': 'stklm'}
    return format_to_extension[_format]


def _make_copies(_phylobuddy):
    copies = deepcopy(_phylobuddy)
    return copies

def _extract_figtree_metadata(_file_path):
    with open(_file_path, "r") as _tree_file:
        filedata = _tree_file.read()
    extract_fig = re.search('(begin figtree;)', filedata)
    if extract_fig is not None:
        end_regex = re.compile('(end;)')
        end_fig = end_regex.search(filedata, extract_fig.start())
        figdata = filedata[extract_fig.start():end_fig.end()+1]
        filedata = filedata[0:extract_fig.start()-1]
    else:
        return None

    return filedata, figdata

def convert_to_ete(_tree):
    tmp_dir = TemporaryDirectory()
    with open("%s/tree.tmp" % tmp_dir.name, "w") as _ofile:
        _ofile.write(re.sub('!color', 'pb_color', _tree.as_string(schema='newick', annotations_as_nhx=True,
                                                                  suppress_annotations=False, suppress_rooting=True)))

    ete_tree = ete3.TreeNode(newick="%s/tree.tmp" % tmp_dir.name)

    for node in ete_tree.traverse():
        if hasattr(node, 'pb_color'):
            style = ete3.NodeStyle()
            style['fgcolor'] = node.pb_color
            style['hz_line_color'] = node.pb_color
            node.set_style(style)

    return ete_tree

# #################################################################################################################### #
class PhyloBuddy:
    def __init__(self, _input, _in_format=None, _out_format=None):
        # ####  IN AND OUT FORMATS  #### #
        # Holders for input type. Used for some error handling below

        in_handle = None
        raw_seq = None
        in_file = None
        self.trees = []
        tree_classes = [Bio.Phylo.Newick.Tree, Bio.Phylo.NeXML.Tree, Bio.Phylo.PhyloXML.Phylogeny, Tree]
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
            in_handle = _input.read()
            _input.seek(0)

        # File paths
        try:
            if os.path.isfile(_input):
                in_file = _input

        except TypeError:  # This happens when testing something other than a string.
            pass


        if not _in_format:
            self.in_format = guess_format(_input)
            self.out_format = str(self.in_format).lower() if not _out_format else str(_out_format).lower()

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
            self.trees = _input.trees

        elif isinstance(_input, list):
            # make sure that the list is actually Bio.Phylo records (just test a few...)
            _sample = _input if len(_input) < 5 else sample(_input, 5)
            for _tree in _sample:
                if type(_tree) not in tree_classes:
                    raise TypeError("Tree list is not populated with Phylo objects.")
            self.trees = _input

        elif str(type(_input)) == "<class '_io.TextIOWrapper'>" or isinstance(_input, StringIO):
            tmp_dir = TemporaryDirectory()
            with open("%s/tree.tmp" % tmp_dir.name, "w") as _ofile:
                _ofile.write(in_handle)

            figtree = _extract_figtree_metadata("%s/tree.tmp" % tmp_dir.name)
            if figtree is not None:
                with open("%s/tree.tmp" % tmp_dir.name, "w") as _ofile:
                    _ofile.write(figtree[0])

            _trees = Tree.yield_from_files(files=["%s/tree.tmp" % tmp_dir.name], schema=self.in_format,
                                           extract_comment_metadata=True)
            for _tree in _trees:
                self.trees.append(_tree)

        elif os.path.isfile(_input):
            figtree = _extract_figtree_metadata(_input)  # FigTree data being discarded
            if figtree is not None:
                tmp_dir = TemporaryDirectory()
                with open("%s/tree.tmp" % tmp_dir.name, "w") as _ofile:
                    _ofile.write(figtree[0])
                _input = "%s/tree.tmp" % tmp_dir.name
            _trees = Tree.yield_from_files(files=[_input], schema=self.in_format, extract_comment_metadata=True)
            for _tree in _trees:
                self.trees.append(_tree)
        else:
            raise GuessError("Not sure what type this is...")

        # Set 1.0 as length if all lengths are zero or None
        all_none = True
        for _tree in self.trees:
            for _node in _tree:
                if _node.edge_length not in [0, 0.0, None]:
                    all_none = False

        if all_none:
            for _tree in self.trees:
                for _node in _tree:
                    _node.edge_length = 1.0


    def print(self):
        print(self)
        return

    def __str__(self):
        if len(self.trees) == 0:
            return "Error: No trees in object.\n"

        tree_list = TreeList()

        if self.out_format in ["newick", "nexus", "nexml"]:
            for _tree in self.trees:
                tree_list.append(_tree)

        else:
            raise TypeError("Error: Unsupported output format.")

        _output = tree_list.as_string(schema=self.out_format, annotations_as_nhx=True, suppress_annotations=False)
        return _output

    def write(self, _file_path):
        with open(_file_path, "w") as _ofile:
            _ofile.write(str(self))
        return


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
        contents = _input.read()
        if contents == "":
            sys.exit("Input file is empty.")
        _input.seek(0)

        if re.search('<nex:nexml', contents, re.IGNORECASE):
            return 'nexml'
        # Maddison, Swofford, and Maddison, 1997 DOI: 10.1093/sysbio/46.4.590
        elif re.search('#nexus', contents, re.IGNORECASE):
            return 'nexus'
        elif re.search('\(', contents):
            return 'newick'
        else:
            return None  # Unable to determine format from file handle

    else:
        raise GuessError("Unsupported _input argument in guess_format(). %s" % _input)
# #################################################### TOOL KIT ###################################################### #


def split_polytomies(_phylobuddy):
    for _tree in _phylobuddy.trees:
        _tree.resolve_polytomies()
    return _phylobuddy


def root(_phylobuddy, position="guess"):
    for _tree in _phylobuddy.trees:
        outgroup = [{'name': taxon_name} for taxon_name in ('E', 'F', 'G')]
        if _tree.is_monophyletic(outgroup):
            _tree.root_with_outgroup(*outgroup)
        else:
            raise ValueError("outgroup is paraphyletic")
    return _phylobuddy

def prune_taxa(_phylobuddy, _patterns):
    for _tree in _phylobuddy.trees:
        taxa_to_prune = []
        _namespace = _tree.update_taxon_namespace()
        for _taxon in _namespace.labels():
            for _pattern in _patterns:
                if re.search(_pattern, _taxon):
                    taxa_to_prune.append(_taxon)
        for _taxon in taxa_to_prune:
            _tree.prune_taxa_with_labels(StringIO(_taxon))


# ################################################# COMMAND LINE UI ################################################## #
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="phylobuddy", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("trees", help="Supply a file path or raw tree string", nargs="*", default=[sys.stdin])

    parser.add_argument("-spts", "--split_polys", action="store_true",
                        help="Create a binary tree by splitting polytomies randomly.")

    parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')
    parser.add_argument('-q', '--quiet', help="Suppress stderr messages", action='store_true')
    parser.add_argument('-t', '--test', action='store_true',
                        help="Run the function and return any stderr/stdout other than sequences.")
    parser.add_argument('-prt', '--prune_taxa', action='append', nargs="+")
    parser.add_argument('-ptr', '--print_trees', action='store_true')
    parser.add_argument('-o', '--out_format', help="If you want a specific format output", action='store')
    parser.add_argument('-f', '--in_format', help="If PhyloBuddy can't guess the file format, just specify it directly.",
                        action='store')

    in_args = parser.parse_args()

    phylobuddy = []
    tree_set = ""

    for tree_set in in_args.trees:
        if isinstance(tree_set, TextIOWrapper) and tree_set.buffer.raw.isatty():
            sys.exit("Warning: No input detected. Process will be aborted.")
        tree_set = PhyloBuddy(tree_set, in_args.in_format, in_args.out_format)
        phylobuddy += tree_set.trees
    phylobuddy = PhyloBuddy(phylobuddy, tree_set.in_format, tree_set.out_format)


# ################################################ INTERNAL FUNCTIONS ################################################ #
    def _print_trees(_phylobuddy):  # TODO: Remove the calls to in_args
        if in_args.test:
            _stderr("*** Test passed ***\n", in_args.quiet)
            pass

        elif in_args.in_place:
            _in_place(str(_phylobuddy), in_args.trees[0])

        else:
            _stdout("{0}\n".format(str(_phylobuddy).rstrip()))


    def _in_place(_output, _path):
        if not os.path.exists(_path):
            _stderr("Warning: The -i flag was passed in, but the positional argument doesn't seem to be a "
                    "file. Nothing was written.\n", in_args.quiet)
            _stderr("%s\n" % _output.strip(), in_args.quiet)
        else:
            with open(os.path.abspath(_path), "w") as _ofile:
                _ofile.write(_output)
            _stderr("File over-written at:\n%s\n" % os.path.abspath(_path), in_args.quiet)

# ############################################## COMMAND LINE LOGIC ############################################## #

    # Split polytomies
    if in_args.split_polys:
        split_polytomies(phylobuddy)
        _print_trees(phylobuddy)
        sys.exit()

    # Print trees
    if in_args.print_trees:
        _print_trees(phylobuddy)
        sys.exit()

    # Prune taxa
    if in_args.prune_taxa:
        prune_taxa(phylobuddy, in_args.prune_taxa[0])
        _print_trees(phylobuddy)
        sys.exit()