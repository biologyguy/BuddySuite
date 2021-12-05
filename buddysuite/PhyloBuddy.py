#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: PhyloBuddy.py
author: Stephen R. Bond
email: biologyguy@gmail.com
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description:
PhyloBuddy is a general wrapper for popular phylogenetic programs, handles format conversion, and manipulates tree files
"""

# ##################################################### IMPORTS ###################################################### #
from __future__ import print_function

# BuddySuite specific
try:
    import buddy_resources as br
    import AlignBuddy as Alb
except ImportError:
    try:
        import buddysuite.buddy_resources as br
        import buddysuite.AlignBuddy as Alb
    except AttributeError:
        from . import buddy_resources as br
        from . import AlignBuddy as Alb

# Standard library
import sys
import os
import random
import string
import re
import shutil
import webbrowser
from math import log, ceil
from io import StringIO, TextIOWrapper
from subprocess import Popen, CalledProcessError, check_output, PIPE
from collections import OrderedDict
from copy import deepcopy

# Third party
import Bio.Phylo
# from Bio.Phylo import PhyloXML, NeXML, Newick
# sys.path.insert(0, "./")  # For stand alone executable, where dependencies are packaged with BuddySuite
try:
    import dendropy
except ImportError:
    print("""\
Dendropy not detected on your system. Try running the following:

    pip install dendropy

Or see https://pythonhosted.org/DendroPy/ for installation details.
""")
    sys.exit()

from dendropy.datamodel.treemodel import Tree, Node
from dendropy.datamodel.treecollectionmodel import TreeList
from dendropy.datamodel.taxonmodel import TaxonNamespace
from dendropy.calculate import treecompare

# Patch a bug in older versions of DendroPy (V4.1.0 and older)
dendropy_ver = int(re.sub("[^0-9]", "", dendropy.__version__))
if dendropy_ver < 411:
    from dendropy.datamodel.basemodel import Annotable


    def deep_copy_annotations_from(self, other, memo=None):
        if hasattr(other, "_annotations"):
            # if not isinstance(self, other.__class__) or not isinstance(other, self.__class__):
            if type(self) is not type(other):
                raise TypeError(
                    "Cannot deep-copy annotations from different type (unable to assume "
                    "object equivalence in dynamic or nested annotations)")
            if memo is None:
                memo = {}
            for a1 in other._annotations:
                a2 = deepcopy(a1, memo=memo)
                memo[id(a1)] = a2
                if a2.is_attribute and a1._value[0] is other:
                    a2._value = (self, a1._value[1])
                self.annotations.add(a2)
            if hasattr(self, "_annotations"):
                memo[id(other._annotations)] = self._annotations

    Annotable.deep_copy_annotations_from = deep_copy_annotations_from


# ##################################################### WISH LIST #################################################### #
"""
def modify_metadata(_trees):
    # Probably just loop through dendropy nodes and modify node.annotations["<attr>"] 
    return _trees

def decode_accessions(phylobuddy):
    # If taxa lables are accessions, reach out to the respective database and resolve them into actual names
    return phylobuddy

"""


# - Compare two trees, and add colour to the nodes that differ. [ ]
# - Implement sum_bootstrap(), but generalize to any value.
# - Regex taxa names
# - 'Clean' a tree, as implemented in phyutility
# - Pull random tree, pull every nth tree, and pull specific tree
# - Pull random tips
# - Prune clade (ie, prune away all tips from the common ancestor of two or more tips)
# See http://cegg.unige.ch/system/files/nwutils_tutorial.pdf for ideas
# - Re-implement many or all of Phyultility commands: https://code.google.com/p/phyutility/

# ##################################################### GLOBALS ###################################################### #
CONFIG = br.config_values()
VERSION = br.Version("PhyloBuddy", 1, "4.0", br.contributor_list, {"year": 2021, "month": 12, "day": 5})
OUTPUT_FORMATS = ["newick", "nexus", "nexml"]
PHYLO_INFERENCE_TOOLS = ["raxml", "phyml", "fasttree"]


# #################################################### PHYLOBUDDY #################################################### #
class PhyloBuddy(object):
    def __init__(self, _input, _in_format=None, _out_format=None):
        # ####  IN AND OUT FORMATS  #### #
        # Holders for input type. Used for some error handling below
        tmp_dir = br.TempDir()
        in_from_handle = None
        raw_seq = None
        in_file = None
        self.trees = []
        self.hash_map = []  # Only used when hash_ids() function is called
        tree_classes = [Tree]  # Dendropy Tree
        # Handles
        if str(type(_input)) == "<class '_io.TextIOWrapper'>":
            if not _input.seekable():  # Deal with input streams (e.g., stdout pipes)
                input_txt = _input.read()
                if re.search("Buddy::.* has crashed with the following traceback", input_txt):
                    print(input_txt)
                    sys.exit()
                temp = StringIO(br.utf_encode(input_txt))
                _input = temp
            _input.seek(0)
            in_from_handle = _input.read()
            _input.seek(0)

        # Raw sequences
        if type(_input) == str and not os.path.isfile(_input):
            raw_seq = br.utf_encode(_input)
            temp = StringIO(raw_seq)
            _input = temp
            in_from_handle = _input.read()
            _input.seek(0)

        # File paths
        try:
            if os.path.isfile(_input):
                in_file = _input
                with open(in_file, "r", encoding="utf-8") as ifile:
                    _input = StringIO(ifile.read())
                    in_from_handle = _input.read()
                    _input.seek(0)
        except TypeError:  # This happens when testing something other than a string.
            pass

        if not _in_format:
            self.in_format = _guess_format(_input)
            self.out_format = str(self.in_format).lower() if not _out_format else str(_out_format).lower()

        else:
            self.in_format = _in_format

        if not self.in_format:
            if in_file:
                raise br.GuessError("Could not automatically determine the format of '{0}'.\n"
                                    "Try explicitly setting it with the -f flag.".format(in_file))
            elif raw_seq:
                raise br.GuessError("Could not automatically determine the format from raw input\n{0} ..."
                                    "Try explicitly setting it with the -f flag.".format(raw_seq[:50]))
            elif in_from_handle:
                raise br.GuessError("Could not automatically determine the format from input file-like object\n{0} ..."
                                    "Try explicitly setting it with the -f flag.".format(in_from_handle[:50]))
            else:
                raise br.GuessError("Unable to determine the format or input type. "
                                    "Please check how PhyloBuddy is being called.")

        self.out_format = self.in_format if not _out_format else _out_format

        # ####  RECORDS  #### #
        if _input.__class__.__name__ == "PhyloBuddy":
            self.trees = _input.trees

        elif isinstance(_input, list):
            # make sure that the list is actually Bio.Phylo records (just test a few...)
            _sample = _input if len(_input) < 5 else random.sample(_input, 5)
            for _tree in _sample:
                if type(_tree) not in tree_classes:
                    raise TypeError("Tree list is not populated with Phylo objects.")
            self.trees = _input

        elif str(type(_input)) == "<class '_io.TextIOWrapper'>" or isinstance(_input, StringIO):
            with open(os.path.join(tmp_dir.path, "tree.tmp"), "w", encoding="utf-8") as _ofile:
                in_from_handle = clean_newick(in_from_handle)
                _ofile.write(in_from_handle)

            # Removes figtree data so parser doesn't die
            figtree = _extract_figtree_metadata(os.path.join(tmp_dir.path, "tree.tmp"))
            if figtree:
                in_from_handle = figtree[0]
#                with open("%s/tree.tmp" % tmp_dir.path, "w", encoding="utf-8") as _ofile:
#                    _ofile.write(figtree[0])

            _trees = dendropy.TreeList()
            if self.in_format != 'nexml':
                _trees.read(data=in_from_handle, schema=self.in_format, extract_comment_metadata=True,
                            preserve_underscores=True)
            else:
                _trees.read(data=in_from_handle, schema=self.in_format)

            for _tree in _trees:
                self.trees.append(_tree)

        else:
            raise br.GuessError("Not sure what type this is...")

        # Set 1.0 as length if all lengths are zero or None and order features
        all_none = True
        for _tree in self.trees:
            for _node in _tree.nodes():
                if _node.has_annotations:
                    _node.annotations._item_list = sorted(_node.annotations._item_list, key=lambda x: x.name)
                    _node.annotations._item_set = set(_node.annotations._item_list)
                if _node.edge_length not in [0, 0.0, None]:
                    all_none = False

        if all_none:
            for _tree in self.trees:
                for _node in _tree.nodes():
                    _node.edge_length = 1.0
        self.memory_footprint = len(str(self))

    def __str__(self):
        if len(self.trees) == 0:
            return "Error: No trees in object.\n"

        tree_list = TreeList()

        if self.out_format in OUTPUT_FORMATS:
            for _tree in self.trees:
                tree_list.append(_tree)

        else:
            raise TypeError("Error: Unsupported output format.")

        if self.out_format != 'nexml':
            _output = tree_list.as_string(schema=self.out_format, annotations_as_nhx=False, suppress_annotations=False)
        else:
            _output = tree_list.as_string(schema='nexml')

        _output = '{0}\n'.format(_output.rstrip())

        return _output

    def __len__(self):
        return len(self.trees)

    def write(self, _file_path, out_format=None):
        orig_out_format = str(self.out_format)
        self.out_format = self.out_format if out_format is None else out_format
        with open(_file_path, "w", encoding="utf-8") as _ofile:
            _ofile.write(str(self))
        self.out_format = orig_out_format
        return


# ################################################# HELPER FUNCTIONS ################################################# #
def _extract_figtree_metadata(_file_path):
    """
    Removes the figtree block from nexus files
    :param _file_path: Specifies the nexus file path
    :return: A length 2 tuple containing the nexus data and the figtree block
    """
    with open(_file_path, "r", encoding="utf-8") as _tree_file:
        filedata = _tree_file.read()
    extract_fig = re.search('(begin figtree;)', filedata)
    if extract_fig is not None:
        end_regex = re.compile('(end;)')
        end_fig = end_regex.search(filedata, extract_fig.start())
        figdata = filedata[extract_fig.start():end_fig.end() + 1]
        filedata = filedata[0:extract_fig.start() - 1]
    else:
        return None
    return filedata, figdata


def _get_tree_binaries(_tool):
    """
    Returns a URL where a tool's binaries can be found.
    :param _tool: Specify the tree building tool to be used.
    :return: A string containing the tool's URL.
    """

    tool_dict = {'raxml': 'http://sco.h-its.org/exelixis/web/software/raxml/index.html',
                 'phyml': 'http://www.atgc-montpellier.fr/phyml/versions.php',
                 'fasttree': 'http://www.microbesonline.org/fasttree/#Install'}
    return tool_dict[_tool]


def _guess_format(_input):
    # If input is just a list, there is no BioPython in-format. Default to Newick.
    if isinstance(_input, list):
        return "newick"

    # Pull value directly from object if appropriate
    if _input.__class__.__name__ == "PhyloBuddy":
        return _input.in_format

    # If input is a handle or path, try to read the file in each format, and assume success if not error and # trees > 0
    if os.path.isfile(str(_input)):
        _input = open(_input, "r", encoding='utf-8')

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
        raise br.GuessError("Unsupported _input argument in guess_format(). %s" % _input)


def clean_newick(tree_string):
    tree_string = re.sub(r':([0-9]\.*[0-9]*)\[([0-9]+)\]', r'\2:\1', tree_string)
    return tree_string


def make_copy(_phylobuddy):
    """
    Returns a copy of the PhyloBuddy object
    :param _phylobuddy: The PhyloBuddy object to be copied
    :return: A copy of the original PhyloBuddy object
    """
    _copy = deepcopy(_phylobuddy)
    return _copy


def num_taxa(phylobuddy, nodes=False, split=False):
    count = [0]
    for indx, tree in enumerate(phylobuddy.trees):
        if split and indx > 0:
            count.append(0)
        for node in tree:
            if nodes and node.label:
                count[-1] += 1
            if node.taxon and node.taxon.label:
                count[-1] += 1
    return count if split else count[0]


# ################################################ MAIN API FUNCTIONS ################################################ #
def add_branch(phylobuddy, new_branch, sister_taxa):
    """
    Add a subtree to an existing tree
    :param phylobuddy: PhyloBuddy object of existing tree
    :param new_branch: Single taxon or PhyloBuddy object of subtree
    :param sister_taxa: Single taxon or list of taxa that the new branch should be placed sister to (regex patterns)
    :return:
    """
    if type(new_branch) == str:
        new_branch = PhyloBuddy("(%s);" % new_branch)
    elif type(new_branch) == Tree:
        new_branch = PhyloBuddy([new_branch])
    elif type(new_branch) != PhyloBuddy:
        raise TypeError("new_branch parameter must be a string, dentropy tree, or PhyloBuddy object.")

    for node in new_branch.trees[0]:
        if node.parent_node is None:
            new_branch = node
            if new_branch.num_child_nodes() == 1:
                new_branch = new_branch.child_nodes()[0]
            break

    sister_taxa = [sister_taxa] if type(sister_taxa) == str else sister_taxa

    for tree in phylobuddy.trees:
        all_nodes = []
        for regex in sister_taxa:
            for indx, id_list in list_ids(PhyloBuddy([tree])).items():
                for next_id in id_list:
                    if re.search(regex, next_id):
                        all_nodes.append(next_id)
        if len(all_nodes) == 0:
            raise AttributeError("Unable to identify any sister taxa in tree.")

        elif len(all_nodes) == 1:
            leaf_node = tree.find_node_with_taxon_label(all_nodes[0])
            mrca = leaf_node.parent_node
            new_edge_len = leaf_node.edge.length / 2
            attachment = Node(edge_length=new_edge_len)
            leaf_node.edge.length = new_edge_len
            new_branch.edge.length = new_edge_len
            mrca.remove_child(leaf_node)
            mrca.add_child(attachment)
            attachment.add_child(new_branch)
            attachment.add_child(leaf_node)

        else:
            mrca = tree.mrca(taxon_labels=all_nodes)
            new_edge_len = mrca.edge.length / 2
            attachment = Node(edge_length=new_edge_len)
            mrca.edge.length = new_edge_len
            new_branch.edge.length = new_edge_len

            mrca_parent = mrca.parent_node
            if mrca_parent is not None:
                mrca_parent.remove_child(mrca)
                mrca_parent.add_child(attachment)
                attachment.add_child(mrca)
            else:
                mrca.add_child(attachment)

            attachment.add_child(new_branch)
    return phylobuddy


def collapse_polytomies(phylobuddy, threshold, mode="support"):
    """
    Remove nodes if their support value or branch length are below the given threshold
    :param phylobuddy: PhyloBuddy objecct
    :param threshold: Value below which the node will be removed
    :type threshold: float
    :param mode: Choose which value to apply threshold to (either 'support' or 'length')
    :return: The modified PhyloBuddy object
    """
    if mode not in ["support", "length"]:
        raise NameError("Mode must be 'support' or 'length'")
    for tree in phylobuddy.trees:
        tree.encode_bipartitions()
        for node in tree.postorder_node_iter():
            try:
                value = node.label if mode == "support" else node.edge_length
                if float(value) < threshold:
                    for child in node.child_nodes():
                        child.parent_node = node.parent_node
                    if node.parent_node:  # Ensure not working from the root
                        node.parent_node.remove_child(node)
                        tree.encode_bipartitions()
            except TypeError:
                continue
    return phylobuddy


def consensus_tree(phylobuddy, frequency=.5):
    """
    Create a consensus tree from two or more trees.
    :param phylobuddy: PhyloBuddy object
    :param frequency: The frequency threshold of a node for it to be included in the new tree.
    :return: The modified PhyloBuddy object
    """
    _trees = TreeList(phylobuddy.trees)
    _consensus = _trees.consensus(min_freq=frequency,
                                  add_edge_length_summaries_as_edge_attributes=False,
                                  add_edge_length_summaries_as_edge_annotations=False)
    phylobuddy.trees = [_consensus]
    return phylobuddy


def display_trees(phylobuddy, program=None):
    """
    Displays trees in a web browser window, one-by-one.
    :param phylobuddy: PhyloBuddy object
    :param program: How should the tree be displayed --> ['figtree', 'system']
    :return: None
    """
    available_programs = ['figtree', 'system']
    if program and program not in available_programs:
        raise AttributeError("Unknown program '%s' selected for display. "
                             "Please select between %s" % (program, available_programs))

    if os.name != "nt" and sys.platform != "darwin" and "DISPLAY" not in os.environ:
        # We just assume that a Windows/Mac machine is graphical
        raise SystemError("This system does not appear to be graphical, "
                          "so display_trees() will not work. Try using trees_to_ascii()")

    program = "system" if program is None else program

    if program == "figtree":
        if not shutil.which("figtree"):
            raise SystemError("Unable to find FigTree on your system. Please install it and ensure it is present in"
                              " your $PATH. http://tree.bio.ed.ac.uk/software/figtree/")

    tmp_dir = br.TempDir()
    tree_file = tmp_dir.subfile('tree.svg')

    if program == "figtree":
        if os.name == "nt":
            print("Sorry! I'm not sure how to launch FigTree on Windows. Given that you've seen this message, you"
                  " probably want the functionality, so leave a message in the GitHub issue tracker.")
            return False

        phylobuddy.write(tree_file, out_format="nexus")
        tty = open("/dev/tty", mode="r")
        os.dup2(tty.fileno(), 0)

        Popen("figtree %s" % tree_file, shell=True, stdout=PIPE, stderr=PIPE)
        sys.stdin = input("Switching to FigTree, hit Enter to exit PhyloBuddy.\n")
        tty.close()

    elif program == "system":
        try:
            br.dummy_func()  # For monkey patching purposes in unit tests
            import pylab
        except RuntimeError as err:
            if "The Mac OS X backend will not be able to function correctly" in str(err):
                message = """\
It looks like you are on a Mac with an incorrect backend
for matplotlib. This can probably be fixed by amending your 
matplotlibrc file. Would you like PhyloBuddy try to do this 
automatically? y/[n]: """
                if br.ask(message, default="no"):
                    home = os.path.expanduser('~')
                    mplibrc = os.path.join(home, ".matplotlib", "matplotlibrc")
                    os.makedirs(os.path.join(home, ".matplotlib"), exist_ok=True)
                    if os.path.isfile(mplibrc):
                        with open(mplibrc, "r") as ifile:
                            matplotlibrc = ifile.read().strip()
                    else:
                        matplotlibrc = ""

                    with open(mplibrc, "w") as ofile:
                        ofile.write("# Backend set by PhyloBuddy\nbackend: TkAgg\n%s" % matplotlibrc)

                    br._stderr("\n%s amended.\nPlease re-run your command.\n" % mplibrc)
                else:
                    br._stderr("\nOriginal error message:\n\n" + str(err) + "\n\nStack Overflow may help you out:\n"
                               "https://stackoverflow.com/questions/21784641/installation-issue-with-"
                               "matplotlib-python\n\n")
                sys.exit()

            else:
                raise err

        label_colors = []
        indx = 0
        for tree in phylobuddy.trees:
            label_colors.append(OrderedDict())
            for node in tree:
                if node.taxon and node.annotations.get_value('!color'):
                    label_colors[indx][node.taxon.label] = node.annotations.get_value('!color')
            indx += 1

        indx = 0
        tty = open("/dev/tty", mode="r")
        os.dup2(tty.fileno(), 0)

        for tree in Bio.Phylo.parse(StringIO(str(phylobuddy)), phylobuddy.out_format):
            Bio.Phylo.draw(tree, label_func=lambda leaf: leaf.name, do_show=False, label_colors=label_colors[indx])
            pylab.axis('off')
            pylab.savefig(tree_file, format='svg', bbox_inches='tight', dpi=300)
            if os.name == "nt":  # File path specification different between operating systems
                file_location = tree_file
            else:
                file_location = "file:///" + tree_file
            webbrowser.open(file_location, new=2)
            input("Switching to system browser. Hit Enter to continue.\n")
            indx += 1
        tty.close()
    return True


def distance(phylobuddy, method='weighted_robinson_foulds'):
    """
    Calculates distance metrics between pairs of trees
    :param phylobuddy: PhyloBuddy object
    :param method: The tree comparison method ([un]weighted_robinson_foulds/euclidean_distance)
    :return: A dictionary of dictonaries containing the distances between tree pairs. dict[tree1][tree2]
    """
    if not len(phylobuddy.trees) > 1:
        raise ValueError("Distance requires at least two trees.")
    method = method.lower()
    if method in ['wrf', 'weighted_robinson_foulds']:
        method = 'wrf'
    elif method in ['uwrf', 'sym', 'unweighted_robinson_foulds', 'symmetric']:
        method = 'uwrf'
    # elif _method in ['mgk', 'mason_gamer_kellogg']:
    #     _method = 'mgk'
    elif method in ['ed', 'euclid', 'euclidean', 'euclidean_distance']:
        method = 'euclid'
    else:
        raise AttributeError('{0} is an invalid comparison method.'.format(method))

    output = OrderedDict()
    for indx1, tree1 in enumerate(phylobuddy.trees):  # Compares all-by-all
        key1 = 'tree_{0}'.format(indx1 + 1) if tree1.label in [None, 'None', ''] else tree1.label
        output.setdefault(key1, OrderedDict())
        for indx2, tree2 in enumerate(list(phylobuddy.trees)[indx1 + 1:]):
            indx2 += indx1 + 1
            key2 = 'tree_{0}'.format(indx2 + 1) if tree2.label in [None, 'None', ''] else tree2.label
            output.setdefault(key2, OrderedDict())

            if method == 'wrf':
                output[key1][key2] = treecompare.weighted_robinson_foulds_distance(tree1, tree2)

            elif method == 'uwrf':
                output[key1][key2] = treecompare.symmetric_difference(tree1, tree2)

            # elif _method == 'mgk':
            #     _output[_key1][_key2] = treecompare.mason_gamer_kellogg_score(_tree1, _tree2)

            else:
                output[key1][key2] = treecompare.euclidean_distance(tree1, tree2)
    return output


def generate_tree(alignbuddy, alias, params=None, keep_temp=None, quiet=False, r_seed=None):
    # ToDo: Break the function up for each program being wrapped. There's WAY too much going on here...
    # ToDo: Multiple alignments are not un-hashing correctly (at least in raxml)
    """
    Calls tree building tools to generate trees
    :param alignbuddy: The AlignBuddy object containing the alignments for building the trees
    :param alias: The tree building tool to be used (raxml/phyml/fasttree/iqtree)
    :param params: Additional parameters to be passed to the tree building tool
    :param keep_temp: Determines if/where the temporary files will be kept
    :param quiet: Suppress all output form alignment programs
    :param r_seed: Set random seed used for name hashing
    :return: A PhyloBuddy object containing the trees produced.
    """

    if params is None:
        params = ''

    if keep_temp:  # Store files in temp dir in a non-temporary directory
        if os.path.exists(keep_temp):
            raise FileExistsError("Execution of %s was halted to prevent files in '%s' from being overwritten. Please "
                                  "choose another location to save temporary files."
                                  % (alias, os.path.split(keep_temp)[-1]))

    # Figure out what tool is being used
    tool = False
    if alias in ['raxml', 'phyml', 'fasttree']:
        tool = alias
    elif "raxml" in alias.lower():
        tool = "raxml"
    elif "phyml" in alias.lower():
        tool = "phyml"
    elif "fasttree" in alias.lower():
        tool = "fasttree"
    elif "iqtree" in alias.lower():
        tool = "iqtree"
    else:
        for prog in [('raxml', " -v", "This is RAxML version"),
                     ('phyml', " --version", "This is PhyML version"),
                     ('fasttree', "", "Usage for FastTree version"),
                     ('iqtree', " -h", "IQ-TREE")]:
            version = Popen("%s%s" % (alias, prog[1]), shell=True, stderr=PIPE, stdout=PIPE).communicate()
            if prog[2] in version[0].decode() or prog[2] in version[1].decode():
                tool = prog[0]
                break
    if not tool:
        raise AttributeError("{0} is not a recognized tree inference tool. "
                             "Please check your spelling (case sensitive)".format(alias))

    if shutil.which(alias) is None:  # Tool must be callable from command line
        raise ProcessLookupError('#### Could not find {0} on your system. ####\n'
                                 'Please check that your spelling is correct (case sensitive) or find installation '
                                 'instructions at {1}.\n'.format(alias, _get_tree_binaries(tool)))

    else:
        tmp_dir = br.TempDir()
        tmp_in = tmp_dir.subfile("pb_input.aln")
        cwd = os.getcwd()

        def remove_invalid_params(_dict):  # Helper method for blacklisting flags
            parameters = params
            for key in _dict:
                if _dict[key] is True:  # Flag has an argument
                    _pattern = "{0} [^-]*".format(key)
                    # Deletes until it reaches another flag.
                    # May cause issues if flag is right before the file path.
                else:
                    _pattern = "{0}".format(key)
                parameters = re.sub(_pattern, '', parameters)
            return parameters

        params = re.split(' ', params, )  # Expands paths in _params to absolute paths
        file_list = []
        for indx, token in enumerate(params):
            if os.path.exists(token):
                params[indx] = os.path.abspath(token)
                if os.path.isfile(token):
                    file_list.append(params[indx])

        params = ' '.join(params)

        phylo_objs = []
        for alignment in alignbuddy.alignments:  # Need to loop through one tree at a time
            sub_alignbuddy = Alb.AlignBuddy([alignment])
            sub_alignbuddy.hash_map = alignbuddy.hash_map
            Alb.hash_ids(sub_alignbuddy, 8, r_seed=r_seed)
            sub_alignbuddy = Alb.clean_seq(sub_alignbuddy)
            sub_alignbuddy.set_format('phylipss') if tool == "phyml" else sub_alignbuddy.set_format('phylipi')

            sub_alignbuddy.write(os.path.join(tmp_dir.path, "pb_input.aln"))  # Most tree builders require an input file

            for next_file in file_list:
                with open(next_file, "r", encoding="utf-8", errors='ignore') as ifile:
                    contents = ifile.read()
                for _hash, _id in sub_alignbuddy.hash_map.items():
                    contents = re.sub(r'%s([,)\'\":])' % _id, r'%s\1' % _hash, contents)
                    contents = re.sub(r'[\'\"]%s[\'\"]' % _hash, r'%s' % _hash, contents)
                file_name = os.path.split(next_file)[1]
                file_name = tmp_dir.subfile(file_name)
                with open(file_name, "w", encoding="utf-8") as ofile:
                    ofile.write(contents)
                params = re.sub(next_file, file_name, params)

            if tool == 'raxml':
                params = remove_invalid_params({'-s': True, '-n': True, '-w': True})
                if '-T' not in params:  # Num threads
                    params += ' -T 2'
                if '-m' not in params:  # An evolutionary model is required
                    if sub_alignbuddy.alpha in ["DNA", "RNA"]:
                        br._stderr("Warning: Using default evolutionary model GTRCAT\n", quiet=quiet)
                        params += " -m GTRCAT"
                    elif sub_alignbuddy.alpha == "protein":
                        br._stderr("Warning: Using default evolutionary model PROTCATLG\n", quiet=quiet)
                        params += " -m PROTCATLG"

                if '-p' not in params:  # RNG seed
                    params += ' -p 12345' if not r_seed else ' -p %s' % r_seed
                if '-#' not in params and '-N' not in params:  # Number of trees to build
                    params += ' -# 1'
                command = '{0} -s {1} {2} -n result -w {3}'.format(alias, tmp_in, params, tmp_dir.path)

            elif tool == 'phyml':
                params = remove_invalid_params({'-q': False, '--sequential': False, '-u': True, '--inputtree': True,
                                                '--run_id': True})
                if sub_alignbuddy.alpha in ["DNA", "RNA"]:
                    if '-d nt' not in params and '--datatype nt' not in params:
                        params += ' -d nt'

                elif sub_alignbuddy.alpha == "protein":
                    if '-d aa' not in params and '--datatype aa' not in params:
                        params += ' -d aa'

                command = '{0} -i {1} {2}'.format(alias, tmp_in, params)

            elif tool == 'fasttree':
                params = remove_invalid_params({'-n': True})
                if sub_alignbuddy.alpha in ["DNA", "RNA"]:
                    command = '{0} {1} -nt {2}'.format(alias, params, tmp_in)  # FastTree must be told what alpha to use
                else:
                    command = '{0} {1} {2}'.format(alias, params, tmp_in)

            elif tool == "iqtree":
                params = remove_invalid_params({"-s": True, "-st": True})
                if sub_alignbuddy.alpha in ["DNA", "RNA"]:
                    params += ' -st DNA'
                else:
                    params += ' -st AA'

                if "-nt" not in params:
                    params += ' -nt AUTO'

                command = '{0} -s {1} {2}'.format(alias, tmp_in, params)
                os.chdir(tmp_dir.path)
            output = ''

            try:
                if tool in ['raxml', 'phyml', 'iqtree']:  # If tool writes to file
                    if quiet:
                        Popen(command, shell=True, universal_newlines=True, stdout=PIPE, stderr=PIPE).communicate()
                    else:
                        Popen(command, shell=True, universal_newlines=True, stdout=sys.__stderr__).wait()
                    file_found = False
                    outfiles = {'raxml': ['RAxML_bestTree.result', 'RAxML_bootstrap.result',
                                          'RAxML_bipartitions.result', 'RAxML_fastTreeSH_Support.result',
                                          'RAxML_MajorityRuleExtendedConsensusTree.result'],
                                'phyml': ['pb_input.aln_phyml_tree', 'pb_input.aln_phyml_tree.txt'],
                                'iqtree': ['pb_input.aln.treefile']}
                    for path in [os.path.join(tmp_dir.path, x) for x in outfiles[tool]]:
                        if os.path.isfile(path):
                            file_found = True
                            break
                    if not file_found:
                        raise FileNotFoundError("Error: {0} failed to generate a tree.".format(alias))
                else:  # If tool outputs to stdout
                    if quiet:
                        output = check_output(command, shell=True, universal_newlines=True, stderr=PIPE)
                    else:
                        output = check_output(command, shell=True, universal_newlines=True)

            except CalledProcessError:  # Haven't been able to find a way to get here. Needs a test.
                raise RuntimeError('{0} threw an error. Scroll up for more info.\n'.format(alias))
            if tool == 'raxml':  # Pull tree from written file
                num_runs = re.search('-[#N] ([0-9]+)', params)
                num_runs = 0 if not num_runs else int(num_runs.group(1))
                if re.search('-b [0-9]+', params) or re.search('-x [0-9]+', params):
                    with open(os.path.join(tmp_dir.path, "RAxML_bootstrap.result"), "r", encoding="utf-8") as result:
                        output += result.read()
                elif os.path.isfile(os.path.join(tmp_dir.path, "RAxML_bipartitions.result")):
                    with open(os.path.join(tmp_dir.path, "RAxML_bipartitions.result"), "r", encoding="utf-8") as result:
                        output += result.read()
                elif os.path.isfile(os.path.join(tmp_dir.path, "RAxML_bestTree.result")):
                    with open(os.path.join(tmp_dir.path, "RAxML_bestTree.result"), "r", encoding="utf-8") as result:
                        output += result.read()
                elif os.path.isfile(os.path.join(tmp_dir.path, "RAxML_fastTreeSH_Support.result")):
                    with open(os.path.join(tmp_dir.path, "RAxML_fastTreeSH_Support.result"), "r", encoding="utf-8") \
                            as result:
                        output += result.read()
                elif os.path.isfile(os.path.join(tmp_dir.path, "RAxML_MajorityRuleExtendedConsensusTree.result")):
                    with open(os.path.join(tmp_dir.path, "RAxML_MajorityRuleExtendedConsensusTree.result"), "r",
                              encoding="utf-8") as result:
                        result = result.read()
                        result = re.sub(":1\.0", "", result)
                        result = re.sub(r'\[([0-9]+)\]', r'\1', result)
                        output += result
                elif num_runs > 1:
                    for tree_indx in range(num_runs):
                        with open(os.path.join(tmp_dir.path, "RAxML_result.result.RUN.{}".format(tree_indx)),
                                  "r", encoding="utf-8") as result:
                            output += result.read()

            elif tool == 'phyml':
                phyml_out = os.path.join(tmp_dir.path, "pb_input.aln_phyml_tree")
                phyml_out += '.txt' if not os.path.isfile(phyml_out) else ''
                with open(phyml_out, "r", encoding="utf-8") as result:
                    output += result.read()
            elif tool == 'iqtree':
                with open(os.path.join(tmp_dir.path, 'pb_input.aln.treefile'), "r") as result:
                    output += result.read()

            os.chdir(cwd)

            if keep_temp:  # Store temp files
                shutil.copytree(tmp_dir.path, keep_temp)

            phylobuddy = PhyloBuddy(output)

            for _hash, _id in sub_alignbuddy.hash_map.items():
                rename(phylobuddy, _hash, _id)

            if keep_temp:
                _root, dirs, files = next(br.walklevel(keep_temp))
                for file in files:
                    with open(os.path.join(_root, file), "r", encoding="utf-8", errors='ignore') as ifile:
                        contents = ifile.read()
                    for _hash, _id in sub_alignbuddy.hash_map.items():
                        contents = re.sub(_hash, _id, contents)
                    with open(os.path.join(_root, file), "w", encoding="utf-8") as ofile:
                        ofile.write(contents)
            phylo_objs += phylobuddy.trees

        phylobuddy = PhyloBuddy(phylo_objs)
        br._stderr("Returning to PhyloBuddy...\n\n", quiet)
        return phylobuddy


def hash_ids(phylobuddy, hash_length=10, nodes=False, r_seed=None):
    """
    Replaces the sequence IDs with random hashes
    :param phylobuddy: PhyloBuddy object
    :param hash_length: Specifies the length of the new hashed IDs
    :param nodes: Also hash node labels
    :param r_seed: Set the random generator seed value
    :return: The modified PhyloBuddy object, with a new attribute `hash_map` added
    """

    class HashFactory(object):
        def __init__(self):
            self.hash_map = []
            # It seems that Dendropy does not create unique labels for each node/tip if the labels are the same, instead
            # it shares the object among everything with the same label name (even between trees). The all_hashes dict
            # allows me to account for this.
            self.all_hashes = {}

        def new_hash(self, label):
            if str(label) in self.all_hashes:
                self.hash_map[-1][label] = self.all_hashes[label]
                return str(label)

            while True:
                output_hash = "".join([rand_gen.choice(string.ascii_letters + string.digits)
                                       for _ in range(hash_length)])
                if output_hash in self.all_hashes:
                    continue
                else:
                    self.hash_map[-1][output_hash] = str(label)
                    self.all_hashes[output_hash] = str(label)
                    return output_hash

        def add_tree(self):
            self.hash_map.append(OrderedDict())
            return

    try:
        hash_length = int(hash_length)
    except ValueError:
        raise TypeError("Hash length argument must be an integer, not %s" % type(hash_length))

    if hash_length < 1:
        raise ValueError("Hash length must be greater than 0")

    if 32 ** hash_length <= num_taxa(phylobuddy, nodes) * 2:
        raise ValueError("Insufficient number of hashes available to cover all sequences. "
                         "Hash length must be increased.")

    rand_gen = random.Random() if not r_seed else random.Random(r_seed)
    hashes = HashFactory()
    for tree in phylobuddy.trees:
        hashes.add_tree()
        for node in tree:
            if nodes and node.label:
                node.label = hashes.new_hash(str(node.label))

            if node.taxon and node.taxon.label:
                node.taxon.label = hashes.new_hash(str(node.taxon.label))

    phylobuddy.hash_map = hashes.hash_map
    return phylobuddy


def ladderize(phylobuddy, ascending=True):
    for tree in phylobuddy.trees:
        tree.ladderize(ascending)
    return phylobuddy


def list_ids(phylobuddy):
    """
    Returns a dictionary of tree names and node labels
    :param phylobuddy: The PhyloBuddy object to be analyzed
    :return: A dictionary of tree names and node labels
    """
    output = OrderedDict()
    for indx, tree in enumerate(phylobuddy.trees):
        namespace = TaxonNamespace()
        for node in tree:
            if node.taxon:
                namespace.add_taxon(node.taxon)
        key = tree.label if tree.label not in [None, ''] else 'tree_{0}'.format(str(indx + 1))
        output[key] = list(namespace.labels())
    return output


def prune_taxa(phylobuddy, *patterns):
    """
    Prunes taxa that match one or more regex patterns
    :param phylobuddy: The PhyloBuddy object whose trees will be pruned.
    :param patterns: One or more regex patterns.
    :return: The same PhyloBuddy object after pruning.
    """
    for indx, tree in enumerate(phylobuddy.trees):
        taxa_to_prune = []
        namespace = TaxonNamespace()
        for node in tree:  # Populate the namespace for easy iteration
            if node.taxon:
                namespace.add_taxon(node.taxon)
        for taxon in namespace.labels():
            for pattern in patterns:
                if re.search(pattern, taxon):  # Sets aside the names of the taxa to be pruned
                    taxa_to_prune.append(taxon)
        try:
            for taxon in taxa_to_prune:  # Removes the nodes from the tree
                tree.prune_taxa_with_labels(StringIO(taxon))
        except AttributeError:  # If all leaves are removed, then delete the tree
            del phylobuddy.trees[indx]


def rename(phylobuddy, query, replace):
    """
    Substitutes matches in node names with a string
    :param phylobuddy: PhyloBuddy object
    :param query: The regex pattern to be searched
    :param replace: The string to replace the matches with
    :return: The modified PhyloBuddy object
    """
    for indx, tree in enumerate(phylobuddy.trees):
        for node in tree:
            if node.label:
                node.label = re.sub(query, replace, node.label)
            if node.taxon and node.taxon.label:
                node.taxon.label = re.sub(query, replace, node.taxon.label)
    return phylobuddy


def root(phylobuddy, *root_nodes):
    """
    Place a new root on trees
    Note: There may be some sort of weird behavior (likely crash) if an internal node is selected as root. Not sure how
    to test for it though, so it isn't currently accounted for. Need a use case.
    :param phylobuddy: PhyloBuddy object
    :param root_nodes: A string with the taxon label to root on, or a list of labels and the most common ancestor
    node will be rooted on
    :return: The modified PhyloBuddy object
    """
    def _root(_tree, _root_nodes=None):
        if _root_nodes:
            all_nodes = []
            for regex in root_nodes:
                for indx, id_list in list_ids(PhyloBuddy([_tree])).items():
                    for next_id in id_list:
                        if re.search(regex, next_id):
                            all_nodes.append(next_id)
            if len(all_nodes) == 0:
                return _tree
            elif len(all_nodes) == 1:
                leaf_node = _tree.find_node_with_taxon_label(all_nodes[0])
                mrca = leaf_node._parent_node

            else:
                mrca = _tree.mrca(taxon_labels=all_nodes)

            _tree.reroot_at_node(mrca, update_bipartitions=True, suppress_unifurcations=False)

        else:
            # WARNING! There is a bug in DendroPy leading to an infinite loop here. The dendropy folks have fixed it
            # in their development branch but it is not yet in the main branch
            _tree.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)

    for tree in phylobuddy.trees:
        _root(tree, root_nodes)
        tree.is_rooted = True

    return phylobuddy


def show_unique(phylobuddy):
    """
    Colors shared leaves green and unique leaves red (between two trees)
    :param phylobuddy: PhyloBuddy object
    :return: The labeled PhyloBuddy object
    """
    if len(phylobuddy.trees) != 2:
        raise AssertionError("PhyloBuddy object should have exactly 2 trees.")

    tree1_leaves = [node.taxon.label for node in phylobuddy.trees[0] if node.taxon]
    tree2_leaves = [node.taxon.label for node in phylobuddy.trees[1] if node.taxon]
    unique = []
    for leaf in tree1_leaves:
        if leaf not in tree2_leaves and leaf not in unique:
            unique.append(leaf)

    for leaf in tree2_leaves:
        if leaf not in tree1_leaves and leaf not in unique:
            unique.append(leaf)

    # Add color
    for node in phylobuddy.trees[0]:
        if node.taxon:
            if node.taxon.label in unique:
                node.annotations['!color'] = '#ff0000'
            else:
                node.annotations['!color'] = '#00ff00'

    for node in phylobuddy.trees[1]:
        if node.taxon:
            if node.taxon.label in unique:
                node.annotations['!color'] = '#ff0000'
            else:
                node.annotations['!color'] = '#00ff00'

    return phylobuddy


def split_polytomies(phylobuddy, r_seed=None):
    """
    Randomly splits polytomies. This function was drawn almost verbatim from the DendroPy Tree.resolve_polytomies()
    method. The main difference is that a very small edge_length is assigned to new nodes when a polytomy is split.
    :param phylobuddy: PhyloBuddy object
    :param r_seed: Specify a random seed so results can be reproducible
    :return: Modified PhyloBuddy object
    """
    rand = random.Random(r_seed)
    for tree in phylobuddy.trees:
        polytomies = []
        for node in tree.postorder_node_iter():
            if len(node._child_nodes) > 2:
                polytomies.append(node)
        for node in polytomies:
            to_attach = rand.sample(node._child_nodes, len(node._child_nodes) - 2)
            for child in to_attach:
                node.remove_child(child)
            attachment_points = list(node._child_nodes)
            while len(to_attach) > 0:
                next_child = to_attach.pop()
                next_sib = rand.choice(attachment_points)
                next_attachment = Node(edge_length=0.000001)
                p = next_sib._parent_node
                p.add_child(next_attachment)
                p.remove_child(next_sib)
                next_attachment.add_child(next_sib)
                next_attachment.add_child(next_child)
                attachment_points.append(next_attachment)
                attachment_points.append(next_child)

    return phylobuddy


def trees_to_ascii(phylobuddy):
    """
    Returns an ascii representation of the tree. Scales to terminal window size.
    :param phylobuddy: The PhyloBuddy object whose trees will be converted.
    :return: A string containing ASCII representations of the trees.
    """
    output = OrderedDict()
    for indx, tree in enumerate(phylobuddy.trees):
        key = 'tree_{0}'.format(indx + 1) if tree.label in [None, ''] else tree.label
        output[key] = tree.as_ascii_plot()
    return output


def unroot(phylobuddy):
    """
    Convert to unrooted trees
    :param phylobuddy: PhyloBuddy object
    :return: The modified PhyloBuddy object
    """
    for tree in phylobuddy.trees:
        tree.is_rooted = False
        tree.update_bipartitions()
    return phylobuddy


# ################################################# COMMAND LINE UI ################################################## #
def argparse_init():
    import argparse

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

    br.flags(parser, ("trees", "Supply file path(s) or raw tree string. If piping trees into PhyloBuddy "
                               "this argument must be left blank."), br.pb_flags, br.pb_modifiers, VERSION)

    in_args = parser.parse_args()
    br.check_garbage_flags(in_args, "PhyloBuddy")

    phylobuddy = []
    tree_set = ""

    if in_args.in_format and in_args.in_format.lower() not in OUTPUT_FORMATS:
        br._stderr("Error: The format '%s' passed in with the -f flag is not recognized. "
                   "Valid options include %s.\n" % (in_args.in_format, OUTPUT_FORMATS))
        sys.exit()

    if in_args.out_format and in_args.out_format.lower() not in OUTPUT_FORMATS:
        br._stderr("Error: Output type %s is not recognized/supported\n" % in_args.out_format)
        sys.exit()

    if not in_args.generate_tree:  # If passing in an alignment, don't want to try and build PhyloBuddy obj
        for tree_set in in_args.trees:
            if isinstance(tree_set, TextIOWrapper) and tree_set.buffer.raw.isatty():
                br._stderr("Warning: No input detected so PhyloBuddy is aborting...\n"
                           "For more information, try:\n%s --help\n" % sys.argv[0])
                sys.exit()
            tree_set = PhyloBuddy(tree_set, in_args.in_format, in_args.out_format)
            phylobuddy += tree_set.trees
        phylobuddy = PhyloBuddy(phylobuddy, tree_set.in_format, tree_set.out_format)

    in_args.random_seed = None if not in_args.random_seed else in_args.random_seed
    return in_args, phylobuddy


def command_line_ui(in_args, phylobuddy, skip_exit=False, pass_through=False):   # ToDo: Convert to a class
    # ############################################## INTERNAL FUNCTIONS ############################################## #
    def _print_trees(_phylobuddy):
        if in_args.test:
            br._stderr("*** Test passed ***\n", in_args.quiet)

        elif in_args.in_place:
            _in_place(str(_phylobuddy), in_args.trees[0])

        else:
            br._stdout("{0}\n".format(str(_phylobuddy).rstrip()))

    def _in_place(_output, file_path):
        if not os.path.isfile(str(file_path)):
            br._stderr("Warning: The -i flag was passed in, but the positional argument doesn't seem to be a "
                       "file. Nothing was written.\n", in_args.quiet)
            br._stderr("%s\n" % _output.strip(), in_args.quiet)
        else:
            with open(os.path.abspath(file_path), "w", encoding="utf-8") as _ofile:
                _ofile.write(_output)
            br._stderr("File overwritten at:\n%s\n" % os.path.abspath(file_path), in_args.quiet)

    def _exit(_tool, skip=skip_exit):
        if skip:
            return
        usage = br.Usage()
        memory_footprint = 0 if type(phylobuddy) != PhyloBuddy else phylobuddy.memory_footprint
        usage.increment("PhyloBuddy", VERSION.short(), _tool, memory_footprint)
        usage.save()
        sys.exit()

    def _raise_error(_err, _tool, check_string=None):
        if pass_through:
            raise _err
        if check_string:
            if type(check_string) == str:
                check_string = [check_string]
            re_raise = True
            for _string in check_string:
                if re.search(_string, str(_err)):
                    re_raise = False
                    break
            if re_raise:
                raise _err
        br._stderr("{0}: {1}\n".format(_err.__class__.__name__, str(_err)), in_args.quiet)
        _exit(_tool)

    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Add branch
    if in_args.add_branch:
        args = in_args.add_branch[0]
        if len(args) < 2:
            _raise_error(AttributeError("Add branch tool requires at least two arguments: new_branch and sister_taxa."),
                         "add_branch")
        try:
            new_branch = PhyloBuddy(args[0])
        except br.GuessError:
            new_branch = args[0]
        sister_taxa = args[1:]

        try:
            _print_trees(add_branch(phylobuddy, new_branch, sister_taxa))
        except AttributeError as err:
            _raise_error(err, "add_branch", "Unable to identify any sister taxa in tree")
        _exit("add_branch")

    # Collapse polytomies
    if in_args.collapse_polytomies:
        args = in_args.collapse_polytomies[0]
        mode = "support"
        threshold = 50
        for arg in args:
            try:
                threshold = float(arg)
            except ValueError:
                mode = arg.lower() if arg.lower() in ["support", "length"] else mode

        phylobuddy = collapse_polytomies(phylobuddy, threshold=threshold, mode=mode)
        _print_trees(phylobuddy)
        _exit("collapse_polytomies")

    # Consensus tree
    if in_args.consensus_tree:
        frequency = in_args.consensus_tree[0]
        frequency = 0.5 if not frequency else frequency
        if not 0 <= frequency <= 1:
            br._stderr("Warning: The frequency value should be between 0 and 1. Defaulting to 0.5.\n\n")
            frequency = 0.5

        _print_trees(consensus_tree(phylobuddy, frequency))
        _exit("consensus_tree")

    # Display trees
    if in_args.display_trees:
        program = None if in_args.display_trees[0] is None else in_args.display_trees[0].lower()
        program = "figtree" if program and "figtree".startswith(program) else program
        try:
            display_trees(phylobuddy, program)
        except SystemError as err:
            print(sys.platform, str(err))
            if sys.platform == "darwin" and "graphical" in str(err):
                br._stderr("Error: Your system does not appear to be graphical. "
                           "Try installing XQuartz (https://www.xquartz.org/), or use print_trees instead.")
            elif "graphical" in str(err):
                br._stderr("Error: Your system does not appear to be graphical, so display_trees can not work. "
                           "Please use print_trees instead.")
            else:
                _raise_error(err, "display_trees", "FigTree")

        except AttributeError as err:
            _raise_error(err, "display_trees", "Unknown program")

        _exit("display_trees")

    # Distance
    if in_args.distance:
        if not len(phylobuddy.trees) > 1:
            _raise_error(ValueError("Distance requires at least two trees."), "distance")
        if in_args.distance[0]:
            output = distance(phylobuddy, in_args.distance[0])
        else:
            output = distance(phylobuddy)

        br._stderr('Tree 1\tTree 2\tValue\n')
        for key1 in output:
            for key2 in output[key1]:
                br._stdout('{0}\t{1}\t{2}\n'.format(key1, key2, output[key1][key2]))
        _exit("distance")

    # Generate Tree
    if in_args.generate_tree:
        # ToDo: The extra arguments parameter probably doesn't need to be dependent on the tool parameter being passed
        args = in_args.generate_tree[0]
        if not args:
            for tool in ['raxml', 'phyml', 'fasttree', 'iqtree']:
                if shutil.which(tool):
                    args = [tool]
                    break
        if not args:
            _raise_error(AttributeError("Unable to identify any supported phylogenetic inference "
                                        "software on your system."), "generate_alignment")

        alignbuddy = []
        align_set = None

        out_format = None if not in_args.out_format else str(in_args.out_format)
        in_args.out_format = None if not in_args.out_format else "phylipsr"

        for align_set in in_args.trees:  # Build an AlignBuddy object
            if isinstance(align_set, TextIOWrapper) and align_set.buffer.raw.isatty():
                br._stderr("Warning: No input detected so PhyloBuddy is aborting...\n"
                           "For more information, try:\n%s --help\n" % sys.argv[0])
                sys.exit()
            align_set = Alb.AlignBuddy(align_set, in_args.in_format, in_args.out_format)
            alignbuddy += align_set.alignments
        if align_set:
            alignbuddy = Alb.AlignBuddy(alignbuddy, align_set.in_format, align_set.out_format)
        else:
            alignbuddy = Alb.AlignBuddy(alignbuddy, in_args.in_format, in_args.out_format)

        params = re.sub("\[(.*)\]", "\1", args[1]) if len(args) > 1 else None

        generated_trees = None
        try:
            generated_trees = generate_tree(alignbuddy, args[0], params, in_args.keep_temp, quiet=in_args.quiet,
                                            r_seed=in_args.random_seed)
        except (FileExistsError, AttributeError, ProcessLookupError, RuntimeError) as e:
            _raise_error(e, "generate_tree")
        except FileNotFoundError as e:
            _raise_error(e, "generate_tree", "Error: {0} failed to generate a tree.".format(args[0]))

        if in_args.out_format:
            generated_trees.out_format = out_format
        phylobuddy = generated_trees
        _print_trees(phylobuddy)
        _exit("generate_tree")

    # Hash sequence ids
    if in_args.hash_ids:
        args = in_args.hash_ids[0]
        hash_length = 10
        hash_nodes = False
        if args:
            for arg in args:
                try:
                    hash_length = int(arg)
                except ValueError:
                    if arg == "nodes":
                        hash_nodes = True

        if hash_length < 1:
            br._stderr("Warning: The hash_length parameter was passed in with the value %s. This is not a positive "
                       "integer, so the hash length as been set to 10.\n\n" % hash_length, quiet=in_args.quiet)
            hash_length = 10

        try:
            hash_ids(phylobuddy, hash_length=hash_length, nodes=hash_nodes, r_seed=in_args.random_seed)
        except ValueError as e:
            if "Insufficient number of hashes available" in str(e):
                holder = ceil(log(num_taxa(phylobuddy) * 2, 32))
                br._stderr("Warning: The hash_length parameter was passed in with the value %s. "
                           "This is too small to properly cover all sequences, so it has been increased to %s.\n\n" %
                           (hash_length, holder), in_args.quiet)
                hash_length = int(holder)
                hash_ids(phylobuddy, hash_length=hash_length, nodes=hash_nodes, r_seed=in_args.random_seed)
            else:
                raise e

        hash_table = "##### Hash table #####\n"
        for indx, tree_map in enumerate(phylobuddy.hash_map):
            if len(phylobuddy.hash_map) > 1:
                hash_table += "# Tree %s\n" % (indx + 1)
            for _hash, orig_id in tree_map.items():
                hash_table += "%s\t%s\n" % (_hash, orig_id)
            hash_table += "\n"
        hash_table = "%s\n######################\n\n" % hash_table.strip()
        br._stderr(hash_table, in_args.quiet)
        _print_trees(phylobuddy)
        _exit("hash_ids")

    # Ladderize
    if in_args.ladderize:
        ascending = True if not in_args.ladderize[0] or "rev" not in in_args.ladderize[0].lower() else False
        phylobuddy = ladderize(phylobuddy, ascending)
        _print_trees(phylobuddy)
        _exit("ladderize")

    # List ids
    if in_args.list_ids:
        listed_ids = list_ids(phylobuddy)
        columns = 1 if not in_args.list_ids[0] or in_args.list_ids[0] <= 0 else abs(in_args.list_ids[0])
        output = ""
        for key in listed_ids:
            count = 0
            output += '#### {0} ####\n'.format(key)
            if len(listed_ids[key]) == 0:
                output += 'None\n'
            else:
                for identifier in listed_ids[key]:
                    if count % columns == 0:
                        output = "%s\n" % output.strip()

                    output += "%s\t" % identifier
                    count += 1

            output = '%s\n\n' % output.strip()
        br._stdout('%s\n\n' % output.strip())
        _exit("list_ids")

    # Number of tips
    if in_args.num_tips:
        counts = num_taxa(phylobuddy, split=True)
        output = ""
        for indx, count in enumerate(counts):
            if len(counts) > 1:
                output += "# Tree %s\n" % (indx + 1)
            output += "%s\n\n" % count
        br._stdout("%s\n" % output.strip())

    # Print trees
    if in_args.print_trees:
        tree_table = trees_to_ascii(phylobuddy)
        output = ''
        for key in tree_table:
            output += '\n#### {0} ####\n'.format(key)
            output += tree_table[key]
        output += '\n'
        br._stdout(output)
        _exit("print_trees")

    # Prune taxa
    if in_args.prune_taxa:
        patterns = br.clean_regex(in_args.prune_taxa[0], in_args.quiet)
        prune_taxa(phylobuddy, *patterns)
        _print_trees(phylobuddy)
        _exit("prune_taxa")

    # Rename IDs
    if in_args.rename_ids:
        patterns = br.clean_regex(in_args.rename_ids[0], in_args.quiet)
        if patterns:
            _print_trees(rename(phylobuddy, patterns[0], in_args.rename_ids[1]))
        else:
            _print_trees(phylobuddy)
        _exit("rename_ids")

    # Root
    if in_args.root:
        root_nodes = br.clean_regex(in_args.root[0], in_args.quiet)
        if root_nodes:
            root(phylobuddy, *root_nodes)
        else:
            root(phylobuddy)
        _print_trees(phylobuddy)
        _exit("root")

    # Screw formats
    if in_args.screw_formats:
        if in_args.screw_formats not in OUTPUT_FORMATS:
            br._stderr("Error: unknown format '%s'\n" % in_args.screw_formats)
        else:
            phylobuddy.out_format = in_args.screw_formats
            if in_args.in_place and os.path.isfile(str(in_args.trees[0])):
                # Need to change the file extension
                _path, ext = os.path.splitext(os.path.abspath(in_args.trees[0]))
                _path = "%s.%s" % (_path, br.format_to_extension[phylobuddy.out_format])

                os.remove(in_args.trees[0])
                in_args.trees[0] = _path
                open(in_args.trees[0], "w", encoding="utf-8").close()

            _print_trees(phylobuddy)
        _exit("screw_formats")

    # Show unique
    if in_args.show_unique:
        try:
            _print_trees(show_unique(phylobuddy))
        except AssertionError as e:
            _raise_error(e, "show_unique", "PhyloBuddy object should have exactly 2 trees.")
        _exit("show_unique")

    # Split polytomies
    if in_args.split_polytomies:
        split_polytomies(phylobuddy, r_seed=in_args.random_seed)
        _print_trees(phylobuddy)
        _exit("split_polytomies")

    # Unroot
    if in_args.unroot:
        unroot(phylobuddy)
        _print_trees(phylobuddy)
        _exit("unroot")


def main():
    br.preparse_flags()
    initiation = []
    try:
        initiation = argparse_init()  # initiation = [in_agrs, phylobuddy]
        command_line_ui(*initiation)
    except (KeyboardInterrupt, br.GuessError) as _e:
        print(_e)
        return False
    except SystemExit:
        return False
    except Exception as _e:
        func = ""
        for next_arg in vars(initiation[0]):
            if getattr(initiation[0], next_arg) and next_arg in br.pb_flags:
                func = next_arg
                break
        br.send_traceback("PhyloBuddy", func, _e, VERSION)
        return False
    return True


if __name__ == '__main__':
    main()
