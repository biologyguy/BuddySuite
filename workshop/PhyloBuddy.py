#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: PhyloBuddy.py
version: 1, beta
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description:
PhyloBuddy is a general wrapper for popular phylogenetic programs, handles format conversion, and manipulates tree files
"""

# ##################################################### IMPORTS ###################################################### #
# Standard library
import sys
import os
import random
import re
import shutil
from io import StringIO, TextIOWrapper
from subprocess import Popen, CalledProcessError, check_output
from collections import OrderedDict
from random import sample
from copy import deepcopy

# Third party
# import Bio.Phylo
# from Bio.Phylo import PhyloXML, NeXML, Newick
sys.path.insert(0, "./")  # For stand alone executable, where dependencies are packaged with BuddySuite
from Bio.Alphabet import IUPAC

try:
    import ete3
    from ete3.coretype.tree import TreeError
except ImportError:
    confirm = input("PhyloBuddy requires ETE v3+, which was not detected on your system. Try to install [y]/n? ")
    if confirm.lower() in ["", "y", "yes"]:
        Popen("pip install --upgrade  https://github.com/jhcepas/ete/archive/3.0.zip", shell=True).wait()
        Popen("pip install six", shell=True).wait()
        try:
            import ete3
            from ete3.coretype.tree import TreeError

        except ImportError:
            sys.exit("Failed to install ETE3, please see http://etetoolkit.org/download/ for further details")
    else:
        sys.exit("Aborting. Please see http://etetoolkit.org/download/ for installation details\n")

try:
    import dendropy
except ImportError:
    confirm = input("PhyloBuddy requires dendropy, which was not detected on your system. Try to install [y]/n? ")
    if confirm.lower() in ["", "y", "yes"]:
        from subprocess import Popen

        Popen("pip install dendropy", shell=True).wait()
        try:
            import dendropy
        except ImportError:
            sys.exit("Failed to install dendropy, please see https://pythonhosted.org/DendroPy/ for further details")
    else:
        sys.exit("Aborting. Please see https://pythonhosted.org/DendroPy/ for installation details\n")

from dendropy.datamodel.treemodel import Tree, Node
from dendropy.datamodel.treecollectionmodel import TreeList
from dendropy.datamodel.taxonmodel import TaxonNamespace
from dendropy.calculate import treecompare

# BuddySuite specific
from MyFuncs import TempDir
import buddy_resources as br


# ##################################################### WISH LIST #################################################### #
def delete_metadata(_trees):
    return _trees


def decode_accessions(_phylobuddy):  # TODO: Implement decode_accessions
    return _phylobuddy


# Compare two trees, and add colour to the nodes that differ. [ ]

# Implement sum_bootstrap(), but generalize to any value.

# Prune taxa [X]

# Regex taxa names

# List all leaf names [X]

# 'Clean' a tree, as implemented in phyutility

# See http://cegg.unige.ch/system/files/nwutils_tutorial.pdf for ideas
# Re-implement many or all of Phyultility commands: https://code.google.com/p/phyutility/


# #################################################### CHANGE LOG #################################################### #
# ##################################################### GLOBALS ###################################################### #
CONFIG = br.config_values()
VERSION = br.Version("PhyloBuddy", 1, 'beta', br.contributors)
OUTPUT_FORMATS = ["newick", "nexus", "nexml"]
PHYLO_INFERENCE_TOOLS = ["raxml", "phyml", "fasttree"]


# #################################################### PHYLOBUDDY #################################################### #
class PhyloBuddy:
    def __init__(self, _input, _in_format=None, _out_format=None):
        # ####  IN AND OUT FORMATS  #### #
        # Holders for input type. Used for some error handling below

        in_handle = None
        raw_seq = None
        in_file = None
        self.trees = []
        tree_classes = [Tree]  # Dendropy Tree
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
            self.in_format = _guess_format(_input)
            self.out_format = str(self.in_format).lower() if not _out_format else str(_out_format).lower()

        else:
            self.in_format = _in_format

        if not self.in_format:
            if in_file:
                raise GuessError("Could not automatically determine the format of '{0}'.\n"
                                 "Try explicitly setting it with the -f flag.".format(in_file))
            elif raw_seq:
                raise GuessError("Could not automatically determine the format from raw input\n{0} ..."
                                 "Try explicitly setting it with the -f flag.".format(raw_seq)[:50])
            elif in_handle:
                raise GuessError("Could not automatically determine the format from input file-like object\n{0} ..."
                                 "Try explicitly setting it with the -f flag.".format(in_handle)[:50])
            else:
                raise GuessError("Unable to determine the format or input type. "
                                 "Please check how PhyloBuddy is being called.")

        self.out_format = self.in_format if not _out_format else _out_format

        # ####  RECORDS  #### #
        if type(_input) == PhyloBuddy:
            self.trees = _input.trees

        elif isinstance(_input, list):
            # make sure that the list is actually Bio.Phylo records (just test a few...)
            _sample = _input if len(_input) < 5 else sample(_input, 5)
            for _tree in _sample:
                if type(_tree) not in tree_classes:
                    raise TypeError("Tree list is not populated with Phylo objects.")
            self.trees = _input

        elif str(type(_input)) == "<class '_io.TextIOWrapper'>" or isinstance(_input, StringIO):
            tmp_dir = TempDir()
            with open("%s/tree.tmp" % tmp_dir.path, "w") as _ofile:
                _ofile.write(in_handle)

            # Removes figtree data so parser doesn't die
            figtree = _extract_figtree_metadata("%s/tree.tmp" % tmp_dir.path)
            if figtree is not None:
                with open("%s/tree.tmp" % tmp_dir.path, "w") as _ofile:
                    _ofile.write(figtree[0])

            if self.in_format != 'nexml':
                _trees = Tree.yield_from_files(files=["%s/tree.tmp" % tmp_dir.path], schema=self.in_format,
                                               extract_comment_metadata=True)
            else:
                _trees = Tree.yield_from_files(files=["%s/tree.tmp" % tmp_dir.path], schema=self.in_format)

            for _tree in _trees:
                self.trees.append(_tree)

        elif os.path.isfile(_input):
            figtree = _extract_figtree_metadata(_input)  # FigTree data being discarded here too
            if figtree is not None:
                tmp_dir = TempDir()
                with open("%s/tree.tmp" % tmp_dir.path, "w") as _ofile:
                    _ofile.write(figtree[0])
                _input = "%s/tree.tmp" % tmp_dir.path
            if self.in_format != 'nexml':
                _trees = Tree.yield_from_files(files=[_input], schema=self.in_format, extract_comment_metadata=True)
            else:
                _trees = Tree.yield_from_files(files=[_input], schema=self.in_format)
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

    def write(self, _file_path):
        with open(_file_path, "w") as _ofile:
            _ofile.write(str(self))
        return


# ################################################# HELPER FUNCTIONS ################################################# #
class GuessError(Exception):
    """Raised when input format cannot be guessed"""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


def _convert_to_ete(_tree, ignore_color=False):
    """
    Converts dendropy trees to ete trees
    :param _tree: A dendropy Tree object
    :param ignore_color: Specifies if figtree color metadata should be turned into ETE NodeStyle objects
    :return: An ETE Tree object
    """
    tmp_dir = TempDir()
    with open("%s/tree.tmp" % tmp_dir.path, "w") as _ofile:
        _ofile.write(re.sub('!color', 'pb_color', _tree.as_string(schema='newick', annotations_as_nhx=True,
                                                                  suppress_annotations=False, suppress_rooting=True)))

    ete_tree = ete3.TreeNode(newick="%s/tree.tmp" % tmp_dir.path)

    if not ignore_color:  # Converts color annotations from figtree into NodeStyle objects.
        for node in ete_tree.traverse():
            if hasattr(node, 'pb_color'):
                style = ete3.NodeStyle()
                style['fgcolor'] = node.pb_color
                style['hz_line_color'] = node.pb_color
                node.set_style(style)
    else:
        for node in ete_tree.traverse():
            node.del_feature('pb_color')

    return ete_tree


def _extract_figtree_metadata(_file_path):
    """
    Removes the figtree block from nexus files
    :param _file_path: Specifies the nexus file path
    :return: A length 2 tuple containing the nexus data and the figtree block
    """
    with open(_file_path, "r") as _tree_file:
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


def _format_to_extension(_format):
    format_to_extension = {'fasta': 'fa', 'fa': 'fa', 'genbank': 'gb', 'gb': 'gb', 'newick': 'nwk', 'nwk': 'nwk',
                           'nexus': 'nex', 'nex': 'nex', 'phylip': 'phy', 'phy': 'phy', 'phylip-relaxed': 'phyr',
                           'phyr': 'phyr', 'stockholm': 'stklm', 'stklm': 'stklm'}
    return format_to_extension[_format]


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


def _order_dendropy_features(node):
    node.annotations._item_list = sorted(node.annotations._item_list, key=lambda x: x.name)
    node.annotations._item_set = set(node.annotations._item_list)
    return node


def _make_copy(_phylobuddy):
    """
    Returns a copy of the PhyloBuddy object
    :param _phylobuddy: The PhyloBuddy object to be copied
    :return: A copy of the original PhyloBuddy object
    """
    _copy = deepcopy(_phylobuddy)
    return _copy


def _stderr(message, quiet=False):
    if not quiet:
        sys.stderr.write(message)
    return


def _stdout(message, quiet=False):
    if not quiet:
        sys.stdout.write(message)
    return


# ################################################ MAIN API FUNCTIONS ################################################ #
def consensus_tree(phylobuddy, frequency=.5):  # ToDo: Remove the length parameters from the tree
    """
    Create a consensus tree from two or more trees.
    :param phylobuddy: PhyloBuddy object
    :param frequency: The frequency threshold of a node for it to be included in the new tree.
    :return: The modified PhyloBuddy object
    """
    _trees = TreeList(phylobuddy.trees)
    _consensus = _trees.consensus(min_freq=frequency)
    phylobuddy.trees = [_consensus]
    return phylobuddy


def display_trees(phylobuddy):
    """
    Displays trees in an ETE GUI window, one-by-one.
    :param phylobuddy: PhyloBuddy object
    :return: None
    """
    if "DISPLAY" not in os.environ:
        raise SystemError("This system is not graphical, so display_trees() will not work. Try using trees_to_ascii()")

    for _tree in phylobuddy.trees:
        _convert_to_ete(_tree).show()
    return


def distance(phylobuddy, method='weighted_robinson_foulds'):
    """
    Calculates distance metrics between pairs of trees
    :param phylobuddy: PhyloBuddy object
    :param method: The tree comparison method ([un]weighted_robinson_foulds/euclidean_distance)
    :return: A dictionary of dictonaries containing the distances between tree pairs. dict[tree1][tree2]
    """
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
    keypairs = []

    for indx1, tree1 in enumerate(phylobuddy.trees):  # Compares all-by-all
        for indx2, tree2 in enumerate(phylobuddy.trees):
            if tree1 is not tree2:  # Will not compare to itself
                key1 = 'tree_{0}'.format(indx1 + 1) if tree1.label not in [None, ''] else tree1.label
                key2 = 'tree_{0}'.format(indx2 + 1) if tree2.label not in [None, ''] else tree2.label
                if key1 not in output.keys():
                    output[key1] = OrderedDict()
                if key2 not in output.keys():
                    output[key2] = OrderedDict()
                if (key2, key1) not in keypairs:  # Prevent repeated (reversed) searches
                    keypairs.append((key1, key2))
                    if method == 'wrf':
                        output[key1][key2] = treecompare.weighted_robinson_foulds_distance(tree1, tree2)
                        output[key2][key1] = output[key1][key2]
                    elif method == 'uwrf':
                        output[key1][key2] = treecompare.symmetric_difference(tree1, tree2)
                        output[key2][key1] = output[key1][key2]
                    # elif _method == 'mgk':
                    #     _output[_key1][_key2] = treecompare.mason_gamer_kellogg_score(_tree1, _tree2)
                    #     _output[_key2][_key1] = _output[_key1][_key2]
                    else:
                        output[key1][key2] = treecompare.euclidean_distance(tree1, tree2)
                        output[key2][key1] = output[key1][key2]
    return output


def generate_tree(alignbuddy, tool, params=None, keep_temp=None):
    # ToDo Check that this works for other versions of RAxML and PhyML
    """
    Calls tree building tools to generate trees
    :param alignbuddy: The AlignBuddy object containing the alignments for building the trees
    :param tool: The tree building tool to be used (raxml/phyml/fasttree)
    :param params: Additional parameters to be passed to the tree building tool
    :param keep_temp: Determines if/where the temporary files will be kept
    :return: A PhyloBuddy object containing the trees produced.
    """

    if params is None:
        params = ''
    tool = tool.lower()

    if keep_temp:  # Store files in temp dir in a non-temporary directory
        if os.path.exists(keep_temp):
            raise FileExistsError()

    if tool not in ['raxml', 'phyml', 'fasttree']:  # Supported tools
        raise AttributeError("{0} is not a valid alignment tool.".format(tool))

    if shutil.which(tool) is None:  # Tool must be callable from command line
        raise ProcessLookupError('#### Could not find {0} in $PATH. ####\nInstallation instructions '
                                 'may be found at {1}.\n'.format(tool, _get_tree_binaries(tool)))

    else:
        tmp_dir = TempDir()
        tmp_in = "{0}/tmp.del".format(tmp_dir.path)

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
        for indx, token in enumerate(params):
            if os.path.exists(token):
                params[indx] = os.path.abspath(token)
        params = ' '.join(params)

        alignbuddy.out_format = 'phylip-interleaved'  # Supported by most tree builders
        with open("{0}/tmp.del".format(tmp_dir.path), 'w') as out_file:
            out_file.write(str(alignbuddy))  # Most tree builders require an input file

        if tool == 'raxml':
            params = remove_invalid_params({'-s': True, '-n': True, '-w': True})
            if '-T' not in params:  # Num threads
                params += ' -T 2'
            if '-m' not in params:  # An evolutionary model is required
                if alignbuddy.alpha in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna, IUPAC.ambiguous_rna,
                                        IUPAC.unambiguous_rna]:
                    _stderr("Warning: Using default evolutionary model GTRCAT\n")
                    params += " -m GTRCAT"
                elif alignbuddy.alpha == IUPAC.protein:
                    _stderr("Warning: Using default evolutionary model PROTCATLG\n")
                    params += " -m PROTCATLG"

            if '-p' not in params:  # RNG seed
                params += ' -p 12345'
            if '-#' not in params and '-N' not in params:  # Number of trees to build
                params += ' -# 1'
            command = '{0} -s {1} {2} -n result -w {3}'.format(tool, tmp_in, params, tmp_dir.path)

        elif tool == 'phyml':
            params = remove_invalid_params({'-q': False, '--sequential': False, '-u': True, '--inputtree': True,
                                            '--run_id': True})
            if alignbuddy.alpha in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna, IUPAC.ambiguous_rna,
                                    IUPAC.unambiguous_rna]:
                if '-d nt' not in params and '--datatype nt' not in params:
                    params += ' -d nt'

            elif alignbuddy.alpha == IUPAC.protein:
                if '-d aa' not in params and '--datatype aa' not in params:
                    params += ' -d aa'

            command = '{0} -i {1} {2}'.format(tool, tmp_in, params)

        elif tool == 'fasttree':
            if '-n ' not in params and '--multiple' not in params and len(alignbuddy.alignments) > 1:
                params += ' -n {0}'.format(len(alignbuddy.alignments))  # Number of alignments to be input
            if alignbuddy.alpha in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna, IUPAC.ambiguous_rna,
                                    IUPAC.unambiguous_rna]:
                command = '{0} {1} -nt {2}'.format(tool, params, tmp_in)  # fasttree must be told what alphabet to use
            else:
                command = '{0} {1} {2}'.format(tool, params, tmp_in)
        else:
            raise AttributeError("'%s' is an unknown phylogenetics tool." % tool)  # Should be unreachable

        output = ''

        try:
            if tool in ['raxml', 'phyml']:  # If tool writes to file
                Popen(command, shell=True, universal_newlines=True, stdout=sys.stderr).wait()
                if not os.path.isfile('{0}/RAxML_bestTree.result'.format(tmp_dir.path)) \
                        and not os.path.isfile('{0}/tmp.del_phyml_tree.txt'.format(tmp_dir.path)):
                    raise FileNotFoundError("Error: {0} failed to generate a tree.".format(tool))
            else:  # If tool outputs to stdout
                output = check_output(command, shell=True, universal_newlines=True)

        except CalledProcessError:  # Haven't been able to find a way to get here. Needs a test.
            raise RuntimeError('\n#### {0} threw an error. Scroll up for more info. ####\n\n'.format(tool))

        if tool == 'raxml':  # Pull tree from written file
            num_runs = re.search('-[#N] ([0-9]+)', params)
            num_runs = int(num_runs.group(1))
            if num_runs > 1:
                for tree_indx in range(num_runs):
                    with open('{0}/RAxML_result.result.RUN.{1}'.format(tmp_dir.path, tree_indx)) as result:
                        output += result.read()

            else:
                with open('{0}/RAxML_bestTree.result'.format(tmp_dir.path)) as result:
                    output += result.read()

        elif tool == 'phyml':
            with open('{0}/tmp.del_phyml_tree.txt'.format(tmp_dir.path)) as result:
                output += result.read()

        if keep_temp:  # Store temp files
            shutil.copytree(tmp_dir.path, keep_temp)

        phylobuddy = PhyloBuddy(output)

        _stderr("Returning to PhyloBuddy...\n\n")

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
            if node.taxon is not None:
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
    for tree in phylobuddy.trees:
        taxa_to_prune = []
        namespace = TaxonNamespace()
        for node in tree:  # Populate the namespace for easy iteration
            if node.taxon is not None:
                namespace.add_taxon(node.taxon)
        for taxon in namespace.labels():
            for pattern in patterns:
                if re.search(pattern, taxon):  # Sets aside the names of the taxa to be pruned
                    taxa_to_prune.append(taxon)
        for taxon in taxa_to_prune:  # Removes the nodes from the tree
            tree.prune_taxa_with_labels(StringIO(taxon))


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
            if node.label is not None:
                node.label = re.sub(query, replace, node.label)
            if node.taxon is not None and node.taxon.label is not None:
                node.taxon.label = re.sub(query, replace, node.taxon.label)
    return phylobuddy


def root(phylobuddy, root_nodes=None):
    """
    Place a new root on trees
    :param phylobuddy: PhyloBuddy object
    :param root_nodes: A string with the taxon label to root on, or a list of labels and the most common ancestor
    node will be rooted on
    :return: The modified PhyloBuddy object
    """
    def _root(_tree, _root_nodes=None):
        if _root_nodes:
            _root_nodes = _root_nodes if type(_root_nodes) == list else [_root_nodes]
            mrca = None
            if len(_root_nodes) == 1:
                leaf_node = _tree.find_node_with_taxon_label(_root_nodes[0])
                if leaf_node:
                    mrca = leaf_node._parent_node

            else:
                mrca = _tree.mrca(taxon_labels=root_nodes)

            if mrca:
                _tree.reroot_at_node(mrca, update_bipartitions=True, suppress_unifurcations=False)
            else:
                _root(_tree)
        else:
            # WARNING! There is a bug in DendroPy leading to an infinite loop here. The dendropy folks have fixed it
            # in their development branch but it is not yet in the main branch
            _tree.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)

    for tree in phylobuddy.trees:
        _root(tree, root_nodes)
        tree.is_rooted = True

    return phylobuddy


def show_diff(phylobuddy):  # Doesn't work.
    sys.exit('show_diff() is not implemented yet.')
    # if len(_phylobuddy.trees) != 2:
    #    raise AssertionError("PhyloBuddy object should have exactly 2 trees.")
    trees = [_convert_to_ete(phylobuddy.trees[0], ignore_color=True),
             _convert_to_ete(phylobuddy.trees[1], ignore_color=True)]

    paths = []

    def get_all_paths(_tree, _path_list=list()):
        new_list = deepcopy(_path_list)
        new_list.append(_tree)
        print(_tree)
        if not _tree.is_leaf():
            children = _tree.child_nodes()
            for child in children:
                new_list.append(get_all_paths(child, _path_list=new_list))
        else:
            paths.append(new_list)

    get_all_paths(phylobuddy.trees[0].seed_node)
    print(paths)

    data = trees[0].robinson_foulds(trees[1])
    tree1_only = data[3] - data[4]
    tree2_only = data[4] - data[3]

    tree1_only_set = set()
    for tup in tree1_only:
        for _name in tup:
            tree1_only_set.add(_name)

    tree2_only_set = set()
    for tup in tree2_only:
        for _name in tup:
            tree2_only_set.add(_name)

    for node in trees[0]:
        if node.name in tree1_only_set:
            node.add_feature('pb_color', '#ff0000')
        else:
            node.add_feature('pb_color', '#00ff00')

    for node in trees[1]:
        if node.name in tree2_only_set:
            node.add_feature('pb_color', '#ff0000')
        else:
            node.add_feature('pb_color', '#00ff00')

    tmp_dir = TempDir()
    with open("%s/tree1.tmp" % tmp_dir.path, "w") as ofile:
        ofile.write(re.sub('pb_color', '!color', trees[0].write(features=[])))
    with open("%s/tree2.tmp" % tmp_dir.path, "w") as ofile:
        ofile.write(re.sub('pb_color', '!color', trees[1].write(features=[])))

    pb1 = PhyloBuddy(_input="%s/tree1.tmp" % tmp_dir.path, _in_format=phylobuddy.in_format,
                     _out_format=phylobuddy.out_format)
    pb2 = PhyloBuddy(_input="%s/tree2.tmp" % tmp_dir.path, _in_format=phylobuddy.in_format,
                     _out_format=phylobuddy.out_format)

    trees = pb1.trees + pb2.trees

    phylobuddy = PhyloBuddy(_input=trees, _in_format=phylobuddy.in_format, _out_format=phylobuddy.out_format)

    return phylobuddy


def show_unique(phylobuddy):
    """
    Colors all of the nodes that aren't common between two trees
    :param phylobuddy: PhyloBuddy object
    :return: The labeled PhyloBuddy object
    """
    if len(phylobuddy.trees) != 2:
        raise AssertionError("PhyloBuddy object should have exactly 2 trees.")

    trees = [_convert_to_ete(phylobuddy.trees[0], ignore_color=True),
             _convert_to_ete(phylobuddy.trees[1], ignore_color=True)]  # Need ETE so we can compare them

    try:
        data = trees[0].robinson_foulds(trees[1])

    except TreeError as _e:
        if "Unrooted tree found!" in str(_e):
            data = trees[0].robinson_foulds(trees[1], unrooted_trees=True)
        else:
            raise _e

    tree1_only = []
    tree2_only = []

    common_leaves = data[2]
    for names in data[3]:
        if len(names) == 1:
            tree1_only.append(names[0])
    for names in data[4]:
        if len(names) == 1:
            tree2_only.append(names[0])

    for node in trees[0]:  # Colors the nodes
        if node.name in common_leaves:
            node.add_feature('pb_color', '#00ff00')
        else:
            node.add_feature('pb_color', '#ff0000')
    for node in trees[1]:
        if node.name in common_leaves:
            node.add_feature('pb_color', '#00ff00')
        else:
            node.add_feature('pb_color', '#ff0000')

    tmp_dir = TempDir()  # Convert back to dendropy

    # Unnamed nodes are still given a 'name' feature, which is breaking downstream stuff
    def delete_inner_names(input_tree):
        output = re.sub('pb_color', '!color', input_tree)
        output = re.sub("name=:", "", output)
        output = re.sub("name=\]", "]", output)
        return output

    with open("%s/tree1.tmp" % tmp_dir.path, "w") as ofile:
        ofile.write(delete_inner_names(trees[0].write(features=[])))

    with open("%s/tree2.tmp" % tmp_dir.path, "w") as ofile:
        ofile.write(delete_inner_names(trees[1].write(features=[])))

    pb1 = PhyloBuddy(_input="%s/tree1.tmp" % tmp_dir.path, _in_format=phylobuddy.in_format,
                     _out_format=phylobuddy.out_format)
    pb2 = PhyloBuddy(_input="%s/tree2.tmp" % tmp_dir.path, _in_format=phylobuddy.in_format,
                     _out_format=phylobuddy.out_format)

    # The 'apply()' method traverses all nodes and executes the provided function with 'node' as the argument
    pb1.trees[0].seed_node.apply(_order_dendropy_features, _order_dendropy_features, _order_dendropy_features)
    pb2.trees[0].seed_node.apply(_order_dendropy_features, _order_dendropy_features, _order_dendropy_features)

    phylobuddy.trees = pb1.trees + pb2.trees

    return phylobuddy


def split_polytomies(phylobuddy):
    """
    Randomly splits polytomies. This function was drawn almost verbatim from the DendroPy Tree.resolve_polytomies()
    method. The main difference is that a very small edge_length is assigned to new nodes when a polytomy is split.
    :param phylobuddy: PhyloBuddy object
    :return: Modified PhyloBuddy object
    """
    rng = random.Random()
    for tree in phylobuddy.trees:
        polytomies = []
        for node in tree.postorder_node_iter():
            if len(node._child_nodes) > 2:
                polytomies.append(node)
        for node in polytomies:
            to_attach = rng.sample(node._child_nodes, len(node._child_nodes) - 2)
            for child in to_attach:
                node.remove_child(child)
            attachment_points = list(node._child_nodes)
            while len(to_attach) > 0:
                next_child = to_attach.pop()
                next_sib = rng.choice(attachment_points)
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
                               "this argument can be left blank."), br.pb_flags, br.pb_modifiers, VERSION)

    in_args = parser.parse_args()

    phylobuddy = []
    tree_set = ""

    if in_args.in_format and in_args.in_format.lower() not in OUTPUT_FORMATS:
        _stderr("Error: The format '%s' passed in with the -f flag is not recognized. "
                "Valid options include %s.\n" % (in_args.in_format, OUTPUT_FORMATS))
        sys.exit()

    if in_args.out_format and in_args.out_format.lower() not in OUTPUT_FORMATS:
        _stderr("Error: Output type %s is not recognized/supported\n" % in_args.out_format)
        sys.exit()

    if not in_args.generate_tree:  # If passing in an alignment, don't want to try and build PhyloBuddy obj
        for tree_set in in_args.trees:
            if isinstance(tree_set, TextIOWrapper) and tree_set.buffer.raw.isatty():
                _stderr("Warning: No input detected. Process will be aborted.")
                sys.exit()
            tree_set = PhyloBuddy(tree_set, in_args.in_format, in_args.out_format)
            phylobuddy += tree_set.trees
        phylobuddy = PhyloBuddy(phylobuddy, tree_set.in_format, tree_set.out_format)

    return in_args, phylobuddy


def command_line_ui(in_args, phylobuddy, skip_exit=False):
    # ############################################## INTERNAL FUNCTIONS ############################################## #
    def _print_trees(_phylobuddy):
        if in_args.test:
            _stderr("*** Test passed ***\n", in_args.quiet)

        elif in_args.in_place:
            _in_place(str(_phylobuddy), in_args.trees[0])

        else:
            _stdout("{0}\n".format(str(_phylobuddy).rstrip()))

    def _in_place(_output, _path):
        if not os.path.isfile(str(_path)):
            _stderr("Warning: The -i flag was passed in, but the positional argument doesn't seem to be a "
                    "file. Nothing was written.\n", in_args.quiet)
            _stderr("%s\n" % _output.strip(), in_args.quiet)
        else:
            with open(os.path.abspath(_path), "w") as _ofile:
                _ofile.write(_output)
            _stderr("File over-written at:\n%s\n" % os.path.abspath(_path), in_args.quiet)

    def _exit(_tool, skip=skip_exit):
        if skip:
            return
        usage = br.Usage()
        usage.increment("PhyloBuddy", VERSION.short(), _tool)
        usage.save()
        sys.exit()

    def _raise_error(_err, _tool, check_string=None):
        if check_string:
            if type(check_string) == str:
                check_string = [check_string]
            re_raise = True
            for _string in check_string:
                if _string in str(_err):
                    re_raise = False
                    break
            if re_raise:
                raise _err
        _stderr("{0}: {1}\n".format(_err.__class__.__name__, str(_err)), in_args.quiet)
        _exit(_tool)

    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Consensus tree
    if in_args.consensus_tree:
        frequency = in_args.consensus_tree[0]
        frequency = 0.5 if not frequency else frequency
        if not 0 <= frequency <= 1:
            _stderr("Warning: The frequency value should be between 0 and 1. Defaulting to 0.5.\n\n")
            frequency = 0.5

        _print_trees(consensus_tree(phylobuddy, frequency))
        _exit("consensus_tree")

    # Display trees
    if in_args.display_trees:
        try:
            display_trees(phylobuddy)
        except SystemError:
            _stderr("Error: Your system is non-graphical, so display_trees can not work. "
                    "Please use print_trees instead.")
        _exit("display_trees")

    # Distance
    if in_args.distance:
        if in_args.distance[0]:
            output = distance(phylobuddy, in_args.distance[0])
        else:
            output = distance(phylobuddy)

        _stderr('Tree 1\tTree 2\tValue\n')
        keypairs = []
        for key1 in output:
            for key2 in output[key1]:
                if (key2, key1) not in keypairs:
                    keypairs.append((key1, key2))
                    _stdout('{0}\t{1}\t{2}\n'.format(key1, key2, output[key1][key2]))
        _exit("distance")

    # Generate Tree
    if in_args.generate_tree:
        in_args.generate_tree = in_args.generate_tree[0]
        alignbuddy = []
        align_set = None
        try:
            import AlignBuddy as Alb
        except ImportError:
            _raise_error(ImportError("AlignBuddy is needed to use generate_msa(). Please re-run the installer and "
                                     "add AlignBuddy to your system."), "generate_tree")

        for align_set in in_args.trees:  # Build an AlignBuddy object
            if isinstance(align_set, TextIOWrapper) and align_set.buffer.raw.isatty():
                _stderr("Warning: No input detected. Process will be aborted.")
                sys.exit()
            align_set = Alb.AlignBuddy(align_set, in_args.in_format, in_args.out_format)
            alignbuddy += align_set.alignments
        if align_set:
            alignbuddy = Alb.AlignBuddy(alignbuddy, align_set.in_format, align_set.out_format)
        else:
            alignbuddy = Alb.AlignBuddy(alignbuddy, in_args.in_format, in_args.out_format)

        # Glean the argument order being passed in
        phylo_program = None

        for tool in PHYLO_INFERENCE_TOOLS:
            breakout = False
            while True:
                if breakout:
                    break
                breakout = True
                for arg in in_args.generate_tree:
                    if tool == arg.lower():
                        phylo_program = tool
                        del in_args.generate_tree[in_args.generate_tree.index(arg)]
                        breakout = False
                        break

        if not phylo_program:
            for tool in PHYLO_INFERENCE_TOOLS:
                if shutil.which(tool):
                    phylo_program = tool
                    break

        if not phylo_program:
            _raise_error(AttributeError("A valid phylogenetic inference program was not detected."), "generate_tree")

        params = None if not in_args.generate_tree else in_args.generate_tree[0]
        generated_trees = None
        try:
            generated_trees = generate_tree(alignbuddy, phylo_program, params, in_args.keep_temp)
        except (FileExistsError, AttributeError, ProcessLookupError, RuntimeError) as _e:
            _raise_error(_e, "generate_tree")
        except FileNotFoundError as _e:
            if str(_e) != "Error: {0} failed to generate a tree.".format(phylo_program):
                raise FileNotFoundError(_e)
            else:
                _raise_error(str(_e), "generate_tree")

        if in_args.out_format:
            generated_trees.out_format = in_args.out_format

        _print_trees(generated_trees)
        _exit("generate_tree")

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
        _stdout('%s\n\n' % output.strip())
        _exit("list_ids")

    # Print trees
    if in_args.print_trees:
        tree_table = trees_to_ascii(phylobuddy)
        output = ''
        for key in tree_table:
            output += '\n#### {0} ####\n'.format(key)
            output += tree_table[key]
        output += '\n'
        _stdout(output)
        _exit("print_trees")

    # Prune taxa
    if in_args.prune_taxa:
        prune_taxa(phylobuddy, *in_args.prune_taxa[0])
        _print_trees(phylobuddy)
        _exit("prune_taxa")

    # Rename IDs
    if in_args.rename_ids:
        _print_trees(rename(phylobuddy, in_args.rename_ids[0], in_args.rename_ids[1]))
        _exit("rename_ids")

    # Root
    if in_args.root:
        root_nodes = in_args.root[0]
        if root_nodes:
            root(phylobuddy, root_nodes)
        else:
            root(phylobuddy)
        _print_trees(phylobuddy)
        _exit("root")

    # Screw formats
    if in_args.screw_formats:
        if in_args.screw_formats not in OUTPUT_FORMATS:
            _stderr("Error: unknown format '%s'\n" % in_args.screw_formats)
        else:
            phylobuddy.out_format = in_args.screw_formats
            if in_args.in_place and os.path.isfile(str(in_args.trees[0])):
                # Need to change the file extension
                os.remove(in_args.trees[0])
                if "." in os.path.abspath(in_args.trees[0]):
                    in_args.trees[0] = ".".join(os.path.abspath(in_args.trees[0]).split(".")[:-1]) + \
                                       "." + _format_to_extension(phylobuddy.out_format)
                else:
                    in_args.trees[0] = "%s.%s" % (os.path.abspath(in_args.trees[0]),
                                                  _format_to_extension(phylobuddy.out_format))

                open(in_args.trees[0], "w").close()
            _print_trees(phylobuddy)
        _exit("screw_formats")

    # Show unique
    if in_args.show_unique:
        try:
            _print_trees(show_unique(phylobuddy))
        except AssertionError as _e:
            if str(_e) == "PhyloBuddy object should have exactly 2 trees.":
                _raise_error(_e, "show_unique")
            else:
                raise _e

        _exit("show_unique")

    # Split polytomies
    if in_args.split_polytomies:
        split_polytomies(phylobuddy)
        _print_trees(phylobuddy)
        _exit("split_polytomies")

    # Unroot
    if in_args.unroot:
        unroot(phylobuddy)
        _print_trees(phylobuddy)
        _exit("unroot")

if __name__ == '__main__':
    try:
        command_line_ui(*argparse_init())
    except (KeyboardInterrupt, GuessError) as e:
        print(e)
    except SystemExit:
        pass
    except Exception as e:
        br.send_traceback("PhyloBuddy", e)
