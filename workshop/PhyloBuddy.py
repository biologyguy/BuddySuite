#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, version 2 of the License (GPLv2).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details at http://www.gnu.org/licenses/.

name: PhyloBuddy.py
date: Dec-6-2014
version: 1, unstable
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: Gnu General Public License, Version 2.0 (http://www.gnu.org/licenses/gpl.html)
derivative work: No

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
except ImportError:
    confirm = input("PhyloBuddy requires ETE v3+, which was not detected on your system. Try to install [y]/n? ")
    if confirm.lower() in ["", "y", "yes"]:
        Popen("pip install --upgrade  https://github.com/jhcepas/ete/archive/3.0.zip", shell=True).wait()
        try:
            import ete3
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


from dendropy.datamodel.treemodel import Tree
from dendropy.datamodel.treecollectionmodel import TreeList
from dendropy.datamodel.taxonmodel import TaxonNamespace
from dendropy.calculate import treecompare

# BuddySuite specific
from MyFuncs import TempDir
import buddy_resources as br


# ##################################################### WISH LIST #################################################### #
def unroot(_trees):
    return _trees


def screw_formats(_phylobuddy, _format):
    return _phylobuddy, _format


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
VERSION = br.Version("PhyloBuddy", 1, 'alpha', br.contributors)


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
            stored_input = _input

        elif isinstance(_input, list):
            # make sure that the list is actually Bio.Phylo records (just test a few...)
            _sample = _input if len(_input) < 5 else sample(_input, 5)
            for _tree in _sample:
                if type(_tree) not in tree_classes:
                    raise TypeError("Tree list is not populated with Phylo objects.")
            self.trees = _input
            stored_input = deepcopy(_input)

        elif str(type(_input)) == "<class '_io.TextIOWrapper'>" or isinstance(_input, StringIO):
            tmp_dir = TempDir()
            stored_input = deepcopy(in_handle)
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
            stored_input = deepcopy(_input)
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

        if self.out_format in ["newick", "nexus", "nexml"]:
            for _tree in self.trees:
                tree_list.append(_tree)

        else:
            raise TypeError("Error: Unsupported output format.")

        if self.out_format != 'nexml':
            _output = tree_list.as_string(schema=self.out_format, annotations_as_nhx=False, suppress_annotations=False)
        else:
            _output = tree_list.as_string(schema=self.out_format)

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
    format_to_extension = {'fasta': 'fa', 'fa': 'fa', 'genbank': 'gb', 'gb': 'gb', 'nexus': 'nex',
                           'nex': 'nex', 'phylip': 'phy', 'phy': 'phy', 'phylip-relaxed': 'phyr', 'phyr': 'phyr',
                           'stockholm': 'stklm', 'stklm': 'stklm'}
    return format_to_extension[_format]


def _get_tree_binaries(_tool):
    """
    Returns a URL where a tool's binaries can be found.
    :param _tool: Specify the tree building tool to be used.
    :return: A string containing the tool's URL.
    """

    tool_dict = {'raxml': 'http://sco.h-its.org/exelixis/web/software/raxml/index.html',
                 'phyml': 'http://www.atgc-montpellier.fr/phyml/versions.php',
                 'fastree': 'http://www.microbesonline.org/fasttree/#Install'}
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


def _make_copy(_phylobuddy):
    """
    Returns a copy of the PhyloBuddy object
    :param _phylobuddy: The PhyloBuddy object to be copied
    :return: A copy of the original PhyloBuddy object
    """
    try:
        _copy = deepcopy(_phylobuddy)
    except AttributeError:
        _stderr("Warning: Deepcopy failed. Attempting workaround. Some metadata may be lost.")
        _copy = PhyloBuddy(str(_phylobuddy), _in_format=_phylobuddy.out_format, _out_format=_phylobuddy.out_format)
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
def calculate_distance(_phylobuddy, _method='weighted_robinson_foulds'):
    """
    Calculates tree distance with various algorithms
    :param _phylobuddy: The PhyloBuddy object containing the trees to be compared
    :param _method: The tree comparison method ([un]weighted_robinson_foulds/euclidean_distance)
    :return: A dictionary of dictonaries containing the distances between tree pairs. dict[tree1][tree2]
    """
    _method = _method.lower()
    if _method in ['wrf', 'weighted_robinson_foulds']:
        _method = 'wrf'
    elif _method in ['uwrf', 'sym', 'unweighted_robinson_foulds', 'symmetric']:
        _method = 'uwrf'
    # elif _method in ['mgk', 'mason_gamer_kellogg']:
    #     _method = 'mgk'
    elif _method in ['ed', 'euclid', 'euclidean']:
        _method = 'euclid'
    else:
        raise AttributeError('{0} is an invalid comparison method.'.format(_method))

    _output = OrderedDict()
    _keypairs = []

    for indx1, _tree1 in enumerate(_phylobuddy.trees):  # Compares all-by-all
        for indx2, _tree2 in enumerate(_phylobuddy.trees):
            if _tree1 is not _tree2:  # Will not compare to itself
                _key1 = 'tree_{0}'.format(indx1 + 1) if _tree1.label not in [None, ''] else _tree1.label
                _key2 = 'tree_{0}'.format(indx2 + 1) if _tree2.label not in [None, ''] else _tree2.label
                if _key1 not in _output.keys():
                    _output[_key1] = OrderedDict()
                if _key2 not in _output.keys():
                    _output[_key2] = OrderedDict()
                if (_key2, _key1) not in _keypairs:  # Prevent repeated (reversed) searches
                    _keypairs.append((_key1, _key2))
                    if _method == 'wrf':
                        _output[_key1][_key2] = treecompare.weighted_robinson_foulds_distance(_tree1, _tree2)
                        _output[_key2][_key1] = _output[_key1][_key2]
                    elif _method == 'uwrf':
                        _output[_key1][_key2] = treecompare.symmetric_difference(_tree1, _tree2)
                        _output[_key2][_key1] = _output[_key1][_key2]
                    # elif _method == 'mgk':
                    #     _output[_key1][_key2] = treecompare.mason_gamer_kellogg_score(_tree1, _tree2)
                    #     _output[_key2][_key1] = _output[_key1][_key2]
                    else:
                        _output[_key1][_key2] = treecompare.euclidean_distance(_tree1, _tree2)
                        _output[_key2][_key1] = _output[_key1][_key2]
    return _output


def consensus_tree(_phylobuddy, _frequency=.5):
    """
    Generates a consensus tree based on all the trees in phylobuddy
    :param _phylobuddy: The PhyloBuddy object to be modified
    :param _frequency: The frequency threshold of a taxa for it to be included
    :return: The modified PhyloBuddy object
    """
    _trees = TreeList(_phylobuddy.trees)
    _consensus = _trees.consensus(_frequency=_frequency)
    _phylobuddy.trees = [_consensus]
    return _phylobuddy


def display_trees(_phylobuddy):
    """
    Displays trees in an ETE GUI window, one-by-one.
    :param _phylobuddy: The PhyloBuddy object whose trees will be displayed.
    :return:
    """
    for _tree in _phylobuddy.trees:
        _convert_to_ete(_tree).show()


def generate_tree(_alignbuddy, _tool, _params=None, _keep_temp=None):
    """
    Calls tree building tools to generate trees
    :param _alignbuddy: The AlignBuddy object containing the alignments for building the trees
    :param _tool: The tree building tool to be used (raxml/phyml/fasttree)
    :param _params: Additional parameters to be passed to the tree building tool
    :param _keep_temp: Determines if/where the temporary files will be kept
    :return: A PhyloBuddy object containing the trees produced.
    """
    # NOTE: FastTree segfaults with protein alignments for an unknown reason (may be OSX only?)

    if _params is None:
        _params = ''
    _tool = _tool.lower()

    if _keep_temp:  # Store files in temp dir in a non-temporary directory
        if os.path.exists(_keep_temp):
            _stderr("Warning: {0} already exists. Please specify a different path.\n".format(_keep_temp))
            sys.exit()

    if _tool not in ['raxml', 'phyml', 'fasttree']:  # Supported tools
        raise AttributeError("{0} is not a valid alignment tool.".format(_tool))
    if shutil.which(_tool) is None:  # Tool must be callable from command line
        _stderr('#### Could not find {0} in $PATH. ####\n'.format(_tool))
        _stderr('Please go to {0} to install {1}.\n'.format(_get_tree_binaries(_tool), _tool))
        sys.exit()
    else:
        tmp_dir = TempDir()
        tmp_in = "{0}/tmp.del".format(tmp_dir.path)

        def remove_invalid_params(_dict):  # Helper method for blacklisting flags
            parameters = _params
            for _key in _dict:
                if _dict[_key] is True:  # Flag has an argument
                    _pattern = "{0} [^-]*".format(_key)
                    # Deletes until it reaches another flag.
                    # May cause issues if flag is right before the file path.
                else:
                    _pattern = "{0}".format(_key)
                parameters = re.sub(_pattern, '', parameters)
            return parameters

        _params = re.split(' ', _params, )  # Expands paths in _params to absolute paths
        for _indx, _token in enumerate(_params):
            if os.path.exists(_token):
                _params[_indx] = os.path.abspath(_token)
        _params = ' '.join(_params)

        _alignbuddy.out_format = 'phylip-interleaved'  # Supported by most tree builders
        with open("{0}/tmp.del".format(tmp_dir.path), 'w') as out_file:
            out_file.write(str(_alignbuddy))  # Most tree builders require an input file
        if _tool == 'raxml':
            _params = remove_invalid_params({'-s': True, '-n': True, '-w': True})
            if '-T' not in _params:  # Num threads
                _params += ' -T 2'
            if '-m' not in _params:  # Tree building model (REQUIRED)
                _stderr("No tree-building method specified! Use the -m flag!\n")
                sys.exit()
            if '-p' not in _params:  # RNG seed
                _params += ' -p 12345'
            if '-#' not in _params and '-N' not in _params:  # Number of trees to build
                _params += ' -# 1'
            command = '{0} -s {1} {2} -n result -w {3}'.format(_tool, tmp_in, _params, tmp_dir.path)
        elif _tool == 'phyml':
            _params = remove_invalid_params({'-q': False, '--sequential': False, '-u': True, '--inputtree': True,
                                             '--run_id': True})
            if '-m' not in _params and '--model' not in _params:  # Tree building model (REQUIRED)
                _stderr("No tree-building method specified! Use the -m flag!\n")
                sys.exit()
            if _alignbuddy.alpha in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna, IUPAC.ambiguous_rna,
                                     IUPAC.unambiguous_rna] and ('-d nt' not in _params or
                                                                 '--datatype nt' not in _params):
                _params += ' -d nt'  # phyml needs to be told which alphabet to use
            elif _alignbuddy.alpha == IUPAC.protein and ('-d aa' not in _params or '--datatype aa' not in _params):
                _params += ' -d aa'
            command = '{0} -i {1} {2}'.format(_tool, tmp_in, _params)
        elif _tool == 'fasttree':
            if '-n ' not in _params and '--multiple' not in _params and len(_alignbuddy.alignments) > 1:
                _params += ' -n {0}'.format(len(_alignbuddy.alignments))  # Number of alignments to be input
            if _alignbuddy.alpha in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna, IUPAC.ambiguous_rna,
                                     IUPAC.unambiguous_rna]:
                command = '{0} {1} -nt {2}'.format(_tool, _params, tmp_in)  # fasttree must be told what alphabet to use
            else:
                command = '{0} {1} {2}'.format(_tool, _params, tmp_in)
        else:
            command = '{0} {1} {2}'.format(_tool, _params, tmp_in)

        _output = ''

        try:
            if _tool in ['raxml', 'phyml']:  # If tool writes to file
                Popen(command, shell=True, universal_newlines=True, stdout=sys.stderr).wait()
            else:  # If tool outputs to stdout
                _output = check_output(command, shell=True, universal_newlines=True)
        except CalledProcessError:
            _stderr('\n#### {0} threw an error. Scroll up for more info. ####\n\n'.format(_tool))
            sys.exit()

        if _tool == 'raxml':  # Pull tree from written file
            num_runs = re.search('-#', _params)
            if _params[num_runs.end() + 1].isdigit():
                num_runs = int(_params[num_runs.end() + 1])
            else:
                num_runs = int(_params[num_runs.end() + 2])
            if num_runs > 1:
                for tree_indx in range(num_runs):
                    with open('{0}/RAxML_result.result.RUN.{1}'.format(tmp_dir.path, tree_indx)) as result:
                        _output += result.read()
            else:
                with open('{0}/RAxML_bestTree.result'.format(tmp_dir.path)) as result:
                    _output += result.read()
        elif _tool == 'phyml':
            with open('{0}/tmp.del_phyml_tree.txt'.format(tmp_dir.path)) as result:
                _output += result.read()

        if _keep_temp:  # Store temp files
            try:
                shutil.copytree(tmp_dir.path, _keep_temp)
            except FileExistsError:
                # Should never get here
                pass

        _phylobuddy = PhyloBuddy(_output)

        _stderr("Returning to PhyloBuddy...\n\n")

        return _phylobuddy


def list_ids(_phylobuddy):
    """
    Returns a dictionary of tree names and node labels
    :param _phylobuddy: The PhyloBuddy object to be analyzed
    :return: A dictionary of tree names and node labels
    """
    _output = OrderedDict()
    for indx, _tree in enumerate(_phylobuddy.trees):
        _namespace = TaxonNamespace()
        for node in _tree:
            if node.taxon is not None:
                _namespace.add_taxon(node.taxon)
        _key = _tree.label if _tree.label not in [None, ''] else 'tree_{0}'.format(str(indx + 1))
        _output[_key] = list(_namespace.labels())
    return _output


def prune_taxa(_phylobuddy, *_patterns):
    """
    Prunes taxa that match one or more regex patterns
    :param _phylobuddy: The PhyloBuddy object whose trees will be pruned.
    :param _patterns: One or more regex patterns.
    :return: The same PhyloBuddy object after pruning.
    """
    for _tree in _phylobuddy.trees:
        taxa_to_prune = []
        _namespace = TaxonNamespace()
        for node in _tree:  # Populate the namespace for easy iteration
            if node.taxon is not None:
                _namespace.add_taxon(node.taxon)
        for _taxon in _namespace.labels():
            for _pattern in _patterns:
                if re.search(_pattern, _taxon):  # Sets aside the names of the taxa to be pruned
                    taxa_to_prune.append(_taxon)
        for _taxon in taxa_to_prune:  # Removes the nodes from the tree
            _tree.prune_taxa_with_labels(StringIO(_taxon))


def rename(_phylobuddy, _query, _replace):
    """
    Substitutes matches in node names with a string
    :param _phylobuddy: The PhyloBuddy object to be modified
    :param _query: The regex pattern to be searched
    :param _replace: The string to replace the matches with
    :return: The modified PhyloBuddy object
    """
    for indx, _tree in enumerate(_phylobuddy.trees):
        for node in _tree:
            if node.label is not None:
                node.label = re.sub(_query, _replace, node.label)
            if node.taxon is not None and node.taxon.label is not None:
                node.taxon.label = re.sub(_query, _replace, node.taxon.label)
    return _phylobuddy


def show_diff(_phylobuddy):  # Doesn't work.
    sys.exit('show_diff() is not implemented yet.')
    # if len(_phylobuddy.trees) != 2:
    #    raise AssertionError("PhyloBuddy object should have exactly 2 trees.")
    _trees = [_convert_to_ete(_phylobuddy.trees[0], ignore_color=True),
              _convert_to_ete(_phylobuddy.trees[1], ignore_color=True)]

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
    get_all_paths(_phylobuddy.trees[0].seed_node)
    print(paths)

    data = _trees[0].robinson_foulds(_trees[1])
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

    for _node in _trees[0]:
        if _node.name in tree1_only_set:
            _node.add_feature('pb_color', '#ff0000')
        else:
            _node.add_feature('pb_color', '#00ff00')

    for _node in _trees[1]:
        if _node.name in tree2_only_set:
            _node.add_feature('pb_color', '#ff0000')
        else:
            _node.add_feature('pb_color', '#00ff00')

    tmp_dir = TempDir()
    with open("%s/tree1.tmp" % tmp_dir.path, "w") as _ofile:
        _ofile.write(re.sub('pb_color', '!color', _trees[0].write(features=[])))
    with open("%s/tree2.tmp" % tmp_dir.path, "w") as _ofile:
        _ofile.write(re.sub('pb_color', '!color', _trees[1].write(features=[])))

    pb1 = PhyloBuddy(_input="%s/tree1.tmp" % tmp_dir.path, _in_format=_phylobuddy.in_format,
                     _out_format=_phylobuddy.out_format)
    pb2 = PhyloBuddy(_input="%s/tree2.tmp" % tmp_dir.path, _in_format=_phylobuddy.in_format,
                     _out_format=_phylobuddy.out_format)

    _trees = pb1.trees + pb2.trees

    _phylobuddy = PhyloBuddy(_input=_trees, _in_format=_phylobuddy.in_format, _out_format=_phylobuddy.out_format)

    return _phylobuddy


def show_unique_nodes(_phylobuddy):
    """
    Colors all of the nodes that aren't common between two trees
    :param _phylobuddy: The PhyloBuddy object to be labeled
    :return: The labeled PhyloBuddy object
    """
    if len(_phylobuddy.trees) != 2:
        raise AssertionError("PhyloBuddy object should have exactly 2 trees.")

    _trees = [_convert_to_ete(_phylobuddy.trees[0], ignore_color=True),
              _convert_to_ete(_phylobuddy.trees[1], ignore_color=True)]  # Need ETE so we can compare them

    data = _trees[0].robinson_foulds(_trees[1])
    tree1_only = []
    tree2_only = []

    common_leaves = data[2]
    for _names in data[3]:
        if len(_names) == 1:
            tree1_only.append(_names[0])
    for _names in data[4]:
        if len(_names) == 1:
            tree2_only.append(_names[0])

    for _node in _trees[0]:  # Colors the nodes
        if _node.name in common_leaves:
            _node.add_feature('pb_color', '#00ff00')
        else:
            _node.add_feature('pb_color', '#ff0000')
    for _node in _trees[1]:
        if _node.name in common_leaves:
            _node.add_feature('pb_color', '#00ff00')
        else:
            _node.add_feature('pb_color', '#ff0000')

    tmp_dir = TempDir()  # Convert back to dendropy
    with open("%s/tree1.tmp" % tmp_dir.path, "w") as _ofile:
        _ofile.write(re.sub('pb_color', '!color', _trees[0].write(features=[])))
    with open("%s/tree2.tmp" % tmp_dir.path, "w") as _ofile:
        _ofile.write(re.sub('pb_color', '!color', _trees[1].write(features=[])))

    pb1 = PhyloBuddy(_input="%s/tree1.tmp" % tmp_dir.path, _in_format=_phylobuddy.in_format,
                     _out_format=_phylobuddy.out_format)
    pb2 = PhyloBuddy(_input="%s/tree2.tmp" % tmp_dir.path, _in_format=_phylobuddy.in_format,
                     _out_format=_phylobuddy.out_format)

    _trees = pb1.trees + pb2.trees

    _phylobuddy = PhyloBuddy(_input=_trees, _in_format=_phylobuddy.in_format, _out_format=_phylobuddy.out_format)

    return _phylobuddy


def split_polytomies(_phylobuddy):
    """
    Randomly splits polytomies.
    :param _phylobuddy: The PhyloBuddy object whose trees will be processed.
    :return: The same PhyloBuddy object after processing.
    """
    for _tree in _phylobuddy.trees:
        _tree.resolve_polytomies(rng=random.Random())
    return _phylobuddy


def trees_to_ascii(_phylobuddy):
    """
    Returns an ascii representation of the tree. Scales to terminal window size.
    :param _phylobuddy: The PhyloBuddy object whose trees will be converted.
    :return: A string containing ASCII representations of the trees.
    """
    _output = OrderedDict()
    for _indx, _tree in enumerate(_phylobuddy.trees):
        _key = 'tree_{0}'.format(_indx + 1) if _tree.label in [None, ''] else _tree.label
        _output[_key] = _tree.as_ascii_plot()
    return _output


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

    br.flags(parser, ("trees", "Supply file path(s) or raw tree string, If piping trees into PhyloBuddy "
                               "this argument can be left blank."), br.pb_flags, br.pb_modifiers, VERSION)

    in_args = parser.parse_args()

    phylobuddy = []
    tree_set = ""

    if not in_args.generate_tree:  # If passing in an alignment, don't want to try and build PhyloBuddy obj
        for tree_set in in_args.trees:
            if isinstance(tree_set, TextIOWrapper) and tree_set.buffer.raw.isatty():
                sys.exit("Warning: No input detected. Process will be aborted.")
            tree_set = PhyloBuddy(tree_set, in_args.in_format, in_args.out_format)
            phylobuddy += tree_set.trees
        phylobuddy = PhyloBuddy(phylobuddy, tree_set.in_format, tree_set.out_format)

    return in_args, phylobuddy


def command_line_ui(in_args, phylobuddy, skip_exit=False):
    # ############################################## INTERNAL FUNCTIONS ############################################## #
    def _print_trees(_phylobuddy):
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

    def _exit(tool, skip=skip_exit):
        if skip:
            return
        usage = br.Usage()
        usage.increment("PhyloBuddy", VERSION.short(), tool)
        usage.save()
        sys.exit()

    def _raise_error(_err, tool):
        _stderr("{0}: {1}\n".format(_err.__class__.__name__, str(_err)))
        _exit(tool)

# ############################################## COMMAND LINE LOGIC ############################################## #
    # Calculate distance
    if in_args.calculate_distance:
        output = calculate_distance(phylobuddy, in_args.calculate_distance)
        _stderr('Tree 1\tTree 2\tValue\n')
        keypairs = []
        for key1 in output:
            for key2 in output[key1]:
                if (key2, key1) not in keypairs:
                    keypairs.append((key1, key2))
                    _stdout('{0}\t{1}\t{2}\n'.format(key1, key2, output[key1][key2]))
        _exit("calculate_distance")

    # Consensus tree
    if in_args.consensus_tree:
        _print_trees(consensus_tree(phylobuddy, in_args.consensus_tree))
        _exit("consensus_tree")

    # Display trees
    if in_args.display_trees:
        display_trees(phylobuddy)
        _exit("display_trees")

    # Generate Tree
    if in_args.generate_tree:
        alignbuddy = []
        align_set = None
        try:
            import AlignBuddy as Alb
        except ImportError:
            _raise_error(ImportError("AlignBuddy is needed to use generate_msa(). Please re-run the installer and "
                                     "add AlignBuddy to your system."), "generate_tree")
            
        for align_set in in_args.trees:  # Build an AlignBuddy object
            if isinstance(align_set, TextIOWrapper) and align_set.buffer.raw.isatty():
                sys.exit("Warning: No input detected. Process will be aborted.")
            align_set = Alb.AlignBuddy(align_set, in_args.in_format, in_args.out_format)
            alignbuddy += align_set.alignments
        if align_set:
            alignbuddy = Alb.AlignBuddy(alignbuddy, align_set.in_format, align_set.out_format)
        else:
            alignbuddy = Alb.AlignBuddy(alignbuddy, in_args.in_format, in_args.out_format)

        params = in_args.params if in_args.params is None else in_args.params[0]
        generated_trees = generate_tree(alignbuddy, in_args.generate_tree[0], params, in_args.keep_temp)
        if in_args.out_format:
            generated_trees.out_format = in_args.out_format
        _stdout(str(generated_trees))
        _exit("generate_tree")

    # List ids
    if in_args.list_ids:
        listed_ids = list_ids(phylobuddy)
        columns = 1 if not in_args.list_ids[0] or in_args.list_ids[0] <= 0 else abs(in_args.list_ids[0])
        output = ""
        for key in listed_ids:
            count = 1
            output += '#### {0} ####\n'.format(key)
            if len(listed_ids[key]) == 0:
                output += 'None\n'
            else:
                for identifier in listed_ids[key]:
                    if count < columns:
                        output += "%s\t" % identifier
                        count += 1
                    else:
                        output += "%s\n" % identifier
                        count = 1
            output += '\n'
        _stdout(output)
        _exit("list_ids")

    # Prune taxa
    if in_args.prune_taxa:
        prune_taxa(phylobuddy, *in_args.prune_taxa[0])
        _print_trees(phylobuddy)
        _exit("prune_taxa")

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

    # Rename IDs
    if in_args.rename_ids:
        _print_trees(rename(phylobuddy, in_args.rename_ids[0], in_args.rename_ids[1]))
        _exit("rename_ids")

    # Split polytomies
    if in_args.split_polytomies:
        split_polytomies(phylobuddy)
        _print_trees(phylobuddy)
        _exit("split_polytomies")


if __name__ == '__main__':
    try:
        command_line_ui(*argparse_init())
    except (KeyboardInterrupt, GuessError) as e:
        print(e)
    except SystemExit:
        pass
    except Exception as e:
        br.send_traceback("PhyloBuddy", e)
