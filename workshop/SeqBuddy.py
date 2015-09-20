#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, version 2 of the License (GPLv2).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details at http://www.gnu.org/licenses/.

name: SeqBuddy.py
date: Nov-20-2014
version: 2, unstable
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: Gnu General Public License, Version 2.0 (http://www.gnu.org/licenses/gpl.html)
derivative work: No

Description:
Collection of functions that do fun stuff with sequences. Pull them into a script, or run as a command line tool.
"""

# Standard library imports
# from pprint import pprint
# import pdb
# import time
import sys
import os
import re
import string
import zipfile
import shutil
from urllib import request, error
from copy import copy, deepcopy
from random import sample, choice, randint, random
from math import floor, ceil, log
from tempfile import TemporaryDirectory
from subprocess import Popen, PIPE
from shutil import which
from hashlib import md5
from io import StringIO, TextIOWrapper
from collections import OrderedDict
from xml.sax import SAXParseException

# Third party package imports
sys.path.insert(0, "./")  # For stand alone executable, where dependencies are packaged with BuddySuite
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import *
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError
from Bio.Data import CodonTable

# My functions
from MyFuncs import run_multicore_function
import buddy_resources as br


# ##################################################### WISH LIST #################################################### #
def sim_ident(matrix):  # Return the pairwise similarity and identity scores among sequences
    x = matrix
    return x


def predict_orfs():
    return


def delete_pattern():
    # remove residues that match a given pattern from all records
    return


def mutate():
    # Apply some model of evolution to generate new sequences from input
    return


def random_aa(_length, number, matrix):
    # create random prot sequences. Not sure exactly how to implement this yet, because it would theoretically not start
    # from input sequences...
    x = [_length, number, matrix]
    return x


def random_dna(_length, number, matrix):
    # create random DNA sequences.
    x = [_length, number, matrix]
    return x


def divergence_value():
    # http://bioinformatics.org/sms/uneven.html
    return


def degenerate_dna():
    # https://github.com/carlosp420/degenerate-dna
    # Degen: http://www.phylotools.com/ptdegenoverview.htm
    return

# - Allow batch calls. E.g., if 6 files are fed in as input, run the SeqBuddy command provided independently on each
# - Add support for selecting individual sequences to modify
# - Add FASTQ support... More generally, support letter annotation mods
# - Add Clustal support
# - Get BuddySuite into PyPi
# - Check on memory requirements before execution
# - Execution timer, for long running jobs
# - Handle all stderr output from private function (to allow quiet execution)
# - Sort out a good way to manage 'lazy' imports
# - Fix --clean_seq for .phy, .phyr, .stklm, .nex

# ################################################# CHANGE LOG for V2 ################################################ #
# - New flag -t/--test, which runs a function but suppresses all stdout (stderr still returned)
# - Standard-in is handled as input now, allowing SeqBuddy to be daisy chained with pipes (|)
# - Remove the the -p flag dependencies for -prr, -li, -btr, -asl, -cts, -hsi, -frp, -drp, -ofa, -ofp, and -oi
# - Add print method to SeqBuddy class that outputs all sequences to string
# - Add molecular_weight() function that calculates molecular weight
# - Add isoelectric_point() function that calculates isoelectric point
# - Add find_CpG() function to predict CpG islands
# - Add find_pattern() function to search sequences for specific pattern
# - Add find_restriction_sites() function to find restriction sites
# - Add split_file() function to separate all seq records into their own SeqBuddy object
# - Add split_by_taxa() function. Writes individual files for groups of sequences based on an identifier in their ids
# - Unit tests
# - New graphical installer
# - Rework argparse output


# ###################################################### GLOBALS ##################################################### #
VERSION = br.Version("SeqBuddy", 2, 'alpha', br.contributors)
FORMATS = ["ids", "accessions", "summary", "full-summary", "clustal", "embl", "fasta", "fastq", "fastq-sanger",
           "fastq-solexa", "fastq-illumina", "genbank", "gb", "imgt", "nexus", "phd", "phylip", "seqxml", "sff",
           "stockholm", "tab", "qual"]


# ##################################################### SEQBUDDY ##################################################### #
class SeqBuddy:  # Open a file or read a handle and parse, or convert raw into a Seq object
    def __init__(self, _input, _in_format=None, _out_format=None, _alpha=None):
        # ####  IN AND OUT FORMATS  #### #
        # Holders for input type. Used for some error handling below
        in_handle = None
        _raw_seq = None
        in_file = None
        self.alpha = _alpha

        # Handles
        if str(type(_input)) == "<class '_io.TextIOWrapper'>":
            if not _input.seekable():  # Deal with input streams (e.g., stdout pipes)
                _temp = StringIO(_input.read())
                _input = _temp
            _input.seek(0)
            in_handle = _input.read()
            _input.seek(0)

        # Raw sequences
        if _in_format == "raw":
            _in_format = "fasta"
            _out_format = "fasta"
            if type(_input) == str:
                _input = [SeqRecord(Seq(_input), id="raw_input", description="")]
            else:
                _input = [SeqRecord(Seq(_input.read()), id="raw_input", description="")]

        # Plain text in a specific format
        if type(_input) == str and not os.path.isfile(_input):
            _raw_seq = _input
            _temp = StringIO(_input)
            _input = _temp
            _input.seek(0)

        # File paths
        try:
            if os.path.isfile(_input):
                in_file = _input
        except TypeError:  # This happens when testing something other than a string.
            pass

        if not _in_format:
            self.in_format = _guess_format(_input)
            self.out_format = str(self.in_format) if not _out_format else _out_format

        else:
            self.in_format = _in_format

        if not self.in_format:
            if in_file:
                raise GuessError("Could not determine format from _input file '{0}'.\n"
                                 "Try explicitly setting with -f flag.".format(in_file))
            elif _raw_seq:
                raise GuessError("File not found, or could not determine format from raw input\n{0} ..."
                                 "Try explicitly setting with -f flag.".format(_raw_seq)[:60])
            elif in_handle:
                raise GuessError("Could not determine format from input file-like object\n{0} ..."
                                 "Try explicitly setting with -f flag.".format(in_handle)[:50])
            else:
                raise GuessError("Unable to determine format or input type. Please check how SeqBuddy is being called.")

        self.out_format = self.in_format if not _out_format else _out_format

        # ####  RECORDS  #### #
        if type(_input) == SeqBuddy:
            _sequences = _input.records

        elif isinstance(_input, list):
            # make sure that the list is actually SeqIO records (just test a few...)
            _sample = _input if len(_input) < 5 else sample(_input, 5)
            for _seq in _sample:
                if type(_seq) != SeqRecord:
                    raise TypeError("Seqlist is not populated with SeqRecords.")
            _sequences = _input

        elif str(type(_input)) == "<class '_io.TextIOWrapper'>" or isinstance(_input, StringIO):
            _sequences = list(SeqIO.parse(_input, self.in_format))

        elif os.path.isfile(_input):
            with open(_input, "r") as _input:
                _sequences = list(SeqIO.parse(_input, self.in_format))
        else:
            _sequences = [SeqRecord(Seq(_input))]  # may be unreachable?

        if self.alpha is None:
            self.alpha = _guess_alphabet(_sequences)
        elif self.alpha in ['protein', 'prot', 'p', 'pep', IUPAC.protein]:
            self.alpha = IUPAC.protein
        elif self.alpha in ['dna', 'd', 'cds', IUPAC.ambiguous_dna]:
            self.alpha = IUPAC.ambiguous_dna
        elif self.alpha in ['rna', 'r', IUPAC.ambiguous_rna]:
            self.alpha = IUPAC.ambiguous_rna
        else:
            _stderr("WARNING: Alphabet not recognized. Correct alphabet will be guessed.\n")
            self.alpha = _guess_alphabet(_sequences)

        for _seq in _sequences:
            _seq.seq.alphabet = self.alpha

        self.records = _sequences

    def to_dict(self):
        _unique, _rep_ids, _rep_seqs, output_str = find_repeats(self)
        if len(_rep_ids) > 0:
            raise RuntimeError("There are repeat IDs in self.records\n%s" % _rep_ids)

        records_dict = {}
        for _rec in self.records:
            records_dict[_rec.id] = _rec
        return records_dict

    def print(self):
        print(self)
        return

    def __str__(self):
        if len(self.records) == 0:
            return "Error: No sequences in object.\n"

        # There is a weird bug in genbank write() that concatenates dots to the organism name (if set).
        # The following is a work around...
        if self.out_format in ["gb", "genbank"]:
            for _rec in self.records:
                try:
                    if re.search("(\. )+", _rec.annotations['organism']):
                        _rec.annotations['organism'] = "."
                except KeyError:
                    pass

        if self.out_format == "phylipi":
            _output = _phylipi(self)

        elif self.out_format == "phylipis":
            _output = _phylipi(self, "strict")

        else:
            tmp_dir = TemporaryDirectory()
            with open("%s/seqs.tmp" % tmp_dir.name, "w") as _ofile:
                SeqIO.write(self.records, _ofile, self.out_format)

            with open("%s/seqs.tmp" % tmp_dir.name, "r") as ifile:
                _output = ifile.read()

        return _output

    def write(self, _file_path):
        with open(_file_path, "w") as _ofile:
            _ofile.write(str(self))
        return


# ################################################# HELPER FUNCTIONS ################################################# #
class GuessError(Exception):
    """Raised when input format cannot be guessed"""
    def __init__(self, _value):
        self.value = _value

    def __str__(self):
        return self.value


def _download_blast_binaries(_blastn=True, _blastp=True, _blastdcmd=True):
    if os.path.exists('.buddysuite'):
        current_path = '.buddysuite/'
    else:
        current_path = os.getcwd()

    binary_source = 'https://raw.github.com/biologyguy/BuddySuite/master/workshop/build_dir/blast_binaries/'
    bins_to_dl = []
    current_os = sys.platform
    if current_os.startswith('darwin'):
        if _blastdcmd:
            bins_to_dl.append('Darwin_blastdbcmd.zip')
        if _blastn:
            bins_to_dl.append('Darwin_blastn.zip')
        if _blastp:
            bins_to_dl.append('Darwin_blastp.zip')
    elif current_os.startswith('linux'):
        if _blastdcmd:
            bins_to_dl.append('Linux_blastdbcmd.zip')
        if _blastn:
            bins_to_dl.append('Linux_blastn.zip')
        if _blastp:
            bins_to_dl.append('Linux_blastp.zip')
    elif current_os.startswith('win'):
        if _blastdcmd:
            bins_to_dl.append('Win32_blastdbcmd.zip')
        if _blastn:
            bins_to_dl.append('Win32_blastn.zip')
        if _blastp:
            bins_to_dl.append('Win32_blastp.zip')
    else:
        return False

    file_to_name = {'Darwin_blastdbcmd.zip': 'blastdbcmd', 'Darwin_blastn.zip': 'blastn',
                    'Darwin_blastp.zip': 'blastp', 'Linux_blastdbcmd.zip': 'blastdbcmd',
                    'Linux_blastn.zip': 'blastn', 'Linux_blastp.zip': 'blastp',
                    'Win32_blastdbcmd.zip': 'blastdbcmd', 'Win32_blastn.zip': 'blastn',
                    'Win32_blastp.zip': 'blastp'}

    os.makedirs("{0}/__tempdir__".format(current_path), exist_ok=True)
    try:
        for blast_bin in bins_to_dl:
            with request.urlopen('{0}{1}'.format(binary_source, blast_bin)) as reader, \
                    open("{0}/__tempdir__/{1}".format(current_path, blast_bin), mode='wb') as writer:
                shutil.copyfileobj(reader, writer)
            zip_file = zipfile.ZipFile("{0}/__tempdir__/{1}".format(current_path, blast_bin))
            zip_file.extractall(path=current_path)
            os.rename('{0}/{1}'.format(current_path, re.sub('\.zip', '', blast_bin)),
                      '{0}/{1}'.format(current_path, file_to_name[blast_bin]))
            os.chmod('{0}/{1}'.format(current_path, file_to_name[blast_bin]), 0o755)
            print("File added: {0}/{1}".format(current_path, file_to_name[blast_bin]))
        shutil.rmtree("{0}/__tempdir__".format(current_path))
    except error.URLError:
        return False

    # TODO account for manual install w/ symlinks
    return True


def _feature_rc(_feature, seq_len):  # BioPython does not properly handle reverse complement of features, so implement..
    if type(_feature.location) == CompoundLocation:
        new_compound_location = []
        for sub_feature in _feature.location.parts:
            sub_feature = _feature_rc(SeqFeature(sub_feature), seq_len)
            new_compound_location.append(sub_feature.location)
        _feature.location = CompoundLocation(new_compound_location, _feature.location.operator)

    elif type(_feature.location) == FeatureLocation:
        _end = seq_len - _feature.location.end
        _shift = _end - _feature.location.start
        _feature = _shift_features(_feature, _shift, seq_len)[0]
        _feature.strand *= -1
    else:
        raise TypeError("_feature_rc requires a feature with either FeatureLocation or CompoundLocation, "
                        "not %s" % type(_feature.location))
    return _feature


def _format_to_extension(_format):
    format_to_extension = {'fasta': 'fa', 'fa': 'fa', 'genbank': 'gb', 'gb': 'gb', 'nexus': 'nex',
                           'nex': 'nex', 'phylip': 'phy', 'phy': 'phy', 'phylip-relaxed': 'phyr', 'phyr': 'phyr',
                           'stockholm': 'stklm', 'stklm': 'stklm'}
    return format_to_extension[_format]


# Does not attempt to explicitly deal with weird cases (e.g., ambiguous residues).
# The user will need to specify an alphabet with the -a flag if using many non-standard characters in their sequences.
def _guess_alphabet(_seqbuddy):
    _seq_list = _seqbuddy if isinstance(_seqbuddy, list) else _seqbuddy.records
    _seq_list = [str(x.seq) for x in _seq_list]
    _sequence = "".join(_seq_list).upper()
    _sequence = re.sub("[NX\-?]", "", _sequence)

    if len(_sequence) == 0:
        return None

    if 'U' in _sequence:  # U is unique to RNA
        return IUPAC.ambiguous_rna

    percent_dna = len(re.findall("[ATCG]", _sequence)) / float(len(_sequence))
    if percent_dna > 0.85:  # odds that a sequence with no Us and such a high ATCG count be anything but DNA is low
        return IUPAC.ambiguous_dna
    else:
        return IUPAC.protein


def _guess_format(_input):  # _input can be list, SeqBuddy object, file handle, or file path.
    # If input is just a list, there is no BioPython in-format. Default to gb.
    if isinstance(_input, list):
        return "gb"

    # Pull value directly from object if appropriate
    if type(_input) == SeqBuddy:
        return _input.in_format

    # If input is a handle or path, try to read the file in each format, and assume success if not error and # seqs > 0
    if os.path.isfile(str(_input)):
        _input = open(_input, "r")

    if str(type(_input)) == "<class '_io.TextIOWrapper'>" or isinstance(_input, StringIO):
        # Die if file is empty
        if _input.read() == "":
            sys.exit("Input file is empty.")
        _input.seek(0)

        possible_formats = ["phylip-relaxed", "stockholm", "fasta", "gb", "fastq", "nexus", "embl", "seqxml"]  # ToDo: Glean CLUSTAL
        for _format in possible_formats:
            try:
                _input.seek(0)
                _seqs = SeqIO.parse(_input, _format)
                if next(_seqs):
                    _input.seek(0)
                    return _format
                else:
                    continue
            except StopIteration:  # ToDo check that other types of error are not possible
                continue
            except ValueError:
                continue
            except SAXParseException:  # Thrown by seqxml parser
                continue
        return None  # Unable to determine format from file handle

    else:
        raise GuessError("Unsupported _input argument in guess_format(). %s" % _input)


def _make_copies(_seqbuddy):
    alphabet_list = [_rec.seq.alphabet for _rec in _seqbuddy.records]
    copies = deepcopy(_seqbuddy)
    copies.alpha = _seqbuddy.alpha
    for _indx, _rec in enumerate(copies.records):
        _rec.seq.alphabet = alphabet_list[_indx]
    return copies


def _phylipi(_seqbuddy, _format="relaxed"):  # _format in ["strict", "relaxed"]
    max_id_length = 0
    max_seq_length = 0
    for _rec in _seqbuddy.records:
        max_id_length = len(_rec.id) if len(_rec.id) > max_id_length else max_id_length
        max_seq_length = len(_rec.seq) if len(_rec.seq) > max_seq_length else max_seq_length

    _output = " %s %s\n" % (len(_seqbuddy.records), max_seq_length)
    for _rec in _seqbuddy.records:
        _seq_id = _rec.id.ljust(max_id_length) if _format == "relaxed" else _rec.id[:10].ljust(10)
        _output += "%s  %s\n" % (_seq_id, _rec.seq)

    return _output


def _shift_features(_features, _shift, full_seq_len):  # shift is an int, how far the new feature should move from 0
    if type(_features) != list:  # Duck type for single feature input
        _features = [_features]

    shifted_features = []
    for _feature in _features:
        if type(_feature.location) == CompoundLocation:  # Recursively call _shift_features() for compound locations
            new_compound_location = []
            for sub_feature in _feature.location.parts:
                sub_feature = _shift_features(SeqFeature(sub_feature), _shift, full_seq_len)
                if not sub_feature:
                    continue
                new_compound_location.append(sub_feature[0].location)

            if not new_compound_location:
                continue

            elif len(new_compound_location) == 1:
                _feature.location = new_compound_location[0]

            else:
                _feature.location = CompoundLocation(new_compound_location, _feature.location.operator)

        elif type(_feature.location) == FeatureLocation:
            _start = _feature.location.start + _shift
            _end = _feature.location.end + _shift
            if _start > full_seq_len or _end < 0:
                continue

            _start = _start if _start >= 0 else 0
            _end = _end if _end <= full_seq_len else full_seq_len

            _feature.location = FeatureLocation(_start, _end, _feature.strand)

        else:
            raise TypeError("_shift_feature requires a feature with either FeatureLocation or CompoundLocation, "
                            "not %s" % type(_feature.location))
        shifted_features.append(_feature)

    return shifted_features


def _stderr(message, quiet=False):
    if not quiet:
        sys.stderr.write(message)
    return


def _stdout(message, quiet=False):
    if not quiet:
        sys.stdout.write(message)
    return


# ################################################ MAIN API FUNCTIONS ################################################ #
def add_feature(_seqbuddy, _type, _location, _strand=None, _qualifiers=None, _pattern=None):
    """
    Adds a feature annotation to all sequences in the SeqBuddy object
    :param _seqbuddy: The SeqBuddy object to be annotated
    :param _type: The type attribute of tha annotation
    :param _location: The location of the annotation - (start, end) or [(start1, end1), (start2, end2)]
    :param _strand: The feature's strand - (+/-/None)
    :param _qualifiers: A dictionary of qualifiers, or a string "foo: bar, fizz: buzz"
    :param _pattern: A regex pattern to specify which sequences to add the feature to
    :return: The annotated SeqBuddy object
    """
    # http://www.insdc.org/files/feature_table.html
    old = _make_copies(_seqbuddy)
    if _pattern is not None:
        recs = pull_recs(_seqbuddy, _pattern).records
    else:
        recs = _seqbuddy.records

    for _rec1 in recs:
        for _rec2 in old.records:
            if _rec1.id == _rec2.id:
                old.records.remove(_rec2)

    if isinstance(_location, FeatureLocation):
        pass
    elif isinstance(_location, CompoundLocation):
        pass
    elif isinstance(_location, list) or isinstance(_location, tuple):
        _locations = []
        if isinstance(_location[0], int):
            _locations.append(FeatureLocation(start=_location[0], end=_location[1]))
        elif isinstance(_location[0], tuple) or isinstance(_location[0], list):
            for _tup in _location:
                _locations.append(FeatureLocation(start=_tup[0], end=_tup[1]))
        elif isinstance(_location[0], str):
            for substr in _location:
                substr = re.sub('[ ()]', '', substr)
                substr = re.sub('-|\.\.', ',', substr)
                _locations.append(FeatureLocation(start=int(re.split(',', substr)[0]),
                                                  end=int(re.split(',', substr)[1])))
        _location = CompoundLocation(sorted(_locations, key=lambda x: x.start), operator='order') \
            if len(_locations) > 1 else _locations[0]
    elif isinstance(_location, str):
        _location = re.sub('[ ()]', '', _location)
        _location = re.split(',', _location)
        _locations = []
        for substr in _location:
            _locations.append(FeatureLocation(start=int(substr.split('-')[0]), end=int(substr.split('-')[1])))
        _location = CompoundLocation(sorted(_locations, key=lambda x: x.start), operator='order') \
            if len(_locations) > 1 else _locations[0]
    else:
        raise TypeError("Input must be list, tuple, or string. Not {0}.".format(type(_location)))

    if _strand in ['+', 'plus', 'sense', 'pos', 'positive', '1', 1]:
        _strand = 1
    elif _strand in ['-', 'minus', 'anti', 'antisense', 'anti-sense', 'neg', 'negative', '-1', -1]:
        _strand = -1
    elif _strand in ['0', 0]:
        _strand = 0
    elif _strand is None:
        pass
    else:
        _strand = None
        _stderr("Warning: _strand input not recognized. Value set to None.")

    if isinstance(_qualifiers, dict):
        pass
    elif isinstance(_qualifiers, str):
        _qualifiers = re.sub(' ', '', _qualifiers)
        _qualifiers = re.sub('=', ':', _qualifiers)
        _qualifiers = re.split(',', _qualifiers)
        qual_dict = {}
        for substr in _qualifiers:
            qual_dict[re.split(':', substr)[0]] = [re.split(':', substr)[1]]
        _qualifiers = qual_dict

    for _rec in recs:
        _rec.features.append(SeqFeature(location=_location, type=_type, strand=_strand, qualifiers=_qualifiers))

    _seqbuddy.records = recs
    _seqbuddy = merge([old, _seqbuddy])
    return _seqbuddy


def ave_seq_length(_seqbuddy, _clean=False):
    """
    Returns the value of the average sequence length
    :param _seqbuddy: The SeqBuddy object to be averaged
    :param _clean: Specifies if non-sequence characters should be counted as well.
    :return: The value of the average sequence length
    """
    if _clean:  # Strip out all gaps and stuff before counting
        clean_seq(_seqbuddy)

    sum_length = 0.
    for _rec in _seqbuddy.records:
        sum_length += len(_rec.seq)
    return sum_length / len(_seqbuddy.records)


def back_translate(_seqbuddy, _mode='random', _species=None):
    """
    Back-translates protein sequences into DNA sequences
    :param _seqbuddy: The SeqBuddy object to be back-translated
    :param _mode: The codon selection mode (random/optimized)
    :param _species: The model to use for optimized codon selection (human/mouse/yeast/ecoli)
    :return: The back-translated SeqBuddy object
    """
    # available modes --> random, optimized
    # available species --> human, mouse, yeast, ecoli
    # codon preference tables derived from the data at http://www.kazusa.or.jp
    # Homo sapiens, species=9606
    if _mode.upper() not in ['RANDOM', 'R', 'OPTIMIZED', 'O']:
        raise AttributeError("Back_translate modes accepted are 'random' or 'r' and 'optimized' or 'o'. "
                             "You entered '%s'" % _mode)

    h_sapi = {'A': (['GCT', 'GCC', 'GCA', 'GCG'], [0.27, 0.40, 0.23, 0.11]),
              'C': (['TGT', 'TGC'], [0.46, 0.54]),
              'D': (['GAT', 'GAC'], [0.46, 0.54]),
              'E': (['GAA', 'GAG'], [0.42, 0.58]),
              'F': (['TTT', 'TTC'], [0.46, 0.54]),
              'G': (['GGT', 'GGC', 'GGA', 'GGG'], [0.16, 0.34, 0.25, 0.25]),
              'H': (['CAT', 'CAC'], [0.42, 0.58]),
              'I': (['ATT', 'ATC', 'ATA'], [0.36, 0.47, 0.17]),
              'K': (['AAA', 'AAG'], [0.43, 0.57]),
              'L': (['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], [0.08, 0.13, 0.13, 0.20, 0.07, 0.40]),
              'M': (['ATG'], [1.00]),
              'N': (['AAT', 'AAC'], [0.47, 0.53]),
              'P': (['CCT', 'CCC', 'CCA', 'CCG'], [0.29, 0.32, 0.28, 0.11]),
              'Q': (['CAA', 'CAG'], [0.27, 0.73]),
              'R': (['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], [0.08, 0.18, 0.11, 0.20, 0.21, 0.21]),
              'S': (['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], [0.19, 0.22, 0.15, 0.05, 0.15, 0.24]),
              '*': (['TAA', 'TGA', 'TAG'], [0.30, 0.47, 0.24]),
              'T': (['ACT', 'ACC', 'ACA', 'ACG'], [0.25, 0.36, 0.28, 0.11]),
              'V': (['GTT', 'GTC', 'GTA', 'GTG'], [0.18, 0.24, 0.12, 0.46]),
              'W': (['TGG'], [1.00]),
              'Y': (['TAT', 'TAC'], [0.44, 0.56]),
              'X': (['NNN'], [1.0])}

    # Mus musculus, species=10090
    m_muscul = {'A': (['GCT', 'GCC', 'GCA', 'GCG'], [0.29, 0.38, 0.23, 0.09]),
                'C': (['TGT', 'TGC'], [0.48, 0.52]),
                'D': (['GAT', 'GAC'], [0.45, 0.55]),
                'E': (['GAA', 'GAG'], [0.41, 0.59]),
                'F': (['TTT', 'TTC'], [0.44, 0.56]),
                'G': (['GGT', 'GGC', 'GGA', 'GGG'], [0.18, 0.33, 0.26, 0.23]),
                'H': (['CAT', 'CAC'], [0.41, 0.59]),
                'I': (['ATT', 'ATC', 'ATA'], [0.34, 0.50, 0.16]),
                'K': (['AAA', 'AAG'], [0.39, 0.61]),
                'L': (['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], [0.07, 0.13, 0.13, 0.20, 0.08, 0.39]),
                'M': (['ATG'], [1.00]),
                'N': (['AAT', 'AAC'], [0.43, 0.57]),
                'P': (['CCT', 'CCC', 'CCA', 'CCG'], [0.31, 0.30, 0.29, 0.10]),
                'Q': (['CAA', 'CAG'], [0.26, 0.74]),
                'R': (['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], [0.08, 0.17, 0.12, 0.19, 0.22, 0.22]),
                'S': (['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], [0.20, 0.22, 0.14, 0.05, 0.15, 0.24]),
                '*': (['TAA', 'TGA', 'TAG'], [0.28, 0.49, 0.23]),
                'T': (['ACT', 'ACC', 'ACA', 'ACG'], [0.25, 0.35, 0.29, 0.10]),
                'V': (['GTT', 'GTC', 'GTA', 'GTG'], [0.17, 0.25, 0.12, 0.46]),
                'W': (['TGG'], [1.00]),
                'Y': (['TAT', 'TAC'], [0.43, 0.57]),
                'X': (['NNN'], [1.0])}

    # Escherichia coli O157:H7 EDL933, species=155864
    e_coli = {'A': (['GCT', 'GCC', 'GCA', 'GCG'], [0.16, 0.27, 0.22, 0.35]),
              'C': (['TGT', 'TGC'], [0.45, 0.55]),
              'D': (['GAT', 'GAC'], [0.63, 0.37]),
              'E': (['GAA', 'GAG'], [0.68, 0.32]),
              'F': (['TTT', 'TTC'], [0.58, 0.42]),
              'G': (['GGT', 'GGC', 'GGA', 'GGG'], [0.33, 0.39, 0.12, 0.16]),
              'H': (['CAT', 'CAC'], [0.58, 0.42]),
              'I': (['ATT', 'ATC', 'ATA'], [0.50, 0.40, 0.09]),
              'K': (['AAA', 'AAG'], [0.76, 0.24]),
              'L': (['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], [0.13, 0.13, 0.11, 0.10, 0.04, 0.49]),
              'M': (['ATG'], [1.00]),
              'N': (['AAT', 'AAC'], [0.47, 0.53]),
              'P': (['CCT', 'CCC', 'CCA', 'CCG'], [0.17, 0.13, 0.19, 0.51]),
              'Q': (['CAA', 'CAG'], [0.33, 0.67]),
              'R': (['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], [0.36, 0.37, 0.07, 0.11, 0.05, 0.03]),
              'S': (['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], [0.14, 0.15, 0.14, 0.15, 0.16, 0.27]),
              '*': (['TAA', 'TGA', 'TAG'], [0.59, 0.33, 0.08]),
              'T': (['ACT', 'ACC', 'ACA', 'ACG'], [0.17, 0.41, 0.15, 0.27]),
              'V': (['GTT', 'GTC', 'GTA', 'GTG'], [0.26, 0.21, 0.16, 0.37]),
              'W': (['TGG'], [1.00]),
              'Y': (['TAT', 'TAC'], [0.57, 0.43]),
              'X': (['NNN'], [1.0])}

    # Saccharomyces cerevisiae, species=4932
    s_cerev = {'A': (['GCT', 'GCC', 'GCA', 'GCG'], [0.38, 0.22, 0.29, 0.11]),
               'C': (['TGT', 'TGC'], [0.63, 0.37]),
               'D': (['GAT', 'GAC'], [0.65, 0.35]),
               'E': (['GAA', 'GAG'], [0.70, 0.30]),
               'F': (['TTT', 'TTC'], [0.59, 0.41]),
               'G': (['GGT', 'GGC', 'GGA', 'GGG'], [0.47, 0.19, 0.22, 0.12]),
               'H': (['CAT', 'CAC'], [0.64, 0.36]),
               'I': (['ATT', 'ATC', 'ATA'], [0.46, 0.26, 0.27]),
               'K': (['AAA', 'AAG'], [0.58, 0.42]),
               'L': (['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], [0.28, 0.29, 0.13, 0.06, 0.14, 0.11]),
               'M': (['ATG'], [1.00]),
               'N': (['AAT', 'AAC'], [0.59, 0.41]),
               'P': (['CCT', 'CCC', 'CCA', 'CCG'], [0.31, 0.15, 0.42, 0.12]),
               'Q': (['CAA', 'CAG'], [0.69, 0.31]),
               'R': (['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], [0.14, 0.06, 0.07, 0.04, 0.48, 0.21]),
               'S': (['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], [0.26, 0.16, 0.21, 0.10, 0.16, 0.11]),
               '*': (['TAA', 'TGA', 'TAG'], [0.47, 0.30, 0.23]),
               'T': (['ACT', 'ACC', 'ACA', 'ACG'], [0.35, 0.22, 0.30, 0.14]),
               'V': (['GTT', 'GTC', 'GTA', 'GTG'], [0.39, 0.21, 0.21, 0.19]),
               'W': (['TGG'], [1.00]),
               'Y': (['TAT', 'TAC'], [0.56, 0.44]),
               'X': (['NNN'], [1.0])}

    # random
    rand_table = {'A': (['GCT', 'GCC', 'GCA', 'GCG'], [0.25, 0.25, 0.25, 0.25]),
                  'C': (['TGT', 'TGC'], [0.5, 0.5]),
                  'D': (['GAT', 'GAC'], [0.5, 0.5]),
                  'E': (['GAA', 'GAG'], [0.5, 0.5]),
                  'F': (['TTT', 'TTC'], [0.5, 0.5]),
                  'G': (['GGT', 'GGC', 'GGA', 'GGG'], [0.25, 0.25, 0.25, 0.25]),
                  'H': (['CAT', 'CAC'], [0.5, 0.5]),
                  'I': (['ATT', 'ATC', 'ATA'], [0.3333, 0.3333, 0.3334]),
                  'K': (['AAA', 'AAG'], [0.5, 0.5]),
                  'L': (['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], [0.167, 0.167, 0.167, 0.167, 0.166, 0.166]),
                  'M': (['ATG'], [1.00]),
                  'N': (['AAT', 'AAC'], [0.5, 0.5]),
                  'P': (['CCT', 'CCC', 'CCA', 'CCG'], [0.25, 0.25, 0.25, 0.25]),
                  'Q': (['CAA', 'CAG'], [0.5, 0.5]),
                  'R': (['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], [0.167, 0.167, 0.167, 0.167, 0.166, 0.166]),
                  'S': (['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], [0.167, 0.167, 0.167, 0.167, 0.166, 0.166]),
                  '*': (['TAA', 'TGA', 'TAG'], [0.3333, 0.3333, 0.3334]),
                  'T': (['ACT', 'ACC', 'ACA', 'ACG'], [0.25, 0.25, 0.25, 0.25]),
                  'V': (['GTT', 'GTC', 'GTA', 'GTG'], [0.25, 0.25, 0.25, 0.25]),
                  'W': (['TGG'], [1.00]),
                  'Y': (['TAT', 'TAC'], [0.5, 0.5]),
                  'X': (['NNN'], [1.0])}

    if _seqbuddy.alpha != IUPAC.protein:
        raise TypeError("The input sequence needs to be IUPAC.protein'>, not %s" %
                        str(type(_seqbuddy.alpha)))

    if not _species:
        lookup_table = rand_table
    elif _species.upper() in ["HUMAN", "H"]:
        lookup_table = h_sapi
    elif _species.upper() in ["MOUSE", "M"]:
        lookup_table = m_muscul
    elif _species.upper() in ["ECOLI", "E"]:
        lookup_table = e_coli
    elif _species.upper() in ["YEAST", "Y"]:
        lookup_table = s_cerev
    else:
        raise AttributeError("The species requested does not match any lookup tables currently implemented. "
                             "Please leave blank or select from human, mouse, ecoli, or yeast.")

    if _mode.upper() in ['OPTIMIZED', 'O']:
        for aa in lookup_table:
            best = ["", 0.]
            for _i in range(len(lookup_table[aa][1])):
                if lookup_table[aa][1][_i] > best[1]:
                    best = [lookup_table[aa][0][_i], lookup_table[aa][1][_i]]
            lookup_table[aa] = ([best[0]], [1.0])
    originals = deepcopy(_seqbuddy)
    for _rec in _seqbuddy.records:
        _rec.features = []
        dna_seq = ""
        for aa in _rec.seq.upper():
            rand_num = random()
            sum_probs = 0.
            for _i in range(len(lookup_table[aa][1])):
                sum_probs += lookup_table[aa][1][_i]
                if sum_probs >= rand_num:
                    dna_seq += lookup_table[aa][0][_i]
                    break
            _rec.seq = Seq(dna_seq, alphabet=IUPAC.ambiguous_dna)

    mapped_features = map_features_prot2dna(originals, _seqbuddy)
    mapped_features.out_format = _seqbuddy.out_format
    return mapped_features


def bl2seq(_seqbuddy):  # TODO do string formatting in command line ui
    """
    Does an all-by-all analysis of the sequences
    :param _seqbuddy: The SeqBuddy object to be BLASTed
    :return: A tuple containing a table of results and a string output
    """
    # Note on blast2seq: Expect (E) values are calculated on an assumed database size of (the rather large) nr, so the
    # threshold may need to be increased quite a bit to return short alignments

    current_dir = os.getcwd()
    script_location = os.path.realpath(__file__)
    script_location = re.sub('SeqBuddy\.py', '', script_location)

    if not which("blastp") and _seqbuddy.alpha not in [IUPAC.protein]:
        _stderr("Blastp binary not found. Would you like to download it? (program will be aborted) [yes]/no\n")
        prompt = input()
        while True:
            if prompt.lower() in ['yes', 'y', '']:
                os.chdir(script_location)
                if _download_blast_binaries(_blastdcmd=False, _blastn=False, _blastp=True):
                    _stderr("Blastp downloaded.\n")
                else:
                    _stderr("Failed to download blastp.\n")
                break
            elif prompt.lower() in ['no', 'n']:
                break
            else:
                _stderr("Input not understood.\n")
                _stderr("Would you like to download blastp? (program will be aborted) [yes]/no\n")
                prompt = input()
        os.chdir(current_dir)
        if not which("blastp"):
            raise RuntimeError("Blastp not present in $PATH or working directory.")
        sys.exit()

    elif not which("blastn") and _seqbuddy.alpha not in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna,
                                                         IUPAC.ambiguous_rna, IUPAC.unambiguous_rna]:
        _stderr("Blastn binary not found. Would you like to download it? (program will be aborted) [yes]/no\n")
        prompt = input()
        while True:
            if prompt.lower() in ['yes', 'y', '']:
                os.chdir(script_location)
                if _download_blast_binaries(_blastdcmd=False, _blastn=False, _blastp=True):
                    _stderr("Blastn downloaded.\n")
                else:
                    _stderr("Failed to download blastn.\n")
                break
            elif prompt.lower() in ['no', 'n']:
                break
            else:
                _stderr("Input not understood.\n")
                _stderr("Would you like to download blastn? (program will be aborted) [yes]/no\n")
                prompt = input()
        os.chdir(current_dir)
        if not which("blastn"):
            raise RuntimeError("Blastn not present in $PATH or working directory.")
        sys.exit()

    def mc_blast(_query, _args):
        _subject_file = _args[0]

        if subject.id == _query.id:
            return

        _blast_res = Popen("echo '%s' | %s -subject %s -outfmt 6" %
                           (_query.format("fasta"), blast_bin, _subject_file), stdout=PIPE, shell=True).communicate()
        _blast_res = _blast_res[0].decode().split("\n")[0].split("\t")
        _result = ""
        while True:
            try:
                if len(_blast_res) == 1:
                    _result = "%s\t%s\t0\t0\t0\t0\n" % (subject.id, _query.id)
                else:
                    # values are: query, subject, %_ident, length, evalue, bit_score
                    if _blast_res[10] == '0.0':
                        _blast_res[10] = '1e-180'
                    _result = "%s\t%s\t%s\t%s\t%s\t%s\n" % (_blast_res[0], _blast_res[1], _blast_res[2],
                                                            _blast_res[3], _blast_res[10], _blast_res[11].strip())
                break
            except ConnectionRefusedError:
                continue
        with lock:
            with open("%s/blast_results.txt" % tmp_dir.name, "a") as _ofile:
                _ofile.write(_result)
        return

    blast_bin = "blastp" if _seqbuddy.alpha == IUPAC.protein else "blastn"
    if not which(blast_bin):
        raise RuntimeError("%s not present in $PATH or working directory." % blast_bin)  # ToDo: Implement -p flag

    from multiprocessing import Lock
    lock = Lock()
    tmp_dir = TemporaryDirectory()

    # Copy the seqbuddy records into new list, so they can be iteratively deleted below
    _seqs_copy = _seqbuddy.records[:]
    subject_file = "%s/subject.fa" % tmp_dir.name
    for subject in _seqbuddy.records:
        with open(subject_file, "w") as ifile:
            SeqIO.write(subject, ifile, "fasta")

        run_multicore_function(_seqs_copy, mc_blast, [subject_file], out_type=sys.stderr, quiet=True)  # Todo Benchmark

        _seqs_copy = _seqs_copy[1:]

    with open("%s/blast_results.txt" % tmp_dir.name, "r") as _ifile:
        output_list = _ifile.read().strip().split("\n")

    # Push output into a dictionary of dictionaries, for more flexible use outside of this function
    output_list = [x.split("\t") for x in output_list]
    output_list = sorted(output_list, key=lambda l: l[0])
    output_dict = {}
    for match in output_list:
        query, subj, _ident, _length, _evalue, _bit_score = match
        if query not in output_dict:
            output_dict[query] = {subj: [float(_ident), int(_length), float(_evalue), float(_bit_score)]}
        else:
            output_dict[query][subj] = [float(_ident), int(_length), float(_evalue), float(_bit_score)]

        if subj not in output_dict:
            output_dict[subj] = {query: [float(_ident), int(_length), float(_evalue), float(_bit_score)]}
        else:
            output_dict[subj][query] = [float(_ident), int(_length), float(_evalue), float(_bit_score)]

    output_str = "#query\tsubject\t%_ident\tlength\tevalue\tbit_score\n"

    output_list = [(_key, _value) for _key, _value in output_dict.items()]
    output_list = sorted(output_list, key=lambda l: l[0])

    ids_already_seen = []
    for query_id, query_values in output_list:
        ids_already_seen.append(query_id)
        query_values = [(_key, _value) for _key, _value in query_values.items()]
        query_values = sorted(query_values, key=lambda l: l[0])
        for subj_id, subj_values in query_values:
            if subj_id in ids_already_seen:
                continue

            ident, length, evalue, bit_score = subj_values
            output_str += "%s\t%s\t%s\t%s\t%s\t%s\n" % (query_id, subj_id, ident, length, evalue, bit_score)

    return [output_dict, output_str]


def blast(_seqbuddy, blast_db, blast_path=None, blastdbcmd=None):  # ToDo: Allow weird binary names to work
    """
    Runs a BLAST search for all of the sequences in the SeqBuddy object against a specified database.
    :param _seqbuddy: The SeqBuddy object containing the sequences to be BLASTed
    :param blast_db: The location of the BLAST database to run the sequences against
    :param blast_path: The location of the blastn/blastp executable
    :param blastdbcmd: The location of the blastdbcmd executable
    :return: A SeqBuddy object containing all of the BLAST database matches
    """
    if not blast_path:
        blast_path = which("blastp") if _seqbuddy.alpha == IUPAC.protein else which("blastn")

    current_dir = os.getcwd()
    script_location = os.path.realpath(__file__)
    script_location = re.sub(str(__file__), '', script_location)

    blast_check = Popen("%s -version" % blast_path, stdout=PIPE, shell=True).communicate()
    blast_check = re.search("([a-z])*[^:]", blast_check[0].decode("utf-8"))
    if blast_check:
        blast_check = blast_check.group(0)

    extensions = {"blastp": ["phr", "pin", "pog", "psd", "psi", "psq"],
                  "blastn": ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]}

    # Try to catch the common variations of the database names that might be given as input
    if blast_db[-2:] in [".p", ".n"]:
        blast_db = blast_db[:-2]

    if blast_db[-3:] in extensions[blast_check]:
        blast_db = blast_db[:-4]

    blast_db = os.path.abspath(blast_db)

    # ToDo Check NCBI++ tools are a conducive version (2.2.29 and above, I think [maybe .28])
    # Check to make sure blast is in $PATH and ensure that the blast_db is present

    if blast_check == "blastp":
        if not which(blast_path):
            _stderr("Blastp binary not found. Would you like to download it? (program will be aborted) [yes]/no\n")
            prompt = input()
            while True:
                if prompt.lower() in ['yes', 'y', '']:
                    os.chdir(script_location)
                    if _download_blast_binaries(_blastdcmd=False, _blastn=False, _blastp=True):
                        _stderr("Blastp downloaded.\n")
                    else:
                        _stderr("Failed to download blastp.\n")
                    break
                elif prompt.lower() in ['no', 'n']:
                    break
                else:
                    _stderr("Input not understood.\n")
                    _stderr("Would you like to download blastp? (program will be aborted) [yes]/no\n")
                    prompt = input()
            os.chdir(current_dir)
            if not which("blastp"):
                raise FileNotFoundError("blastp binary not found")
            return

        if not os.path.isfile("%s.pin" % blast_db) or not os.path.isfile("%s.phr" % blast_db) \
                or not os.path.isfile("%s.psq" % blast_db):
            raise RuntimeError("Blastp database not found at '%s'" % blast_db)
    elif blast_check == "blastn":
        if not which(blast_path):
            _stderr("Blastn binary not found. Would you like to download it? (program will be aborted) [yes]/no\n")
            prompt = input()
            while True:
                if prompt.lower() in ['yes', 'y', '']:
                    os.chdir(script_location)
                    if _download_blast_binaries(_blastdcmd=False, _blastn=True, _blastp=False):
                        _stderr("Blastn downloaded.\n")
                    else:
                        _stderr("Failed to download blastn.\n")
                    break
                elif prompt.lower() in ['no', 'n']:
                    break
                else:
                    _stderr("Input not understood.\n")
                    _stderr("Would you like to download blastn? (program will be aborted) [yes]/no\n")
                    prompt = input()
            os.chdir(current_dir)
            if not which("blastn"):
                raise FileNotFoundError("blastn binary not found")
            return

        if not os.path.isfile("%s.nin" % blast_db) or not os.path.isfile("%s.nhr" % blast_db) \
                or not os.path.isfile("%s.nsq" % blast_db):
            raise RuntimeError("Blastn database not found at '%s'" % blast_db)
    else:
        raise RuntimeError("Blast binary doesn't seem to work, at %s" % blast_path)

    if not blastdbcmd:
        blastdbcmd = "blastdbcmd"

    if not which(blastdbcmd):
        _stderr("Blastdbcmd binary not found. Would you like to download it? (program will be aborted) [yes]/no\n")
        prompt = input()
        while True:
            if prompt.lower() in ['yes', 'y', '']:
                os.chdir(script_location)
                if _download_blast_binaries(_blastdcmd=True, _blastn=False, _blastp=False):
                    _stderr("Blastdbcmd downloaded.\n")
                else:
                    _stderr("Failed to download blastdbcmd.\n")
                break
            elif prompt.lower() in ['no', 'n']:
                break
            else:
                _stderr("Input not understood.\n")
                _stderr("Would you like to download blastdbcmd? (program will be aborted) [yes]/no\n")
                prompt = input()
        os.chdir(current_dir)
        if not which("blastdbcmd"):
            raise FileNotFoundError("blastdbcmd")
        return

    # Check that compelte blastdb is present and was made with the -parse_seqids flag
    for extension in extensions[blast_check]:
        if not os.path.isfile("%s.%s" % (blast_db, extension)):
            raise RuntimeError("The .%s file of your blast database was not found. Ensure the -parse_seqids flag was "
                               "used with makeblastdb." % extension)

    _seqbuddy = clean_seq(_seqbuddy)  # in case there are gaps or something in the sequences

    tmp_dir = TemporaryDirectory()
    with open("%s/tmp.fa" % tmp_dir.name, "w") as _ofile:
        SeqIO.write(_seqbuddy.records, _ofile, "fasta")

    Popen("%s -db %s -query %s/tmp.fa -out %s/out.txt -num_threads 4 -evalue 0.01 -outfmt 6" %
          (blast_path, blast_db, tmp_dir.name, tmp_dir.name), shell=True).wait()

    with open("%s/out.txt" % tmp_dir.name, "r") as ifile:
        blast_results = ifile.read()
        _records = blast_results.split("\n")

    hit_ids = []
    for record in _records:
        record = record.split("\t")
        if len(record) == 1:
            continue
        hit_id = record[1].strip()
        if hit_id in hit_ids:
            continue

        hit_ids.append(hit_id)

    if len(hit_ids) == 0:
        sys.stderr.write("No matches identified.\n")
        return None

    _ofile = open("%s/seqs.fa" % tmp_dir.name, "w")
    for hit_id in hit_ids:
        hit = Popen("blastdbcmd -db %s -entry 'lcl|%s'" % (blast_db, hit_id), stdout=PIPE, shell=True).communicate()
        hit = hit[0].decode("utf-8")
        hit = re.sub("lcl\|", "", hit)
        _ofile.write("%s\n" % hit)

    _ofile.close()

    with open("%s/seqs.fa" % tmp_dir.name, "r") as ifile:
        _new_seqs = SeqBuddy(ifile)

    return _new_seqs


def clean_seq(_seqbuddy, skip_list=None, ambiguous=True):
    """
    Removes all non-sequence characters from the sequences
    :param _seqbuddy: The SeqBuddy object to be cleaned
    :param skip_list: A list of characters to be left alone
    :param ambiguous: Specifies whether ambiguous characters should be kept or not
    :return: The cleaned SeqBuddy object
    """
    skip_list = "" if not skip_list else "".join(skip_list)
    for _rec in _seqbuddy.records:
        if _seqbuddy.alpha == IUPAC.protein:
            full_skip = "ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy%s" % skip_list
            _rec.seq = Seq(re.sub("[^%s]" % full_skip, "", str(_rec.seq)),
                           alphabet=_seqbuddy.alpha)
        else:
            if ambiguous:
                full_skip = "ATGCURYWSMKHBVDNXatgcurywsmkhbvdnx%s" % skip_list
                _rec.seq = Seq(re.sub("[^%s]" % full_skip, "", str(_rec.seq)),
                               alphabet=_seqbuddy.alpha)
            else:
                full_skip = "ATGCUatgcu%s" % skip_list
                _rec.seq = Seq(re.sub("[^%s]" % full_skip, "", str(_rec.seq)), alphabet=_seqbuddy.alpha)

    return _seqbuddy


def combine_features(_seqbuddy1, _seqbuddy2):  # ToDo: rewrite this to accept any number of input files.
    """
    Merges the feature lists of two SeqBuddy objects
    :param _seqbuddy1: The first SeqBuddy object
    :param _seqbuddy2: The second SeqBuddy object
    :return: A SeqBuddy object with merged features
    """
    # make sure there are no repeat ids
    _unique, _rep_ids, _rep_seqs, output_str = find_repeats(_seqbuddy1)
    if len(_rep_ids) > 0:
        raise RuntimeError("There are repeat IDs in the first file provided\n%s" % _rep_ids)

    _unique, _rep_ids, _rep_seqs, output_str = find_repeats(_seqbuddy2)
    if len(_rep_ids) > 0:
        raise RuntimeError("There are repeat IDs in the second file provided\n%s" % _rep_ids)

    seq_dict1 = {}
    seq_dict2 = {}
    seq_order = []

    for _rec in _seqbuddy1.records:
        seq_dict1[_rec.id] = _rec
        seq_order.append(_rec.id)

    for _rec in _seqbuddy2.records:
        seq_dict2[_rec.id] = _rec
        if _rec.id not in seq_order:
            seq_order.append(_rec.id)

    # make sure that we're comparing apples to apples across all sequences (i.e., same alphabet)
    reference_alphabet = sample(seq_dict1.items(), 1)[0][1].seq.alphabet
    for _seq_id in seq_dict1:
        if type(seq_dict1[_seq_id].seq.alphabet) != type(reference_alphabet):
            raise RuntimeError("You have mixed multiple alphabets into your sequences. Make sure everything is the same"
                               "\n\t%s in first set\n\tOffending alphabet: %s\n\tReference alphabet: %s"
                               % (_seq_id, seq_dict1[_seq_id].seq.alphabet, reference_alphabet))

    for _seq_id in seq_dict2:
        if type(seq_dict2[_seq_id].seq.alphabet) != type(reference_alphabet):
            raise RuntimeError("You have mixed multiple alphabets into your sequences. Make sure everything is the same"
                               "\n\t%s in first set\n\tOffending alphabet: %s\n\tReference alphabet: %s"
                               % (_seq_id, seq_dict2[_seq_id].seq.alphabet, reference_alphabet))

    _new_seqs = {}
    warning_used = False
    for _seq_id in seq_dict1:
        if _seq_id in seq_dict2:
            _seq_feats1 = []  # Test list so features common to both records are not duplicated
            for feature in seq_dict1[_seq_id].features:
                _seq_feats1.append("%s-%s-%s" % (feature.location.start, feature.location.end, feature.type))
            for feature in seq_dict2[_seq_id].features:
                feature_check = "%s-%s-%s" % (feature.location.start, feature.location.end, feature.type)
                if feature_check in _seq_feats1:
                    continue
                else:
                    seq_dict1[_seq_id].features.append(feature)
        else:
            warning_used = True
            sys.stderr.write("Warning: %s is only in the first set of sequences\n" % _seq_id)

        _new_seqs[_seq_id] = seq_dict1[_seq_id]

    for _seq_id in seq_dict2:
        if _seq_id not in seq_dict1:
            warning_used = True
            sys.stderr.write("Warning: %s is only in the first set of sequences\n" % _seq_id)
            _new_seqs[_seq_id] = seq_dict2[_seq_id]

    if warning_used:
        sys.stderr.write("\n")

    _seqbuddy = SeqBuddy([_new_seqs[_seq_id] for _seq_id in seq_order], _out_format=_seqbuddy1.in_format)
    _seqbuddy = order_features_by_position(_seqbuddy)
    return _seqbuddy


def complement(_seqbuddy):
    """
    Converts DNA/RNA sequences to their complementary sequence
    :param _seqbuddy: The SeqBuddy object to be modified
    :return: The modified SeqBuddy object
    """
    if _seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for _rec in _seqbuddy.records:
        _rec.seq = _rec.seq.complement()
    return _seqbuddy


def concat_seqs(_seqbuddy, _clean=False):
    """
    Concatenates all of the sequences in the SeqBuddy object into one
    :param _seqbuddy: The SeqBuddy object to be concatenated
    :param _clean: Specifies whether non-sequence characters should be cleaned
    :return: The concatenated SeqBuddy object
    """
    if _clean:
        clean_seq(_seqbuddy)

    _new_seq = ""
    concat_ids = []
    features = []
    for _rec in _seqbuddy.records:
        _shift = len(_new_seq)
        full_seq_len = len(_new_seq) + len(str(_rec.seq))
        _rec.features = _shift_features(_rec.features, _shift, full_seq_len)

        _location = FeatureLocation(len(_new_seq), len(_new_seq) + len(str(_rec.seq)))
        feature = SeqFeature(location=_location, id=_rec.id, type=_rec.id[:15])
        features.append(feature)
        features += _rec.features
        concat_ids.append(_rec.id)
        _new_seq += str(_rec.seq)

    _new_seq = [SeqRecord(Seq(_new_seq, alphabet=_seqbuddy.alpha),
                          description="", id="concatination", features=features)]
    _seqbuddy = SeqBuddy(_new_seq)
    _seqbuddy.out_format = "gb"
    return _seqbuddy


def count_codons(_seqbuddy):
    """
    Generate frequency statistics for codon composition
    :param _seqbuddy: The SeqBuddy object to be analyzed
    :return: A tuple containing the original SeqBuddy object and a dictionary - dict[id][codon] = (Amino acid, num, %)
    """
    if _seqbuddy.alpha not in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna, IUPAC.ambiguous_rna, IUPAC.unambiguous_rna]:
        raise TypeError("Nucleic acid sequence required, not protein or other.")
    if _seqbuddy.alpha in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna]:
        codontable = CodonTable.ambiguous_dna_by_name['Standard'].forward_table
    else:
        codontable = CodonTable.ambiguous_rna_by_name['Standard'].forward_table
    _output = OrderedDict()
    for _rec in _seqbuddy.records:
        _sequence = _rec.seq
        if len(_sequence) % 3 != 0:
            _stderr("Warning: {0} length not a multiple of 3. Sequence will be truncated.\n".format(_rec.id))
            while len(_sequence) % 3 != 0:
                _sequence = _sequence[:-1]
        data_table = OrderedDict()
        num_codons = len(_sequence) / 3
        while len(_sequence) > 0:
            _codon = str(_sequence[:3]).upper()
            if _codon in data_table.keys():
                data_table[_codon][1] += 1
            else:
                if _codon.upper() in ['ATG', 'AUG']:
                    data_table[_codon] = ['M', 1, 0.0]
                elif _codon.upper() == 'NNN':
                    data_table[_codon] = ['X', 1, 0.0]
                elif _codon.upper() in ['TAA', 'TAG', 'TGA', 'UAA', 'UAG', 'UGA']:
                    data_table[_codon] = ['*', 1, 0.0]
                else:
                    try:
                        data_table[_codon] = [codontable[_codon.upper()], 1, 0.0]
                    except KeyError:
                        _stderr("Warning: Codon '{0}' is invalid. Codon will be skipped.\n".format(_codon))
            _sequence = _sequence[3:]
        for _codon in data_table:
            data_table[_codon][2] = round(data_table[_codon][1] / float(num_codons) * 100, 3)
        _output[_rec.id] = OrderedDict(sorted(data_table.items(), key=lambda x: x[0]))
    for _rec in _seqbuddy.records:
        try:
            _rec.buddy_data['Codon_frequency'] = _output[_rec.id]
        except AttributeError:
            _rec.buddy_data = {'Codon_frequency': _output[_rec.id]}
    return _seqbuddy, _output


def count_residues(_seqbuddy):
    """
    Count the number of each type of residue in the SeqBuddy object
    :param _seqbuddy: The SeqBuddy object to be analyzed
    :return: A tuple containing an annotated SeqBuddy object and a dictionary - dict[id][residue]
    """
    _output = OrderedDict()
    for _rec in _seqbuddy.records:
        resid_count = {}
        _seq = str(_rec.seq).upper()
        for char in _seq:
            resid_count.setdefault(char, 0)
            resid_count[char] += 1

        seq_len = len(_rec)
        for _residue, _count in resid_count.items():
            resid_count[_residue] = [_count, _count / seq_len]

        if _seqbuddy.alpha is IUPAC.protein:
            ambig = len(re.findall("X", _seq))
            if ambig > 0:
                resid_count['% Ambiguous'] = round(100 * ambig / seq_len, 2)

            pos = len(re.findall("[HKR]", _seq))
            resid_count['% Positive'] = round(100 * pos / seq_len, 2)

            neg = len(re.findall("[DEC]", _seq))
            resid_count['% Negative'] = round(100 * neg / seq_len, 2)

            neut = len(re.findall("[GAVLIPFYWSTNQM]", _seq))
            resid_count['% Uncharged'] = round(100 * neut / seq_len, 2)

            hyrdophobic = len(re.findall("[AVLIPYFWMC]", _seq))
            resid_count['% Hyrdophobic'] = round(100 * hyrdophobic / seq_len, 2)

            hyrdophilic = len(re.findall("[NQSTKRHDE]", _seq))
            resid_count['% Hyrdophilic'] = round(100 * hyrdophilic / seq_len, 2)

            for _residue in ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
                             "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]:
                resid_count.setdefault(_residue, [0, 0])

        else:
            ambig = len(re.findall("[^ATCGU]", _seq))
            if ambig > 0:
                resid_count['% Ambiguous'] = round(100 * ambig / seq_len, 2)

            for _residue in ["A", "G", "C"]:
                resid_count.setdefault(_residue, [0, 0])

            if "T" not in resid_count and "U" not in resid_count:
                resid_count["T"] = [0, 0]

        try:
            _rec.buddy_data['Residue_frequency'] = resid_count
        except AttributeError:
            _rec.buddy_data = {'Residue_frequency': resid_count}

        _output[_rec.id] = resid_count
    return _seqbuddy, _output


def delete_features(_seqbuddy, _pattern):
    """
    Deletes features with IDs matching a regex pattern
    :param _seqbuddy: The SeqBuddy object to be modified
    :param _pattern: The regex pattern to search with
    :return: The modified SeqBuddy object
    """
    for _rec in _seqbuddy.records:
        retained_features = []
        for _feature in _rec.features:
            if not re.search(_pattern, _feature.type):
                retained_features.append(_feature)
        _rec.features = retained_features
    return _seqbuddy


def delete_large(_seqbuddy, max_value):
    """
    Deletes records larger than a certain size
    :param _seqbuddy: The SeqBuddy object to be modified
    :param max_value: The maximum threshold for sequence length
    :return: The modified SeqBuddy object
    """
    retained_records = []
    for _rec in _seqbuddy.records:
        if len(str(_rec.seq)) <= max_value:
            retained_records.append(_rec)
    _seqbuddy.records = retained_records
    return _seqbuddy


def delete_metadata(_seqbuddy):
    """
    Removes all of the metadata from the records
    :param _seqbuddy: The SeqBuddy object to be stripped of its metadata
    :return: The stripped SeqBuddy object
    """
    _new_seqs = []
    for _rec in _seqbuddy.records:
        _new_seqs.append(SeqRecord(Seq(str(_rec.seq), alphabet=_seqbuddy.alpha), id=_rec.id, name='', description=''))
    _seqbuddy.records = _new_seqs
    return _seqbuddy


def delete_records(_seqbuddy, search_str):
    """
    Deletes records with IDs matching a regex pattern
    :param _seqbuddy: The SeqBuddy object to be modified
    :param search_str: The regex pattern to search with
    :return: The modified SeqBuddy object
    """
    retained_records = []
    _deleted = pull_recs(copy(_seqbuddy), search_str).records
    for _rec in _seqbuddy.records:
        if _rec in _deleted:
            continue
        else:
            retained_records.append(_rec)
    _seqbuddy.records = retained_records
    return _seqbuddy


def delete_repeats(_seqbuddy, scope='all'):  # scope in ['all', 'ids', 'seqs']
    """
    Deletes records with repeated IDs/seqs
    :param _seqbuddy: The SeqBuddy object to be modified
    :param scope: Specifies if deleting repeat seqs, ids, or all
    :return: The modified SeqBuddy object
    """
    # First, remove duplicate IDs
    if scope in ['all', 'ids']:
        _unique, _rep_ids, _rep_seqs, output_str = find_repeats(_seqbuddy)
        if len(_rep_ids) > 0:
            for _rep_id in _rep_ids:
                store_one_copy = pull_recs(copy(_seqbuddy), "^%s$" % _rep_id).records[0]
                delete_records(_seqbuddy, "^%s$" % _rep_id)
                _seqbuddy.records.append(store_one_copy)

    # Then remove duplicate sequences
    if scope in ['all', 'seqs']:
        _unique, _rep_ids, _rep_seqs, output_str = find_repeats(_seqbuddy)
        if len(_rep_seqs) > 0:
            _rep_seq_ids = []
            for _seq in _rep_seqs:
                _rep_seq_ids.append([])
                for _rep_seq_id in _rep_seqs[_seq]:
                    _rep_seq_ids[-1].append(_rep_seq_id)

            repeat_regex = ""

            for _rep_seqs in _rep_seq_ids:
                for _rep_seq in _rep_seqs[1:]:
                    _rep_seq = re.sub("([|.*?^\[\]()])", r"\\\1", _rep_seq)
                    repeat_regex += "^%s$|" % _rep_seq

            repeat_regex = repeat_regex[:-1]
            delete_records(_seqbuddy, repeat_regex)

    return _seqbuddy


def delete_small(_seqbuddy, min_value):
    """
    Deletes records smaller than a certain size
    :param _seqbuddy: The SeqBuddy object to be modified
    :param min_value: The minimum threshold for sequence length
    :return: The modified SeqBuddy object
    """
    retained_records = []
    for _rec in _seqbuddy.records:
        if len(str(_rec.seq)) >= min_value:
            retained_records.append(_rec)
    _seqbuddy.records = retained_records
    return _seqbuddy


def dna2rna(_seqbuddy):
    """
    Back-transcribes DNA into RNA sequences
    :param _seqbuddy: The SeqBuddy object to be back-transcribed
    :return: The back-transcribed SeqBuddy object
    """
    if _seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq.transcribe()), alphabet=IUPAC.ambiguous_rna)
    _seqbuddy.alpha = IUPAC.ambiguous_rna
    return _seqbuddy


def extract_range(_seqbuddy, _start, _end):
    """
    Retrieves subsequences in a specified range
    :param _seqbuddy: The SeqBuddy object to be pulled from
    :param _start: The starting point
    :param _end: The end point
    :return: The modified SeqBuddy object
    """
    _start = 1 if int(_start) < 1 else _start
    # Don't use the standard index-starts-at-0... _end must be left for the range to be inclusive
    _start, _end = int(_start) - 1, int(_end)
    if _end < _start:
        raise ValueError("Error at extract range: The value given for end of range is smaller than for the start "
                         "of range.")

    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq)[_start:_end], alphabet=_rec.seq.alphabet)
        _rec.description += " Sub-sequence extraction, from residue %s to %s" % (_start + 1, _end)
        _features = []
        for _feature in _rec.features:
            if _feature.location.end < _start:
                continue
            if _feature.location.start > _end:
                continue

            feat_start = _feature.location.start - _start
            if feat_start < 0:
                feat_start = 0

            feat_end = _feature.location.end - _start
            if feat_end > len(str(_rec.seq)):
                feat_end = len(str(_rec.seq))

            new_location = FeatureLocation(feat_start, feat_end)
            _feature.location = new_location
            _features.append(_feature)
        _rec.features = _features
    return _seqbuddy


def find_cpg(_seqbuddy):
    """
    Predicts locations of CpG islands in DNA sequences
    :param _seqbuddy: The SeqBuddy object to be analyzed
    :return: A tuple containing an SeqBuddy object and a dict of indices - dict[id]
    """
    _seqbuddy = clean_seq(_seqbuddy)
    if _seqbuddy.alpha not in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna]:
        raise TypeError("DNA sequence required, not protein or RNA.")

    _output = OrderedDict()
    records = []

    def cpg_calc(in_seq):  # Returns observed/expected value of a sequence
        in_seq = in_seq.upper()
        observed_cpg = len(re.findall("CG", in_seq)) * len(in_seq)
        expected = (len(re.findall("[CG]", in_seq)) / 2) ** 2
        return observed_cpg / expected

    def cg_percent(in_seq):  # Returns the CG % of a sequence
        in_seq = in_seq.upper()
        return len(re.findall("[CG]", in_seq)) / len(in_seq)

    def find_islands(cg_percents, oe_values):  # Returns a list of tuples containing the start and end of an island
        _indx = 0
        out_list = []
        while _indx < len(cg_percents):
            if cg_percents[_indx] > .5 and oe_values[_indx] > .6:
                start = _indx
                while _indx < len(cg_percents) and cg_percents[_indx] > .5 and oe_values[_indx] > .6:
                    _indx += 1
                end = _indx
                out_list.append((start, end))
            _indx += 1
        return out_list

    def map_cpg(in_seq, island_ranges):  # Maps CpG islands onto a sequence as capital letters
        cpg_seq = in_seq.lower()
        for pair in island_ranges:
            cpg_seq = cpg_seq[:pair[0]] + cpg_seq[pair[0]:pair[1] + 1].upper() + cpg_seq[pair[1] + 1:]
        return cpg_seq

    for rec in _seqbuddy.records:
        _seq = rec.seq

        oe_vals_list = [0 for _ in range(len(_seq))]
        cg_percent_list = [0 for _ in range(len(_seq))]
        window_size = len(_seq) if len(_seq) < 200 else 200

        for _index, window in enumerate([_seq[i:i + window_size] for i in range(len(_seq) - window_size + 1)]):
            cpg = cpg_calc(str(window))
            for i in range(window_size):
                oe_vals_list[_index + i] += cpg

            cg_perc = cg_percent(str(window))
            for i in range(window_size):
                cg_percent_list[_index + i] += cg_perc

        for _index in range(len(oe_vals_list)):
            if _index + 1 <= window_size:
                oe_vals_list[_index] /= (_index + 1)
                cg_percent_list[_index] /= (_index + 1)

            elif (len(oe_vals_list) - window_size) - (_index + 1) < 0:
                oe_vals_list[_index] /= (len(oe_vals_list) - _index)
                cg_percent_list[_index] /= (len(cg_percent_list) - _index)

            else:
                oe_vals_list[_index] /= window_size
                cg_percent_list[_index] /= window_size

        indices = find_islands(cg_percent_list, oe_vals_list)
        cpg_features = [SeqFeature(location=FeatureLocation(start, end), type="CpG_island",
                                   qualifiers={'created_by': 'SeqBuddy'}) for (start, end) in indices]
        for feature in rec.features:
            cpg_features.append(feature)
        _rec = SeqRecord(map_cpg(_seq, indices), id=rec.id, name=rec.name, description=rec.description,
                         dbxrefs=rec.dbxrefs, features=cpg_features, annotations=rec.annotations,
                         letter_annotations=rec.letter_annotations)

        records.append(_rec)
        _output[rec.id] = indices
    _seqbuddy.records = records
    return _seqbuddy, _output


def find_pattern(_seqbuddy, _pattern):  # TODO ambiguous letters mode
    """
    Finds occurences of a pattern in a SeqBuddy object
    :param _seqbuddy: The SeqBuddy object to be searched
    :param _pattern: The regex pattern to search with
    :return: A tuple containing an annotated SeqBuddy object and a dictionary of matches dict[id]
    """
    # search through sequences for regex matches. For example, to find micro-RNAs
    _pattern = _pattern.upper()
    _output = OrderedDict()
    for _rec in _seqbuddy.records:
        indices = []
        matches = re.finditer(_pattern, str(_rec.seq).upper())
        for match in matches:
            indices.append(match.start())
            _rec.features.append(SeqFeature(location=FeatureLocation(start=match.start(), end=match.end()),
                                            type='match', qualifiers={'regex': _pattern, 'added_by': 'SeqBuddy'}))
        _output[_rec.id] = indices
    return _seqbuddy, _output


# TODO do string formatting in command line ui
def find_repeats(_seqbuddy, _columns=1):
    """
    Finds sequences with identical IDs or sequences
    :param _seqbuddy: The SeqBuddy object to be searched
    :param _columns: The number of columns to be output
    :return: A tuple containing the unique records, the repeat IDs, the repeat sequences, and the string output
    """
    _columns = 1 if _columns == 0 else abs(_columns)
    unique_seqs = {}
    repeat_ids = {}
    repeat_seqs = {}

    # First find replicate IDs
    # MD5 hash all sequences as we go for memory efficiency when looking for replicate sequences (below)
    # Need to work from a copy though, so sequences aren't overwritten
    _seqbuddy = deepcopy(_seqbuddy)
    for _rec in _seqbuddy.records:
        _seq = str(_rec.seq).encode()
        _seq = md5(_seq).hexdigest()
        _rec.seq = Seq(_seq)
        if _rec.id in repeat_ids:
            repeat_ids[_rec.id].append(_rec)
        elif _rec.id in unique_seqs:
            repeat_ids[_rec.id] = [_rec]
            repeat_ids[_rec.id].append(unique_seqs[_rec.id])
            del(unique_seqs[_rec.id])
        else:
            unique_seqs[_rec.id] = _rec

    # Then look for replicate sequences
    flip_uniqe = {}
    del_keys = []
    for _key, _value in unique_seqs.items():  # find and remove duplicates in/from the unique list
        _value = str(_value.seq)
        if _value not in flip_uniqe:
            flip_uniqe[_value] = [_key]
        else:
            if _value not in repeat_seqs:
                repeat_seqs[_value] = [_key]
                repeat_seqs[_value] += flip_uniqe[_value]
                if flip_uniqe[_value][0] in unique_seqs:
                    del_keys.append(flip_uniqe[_value][0])
            else:
                repeat_seqs[_value].append(_key)
            del_keys.append(unique_seqs[_key].id)

    for _key in del_keys:
        if _key in unique_seqs:
            del(unique_seqs[_key])

    for _key, _value in repeat_ids.items():  # find duplicates in the repeat ID list
        for _rep_seq in _value:
            _rep_seq = str(_rep_seq.seq)
            if _rep_seq not in flip_uniqe:
                flip_uniqe[_rep_seq] = [_key]
            else:
                if _rep_seq not in repeat_seqs:
                    repeat_seqs[_rep_seq] = [_key]
                    repeat_seqs[_rep_seq] += flip_uniqe[_rep_seq]

                else:
                    repeat_seqs[_rep_seq].append(_key)

    output_str = ""
    if len(repeat_ids) > 0:
        output_str += "#### Records with duplicate IDs: ####\n"
        _counter = 1
        for _next_id in repeat_ids:
            output_str += "%s\t" % _next_id
            if _counter % _columns == 0:
                output_str = "%s\n" % output_str.strip()
            _counter += 1

        output_str = "%s\n\n" % output_str.strip()

    else:
        output_str += "#### No records with duplicate IDs ####\n\n"

    if len(repeat_seqs) > 0:
        output_str += "#### Records with duplicate sequences: ####\n"
        _counter = 1
        for _next_id in repeat_seqs:
            output_str += "["
            for seq_id in repeat_seqs[_next_id]:
                output_str += "%s, " % seq_id
            output_str = "%s], " % output_str.strip(", ")

            if _counter % _columns == 0:
                output_str = "%s\n" % output_str.strip(", ")

            _counter += 1

        output_str = "%s\n\n" % output_str.strip(", ")
    else:
        output_str += "#### No records with duplicate sequences ####\n\n"

    output_str = "{0}\n".format(output_str.strip())
    return [unique_seqs, repeat_ids, repeat_seqs, output_str]


def find_restriction_sites(_seqbuddy, _enzymes="commercial", _min_cuts=1, _max_cuts=None):  # ToDo: Make sure cut sites are not already in the features list
    """
    Finds the restriction sites in the sequences in the SeqBuddy object
    :param _seqbuddy: The SeqBuddy object to be analyzed
    :param _enzymes: "commercial", "all", or a list of specific enzyme names
    :param _min_cuts: The minimum cut threshold
    :param _max_cuts: The maximum cut threshold
    :return: Returns a tuple containing an annotated SeqBuddy object, and a dictionary of restriction sites dict[id][re]
    """
    if _seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Unable to identify restriction sites in protein sequences.")
    if _max_cuts and _min_cuts > _max_cuts:
        raise ValueError("min_cuts parameter has been set higher than max_cuts.")
    _max_cuts = 1000000000 if not _max_cuts else _max_cuts

    _enzymes = _enzymes if type(_enzymes) == list else [_enzymes]

    blacklist = ["AbaSI", "FspEI", "MspJI", "SgeI", "AspBHI", "SgrTI", "YkrI", "BmeDI"]  # highly nonspecific
    blacklist += ["AjuI", "AlfI", "AloI", "ArsI", "BaeI", "BarI", "BcgI", "BdaI", "BplI", "BsaXI", "Bsp24I", "CjeI",
                  "CjePI", "CspCI", "FalI", "Hin4I", "NgoAVIII", "NmeDI", "PpiI", "PsrI", "R2_BceSIV", "RdeGBIII",
                  "SdeOSI", "TstI", "UcoMSI"]  # two-cutting
    blacklist += ["AlwFI", "AvaIII", "BmgI", "BscGI", "BspGI", "BspNCI", "Cdi630V", "Cgl13032I", "Cgl13032II",
                  "CjeFIII", "CjeFV", "CjeNII", "CjeP659IV", "CjuI", "CjuII", "DrdII", "EsaSSI", "FinI", "GauT27I",
                  "HgiEII", "Hpy99XIII", "Hpy99XIV", "Jma19592I", "MjaIV", "MkaDII", "NhaXI", "PenI", "Pfl1108I",
                  "RdeGBI", "RflFIII", "RlaI", "RpaTI", "SnaI", "Sno506I", "SpoDI", "TssI", "TsuI", "UbaF11I",
                  "UbaF12I", "UbaF13I", "UbaF14I", "UbaF9I", "UbaPI"]  # non-cutters

    batch = RestrictionBatch([])
    for enzyme in _enzymes:
        if enzyme == "commercial":
            for res in CommOnly:
                if str(res) not in blacklist:
                    batch.add(res)

        elif enzyme == "all":
            for res in AllEnzymes:
                if str(res) not in blacklist:
                    batch.add(res)

        else:
            try:
                batch.add(enzyme)
            except ValueError:
                _stderr("Warning: %s not a known enzyme\n" % enzyme)

    sites = []
    for _rec in _seqbuddy.records:
        _rec.res_sites = {}
        analysis = Analysis(batch, _rec.seq)
        result = analysis.with_sites()
        for _key, _value in result.items():
            if _key.cut_twice():
                _stderr("Warning: Double-cutters not supported.\n")
                pass
            elif _min_cuts <= len(_value) <= _max_cuts:
                try:
                    for zyme in _value:
                        cut_start = zyme + _key.fst3 - 1
                        cut_end = zyme + _key.fst5 + abs(_key.ovhg) - 1
                        _rec.features.append(SeqFeature(FeatureLocation(start=cut_start, end=cut_end), type=str(_key)))
                except TypeError:
                    _stderr("Warning: No-cutters not supported.\n")
                    pass
                _rec.res_sites[_key] = _value
        _rec.res_sites = OrderedDict(sorted(_rec.res_sites.items(), key=lambda x: x[0]))
        sites.append((_rec.id, _rec.res_sites))
    order_features_alphabetically(_seqbuddy)

    return [_seqbuddy, sites]


# TODO do string formatting in command line ui
def hash_sequence_ids(_seqbuddy, _hash_length=10):
    """
    Replaces the sequence IDs with random hashes
    :param _seqbuddy: The SeqBuddy to be hashed
    :param _hash_length: Specifies the length of the random hashes
    :return: A tuple containing: the hashed SeqBuddy, a dictionary mapping hashes to IDs, and a string representation
    """
    hash_list = []
    seq_ids = []
    if type(_hash_length) != int or _hash_length < 1:
        sys.stderr.write("Warning: The _hash_length parameter was passed in with the value %s. This is not an integer"
                         " greater than 0, so the hash length as been set to 10.\n\n" % _hash_length)
        _hash_length = 10

    if 32 ** _hash_length <= len(_seqbuddy.records) * 2:
        holder = ceil(log(len(_seqbuddy.records) * 2, 32))
        sys.stderr.write("Warning: The _hash_length parameter was passed in with the value %s. This is too small to "
                         "properly cover all sequences, so it has been increased to %s.\n\n" % (_hash_length, holder))
        _hash_length = holder

    for i in range(len(_seqbuddy.records)):
        new_hash = ""
        seq_ids.append(_seqbuddy.records[i].id)
        while True:
            new_hash = "".join([choice(string.ascii_letters + string.digits) for _ in range(_hash_length)])
            if new_hash in hash_list:
                continue
            else:
                hash_list.append(new_hash)
                break
        _seqbuddy.records[i].id = new_hash
        _seqbuddy.records[i].name = new_hash

    _hash_map = OrderedDict()
    for i in range(len(hash_list)):
        _hash_map[hash_list[i]] = seq_ids[i]

    _hash_table = "# Hash table\n"
    for _seq in _hash_map:
        _hash_table += "%s,%s\n" % (_seq, _hash_map[_seq])

    return [_seqbuddy, _hash_map, _hash_table]


def insert_sequence(_seqbuddy, _sequence, _location):
    """
    Add a specific sequence at a defined location in all records. E.g., adding a barcode (zero-indexed)
    :param _seqbuddy: The SeqBuddy object to be modified
    :param _sequence: The sequence to be inserted
    :param _location: The location to insert the sequences at
    :return: The modified SeqBuddy object
    """
    if type(_location) is int:
        for _rec in _seqbuddy.records:
            new_seq = _rec.seq[:_location] + _sequence + _rec.seq[_location:]
            _rec.seq = new_seq
    elif 'start' in _location:
        for _rec in _seqbuddy.records:
            new_seq = _sequence + _rec.seq
            _rec.seq = new_seq
    elif 'end' in _location:
        for _rec in _seqbuddy.records:
            new_seq = _rec.seq + _sequence
            _rec.seq = new_seq
    else:
        raise TypeError("Location must be 'start', 'end', or int.")
    return _seqbuddy


def isoelectric_point(_seqbuddy):
    """
    Calculate the isoelectric point of each sequence
    :param _seqbuddy: The SeqBuddy object to be analyzed
    :return: A tuple containing an annotated SeqBuddy object and a dictionary of isoelectric point values - dict[id]
    """
    if _seqbuddy.alpha is not IUPAC.protein:
        raise TypeError("Protein sequence required, not nucleic acid.")
    _isoelectric_points = OrderedDict()
    for _rec in _seqbuddy.records:
        _pI = ProteinAnalysis(str(_rec.seq))
        _pI = round(_pI.isoelectric_point(), 10)
        _isoelectric_points[_rec.id] = _pI
        _rec.features.append(SeqFeature(location=FeatureLocation(start=1, end=len(_rec.seq)), type='pI',
                                        qualifiers={'value': _pI}))
    return _seqbuddy, _isoelectric_points


def list_features(_seqbuddy):
    """
    Returns a dictionary of sequence features
    :param _seqbuddy: The SeqBuddy object to be analyzed
    :return: A dictionary of sequence features - dict[id] = list(features)
    """
    _output = OrderedDict()
    for _rec in _seqbuddy.records:
        _output[_rec.id] = _rec.features if len(_rec.features) > 0 else None
    return _output


def list_ids(_seqbuddy, _columns=1):  # TODO Make this return a list
    """
    Returns a list of sequence IDs
    :param _seqbuddy: The SeqBuddy object to be analyzed
    :return: A string listing sequence IDs
    """
    _columns = 1 if _columns == 0 else abs(_columns)
    _output = ""
    _counter = 1
    for rec in _seqbuddy.records:
        _output += "%s\t" % rec.id
        if _counter % _columns == 0:
            _output = "%s\n" % _output.strip()
        _counter += 1
    return "%s\n" % _output.strip()


def lowercase(_seqbuddy):
    """
    Converts all sequence residues to lowercase.
    :param _seqbuddy: The SeqBuddy object to be modified.
    :return: The modified SeqBuddy object
    """
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq).lower(), alphabet=_rec.seq.alphabet)
    return _seqbuddy


def map_features_dna2prot(dna_seqbuddy, prot_seqbuddy):
    """
    Applies DNA features to protein sequences
    :param dna_seqbuddy: A DNA SeqBuddy object with features to map
    :param prot_seqbuddy: A protein SeqBuddy
    :return: A protein SeqBuddy with the DNA SeqBuddy's features
    """
    def _feature_map(_feature):
        if type(_feature.location) == CompoundLocation:
            new_compound_location = []
            for sub_feature in _feature.location.parts:
                sub_feature = _feature_map(SeqFeature(sub_feature))
                new_compound_location.append(sub_feature.location)
            _feature.location = CompoundLocation(new_compound_location, _feature.location.operator)

        elif type(_feature.location) == FeatureLocation:
            _start = _feature.location.start / 3
            _end = _feature.location.end / 3
            _location = FeatureLocation(floor(_start), floor(_end))
            _feature.location = _location

        else:
            raise TypeError("_feature_map requires a feature with either FeatureLocation or CompoundLocation, "
                            "not %s" % type(_feature.location))
        return _feature

    prot_seqbuddy = clean_seq(prot_seqbuddy, "*")
    dna_seqbuddy = clean_seq(dna_seqbuddy)
    prot_dict = SeqIO.to_dict(prot_seqbuddy.records)
    dna_dict = SeqIO.to_dict(dna_seqbuddy.records)
    _new_seqs = {}
    stderr_written = False
    for _seq_id, dna_rec in dna_dict.items():
        if _seq_id not in prot_dict:
            stderr_written = True
            sys.stderr.write("Warning: %s is in the cDNA file, but not in the protein file\n" % _seq_id)
            continue

        if len(prot_dict[_seq_id].seq) * 3 not in [len(dna_rec.seq), len(dna_rec.seq) - 3]:  # len(cds) or len(cds minus stop)
            sys.stderr.write("Warning: size mismatch between aa and nucl seqs for %s --> %s, %s\n" %
                             (_seq_id, len(dna_rec.seq), len(prot_dict[_seq_id].seq)))
        _new_seqs[_seq_id] = prot_dict[_seq_id]
        prot_feature_hashes = []
        for feature in prot_dict[_seq_id].features:
            prot_feature_hashes.append(md5(str(feature).encode()).hexdigest())

        for feature in dna_rec.features:
            feature = _feature_map(feature)
            if md5(str(feature).encode()).hexdigest() not in prot_feature_hashes:
                prot_dict[_seq_id].features.append(feature)

    for _seq_id, prot_rec in prot_dict.items():
        if _seq_id not in dna_dict:
            stderr_written = True
            sys.stderr.write("Warning: %s is in the protein file, but not in the cDNA file\n" % _seq_id)
            _new_seqs[_seq_id] = prot_rec

    if stderr_written:
        sys.stderr.write("\n")

    _seqs_list = [_new_seqs[_rec.id] for _rec in prot_seqbuddy.records]
    _seqbuddy = SeqBuddy(_seqs_list)
    _seqbuddy.out_format = "gb"
    return _seqbuddy


def map_features_prot2dna(prot_seqbuddy, dna_seqbuddy):
    """
    Applies protein features to DNA sequences
    :param prot_seqbuddy: A protein SeqBuddy object with features to map
    :param dna_seqbuddy: A DNA SeqBuddy
    :return: A DNA SeqBuddy with the protein SeqBuddy's features
    """
    def _feature_map(_feature):
        if type(_feature.location) == CompoundLocation:
            new_compound_location = []
            for sub_feature in _feature.location.parts:
                sub_feature = _feature_map(SeqFeature(sub_feature))
                new_compound_location.append(sub_feature.location)
            _feature.location = CompoundLocation(new_compound_location, _feature.location.operator)

        elif type(_feature.location) == FeatureLocation:
            _start = feature.location.start * 3
            _end = feature.location.end * 3
            _location = FeatureLocation(_start, _end)
            _feature.location = _location

        else:
            raise TypeError("_feature_map requires a feature with either FeatureLocation or CompoundLocation, "
                            "not %s" % type(_feature.location))
        return _feature

    prot_seqbuddy = clean_seq(prot_seqbuddy, "*")
    dna_seqbuddy = clean_seq(dna_seqbuddy)
    prot_dict = SeqIO.to_dict(prot_seqbuddy.records)
    dna_dict = SeqIO.to_dict(dna_seqbuddy.records)
    _new_seqs = {}
    stderr_written = False
    for _seq_id, prot_rec in prot_dict.items():
        if _seq_id not in dna_dict:
            stderr_written = True
            sys.stderr.write("Warning: %s is in the protein file, but not in the cDNA file\n" % _seq_id)
            continue

        if len(prot_rec.seq) * 3 not in [len(dna_dict[_seq_id].seq), len(dna_dict[_seq_id].seq) - 3]:  # len(cds) or len(cds minus stop)
            sys.stderr.write("Warning: size mismatch between aa and nucl seqs for %s --> %s, %s\n" %
                             (_seq_id, len(prot_rec.seq), len(dna_dict[_seq_id].seq)))
        _new_seqs[_seq_id] = dna_dict[_seq_id]
        dna_feature_hashes = []
        for feature in dna_dict[_seq_id].features:
            dna_feature_hashes.append(md5(str(feature).encode()).hexdigest())

        for feature in prot_rec.features:
            feature = _feature_map(feature)
            prot_feature_hashes = [md5(str(feature).encode()).hexdigest()]
            # Need to account for strand orientation
            feature.strand = 0
            prot_feature_hashes.append(md5(str(feature).encode()).hexdigest())
            feature.strand = 1
            prot_feature_hashes.append(md5(str(feature).encode()).hexdigest())
            if not set(prot_feature_hashes) & set(dna_feature_hashes):
                dna_dict[_seq_id].features.append(feature)

    for _seq_id, dna_rec in dna_dict.items():
        if _seq_id not in prot_dict:
            stderr_written = True
            sys.stderr.write("Warning: %s is in the cDNA file, but not in the protein file\n" % _seq_id)
            _new_seqs[_seq_id] = dna_rec

    if stderr_written:
        sys.stderr.write("\n")

    _seqs_list = [_new_seqs[_rec.id] for _rec in dna_seqbuddy.records]
    _seqbuddy = SeqBuddy(_seqs_list)
    _seqbuddy.out_format = "gb"
    return _seqbuddy


def merge(_seqbuddy_list):
    """
    Combines two or more SeqBuddy objects
    :param _seqbuddy_list: The list of SeqBuddy objects to be merged
    :return: A single, merged SeqBuddy object
    """
    _output = _seqbuddy_list[0]
    for _seqbuddy in _seqbuddy_list[1:]:
        _output.records += _seqbuddy.records
    return _output


def molecular_weight(_seqbuddy):
    """
    Calculates the mass of each sequence in daltons
    :param _seqbuddy: The SeqBuddy object to be analyzed
    :return: A tuple containing an annotated SeqBuddy object and a dictionary of molecular weight values -
    dict[id][(ssRNA_value/ssDNA_value, dsDNA_value/peptide_value)]
    """

    amino_acid_weights = {'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09, 'C': 103.14, 'Q': 128.13, 'E': 129.12,
                          'G': 57.05, 'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17, 'M': 131.19, 'F': 147.18,
                          'P': 97.12, 'S': 87.08, 'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13, '-': 0, '*': 0,
                          'X': 110}
    deoxynucleotide_weights = {'A': 313.2, 'G': 329.2, 'C': 289.2, 'T': 304.2, 'Y': 296.7, 'R': 321.2, 'W': 308.7,
                               'S': 309.2, 'K': 316.7, 'M': 301.2, 'D': 315.53, 'V': 310.53, 'H': 302.2, 'B': 307.53,
                               'X': 308.95, 'N': 308.95, '-': 0, '.': 0}
    deoxyribonucleotide_weights = {'A': 329.2, 'G': 306.2, 'C': 305.2, 'U': 345.2, 'Y': 325.2, 'R': 317.7, 'W': 337.2,
                                   'S': 305.7, 'K': 325.7, 'M': 317.2, 'D': 326.87, 'V': 313.53, 'H': 326.53,
                                   'B': 318.87, 'X': 321.45, 'N': 321.45, '-': 0, '.': 0}
    deoxynucleotide_compliments = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A', 'Y': 'R', 'R': 'Y', 'W': 'W',
                                   'S': 'S', 'K': 'M', 'M': 'K', 'D': 'H', 'V': 'B', 'H': 'D', 'B': 'V',
                                   'X': 'X', 'N': 'N', '-': '-', '.': '.'}
    _dna = False
    _output = {'masses_ss': [], 'masses_ds': [], 'ids': []}
    _dict = amino_acid_weights
    if _seqbuddy.alpha == IUPAC.protein:
        _dict = amino_acid_weights
    elif _seqbuddy.alpha in [IUPAC.ambiguous_dna or IUPAC.unambiguous_dna]:
        _dict = deoxynucleotide_weights
        _dna = True
    elif _seqbuddy.alpha in [IUPAC.ambiguous_rna or IUPAC.unambiguous_rna]:
        _dict = deoxyribonucleotide_weights
    for _rec in _seqbuddy.records:
        _rec.mass_ds = 0
        _rec.mass_ss = 0
        if _seqbuddy.alpha == IUPAC.protein:
            _rec.mass_ss += 18.02  # molecular weight of a water molecule
        else:
            if _dna:
                _rec.mass_ss += 79.0  # molecular weight of 5' monophosphate in ssDNA
                _rec.mass_ds += 157.9  # molecular weight of the 5' triphosphate in dsDNA
            else:
                _rec.mass_ss += 159.0  # molecular weight of a 5' triphosphate in ssRNA
        for _indx, _value in enumerate(str(_rec.seq).upper()):
            _rec.mass_ss += _dict[_value]
            if _dna:
                _rec.mass_ds += _dict[_value] + deoxynucleotide_weights[deoxynucleotide_compliments[_value]]
        _output['masses_ss'].append(round(_rec.mass_ss, 3))

        _qualifiers = {}
        if _seqbuddy.alpha == IUPAC.protein:
            _qualifiers["peptide_value"] = round(_rec.mass_ss, 3)
        elif _dna:
            _qualifiers["ssDNA_value"] = round(_rec.mass_ss, 3)
            _qualifiers["dsDNA_value"] = round(_rec.mass_ds, 3)
            _output['masses_ds'].append(round(_rec.mass_ds, 3))
        elif _seqbuddy.alpha in [IUPAC.ambiguous_rna or IUPAC.unambiguous_rna]:
            _qualifiers["ssRNA_value"] = round(_rec.mass_ss, 3)
        _output['ids'].append(_rec.id)
        mw_feature = SeqFeature(location=FeatureLocation(start=1, end=len(_rec.seq)), type='mw', qualifiers=_qualifiers)
        _rec.features.append(mw_feature)
    return _seqbuddy, _output


def num_seqs(_seqbuddy):
    """
    Counts the number of sequences in the SeqBuddy object
    :param _seqbuddy: The SeqBuddy object to be counted
    :return: The int number of sequences
    """
    return len(_seqbuddy.records)


def order_features_alphabetically(_seqbuddy, _reverse=False):
    """
    Sorts features in alphabetical order
    :param _seqbuddy: The SeqBuddy object to have its features sorted
    :param _reverse: Specifies if the features should be sorted backwards
    :return: The SeqBuddy object with sorted features
    """
    for _rec in _seqbuddy.records:
        new_feature_list = [(_feature.type, _feature) for _feature in _rec.features]
        new_feature_list = sorted(new_feature_list, key=lambda x: x[0], reverse=_reverse)
        new_feature_list = [_feature[1] for _feature in new_feature_list]
        _rec.features = new_feature_list
    return _seqbuddy


def order_features_by_position(_seqbuddy, _reverse=False):
    """
    Sorts features by the order in which they appear in the sequence
    :param _seqbuddy: The SeqBuddy object to have its features sorted
    :param _reverse: Specifies if the features should be sorted backwards
    :return: The SeqBuddy object with sorted features
    """
    for _rec in _seqbuddy.records:
        new_feature_list = [(int(_feature.location.start), _feature) for _feature in _rec.features]
        new_feature_list = sorted(new_feature_list, key=lambda x: x[0], reverse=_reverse)
        new_feature_list = [_feature[1] for _feature in new_feature_list]
        _rec.features = new_feature_list
    return _seqbuddy


def order_ids(_seqbuddy, _reverse=False):
    """
    Sorts the sequences by ID, alphabetically
    :param _seqbuddy: The SeqBuddy object to be sorted
    :param _reverse: Reverses the sequence order
    :return: The sorted SeqBuddy object
    """
    _output = [(_rec.id, _rec) for _rec in _seqbuddy.records]
    _output = sorted(_output, key=lambda x: x[0], reverse=_reverse)
    _output = [_rec[1] for _rec in _output]
    _seqbuddy.records = _output
    return _seqbuddy


def order_ids_randomly(_seqbuddy):
    """
    Randomly reorders the sequences
    :param _seqbuddy: The SeqBuddy object to be shuffled
    :return: The shuffled SeqBuddy object
    """
    _output = []
    for _ in range(len(_seqbuddy.records)):
        random_index = randint(1, len(_seqbuddy.records)) - 1
        _output.append(_seqbuddy.records.pop(random_index))
    _seqbuddy.records = _output
    return _seqbuddy


def pull_random_recs(_seqbuddy, _count=1):  # Return a random set of sequences (without replacement)
    """
    Randomly retrieves sequences
    :param _seqbuddy: The SeqBuddy object to be pulled from
    :param _count: The number of records to pull
    :return: The modified SeqBuddy object
    """
    random_recs = []
    _count = abs(_count) if abs(_count) <= len(_seqbuddy.records) else len(_seqbuddy.records)
    for i in range(_count):
        rand_index = randint(0, len(_seqbuddy.records) - 1)
        random_recs.append(_seqbuddy.records.pop(rand_index))

    _seqbuddy.records = random_recs
    return _seqbuddy


def pull_record_ends(_seqbuddy, _amount, _which_end):
    """
    Retrieves subsequences from the ends of the sequences
    :param _seqbuddy: The SeqBuddy object to be pulled from
    :param _amount: The number of residues to be pulled
    :param _which_end: Which end to pull from (front/rear)
    :return: The modified SeqBuddy object
    """
    _amount = int(_amount)
    if _amount < 0:
        raise ValueError("Positive integer required for '_amount' argument in pull_record_ends.")

    seq_ends = []
    for _rec in _seqbuddy.records:
        if _which_end == 'front':
            _rec.seq = Seq(str(_rec.seq)[:_amount], alphabet=_rec.seq.alphabet)
            _rec.features = _shift_features(_rec.features, 0, len(str(_rec.seq)))

        elif _which_end == "rear":
            _shift = -1 * (len(str(_rec.seq)) - _amount)
            _rec.features = _shift_features(_rec.features, _shift, len(str(_rec.seq)))
            _rec.seq = _rec.seq[-1 * _amount:]

        else:
            raise AttributeError("You must pick 'front' or 'rear' for the '_which_end' argument in pull_record_ends.")

        seq_ends.append(_rec)

    _seqbuddy.records = seq_ends
    return _seqbuddy


def pull_recs(_seqbuddy, _search):  # _search can be a list of regex expressions or single string
    """
    Retrieves sequences with names/IDs matching a search pattern
    :param _seqbuddy: The SeqBuddy object to be pulled from
    :param _search: The regex pattern to search with
    :return: The modified SeqBuddy object
    """
    _search = "|".join(_search) if type(_search) == list else _search
    matched_records = []
    for _rec in _seqbuddy.records:
        if re.search(_search, _rec.description) or re.search(_search, _rec.id) or re.search(_search, _rec.name):
            matched_records.append(_rec)
    _seqbuddy.records = matched_records
    return _seqbuddy


def purge(_seqbuddy, threshold):  # ToDo: Implement a way to return a certain # of seqs (i.e. auto-determine threshold)
    """
    Deletes highly similar sequences
    :param _seqbuddy: The SeqBuddy object to be purged
    :param threshold: Sets the similarity threshold
    :return: The purged SeqBuddy object
    """
    keep_set = {}
    purged = []
    _blast_res = bl2seq(_seqbuddy)[0]
    _blast_res = [(_key, _value) for _key, _value in _blast_res.items()]
    _blast_res = sorted(_blast_res, key=lambda l: l[0])
    for _query_id, match_list in _blast_res:
        if _query_id in purged:
            continue
        else:
            keep_set[_query_id] = []
            for _subj_id in match_list:
                _ident, _length, _evalue, _bit_score = match_list[_subj_id]

                if _bit_score >= threshold:
                    purged.append(_subj_id)
                    keep_set[_query_id].append(_subj_id)

    _output = []
    for _rec in _seqbuddy.records:
        if _rec.id in keep_set:
            _output.append(_rec)

    _seqbuddy.records = _output
    # TODO do string formatting in command line ui
    _record_map = "### Deleted record mapping ###\n"
    keep_set = [(_key, sorted(_value)) for _key, _value in keep_set.items()]
    keep_set = sorted(keep_set, key=lambda l: l[0])
    for _seq_id, seq_list in keep_set:
        _record_map += "%s\n" % _seq_id
        for del_seq_id in seq_list:
            _record_map += "%s, " % del_seq_id
        _record_map = _record_map.strip(", ") + "\n\n"

    _record_map = _record_map.strip() + "\n##############################\n\n"

    return [_seqbuddy, purged, _record_map]


def raw_seq(_seqbuddy):  # TODO Make this return a dict
    """
    Returns the raw sequence data from the SeqBuddy object
    :param _seqbuddy: The SeqBuddy object to be analyzed
    :return: A string containing the raw sequences
    """
    _seqbuddy = clean_seq(_seqbuddy)
    _output = ""
    for _rec in _seqbuddy.records:
        _output += "%s\n\n" % _rec.seq

    return "%s\n" % _output.strip()


def rename(_seqbuddy, query, replace="", _num=0):  # TODO Allow a replacement pattern increment (like numbers)
    """
    Rename sequence IDs
    :param _seqbuddy: The SeqBuddy object to be modified
    :param query: The pattern to be searched for
    :param replace: The string to be substituted
    :param _num: The maximum number of substitutions to make
    :return: The modified SeqBuddy object
    """
    for _rec in _seqbuddy.records:
        new_name = re.sub(query, replace, _rec.id, _num)
        _rec.id = new_name
        _rec.name = new_name
    return _seqbuddy


def reverse_complement(_seqbuddy):
    """
    Converts DNA/RNA sequences to their reverse complementary sequence
    :param _seqbuddy: The SeqBuddy object to be modified
    :return: The modified SeqBuddy object
    """
    if _seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for _rec in _seqbuddy.records:
        _rec.seq = _rec.seq.reverse_complement()
        seq_len = len(_rec.seq)
        shifted_features = [_feature_rc(_feature, seq_len) for _feature in _rec.features]
        _rec.features = shifted_features
    return _seqbuddy


def rna2dna(_seqbuddy):
    """
    Transcribes RNA into DNA sequences.
    :param _seqbuddy: The SeqBuddy object to be transcribed
    :return: The transcribed SeqBuddy object
    """
    if _seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq.back_transcribe()), alphabet=IUPAC.ambiguous_dna)
    _seqbuddy.alpha = IUPAC.ambiguous_dna
    return _seqbuddy


def select_frame(_seqbuddy, frame):  # ToDo: record the deleted residues so the earlier frame can be returned to.
    """
    Changes the reading frame of the sequences
    :param _seqbuddy: The SeqBuddy object to be shifted
    :param frame: The reading frame to shift to
    :return: The shifted SeqBuddy object
    """
    if _seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Select frame requires nucleic acid, not protein.")
    for _rec in _seqbuddy.records:
        _rec.features = _shift_features(_rec.features, (frame - 1) * -1, len(_rec.seq))
        _rec.seq = Seq(str(_rec.seq)[frame - 1:], alphabet=_rec.seq.alphabet)
    return _seqbuddy


def shuffle_seqs(_seqbuddy):
    """
    Randomly reorder the residues in each sequence
    :param _seqbuddy: The SeqBuddy object to be shuffled
    :return: The shuffled SeqBuddy object
    """
    for _rec in _seqbuddy.records:
        tokens = []
        for letter in _rec.seq:
            tokens.append(letter)
        new_seq = ''
        while len(tokens) > 0:
            rand_indx = randint(0, len(tokens) - 1)
            new_seq += tokens.pop(rand_indx)
        _rec.seq = Seq(data=new_seq, alphabet=_seqbuddy.alpha)
    return _seqbuddy


def split_by_taxa(_seqbuddy, split_pattern):
    """
    Splits a SeqBuddy object by a specified pattern
    :param _seqbuddy: The SeqBuddy object to be split
    :param split_pattern: The regex pattern to split with
    :return: A dictionary of SeqRecords
    """
    recs_by_taxa = {}
    for _rec in _seqbuddy.records:
        split = re.split(split_pattern, _rec.id)
        recs_by_taxa.setdefault(split[1], []).append(_rec)
    return recs_by_taxa


def split_file(_seqbuddy):
    """
    Split the records in a SeqBuddy object up into a collection of new SeqBuddy objects
    :param _seqbuddy: The SeqBuddy object to be split
    :return: A list of SeqBuddy objects
    """
    sb_objs_list = []
    for _rec in _seqbuddy.records:
        _sb = SeqBuddy([_rec])
        _sb.in_format = _seqbuddy.in_format
        _sb.out_format = _seqbuddy.out_format
        sb_objs_list.append(_sb)
    return sb_objs_list


def translate6frames(_seqbuddy):
    """
    Translates a nucleotide sequence into a protein sequence across all six reading frames.
    :param _seqbuddy: The SeqBuddy object to be translated
    :return: The translated SeqBuddy object
    """
    frame1, frame2, frame3 = deepcopy(_seqbuddy), deepcopy(_seqbuddy), deepcopy(_seqbuddy)
    _seqbuddy = reverse_complement(_seqbuddy)

    rframe1, rframe2, rframe3 = deepcopy(_seqbuddy), deepcopy(_seqbuddy), deepcopy(_seqbuddy)

    frame2 = select_frame(frame2, 2)
    frame3 = select_frame(frame3, 3)
    rframe2 = select_frame(rframe2, 2)
    rframe3 = select_frame(rframe3, 3)

    frame1 = translate_cds(frame1, quiet=True)
    frame2 = translate_cds(frame2, quiet=True)
    frame3 = translate_cds(frame3, quiet=True)
    rframe1 = translate_cds(rframe1, quiet=True)
    rframe2 = translate_cds(rframe2, quiet=True)
    rframe3 = translate_cds(rframe3, quiet=True)

    _output = []
    for _i in range(len(frame1.records)):
        frame1.records[_i].id = "%s_f1" % frame1.records[_i].id
        frame2.records[_i].id = "%s_f2" % frame2.records[_i].id
        frame3.records[_i].id = "%s_f3" % frame3.records[_i].id
        rframe1.records[_i].id = "%s_rf1" % rframe1.records[_i].id
        rframe2.records[_i].id = "%s_rf2" % rframe2.records[_i].id
        rframe3.records[_i].id = "%s_rf3" % rframe3.records[_i].id

        _output += [frame1.records[_i], frame2.records[_i], frame3.records[_i],
                    rframe1.records[_i], rframe2.records[_i], rframe3.records[_i]]

    _seqbuddy = SeqBuddy(_output, _out_format=_seqbuddy.out_format)
    return _seqbuddy


# ToDo: Deal with alignments...
def translate_cds(_seqbuddy, quiet=False):  # adding 'quiet' will suppress the errors thrown by translate(cds=True)
    """
    Translates a nucleotide sequence into a protein sequence.
    :param _seqbuddy: The SeqBuddy object to be translated
    :param quiet: Suppress error messages and warnings
    :return: The translated SeqBuddy object
    """
    def trans(in_seq):
        try:
            in_seq.seq = in_seq.seq.translate(cds=True, to_stop=True)
            return in_seq

        except TranslationError as _e1:
            if not quiet:
                sys.stderr.write("Warning: %s in %s\n" % (_e1, in_seq.id))
            return _e1

        except ValueError:
            raise TypeError("Nucleic acid sequence required, not protein.")

    _translation = deepcopy(_seqbuddy)
    _translation.alpha = IUPAC.protein
    for _rec in _translation.records:
        _rec.features = []
        temp_seq = deepcopy(_rec)
        while True:  # Modify a copy of the sequence as needed to complete the cds translation
            test_trans = trans(temp_seq)
            # success
            if str(type(test_trans)) == "<class 'Bio.SeqRecord.SeqRecord'>":
                break

            # not standard length
            if re.search("Sequence length [0-9]+ is not a multiple of three", str(test_trans)):
                temp_seq.seq = Seq(str(temp_seq.seq)[:(len(str(temp_seq.seq)) - len(str(temp_seq.seq)) % 3)],
                                   alphabet=temp_seq.seq.alphabet)
                _rec.seq = Seq(str(_rec.seq)[:(len(str(_rec.seq)) - len(str(_rec.seq)) % 3)],
                               alphabet=_rec.seq.alphabet)
                continue

            # not a start codon
            if re.search("First codon '[A-Za-z]{3}' is not a start codon", str(test_trans)):
                temp_seq.seq = Seq("ATG" + str(temp_seq.seq)[3:], alphabet=temp_seq.seq.alphabet)
                continue

            # not a stop codon
            if re.search("Final codon '[A-Za-z]{3}' is not a stop codon", str(test_trans)):
                temp_seq.seq = Seq(str(temp_seq.seq) + "TGA", alphabet=temp_seq.seq.alphabet)
                continue

            # non-standard characters
            if re.search("Codon '[A-Za-z]{3}' is invalid", str(test_trans)):
                regex = re.findall("Codon '([A-Za-z]{3})' is invalid", str(test_trans))
                regex = "(?i)%s" % regex[0]
                temp_seq.seq = Seq(re.sub(regex, "NNN", str(temp_seq.seq), count=1), alphabet=temp_seq.seq.alphabet)
                _rec.seq = Seq(re.sub(regex, "NNN", str(_rec.seq), count=1), alphabet=_rec.seq.alphabet)
                continue

            # internal stop codon(s) found
            if re.search("Extra in frame stop codon found", str(test_trans)):
                for _i in range(round(len(str(temp_seq.seq)) / 3) - 1):
                    _codon = str(temp_seq.seq)[(_i * 3):(_i * 3 + 3)]
                    if _codon.upper() in ["TGA", "TAG", "TAA"]:
                        new_seq = str(temp_seq.seq)[:(_i * 3)] + "NNN" + str(temp_seq.seq)[(_i * 3 + 3):]
                        temp_seq.seq = Seq(new_seq, alphabet=temp_seq.seq.alphabet)
                continue

            break

        try:
            _rec.seq = _rec.seq.translate()
            _rec.seq.alphabet = IUPAC.protein

        except TranslationError as e1:
            raise TranslationError("%s failed to translate  --> %s\n" % (_rec.id, e1))

    _output = map_features_dna2prot(_seqbuddy, _translation)
    _output.out_format = _seqbuddy.out_format
    _seqbuddy = _output
    return _seqbuddy


def uppercase(_seqbuddy):
    """
    Converts all sequence residues to uppercase.
    :param _seqbuddy: The SeqBuddy object to be modified.
    :return: The modified SeqBuddy object
    """
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq).upper(), alphabet=_rec.seq.alphabet)
    return _seqbuddy

def degenerate_sequence(_seqbuddy, table=1, reading_frame =1 ): 
    """
    Generate degenerate codon sequence 
    :param _seqbuddy: The SeqBuddy object to be analyzed
    :param _table: The degenerate codon table to use
    :param _reading_frame: selects reading frame to start creating degenerate codon
    :return: A SeqBuddy object containing a degenerate nucleotide seqeuence
   
    The method is developed based on the Perl script Degen v1.4.
    http://www.phylotools.com/ptdegenoverview.htm
    
    Zwick, A., Regier, J.C. & Zwickl, D.J. (2012). "Resolving Discrepancy between Nucleotides and Amino Acids in 
    Deep-Level Arthropod Phylogenomics: Differentiating Serine Codons in 21-Amino-Acid Models". PLoS ONE 7(11): e47450.

    Regier, J.C., Shultz, J.W., Zwick, A., Hussey, A., Ball, B., Wetzer, R. Martin, J.W. & Cunningham, C.W. (2010).
     "Arthropod relationships revealed by phylogenomic analysis of nuclear protein-coding sequences". Nature 463: 1079-1083.
    """
    #degenerate codons that are  universal to each all codon tables
    base_dict = {'ACA': 'ACN','ACC': 'ACN','ACG': 'ACN','ACT': 'ACN','CAC': 'CAY','CAT': 'CAY','CCA': 'CCN','CCC': 'CCN', 
                 'CCG': 'CCN','CCT': 'CCN','GAA': 'GAR','GAC': 'GAY','GAG': 'GAR','GAT': 'GAY','GCA': 'GCN','GCC': 'GCN',
                 'GCG': 'GCN','GCT': 'GCN','GTA': 'GTN','GTC': 'GTN','GTG': 'GTN','GTT': 'GTN','TCA': 'TCN','TCC': 'TCN', 
                 'TCG': 'TCN','TCT': 'TCN','TTC': 'TTY','TTT': 'TTY','NNN': 'NNN','???': 'NNN','---': '---'} 
    
    #Standard Genetic Code codons.  It is the default table
    dgn_dict_1 = {'AAA': 'AAR','AAC': 'AAY','AAG': 'AAR','AAT': 'AAY','AGA': 'MGN','AGC': 'TCN','AGG': 'MGN','AGT': 'AGY',
                  'ATA': 'ATH','ATC': 'ATH','ATG': 'ATG','ATT': 'ATH','CAA': 'CAR','CAG': 'CAR','CGA': 'MGN','CGC': 'MGN',
                  'CGG': 'MGN','CGT': 'MGN','CTA': 'YTN','CTC': 'YTN','CTG': 'YTN','CTT': 'YTN','GGA': 'GGN','GGC': 'GGN',
                  'GGG': 'GGN','GGT': 'GGN','TAA': 'TAA','TAC': 'TAY','TAG': 'TAG','TAT': 'TAY','TGA': 'TGA','TGC': 'TGY',
                  'TGG': 'TGG','TGT': 'TGY','TTA': 'YTN','TTG': 'YTN'}
    
    #Vertebrate Mitochondiral Code codons 
    dgn_dict_2 = {'AAA': 'AAR','AAC': 'AAY','AAG': 'AAR','AAT': 'AAY','ATA': 'ATR','ATC': 'ATY','ATG': 'ATR','ATT': 'ATY',
                  'AGA': 'AGA','AGG': 'AGG','AGC': 'AGY','AGT': 'AGY','CAA': 'CAR','CAG': 'CAR','CGA': 'CGN','CGC': 'CGN',
                  'CGG': 'CGN','CGT': 'CGN','CTA': 'YTN','CTC': 'YTN','CTG': 'YTN','CTT': 'YTN','GGA': 'GGN','GGC': 'GGN',
                  'GGG': 'GGN','GGT': 'GGN','TAA': 'TAA','TAC': 'TAY','TAG': 'TAG','TAT': 'TAY','TGA': 'TGR','TGC': 'TGY',
                  'TGG': 'TGR','TGT': 'TGY','TTA': 'YTN','TTG': 'YTN'}
    
    #Yeast Mitochondiral Code codons
    dgn_dict_3 = {'AAA': 'AAR','AAC': 'AAY','AAG': 'AAR','AAT': 'AAY','AGA': 'MGN','AGC': 'AGY','AGG': 'MGN','AGT': 'AGY',
                  'ATA': 'ATR','ATC': 'ATY','ATG': 'ATR','ATT': 'ATY','CAA': 'CAR','CAG': 'CAR','CGA': 'CGA','CGC': 'CGC',
                  'CGG': 'MGN','CGT': 'MGN','CTA': 'CTN','CTC': 'CTN','CTG': 'CTN','CTT': 'CTN','GGA': 'GGN','GGC': 'GGN',
                  'GGG': 'GGN','GGT': 'GGN','TAA': 'TAA','TAC': 'TAY','TAG': 'TAG','TAT': 'TAY','TGA': 'TGR','TGC': 'TGY',
                  'TGG': 'TGR','TGT': 'TGY','TTA': 'TTR','TTG': 'TTR'}
    
    #Mold / Protozoan / Coelenterate Mitochondrial Code & Mycoplasma / Spiroplasma Code Codons
    dgn_dict_4 = {'AAA': 'AAR','AAC': 'AAY','AAG': 'AAR','AAT': 'AAY','AGA': 'MGN','AGC': 'AGY','AGG': 'MGN','AGT': 'AGY',
                  'ATA': 'ATH','ATC': 'ATH','ATG': 'ATG','ATT': 'ATH','CAA': 'CAR','CAG': 'CAR','CGA': 'MGN','CGC': 'MGN',
                  'CGG': 'MGN','CGT': 'MGN','CTA': 'YTN','CTC': 'YTN','CTG': 'YTN','CTT': 'YTN','GGA': 'GGN','GGC': 'GGN',
                  'GGG': 'GGN','GGT': 'GGN','TAA': 'TAA','TAC': 'TAY','TAG': 'TAG','TAT': 'TAY','TGA': 'TGR','TGC': 'TGY',
                  'TGG': 'TGR','TGT': 'TGY','TTA': 'YTN','TTG': 'YTN'}
    
    #Invertebrate Mitochondrial Code codons
    dgn_dict_5 = {'AAA': 'AAR','AAC': 'AAY','AAG': 'AAR','AAT': 'AAY','AGA': 'AGN','AGC': 'AGN','AGG': 'AGN','AGT': 'AGN',
                  'ATA': 'ATR','ATC': 'ATY','ATG': 'ATR','ATT': 'ATY','CAA': 'CAR','CAG': 'CAR','CGA': 'CGN','CGC': 'CGN',
                  'CGG': 'CGN','CGT': 'CGN','CTA': 'YTN','CTC': 'YTN','CTG': 'YTN','CTT': 'YTN','GGA': 'GGN','GGC': 'GGN',
                  'GGG': 'GGN','GGT': 'GGN','TAA': 'TAA','TAC': 'TAY','TAG': 'TAG','TAT': 'TAY','TGA': 'TGR','TGC': 'TGY',
                  'TGG': 'TGR','TGT': 'TGY','TTA': 'YTN','TTG': 'YTN'}
    
    #Ciliate / Dasycladacean / Hexamita Nuclear Code codons
    dgn_dict_6 = {'AAA': 'AAR','AAC': 'AAY','AAG': 'AAR','AAT': 'AAY','AGA': 'MGN','AGC': 'AGY','AGG': 'MGN','AGT': 'AGY',
                  'ATA': 'ATH','ATC': 'ATH','ATG': 'ATG','ATT': 'ATH','CAA': 'YAR','CAG': 'YAR','CGA': 'MGN','CGC': 'MGN',
                  'CGG': 'MGN','CGT': 'MGN','CTA': 'YTN','CTC': 'YTN','CTG': 'YTN','CTT': 'YTN','GGA': 'GGN','GGC': 'GGN',
                  'GGG': 'GGN','GGT': 'GGN','TAA': 'YAR','TAC': 'TAY','TAG': 'YAR','TAT': 'TAY','TGA': 'TGA','TGC': 'TGY',
                  'TGG': 'TGG','TGT': 'TGY','TTA': 'YTN','TTG': 'YTN'}
    
    #Echinoderm and Flatworm Mitochondrial Code codons
    dgn_dict_9 = {'AAA': 'AAH','AAC': 'AAH','AAG': 'AAG','AAT': 'AAH','AGA': 'AGN','AGC': 'AGN','AGG': 'AGN','AGT': 'AGN',
                  'ATA': 'ATH','ATC': 'ATH','ATG': 'ATG','ATT': 'ATH','CAA': 'CAR','CAG': 'CAR','CGA': 'CGN','CGC': 'CGN',
                  'CGG': 'CGN','CGT': 'CGN','CTA': 'YTN','CTC': 'YTN','CTG': 'YTN','CTT': 'YTN','GGA': 'GGN','GGC': 'GGN',
                  'GGG': 'GGN','GGT': 'GGN','TAA': 'TAA','TAC': 'TAY','TAG': 'TAG','TAT': 'TAY','TGA': 'TGR','TGC': 'TGY',
                  'TGG': 'TGR','TGT': 'TGY','TTA': 'YTN','TTG': 'YTN'}
    
    #Euplotid Nuclear Code codons
    dgn_dict_10 = {'AAA': 'AAR','AAC': 'AAY','AAG': 'AAR','AAT': 'AAY','AGA': 'MGN','AGC': 'AGY','AGG': 'MGN','AGT': 'AGY',
                   'ATA': 'ATH','ATC': 'ATH','ATG': 'ATG','ATT': 'ATH','CAA': 'CAR','CAG': 'CAR','CGA': 'MGN','CGC': 'MGN',
                   'CGG': 'MGN','CGT': 'MGN','CTA': 'YTN','CTC': 'YTN','CTG': 'YTN','CTT': 'YTN','GGA': 'GGN','GGC': 'GGN',
                   'GGG': 'GGN','GGT': 'GGN','TAA': 'TAA','TAC': 'TAY','TAG': 'TAG','TAT': 'TAY','TGA': 'TGH','TGC': 'TGH',
                   'TGG': 'TGG','TGT': 'TGH','TTA': 'YTN','TTG': 'YTN'}            
    
    #Bacterial / Archaeal / Plant Plastid Code condons
    dgn_dict_11 = {'AAA': 'AAR','AAC': 'AAY','AAG': 'AAR','AAT': 'AAY','AGA': 'MGN','AGC': 'AGY','AGG': 'MGN','AGT': 'AGY',
                   'ATA': 'ATH','ATC': 'ATH','ATG': 'ATG','ATT': 'ATH','CAA': 'CAR','CAG': 'CAR','CGA': 'MGN','CGC': 'MGN',
                   'CGG': 'MGN','CGT': 'MGN','CTA': 'YTN','CTC': 'YTN','CTG': 'YTN','CTT': 'YTN','GGA': 'GGN','GGC': 'GGN',
                   'GGG': 'GGN','GGT': 'GGN','TAA': 'TAA','TAC': 'TAY','TAG': 'TAG','TAT': 'TAY','TGA': 'TGA','TGC': 'TGY',
                   'TGG': 'TGG','TGT': 'TGY','TTA': 'YTN','TTG': 'YTN'}
    
    #Alternative Yeast Nuclear Code codons
    dgn_dict_12 = {'AAA': 'AAR','AAC': 'AAY','AAG': 'AAR','AAT': 'AAY','AGA': 'MGN','AGC': 'AGY','AGG': 'MGN','AGT': 'AGY',
                   'ATA': 'ATH','ATC': 'ATH','ATG': 'ATG','ATT': 'ATH','CAA': 'CAR','CAG': 'CAR','CGA': 'MGN','CGC': 'MGN',
                   'CGG': 'MGN','CGT': 'MGN','CTA': 'YTN','CTC': 'YTN','CTG': 'CTG','CTT': 'YTN','GGA': 'GGN','GGC': 'GGN',
                   'GGG': 'GGN','GGT': 'GGN','TAA': 'TAA','TAC': 'TAY','TAG': 'TAG','TAT': 'TAY','TGA': 'TGA','TGC': 'TGY',
                   'TGG': 'TGG','TGT': 'TGY','TTA': 'YTN','TTG': 'YTN'}
    
    #Ascidian Mitochondrial Code codons
    dgn_dict_13 = {'AAA': 'AAR','AAC': 'AAY','AAG': 'AAR','AAT': 'AAY','AGA': 'RGN','AGC': 'AGY','AGG': 'RGN','AGT': 'AGY',
                   'ATA': 'ATR','ATC': 'ATY','ATG': 'ATR','ATT': 'ATY','CAA': 'CAR','CAG': 'CAR','CGA': 'CGN','CGC': 'CGN',
                   'CGG': 'CGN','CGT': 'CGN','CTA': 'YTN','CTC': 'YTN','CTG': 'YTN','CTT': 'YTN','GGA': 'RGN','GGC': 'RGN',
                   'GGG': 'RGN','GGT': 'RGN','TAA': 'TAA','TAC': 'TAY','TAG': 'TAG','TAT': 'TAY','TGA': 'TGR','TGC': 'TGY',
                   'TGG': 'TGR','TGT': 'TGY','TTA': 'YTN','TTG': 'YTN'}
    
    #Alternative Flatworm Mitochondrial Code condons
    dgn_dict_14 = {'AAA': 'AAH','AAC': 'AAH','AAG': 'AAG','AAT': 'AAH','AGA': 'AGN','AGC': 'AGN','AGG': 'AGN','AGT': 'AGN',
                   'ATA': 'ATH','ATC': 'ATH','ATG': 'ATG','ATT': 'ATH','CAA': 'CAR','CAG': 'CAR','CGA': 'CGN','CGC': 'CGN',
                   'CGG': 'CGN','CGT': 'CGN','CTA': 'YTN','CTC': 'YTN','CTG': 'YTN','CTT': 'YTN','GGA': 'GGN','GGC': 'GGN',
                   'GGG': 'GGN','GGT': 'GGN','TAA': 'TAH','TAC': 'TAH','TAG': 'TAG','TAT': 'TAH','TGA': 'TGR','TGC': 'TGY',
                   'TGG': 'TGR','TGT': 'TGY','TTA': 'YTN','TTG': 'YTN'}  
    
    dgn_tables = {1: dgn_dict_1, 2: dgn_dict_2, 3: dgn_dict_3, 4: dgn_dict_4, 5: dgn_dict_5, 6:
              dgn_dict_6, 9: dgn_dict_9, 10: dgn_dict_10, 11: dgn_dict_11, 12: dgn_dict_12, 13: dgn_dict_13}

    #add variable codons to working dictionary 
    working_dict = base_dict.copy()

    #Handle if codon table not supported. 
    try:
        working_dict.update(dgn_tables[table])
    except KeyError:
        print("Could not locate codon dictionary. Supported codon tables are 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, and 13")
        sys.exit(0)

    if str(_seqbuddy.alpha) == str(IUPAC.protein):
        raise TypeError("DNA sequence required, not protein.")
    if str(_seqbuddy.alpha) == str(IUPAC.unambiguous_rna) or str(_seqbuddy.alpha) == str(IUPAC.unambiguous_rna):
        raise TypeError("Please use a DNA seqeunce instead of an RNA sequence.")

    _seqbuddy = clean_seq(_seqbuddy)
    _seqbuddy = uppercase(_seqbuddy)
    for _rec in _seqbuddy.records:
        # shift reading frame
        _rec.seq = _rec.seq[reading_frame-1:]
        seq_length = len(str(_rec.seq))
        i=0
        degen_string=""
        while i < seq_length:
            codon = str(_rec.seq[i:i+3])           
            degen_string += working_dict[codon] if codon in working_dict else codon
            i = i+3
        _rec.seq = Seq(str(degen_string), alphabet=IUPAC.ambiguous_dna)
    return _seqbuddy




# ################################################# COMMAND LINE UI ################################################## #
def argparse_init():
    import argparse

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

    seqbuddy = []
    seq_set = ""

    try:
        for seq_set in in_args.sequence:
            if isinstance(seq_set, TextIOWrapper) and seq_set.buffer.raw.isatty():
                sys.exit("Warning: No input detected. Process will be aborted.")
            seq_set = SeqBuddy(seq_set, in_args.in_format, in_args.out_format, in_args.alpha)
            seqbuddy += seq_set.records

        seqbuddy = SeqBuddy(seqbuddy, seq_set.in_format, seq_set.out_format, seq_set.alpha)
    except GuessError:
        sys.exit("Error: SeqBuddy could not understand your input. "
                 "Check the file path or try specifying an input type with -f")

    return in_args, seqbuddy


def command_line_ui(in_args, seqbuddy, skip_exit=False):
    # ############################################# INTERNAL FUNCTION ################################################ #
    def _print_recs(_seqbuddy):
        if in_args.test:
            _stderr("*** Test passed ***\n", in_args.quiet)
            pass

        elif in_args.in_place:
            _in_place(str(_seqbuddy), in_args.sequence[0])

        else:
            _stdout("{0}\n".format(str(_seqbuddy).rstrip()))

    def _in_place(_output, _path):
        if not os.path.exists(_path):
            _stderr("Warning: The -i flag was passed in, but the positional argument doesn't seem to be a "
                    "file. Nothing was written.\n", in_args.quiet)
            _stderr("%s\n" % _output.strip(), in_args.quiet)
        else:
            with open(os.path.abspath(_path), "w") as _ofile:
                _ofile.write(_output)
            _stderr("File over-written at:\n%s\n" % os.path.abspath(_path), in_args.quiet)

    def _get_blast_binaries():
        blastp = None
        blastn = None
        blastdbcmd = None
        if in_args.params:
            for _param in in_args.params:
                binary = Popen("%s -version" % _param, stdout=PIPE, shell=True).communicate()
                binary = re.search("([a-z])*[^:]", binary[0].decode("utf-8"))
                binary = binary.group(0)
                if binary == "blastp":
                    blastp = _param
                elif binary == "blastn":
                    blastn = _param
                elif binary == "blastdbcmd":
                    blastdbcmd = _param

        blastp = blastp if blastp else which("blastp")
        blastn = blastn if blastn else which("blastn")
        blastdbcmd = blastdbcmd if blastdbcmd else which("blastdbcmd")

        return {"blastdbcmd": blastdbcmd, "blastp": blastp, "blastn": blastn}

    def _raise_error(_err):
        sys.exit("{0}: {1}\n".format(_err.__class__.__name__, str(_err)))

    def _exit(tool, skip=skip_exit):
        if skip:
            return
        usage = br.Usage()
        usage.increment("SeqBuddy", VERSION.short(), tool)
        usage.save()
        sys.exit()

    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Add feature
    if in_args.add_feature:
        # _type, _location, _strand=None, _qualifiers=None, _pattern=None
        strand = None
        qualifiers = None
        pattern = None
        if len(in_args.add_feature) < 2:
            raise AttributeError("Too few parameters provided. Must provide at least a feature type and location.")
        elif len(in_args.add_feature) == 5:
            strand = in_args.add_feature[2]
            qualifiers = in_args.add_feature[3]
            pattern = in_args.add_feature[4]
        elif len(in_args.add_feature) == 4:
            if in_args.add_feature[2] in ['+', 'plus', 'sense', 'pos', 'positive', '1', 1, '-', 'minus', 'anti',
                                          'antisense', 'anti-sense', 'neg', 'negative', '-1', -1, '0', 0]:
                strand = in_args.add_feature[2]
                if '=' in in_args.add_feature[3] or ':' in in_args.add_feature[3]:
                    qualifiers = in_args.add_feature[3]
                else:
                    pattern = in_args.add_feature[3]
            else:
                qualifiers = in_args.add_feature[3]
                pattern = in_args.add_feature[3]

        elif len(in_args.add_feature) == 3:
            if in_args.add_feature[2] in ['+', 'plus', 'sense', 'pos', 'positive', '1', 1, '-', 'minus', 'anti',
                                          'antisense', 'anti-sense', 'neg', 'negative', '-1', -1, '0', 0]:
                strand = in_args.add_feature[2]
            elif '=' in in_args.add_feature[2] or ':' in in_args.add_feature[2]:
                qualifiers = in_args.add_feature[2]
            else:
                pattern = in_args.add_feature[2]
        elif len(in_args.add_feature) == 2:
            pass
        else:
            raise AttributeError("Invalid parameters were provided.")
        ftype = in_args.add_feature[0]
        flocation = in_args.add_feature[1]
        _print_recs(add_feature(seqbuddy, ftype, flocation, _strand=strand, _qualifiers=qualifiers, _pattern=pattern))
        _exit("add_feature")

    # Average length of sequences
    if in_args.ave_seq_length:
        clean = False if not in_args.ave_seq_length[0] or in_args.ave_seq_length[0] != "clean" else True
        _stdout("%s\n" % round(ave_seq_length(seqbuddy, clean), 2))
        _exit("ave_seq_length")

    # Back Transcribe
    if in_args.back_transcribe:
        if seqbuddy.alpha != IUPAC.ambiguous_rna:
            raise ValueError("You need to provide an RNA sequence.")
        _print_recs(rna2dna(seqbuddy))
        _exit("back_transcribe")

    # Back translate CDS
    if in_args.back_translate:
        if in_args.back_translate[0]:  # All this logic is to determine what mode is being used by the UI
            in_args.back_translate = in_args.back_translate[0]
            in_args.back_translate = [i.upper() for i in in_args.back_translate]
            mode = [i for i in in_args.back_translate if i in ['RANDOM', 'R', "OPTIMIZED", "O"]]
            mode = "RANDOM" if len(mode) == 0 else mode[0]
            species = [i for i in in_args.back_translate if i in ['HUMAN', 'H', "MOUSE", "M",
                                                                  "YEAST", "Y", "ECOLI", "E"]]
            species = None if len(species) == 0 else species[0]
        else:
            mode = "RANDOM"
            species = None
        _print_recs(back_translate(seqbuddy, mode, species))
        _exit("back_translate")

    # BL2SEQ
    if in_args.bl2seq:
        output = bl2seq(seqbuddy)
        _stdout(output[1])
        _exit("bl2seq")

    # BLAST  ToDo: Determine if this can be refactored
    if in_args.blast:
        blast_binaries = _get_blast_binaries()
        blast_binary_path = blast_binaries["blastp"] if seqbuddy.alpha == IUPAC.protein else blast_binaries["blastn"]
        try:
            blast_res = blast(seqbuddy, in_args.blast, blast_path=blast_binary_path,
                              blastdbcmd=blast_binaries["blastdbcmd"])

        except FileNotFoundError as e:
            raise FileNotFoundError("%s binary not found, explicitly set with the -p flag.\nTo pass in the path to "
                                    "both blast(p/n) and blastdbcmd, separate them with a space." % e)

        if blast_res:
            _print_recs(blast_res)
        _exit("blast")

    # Clean Seq
    if in_args.clean_seq:
        if in_args.clean_seq[0] == "strict":
            _print_recs(clean_seq(seqbuddy, ambiguous=False))
        else:
            _print_recs(clean_seq(seqbuddy, ambiguous=True))
        _exit("clean_seq")

    # Combine feature sets from two files into one
    if in_args.combine_features:
        file1, file2 = in_args.sequence[:2]
        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)
        _print_recs(combine_features(file1, file2))
        _exit("combine_features")

    # Complement
    if in_args.complement:
        _print_recs(complement(seqbuddy))
        _exit("complement")

    # Concatenate sequences
    if in_args.concat_seqs:
        clean = False if not in_args.concat_seqs[0] or in_args.concat_seqs[0] != "clean" else True
        seqbuddy = concat_seqs(seqbuddy, clean)
        if in_args.out_format:
            seqbuddy.out_format = in_args.out_format
        _print_recs(seqbuddy)
        _exit("concat_seqs")

    # Codon counter
    if in_args.count_codons:
        try:
            codon_table = count_codons(seqbuddy)[1]
            for sequence_id in codon_table:
                _stdout('#### {0} ####\n'.format(sequence_id))
                _stdout('Codon\tAmino Acid\tNum\tPercent\n')
                for codon in codon_table[sequence_id]:
                    data = codon_table[sequence_id][codon]
                    _stdout('{0}\t{1}\t{2}\t{3}\n'.format(codon, data[0], data[1], data[2]))
                _stdout('\n')
        except TypeError as e:
            _raise_error(e)
        _exit("count_codons")

    # Count residues
    if in_args.count_residues:
        output = count_residues(seqbuddy)[1]
        for sequence in output:
            print(sequence)
            sorted_residues = OrderedDict(sorted(output[sequence].items()))
            for residue in sorted_residues:
                try:
                    print("{0}:\t{1}\t{2} %".format(residue, output[sequence][residue][0],
                                                    round(output[sequence][residue][1] * 100, 2)))
                except TypeError:
                    print("{0}:\t{1}".format(residue, output[sequence][residue]))
            print()
        _exit("count_residues")

    # Delete features
    if in_args.delete_features:
        for next_pattern in in_args.delete_features:
            delete_features(seqbuddy, next_pattern)
        _print_recs(seqbuddy)
        _exit("delete_features")

    # Delete sequences above threshold
    if in_args.delete_large:
        _print_recs(delete_large(seqbuddy, in_args.delete_large))
        _exit("delete_large")

    # Delete metadata
    if in_args.delete_metadata:
        _print_recs(delete_metadata(seqbuddy))
        _exit("delete_metadata")

    # Delete records
    if in_args.delete_records:
        if in_args.params:
            columns = int(in_args.params[0])
        else:
            columns = 1

        new_list = SeqBuddy(list(seqbuddy.records))
        deleted_seqs = []
        for next_pattern in in_args.delete_records:
            deleted_seqs += pull_recs(copy(new_list), next_pattern).records
            delete_records(new_list, next_pattern)

        if len(deleted_seqs) > 0 and not in_args.quiet:
                counter = 1
                output = "# ####################### Deleted records ######################## #\n"
                for seq in deleted_seqs:
                    output += "%s\t" % seq.id
                    if counter % columns == 0:
                        output = "%s\n" % output.strip()
                    counter += 1
                output = "%s\n# ################################################################ #\n" % output.strip()
                _stderr(output)

        if len(deleted_seqs) == 0:
            _stderr("# ################################################################ #\n")
            _stderr("# No sequence identifiers match %s\n" % ", ".join(in_args.delete_records))
            _stderr("# ################################################################ #\n")

        new_list.out_format = in_args.out_format if in_args.out_format else seqbuddy.out_format
        _print_recs(new_list)
        _exit("delete_records")

    # Delete repeats
    if in_args.delete_repeats:
        if in_args.delete_repeats[0]:
            columns = int(in_args.delete_repeats[0])
        else:
            columns = 1

        unique, rep_ids, rep_seqs, out_string = find_repeats(seqbuddy)
        stderr_output = ""
        if len(rep_ids) > 0:
            stderr_output += "# Records with duplicate ids deleted (first instance retained)\n"
            counter = 1
            for seq in rep_ids:
                stderr_output += "%s\t" % seq
                if counter % columns == 0:
                    stderr_output = "%s\n" % stderr_output.strip()
                counter += 1
            stderr_output = "%s\n\n" % stderr_output.strip()
            seqbuddy = delete_repeats(seqbuddy, 'ids')
            unique, rep_ids, rep_seqs, out_string = find_repeats(seqbuddy)

        rep_seq_ids = []
        for seq in rep_seqs:
            rep_seq_ids.append([])
            for rep_seq_id in rep_seqs[seq]:
                rep_seq_ids[-1].append(rep_seq_id)

        if len(rep_seq_ids) > 0:
            stderr_output += "# Records with duplicate sequence deleted (first instance retained)\n"
            counter = 1
            for rep_seqs in rep_seq_ids:
                stderr_output += "["
                for rep_seq in rep_seqs:
                    stderr_output += "%s, " % rep_seq
                stderr_output = "%s], " % stderr_output.strip(", ")
                if counter % columns == 0:
                    stderr_output = "%s\n" % stderr_output.strip(", ")
                counter += 1
            stderr_output = "%s\n" % stderr_output.strip(", ")

        if stderr_output != "" and in_args.quiet:
            _print_recs(delete_repeats(seqbuddy, 'seqs'))

        elif stderr_output != "":
            _stderr("# ################################################################ #\n")
            _stderr("%s\n" % stderr_output.strip())
            _stderr("# ################################################################ #\n\n")

            _print_recs(delete_repeats(seqbuddy, 'seqs'))

        else:
            _stderr("No duplicate records found\n")
        _exit("delete_repeats")

    # Delete sequences below threshold
    if in_args.delete_small:
        _print_recs(delete_small(seqbuddy, in_args.delete_small))
        _exit("delete_small")

    # degenerate_sequence
    if in_args.degenerate_sequence:
        table, reading_frame = 1, 1
        in_args.degenerate_sequence = in_args.degenerate_sequence[0]
        
        # check if to make sure letters are not in argument
        check_numbers = [n for n in in_args.degenerate_sequence if n.isdigit()]
        if len(check_numbers) != len(in_args.degenerate_sequence):    
                raise AttributeError('Please use integers not strings')
        
        # notify user only need two arguments
        if len(in_args.degenerate_sequence) > 2:
            raise AttributeError(
                'Too many attributes provided please only provided 1 or 2 parameters (table or table reading frame') 
       
       # if no argument provided will use table 1 first reading frame as default(set above)
        if not in_args.degenerate_sequence:
            pass
        # if one argument provided will set the given argument as the codon 
        # table.
        elif len(in_args.degenerate_sequence) == 1:
            print(
                "Only one parameter detected, will use the given parameter as a codon table and start at the first reading frame")
            table = int(in_args.degenerate_sequence[0])
            read_frame = 1
        else:
            table = int(in_args.degenerate_sequence[0])
            reading_frame = int(in_args.degenerate_sequence[1])
        _print_recs(degenerate_sequence(seqbuddy, table, reading_frame))

    # Extract regions
    if in_args.extract_region:
        _print_recs(extract_range(seqbuddy, *in_args.extract_region))
        _exit("extract_region")

    # Find CpG
    if in_args.find_CpG:
        output = find_cpg(seqbuddy)
        if output[1]:
            out_string = ""
            for key, value in output[1].items():
                if value:
                    value = ["%s-%s" % (x[0], x[1]) for x in value]
                    out_string += "{0}: {1}\n".format(key, ", ".join(value))
            _stderr('########### Islands identified ###########\n%s\n##########################################\n\n' %
                    out_string.strip(), in_args.quiet)
        else:
            _stderr("# No Islands identified\n\n", in_args.quiet)
        _print_recs(seqbuddy)
        _exit("find_CpG")

    # Find pattern
    if in_args.find_pattern:
        for pattern in in_args.find_pattern:
            output = find_pattern(seqbuddy, pattern)
            out_string = ""
            num_matches = 0
            for key in output[1]:
                out_string += "{0}: {1}\n".format(key, ", ".join([str(x) for x in output[1][key]]))
                num_matches += len(output[1][key])
            _stderr("#### {0} matches found across {1} sequences for "
                    "pattern '{2}' ####\n".format(num_matches, len(output[1]), pattern), in_args.quiet)
            _stderr("%s\n" % out_string, in_args.quiet)
        _print_recs(seqbuddy)
        _exit("find_pattern")

    # Find repeat sequences or ids
    if in_args.find_repeats:
        columns = 1 if not in_args.find_repeats[0] else in_args.find_repeats[0]
        unique, rep_ids, rep_seqs, out_string = find_repeats(seqbuddy, columns)
        _stdout(out_string)
        _exit("find_repeats")

    # Find restriction sites
    if in_args.find_restriction_sites:
        min_cuts, max_cuts, enzymes, order = None, None, [], 'position'
        in_args.find_restriction_sites = in_args.find_restriction_sites[0]
        for param in in_args.find_restriction_sites:
            try:
                param = int(param)
            except ValueError:
                pass
            if type(param) == int:
                if not min_cuts:
                    min_cuts = param
                elif not max_cuts:
                    max_cuts = param
                else:
                    _raise_error(ValueError("To many integers in parameter list for find_restriction_sites."))
            elif param in ['alpha', 'position']:
                order = param
            else:
                enzymes.append(param)

        enzymes = ["commercial"] if len(enzymes) == 0 else enzymes
        max_cuts = int(min_cuts) if min_cuts and not max_cuts else max_cuts
        min_cuts = 1 if not min_cuts else min_cuts
        if max_cuts and min_cuts > max_cuts:
            temp = int(max_cuts)
            max_cuts = int(min_cuts)
            min_cuts = temp

        clean_seq(seqbuddy)
        try:
            output = find_restriction_sites(seqbuddy, enzymes, min_cuts, max_cuts)
        except TypeError as e:
            _raise_error(e)
            sys.exit()

        out_string = ''
        for tup in output[1]:
            out_string += "{0}\n".format(tup[0])
            restriction_list = tup[1]
            restriction_list = [[key, value] for key, value in restriction_list.items()]
            restriction_list = sorted(restriction_list, key=lambda l: str(l[0])) if order == 'alpha' else \
                sorted(restriction_list, key=lambda l: l[1])

            for _enzyme in restriction_list:
                cut_sites = [str(x) for x in _enzyme[1]]
                out_string += "{0}:\t{1}\n".format(_enzyme[0], ", ".join(cut_sites))
            out_string += "\n"
        _stdout(out_string)
        _exit("find_restriction_sites")

    # Guess alphabet
    if in_args.guess_alphabet:
        for seq_set in in_args.sequence:
            if str(type(seq_set)) != "<class '_io.TextIOWrapper'>":
                seqbuddy = SeqBuddy(seq_set)
                _stdout("%s\t-->\t" % seq_set)
            if seqbuddy.alpha == IUPAC.protein:
                _stdout("prot\n")
            elif seqbuddy.alpha == IUPAC.ambiguous_dna:
                _stdout("dna\n")
            elif seqbuddy.alpha == IUPAC.ambiguous_rna:
                _stdout("rna\n")
            else:
                _stdout("Undetermined\n")
        _exit("guess_alphabet")

    # Guess format
    if in_args.guess_format:
        if str(type(in_args.sequence[0])) == "<class '_io.TextIOWrapper'>":
            _stdout("{0}\n".format(seqbuddy.in_format))
        else:
            for seq_set in in_args.sequence:
                _stdout("%s\t-->\t%s\n" % (seq_set, SeqBuddy(seq_set).in_format))
        _exit("guess_format")

    # Hash sequence ids
    if in_args.hash_seq_ids:
        hash_length = in_args.hash_seq_ids[0] if in_args.hash_seq_ids[0] else 10
        seqbuddy, hash_map, hash_table = hash_sequence_ids(seqbuddy, hash_length)
        _stderr(hash_table, in_args.quiet)
        _print_recs(seqbuddy)
        _exit("hash_seq_ids")

    # Insert Seq
    if in_args.insert_seq:
        location = in_args.insert_seq[1]
        if location not in ['start', 'end']:
            try:
                location = int(location)
            except ValueError("Location must be start, end, or integer index") as e:
                _raise_error(e)
        _print_recs(insert_sequence(seqbuddy, in_args.insert_seq[0], location))
        _exit("insert_seq")

    # Calculate Isoelectric Point
    if in_args.isoelectric_point:
        try:
            isoelectric_points = isoelectric_point(seqbuddy)[1]
            _stderr("ID\tpI\n")
            for rec_id in isoelectric_points:
                print("{0}\t{1}".format(rec_id, isoelectric_points[rec_id]))
        except ValueError as e:
            _raise_error(e)
        _exit("isoelectric_point")

    # List features
    if in_args.list_features:
        feature_table = list_features(seqbuddy)
        for rec_id in feature_table:
            _stdout('#### {0} ####\n'.format(rec_id))
            out_string = ''
            if feature_table[rec_id] is not None:
                for feat in feature_table[rec_id]:
                    out_string += '{0}\n'.format(feat.type)
                    if isinstance(feat.location, CompoundLocation):
                        out_string += '\tLocations:\n'
                        for part in feat.location.parts:
                            out_string += '\t\t{0}-{1}\n'.format(part.start, part.end)
                    elif isinstance(feat.location, FeatureLocation):
                        out_string += '\tLocation:\n'
                        out_string += '\t\t{0}-{1}\n'.format(feat.location.start, feat.location.end)
                    if feat.strand == 1:
                        out_string += '\tStrand: Sense(+)\n'
                    elif feat.strand == -1:
                        out_string += '\tStrand: Antisense(-)\n'
                    if str(feat.id) != '<unknown id>':
                        out_string += '\tID: {0}\n'.format(feat.id)
                    if len(feat.qualifiers) > 0:
                        out_string += '\tQualifiers:\n'
                        for key in feat.qualifiers:
                            out_string += '\t\t{0}: {1}\n'.format(key, str(feat.qualifiers[key]).strip("[]'"))
                    if feat.ref is not None:
                        out_string += '\ref: {0}'.format(feat.ref)
            else:
                out_string = 'None\n'
            _stdout('%s\n' % out_string)
        _exit("list_features")

    # List identifiers
    if in_args.list_ids:
        columns = 1 if not in_args.list_ids[0] else in_args.list_ids[0]
        _stdout(list_ids(seqbuddy, columns))
        _exit("list_ids")

    # Lowercase
    if in_args.lowercase:
        _print_recs(lowercase(seqbuddy))
        _exit("lowercase")

    # Map features from cDNA over to protein
    if in_args.map_features_dna2prot:
        file1, file2 = in_args.sequence[:2]
        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)
        if file1.alpha == file2.alpha:
            raise ValueError("You must provide one DNA file and one protein file")
        if file1.alpha == IUPAC.protein:
            prot = file1
            dna = file2
        else:
            in_args.sequence[0] = in_args.sequence[1]  # in case the -i flag is thrown
            prot = file2
            dna = file1
        _print_recs(map_features_dna2prot(dna, prot))
        _exit("map_features_dna2prot")

    # Map features from protein over to cDNA
    if in_args.map_features_prot2dna:
        file1, file2 = in_args.sequence[:2]
        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)
        if file1.alpha == file2.alpha:
            raise ValueError("You must provide one DNA file and one protein file")
        if file1.alpha != IUPAC.protein:
            dna = file1
            prot = file2
        else:
            in_args.sequence[0] = in_args.sequence[1]  # in case the -i flag is thrown
            dna = file2
            prot = file1
        _print_recs(map_features_prot2dna(prot, dna))
        _exit("map_features_prot2dna")

    # Merge
    if in_args.merge:
        _print_recs(seqbuddy)
        _exit("merge")

    # Calculate Molecular Weight
    if in_args.molecular_weight:
        lists = molecular_weight(seqbuddy)[1]
        if seqbuddy.alpha == (IUPAC.ambiguous_dna or IUPAC.unambiguous_dna):
            _stderr("ID\tssDNA\tdsDNA\n")
        elif seqbuddy.alpha == (IUPAC.ambiguous_rna or IUPAC.unambiguous_rna):
            _stderr("ID\tssRNA\n")
        else:
            _stderr("ID\tProtein\n")
        for indx, value in enumerate(lists['ids']):
            if len(lists['masses_ds']) != 0:
                print("{0}\t{1}\t{2}".format(value, lists['masses_ss'][indx], lists['masses_ds'][indx]))
            else:
                print("{0}\t{1}".format(value, lists['masses_ss'][indx]))
        _exit("molecular_weight")

    # Count number of sequences in a file
    if in_args.num_seqs:
        _stdout("%s\n" % num_seqs(seqbuddy))
        _exit("num_seqs")

    # Order sequence features alphabetically
    if in_args.order_features_alphabetically:
        reverse = True if in_args.order_features_alphabetically[0] \
            and in_args.order_features_alphabetically[0] == "rev" else False
        _print_recs(order_features_alphabetically(seqbuddy, reverse))
        _exit("order_features_alphabetically")

    # Order sequence features by their position in the sequence
    if in_args.order_features_by_position:
        reverse = True if in_args.order_features_by_position[0] \
            and in_args.order_features_by_position[0] == "rev" else False
        _print_recs(order_features_by_position(seqbuddy, reverse))
        _exit("order_features_by_position")

    # Order ids
    if in_args.order_ids:
        reverse = True if in_args.order_ids[0] and in_args.order_ids[0] == "rev" else False
        _print_recs(order_ids(seqbuddy, _reverse=reverse))
        _exit("order_ids")

    # Order ids randomly
    if in_args.order_ids_randomly:
        _print_recs(order_ids_randomly(seqbuddy))
        _exit("order_ids_randomly")

    # Pull random records
    if in_args.pull_random_record:
        count = 1 if not in_args.pull_random_record[0] else in_args.pull_random_record[0]
        _print_recs(pull_random_recs(seqbuddy, count))
        _exit("pull_random_record")

    # Pull sequence ends
    if in_args.pull_record_ends:
        try:
            _print_recs(pull_record_ends(seqbuddy, *in_args.pull_record_ends))
            sys.exit()
        except ValueError:
            pass
        except AttributeError:
            pass
        try:
            _print_recs(pull_record_ends(seqbuddy, in_args.pull_record_ends[1], in_args.pull_record_ends[0]))
        except ValueError:
            _raise_error(ValueError("Arguments are <amount (int)> <front|rear>"))
        except AttributeError:
            _raise_error(AttributeError("Choose 'front' or 'rear' to specify where the sequence "
                                        "should come from in pull_record_ends"))
        _exit("pull_record_ends")

    # Pull records
    if in_args.pull_records:
        _print_recs(pull_recs(seqbuddy, in_args.pull_records))
        _exit("pull_records")

    # Purge
    if in_args.purge:
        purged_seqs, deleted, record_map = purge(seqbuddy, in_args.purge)
        _stderr(record_map, in_args.quiet)
        _print_recs(purged_seqs)
        _exit("purge")

    # Raw Seq
    if in_args.raw_seq:
        output = raw_seq(seqbuddy)
        if in_args.in_place:
            _in_place(output, in_args.sequence[0])
        else:
            _stdout(output)
        _exit("raw_seq")

    # Renaming
    if in_args.rename_ids:
        num = 0 if not in_args.params else int(in_args.params[0])
        _print_recs(rename(seqbuddy, in_args.rename_ids[0], in_args.rename_ids[1], num))
        _exit("rename_ids")

    # Reverse complement
    if in_args.reverse_complement:
        _print_recs(reverse_complement(seqbuddy))
        _exit("reverse_complement")

    # Screw formats
    if in_args.screw_formats:
        if in_args.screw_formats not in FORMATS:
            _stderr("Error: unknown format '%s'\n" % in_args.screw_formats)
        else:
            seqbuddy.out_format = in_args.screw_formats
            if in_args.in_place:  # Need to change the file extension
                os.remove(in_args.sequence[0])
                in_args.sequence[0] = ".".join(os.path.abspath(in_args.sequence[0]).split(".")[:-1]) + \
                                      "." + seqbuddy.out_format
                open(in_args.sequence[0], "w").close()
            _print_recs(seqbuddy)
        _exit("screw_formats")

    # Shift reading frame
    if in_args.select_frame:
        _print_recs(select_frame(seqbuddy, in_args.select_frame))
        _exit("select_frame")

    # Shuffle Seqs
    if in_args.shuffle_seqs:
        _print_recs(shuffle_seqs(seqbuddy))
        _exit("shuffle_seqs")

    # Split sequences by taxa.
    if in_args.split_by_taxa:
        in_args.in_place = True
        out_dir = os.path.abspath(in_args.split_by_taxa[1])
        os.makedirs(out_dir, exist_ok=True)
        taxa_groups = split_by_taxa(seqbuddy, in_args.split_by_taxa[0])
        check_quiet = in_args.quiet  # 'quiet' must be toggled to 'on' _print_recs() here.
        in_args.quiet = True
        for taxa_heading in taxa_groups:
            seqbuddy.records = taxa_groups[taxa_heading]
            in_args.sequence[0] = "%s/%s.%s" % (out_dir, taxa_heading, _format_to_extension(seqbuddy.out_format))
            _stderr("New file: %s\n" % in_args.sequence[0], check_quiet)
            open(in_args.sequence[0], "w").close()
            _print_recs(seqbuddy)
        _exit("split_by_taxa")

    # Split sequences into files
    if in_args.split_to_files:
        in_args.in_place = True
        out_dir = os.path.abspath(in_args.split_to_files)
        os.makedirs(out_dir, exist_ok=True)
        check_quiet = in_args.quiet  # 'quiet' must be toggled to 'on' _print_recs() here.
        in_args.quiet = True
        for buddy in split_file(seqbuddy):
            seqbuddy.records = buddy.records
            ext = _format_to_extension(seqbuddy.out_format)
            in_args.sequence[0] = "%s/%s.%s" % (out_dir, seqbuddy.records[0].id, ext)
            _stderr("New file: %s\n" % in_args.sequence[0], check_quiet)
            open(in_args.sequence[0], "w").close()
            _print_recs(seqbuddy)
        _exit("split_to_files")

    # Transcribe
    if in_args.transcribe:
        if seqbuddy.alpha != IUPAC.ambiguous_dna:
            raise ValueError("You need to provide a DNA sequence.")
        _print_recs(dna2rna(seqbuddy))
        _exit("transcribe")

    # Translate CDS
    if in_args.translate:
        if seqbuddy.alpha == IUPAC.protein:
            raise ValueError("You need to supply DNA or RNA sequences to translate")
        _print_recs(translate_cds(seqbuddy, quiet=in_args.quiet))
        _exit("translate")

    # Translate 6 reading frames
    if in_args.translate6frames:
        if seqbuddy.alpha == IUPAC.protein:
            raise ValueError("You need to supply DNA or RNA sequences to translate")
        seqbuddy = translate6frames(seqbuddy)
        if in_args.out_format:
            seqbuddy.out_format = in_args.out_format
        _print_recs(seqbuddy)
        _exit("translate6frames")

    # Uppercase
    if in_args.uppercase:
        _print_recs(uppercase(seqbuddy))
        _exit("uppercase")

if __name__ == '__main__':
    try:
        command_line_ui(*argparse_init())
    except (KeyboardInterrupt, GuessError) as _e:
        print(_e)
    except SystemExit:
        pass
    except Exception as _e:
        br.send_traceback("SeqBuddy", _e)
