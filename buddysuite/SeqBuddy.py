#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: SeqBuddy.py
author: Stephen R. Bond
email: biologyguy@gmail.com
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description:
Collection of functions that do fun stuff with sequences. Pull them into a script, or run as a command line tool.
"""

# ##################################################### IMPORTS ###################################################### #
# Standard library
# from pprint import pprint
# import pdb
# import time
from __future__ import print_function

# BuddySuite specific
import warnings

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
from os.path import join
import re
import string
import zipfile
import shutil
import time
import urllib.parse
import urllib.request
import urllib.error
from copy import deepcopy
from random import sample, Random
from math import floor, ceil, log
from subprocess import Popen, PIPE
from multiprocessing import Lock
from shutil import which
from hashlib import md5
from io import StringIO, TextIOWrapper
from collections import OrderedDict
from xml.sax import SAXParseException

# Third party
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import RestrictionBatch, CommOnly, AllEnzymes, Analysis
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Nexus.Trees import TreeError

# ##################################################### WISH LIST #################################################### #
'''
def sim_ident(matrix):  # Return the pairwise similarity and identity scores among sequences
    """
    (This may be more appropriate for AlignBuddy)
    :param matrix:
    """
    x = matrix
    return x


def cd_hit(seqbuddy, threshold):
    """
    :param seqbuddy: SeqBuddy object
    :param threshold: Maximum sequence identity
    """
    return seqbuddy


def auto_annotate():
    """
    Find common plasmid features in sequences
    """
    return


def mutate():
    """
    Apply some model of evolution to generate new sequences from input
    :return:
    """
    return


def random_aa(_length, number, matrix):
    """
    create random prot sequences. Not sure exactly how to implement this yet, because it would theoretically not start
    from input sequences...
    :param _length:
    :param number:
    :param matrix:
    """
    x = [_length, number, matrix]
    return x


def random_dna(_length, number, matrix):
    """
    create random DNA sequences.
    :param _length:
    :param number:
    :param matrix:
    :return:
    """
    x = [_length, number, matrix]
    return x


def divergence_value():
    """
    http://bioinformatics.org/sms/uneven.html
    """
    return


def incremental_rename(query, replace):
    """
    Append a number to the end of each replacement to ensure unique ids
    :param query:
    :param replace:
    """
    x = (query, replace)
    return x


def keep_features(seqbuddy, regex(s)):
    """
    Opposite of delete features
    """
    return
'''
# - Allow batch calls. E.g., if 6 files are fed in as input, run the SeqBuddy command independently on each
# - Check on memory requirements before execution
# - Execution timer, for long running jobs
# - Sort out a good way to manage 'lazy' imports (might not be that important)
# - Try to speed things up by reading in all sequence data only when necessary

# ###################################################### GLOBALS ##################################################### #
VERSION = br.Version("SeqBuddy", 1, "4.0", br.contributor_list, {"year": 2021, "month": 12, "day": 5})
OUTPUT_FORMATS = ["ids", "accessions", "summary", "full-summary", "clustal", "embl", "fasta", "fastq", "fastq-sanger",
                  "fastq-solexa", "fastq-illumina", "genbank", "gb", "imgt", "nexus", "nexuss", "nexusi", "phylip",
                  "phd", "phylip-relaxed", "phylipss", "phylipsr", "raw", "sff", "stockholm", "tab", "qual"]


# ##################################################### SEQBUDDY ##################################################### #
class SeqBuddy(object):
    """
    Core class.
    Open a file or read a handle and parse, or convert raw into a Seq object
    """
    def __init__(self, sb_input, in_format=None, out_format=None, alpha=None):
        # ####  IN AND OUT FORMATS  #### #
        # Holders for input type. Used for some error handling below
        in_handle = None
        raw_sequence = None
        in_file = None
        self.alpha = None if alpha is None else alpha.lower()
        self.hash_map = OrderedDict()  # This is only used by functions that use hash_id()

        # SeqBuddy obj
        if sb_input.__class__.__name__ == "SeqBuddy":
            sb_input = make_copy(sb_input)

        # Handles
        if str(type(sb_input)) == "<class '_io.TextIOWrapper'>":
            if not sb_input.seekable():  # Deal with input streams (e.g., stdout pipes)
                input_txt = sb_input.read()
                if re.search("Buddy::.+ has crashed with the following traceback", input_txt):
                    print(input_txt)
                    sys.exit()
                temp = StringIO(br.utf_encode(input_txt))
                sb_input = temp
            sb_input.seek(0)
            in_handle = sb_input.read()
            sb_input.seek(0)

        # Raw sequences
        if in_format == "raw":
            in_format = "fasta"
            out_format = "fasta" if not out_format else out_format
            if type(sb_input) == str:
                sb_input = br.utf_encode(sb_input)
            else:
                sb_input = br.utf_encode(sb_input.read())
            sb_input = [SeqRecord(Seq(sb_input), id="raw_input", description="")]

        # Plain text in a specific format
        if type(sb_input) == str and not os.path.isfile(sb_input):
            raw_sequence = br.utf_encode(sb_input)
            temp = StringIO(sb_input)
            sb_input = temp
            sb_input.seek(0)

        # File paths
        try:
            if os.path.isfile(sb_input):
                in_file = sb_input
                with open(sb_input, "r", encoding="utf-8") as ifile:
                    sb_input = StringIO(ifile.read())

        except TypeError:  # This happens when testing something other than a string.
            pass

        if in_format:
            self.in_format = br.parse_format(in_format)

        elif sb_input.__class__.__name__ == "SeqBuddy":
            self.in_format = sb_input.in_format

        else:
            self.in_format = br.guess_format(sb_input)
            if self.in_format == "empty file":
                self.in_format = "fasta"
            self.out_format = str(self.in_format) if not out_format else br.parse_format(out_format)

        if not self.in_format:
            if in_file:
                raise br.GuessError("Could not determine format from sb_input file '%s'.\n"
                                    "Try explicitly setting with -f flag." % in_file)
            elif raw_sequence:
                raise br.GuessError("File not found, or could not determine format from raw input\n --> %s ...\n"
                                    "Try explicitly setting with -f flag." % raw_sequence[:60])
            elif in_handle:
                raise br.GuessError("Could not determine format from input file-like object\n --> %s ...\n"
                                    "Try explicitly setting with -f flag." % in_handle[:50])
            else:
                raise br.GuessError("Unable to determine format or input type. "
                                    "Please check how SeqBuddy is being called.")

        self.out_format = self.in_format if not out_format else br.parse_format(out_format)
        self.out_format = "clustal" if self.out_format == "blast-aln" else self.out_format  # Cant write blast-aln

        # ####  RECORDS  #### #
        if sb_input.__class__.__name__ == "SeqBuddy":
            sequences = sb_input.records

        elif isinstance(sb_input, list):
            # make sure that the list is actually SeqIO records (just test a few...)
            rand_sample = sb_input if len(sb_input) < 5 else sample(sb_input, 5)
            for seq in rand_sample:
                if type(seq) != SeqRecord:
                    raise TypeError("Seqlist is not populated with SeqRecords.")
            sequences = sb_input

        elif str(type(sb_input)) == "<class '_io.TextIOWrapper'>" or isinstance(sb_input, StringIO):
            if self.in_format in ["phylipss", "phylipsr"]:
                relaxed = False if self.in_format == "phylipss" else True
                aligns = br.phylip_sequential_read(sb_input.read(), relaxed=relaxed)
                sequences = []
                for align in aligns:
                    sequences += [rec for rec in align]
            else:
                sequences = list(SeqIO.parse(sb_input, self.in_format))

        elif os.path.isfile(sb_input):
            with open(sb_input, "r", encoding="utf-8") as sb_input:
                if self.in_format in ["phylipss", "phylipsr"]:
                    relaxed = False if self.in_format == "phylipss" else True
                    aligns = br.phylip_sequential_read(sb_input.read(), relaxed=relaxed)
                    sequences = []
                    for align in aligns:
                        sequences += [rec for rec in align]
                else:
                    sequences = list(SeqIO.parse(sb_input, self.in_format))
        else:
            sequences = [SeqRecord(Seq(sb_input))]  # may be unreachable?

        if self.alpha is None:
            self.alpha = guess_alphabet(sequences)
        elif "protein".startswith(str(self.alpha)) or self.alpha in ['pep', "protein"]:
            self.alpha = "protein"
        elif "dna".startswith(str(self.alpha)) or self.alpha in ['cds', "dna"]:
            self.alpha = "DNA"
        elif "rna".startswith(str(self.alpha)) or self.alpha == "rna":
            self.alpha = "RNA"
        else:
            br._stderr("WARNING: Alphabet not recognized. Correct alphabet will be guessed.\n")
            self.alpha = guess_alphabet(sequences)

        for seq in sequences:
            seq.annotations["molecule_type"] = self.alpha

        # The NEXUS parser adds '.copy' to any repeat taxa, strip that off...
        if self.in_format == "nexus":
            for rec in sequences:
                rec.id = re.sub(r"\.copy[0-9]*$", "", rec.id)

        self.records = sequences
        self.memory_footprint = sum([len(rec) for rec in sequences])
        self.hidden_recs = []  # If anything is in here it will be printed with __str__() or write()

    def __str__(self):
        records = self.records + self.hidden_recs
        if len(records) == 0:
            return "Error: No sequences in object.\n"

        # There is a weird bug in genbank write() that concatenates dots to the organism name (if set).
        # The following is a work around...
        self.out_format = self.out_format.lower()
        if self.out_format in ["gb", "genbank"]:
            for rec in records:
                try:
                    if re.search(r"(\. )+", rec.annotations['organism']):
                        rec.annotations['organism'] = "."
                except KeyError:
                    pass

        output = ""
        if self.out_format == "phylipsr":
            output = br.phylip_sequential_out(self)

        elif self.out_format == "phylipss":
            output = br.phylip_sequential_out(self, relaxed=False)

        elif self.out_format.startswith("nexus"):
            try:
                output = br.nexus_out(self, self.out_format)
            except ValueError as e:
                if "Sequences must all be the same length" in str(e):
                    br._stderr("Warning: Alignment format detected but sequences are different lengths. "
                               "Format changed to fasta to accommodate proper printing of records.\n\n")
                    self.out_format = "fasta"
                else:
                    raise e
        elif self.out_format == "raw":
            output = "\n\n".join([str(rec.seq) for rec in records])

        if not output:
            tmp_dir = br.TempDir()
            with open("%s%sseqs.tmp" % (tmp_dir.path, os.path.sep), "w", encoding="utf-8") as _ofile:
                try:
                    SeqIO.write(records, _ofile, self.out_format)
                except ValueError as e:
                    if "Sequences must all be the same length" in str(e):
                        br._stderr("Warning: Alignment format detected but sequences are different lengths. "
                                   "Format changed to fasta to accommodate proper printing of records.\n\n")
                        SeqIO.write(records, _ofile, "fasta")
                    elif "Repeated name" in str(e) and self.out_format == "phylip":
                        br._stderr("Warning: Phylip format returned a 'repeat name' error, probably due to truncation. "
                                   "Attempting phylip-relaxed.\n")
                        SeqIO.write(records, _ofile, "phylip-relaxed")
                    else:
                        raise e

            with open("%s%sseqs.tmp" % (tmp_dir.path, os.path.sep), "r", encoding="utf-8") as ifile:
                output = ifile.read()

        return "%s\n" % output.rstrip()

    def __len__(self):
        return len(self.records)

    def to_dict(self):
        sb_copy = find_repeats(make_copy(self))
        if len(sb_copy.repeat_ids) > 0:
            raise RuntimeError("There are repeat IDs in self.records\n%s" %
                               ", ".join([key for key, recs in sb_copy.repeat_ids.items()]))

        records_dict = OrderedDict()
        for rec in self.records:
            records_dict[rec.id] = rec
        return records_dict

    def write(self, file_path, out_format=None):
        with open(file_path, "w", encoding="utf-8") as ofile:
            if out_format:
                out_format_save = str(self.out_format)
                self.out_format = out_format
                ofile.write(str(self))
                self.out_format = out_format_save
            else:
                ofile.write(str(self))
        return

    def print_hashmap(self):
        output = ""
        if self.hash_map:
            for _hash, orig_id in self.hash_map.items():
                output += "%s\t%s\n" % (_hash, orig_id)
        return output

    def reverse_hashmap(self):
        if self.hash_map:
            for _hash, seq_id in self.hash_map.items():
                rename(self, _hash, seq_id)
        return


# ################################################# HELPER FUNCTIONS ################################################# #
def _add_buddy_data(rec, key=None, data=None):
    """
    Append the buddy_data attribute (an OrderedDict()) to a BioPython.SeqRecord object
    :param rec: The SeqRecord object
    :param key: New Dict key
    :param data: Data to be put in the key location
    :return: Modified SeqRecord obj
    """
    if not hasattr(rec, 'buddy_data'):
        rec.buddy_data = OrderedDict()

    if key:
        if data:
            rec.buddy_data[key] = data
        else:
            if key not in rec.buddy_data:
                rec.buddy_data[key] = None
    return rec


def _check_for_blast_bin(blast_bin):
    """
    Check the user's system for the blast bin in $PATH, try to download if not.
    :param blast_bin: the name of the binary to look for (str)
    :return: success or failure (bool)
    """
    if which(blast_bin):
        return True
    else:
        br._stderr("%s binary not found. Please install BLAST+ executables.\n" % blast_bin)
        return False


def _feature_rc(feature, seq_len):
    """
    BioPython does not properly handle reverse complement of features, so implement it...
    :param feature: BioPython feature with either a CompoundLocation or FeatureLocation
    :param seq_len: The full length of the original sequence
    :return: SeqFeature object
    """
    if type(feature.location) == CompoundLocation:
        new_compound_location = []
        for sub_feature in feature.location.parts:
            sub_feature = _feature_rc(SeqFeature(sub_feature), seq_len)
            new_compound_location.append(sub_feature.location)
        feature.location = CompoundLocation(new_compound_location, feature.location.operator)

    elif type(feature.location) == FeatureLocation:
        end = seq_len - feature.location.end
        shift = end - feature.location.start
        feature = br.shift_features(feature, shift, seq_len)[0]
        if feature.strand is not None:
            feature.strand *= -1
    else:
        raise TypeError("_feature_rc requires a feature with either FeatureLocation or CompoundLocation, "
                        "not %s" % type(feature.location))
    return feature


class FeatureReMapper(object):
    """
    Build a list that maps original residues to new positions if residues have been removed
    This will not work if new columns are being added.
    :usage: Instantiate a new object, and for each position in the original sequence, call the 'extend' method,
            specifying whether that residue exists in the new alignment or not. Remap the features on the new sequence
            by calling the remap_features method.
    """
    def __init__(self, old_seq):
        self.old_seq = old_seq  # SeqRecord
        self.position_map = []
        self.starting_position_filled = False

    def extend(self, exists=True):
        """
        Iterates the position map, adding an index and whether the new sequence contains the residue.
        Extend() must be called exactly len(self.old_seq) times.
        :param exists: Specify whether the next residue exists or not
        """
        if len(self.old_seq.seq) < len(self.position_map) + 1:
            raise AttributeError("The position map has already been fully populated.")

        if not self.position_map:
            if exists:
                self.starting_position_filled = True
            self.position_map.append((0, exists))

        else:
            if exists and not self.starting_position_filled:
                self.position_map.append((0, True))
                self.starting_position_filled = True
            elif exists and self.starting_position_filled:
                self.position_map.append((self.position_map[-1][0] + 1, True))
            else:
                self.position_map.append((self.position_map[-1][0], False))
        return

    def remap_features(self, new_seq):
        """
        Add all the features from self.old_seq that still exist onto new_seq
        :param new_seq: SeqRecord containing the new sequence that was used to build self.position_map
        """
        if len(self.old_seq.seq) != len(self.position_map):
            raise AttributeError("The position map has not been fully populated.")

        new_features = []
        for feature in self.old_seq.features:
            feature = self._remap(feature)
            if feature:
                new_features.append(feature)
        new_seq.features = new_features
        return new_seq

    def _remap(self, feature):
        """
        Deal with the weirdness that is SeqRecord features... Compares self.old_seq features against self.position_map
        :param feature: A feature from self.old_seq
        """
        if type(feature.location) == FeatureLocation:
            for pos, present in self.position_map[feature.location.start:feature.location.end]:
                if present:
                    start = pos
                    end = self.position_map[feature.location.end - 1][0] + 1
                    location = FeatureLocation(start, end, strand=feature.location.strand)
                    feature.location = location
                    return feature
            else:
                return None

        else:  # CompoundLocation
            parts = []
            for sub_feature in feature.location.parts:
                sub_feature = self._remap(SeqFeature(sub_feature))
                if sub_feature:
                    parts.append(sub_feature.location)
            if len(parts) > 1:
                feature.location = CompoundLocation(parts, operator='order')
            elif len(parts) == 1:
                feature.location = FeatureLocation(parts[0].start, parts[0].end, strand=parts[0].strand)
            else:
                feature = None
            return feature


def guess_alphabet(seqbuddy):
    """
    Looks through the characters in the SeqBuddy records to determine the most likely alphabet
    Does not attempt to explicitly deal with weird cases (e.g., ambiguous residues).
    The user will need to specify an alphabet with the -a flag if using many non-standard characters in their sequences.
    :param seqbuddy: SeqBuddy object
    :return: DNA/RNA/protein
    """
    seq_list = seqbuddy if isinstance(seqbuddy, list) else seqbuddy.records
    seq_list = [str(x.seq) for x in seq_list]
    sequence = "".join(seq_list).upper()
    sequence = re.sub(r"[NX\-?]", "", sequence)

    if len(sequence) == 0:
        return None

    if 'U' in sequence:  # U is unique to RNA
        return "RNA"

    percent_dna = len(re.findall("[ATCG]", sequence)) / float(len(sequence))
    percent_protein = len(re.findall("[ACDEFGHIKLMNPQRSTVWXY]", sequence)) / float(len(sequence))
    if percent_dna > 0.85:  # odds that a sequence with no Us and such a high ATCG count be anything but DNA is low
        return "DNA"
    elif percent_protein > 0.85:
        return "protein"
    else:
        return None


def make_copy(seqbuddy):
    """
    Deepcopy a SeqBuddy object. The alphabet objects are not handled properly when deepcopy is called,
    so need to wrap it.
    :param seqbuddy: SeqBuddy object
    :return: SeqBuddy object
    """
    alphabet_list = [rec.annotations["molecule_type"] for rec in seqbuddy.records]
    hidden_alphabet_list = [rec.annotations["molecule_type"] for rec in seqbuddy.hidden_recs]
    _copy = deepcopy(seqbuddy)
    _copy.alpha = seqbuddy.alpha
    for indx, rec in enumerate(_copy.records):
        rec.annotations["molecule_type"] = alphabet_list[indx]

    for indx, rec in enumerate(_copy.hidden_recs):
        rec.annotations["molecule_type"] = hidden_alphabet_list[indx]
    return _copy


def _prepare_restriction_sites(parameters):
    min_cuts, max_cuts, _enzymes, order, topology = None, None, [], 'position', None
    for param in parameters:
        try:
            param = int(param)
        except ValueError:
            pass
        if type(param) == int:
            if not min_cuts:
                min_cuts = param
            elif not max_cuts:
                max_cuts = param if param >= min_cuts else min_cuts
                min_cuts = param if param <= min_cuts else min_cuts
            elif param > max_cuts:
                max_cuts = param
            elif param < min_cuts:
                min_cuts = param

        elif param in ['alpha', 'position']:
            order = param
        elif "circular".startswith(param):
            topology = "circular"
        elif "linear".startswith(param):
            topology = "linear"
        else:
            _enzymes.append(param)

    _enzymes = ["commercial"] if len(_enzymes) == 0 else _enzymes
    max_cuts = min_cuts if min_cuts and not max_cuts else max_cuts
    min_cuts = 1 if not min_cuts else min_cuts
    return _enzymes, order, min_cuts, max_cuts, topology


# ################################################ MAIN API FUNCTIONS ################################################ #
def amend_metadata(seqbuddy, attr, sub_value, regex):
    """
    Modify, delete, or set metadata in SeqRecord objects
    :param seqbuddy: SeqBuddy object
    :param attr: The attribute to modify. See below for list of supported attributes
    :param sub_value: The string to replace the regex match with
    :param regex: The regular expression that will be fed into re.sub()

    Supported attr types:
        Strings: description, topology, data_file_division, date, source, organism, comment
        OrderedDict: structured_comment
        Lists: accessions, keywords, taxonomy, references, dbxrefs
        Ints: sequence_version
    """
    regex = "|".join(br.clean_regex(regex))
    for rec in seqbuddy.records:
        if attr in ["description", "definition"]:
            rec.description = re.sub(regex, sub_value, rec.description, flags=re.DOTALL)
        elif attr == "topology":
            if sub_value not in ['', 'linear', 'circular']:
                raise ValueError("Topology values are limited to ['', 'linear', 'circular']")
            rec.annotations.setdefault(attr, " ")
            rec.annotations[attr] = sub_value
        elif attr in ["data_file_division", "date", "source", "organism", "comment"]:
            rec.annotations.setdefault(attr, " ")
            rec.annotations[attr] = re.sub(regex, sub_value, rec.annotations[attr], flags=re.DOTALL)
            if attr == "comment" and "structured_comment" in rec.annotations:
                new_structured_comment = OrderedDict()
                for outer_dict_key, inner_dict in rec.annotations["structured_comment"].items():
                    new_key = re.sub(regex, sub_value, outer_dict_key, flags=re.DOTALL)
                    if new_key:
                        new_structured_comment[new_key] = OrderedDict()
                        for inner_dict_key, value in inner_dict.items():
                            new_inner_key = re.sub(regex, sub_value, inner_dict_key, flags=re.DOTALL)
                            if new_inner_key:
                                new_structured_comment[new_key][new_inner_key] = re.sub(regex, sub_value, value,
                                                                                        flags=re.DOTALL)
                rec.annotations["structured_comment"] = new_structured_comment
        elif attr == "references":
            if not sub_value and regex == ".+":
                rec.annotations["references"] = []
            elif regex and regex != ".+":
                rec.annotations.setdefault("references", [])
                new_refs = []
                for ref in rec.annotations["references"]:
                    for ref_attrib in [x for x in dir(ref) if not x.startswith("_") and x != "location"]:
                        setattr(ref, ref_attrib, re.sub(regex, sub_value, getattr(ref, ref_attrib), flags=re.DOTALL))

                    new_refs.append(ref)
                rec.annotations["references"] = new_refs
        elif attr in ["taxonomy", "accessions", "keywords"]:
            # Check to see if the user is trying to send in a whole new list
            if type(sub_value) == str:
                sub_value = sub_value.split() if sub_value else [""]
            if len(sub_value) == 1 and regex != ".+":
                sub_value = sub_value[0]
                rec.annotations.setdefault(attr, [])
                if not rec.annotations[attr]:
                    new_list = [sub_value]
                else:
                    new_list = []
                    for value in rec.annotations[attr]:
                        new_list.append(re.sub(regex, sub_value, value, flags=re.DOTALL))
                new_list = [x for x in new_list if x]
            else:
                new_list = sub_value
            rec.annotations[attr] = new_list
        elif attr in ["dbxrefs", "dblink"]:
            # Check to see if the user is trying to send in a whole new list
            if type(sub_value) == str:
                sub_value = sub_value.split() if sub_value else [""]
            if len(sub_value) == 1 and regex != ".+":
                sub_value = sub_value[0]
                if not rec.dbxrefs:
                    new_list = [sub_value]
                else:
                    new_list = []
                    for value in rec.dbxrefs:
                        new_list.append(re.sub(regex, sub_value, value, flags=re.DOTALL))
                new_list = [x for x in new_list if x]
            else:
                new_list = sub_value
            rec.dbxrefs = new_list
        elif attr in ["version", "sequence_version"]:
            try:
                rec.annotations[attr] = int(sub_value)
                # rec.id must be updated manually, because it is built during rec instantiation
                rec_id_split = rec.id.split(".")
                rec.id = ".".join(rec_id_split[:-1]) if len(rec_id_split) > 1 else rec_id_split[0]
                rec.id += ".%i" % rec.annotations[attr]
            except ValueError:
                pass
        elif attr in dir(rec):
            if type(getattr(rec, attr)) == str:
                setattr(rec, attr, re.sub(regex, sub_value, getattr(rec, attr)))
        else:
            setattr(rec, attr, sub_value)
    return seqbuddy


def annotate(seqbuddy, _type, location, strand=None, qualifiers=None, pattern=None):
    """
    Adds a feature annotation to sequences in the SeqBuddy object
    :param seqbuddy: SeqBuddy object
    :param _type: Name/designation of the annotation
    :param location: The location of the new annotation in the sequence.
    If a single SeqFeature, use a tuple (start, end) or FeatureLocation object
    If a CompoundFeature, us a list of tuples [(start1, end1), (start2, end2)] or CompoundFeature object
    NOTE!!! If feeding in tuples, the 'start' index begins at 1, while Feature objects start at 0.
    :type location: str list tuple FeatureLocation CompoundFeature
    :param strand: The feature's orientation (+/-/None)
    :param qualifiers: Further information to append to the new feature
    The argument can be a dictionary or a list ["foo: bar", "fizz: buzz"]
    :param pattern: List of regex patterns to specify which sequences to add the feature to
    :return: The updated SeqBuddy object
    """
    # http://www.insdc.org/files/feature_table.html
    old = make_copy(seqbuddy)
    if pattern:
        recs = pull_recs(seqbuddy, pattern).records
    else:
        recs = seqbuddy.records

    for rec1 in recs:
        while True:
            breakout = True
            for indx2, rec2 in enumerate(old.records):
                if rec1.id == rec2.id:
                    del old.records[indx2]
                    breakout = False
                    break
            if breakout:
                break

    if isinstance(location, list) or isinstance(location, tuple):
        locations = []
        if isinstance(location[0], int):
            locations.append(FeatureLocation(start=location[0] - 1, end=location[1]))
        elif isinstance(location[0], tuple) or isinstance(location[0], list):
            for _tup in location:
                locations.append(FeatureLocation(start=_tup[0] - 1, end=_tup[1]))
        elif isinstance(location[0], str):
            try:
                for substr in location:
                    substr = re.sub('[ ()]', '', substr)
                    substr = re.sub(r'-|\.\.', ',', substr)
                    locations.append(FeatureLocation(start=int(re.split(',', substr)[0]) - 1,
                                                     end=int(re.split(',', substr)[1])))
            except ValueError:
                raise AttributeError("The provided location string is invalid")
        location = CompoundLocation(sorted(locations, key=lambda x: x.start), operator='order') \
            if len(locations) > 1 else locations[0]
    elif isinstance(location, str):
        try:
            location = re.sub('[ ()]', '', location)
            location = re.split(',', location)
            locations = []
            for substr in location:
                locations.append(FeatureLocation(start=int(substr.split('-')[0]) - 1, end=int(substr.split('-')[1])))
            location = CompoundLocation(sorted(locations, key=lambda x: x.start), operator='order') \
                if len(locations) > 1 else locations[0]
        except ValueError:
                raise AttributeError("The provided location string is invalid")
    elif isinstance(location, FeatureLocation) or isinstance(location, CompoundLocation):
        pass
    else:
        raise TypeError("Input must be list, tuple, or string. Not {0}.".format(type(location)))

    if location.start < 0:
        if isinstance(location, FeatureLocation):
            location = FeatureLocation(0, location.end, location.strand)
        else:
            first_part = FeatureLocation(0, location.parts[0].end, location.parts[0].strand)
            location = CompoundLocation([first_part] + location.parts[1:], operator='order')

    if seqbuddy.alpha == "DNA":
        if strand in ['+', 'plus', 'sense', 'pos', 'positive', '1', 1]:
            strand = 1
        elif strand in ['-', 'minus', 'anti', 'antisense', 'anti-sense', 'neg', 'negative', '-1', -1]:
            strand = -1
        elif strand in ['0', 0]:
            strand = 0
        elif strand is None:
            pass
        else:
            strand = None
            br._stderr("Warning: strand input not recognized. Value set to None.")
    else:
        strand = None
    qualifiers = [qualifiers] if isinstance(qualifiers, str) else qualifiers

    if isinstance(qualifiers, list):
        qual_dict = OrderedDict()
        for qual in qualifiers:
            qual = qual.split("=")
            qual_dict[qual[0]] = "=".join(qual[1:])
        qualifiers = qual_dict

    elif qualifiers and not isinstance(qualifiers, dict):
        raise TypeError("Qualifiers must a list or dict")

    for _rec in recs:
        copy_loc = deepcopy(location)
        if len(_rec) < location.end:
            if isinstance(location, FeatureLocation):
                copy_loc = FeatureLocation(start=location.start, end=len(_rec), strand=location.strand)
            else:  # must be CompoundLocation
                indx = 0
                for indx, part in enumerate(location.parts):
                    if part.end > len(_rec):
                        break
                if indx > 0:
                    copy_loc = CompoundLocation(location.parts[:indx + 1], operator='order')
                    copy_loc.parts[-1] = FeatureLocation(start=location.parts[-1].start, end=len(_rec),
                                                         strand=location.strand)
                else:
                    copy_loc = FeatureLocation(start=location.start, end=len(_rec), strand=location.strand)

        _rec.features.append(SeqFeature(location=copy_loc, type=_type, strand=strand, qualifiers=qualifiers))

    seqbuddy.records = recs
    order_features_by_position(seqbuddy)
    seqbuddy = merge(old, seqbuddy)
    return seqbuddy


def ave_seq_length(seqbuddy, clean=False):
    """
    Calculate the the average length of sequences in the seqbuddy object
    :param seqbuddy: SeqBuddy object
    :param clean: Specifies if non-sequence characters should be counted as well.
    :return: average sequence length (float)
    """
    if clean:  # Strip out all gaps and stuff before counting
        clean_seq(seqbuddy)

    sum_length = 0.
    for rec in seqbuddy.records:
        sum_length += len(rec.seq)
    return sum_length / len(seqbuddy)


def back_translate(seqbuddy, mode='random', species=None, r_seed=None):
    """
    Back-translates protein sequences into DNA sequences
    :param seqbuddy: SeqBuddy object
    :param mode: The codon selection mode (random/optimized)
    :param species: The model to use for optimized codon selection (human/mouse/yeast/ecoli)
    codon preference tables derived from the data at http://www.kazusa.or.jp
    :param r_seed: Set the random generator seed value
    :return: Modified SeqBuddy object
    """
    rand_gen = Random(r_seed)

    # Homo sapiens, species=9606
    if mode.upper() not in ['RANDOM', 'R', 'OPTIMIZED', 'O']:
        raise AttributeError("Back_translate modes accepted are 'random' or 'r' and 'optimized' or 'o'. "
                             "You entered '%s'" % mode)

    h_sapi = {'A': (['GCT', 'GCC', 'GCA', 'GCG'], [0.2675, 0.3975, 0.2275, 0.1075]),
              'C': (['TGT', 'TGC'], [0.46, 0.54]),
              'D': (['GAT', 'GAC'], [0.46, 0.54]),
              'E': (['GAA', 'GAG'], [0.42, 0.58]),
              'F': (['TTT', 'TTC'], [0.46, 0.54]),
              'G': (['GGT', 'GGC', 'GGA', 'GGG'], [0.16, 0.34, 0.25, 0.25]),
              'H': (['CAT', 'CAC'], [0.42, 0.58]),
              'I': (['ATT', 'ATC', 'ATA'], [0.36, 0.47, 0.17]),
              'K': (['AAA', 'AAG'], [0.43, 0.57]),
              'L': (['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], [0.08, 0.1275, 0.1275, 0.1975, 0.07, 0.3975]),
              'M': (['ATG'], [1.00]),
              'N': (['AAT', 'AAC'], [0.47, 0.53]),
              'P': (['CCT', 'CCC', 'CCA', 'CCG'], [0.29, 0.32, 0.28, 0.11]),
              'Q': (['CAA', 'CAG'], [0.27, 0.73]),
              'R': (['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], [0.08, 0.18, 0.1125, 0.2025, 0.2125, 0.2125]),
              'S': (['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], [0.19, 0.22, 0.15, 0.05, 0.15, 0.24]),
              '*': (['TAA', 'TGA', 'TAG'], [0.30, 0.46, 0.24]),
              'T': (['ACT', 'ACC', 'ACA', 'ACG'], [0.25, 0.36, 0.28, 0.11]),
              'V': (['GTT', 'GTC', 'GTA', 'GTG'], [0.18, 0.24, 0.12, 0.46]),
              'W': (['TGG'], [1.00]),
              'Y': (['TAT', 'TAC'], [0.44, 0.56]),
              'X': (['NNN'], [1.0]),
              '-': (['---'], [1.0])}

    # Mus musculus, species=10090
    m_muscul = {'A': (['GCT', 'GCC', 'GCA', 'GCG'], [0.2925, 0.3825, 0.2325, 0.0925]),
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
                'T': (['ACT', 'ACC', 'ACA', 'ACG'], [0.2525, 0.3525, 0.2925, 0.1025]),
                'V': (['GTT', 'GTC', 'GTA', 'GTG'], [0.17, 0.25, 0.12, 0.46]),
                'W': (['TGG'], [1.00]),
                'Y': (['TAT', 'TAC'], [0.43, 0.57]),
                'X': (['NNN'], [1.0]),
                '-': (['---'], [1.0])}

    # Escherichia coli O157:H7 EDL933, species=155864
    e_coli = {'A': (['GCT', 'GCC', 'GCA', 'GCG'], [0.16, 0.27, 0.22, 0.35]),
              'C': (['TGT', 'TGC'], [0.45, 0.55]),
              'D': (['GAT', 'GAC'], [0.63, 0.37]),
              'E': (['GAA', 'GAG'], [0.68, 0.32]),
              'F': (['TTT', 'TTC'], [0.58, 0.42]),
              'G': (['GGT', 'GGC', 'GGA', 'GGG'], [0.33, 0.39, 0.12, 0.16]),
              'H': (['CAT', 'CAC'], [0.58, 0.42]),
              'I': (['ATT', 'ATC', 'ATA'], [0.505, 0.404, 0.091]),
              'K': (['AAA', 'AAG'], [0.76, 0.24]),
              'L': (['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], [0.13, 0.13, 0.11, 0.10, 0.04, 0.49]),
              'M': (['ATG'], [1.00]),
              'N': (['AAT', 'AAC'], [0.47, 0.53]),
              'P': (['CCT', 'CCC', 'CCA', 'CCG'], [0.17, 0.13, 0.19, 0.51]),
              'Q': (['CAA', 'CAG'], [0.33, 0.67]),
              'R': (['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], [0.3625, 0.3725, 0.0725, 0.1125, 0.05, 0.03]),
              'S': (['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], [0.14, 0.1475, 0.14, 0.1475, 0.1575, 0.2675]),
              '*': (['TAA', 'TGA', 'TAG'], [0.59, 0.33, 0.08]),
              'T': (['ACT', 'ACC', 'ACA', 'ACG'], [0.17, 0.41, 0.15, 0.27]),
              'V': (['GTT', 'GTC', 'GTA', 'GTG'], [0.26, 0.21, 0.16, 0.37]),
              'W': (['TGG'], [1.00]),
              'Y': (['TAT', 'TAC'], [0.57, 0.43]),
              'X': (['NNN'], [1.0]),
              '-': (['---'], [1.0])}

    # Saccharomyces cerevisiae, species=4932
    s_cerev = {'A': (['GCT', 'GCC', 'GCA', 'GCG'], [0.38, 0.22, 0.29, 0.11]),
               'C': (['TGT', 'TGC'], [0.63, 0.37]),
               'D': (['GAT', 'GAC'], [0.65, 0.35]),
               'E': (['GAA', 'GAG'], [0.70, 0.30]),
               'F': (['TTT', 'TTC'], [0.59, 0.41]),
               'G': (['GGT', 'GGC', 'GGA', 'GGG'], [0.47, 0.19, 0.22, 0.12]),
               'H': (['CAT', 'CAC'], [0.64, 0.36]),
               'I': (['ATT', 'ATC', 'ATA'], [0.464, 0.263, 0.273]),
               'K': (['AAA', 'AAG'], [0.58, 0.42]),
               'L': (['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], [0.2775, 0.2875, 0.1275, 0.06, 0.1375, 0.11]),
               'M': (['ATG'], [1.00]),
               'N': (['AAT', 'AAC'], [0.59, 0.41]),
               'P': (['CCT', 'CCC', 'CCA', 'CCG'], [0.31, 0.15, 0.42, 0.12]),
               'Q': (['CAA', 'CAG'], [0.69, 0.31]),
               'R': (['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], [0.14, 0.06, 0.07, 0.04, 0.48, 0.21]),
               'S': (['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], [0.26, 0.16, 0.21, 0.10, 0.16, 0.11]),
               '*': (['TAA', 'TGA', 'TAG'], [0.47, 0.30, 0.23]),
               'T': (['ACT', 'ACC', 'ACA', 'ACG'], [0.3475, 0.2175, 0.2975, 0.1375]),
               'V': (['GTT', 'GTC', 'GTA', 'GTG'], [0.39, 0.21, 0.21, 0.19]),
               'W': (['TGG'], [1.00]),
               'Y': (['TAT', 'TAC'], [0.56, 0.44]),
               'X': (['NNN'], [1.0]),
               '-': (['---'], [1.0])}

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
                  'X': (['NNN'], [1.0]),
                  '-': (['---'], [1.0])}

    if seqbuddy.alpha != "protein":
        raise TypeError("The input sequence needs to be 'protein', not '%s'" % str(type(seqbuddy.alpha)))

    if not species:
        lookup_table = rand_table
    elif species.upper() in ["HUMAN", "H"]:
        lookup_table = h_sapi
    elif species.upper() in ["MOUSE", "M"]:
        lookup_table = m_muscul
    elif species.upper() in ["ECOLI", "E"]:
        lookup_table = e_coli
    elif species.upper() in ["YEAST", "Y"]:
        lookup_table = s_cerev
    else:
        raise AttributeError("The species requested does not match any lookup tables currently implemented. "
                             "Please leave blank or select from human, mouse, ecoli, or yeast.")

    # Modify the lookup table so each amino acid maps to a single 'optimal' codon
    if mode.upper() in ['OPTIMIZED', 'O']:
        for aa in lookup_table:
            best = ["", 0.]
            for i in range(len(lookup_table[aa][1])):
                if lookup_table[aa][1][i] > best[1]:
                    best = [lookup_table[aa][0][i], lookup_table[aa][1][i]]
            lookup_table[aa] = ([best[0]], [1.0])

    clean_seq(seqbuddy, skip_list=r"\-*")
    originals = make_copy(seqbuddy)
    for rec in seqbuddy.records:
        rec.features = []
        dna_seq = ["" for _ in range(len(rec))]
        for indx, aa in enumerate(rec.seq.upper()):
            rand_num = rand_gen.random()
            sum_probs = 0.
            for i in range(len(lookup_table[aa][1])):
                sum_probs += lookup_table[aa][1][i]
                if sum_probs >= rand_num:
                    dna_seq[indx] = lookup_table[aa][0][i]
                    break
            rec.seq = Seq("".join(dna_seq))
            rec.annotations["molecule_type"] = "DNA"

    seqbuddy.alpha = "DNA"
    map_features_prot2nucl(originals, seqbuddy, mode="list")
    return seqbuddy


# The following function is used by bl2seq() below
def _mc_blast(_query, _args):
    subject_file, lock, tmp_dir, blast_bin = _args
    subject_sb = SeqBuddy(subject_file)
    assert len(subject_sb) == 1
    _query_file = br.TempFile()
    query_sb = SeqBuddy(_query)
    query_sb.write(_query_file.path, out_format="fasta")

    blast_res = Popen('%s -query "%s" -subject "%s" -outfmt 6' %
                      (blast_bin, _query_file.path, subject_file), stdout=PIPE, shell=True).communicate()
    blast_res = blast_res[0].decode().strip().split("\n")
    result_file = br.TempFile()
    hit_ids = {}
    current_id = ""
    for res in blast_res:
        res = res.strip().split("\t")
        if current_id == res[0] or res[0] == res[1]:
            continue
        current_id = res[0]
        hit_ids[res[0]] = None
        if res[10] == '0.0':
            res[10] = '1e-180'

        res = "%s\t%s\t%s\t%s\t%s\t%s\n" % (res[0], res[1], res[2], res[3], res[10], res[11])
        result_file.write(res)

    for _rec in query_sb.records:
        if _rec.id not in hit_ids and _rec.id != subject_sb.records[0].id:
            result_file.write("%s\t%s\t0\t0\t0\t0\n" % (_rec.id, subject_sb.records[0].id))

    _result = result_file.read()
    with lock:
        with open(join(tmp_dir.path, "blast_results.txt"), "a", encoding="utf-8") as _ofile:
            _ofile.write(_result)
    return


def bl2seq(seqbuddy):
    """
    Does an all-by-all analysis of the sequences
    :param seqbuddy: SeqBuddy object
    :return: OrderedDict of results dict[key][matches]
    """
    # Note on blast2seq: Expect (E) values are calculated on an assumed database size of (the rather large) nr, so the
    # threshold may need to be increased quite a bit to return short alignments
    if seqbuddy.alpha == "protein" and not _check_for_blast_bin("blastp"):
        raise RuntimeError("Blastp not present in $PATH or working directory.")

    elif (seqbuddy.alpha == "DNA" or seqbuddy.alpha == "RNA") \
            and not _check_for_blast_bin("blastn"):
        raise RuntimeError("Blastn not present in $PATH or working directory.")

    blast_bin = "blastp" if seqbuddy.alpha == "protein" else "blastn"

    from multiprocessing import Lock
    lock = Lock()
    tmp_dir = br.TempDir()

    # Remove any gaps
    seqbuddy = clean_seq(seqbuddy, skip_list=["*"])

    # Copy the seqbuddy records into new list, so they can be iteratively deleted below
    make_ids_unique(seqbuddy, sep="-")
    seqs_copy = make_copy(seqbuddy).records
    subject_file = join(tmp_dir.path, "subject.fa")
    cpus = br.usable_cpu_count()
    for subject in seqbuddy.records:
        with open(subject_file, "w", encoding="utf-8") as ifile:
            SeqIO.write(subject, ifile, "fasta")

        br.run_multicore_function(br.chunk_list(seqs_copy, cpus), _mc_blast, [subject_file, lock, tmp_dir, blast_bin],
                                  out_type=sys.stderr, quiet=True)

        seqs_copy = seqs_copy[1:]

    with open(join(tmp_dir.path, "blast_results.txt"), "r", encoding="utf-8") as _ifile:
        output_list = _ifile.read().strip().split("\n")

    # Push output into a dictionary of dictionaries, for more flexible use outside of this function
    output_list = [x.split("\t") for x in output_list]
    output_list = sorted(output_list, key=lambda l: l[0])
    output_dict = {}
    for match in output_list:
        query, subj, ident, length, evalue, bit_score = match
        output_dict.setdefault(query, {})
        output_dict[query][subj] = [float(ident), int(length), float(evalue), float(bit_score)]

        output_dict.setdefault(subj, {})
        output_dict[subj][query] = [float(ident), int(length), float(evalue), float(bit_score)]

    for key, value in output_dict.items():
        output_dict[key] = [(x, y) for x, y in output_dict[key].items()]
        output_dict[key] = OrderedDict(sorted(output_dict[key], key=lambda l: l[0]))

    output_dict = [(key, value) for key, value in output_dict.items()]
    output_dict = OrderedDict(sorted(output_dict, key=lambda l: l[0]))

    return output_dict


def blast(subject, query, **kwargs):
    """
    Runs a BLAST search against a specified database or query SeqBuddy obj, returning all significant matches.
    :param subject: SeqBuddy object
    :param query: Another SeqBuddy object, or the location of the BLAST database to run the sequences against
    :param kwargs:  - quiet -> [True|False]
                    - blast_args -> [<any extra blast command line arguments>]
    :return: A SeqBuddy object containing all of the BLAST database matches
    """
    # ToDo: Allow mixed sequence types (blastx?)
    # ToDo: Implement the 'keep' flag so the Databases can be saved.

    kwargs["quiet"] = False if "quiet" not in kwargs or not kwargs["quiet"] else True

    blast_bin = "blastp" if subject.alpha == "protein" else "blastn"
    if not _check_for_blast_bin(blast_bin):
        raise SystemError("%s not found in system path." % blast_bin)

    if not _check_for_blast_bin("blastdbcmd"):
        raise SystemError("blastdbcmd not found in system path.")

    extensions = {"blastp": ["phr", "pin", "pog", "psd", "psi", "psq"],
                  "blastn": ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]}
    tmp_dir = br.TempDir()

    if query.__class__.__name__ == "SeqBuddy":
        if not _check_for_blast_bin("makeblastdb"):
            raise SystemError("makeblastdb not found in system path.")
        if subject.alpha != query.alpha:
            raise ValueError("Trying to compare protein to nucleotide.")
        query_sb = hash_ids(query, r_seed=12345)
        query_sb = clean_seq(query_sb, skip_list="*")
        query_sb.write(os.path.join(tmp_dir.path, "query.fa"), out_format="fasta")
        dbtype = "prot" if subject.alpha == "protein" else "nucl"
        query_path = os.path.join(tmp_dir.path, "query.fa")
        query_db_path = os.path.join(tmp_dir.path, "query_db")
        makeblastdb = Popen('makeblastdb -dbtype {0} -in "{1}" -out "{2}" '
                            '-parse_seqids'.format(dbtype, query_path, query_db_path), shell=True,
                            stdout=PIPE).communicate()[0].decode("utf-8")
        makeblastdb = re.sub("New DB .*\n", "", makeblastdb.strip())
        makeblastdb = re.sub("Building a new DB", "Building a new DB with makeblastdb", makeblastdb)
        br._stderr("%s\n\n" % makeblastdb, quiet=kwargs["quiet"])
        query = "%s%squery_db" % (tmp_dir.path, os.path.sep)

    else:
        query_sb = None

    # Try to catch the common variations of the database names that might be given as input
    if query[-2:] in [".p", ".n"]:
        query = query[:-2]

    if query[-3:] in extensions[blast_bin]:
        query = query[:-4]

    query = os.path.abspath(query)

    # Check that complete blastdb is present and was made with the -parse_seqids flag
    for extension in extensions[blast_bin]:
        if not os.path.isfile("%s.%s" % (query, extension)):
            raise RuntimeError("The .%s file of your %s database was not found. Ensure the -parse_seqids flag was "
                               "used with makeblastdb." % (extension, blast_bin.upper()))

    subject = clean_seq(subject)  # in case there are gaps or something in the sequences

    with open("%s%stmp.fa" % (tmp_dir.path, os.path.sep), "w", encoding="utf-8") as ofile:
        SeqIO.write(subject.records, ofile, "fasta")

    num_threads = 4
    evalue = 0.01
    if "blast_args" in kwargs and kwargs["blast_args"]:
        # Clean up any parameters that the user may be including that will break SeqBuddy
        black_list = ["db", "query", "subject", "out", "outfmt"]
        for no_no in black_list:
            if re.search("-%s" % no_no, kwargs["blast_args"]):
                br._stderr("Warning: Explicitly setting the blast -%s parameter is not supported in SeqBuddy.\n" %
                           no_no, quiet=kwargs["quiet"])
                kwargs["blast_args"] = re.sub("-%s .*-" % no_no, "-", kwargs["blast_args"])
                kwargs["blast_args"] = re.sub("-%s .*$" % no_no, "", kwargs["blast_args"])

        if "-num_threads" in kwargs["blast_args"]:
            try:
                num_threads_search = re.search("-num_threads ([^-]*)", kwargs["blast_args"])
                num_threads = int(num_threads_search.group(1))
                kwargs["blast_args"] = re.sub("-num_threads .*-", "-", kwargs["blast_args"])
                kwargs["blast_args"] = re.sub("-num_threads .*$", "", kwargs["blast_args"])

            except ValueError:
                raise ValueError("-num_threads expects an integer.")

        if "-evalue" in kwargs["blast_args"]:
            try:
                evalue_search = re.search("-evalue ([^-]*)", kwargs["blast_args"])
                evalue = float(evalue_search.group(1))
                kwargs["blast_args"] = re.sub("-evalue .*-", "-", kwargs["blast_args"])
                kwargs["blast_args"] = re.sub("-evalue .*$", "", kwargs["blast_args"])

            except ValueError:
                raise ValueError("-evalue expects a number.")

    else:
        kwargs["blast_args"] = ""

    blast_command = '{0} -db "{1}" -query {2}tmp.fa -out {2}out.txt ' \
                    '-outfmt 6 -num_threads {3} -evalue {4} {5}'.format(blast_bin, query, tmp_dir.path + os.path.sep,
                                                                        num_threads, evalue, kwargs["blast_args"])

    br._stderr("Running...\n%s\n\n" %
               re.sub("-db.*-num_threads", "-num_threads", blast_command), quiet=kwargs["quiet"])
    blast_output = Popen(blast_command, shell=True, stdout=PIPE, stderr=PIPE).communicate()
    blast_output = blast_output[1].decode("utf-8")

    if "Error" in blast_output:
        raise RuntimeError(blast_output)

    with open("%s%sout.txt" % (tmp_dir.path, os.path.sep), "r", encoding="utf-8") as ifile:
        blast_results = ifile.read()
        records = blast_results.split("\n")

    hit_ids = []
    for record in records:
        record = record.split("\t")
        if len(record) == 1:
            continue
        hit_id = record[1].strip()
        if hit_id in hit_ids:
            continue

        hit_ids.append(hit_id)

    with open("%s%sseqs.fa" % (tmp_dir.path, os.path.sep), "w", encoding="utf-8") as ofile:
        for hit_id in hit_ids:
            hit = Popen('blastdbcmd -db "%s" -entry \"lcl|%s\"' % (query, hit_id),
                        stdout=PIPE, shell=True).communicate()
            hit = hit[0].decode("utf-8")
            hit = re.sub(r"lcl\|", "", hit)
            ofile.write("%s\n" % hit)

    new_seqs = SeqBuddy("%s%sseqs.fa" % (tmp_dir.path, os.path.sep))
    new_seqs.out_format = subject.out_format
    if query_sb:
        new_seqs.hash_map = query_sb.hash_map
        new_seqs.reverse_hashmap()

        for key, value in new_seqs.hash_map.items():
            blast_results = re.sub(key, value, blast_results)

    br._stderr("# ######################## BLAST results ######################## #\n%s"
               "# ############################################################### #\n\n" % blast_results,
               quiet=kwargs["quiet"])
    return new_seqs


def clean_seq(seqbuddy, ambiguous=True, rep_char="N", skip_list=None):
    """
    Removes all non-sequence characters, and converts ambiguous characters to 'X' if ambiguous=False
    :param seqbuddy: SeqBuddy object
    :param rep_char: What character should be used to replace ambiguous characters
    :param ambiguous: Specifies whether ambiguous characters should be kept or not
    :param skip_list: Optional list of characters to be left alone
    :return: The cleaned SeqBuddy object
    """
    seqbuddy_copy = make_copy(seqbuddy)
    skip_list = "" if not skip_list else "".join(skip_list)
    for rec, rec_copy in zip(seqbuddy.records, seqbuddy_copy.records):
        if rec.annotations["molecule_type"] == "protein":
            full_skip = "ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy%s" % skip_list
            rec_copy.seq = Seq(re.sub("[^%s]" % full_skip, "-", str(rec.seq)))
            rec_copy.annotations["molecule_type"] = "protein"
            rec.letter_annotations = {}
            rec.seq = Seq(re.sub("[^%s]" % full_skip, "", str(rec.seq)))
            rec.annotations["molecule_type"] = "protein"
        else:
            full_skip = "ATGCURYWSMKHBVDNXatgcurywsmkhbvdnx%s" % skip_list
            rec_copy.seq = Seq(re.sub("[^%s]" % full_skip, "-", str(rec.seq)))
            rec_copy.annotations["molecule_type"] = rec.annotations["molecule_type"]
            rec.seq = Seq(re.sub("[^%s]" % full_skip, "", str(rec.seq)))
            if not ambiguous:
                full_skip = "ATGCUatgcu%s" % skip_list
                rec.seq = Seq(re.sub("[^%s]" % full_skip, rep_char, str(rec.seq)))
                rec_copy.seq = Seq(re.sub("[^%s]" % full_skip, rep_char, str(rec.seq)))
                rec_copy.annotations["molecule_type"] = rec.annotations["molecule_type"]

    br.remap_gapped_features(seqbuddy_copy.records, seqbuddy.records)
    return seqbuddy


def complement(seqbuddy):
    """
    Converts DNA/RNA sequences to their complementary sequence
    :param seqbuddy: SeqBuddy object
    :return: The modified SeqBuddy object
    """
    if seqbuddy.alpha == "protein":
        raise TypeError("Nucleic acid sequence required, not protein.")
    for rec in seqbuddy.records:
        rec.seq = rec.seq.complement()
    return seqbuddy


def concat_seqs(seqbuddy, clean=False):
    """
    Concatenates sequences into a single record
    :param seqbuddy: SeqBuddy object
    :param clean: Specifies whether non-sequence characters should be removed
    :return: modified SeqBuddy object
    """
    if clean:
        clean_seq(seqbuddy)

    new_seq = ""
    concat_ids = []
    features = []
    letter_annotations = {}
    for rec in seqbuddy.records:
        shift = len(new_seq)
        full_seq_len = len(new_seq) + len(str(rec.seq))
        rec.features = br.shift_features(rec.features, shift, full_seq_len)

        location = FeatureLocation(len(new_seq), len(new_seq) + len(str(rec.seq)))
        feature = SeqFeature(location=location, id=rec.id, type=rec.id[:15])
        features.append(feature)
        features += rec.features

        for qual_type, scores in rec.letter_annotations.items():
            letter_annotations.setdefault(qual_type, [])
            letter_annotations[qual_type] += scores

        concat_ids.append(rec.id)
        new_seq += str(rec.seq)

    new_seq = [SeqRecord(Seq(new_seq),
                         description="", id="concatination", name="concatination", features=features,
                         annotations={"molecule_type": seqbuddy.alpha}, letter_annotations=letter_annotations)]
    seqbuddy.records = new_seq
    seqbuddy.out_format = "gb"
    return seqbuddy


def count_codons(seqbuddy):
    """
    Generate frequency statistics for codon composition
    :param seqbuddy: SeqBuddy object
    :return: A tuple containing the original SeqBuddy object and a dictionary - dict[id][codon] = (Amino acid, num, %)
    """
    if seqbuddy.alpha not in ["RNA", "DNA"]:
        raise TypeError("Nucleic acid sequence required, not protein or other.")
    if seqbuddy.alpha == "DNA":
        codontable = CodonTable.ambiguous_dna_by_name['Standard'].forward_table
    else:
        codontable = CodonTable.ambiguous_rna_by_name['Standard'].forward_table
    output = OrderedDict()
    sb_copy = make_copy(seqbuddy)
    replace_subsequence(sb_copy, "[-.]")
    uppercase(sb_copy)
    for rec in sb_copy.records:
        sequence = str(rec.seq)
        if len(sequence) % 3 != 0:
            while len(sequence) % 3 != 0:
                sequence = sequence[:-1]
        data_table = {}
        num_codons = len(sequence) / 3
        while len(sequence) > 0:
            codon = str(sequence[:3])
            if codon in data_table:
                data_table[codon][1] += 1
            else:
                if codon in ['ATG', 'AUG']:
                    data_table[codon] = ['M', 1, 0.0]
                elif codon == 'NNN':
                    data_table[codon] = ['X', 1, 0.0]
                elif codon in ['TAA', 'TAG', 'TGA', 'UAA', 'UAG', 'UGA']:
                    data_table[codon] = ['*', 1, 0.0]
                else:
                    try:
                        data_table[codon] = [codontable[codon], 1, 0.0]
                    except (KeyError, CodonTable.TranslationError):
                        br._stderr("Warning: Codon '{0}' is invalid. Codon will be skipped.\n".format(codon))
            sequence = sequence[3:]
        for codon in data_table:
            data_table[codon][2] = round(data_table[codon][1] / float(num_codons) * 100, 3)
        output[rec.id] = OrderedDict(sorted(data_table.items(), key=lambda x: x[0]))
    for rec in seqbuddy.records:
        try:
            rec.buddy_data['Codon_frequency'] = output[rec.id]
        except AttributeError:
            rec.buddy_data = OrderedDict({'Codon_frequency': output[rec.id]})
    return seqbuddy, output


def count_residues(seqbuddy):
    """
    Generate frequency statistics for residue composition
    :param seqbuddy: SeqBuddy object
    :return: annotated SeqBuddy object. Residue counts are appended to buddy_data in the SeqRecord obects
    """
    for rec in seqbuddy.records:
        resid_count = {}
        seq = str(rec.seq).upper()
        for char in seq:
            resid_count.setdefault(char, 0)
            resid_count[char] += 1

        seq_len = len(rec)
        for residue, count in resid_count.items():
            resid_count[residue] = [count, count / seq_len]

        if seqbuddy.alpha == "protein":
            ambig = len(re.findall("X", seq))
            if ambig > 0:
                resid_count['% Ambiguous'] = round(100 * ambig / seq_len, 2)

            pos = len(re.findall("[HKR]", seq))
            resid_count['% Positive'] = round(100 * pos / seq_len, 2)

            neg = len(re.findall("[DEC]", seq))
            resid_count['% Negative'] = round(100 * neg / seq_len, 2)

            neut = len(re.findall("[GAVLIPFYWSTNQM]", seq))
            resid_count['% Uncharged'] = round(100 * neut / seq_len, 2)

            hydrophobic = len(re.findall("[AVLIPYFWMC]", seq))
            resid_count['% Hydrophobic'] = round(100 * hydrophobic / seq_len, 2)

            hydrophilic = len(re.findall("[NQSTKRHDE]", seq))
            resid_count['% Hydrophilic'] = round(100 * hydrophilic / seq_len, 2)

            for residue in ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
                            "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]:
                resid_count.setdefault(residue, [0, 0])

        else:
            ambig = len(re.findall("[^ATCGU]", seq))
            if ambig > 0:
                resid_count['% Ambiguous'] = round(100 * ambig / seq_len, 2)

            for residue in ["A", "G", "C"]:
                resid_count.setdefault(residue, [0, 0])

            if "T" not in resid_count and "U" not in resid_count:
                resid_count["T"] = [0, 0]

        _add_buddy_data(rec, "res_count", OrderedDict(sorted(resid_count.items())))
    return seqbuddy


def degenerate_sequence(seqbuddy, table=1):
    """
    Generate degenerate codon sequence
    :param seqbuddy: The SeqBuddy object to be analyzed
    :param table: The degenerate codon table to use
    :return: A SeqBuddy object containing a degenerate nucleotide sequence

    Contributed by Jeremy Labarge (https://github.com/biojerm)

    The method is developed based on the Perl script Degen v1.4.
    http://www.phylotools.com/ptdegenoverview.htm

    Other inspiration
    https://github.com/carlosp420/degenerate-dna

    Zwick, A., Regier, J.C. & Zwickl, D.J. (2012). "Resolving Discrepancy between Nucleotides and Amino Acids in
    Deep-Level Arthropod Phylogenomics: Differentiating Serine Codons in 21-Amino-Acid Models". PLoS ONE 7(11): e47450.

    Regier, J.C., Shultz, J.W., Zwick, A., Hussey, A., Ball, B., Wetzer, R. Martin, J.W. & Cunningham, C.W. (2010).
     "Arthropod relationships revealed by phylogenomic analysis of nuclear protein-coding sequences".
     Nature 463: 1079-1083.
    """
    dgndict = {}
    # degenerate codons that are  universal to each all codon tables
    base_dict = {'ACA': 'ACN', 'ACC': 'ACN', 'ACG': 'ACN', 'ACT': 'ACN', 'CAC': 'CAY', 'CAT': 'CAY', 'CCA': 'CCN',
                 'CCC': 'CCN', 'CCG': 'CCN', 'CCT': 'CCN', 'GAA': 'GAR', 'GAC': 'GAY', 'GAG': 'GAR', 'GAT': 'GAY',
                 'GCA': 'GCN', 'GCC': 'GCN', 'GCG': 'GCN', 'GCT': 'GCN', 'GTA': 'GTN', 'GTC': 'GTN', 'GTG': 'GTN',
                 'GTT': 'GTN', 'TCA': 'TCN', 'TCC': 'TCN', 'TCG': 'TCN', 'TCT': 'TCN', 'TTC': 'TTY', 'TTT': 'TTY',
                 'NNN': 'NNN', '???': 'NNN', '---': '---'}

    # Standard Genetic Code codons.  It is the default table
    dgndict[1] = {'AAA': 'AAR', 'AAC': 'AAY', 'AAG': 'AAR', 'AAT': 'AAY', 'AGA': 'MGN', 'AGC': 'TCN', 'AGG': 'MGN',
                  'AGT': 'AGY', 'ATA': 'ATH', 'ATC': 'ATH', 'ATG': 'ATG', 'ATT': 'ATH', 'CAA': 'CAR', 'CAG': 'CAR',
                  'CGA': 'MGN', 'CGC': 'MGN', 'CGG': 'MGN', 'CGT': 'MGN', 'CTA': 'YTN', 'CTC': 'YTN', 'CTG': 'YTN',
                  'CTT': 'YTN', 'GGA': 'GGN', 'GGC': 'GGN', 'GGG': 'GGN', 'GGT': 'GGN', 'TAA': 'TAA', 'TAC': 'TAY',
                  'TAG': 'TAG', 'TAT': 'TAY', 'TGA': 'TGA', 'TGC': 'TGY', 'TGG': 'TGG', 'TGT': 'TGY', 'TTA': 'YTN',
                  'TTG': 'YTN'}

    # Vertebrate Mitochondrial Code codons
    dgndict[2] = {'AAA': 'AAR', 'AAC': 'AAY', 'AAG': 'AAR', 'AAT': 'AAY', 'ATA': 'ATR', 'ATC': 'ATY', 'ATG': 'ATR',
                  'ATT': 'ATY', 'AGA': 'AGA', 'AGG': 'AGG', 'AGC': 'AGY', 'AGT': 'AGY', 'CAA': 'CAR', 'CAG': 'CAR',
                  'CGA': 'CGN', 'CGC': 'CGN', 'CGG': 'CGN', 'CGT': 'CGN', 'CTA': 'YTN', 'CTC': 'YTN', 'CTG': 'YTN',
                  'CTT': 'YTN', 'GGA': 'GGN', 'GGC': 'GGN', 'GGG': 'GGN', 'GGT': 'GGN', 'TAA': 'TAA', 'TAC': 'TAY',
                  'TAG': 'TAG', 'TAT': 'TAY', 'TGA': 'TGR', 'TGC': 'TGY', 'TGG': 'TGR', 'TGT': 'TGY', 'TTA': 'YTN',
                  'TTG': 'YTN'}

    # Yeast Mitochondrial Code codons
    dgndict[3] = {'AAA': 'AAR', 'AAC': 'AAY', 'AAG': 'AAR', 'AAT': 'AAY', 'AGA': 'MGN', 'AGC': 'AGY', 'AGG': 'MGN',
                  'AGT': 'AGY', 'ATA': 'ATR', 'ATC': 'ATY', 'ATG': 'ATR', 'ATT': 'ATY', 'CAA': 'CAR', 'CAG': 'CAR',
                  'CGA': 'CGA', 'CGC': 'CGC', 'CGG': 'MGN', 'CGT': 'MGN', 'CTA': 'CTN', 'CTC': 'CTN', 'CTG': 'CTN',
                  'CTT': 'CTN', 'GGA': 'GGN', 'GGC': 'GGN', 'GGG': 'GGN', 'GGT': 'GGN', 'TAA': 'TAA', 'TAC': 'TAY',
                  'TAG': 'TAG', 'TAT': 'TAY', 'TGA': 'TGR', 'TGC': 'TGY', 'TGG': 'TGR', 'TGT': 'TGY', 'TTA': 'TTR',
                  'TTG': 'TTR'}

    # Mold / Protozoan / Coelenterate Mitochondrial Code & Mycoplasma / Spiroplasma Code Codons
    dgndict[4] = {'AAA': 'AAR', 'AAC': 'AAY', 'AAG': 'AAR', 'AAT': 'AAY', 'AGA': 'MGN', 'AGC': 'AGY', 'AGG': 'MGN',
                  'AGT': 'AGY', 'ATA': 'ATH', 'ATC': 'ATH', 'ATG': 'ATG', 'ATT': 'ATH', 'CAA': 'CAR', 'CAG': 'CAR',
                  'CGA': 'MGN', 'CGC': 'MGN', 'CGG': 'MGN', 'CGT': 'MGN', 'CTA': 'YTN', 'CTC': 'YTN', 'CTG': 'YTN',
                  'CTT': 'YTN', 'GGA': 'GGN', 'GGC': 'GGN', 'GGG': 'GGN', 'GGT': 'GGN', 'TAA': 'TAA', 'TAC': 'TAY',
                  'TAG': 'TAG', 'TAT': 'TAY', 'TGA': 'TGR', 'TGC': 'TGY', 'TGG': 'TGR', 'TGT': 'TGY', 'TTA': 'YTN',
                  'TTG': 'YTN'}

    # Invertebrate Mitochondrial Code codons
    dgndict[5] = {'AAA': 'AAR', 'AAC': 'AAY', 'AAG': 'AAR', 'AAT': 'AAY', 'AGA': 'AGN', 'AGC': 'AGN', 'AGG': 'AGN',
                  'AGT': 'AGN', 'ATA': 'ATR', 'ATC': 'ATY', 'ATG': 'ATR', 'ATT': 'ATY', 'CAA': 'CAR', 'CAG': 'CAR',
                  'CGA': 'CGN', 'CGC': 'CGN', 'CGG': 'CGN', 'CGT': 'CGN', 'CTA': 'YTN', 'CTC': 'YTN', 'CTG': 'YTN',
                  'CTT': 'YTN', 'GGA': 'GGN', 'GGC': 'GGN', 'GGG': 'GGN', 'GGT': 'GGN', 'TAA': 'TAA', 'TAC': 'TAY',
                  'TAG': 'TAG', 'TAT': 'TAY', 'TGA': 'TGR', 'TGC': 'TGY', 'TGG': 'TGR', 'TGT': 'TGY', 'TTA': 'YTN',
                  'TTG': 'YTN'}

    # Ciliate / Dasycladacean / Hexamita Nuclear Code codons
    dgndict[6] = {'AAA': 'AAR', 'AAC': 'AAY', 'AAG': 'AAR', 'AAT': 'AAY', 'AGA': 'MGN', 'AGC': 'AGY', 'AGG': 'MGN',
                  'AGT': 'AGY', 'ATA': 'ATH', 'ATC': 'ATH', 'ATG': 'ATG', 'ATT': 'ATH', 'CAA': 'YAR', 'CAG': 'YAR',
                  'CGA': 'MGN', 'CGC': 'MGN', 'CGG': 'MGN', 'CGT': 'MGN', 'CTA': 'YTN', 'CTC': 'YTN', 'CTG': 'YTN',
                  'CTT': 'YTN', 'GGA': 'GGN', 'GGC': 'GGN', 'GGG': 'GGN', 'GGT': 'GGN', 'TAA': 'YAR', 'TAC': 'TAY',
                  'TAG': 'YAR', 'TAT': 'TAY', 'TGA': 'TGA', 'TGC': 'TGY', 'TGG': 'TGG', 'TGT': 'TGY', 'TTA': 'YTN',
                  'TTG': 'YTN'}

    # Echinoderm and Flatworm Mitochondrial Code codons
    dgndict[7] = {'AAA': 'AAH', 'AAC': 'AAH', 'AAG': 'AAG', 'AAT': 'AAH', 'AGA': 'AGN', 'AGC': 'AGN', 'AGG': 'AGN',
                  'AGT': 'AGN', 'ATA': 'ATH', 'ATC': 'ATH', 'ATG': 'ATG', 'ATT': 'ATH', 'CAA': 'CAR', 'CAG': 'CAR',
                  'CGA': 'CGN', 'CGC': 'CGN', 'CGG': 'CGN', 'CGT': 'CGN', 'CTA': 'YTN', 'CTC': 'YTN', 'CTG': 'YTN',
                  'CTT': 'YTN', 'GGA': 'GGN', 'GGC': 'GGN', 'GGG': 'GGN', 'GGT': 'GGN', 'TAA': 'TAA', 'TAC': 'TAY',
                  'TAG': 'TAG', 'TAT': 'TAY', 'TGA': 'TGR', 'TGC': 'TGY', 'TGG': 'TGR', 'TGT': 'TGY', 'TTA': 'YTN',
                  'TTG': 'YTN'}

    # Euplotid Nuclear Code codons
    dgndict[8] = {'AAA': 'AAR', 'AAC': 'AAY', 'AAG': 'AAR', 'AAT': 'AAY', 'AGA': 'MGN', 'AGC': 'AGY', 'AGG': 'MGN',
                  'AGT': 'AGY', 'ATA': 'ATH', 'ATC': 'ATH', 'ATG': 'ATG', 'ATT': 'ATH', 'CAA': 'CAR', 'CAG': 'CAR',
                  'CGA': 'MGN', 'CGC': 'MGN', 'CGG': 'MGN', 'CGT': 'MGN', 'CTA': 'YTN', 'CTC': 'YTN', 'CTG': 'YTN',
                  'CTT': 'YTN', 'GGA': 'GGN', 'GGC': 'GGN', 'GGG': 'GGN', 'GGT': 'GGN', 'TAA': 'TAA', 'TAC': 'TAY',
                  'TAG': 'TAG', 'TAT': 'TAY', 'TGA': 'TGH', 'TGC': 'TGH', 'TGG': 'TGG', 'TGT': 'TGH', 'TTA': 'YTN',
                  'TTG': 'YTN'}

    # Bacterial / Archaeal / Plant Plastid Code codons
    dgndict[9] = {'AAA': 'AAR', 'AAC': 'AAY', 'AAG': 'AAR', 'AAT': 'AAY', 'AGA': 'MGN', 'AGC': 'AGY', 'AGG': 'MGN',
                  'AGT': 'AGY', 'ATA': 'ATH', 'ATC': 'ATH', 'ATG': 'ATG', 'ATT': 'ATH', 'CAA': 'CAR', 'CAG': 'CAR',
                  'CGA': 'MGN', 'CGC': 'MGN', 'CGG': 'MGN', 'CGT': 'MGN', 'CTA': 'YTN', 'CTC': 'YTN', 'CTG': 'YTN',
                  'CTT': 'YTN', 'GGA': 'GGN', 'GGC': 'GGN', 'GGG': 'GGN', 'GGT': 'GGN', 'TAA': 'TAA', 'TAC': 'TAY',
                  'TAG': 'TAG', 'TAT': 'TAY', 'TGA': 'TGA', 'TGC': 'TGY', 'TGG': 'TGG', 'TGT': 'TGY', 'TTA': 'YTN',
                  'TTG': 'YTN'}

    # Alternative Yeast Nuclear Code codons
    dgndict[10] = {'AAA': 'AAR', 'AAC': 'AAY', 'AAG': 'AAR', 'AAT': 'AAY', 'AGA': 'MGN', 'AGC': 'AGY', 'AGG': 'MGN',
                   'AGT': 'AGY', 'ATA': 'ATH', 'ATC': 'ATH', 'ATG': 'ATG', 'ATT': 'ATH', 'CAA': 'CAR', 'CAG': 'CAR',
                   'CGA': 'MGN', 'CGC': 'MGN', 'CGG': 'MGN', 'CGT': 'MGN', 'CTA': 'YTN', 'CTC': 'YTN', 'CTG': 'CTG',
                   'CTT': 'YTN', 'GGA': 'GGN', 'GGC': 'GGN', 'GGG': 'GGN', 'GGT': 'GGN', 'TAA': 'TAA', 'TAC': 'TAY',
                   'TAG': 'TAG', 'TAT': 'TAY', 'TGA': 'TGA', 'TGC': 'TGY', 'TGG': 'TGG', 'TGT': 'TGY', 'TTA': 'YTN',
                   'TTG': 'YTN'}

    # Ascidian Mitochondrial Code codons
    dgndict[11] = {'AAA': 'AAR', 'AAC': 'AAY', 'AAG': 'AAR', 'AAT': 'AAY', 'AGA': 'RGN', 'AGC': 'AGY', 'AGG': 'RGN',
                   'AGT': 'AGY', 'ATA': 'ATR', 'ATC': 'ATY', 'ATG': 'ATR', 'ATT': 'ATY', 'CAA': 'CAR', 'CAG': 'CAR',
                   'CGA': 'CGN', 'CGC': 'CGN', 'CGG': 'CGN', 'CGT': 'CGN', 'CTA': 'YTN', 'CTC': 'YTN', 'CTG': 'YTN',
                   'CTT': 'YTN', 'GGA': 'RGN', 'GGC': 'RGN', 'GGG': 'RGN', 'GGT': 'RGN', 'TAA': 'TAA', 'TAC': 'TAY',
                   'TAG': 'TAG', 'TAT': 'TAY', 'TGA': 'TGR', 'TGC': 'TGY', 'TGG': 'TGR', 'TGT': 'TGY', 'TTA': 'YTN',
                   'TTG': 'YTN'}

    # Alternative Flatworm Mitochondrial Code codons
    dgndict[12] = {'AAA': 'AAH', 'AAC': 'AAH', 'AAG': 'AAG', 'AAT': 'AAH', 'AGA': 'AGN', 'AGC': 'AGN', 'AGG': 'AGN',
                   'AGT': 'AGN', 'ATA': 'ATH', 'ATC': 'ATH', 'ATG': 'ATG', 'ATT': 'ATH', 'CAA': 'CAR', 'CAG': 'CAR',
                   'CGA': 'CGN', 'CGC': 'CGN', 'CGG': 'CGN', 'CGT': 'CGN', 'CTA': 'YTN', 'CTC': 'YTN', 'CTG': 'YTN',
                   'CTT': 'YTN', 'GGA': 'GGN', 'GGC': 'GGN', 'GGG': 'GGN', 'GGT': 'GGN', 'TAA': 'TAH', 'TAC': 'TAH',
                   'TAG': 'TAG', 'TAT': 'TAH', 'TGA': 'TGR', 'TGC': 'TGY', 'TGG': 'TGR', 'TGT': 'TGY', 'TTA': 'YTN',
                   'TTG': 'YTN'}

    # Add variable codons to working dictionary
    try:
        base_dict.update(dgndict[table])
    except KeyError:
        raise KeyError('Could not locate codon dictionary. Supported codon tables are numbered 1 through 12')

    if seqbuddy.alpha == "protein":
        raise TypeError("Nucleic acid sequence required, not protein.")

    if seqbuddy.alpha == "RNA":
        rna2dna(seqbuddy)

    clean_seq(seqbuddy)
    uppercase(seqbuddy)
    for _rec in seqbuddy.records:
        seq_length = len(str(_rec.seq))
        i = 0
        degen_string = ""
        while i < seq_length:
            _codon = str(_rec.seq[i:i + 3])
            degen_string += base_dict[_codon] if _codon in base_dict else _codon
            i += 3
        _rec.seq = Seq(str(degen_string))
        _rec.annotations["molecule_type"] = "DNA"
    return seqbuddy


def delete_features(seqbuddy, pattern):
    """
    Deletes features with IDs matching a regex pattern
    :param seqbuddy: SeqBuddy object
    :param pattern: The regex pattern to search with
    :return: The modified SeqBuddy object
    """
    for rec in seqbuddy.records:
        retained_features = []
        for feat in rec.features:
            keep = True
            if re.search(pattern, feat.type):
                keep = False
            else:
                for key, quals in feat.qualifiers.items():
                    if re.search(pattern, key):
                        keep = False
                        break
                    breakout = False
                    for qual in quals:
                        if re.search(pattern, qual):
                            keep = False
                            breakout = True
                            break
                    if breakout:
                        break
            if keep:
                retained_features.append(feat)
        rec.features = retained_features
    return seqbuddy


def delete_large(seqbuddy, max_value):
    """
    Deletes records with sequences larger than a certain size
    :param seqbuddy: SeqBuddy object
    :param max_value: The maximum threshold for sequence length
    :return: The modified SeqBuddy object
    """
    retained_records = []
    for rec in seqbuddy.records:
        if len(str(rec.seq)) <= max_value:
            retained_records.append(rec)
    seqbuddy.records = retained_records
    return seqbuddy


def delete_metadata(seqbuddy):
    """
    Removes all metadata from records
    :param seqbuddy: SeqBuddy object
    :return: The modified SeqBuddy object
    """
    if seqbuddy.alpha == "protein":
        temp_seqbuddy = SeqBuddy(">seq1\nAPPTRMWQQQRQQQPC", in_format="fasta")
    elif seqbuddy.alpha == "RNA":
        temp_seqbuddy = SeqBuddy(">seq1\nAUGCGCUAGCAUGUCA", in_format="fasta")
    else:
        temp_seqbuddy = SeqBuddy(">seq1\nATGCGCTAGCATGTCA", in_format="fasta")
    for rec in seqbuddy.records:
        rec.name = ''
        rec.description = ''
        rec.annotations = temp_seqbuddy.records[0].annotations
        rec.features = temp_seqbuddy.records[0].features
    return seqbuddy


def delete_records(seqbuddy, patterns, description=False):
    """
    Deletes records with IDs matching a regex pattern
    :param seqbuddy: SeqBuddy object
    :param patterns: A single regex pattern, or list of patterns, to search with
    :type patterns: list str
    :param description: Allow search in description string
    :return: The modified SeqBuddy object
    """
    if type(patterns) == str:
        patterns = [patterns]
    if type(patterns) != list:
        raise ValueError("'patterns' must be a list or a string.")

    for indx, pattern in enumerate(patterns):
        patterns[indx] = ".*" if pattern == "*" else pattern

    patterns = "|".join(patterns)

    retained_records = [None for _ in range(len(seqbuddy))]
    deleted = [rec.id for rec in pull_recs(make_copy(seqbuddy), patterns, description=description).records]
    counter = 0
    for rec in seqbuddy.records:
        if rec.id in deleted:
            continue
        else:
            retained_records[counter] = rec
            counter += 1
    seqbuddy.records = retained_records[:counter]
    return seqbuddy


def delete_recs_with_feature(seqbuddy, regex):
    """
    Remove sequences with feature names/IDs matching a search pattern
    :param seqbuddy: SeqBuddy object
    :param regex: List of regex expressions or single regex
    :type regex: str list
    :return: The modified SeqBuddy object
    """
    if type(regex) == str:
        regex = [regex]
    for indx, pattern in enumerate(regex):
        regex[indx] = ".*" if pattern == "*" else pattern

    regex = "|".join(regex)
    keep_records = []
    for rec in seqbuddy.records:
        if re.search(regex, rec.description):
            continue
        else:
            retain = True

        for feat in rec.features:
            if re.search(regex, feat.type) or re.search(regex, feat.id):
                retain = False
                break
            breakout = False
            for key, quals in feat.qualifiers.items():
                if breakout:
                    break
                if re.search(regex, key):
                    retain = False
                    breakout = True
                    break
                for qual in quals:
                    if re.search(regex, qual):
                        retain = False
                        breakout = True
                        break
            if breakout:
                break

        if retain:
            keep_records.append(rec)
    seqbuddy.records = keep_records
    return seqbuddy


def delete_repeats(seqbuddy, scope='all'):  # scope in ['all', 'ids', 'seqs']
    """
    Deletes records with repeated IDs/seqs
    :param seqbuddy: SeqBuddy object
    :param scope: Specifies if deleting repeat seqs, ids, or all
    :return: The modified SeqBuddy object
    """
    # First, remove duplicate IDs
    if scope in ['all', 'ids']:
        find_repeats(seqbuddy)
        if len(seqbuddy.repeat_ids) > 0:
            stored_recs = [None for _ in range(len(seqbuddy.repeat_ids))]
            for indx, rep_id in enumerate(seqbuddy.repeat_ids):
                stored_recs[indx] = pull_recs(make_copy(seqbuddy), ["^%s$" % rep_id]).records[0]
            regex = "^%s$" % "$|^".join(seqbuddy.repeat_ids)
            delete_records(seqbuddy, regex)
            seqbuddy.records += stored_recs

    # Then remove duplicate sequences
    if scope in ['all', 'seqs']:
        find_repeats(seqbuddy)
        if len(seqbuddy.repeat_seqs) > 0:
            rep_seq_ids = []
            for seq in seqbuddy.repeat_seqs:
                rep_seq_ids.append([])
                for rep_seq_id in seqbuddy.repeat_seqs[seq]:
                    rep_seq_ids[-1].append(rep_seq_id)

            repeat_regex = ""

            for repeat_seqs in rep_seq_ids:
                for rep_seq in repeat_seqs[1:]:
                    rep_seq = re.sub(r"([|.*?^\[\]()])", r"\\\1", rep_seq)
                    repeat_regex += "^%s$|" % rep_seq

            repeat_regex = repeat_regex[:-1]
            delete_records(seqbuddy, repeat_regex)

    seqbuddy.repeat_seqs = OrderedDict()
    seqbuddy.repeat_ids = OrderedDict()
    seqbuddy.unique_seqs = OrderedDict([(x.id, x) for x in seqbuddy.records])
    return seqbuddy


def delete_small(seqbuddy, min_value):
    """
    Deletes records with sequence smaller than a certain size
    :param seqbuddy: SeqBuddy object
    :param min_value: The minimum threshold for sequence length
    :return: The modified SeqBuddy object
    """
    retained_records = []
    for rec in seqbuddy.records:
        if len(str(rec.seq)) >= min_value:
            retained_records.append(rec)
    seqbuddy.records = retained_records
    return seqbuddy


def delete_taxa(seqbuddy, taxa):
    """
    Pull out records that are annotated with a particular taxon
    :param seqbuddy:
    :param taxa: List of exact matches for taxonomic designation (will sift through entire taxon tree)
    :return:
    """
    if type(taxa) == str:
        taxa = [taxa]

    taxa = [t.lower() for t in taxa]
    keep_list = []
    for rec in seqbuddy.records:
        delete = False
        if 'taxonomy' in rec.annotations:
            taxonomy = [x.lower() for x in rec.annotations['taxonomy']]
            breakout = False
            for taxon in taxa:
                if taxon in taxonomy:
                    delete = True
                    breakout = True
                    break
            if breakout:
                continue

        if 'organism' in rec.annotations:
            for taxon in taxa:
                organism = rec.annotations['organism'].lower().split()
                if taxon in organism:
                    delete = True
                    break
        if not delete:
            keep_list.append(rec)

    seqbuddy.records = keep_list
    return seqbuddy


def dna2rna(seqbuddy):
    """
    Transcribes DNA into RNA
    :param seqbuddy: SeqBuddy object
    :return: Modified SeqBuddy object
    """
    if seqbuddy.alpha != "DNA":
        raise TypeError("DNA sequence required, not %s." % seqbuddy.alpha)
    for rec in seqbuddy.records:
        rec.seq = Seq(str(rec.seq.transcribe()))
        rec.annotations["molecule_type"] = "RNA"
    seqbuddy.alpha = "RNA"
    return seqbuddy


def extract_feature_sequences(seqbuddy, patterns):
    """
    Pull out specific features from annotated sequences
    :param seqbuddy: SeqBuddy object
    :type seqbuddy: SeqBuddy
    :param patterns: The feature(s) to be extracted
    :type patterns: list str
    :return: Modified SeqBuddy object
    :rtype: SeqBuddy
    """
    def check_pattern(_feature, _pattern):
        if re.search(_pattern, _feature.type):
            return True
        else:
            for qual_type, quals in feature.qualifiers.items():
                for qual in quals:
                    if re.search(_pattern, qual):
                        return True
        return False

    if type(patterns) == str:
        patterns = [patterns]

    range_patterns = []
    single_patterns = []
    for pattern in patterns:
        if ":" in pattern:
            range_patterns.append(pattern.split(":"))
        else:
            single_patterns.append(pattern)

    new_recs = []
    for rec in seqbuddy.records:
        keep_ranges = []
        for feature in rec.features:
            check_single_patterns = [check_pattern(feature, pat) for pat in single_patterns]
            if True in check_single_patterns:
                if type(feature.location) == CompoundLocation:
                    keep_ranges += [[int(x.start), int(x.end)] for x in feature.location.parts]
                else:
                    keep_ranges.append([int(feature.location.start), int(feature.location.end)])
        for rang_pat in range_patterns:
            start, end = len(rec.seq), 0
            pat1, pat2 = False, False
            for feature in rec.features:
                if check_pattern(feature, rang_pat[0]):
                    start = int(feature.location.start) if int(feature.location.start) < start else start
                    end = int(feature.location.end) if int(feature.location.end) > end else end
                    pat1 = True
                if check_pattern(feature, rang_pat[1]):
                    start = int(feature.location.start) if int(feature.location.start) < start else start
                    end = int(feature.location.end) if int(feature.location.end) > end else end
                    pat2 = True
            if pat1 and pat2:
                keep_ranges.append([start, end])

        if not keep_ranges:
            rec.seq = Seq("")
            rec.features = []
            new_recs.append(rec)
        else:
            keep_ranges = sorted(keep_ranges, key=lambda x: x[0])
            final_positions = ""
            active_range = keep_ranges[0]
            for _range in keep_ranges[1:]:
                if active_range[1] >= _range[0]:
                    active_range[1] = max(active_range[1], _range[1])
                else:
                    final_positions += "%s:%s," % (active_range[0] + 1, active_range[1])
                    active_range = [_range[0], _range[1]]

            final_positions += "%s:%s" % (active_range[0] + 1, active_range[1])
            single_rec = SeqBuddy([rec])
            single_rec = extract_regions(single_rec, final_positions)
            new_recs.append(single_rec.records[0])

    seqbuddy.records = new_recs
    return seqbuddy


def extract_regions(seqbuddy, positions):
    """
    Fine grained control of what residues to pull out of the sequences
    :param seqbuddy: SeqBuddy object
    :param positions: Position code describing which residues to pull (str)

    Position Code:  - Always a string
                    - Comma-separated
                    - Three types of extraction:
                        - Singlets: "2,5,9,-5"
                        - Ranges: "40:75,89:100,432:-45"
                        - mth of nth: "1/5,3/5"
    """
    positions = re.sub(r"\s|[,/-]$|^[,/]", "", positions).split(",")

    def process_single(num, max_len):
        if num == 0:
            num = 1
        elif num < 0:
            if abs(num) > max_len:
                num = 1
            else:
                num += max_len + 1
        elif num > max_len:
            num = max_len
        return num

    def create_residue_list(_rec, _positions):
        rec_len = len(_rec.seq)
        singlets = []
        for _position in _positions:
            # Singlets
            try:
                single = process_single(int(_position), rec_len)
                singlets.append(single - 1)
                continue
            except ValueError as e:
                if "invalid literal for int() with base 10" in str(e):
                    pass
                else:
                    raise e
            try:
                # mth of nth
                if "/" in _position:
                    start, end = _position.split("/")
                    # NOTE: if end is negative, it is converted to the positive from the back of the sequence (this may
                    # be confusing behaviour).
                    end = process_single(int(end), rec_len)
                    if ":" in start:
                        range_start, range_end = start.split(":")
                        range_end = process_single(-1, end) if not range_end else process_single(int(range_end), end)
                        range_start = 1 if not range_start else process_single(int(range_start), end)

                        if range_start > range_end:
                            raise ValueError

                        for i in range(range_start, range_end + 1):
                            singlets += range(i - 1, rec_len, end)

                    else:
                        start = process_single(int(start), end)
                        singlets += range(start - 1, rec_len, end)

                # Ranges
                elif ":" in _position:
                    start, end = _position.split(":")
                    start = 1 if not start else process_single(int(start), rec_len)
                    end = process_single(-1, rec_len) if not end else process_single(int(end), rec_len)
                    start, end = sorted([start, end])
                    singlets += range(start - 1, end)

                # Fail...
                else:
                    raise ValueError()

            except ValueError:
                raise ValueError("Unable to decode the positions string '%s'." % _position)

        singlets = list(set(singlets))
        singlets = sorted(singlets)
        return singlets

    new_records = []
    for rec in seqbuddy.records:
        new_rec_positions = create_residue_list(rec, positions)
        letter_annotations = {}
        for anno_type in rec.letter_annotations:
            letter_annotations[anno_type] = [None for _ in range(len(new_rec_positions))]
            for indx, anno_pos in enumerate(new_rec_positions):
                letter_annotations[anno_type][indx] = rec.letter_annotations[anno_type][anno_pos]
        new_seq = []
        if rec.features:  # This is super slow for large records...
            remapper = FeatureReMapper(rec)
            for indx, residue in enumerate(str(rec.seq)):
                if indx in new_rec_positions:
                    remapper.extend(True)
                    new_seq.append(residue)
                else:
                    remapper.extend(False)
            new_seq = ''.join(new_seq)
            new_seq = Seq(new_seq)
            new_seq = SeqRecord(new_seq, id=rec.id, name=rec.name, description=rec.description, dbxrefs=rec.dbxrefs,
                                annotations=rec.annotations, letter_annotations=letter_annotations)
            new_seq = remapper.remap_features(new_seq)
        else:
            seq = str(rec.seq)
            for indx in new_rec_positions:
                new_seq.append(seq[indx])
            new_seq = ''.join(new_seq)
            new_seq = Seq(new_seq)
            new_seq = SeqRecord(new_seq, id=rec.id, name=rec.name, description=rec.description, dbxrefs=rec.dbxrefs,
                                annotations=rec.annotations, letter_annotations=letter_annotations)

        new_records.append(new_seq)

    seqbuddy.records = new_records
    return seqbuddy


def find_cpg(seqbuddy):
    """
    Predicts locations of CpG islands in DNA sequences
    :param seqbuddy: SeqBuddy object
    :return: Modified SeqBuddy object (buddy_data["cpgs"] appended to all records)
    """
    seqbuddy = clean_seq(seqbuddy)
    if seqbuddy.alpha != "DNA":
        raise TypeError("DNA sequence required, not protein or RNA.")

    records = []

    def cpg_calc(in_seq):  # Returns observed/expected value of a sequence
        in_seq = in_seq.upper()
        observed_cpg = len(re.findall("CG", in_seq)) * len(in_seq)
        expected = (len(re.findall("[CG]", in_seq)) / 2) ** 2
        expected = 1 if not expected else expected  # Prevent DivByZero
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

    for rec in seqbuddy.records:
        seq = rec.seq

        oe_vals_list = [0 for _ in range(len(seq))]
        cg_percent_list = [0 for _ in range(len(seq))]
        window_size = len(seq) if len(seq) < 200 else 200

        for indx, window in enumerate([seq[i:i + window_size] for i in range(len(seq) - window_size + 1)]):
            cpg = cpg_calc(str(window))
            for i in range(window_size):
                oe_vals_list[indx + i] += cpg

            cg_perc = cg_percent(str(window))
            for i in range(window_size):
                cg_percent_list[indx + i] += cg_perc

        for indx in range(len(oe_vals_list)):
            if indx + 1 <= window_size:
                oe_vals_list[indx] /= (indx + 1)
                cg_percent_list[indx] /= (indx + 1)

            elif (len(oe_vals_list) - window_size) - (indx + 1) < 0:
                oe_vals_list[indx] /= (len(oe_vals_list) - indx)
                cg_percent_list[indx] /= (len(cg_percent_list) - indx)

            else:
                oe_vals_list[indx] /= window_size
                cg_percent_list[indx] /= window_size

        indices = find_islands(cg_percent_list, oe_vals_list)
        cpg_features = [SeqFeature(location=FeatureLocation(start, end), type="CpG_island",
                                   qualifiers={'created_by': 'SeqBuddy'}) for (start, end) in indices]
        for feature in rec.features:
            cpg_features.append(feature)
        rec = SeqRecord(map_cpg(seq, indices), id=rec.id, name=rec.name, description=rec.description,
                        dbxrefs=rec.dbxrefs, features=cpg_features, annotations=rec.annotations,
                        letter_annotations=rec.letter_annotations)

        records.append(rec)
        _add_buddy_data(rec, "cpgs", indices)
    seqbuddy.records = records
    return seqbuddy


def find_orfs(seqbuddy, include_feature=True, include_buddy_data=True, min_size=None, rev_comp=True):
    """
    Finds all the open reading frames in the sequences and their reverse complements.
    :param seqbuddy: SeqBuddy object
    :param include_feature: Add a new 'orf' feature to records
    :param include_buddy_data: Append information directly to records
    :param min_size: Smallest ORF that is accepted
    :param rev_comp: Also include ORFS in the reverse complement sequence
    :return: Annotated SeqBuddy object. The match indices are also stored in rec.buddy_data["find_orfs"].
    """
    if seqbuddy.alpha == "protein":
        raise TypeError("Nucleic acid sequence required, not protein.")

    if min_size and min_size < 6:
        raise ValueError("Open reading frames cannot be smaller than 6 residues.")

    pattern = "a[tu]g(...)*?([tu]aa|[tu]ag|[tu]ga)"

    clean_seq(seqbuddy)
    lowercase(seqbuddy)

    find_pattern(seqbuddy, pattern, ambig=True, include_feature=True, include_buddy_data=False)
    lowercase(seqbuddy)

    if rev_comp:
        reverse_complement(seqbuddy)
        find_pattern(seqbuddy, pattern, ambig=True, include_feature=True, include_buddy_data=False)
        lowercase(seqbuddy)

        reverse_complement(seqbuddy)

    for rec in seqbuddy.records:
        buddy_data = {'+': [], '-': []}
        feature_indicies = []
        orf_list = []
        min_size = min_size if min_size else 9

        for indx, feature in enumerate(rec.features):
            if 'regex' in feature.qualifiers.keys() and feature.qualifiers['regex'] == pattern:
                feature.type = "orf" if feature.strand == +1 else "comp_orf"
                strand = "+" if feature.strand == +1 else "-"
                feature.qualifiers.pop("regex")
                feature_indicies.append(indx)
                if len(feature) >= min_size:
                    buddy_data[strand].append((int(feature.location.start), int(feature.location.end)))
                    orf_list.append(feature)

        for indx in sorted(feature_indicies, reverse=True):
            del rec.features[indx]

        if include_feature:
            orf_list = [(int(feature.location.end) - int(feature.location.start), feature) for feature in orf_list]
            orf_list = sorted(orf_list, key=lambda x: x[0], reverse=True)
            orf_list = [feature[1] for feature in orf_list]
            orf_count = 1
            comp_orf_count = 1
            for orf in orf_list:
                if orf.type == "orf":
                    orf.type = "orf%s" % orf_count
                    orf_count += 1
                else:
                    orf.type = "comp_orf%s" % comp_orf_count
                    comp_orf_count += 1

            rec.features += orf_list

        if include_buddy_data:
            _add_buddy_data(rec, 'find_orfs')
            rec.buddy_data['find_orfs'] = buddy_data
    return seqbuddy


def find_pattern(seqbuddy, *patterns, ambig=False, include_feature=True, include_buddy_data=True):
    """
    Finds ﻿occurrences of a sequence pattern
    :param seqbuddy: SeqBuddy object
    :param patterns: regex patterns
    :param ambig: Convert any ambiguous letter codes in the search pattern into regex
    :param include_feature: Add a new 'match' feature to records
    :param include_buddy_data: Append information directly to records
    :return: Annotated SeqBuddy object. The match indices are also stored in rec.buddy_data["find_patterns"].
    """
    # search through sequences for regex matches. For example, to find micro-RNAs
    lowercase(seqbuddy)
    for pattern in patterns:
        pattern_backup = str(pattern)
        if ambig and seqbuddy.alpha == "protein":
            pattern = re.sub("[xX]", "[ARNDCQEGHILKMFPSTWYVX]", pattern)
            pattern = re.sub("[bB]", "[NDB]", pattern)
            pattern = re.sub("[zZ]", "[QEZ]", pattern)

        elif ambig and seqbuddy.alpha == "DNA":
            pattern = re.sub("[kK]", "[GT]", pattern)
            pattern = re.sub("[mM]", "[AC]", pattern)
            pattern = re.sub("[rR]", "[AG]", pattern)
            pattern = re.sub("[yY]", "[CT]", pattern)
            pattern = re.sub("[sS]", "[CG]", pattern)
            pattern = re.sub("[wW]", "[AT]", pattern)
            pattern = re.sub("[bB]", "[CGT]", pattern)
            pattern = re.sub("[vV]", "[CGA]", pattern)
            pattern = re.sub("[hH]", "[ACT]", pattern)
            pattern = re.sub("[dD]", "[AGT]", pattern)
            pattern = re.sub("[xnXN]", "[ATCG]", pattern)

        elif ambig and seqbuddy.alpha == "RNA":
            pattern = re.sub("[kK]", "[GU]", pattern)
            pattern = re.sub("[mM]", "[AC]", pattern)
            pattern = re.sub("[rR]", "[AG]", pattern)
            pattern = re.sub("[yY]", "[CU]", pattern)
            pattern = re.sub("[sS]", "[CG]", pattern)
            pattern = re.sub("[wW]", "[AU]", pattern)
            pattern = re.sub("[bB]", "[CGU]", pattern)
            pattern = re.sub("[vV]", "[CGA]", pattern)
            pattern = re.sub("[hH]", "[ACU]", pattern)
            pattern = re.sub("[dD]", "[AGU]", pattern)
            pattern = re.sub("[xnXN]", "[AUCG]", pattern)

        if ambig:
            safety_valve = br.SafetyValve()
            # Strip out any double square brackets
            while re.search(r"\[[^[\]]*?\[[^]]*]", pattern):
                safety_valve.step("Ambiguous %s regular expression '%s' failed compile." %
                                  (seqbuddy.alpha, pattern_backup))
                pattern = re.sub(r"(\[[^[\]]*?)\[([^]]*)]", r"\1\2", pattern, count=1)

        for rec in seqbuddy.records:
            if include_buddy_data:
                _add_buddy_data(rec, 'find_patterns')
            indices = []
            matches = re.finditer(pattern, str(rec.seq), flags=re.IGNORECASE)
            new_seq = ""
            last_match = 0
            for match in matches:
                indices.append([match.start(), match.end()])
                if include_feature:
                    rec.features.append(SeqFeature(location=FeatureLocation(start=match.start(), end=match.end()),
                                                   type='match', strand=+1,
                                                   qualifiers=OrderedDict([('regex', pattern_backup),
                                                                           ('added_by', 'SeqBuddy')])))
                if match.start() > 0:
                    new_seq += str(rec.seq[last_match:match.start()])
                new_seq += str(rec.seq[match.start():match.end()]).upper()
                last_match = match.end()
            new_seq += str(rec.seq[last_match:])
            rec.seq = Seq(new_seq)

            if include_buddy_data:
                if not rec.buddy_data['find_patterns']:
                    rec.buddy_data['find_patterns'] = OrderedDict({pattern_backup: indices})
                else:
                    rec.buddy_data['find_patterns'][pattern_backup] = indices
    return seqbuddy


def find_repeats(seqbuddy):
    """
    Finds sequences with identical IDs or sequences
    :param seqbuddy: SeqBuddy object
    :return: modified seqbuddy object with three new attributes --> unique_seqs, repeat_ids, and repeat_seqs
    """
    unique_seqs = OrderedDict()
    repeat_ids = OrderedDict()
    repeat_seqs = OrderedDict()

    # First find replicate IDs
    # MD5 hash all sequences as we go for memory efficiency when looking for replicate sequences (below)
    # Need to work from a copy though, so sequences aren't overwritten
    seqbuddy_copy = make_copy(seqbuddy)
    for rec in seqbuddy_copy.records:
        seq = str(rec.seq).encode("utf-8")
        seq = md5(seq).hexdigest()
        rec.seq = Seq(seq)
        if rec.id in repeat_ids:
            repeat_ids[rec.id].append(rec)
        elif rec.id in unique_seqs:
            repeat_ids[rec.id] = [rec]
            repeat_ids[rec.id].append(unique_seqs[rec.id])
            del (unique_seqs[rec.id])
        else:
            unique_seqs[rec.id] = rec

    # Then look for replicate sequences
    flip_uniqe = {}
    del_keys = []
    for key, value in unique_seqs.items():  # find and remove duplicates in/from the unique list
        value = str(value.seq)
        if value not in flip_uniqe:
            flip_uniqe[value] = [key]
        else:
            if value not in repeat_seqs:
                repeat_seqs[value] = [key]
                repeat_seqs[value] += flip_uniqe[value]
                if flip_uniqe[value][0] in unique_seqs:
                    del_keys.append(flip_uniqe[value][0])
            else:
                repeat_seqs[value].append(key)
            del_keys.append(unique_seqs[key].id)

    for key in del_keys:
        if key in unique_seqs:
            del (unique_seqs[key])

    for key, value in repeat_ids.items():  # find duplicates in the repeat ID list
        for rep_seq in value:
            rep_seq = str(rep_seq.seq)
            if rep_seq not in flip_uniqe:
                flip_uniqe[rep_seq] = [key]
            else:
                if rep_seq not in repeat_seqs:
                    repeat_seqs[rep_seq] = [key]
                    repeat_seqs[rep_seq] += flip_uniqe[rep_seq]

                else:
                    repeat_seqs[rep_seq].append(key)

    seqbuddy.unique_seqs = unique_seqs
    seqbuddy.repeat_ids = repeat_ids
    seqbuddy.repeat_seqs = repeat_seqs
    return seqbuddy


def find_restriction_sites(seqbuddy, enzyme_group=(), min_cuts=1, max_cuts=None, topology=None, quiet=False):
    """
    Finds the restriction sites in the sequences in the SeqBuddy object
    :param seqbuddy: SeqBuddy object
    :param enzyme_group: "commercial", "all", or a list of specific enzyme names
    :param min_cuts: The minimum cut threshold
    :param max_cuts: The maximum cut threshold
    :param topology: linear or circular sequence
    :param quiet: Suppress stderr
    :return: annotated SeqBuddy object, and a dictionary of restriction sites added as the `restriction_sites` attribute
    """

    if seqbuddy.alpha == "protein":
        raise TypeError("Unable to identify restriction sites in protein sequences.")

    convert_rna = False
    if seqbuddy.alpha == "RNA":
        convert_rna = True
        rna2dna(seqbuddy)

    if max_cuts and min_cuts > max_cuts:
        raise ValueError("min_cuts parameter has been set higher than max_cuts.")
    max_cuts = 1000000000 if not max_cuts else max_cuts

    enzyme_group = list(enzyme_group) if enzyme_group else ["commercial"]

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
    for enzyme in enzyme_group:
        if enzyme == "commercial":
            for res in CommOnly:
                if str(res) not in blacklist:
                    batch.add(res)

        elif enzyme == "all":
            for res in sorted([x for x in AllEnzymes], key=lambda x: str(x)):
                if str(res) not in blacklist:
                    batch.add(res)

        else:
            try:
                batch.add(enzyme)
            except ValueError:
                br._stderr("Warning: %s not a known enzyme\n" % enzyme, quiet=quiet)

    sites = []
    no_cutters_found = False
    double_cutters_found = False
    for rec in seqbuddy.records:
        linear = True
        if not topology and 'topology' in rec.annotations:
            linear = False if rec.annotations['topology'] == 'circular' else True
        elif topology == "circular":
            linear = False
        features = ["%s%s%s" % (feat.type, feat.location.start, feat.location.end) for feat in rec.features]
        rec.res_sites = {}
        analysis = Analysis(batch, rec.seq, linear=linear)
        result = analysis.with_sites()
        for key, value in result.items():
            if key.cut_twice():
                if not double_cutters_found:
                    br._stderr("Warning: Double-cutters not supported.\n", quiet=quiet)
                    double_cutters_found = True
                pass
            elif min_cuts <= len(value) <= max_cuts:
                try:
                    for zyme in value:
                        cut_start = zyme + key.fst3 - 1
                        cut_end = zyme + key.fst5 + abs(key.ovhg) - 1
                        if "%s%s%s" % (str(key), cut_start, cut_end) not in features:
                            rec.features.append(SeqFeature(FeatureLocation(start=cut_start, end=cut_end),
                                                           type=str(key)))
                except TypeError:
                    if not no_cutters_found:
                        br._stderr("Warning: No-cutters not supported.\n", quiet=quiet)
                        no_cutters_found = True
                    pass
                rec.res_sites[key] = value
        rec.res_sites = OrderedDict(sorted(rec.res_sites.items(), key=lambda x: x[0]))
        sites.append((rec.id, rec.res_sites))

    order_features_alphabetically(seqbuddy)
    seqbuddy.restriction_sites = sites
    if convert_rna:
        dna2rna(seqbuddy)
    return seqbuddy


def hash_ids(seqbuddy, hash_length=10, r_seed=None):
    """
    Replaces the sequence IDs with random hashes
    :param seqbuddy: SeqBuddy object
    :param hash_length: Specifies the length of the new hashed IDs
    :param r_seed: Set the random generator seed value
    :return: The modified SeqBuddy object, with a new attribute `hash_map` added
    """
    hash_list = []
    seq_ids = []

    try:
        hash_length = int(hash_length)
    except ValueError:
        raise TypeError("Hash length argument must be an integer, not %s" % type(hash_length))

    # If a hash_map already exists and fits all the specs, re-apply it.
    if seqbuddy.hash_map and len(seqbuddy.hash_map) == len(seqbuddy):
        # work from a copy, just in case we find an id that doesn't match and we need to start from scratch
        seqbuddy_copy = make_copy(seqbuddy)
        re_apply_hash_map = True
        for rec, _hash_id in zip(seqbuddy_copy.records, list(seqbuddy_copy.hash_map.items())):
            if rec.id != _hash_id[1]:
                re_apply_hash_map = False
                break
            rec.id = _hash_id[0]
            rec.name = _hash_id[0]

        if re_apply_hash_map:
            seqbuddy.records = seqbuddy_copy.records
            return seqbuddy

    if hash_length < 1:
        raise ValueError("Hash length must be greater than 0")

    if 32 ** hash_length <= len(seqbuddy) * 2:
        raise ValueError("Insufficient number of hashes available to cover all sequences. "
                         "Hash length must be increased.")

    rand_gen = Random() if not r_seed else Random(r_seed)
    for i in range(len(seqbuddy)):
        seq_ids.append(seqbuddy.records[i].id)
        while True:
            chars = string.ascii_letters + string.digits
            new_hash = "".join([rand_gen.choice(chars) for _ in range(hash_length)])
            if new_hash in hash_list:
                continue
            else:
                hash_list.append(new_hash)
                break
        if re.match(seqbuddy.records[i].id, seqbuddy.records[i].description):
            seqbuddy.records[i].description = seqbuddy.records[i].description[len(seqbuddy.records[i].id) + 1:]

        seqbuddy.records[i].id = new_hash
        seqbuddy.records[i].name = new_hash

    hash_map = OrderedDict()
    for indx, _hash in enumerate(hash_list):
        hash_map[hash_list[indx]] = seq_ids[indx]
    seqbuddy.hash_map = hash_map
    return seqbuddy


def head(seqbuddy, count=1):
    """
    Pull the first record(s) in a SeqBuddy object
    :param seqbuddy: SeqBuddy object
    :param count: Number of records to return
    :return: The modified SeqBuddy object
    """
    seqbuddy.records = seqbuddy.records[:count]
    return seqbuddy


def insert_sequence(seqbuddy, sequence, location=0, regexes=None):
    """
    Add a specific sequence at a defined location in all records. E.g., adding a barcode (zero-indexed)
    :param seqbuddy: SeqBuddy object
    :param sequence: The sequence to be inserted
    :param location: The location to insert the sequences at
    :param regexes: Patterns to search for
    :return: The modified SeqBuddy object
    """
    # ToDo: Features... Move and add.
    if regexes:
        recs_to_update = pull_recs(make_copy(seqbuddy), regexes).to_dict()
    else:
        recs_to_update = seqbuddy.to_dict()

    for seq_id, rec in recs_to_update.items():
        if location >= 0:
            new_seq = str(rec.seq)[:location] + sequence + str(rec.seq)[location:]
        elif location == -1:
            new_seq = str(rec.seq) + sequence
        else:
            new_seq = str(rec.seq)[:location + 1] + sequence + str(rec.seq)[location + 1:]

        rec.seq = Seq(new_seq)

    for indx, rec in enumerate(seqbuddy.records):
        if rec.id in recs_to_update:
            seqbuddy.records[indx] = recs_to_update[rec.id]

    return seqbuddy


def in_silico_digest(seqbuddy, enzyme_group=(), quiet=False, topology=None):
    """
    Find restriction sites and break up sequences accordingly
    :param seqbuddy: SeqBuddy object
    :param enzyme_group: list of specific enzyme names
    :param quiet: Suppress stderr
    :param topology: circular or linear sequence
    :return: New seqbuddy object with sequences fragmented
    """
    if seqbuddy.alpha == "protein":
        raise TypeError("Unable to identify restriction sites in protein sequences.")

    seqbuddy_rs_lin = find_restriction_sites(make_copy(seqbuddy), enzyme_group, topology="linear", quiet=quiet)
    seqbuddy_rs_circ = find_restriction_sites(make_copy(seqbuddy), enzyme_group, topology="circular", quiet=quiet)
    new_records = []
    for indx, rec in enumerate(seqbuddy.records):
        sub_seqbuddy = SeqBuddy([rec])
        res_sites_lin = [cut_sites for enzym, cut_sites in seqbuddy_rs_lin.restriction_sites[indx][1].items()]
        res_sites_circ = [cut_sites for enzym, cut_sites in seqbuddy_rs_circ.restriction_sites[indx][1].items()]
        res_sites_lin = sorted([cut_site for sublist in res_sites_lin for cut_site in sublist])
        res_sites_circ = sorted([cut_site for sublist in res_sites_circ for cut_site in sublist])

        if topology == "circular":
            cut_type = "circular"
        elif topology == "linear":
            cut_type = "linear"
        elif not topology:
            if "topology" in rec.annotations:
                cut_type = rec.annotations["topology"]
            else:
                cut_type = "linear"
        else:
            raise ValueError("Invalid topology. Accepted values are None, 'circular' and 'linear' ")

        seq_pointer = 0
        new_fragments = []
        if cut_type == "circular":
            if not res_sites_circ:
                if "topology" in rec.annotations and rec.annotations["topology"] == "linear":
                    rec.annotations["topology"] = "circular"
                new_records.append(sub_seqbuddy.records[0])

            elif res_sites_circ[0] != 1:
                for cut in res_sites_circ:
                    sb_copy = extract_regions(make_copy(sub_seqbuddy), "%s:%s" % (seq_pointer, cut - 1))
                    fragment = amend_metadata(sb_copy, "topology", "linear", "")
                    new_fragments.append(fragment.records[0])
                    seq_pointer = cut
                sb_copy = extract_regions(make_copy(sub_seqbuddy), "%s:%s" % (seq_pointer, ""))
                new_fragments[0] = amend_metadata(sb_copy, "topology", "linear", "").records[0] + new_fragments[0]
                new_records.extend(new_fragments)

            elif res_sites_circ[0] == 1:
                for cut in res_sites_lin:
                    sb_copy = extract_regions(make_copy(sub_seqbuddy), "%s:%s" % (seq_pointer, cut - 1))
                    fragment = amend_metadata(sb_copy, "topology", "linear", "")
                    new_records.append(fragment.records[0])
                    seq_pointer = cut
                    sb_copy = extract_regions(make_copy(sub_seqbuddy), "%s:%s" % (seq_pointer, ""))
                    new_records.append(amend_metadata(sb_copy, "topology", "linear", "").records[0])

        elif cut_type == "linear":
            for cut in res_sites_lin:
                sb_copy = extract_regions(make_copy(sub_seqbuddy), "%s:%s" % (seq_pointer, cut - 1))
                fragment = amend_metadata(sb_copy, "topology", "linear", "")
                new_records.append(fragment.records[0])
                seq_pointer = cut
            final_fragment = extract_regions(make_copy(sub_seqbuddy), "%s:%s" % (seq_pointer, ""))
            if res_sites_lin:
                final_fragment = amend_metadata(final_fragment, "topology", "linear", "")
            new_records.append(final_fragment.records[0])

    seqbuddy.records = new_records
    return seqbuddy


def isoelectric_point(seqbuddy):
    """
    Calculate the isoelectric points
    :param seqbuddy: SeqBuddy object
    :return: SeqBuddy object with isoelectric point appended to each record as the last feature in the feature list
    """
    if seqbuddy.alpha != "protein":
        raise TypeError("Protein sequence required, not nucleic acid.")
    isoelectric_points = OrderedDict()
    for rec in seqbuddy.records:
        iso_point = ProteinAnalysis(str(rec.seq))
        iso_point = round(iso_point.isoelectric_point(), 10)
        isoelectric_points[rec.id] = iso_point
        rec.features.append(SeqFeature(location=FeatureLocation(start=0, end=len(rec.seq)), type='pI',
                                       qualifiers={'value': iso_point}))
    return seqbuddy


def keep_taxa(seqbuddy, taxa, match_all=False):
    """
    Pull out records that are annotated with a particular taxon
    :param seqbuddy:
    :param taxa: List of exact matches for taxonomic designation (will sift through entire taxon tree)
    :param match_all: When True, only accept record if all taxanomic levels are met
    return:
    """
    if type(taxa) == str:
        taxa = [taxa]

    taxa = [t.lower() for t in taxa]
    keep_list = []
    for rec in seqbuddy.records:
        if 'taxonomy' in rec.annotations:
            taxonomy = [x.lower() for x in rec.annotations['taxonomy']]
            if 'organism' in rec.annotations:
                organism = rec.annotations['organism'].split()[1:]
                taxonomy += organism

            append2keeplist = True if match_all else False
            for taxon in taxa:
                if match_all:
                    if taxon not in taxonomy:
                        append2keeplist = False
                        break
                else:
                    if taxon in taxonomy:
                        append2keeplist = True
                        break

            if append2keeplist:
                keep_list.append(rec)

    seqbuddy.records = keep_list
    return seqbuddy


def lowercase(seqbuddy):
    """
    Converts all sequence characters to lowercase.
    :param seqbuddy: SeqBuddy object
    :return: The modified SeqBuddy object
    """
    for rec in seqbuddy.records:
        rec.seq = Seq(str(rec.seq).lower())
    return seqbuddy


def make_groups(seqbuddy, split_patterns=(), num_chars=None, regex=None):
    """
    Splits a SeqBuddy object into new object based on an identifying suffix or regular expression
    :param seqbuddy: SeqBuddy object
    :param split_patterns: The regex pattern(s) to split with
    :param num_chars: Restrict the size of the identifier to a specific number of characters
    :param regex: Uses a regular expression to create group identifiers (downstream of split_patterns and num_chars).
    Duck typed --> str, list, tuple, None
    :return: A list of SeqBuddy objects
    """
    recs_by_identifier = OrderedDict()
    recs_by_identifier["Unknown"] = []
    split_prefix = False
    if split_patterns:
        split_patterns = [split_patterns] if type(split_patterns) == str else split_patterns
        for rec in seqbuddy.records:
            if len(re.split("|".join(split_patterns), rec.id)) > 1:
                split_prefix = True
                break

    for rec in seqbuddy.records:
        if split_prefix:
            split = re.split("|".join(split_patterns), rec.id)
            if len(split) == 1:
                recs_by_identifier["Unknown"].append(rec)
                continue
            else:
                split = split[0]
        else:
            split = rec.id
        split = split if not num_chars else split[:num_chars]
        if regex:
            regex = regex if type(regex) == str else "|".join(regex)
            split = re.search(regex, split)
            if not split:
                recs_by_identifier["Unknown"].append(rec)
                continue
            else:
                if re.search(r"\([^()]*\)", regex):
                    split = "".join(list(split.groups()))
                else:
                    split = split.group(0)
        if split != "":
            recs_by_identifier.setdefault(split, []).append(rec)
        else:
            recs_by_identifier["Unknown"].append(rec)

    if not recs_by_identifier["Unknown"]:
        del recs_by_identifier["Unknown"]

    new_seqbuddies = [(identifier, make_copy(seqbuddy)) for identifier in recs_by_identifier]
    for identifier, sb in new_seqbuddies:
        sb.records = recs_by_identifier[identifier]
        sb.identifier = identifier
    new_seqbuddies = [sb for identifier, sb in new_seqbuddies]
    return new_seqbuddies


def make_ids_unique(seqbuddy, sep="", padding=0):
    """
    Rename all repeat IDs
    Note: the edge case where a new ID creates a new conflict is not handled
    :param seqbuddy: SeqBuddy object
    :param sep: Pattern to place between id and number
    :param padding: Pad number with zeros
    :return: The modified SeqBuddy object
    """
    ids = OrderedDict()
    for rec in seqbuddy.records:
        ids.setdefault(rec.id, [])
        ids[rec.id].append(rec)

    records = []
    for key, recs in ids.items():
        if len(recs) > 1:
            for indx, rec in enumerate(recs):
                if re.match(rec.id, rec.description):
                    rec.description = rec.description[len(rec.id) + 1:]
                rec.id = "%s%s%s" % (rec.id, sep, str(indx + 1).zfill(padding))
        records += recs
    seqbuddy.records = records
    return seqbuddy


def map_features_nucl2prot(nuclseqbuddy, protseqbuddy, mode="key", quiet=False):
    """
    Applies cDNA/mRNA features to protein sequences
    :param nuclseqbuddy: A DNA SeqBuddy object with features to map
    :param protseqbuddy: A protein SeqBuddy
    :param mode: Specify how sequences should be matched up {list, key}
            - list = records are mapped in index order
            - key = records are converted to a dict and matched by key
    :param quiet: Suppress br._stderr messages
    :return: A protein SeqBuddy with the DNA SeqBuddy's features
    """
    def _feature_map(feature):
        if type(feature.location) == CompoundLocation:
            new_compound_location = []
            for sub_feature in feature.location.parts:
                sub_feature = _feature_map(SeqFeature(sub_feature))
                new_compound_location.append(sub_feature.location)
            feature.location = CompoundLocation(new_compound_location, feature.location.operator)

        elif type(feature.location) == FeatureLocation:
            start = feature.location.start / 3
            end = feature.location.end / 3
            location = FeatureLocation(floor(start), floor(end))
            feature.location = location

        else:
            raise TypeError("_feature_map requires a feature with either FeatureLocation or CompoundLocation, "
                            "not %s" % type(feature.location))  # This should be un-reachable because of clean_seq call
        return feature

    prot_copy, nucl_copy = make_copy(protseqbuddy), make_copy(nuclseqbuddy)
    clean_seq(prot_copy, skip_list="*")
    clean_seq(nucl_copy)

    stderr_written = False
    if mode == "list":
        if len(prot_copy) != len(nucl_copy):
            raise ValueError("The two input files do not contain the same number of sequences")

        record_map = list(zip(nucl_copy.records, prot_copy.records))

    elif mode == "key":
        prot_dict = prot_copy.to_dict()
        nucl_dict = nucl_copy.to_dict()

        record_map = []
        for seq_id, nucl_rec in nucl_dict.items():
            if seq_id not in prot_dict:
                stderr_written = True
                br._stderr("Warning: %s is in the cDNA file, but not in the protein file\n" % seq_id, quiet)
                continue
            else:
                record_map.append((nucl_rec, prot_dict[seq_id]))

        for seq_id, prot_rec in prot_dict.items():
            if seq_id not in nucl_dict:
                stderr_written = True
                br._stderr("Warning: %s is in the protein file, but not in the cDNA file\n" % seq_id, quiet)
                record_map.append((None, prot_rec))
    else:
        raise ValueError("'mode' must be either 'key' or 'position'.")

    for nucl_rec, prot_rec in record_map:
        # len(cds) or len(cds minus stop)
        if not nucl_rec:
            continue

        if len(prot_rec.seq) * 3 not in [len(nucl_rec.seq), len(nucl_rec.seq) - 3]:
            br._stderr("Warning: size mismatch between aa and nucl seqs for %s --> %s, %s\n" %
                       (nucl_rec.id, len(nucl_rec.seq), len(prot_rec.seq)), quiet)

        prot_feature_hashes = []
        for feat in prot_rec.features:
            prot_feature_hashes.append(md5(str(feat).encode()).hexdigest())

        for feat in nucl_rec.features:
            feat = _feature_map(feat)
            if md5(str(feat).encode()).hexdigest() not in prot_feature_hashes:
                prot_rec.features.append(feat)

    protseqbuddy.records = br.remap_gapped_features(prot_copy.records, protseqbuddy.records)
    if stderr_written:
        br._stderr("\n", quiet)
    return protseqbuddy


def map_features_prot2nucl(protseqbuddy, nuclseqbuddy, mode="key", quiet=False):
    """
    Applies protein features to cDNA/mRNA sequences
    :param protseqbuddy: A protein SeqBuddy object with features to map
    :param nuclseqbuddy: A DNA SeqBuddy
    :param mode: Specify how sequences should be matched up {list, key}
            - list = records are mapped in index order
            - key = records are converted to a dict and matched by key
    :param quiet: Suppress br._stderr messages
    :return: A DNA SeqBuddy with the protein SeqBuddy's features
    """
    def _feature_map(feature):
        if type(feature.location) == CompoundLocation:
            new_compound_location = []
            for sub_feature in feature.location.parts:
                sub_feature = _feature_map(SeqFeature(sub_feature))
                new_compound_location.append(sub_feature.location)
            feature.location = CompoundLocation(new_compound_location, feature.location.operator)

        elif type(feature.location) == FeatureLocation:
            start = feat.location.start * 3
            end = feat.location.end * 3
            location = FeatureLocation(start, end)
            feature.location = location

        else:
            raise TypeError("_feature_map requires a feature with either FeatureLocation or CompoundLocation, "
                            "not %s" % type(feature.location))  # This should be un-reachable because of clean_seq call
        return feature

    prot_copy, nucl_copy = make_copy(protseqbuddy), make_copy(nuclseqbuddy)
    clean_seq(prot_copy, skip_list="*")
    clean_seq(nucl_copy)

    stderr_written = False
    if mode == "list":
        if len(prot_copy) != len(nucl_copy):
            raise ValueError("The two input files do not contain the same number of sequences, try using 'key' mode.")

        record_map = list(zip(prot_copy.records, nucl_copy.records))

    elif mode == "key":
        prot_dict = prot_copy.to_dict()
        nucl_dict = nucl_copy.to_dict()

        record_map = []
        for seq_id, prot_rec in prot_dict.items():
            if seq_id not in nucl_dict:
                stderr_written = True
                br._stderr("Warning: %s is in the protein file, but not in the cDNA file\n" % seq_id, quiet)
                continue
            else:
                record_map.append((prot_rec, nucl_dict[seq_id]))

        for seq_id, dna_rec in nucl_dict.items():
            if seq_id not in prot_dict:
                stderr_written = True
                br._stderr("Warning: %s is in the cDNA file, but not in the protein file\n" % seq_id, quiet)
                record_map.append((None, dna_rec))

    else:
        raise ValueError("'mode' must be either 'key' or 'position'.")

    for prot_rec, nucl_rec in record_map:
        # len(cds) or len(cds minus stop)
        if not prot_rec:
            continue

        if len(prot_rec.seq) * 3 not in [len(nucl_rec.seq), len(nucl_rec.seq) - 3]:
            br._stderr("Warning: size mismatch between aa and nucl seqs for %s --> %s, %s\n" %
                       (prot_rec.id, len(prot_rec.seq), len(nucl_rec.seq)), quiet)

        dna_feature_hashes = []
        for feat in nucl_rec.features:
            dna_feature_hashes.append(md5(str(feat).encode()).hexdigest())

        for feat in prot_rec.features:
            feat = _feature_map(feat)
            prot_feature_hashes = [md5(str(feat).encode()).hexdigest()]
            # Need to account for strand orientation
            feat.strand = 0
            prot_feature_hashes.append(md5(str(feat).encode()).hexdigest())
            feat.strand = 1
            prot_feature_hashes.append(md5(str(feat).encode()).hexdigest())
            if not set(prot_feature_hashes) & set(dna_feature_hashes):
                nucl_rec.features.append(feat)

    nuclseqbuddy.records = br.remap_gapped_features(nucl_copy.records, nuclseqbuddy.records)
    if stderr_written:
        br._stderr("\n", quiet)
    return nuclseqbuddy


def max_records(seqbuddy, number=1):
    """
    Removes all sequences of length less than the maximum sequence length
    :param seqbuddy:
    :param number: How many records to return?
    :return:
    """
    seqbuddy = order_recs_by_len(seqbuddy, rev=True)
    seqbuddy.records = seqbuddy.records[:number]
    return seqbuddy


def merge(*seqbuddy):
    """
    Merges the feature lists of SeqBuddy objects
    :param seqbuddy: SeqBuddy objects to be combined
    :return: A new SeqBuddy object
    """

    def merge_records(rec1, rec2):
        # Deal with case issues and trailing stop codons
        rec1_seq = str(rec1.seq).lower()
        rec1_seq = re.sub(r"\*$", "", rec1_seq)
        rec2_seq = str(rec2.seq).lower()
        rec2_seq = re.sub(r"\*$", "", rec2_seq)
        if rec1_seq != rec2_seq:
            raise RuntimeError("Sequence mismatch for record '%s'" % rec1.id)
        already_present = False
        for feat2 in rec2.features:
            for feat1 in rec1.features:
                if str(feat2) == str(feat1):
                    already_present = True
                    break
            if not already_present:
                rec1.features.append(feat2)
            already_present = False
        return rec1

    seq_dict = {}
    for sb in seqbuddy:
        for seq_id, rec in sb.to_dict().items():
            if seq_id not in seq_dict:
                seq_dict[seq_id] = rec
            else:
                seq_dict[seq_id] = merge_records(seq_dict[seq_id], rec)

    seqbuddy = SeqBuddy([rec for _id, rec in seq_dict.items()], in_format=seqbuddy[0].in_format,
                        out_format=seqbuddy[0].out_format, alpha=seqbuddy[0].alpha)
    seqbuddy = order_ids(seqbuddy)
    seqbuddy = order_features_by_position(seqbuddy)
    return seqbuddy


def min_records(seqbuddy, number=1):
    """
    Removes all sequences of length greater than the minimum sequence length
    :param seqbuddy:
    :param number: How many records to return?
    :return:
    """
    seqbuddy = order_recs_by_len(seqbuddy)
    seqbuddy.records = seqbuddy.records[:number]
    return seqbuddy


def molecular_weight(seqbuddy):
    """
    Calculates the mass of each sequence in daltons
    :param seqbuddy: SeqBuddy object
    :return: SeqBuddy object with appended molecular_weights dictionary -
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
    dna = False
    output = OrderedDict([('masses_ss', []), ('masses_ds', []), ('ids', [])])
    aa_dict = amino_acid_weights
    if seqbuddy.alpha == "protein":
        aa_dict = amino_acid_weights
    elif seqbuddy.alpha == "DNA":
        aa_dict = deoxynucleotide_weights
        dna = True
    elif seqbuddy.alpha == "RNA":
        aa_dict = deoxyribonucleotide_weights
    for rec in seqbuddy.records:
        rec.mass_ds = 0
        rec.mass_ss = 0
        if seqbuddy.alpha == "protein":
            rec.mass_ss += 18.02  # molecular weight of a water molecule
        else:
            if dna:
                rec.mass_ss += 79.0  # molecular weight of 5' monophosphate in ssDNA
                rec.mass_ds += 157.9  # molecular weight of the 5' triphosphate in dsDNA
            else:
                rec.mass_ss += 159.0  # molecular weight of a 5' triphosphate in ssRNA
        for indx, value in enumerate(str(rec.seq).upper()):
            try:
                rec.mass_ss += aa_dict[value]
                if dna:
                    rec.mass_ds += aa_dict[value] + deoxynucleotide_weights[deoxynucleotide_compliments[value]]
            except KeyError:
                raise KeyError("Invalid residue '{0}' in record {1}. '{0}' is not valid a valid character in "
                               "{2}.".format(value, rec.id, str(seqbuddy.alpha)))
        output['masses_ss'].append(round(rec.mass_ss, 3))

        qualifiers = OrderedDict()
        if seqbuddy.alpha == "protein":
            qualifiers["peptide_value"] = round(rec.mass_ss, 3)
        elif dna:
            qualifiers["ssDNA_value"] = round(rec.mass_ss, 3)
            qualifiers["dsDNA_value"] = round(rec.mass_ds, 3)
            output['masses_ds'].append(round(rec.mass_ds, 3))
        elif seqbuddy.alpha == "RNA":
            qualifiers["ssRNA_value"] = round(rec.mass_ss, 3)
        output['ids'].append(rec.id)
        mw_feature = SeqFeature(location=FeatureLocation(start=1, end=len(rec.seq)), type='mw', qualifiers=qualifiers)
        rec.features.append(mw_feature)
    seqbuddy.molecular_weights = output
    return seqbuddy


def num_seqs(seqbuddy):
    """
    Counts the number of sequences in the SeqBuddy object
    :param seqbuddy: SeqBuddy object
    :return: The int number of sequences
    """
    return len(seqbuddy)


def order_features_alphabetically(seqbuddy, reverse=False):
    """
    Sorts features in alphabetical order
    :param seqbuddy: SeqBuddy object
    :param reverse: Specifies if the features should be sorted backwards
    :return: The SeqBuddy object with sorted features
    """
    for rec in seqbuddy.records:
        new_feature_list = [(feature.type, feature) for feature in rec.features]
        new_feature_list = sorted(new_feature_list, key=lambda x: x[0], reverse=reverse)
        new_feature_list = [feature[1] for feature in new_feature_list]
        rec.features = new_feature_list
    return seqbuddy


def order_features_by_position(seqbuddy, reverse=False):
    """
    Sorts features by the order in which they appear in the sequence
    :param seqbuddy: SeqBuddy object
    :param reverse: Specifies if the features should be sorted backwards
    :return: The SeqBuddy object with sorted features
    """
    for rec in seqbuddy.records:
        new_feature_list = [(int(feature.location.start), feature) for feature in rec.features]
        new_feature_list = sorted(new_feature_list, key=lambda x: x[0], reverse=reverse)
        new_feature_list = [feature[1] for feature in new_feature_list]
        rec.features = new_feature_list
    return seqbuddy


def order_ids(seqbuddy, reverse=False):
    """
    Sorts the sequences by ID, alphabetically
    :param seqbuddy: SeqBuddy object
    :param reverse: Reverses the sequence order
    :return: The sorted SeqBuddy object
    """
    def process_current_group(_group):
        if _group:
            _group = [(int(re.search("[0-9]+", x.id).group(0)), x) for x in _group]
            _group = sorted(_group, key=lambda l: l[0], reverse=reverse)
            _group = [x[1] for x in _group]
        return _group

    # Do an initial sort
    records = [(_rec.id, _rec) for _rec in seqbuddy.records]
    records = sorted(records, key=lambda x: x[0], reverse=reverse)
    records = [rec[1] for rec in records]

    # Now take into account numeric-order
    new_order = []
    current_group = []
    for rec in records:
        if re.search("[0-9]", rec.id):
            if current_group and re.match("([^0-9]*)", rec.id).group(0) != re.match("([^0-9]*)",
                                                                                    current_group[-1].id).group(0):
                new_order += process_current_group(current_group)
                current_group = []
            current_group.append(rec)
        else:
            new_order += process_current_group(current_group)
            current_group = []
            new_order.append(rec)
    new_order += process_current_group(current_group)
    seqbuddy.records = new_order
    return seqbuddy


def order_ids_randomly(seqbuddy, r_seed=None):
    """
    Reorders seqbuddy.records. The order will always be changed if more than 2 recs are fed in.
    :param seqbuddy: SeqBuddy object
    :param r_seed: Set the random generator seed value
    :return: The reordered SeqBuddy object
    """
    rand_gen = Random(r_seed)

    if len(seqbuddy) < 2:
        return seqbuddy
    elif len(seqbuddy) == 2:
        seqbuddy.records.reverse()
        return seqbuddy

    # make sure that every record isn't identical
    differences = False
    for indx, rec in enumerate(seqbuddy.records[1:]):
        if "%s%s" % (rec.id, rec.seq) != "%s%s" % (seqbuddy.records[indx - 1].id, seqbuddy.records[indx - 1].seq):
            differences = True
            break

    if not differences:
        return seqbuddy

    output = [None for _ in range(len(seqbuddy))]
    valve = br.SafetyValve(global_reps=1000)
    while valve.step("order_ids_randomly() was unable to reorder your sequences. This shouldn't happen, so please"
                     "contact the developers to let then know about this error."):
        sb_copy = make_copy(seqbuddy)
        for indx in range(len(sb_copy)):
            random_index = rand_gen.randint(1, len(sb_copy)) - 1
            output[indx] = (sb_copy.records.pop(random_index))
        if ["%s%s" % (rec.id, rec.seq) for rec in seqbuddy.records] != ["%s%s" % (rec.id, rec.seq) for rec in output]:
            break
        output = []

    seqbuddy.records = output
    return seqbuddy


def order_recs_by_len(seqbuddy, rev=False):
    """
    Sort the records by sequence size (smallest to largest)
    :param seqbuddy: SeqBuddy object
    :param rev: Sort from largest to smallest
    :return: The reordered SeqBuddy object
    """
    seqbuddy.records = sorted(seqbuddy.records, key=lambda x: len(x.seq), reverse=rev)
    return seqbuddy


def prepend_organism(seqbuddy, length=4):
    length = round(length)
    if length < 2:
        raise ValueError("Prefix length must be > 2")
    prefix_map = OrderedDict()
    new_prefix_map = {}

    for rec in seqbuddy.records:
        if "organism" in rec.annotations and rec.annotations['organism'] \
                and re.sub("[. ]", "", rec.annotations['organism']):
            organism = rec.annotations['organism'].split()[:2]
            org_name = " ".join(organism)
            if len(organism) > 1:
                prefix = organism[0][:1] + organism[1][:length-1]
            elif len(organism) == 1:
                prefix = organism[0][:length]
            if prefix in prefix_map:
                if org_name not in prefix_map[prefix]:
                    prefix_map[prefix].append(org_name)
                    continue
            else:
                prefix_map[prefix] = [org_name]

    for prefix in prefix_map:
        species = prefix_map[prefix]
        if len(species) > 1:
            for species_name in species:
                new_prefix = prefix + str(species.index(species_name) + 1)
                new_prefix_map[species_name] = new_prefix
        else:
            species_name = species[0]
            new_prefix_map[species_name] = prefix

    for rec in seqbuddy.records:
        if "organism" in rec.annotations and rec.annotations['organism'] \
                and re.sub("[. ]", "", rec.annotations['organism']):
            organism_name = " ".join(rec.annotations['organism'].split()[:2])
            if organism_name in new_prefix_map:
                final_prefix = new_prefix_map[organism_name] + "-"
                rec.id = final_prefix + rec.id
                rec.name = final_prefix + rec.name
        else:
            rec.id = "Unkn-" + rec.id
            rec.name = "Unkn-" + rec.name

    new_prefix_map = {v: k for k, v in new_prefix_map.items()}
    seqbuddy.prefix_map = new_prefix_map
    return seqbuddy


class PrositeScan(object):
    """
    Search for PROSITE scan motifs in sequences (via REST service)
    :param seqbuddy: Input seqbuddy object
    :param quiet: Suppress all stderr
    :return:
    """
    def __init__(self, seqbuddy, quiet=False, r_seed=None):
        self.seqbuddy = seqbuddy
        self.quiet = quiet
        self.r_seed = r_seed
        self.base_url = 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5'
        self.check_interval = 30

        # Prepend client specific agent string.
        user_agent = "SeqBuddy v%s" % VERSION.short()
        self.http_headers = {'User-Agent': user_agent}
        self.user_deets = br.config_values()

    def _rest_request(self, url, request_data=None):
        attempt = 1
        while True:
            try:
                # Set the User-agent.
                req = urllib.request.Request(url, None, self.http_headers)
                if request_data:
                    # Make the submission (HTTP POST).
                    req_h = urllib.request.urlopen(req, request_data)
                else:
                    # Make the request (HTTP GET).
                    req_h = urllib.request.urlopen(req)
                result = req_h.read().decode("utf-8")
                req_h.close()
                return result
            except urllib.error.HTTPError as err:
                if attempt > 5:
                    br._stderr("Error: Bad request at %s\n" % url, quiet=self.quiet)
                    raise err
                elif err.code in [400, 500]:
                    attempt += 1
                    time.sleep(attempt * 10)
                else:
                    raise err
            except urllib.request.http.client.RemoteDisconnected:
                if attempt > 5:
                    br._stderr("Error: Bad request at %s\n" % url, quiet=self.quiet)
                    raise
                else:
                    attempt += 1
                    time.sleep(attempt * 10)

    def _mc_run_prosite(self, _rec, args):
        out_file_path, lock = args
        try:
            if not self.user_deets["email"] or not re.search(r".+@.+\..+", self.user_deets["email"]):
                email = "buddysuite@gmail.com"
            else:
                email = self.user_deets["email"]

            params = {'sequence': str(_rec.seq).upper(),
                      'email': email,
                      'goterms': False,
                      'pathways': False,
                      'appl': "PrositeProfiles"}

            # Submit the job
            request_data = urllib.parse.urlencode(params)
            request_data = request_data.encode("utf-8")
            job_id = self._rest_request('%s/run/' % self.base_url, request_data)
            # ToDo: Consider including a timeout mechanism?
            result = 'PENDING'
            while result == 'RUNNING' or result == 'PENDING':
                result = self._rest_request('%s/status/%s' % (self.base_url, job_id))
                if result == 'RUNNING' or result == 'PENDING':
                    time.sleep(self.check_interval)

            try:
                result = self._rest_request('%s/result/%s/tsv' % (self.base_url, job_id)).strip()
            except urllib.error.HTTPError as err:
                if err.code in [400, 500]:
                    br._stderr("Error: Failed to retrieve %s\n" % _rec.id)
                    return
                else:
                    raise err

            feature_list = []
            if result:
                for feature in result.split("\n"):
                    feature = feature.split("\t")
                    start = feature[6]
                    end = feature[7]
                    feature_description = feature[-1]
                    feature = SeqFeature(FeatureLocation(int(start) - 1, int(end)), type="Region",
                                         qualifiers={"note": feature_description})
                    feature_list.append(feature)

            temp_seq = SeqBuddy([_rec], out_format="gb")
            temp_seq.records[0].features = feature_list
            temp_seq = order_features_by_position(temp_seq)

            with lock:
                with open(out_file_path, "a", encoding="utf-8") as out_file:
                    out_file.write("%s\n" % str(temp_seq))
        except KeyboardInterrupt:
            with lock:
                br._stderr("%s killed with KeyboardInterrupt\n" % _rec.id, quiet=self.quiet)
        return

    def run(self):
        self._rest_request(self.base_url)  # Confirm internet connection prior to multi-core loop

        temp_file = br.TempFile()
        hash_ids(self.seqbuddy, r_seed=self.r_seed)
        clean_seq(self.seqbuddy, skip_list="*")  # Clean once to make sure no wonky characters (no alignments)
        seqbuddy_nucl_copy = None

        if self.seqbuddy.alpha != "protein":
            seqbuddy_nucl_copy = make_copy(self.seqbuddy)
            translate_cds(self.seqbuddy)

        seqbuddy_prot_copy = make_copy(self.seqbuddy)
        clean_seq(self.seqbuddy)  # Clean again to strip *'s (added back in later)

        lock = Lock()
        try:
            br.run_multicore_function(self.seqbuddy.records, self._mc_run_prosite, [temp_file.path, lock],
                                      max_processes=30, out_type=sys.stderr, quiet=self.quiet)
        except KeyboardInterrupt:
            # Allow the user to grab whatever is finished and move on
            pass

        self.seqbuddy.records = SeqBuddy(temp_file.path).records
        new_records = []
        for rec in seqbuddy_prot_copy.records:
            for indx, rec2 in enumerate(self.seqbuddy.records):
                if rec.id == rec2.id:
                    new_records.append(rec2)
                    del self.seqbuddy.records[indx]
                    break

        self.seqbuddy.records = new_records

        find_pattern(seqbuddy_prot_copy, r"\*", include_feature=False)
        for rec in seqbuddy_prot_copy.records:
            for rec2 in self.seqbuddy.records:
                if rec.id == rec2.id:
                    for match in rec.buddy_data['find_patterns'][r"\*"]:
                        new_seq = str(rec2.seq)[:match[0]] + "*" + str(rec2.seq)[match[0]:]
                        rec2.seq = Seq(new_seq)
                        for f_indx, feat in enumerate(rec2.features):
                            if feat.location.start >= match[0]:
                                feat_loc = FeatureLocation(start=feat.location.start + 1, end=feat.location.end + 1)
                            elif feat.location.end > match[0]:
                                feat_loc = FeatureLocation(start=feat.location.start, end=feat.location.end + 1)
                            else:
                                feat_loc = feat.location

                            rec2.features[f_indx] = SeqFeature(location=feat_loc,
                                                               type=feat.type,
                                                               strand=feat.strand,
                                                               qualifiers=feat.qualifiers)

        if seqbuddy_nucl_copy:
            self.seqbuddy = map_features_prot2nucl(self.seqbuddy, seqbuddy_nucl_copy)
        else:
            self.seqbuddy = merge(seqbuddy_prot_copy, self.seqbuddy)

        self.seqbuddy.hash_map = seqbuddy_prot_copy.hash_map
        self.seqbuddy.reverse_hashmap()
        self.seqbuddy = order_ids(self.seqbuddy)
        return self.seqbuddy


def pull_random_recs(seqbuddy, count=1, r_seed=None):
    """
    Return a random record or subset of records (without replacement)
    :param seqbuddy: SeqBuddy object
    :param count: The number of random records to pull (int)
    :param r_seed: Set the random generator seed value
    :return: The original SeqBuddy object with only the selected records remaining
    """
    rand_gen = Random(r_seed)
    count = abs(count) if abs(count) <= len(seqbuddy) else len(seqbuddy)
    random_recs = []
    for _ in range(count):
        rand_index = rand_gen.randint(0, len(seqbuddy) - 1)
        random_recs.append(seqbuddy.records.pop(rand_index))
    seqbuddy.records = random_recs
    return seqbuddy


def pull_record_ends(seqbuddy, amount):
    """
    Retrieves subsequences from the ends of the sequences
    :param seqbuddy: SeqBuddy object
    :param amount: The number of residues to be pulled (negative numbers from rear)
    :return: The modified SeqBuddy object
    """
    amount = int(amount)
    seq_ends = []
    for rec in seqbuddy.records:
        if amount >= 0:
            rec.seq = Seq(str(rec.seq)[:amount])
            rec.features = br.shift_features(rec.features, 0, len(str(rec.seq)))

        else:
            shift = -1 * (len(str(rec.seq)) + amount) if abs(amount) <= len(str(rec.seq)) else 0
            rec.features = br.shift_features(rec.features, shift, len(str(rec.seq)))
            rec.seq = rec.seq[amount:]

        seq_ends.append(rec)

    seqbuddy.records = seq_ends
    return seqbuddy


def pull_recs(seqbuddy, regex, description=False):
    """
    Retrieves sequences with names/IDs matching a search pattern
    :param seqbuddy: SeqBuddy object
    :param regex: List of regex expressions or single regex
    :type regex: str list
    :param description: Allow search in description string
    :return: The modified SeqBuddy object
    """
    if type(regex) == str:
        regex = [regex]
    for indx, pattern in enumerate(regex):
        regex[indx] = ".*" if pattern == "*" else pattern

    regex = "|".join(regex)
    matched_records = []
    for rec in seqbuddy.records:
        if re.search(regex, rec.id) or re.search(regex, rec.name):
            matched_records.append(rec)
            continue
        if description and (re.search(regex, rec.description) or re.search(regex, str(rec.annotations))):
            matched_records.append(rec)
    seqbuddy.records = matched_records
    return seqbuddy


def pull_recs_with_feature(seqbuddy, regex):
    """
    Retrieves sequences with feature names/IDs matching a search pattern
    :param seqbuddy: SeqBuddy object
    :param regex: List of regex expressions or single regex
    :type regex: str list
    :return: The modified SeqBuddy object
    """
    if type(regex) == str:
        regex = [regex]
    for indx, pattern in enumerate(regex):
        regex[indx] = ".*" if pattern == "*" else pattern

    regex = "|".join(regex)
    matched_records = []
    for rec in seqbuddy.records:
        if re.search(regex, rec.description):
            matched_records.append(rec)
            continue
        for feat in rec.features:
            if re.search(regex, feat.type) or re.search(regex, feat.id):
                matched_records.append(rec)
                break
            breakout = False
            for key, quals in feat.qualifiers.items():
                if breakout:
                    break
                if re.search(regex, key):
                    matched_records.append(rec)
                    breakout = True
                    break
                for qual in quals:
                    if re.search(regex, qual):
                        matched_records.append(rec)
                        breakout = True
                        break
            if breakout:
                break

    seqbuddy.records = matched_records
    return seqbuddy


def purge(seqbuddy, threshold):
    """
    Deletes highly similar sequences
    :param seqbuddy: SeqBuddy object
    :param threshold: Sets the similarity threshold
    :return: The purged SeqBuddy object
    """
    # ToDo: Implement a way to return a certain # of seqs (i.e. auto-determine threshold)
    #    - This would probably be a different flag in the UI

    keep_dict = {}
    purged = []
    for query_id, match_list in bl2seq(seqbuddy).items():
        if query_id in purged:
            continue
        else:
            keep_dict[query_id] = []
            for subj_id in match_list:
                ident, length, evalue, bit_score = match_list[subj_id]

                if bit_score >= threshold:
                    purged.append(subj_id)
                    keep_dict[query_id].append(subj_id)
    new_records = []
    for rec in seqbuddy.records:
        if rec.id in keep_dict:
            _add_buddy_data(rec, "purge_set", keep_dict[rec.id])
            new_records.append(rec)

    seqbuddy.records = new_records
    return seqbuddy


def rename(seqbuddy, query, replace="", num=0, store_old_id=False):
    """
    Rename sequence IDs
    :param seqbuddy: SeqBuddy object
    :param query: The pattern to be searched for
    :param replace: The string to be substituted
    :param num: The maximum number of substitutions to make
    :param store_old_id: Keep a copy of the original ID in the description line
    :return: The modified SeqBuddy object
    """
    for rec in seqbuddy.records:
        new_name = br.replacements(rec.id, query, replace, num)
        if re.match(rec.id, rec.description):
            rec.description = rec.description[len(rec.id) + 1:]
        if store_old_id:
            rec.description = "%s %s" % (rec.id, rec.description)
        rec.id = new_name
        rec.name = new_name
    return seqbuddy


def replace_subsequence(seqbuddy, query, replacement=""):
    # ToDo: - Allow variable number of replacements
    #       - Add features to SeqRecords to denote substitutions
    #       - Consider adding case sensitive support.
    #       - Catch regex errors
    """
    :param seqbuddy: SeqBuddy object
    :param query: Regular expression (str)
    :param replacement: string
    :return:
    """
    query = "(?i)%s" % query
    for rec in seqbuddy.records:
        new_seq = re.sub(query, replacement, str(rec.seq))
        rec.seq = Seq(new_seq)
    return seqbuddy


def reverse_complement(seqbuddy):
    """
    Converts DNA/RNA sequences to their reverse complementary sequence.
    If protein sequence is provided, it is reverse transcribed, reverse complemented, and then transcribed again.
    :param seqbuddy: SeqBuddy object
    :return: The modified SeqBuddy object
    """
    for rec in seqbuddy.records:
        if rec.annotations["molecule_type"] == "RNA":
            rec.seq = rec.seq.reverse_complement_rna()
        else:
            rec.seq = rec.seq.reverse_complement()
        seq_len = len(rec.seq)
        shifted_features = [_feature_rc(feature, seq_len) for feature in rec.features]
        rec.features = shifted_features
    return seqbuddy


def rna2dna(seqbuddy):
    """
    Reverse-transcribes RNA into cDNA
    :param seqbuddy: SeqBuddy object
    :return: Modified SeqBuddy object
    """
    if seqbuddy.alpha != "RNA":
        raise TypeError("RNA sequence required, not %s." % seqbuddy.alpha)
    for rec in seqbuddy.records:
        rec.seq = Seq(str(rec.seq.back_transcribe()))
        rec.annotations["molecule_type"] = "DNA"
    seqbuddy.alpha = "DNA"
    return seqbuddy


def select_frame(seqbuddy, frame, add_metadata=True):
    """
    Changes the reading frame of the sequences
    :param seqbuddy: SeqBuddy object
    :param frame: The reading frame to shift to
    :param add_metadata: Appends the frame shift to a feature
    :return: The shifted SeqBuddy object
    """
    def reset_frame(_rec, _residues):
        _rec.seq = Seq("%s%s" % (_residues, str(_rec.seq)))
        for _feature in _rec.features:
            if "shift" in _feature.qualifiers:
                if type(_feature.location) != CompoundLocation:
                    _feature.location = FeatureLocation(_feature.location.start + int(_feature.qualifiers["shift"][0]),
                                                        _feature.location.end, _feature.location.strand)
                else:
                    _feature.location.parts[0] = FeatureLocation(_feature.location.start +
                                                                 int(_feature.qualifiers["shift"][0]),
                                                                 _feature.location.parts[0].end,
                                                                 _feature.location.strand)
                if frame == 1:
                    del _feature.qualifiers["shift"]

        _rec.features = br.shift_features(_rec.features, len(_residues), len(_rec.seq))
        _rec.description = str(re.sub(r"\(frame[23][A-Za-z]{1,2}\)", "", _rec.description)).strip()
        return _rec

    if seqbuddy.alpha == "protein":
        raise TypeError("Select frame requires nucleic acid, not protein.")

    for rec in seqbuddy.records:
        for indx, feature in enumerate(rec.features):
            if feature.location.start + 1 < frame and add_metadata:
                feature.qualifiers["shift"] = [str(feature.location.start + 1 - frame)]

        check_description = re.search(r"\(frame[23]([A-Za-z]{1,2})\)", rec.description)
        if check_description:
            rec = reset_frame(rec, check_description.group(1))

        rec.features = br.shift_features(rec.features, (frame - 1) * -1, len(rec.seq))
        if frame in [2, 3] and add_metadata:
            residues = str(rec.seq)[:frame - 1]
            rec.annotations["frame_shift"] = residues
            rec.description += " (frame%s%s)" % (frame, residues)
        rec.seq = Seq(str(rec.seq)[frame - 1:])
    return seqbuddy


def shuffle_seqs(seqbuddy, r_seed=None):
    """
    Randomly reorder the residues in each sequence
    :param seqbuddy: SeqBuddy object
    :param r_seed: Set the random generator seed value
    :return: The shuffled SeqBuddy object
    """
    rand_gen = Random(r_seed)
    for rec in seqbuddy.records:
        new_seq = list(str(rec.seq))
        rand_gen.shuffle(new_seq)
        rec.seq = Seq("".join(new_seq))
    return seqbuddy


def split_by_x_files(seqbuddy, file_number):
    """
    Splits SeqBuddy object with multiple records into set number of SeqBuddy objects
    :param seqbuddy: SeqBuddy object
    :param file_number: number of output Seqbuddy objects
    :return: A list of seqbuddy objects
    """
    counter = 1
    indx = 0
    uneven_counter = 1
    new_seqbuddies = dict()
    new_seqbuddies[indx] = []
    quotient = len(seqbuddy.records) // file_number
    remainder = len(seqbuddy.records) % file_number
    for rec in seqbuddy.records:
        seq_number = quotient + 1 if uneven_counter <= remainder else quotient
        if counter <= seq_number:
            new_seqbuddies[indx].append(rec)
            counter += 1
        elif counter > seq_number:
            indx += 1
            uneven_counter += 1
            new_seqbuddies[indx] = []
            new_seqbuddies[indx].append(rec)
            counter = 2

    split_buddies = []
    for key in new_seqbuddies:
        new_buddy = make_copy(seqbuddy)
        new_buddy.records = new_seqbuddies[key]
        split_buddies.append(new_buddy)
    return split_buddies


def split_by_x_seqs(seqbuddy, seq_number):
    """
    Splits SeqBuddy object with multiple records into multiple SeqBuddy objects containing set number of records
    :param seqbuddy: SeqBuddy object
    :param seq_number: number of records per output seqbuddy object
    :return: A list of seqbuddy objects
    """
    counter = 1
    indx = 0
    new_seqbuddies = dict()
    new_seqbuddies[indx] = []
    for rec in seqbuddy.records:
        if counter <= seq_number:
            new_seqbuddies[indx].append(rec)
            counter += 1
        elif counter > seq_number:
            indx += 1
            new_seqbuddies[indx] = []
            new_seqbuddies[indx].append(rec)
            counter = 2

    split_buddies = []
    for key in new_seqbuddies:
        new_buddy = make_copy(seqbuddy)
        new_buddy.records = new_seqbuddies[key]
        split_buddies.append(new_buddy)
    return split_buddies


def tail(seqbuddy, count=1):
    """
    Pull the last record(s) in a SeqBuddy object
    :param seqbuddy: SeqBuddy object
    :param count: Number of records to return
    :return: The modified SeqBuddy object
    """
    seqbuddy.records = seqbuddy.records[-1*count:]
    return seqbuddy


def taxonomic_breakdown(seqbuddy, max_depth=5):
    max_depth = abs(max_depth)
    taxonomy_dict = OrderedDict([("rec_ids", [])])
    seqbuddy = make_copy(seqbuddy)
    for rec_id, rec in seqbuddy.to_dict().items():
        if "taxonomy" not in rec.annotations or not rec.annotations['taxonomy']:
            taxonomy_dict.setdefault("Unknown", OrderedDict([("rec_ids", [])]))
            taxonomy_dict["rec_ids"].append(rec.id)
            taxonomy_dict["Unknown"]["rec_ids"].append(rec.id)
        else:
            if "organism" in rec.annotations:
                rec.annotations['organism'] = re.sub("[ ]+", " ", rec.annotations['organism'])
                if rec.annotations['organism']:
                    organism = rec.annotations['organism'].split()[1:]
                    organism = " ".join(organism)
                    rec.annotations['taxonomy'].append(organism)

            taxonomy_dict = _tb_add_tax(rec_id, rec.annotations['taxonomy'], taxonomy_dict)
    return _tb_to_string(taxonomy_dict, depth=0, max_depth=max_depth)


def _tb_add_tax(rec_id, taxonomy, tax_dict):
    tax_dict.setdefault("rec_ids", [])
    tax_dict["rec_ids"].append(rec_id)
    if len(taxonomy) >= 1:
        tax_dict.setdefault(taxonomy[0], OrderedDict())
        tax_dict[taxonomy[0]] = _tb_add_tax(rec_id, taxonomy[1:], tax_dict[taxonomy[0]])
    return tax_dict


def _tb_to_string(tax_dict, depth, max_depth):
    output = "Total: %s\n\n" % len(tax_dict["rec_ids"]) if depth == 0 else ""
    if not max_depth or depth < max_depth:
        del tax_dict["rec_ids"]
        tax_dict = sorted([(taxon, sub_dict) for taxon, sub_dict in tax_dict.items()],
                          key=lambda x: len(x[1]))
        for taxon, sub_dict in tax_dict:
            output += "{0}{1}    {2}\n".format(" |" * depth, taxon, len(sub_dict['rec_ids']))
            output += _tb_to_string(sub_dict, depth + 1, max_depth)
    return output


def translate6frames(seqbuddy):
    """
    Translates a nucleotide sequence into a protein sequence across all six reading frames.
    :param seqbuddy: SeqBuddy object
    :return: The translated SeqBuddy object
    """
    frame1, frame2, frame3 = make_copy(seqbuddy), make_copy(seqbuddy), make_copy(seqbuddy)
    seqbuddy = reverse_complement(seqbuddy)

    rframe1, rframe2, rframe3 = make_copy(seqbuddy), make_copy(seqbuddy), make_copy(seqbuddy)

    frame2 = select_frame(frame2, 2, add_metadata=False)
    frame3 = select_frame(frame3, 3, add_metadata=False)
    rframe2 = select_frame(rframe2, 2, add_metadata=False)
    rframe3 = select_frame(rframe3, 3, add_metadata=False)

    frame1 = translate_cds(frame1, quiet=True)
    frame2 = translate_cds(frame2, quiet=True)
    frame3 = translate_cds(frame3, quiet=True)
    rframe1 = translate_cds(rframe1, quiet=True)
    rframe2 = translate_cds(rframe2, quiet=True)
    rframe3 = translate_cds(rframe3, quiet=True)

    output = []
    for i in range(len(frame1)):
        frame1.records[i].id = "%s_f1" % frame1.records[i].id
        frame2.records[i].id = "%s_f2" % frame2.records[i].id
        frame3.records[i].id = "%s_f3" % frame3.records[i].id
        rframe1.records[i].id = "%s_rf1" % rframe1.records[i].id
        rframe2.records[i].id = "%s_rf2" % rframe2.records[i].id
        rframe3.records[i].id = "%s_rf3" % rframe3.records[i].id

        output += [frame1.records[i], frame2.records[i], frame3.records[i],
                   rframe1.records[i], rframe2.records[i], rframe3.records[i]]

    seqbuddy.records = SeqBuddy(output).records
    return seqbuddy


def translate_cds(seqbuddy, quiet=False, alignment=False):
    """
    Translates a nucleotide sequence into a protein sequence.
    :param seqbuddy: SeqBuddy object
    :param quiet: Suppress the errors thrown by translate(cds=True)
    :param alignment: If the incoming sequence has gaps you want maintained, set to True. Otherwise they will be cleaned
    :return: The translated SeqBuddy object
    """
    if seqbuddy.alpha == "protein":
        raise TypeError("Protein sequence cannot be translated.")

    codon_dict = {'---': '-', '--A': '-', '--C': '-', '--G': '-', '--T': '-', '-A-': '-', '-C-': '-', '-G-': '-',
                  '-T-': '-', 'A--': '-', 'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T',
                  'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S', 'ATA': 'I', 'ATC': 'I',
                  'ATG': 'M', 'ATT': 'I', 'C--': '-', 'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 'CCA': 'P',
                  'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L',
                  'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'G--': '-', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
                  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                  'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'T--': '-', 'TAA': '*', 'TAC': 'Y', 'TAG': '*',
                  'TAT': 'Y', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TGA': '*', 'TGC': 'C', 'TGG': 'W',
                  'TGT': 'C', 'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}

    if not alignment:
        clean_seq(seqbuddy)

    if seqbuddy.alpha == "RNA":
        rna2dna(seqbuddy)

    translated_sb = make_copy(seqbuddy)
    uppercase(translated_sb)
    for rec in translated_sb.records:
        if rec.annotations["molecule_type"] == "protein":
            warnings.warn("Record %s is protein." % rec.id)

        new_seq = [""] * int(ceil(len(rec) / 3))
        new_seq_indx = 0
        old_seq = str(rec.seq)
        for codon in [old_seq[::1][i:i + 3][::1] for i in range(0, len(old_seq), 3)][::1]:
            if len(codon) == 3:
                if codon in codon_dict:
                    new_seq[new_seq_indx] = codon_dict[codon]
                else:
                    new_seq[new_seq_indx] = "N"
                new_seq_indx += 1

        new_seq = "".join(new_seq) if new_seq and new_seq[-1] != "" else "".join(new_seq[:-1])
        new_seq = Seq(new_seq)
        rec.seq = new_seq
        rec.features = []

    # Can't updated records to protein until everything is converted because Nexus format will explode for some reason
    for rec in translated_sb.records:
        rec.annotations["molecule_type"] = "protein"

    map_features_nucl2prot(seqbuddy, translated_sb, mode="list", quiet=quiet)
    for indx, rec in enumerate(translated_sb.records):
        seqbuddy.records[indx] = rec
    seqbuddy.alpha = "protein"
    return seqbuddy


def transmembrane_domains(seqbuddy, job_ids=None, quiet=False, keep_temp=None):
    """
    Access the TOPCONS service to annotate transmembrane domains
    :param seqbuddy: SeqBuddy object
    :param job_ids: If the sequences in SeqBuddy object have previously been run, pickup from the download step by
                  supplying the TOPCONS reference job id
    :type job_ids: str list
    :param quiet: Suppress all output
    :param keep_temp: Save all output files generated by TOPCONS
    """
    # ToDo: Write hashmaps to the current working dir instead of site-packages
    try:
        from suds.client import Client
    except ImportError:
        raise ImportError("Please install the 'suds' package to run transmembrane_domains:\n\n$ pip install suds-py3")

    def dl_progress(count, block_size, total_size):
        percent = count * block_size * 100 / total_size
        valve.test(percent)
        printer.write("Retrieving job %s of %s: ... %d%%" % (len(results) + len(failed) + 1, len(jobs), int(percent)))

    wsdl_url = "http://v2.topcons.net/pred/api_submitseq/?wsdl"
    max_seqsize = 9 * 1024 * 1024
    max_filesize = 1024 * 1024

    printer = br.DynamicPrint(out_type="stderr", quiet=quiet)
    temp_dir = br.TempDir()
    config = br.config_values()

    job_dir = "{0}{1}topcons".format(temp_dir.path, os.path.sep) if not config["data_dir"] \
        else "%s%stopcons" % (config["data_dir"], os.path.sep)
    try:
        os.makedirs(job_dir, exist_ok=True)
        open("%s%sdel" % (job_dir, os.path.sep), "w").close()
        os.remove("%s%sdel" % (job_dir, os.path.sep))
    except PermissionError:
        job_dir = "{0}{1}topcons".format(temp_dir.path, os.path.sep)
        os.makedirs(job_dir, exist_ok=True)

    printer.write("Cleaning sequences")
    clean_seq(seqbuddy, skip_list="*")

    job_ids = [] if not job_ids else job_ids
    job_ids = [job_ids] if type(job_ids) == str else job_ids

    printer.write("Hashing sequence IDs")

    hash_map = OrderedDict()
    seqbuddy_copy = make_copy(seqbuddy)
    seqbuddy.out_format = "fasta"

    printer.write("Stripping meta data")
    delete_metadata(seqbuddy)

    if seqbuddy.alpha != "protein":
        printer.write("Translating to protein")
        translate_cds(seqbuddy)

    jobs = []

    if job_ids:
        seqbuddy_recs = []
        for jobid in job_ids:
            if not os.path.isfile("%s%s%s.hashmap" % (job_dir, os.path.sep, jobid)):
                printer.clear()
                error_message = "SeqBuddy does not have the necessary hash-map to process job id '%s'. This could be" \
                                " a job id typo, a configuration issue, or you may be attempting to access a job" \
                                " submitted by a different computer. See the GitHub wiki for further details" \
                                " https://github.com/biologyguy/BuddySuite/wiki/SB-Transmembrane-domains" % jobid
                raise FileNotFoundError(error_message)

            with open("%s%s%s.hashmap" % (job_dir, os.path.sep, jobid), "r") as ifile:
                jobs.append({"type": "previous", "hash_map": OrderedDict(), "records": []})
                for line in ifile:
                    line = line.strip().split("\t")
                    for indx, rec in enumerate(seqbuddy.records):
                        if rec.id == line[1]:
                            hash_map[line[0]] = line[1]
                            rec = SeqBuddy([rec])
                            rename(rec, line[1], line[0])
                            jobs[-1]["records"].append(rec)
                            jobs[-1]["hash_map"][line[0]] = line[1]
                            seqbuddy_recs.append(rec.records[0])
                            del seqbuddy.records[indx]
                            break

    if len(seqbuddy):  # len() is in case job ids are passed in, but not all sequences are covered
        hash_ids(seqbuddy)
        for _hash, rec_id in seqbuddy.hash_map.items():
            hash_map[_hash] = rec_id

        jobs.append({"type": "new", "hash_map": OrderedDict(), "records": []})
        seqbuddy_size = len(seqbuddy.records)
        printer.write("Preparing jobs for upload (0 of %s records processed)" % seqbuddy_size)
        rec_string = ""
        for indx, rec in enumerate(seqbuddy.records):
            printer.write("Preparing jobs for upload (%s of %s records processed)" % (indx + 1, seqbuddy_size))
            rec_string += rec.format("fasta")
            if len(rec.format("fasta") + rec_string) <= max_filesize:
                jobs[-1]["records"].append(rec)
                jobs[-1]["hash_map"][rec.id] = seqbuddy.hash_map[rec.id]
            else:
                if len(rec.format("fasta")) > max_seqsize:
                    printer.clear()
                    raise ValueError("Record '%s' is too large to send to TOPCONS. Max record size is 9Mb" %
                                     seqbuddy.hash_map[rec.id])
                jobs.append({"type": "new", "hash_map": OrderedDict({rec.id: seqbuddy.hash_map[rec.id]}),
                             "records": [rec]})
                rec_string = ""

        for job in jobs:
            if job["type"] == "new":
                job["records"] = SeqBuddy(job["records"], out_format="fasta")
                job["records"].hash_map = job["hash_map"]

        for indx, job in enumerate(jobs):
            if job["type"] == "new":
                printer.write("Uploading job %s of %s" % (indx + 1, len(jobs)))
                myclient = Client(wsdl_url, cache=None)
                ret_value = myclient.service.submitjob(str(job["records"]), "", "", "")
                if len(ret_value) >= 1:
                    jobid, result_url, numseq_str, errinfo, warninfo = ret_value[0][:5]
                    if jobid not in ["None", ""]:
                        printer.clear()
                        br._stderr("Job '%s' submitted\n" % jobid, quiet=quiet)
                        job_ids.append(jobid)
                        temp_dir.subdir(jobid)
                        with open("%s%s%s.hashmap" % (job_dir, os.path.sep, jobid), "w", encoding="utf-8") as ofile:
                            ofile.write(job["records"].print_hashmap())
                    else:
                        printer.clear()
                        raise ConnectionError("Failed to submit TOPCONS job.\n%s" % errinfo)
                else:
                    printer.clear()
                    raise ConnectionError("Failed to submit TOPCONS job. Are you connected to the internet?")

    # Need to match up all hashed ids in seqbuddy_copy for downstream stuff
    records = []
    for _hash, rec_id in hash_map.items():
        for indx, rec in enumerate(seqbuddy_copy.records):
            if rec.id == rec_id:
                rec.id = _hash
                records.append(rec)
                del seqbuddy_copy.records[indx]
                break
    seqbuddy_copy.records = records

    # Stops are converted to Xs by TOPCONS, so find them now for later replacement
    stop_positions = {}
    if seqbuddy_copy.alpha == "protein":
        printer.write("Identifying stop codons")
        seqbuddy_copy = find_pattern(seqbuddy_copy, r"\*", include_feature=False)
        stop_positions = {rec.id: rec.buddy_data['find_patterns'][r'\*'] for rec in seqbuddy_copy.records}

    results = []
    failed = []
    wait = True
    delay = 1
    while len(results) + len(failed) != len(jobs):
        if wait:
            delay = 1.5 * delay if delay < 300 else delay
            for i in range(round(delay)):
                slash = ["/", "—", "\\", "|"]
                printer.write("Waiting for TOPCONS results (%s of %s jobs complete) %s " %
                              (len(results) + len(failed), len(jobs), slash[i % 4]))
                time.sleep(1)

        wait = True
        for indx, jobid in enumerate(job_ids):
            printer.write("Checking job %s of %s (remaining)" % (indx + 1, len(jobs) - len(results) - len(failed)))
            myclient = Client(wsdl_url, cache=None)
            ret_value = myclient.service.checkjob(jobid)
            if len(ret_value) >= 1:
                status, result_url, errinfo = ret_value[0][:3]
                if status == "Failed":
                    printer.clear()
                    raise ConnectionError("Job failed...\nServer message: %s" % errinfo)
                elif status == "Finished":
                    outfile = "%s%s%s.zip" % (temp_dir.path, os.path.sep, jobid)
                    printer.write("Retrieving job %s of %s" % (len(results) + len(failed) + 1, len(jobs)))
                    tries = 1
                    while True:
                        valve = br.SafetyValve(state_reps=25)
                        try:
                            urllib.request.urlretrieve(result_url, filename=outfile, reporthook=dl_progress)
                            break

                        except RuntimeError:
                            printer.write("Download stalled, restarting... try %s of 5" % tries)
                            time.sleep(5 * tries)
                            tries += 1

                        except urllib.error.ContentTooShortError:
                            printer.write("Download file wrong size, restarting... try %s of 5" % tries)
                            time.sleep(5 * tries)
                            tries += 1

                        except urllib.error.HTTPError:
                            printer.write("HTTPError reported, restarting... try %s of 5" % tries)
                            time.sleep(150 * tries)
                            tries += 1

                        if tries >= 5:
                            break

                    if os.path.exists(outfile) and jobid not in results:
                        results.append(jobid)

                    else:
                        br._stderr("\nError: Failed to download TOPCONS job {0} after 5 attempts. "
                                   "The data will be saved on the server for manual retrieval.\n"
                                   "A sequence name hash-map has been saved to {0}.hashmap".format(jobid), quiet=quiet)
                        with open("%s.hashmap" % jobid, "w", encoding="utf-8") as ofile:
                            ofile.write(seqbuddy.print_hashmap())
                        failed.append(jobid)

                    del job_ids[indx]
                    wait = False
                    break

                elif status == "None":
                    printer.clear()
                    raise ConnectionError("The job seems to have been lost by the server.\n%s" % errinfo)

    for indx, jobid in enumerate(results):
        printer.write("Extracting results %s of %s (%s)" % (indx + 1, len(results), jobid))
        with zipfile.ZipFile("%s%s%s.zip" % (temp_dir.path, os.path.sep, jobid)) as zf:
            zf.extractall(temp_dir.path)
        os.remove("%s%s%s.zip" % (temp_dir.path, os.path.sep, jobid))

    printer.write("Processing results...")
    records = []
    for jobid in results:
        with open("{0}{2}{1}{2}query.result.txt".format(temp_dir.path, jobid, os.path.sep),
                  "r", encoding="utf-8") as ifile:
            topcons = ifile.read()

        topcons = topcons.split(
            "##############################################################################")[2:-1]

        for rec in topcons:
            printer.write("Processing results... %s" % len(records))
            seq_id = re.search("Sequence name: (.+)", rec).group(1).strip()
            if seq_id not in hash_map:
                continue
            seq = re.search("Sequence:\n([A-Z]+)", rec).group(1).strip()
            alignment = ""
            for algorithm in ["TOPCONS", "OCTOPUS", "Philius", "PolyPhobius", "SCAMPI", "SPOCTOPUS"]:
                if re.search(r"%s predicted topology:\n\*\*\*No topology could be "
                             r"produced with this method\*\*\*" % algorithm, rec):
                    continue
                top_file = re.search("%s predicted topology:\n([ioMSs]+)" % algorithm, rec).group(1).strip()
                top_file = re.sub("[^M]", "i", top_file)
                alignment += ">%s\n%s\n\n" % (algorithm, top_file)

            cons_seq = SeqBuddy(">%s\n%s\n" % (seq_id, seq), out_format="genbank")
            if alignment:
                alignment = Alb.AlignBuddy(alignment)
                Alb.consensus_sequence(alignment)
                counter = 1
                for tmd in re.finditer("([MX]+)", str(alignment.records()[0].seq)):
                    annotate(cons_seq, "TMD%s" % counter, "%s-%s" % (tmd.start(), tmd.end()))
                    counter += 1
            records.append(cons_seq.records[0])

    printer.write("Creating new SeqBuddy object")
    seqbuddy.records = SeqBuddy(records).records

    if keep_temp:
        printer.write("Preparing TOPCONS files to be saved")
        for _root, dirs, files in br.walklevel(temp_dir.path):
            for file in files:
                with open("%s%s%s" % (_root, os.path.sep, file), "r", encoding="utf-8") as ifile:
                    contents = ifile.read()
                for _hash, _id in hash_map.items():
                    contents = re.sub(_hash, _id, contents)
                with open("%s%s%s" % (_root, os.path.sep, file), "w", encoding="utf-8") as ofile:
                    ofile.write(contents)

        printer.write("Saving TOPCONS files")
        os.makedirs(keep_temp, exist_ok=True)
        _root, dirs, files = next(br.walklevel(temp_dir.path))
        for file in files:
            shutil.copyfile("%s%s%s" % (_root, os.path.sep, file), "%s%s%s" % (keep_temp, os.path.sep, file))
        for _dir in dirs:
            shutil.copytree("%s%s%s" % (_root, os.path.sep, _dir), "%s%s%s" % (keep_temp, os.path.sep, _dir))

    printer.write("Merging sequence features")
    if seqbuddy_copy.alpha != "protein":
        seqbuddy = map_features_prot2nucl(seqbuddy, seqbuddy_copy)
    else:
        for indx, rec in enumerate(seqbuddy.records):
            printer.write("Merging sequence features --> %s of %s" % (indx, len(seqbuddy)))
            matches = stop_positions[rec.id]
            for match in matches:
                new_seq = str(rec.seq)[:match[0]] + "*" + str(rec.seq)[match[0] + 1:]
                rec.seq = Seq(new_seq)
        printer.write("Merging sequence features --> Executing merge()")
        seqbuddy = merge(seqbuddy_copy, seqbuddy)

    printer.write("Restoring original sequence names")
    for _hash, seq_id in hash_map.items():
        rename(seqbuddy, _hash, seq_id)

    printer.write("Sorting sequences by ID")
    seqbuddy = order_ids(seqbuddy)

    printer.write("************** Complete **************")
    printer.new_line(2)
    return seqbuddy


def uppercase(seqbuddy):
    """
    Converts all sequence characters to uppercase.
    :param seqbuddy: SeqBuddy object
    :return: The modified SeqBuddy object
    """
    for rec in seqbuddy.records:
        rec.seq = Seq(str(rec.seq).upper())
    return seqbuddy


# ################################################# COMMAND LINE UI ################################################## #
def argparse_init():
    # Catching params to prevent weird collisions with 3rd party arguments
    if '--blast' in sys.argv:  # Only blast at the moment, but other flags may come up in the future.
        sys.argv[sys.argv.index('--blast')] = '-bl'
    if '-bl' in sys.argv:
        sb_flag_indx = sys.argv.index('-bl')
        if len(sys.argv) > sb_flag_indx + 1:
            extra_args = None
            for indx, param in enumerate(sys.argv[sb_flag_indx + 1:]):
                if param in ["-a", "--alpha", "-f", "--in_format", "-i", "--in_place", "-k", "--keep_temp", "-o",
                             "--out_format", "-q", "--quiet", "-t", "--test"]:
                    extra_args = sb_flag_indx + 1 + indx
                    break

            if extra_args == sb_flag_indx + 1 or extra_args == sb_flag_indx + 2:
                # No conflicts possible, continue on your way
                pass

            elif len(sys.argv) > sb_flag_indx + 2 or extra_args:
                # There must be optional arguments being passed into the blast tool
                sys.argv[sb_flag_indx + 2] = " %s" % sys.argv[sb_flag_indx + 2].rstrip()

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
''',)

    br.flags(parser, ("sequence", "Supply file path(s) or raw sequence. If piping sequences "
                                  "into SeqBuddy this argument must be left blank."),
             br.sb_flags, br.sb_modifiers, VERSION)

    in_args = parser.parse_args()
    br.check_garbage_flags(in_args, "SeqBuddy")

    seqbuddy = []
    seq_set = ""

    if in_args.out_format and in_args.out_format.lower() not in OUTPUT_FORMATS:
        br._stderr("Error: Output type '%s' is not recognized/supported\n" % in_args.out_format)
        sys.exit()

    if in_args.guess_alphabet or in_args.guess_format:
        return in_args, SeqBuddy

    try:
        for seq_set in in_args.sequence:
            if isinstance(seq_set, TextIOWrapper) and seq_set.buffer.raw.isatty():
                br._stderr("Warning: No input detected so SeqBuddy is aborting...\n"
                           "For more information, try:\n%s --help\n" % sys.argv[0])
                sys.exit()
            seq_set = SeqBuddy(seq_set, in_args.in_format, in_args.out_format, in_args.alpha)
            seqbuddy += seq_set.records

        seqbuddy = SeqBuddy(seqbuddy, seq_set.in_format, seq_set.out_format, seq_set.alpha)
    except br.GuessError as err:
        br._stderr("GuessError: %s\n" % err, in_args.quiet)
        sys.exit()
    except Exception as err:
        br._stderr("Error: Unable to process input file(s)\n%s\n" % err)
        sys.exit()

    if in_args.restrict:
        seqbuddy.hidden_recs = delete_records(make_copy(seqbuddy), in_args.restrict[0]).records
        seqbuddy = pull_recs(seqbuddy, in_args.restrict[0])

    in_args.random_seed = None if not in_args.random_seed else in_args.random_seed
    return in_args, seqbuddy


def command_line_ui(in_args, seqbuddy, skip_exit=False, pass_through=False):  # ToDo: Convert to a class
    # ############################################ INTERNAL FUNCTIONS ################################################ #
    def _print_recs(_seqbuddy):
        try:
            _output = str(_seqbuddy)
        except TypeError as _err:  # Catching biopython/nexus issue
            if " unexpected keyword argument 'interleave'" in str(_err) and _seqbuddy.out_format.startswith("nexus"):
                br._stderr("WARNING: Attempting to write NEXUS format '%s' with an old version of "
                           "BioPython. Switching to default nexus.\n" % _seqbuddy.out_format, in_args.quiet)
                _seqbuddy.out_format = "nexus"
                _output = str(_seqbuddy)
            else:
                raise _err

        if in_args.test:
            br._stderr("*** Test passed ***\n", in_args.quiet)
            pass

        elif in_args.in_place:
            _in_place(_output, in_args.sequence[0])

        else:
            br._stdout("{0}\n".format(_output.rstrip()))

    def _in_place(_output, file_path):
        if not os.path.exists(file_path):
            br._stderr("Warning: The -i flag was passed in, but the positional argument doesn't seem to be a "
                       "file. Nothing was written.\n", in_args.quiet)
            br._stderr("%s\n" % _output.strip(), in_args.quiet)
        else:
            with open(os.path.abspath(file_path), "w", encoding="utf-8") as _ofile:
                _ofile.write(_output)
            br._stderr("File overwritten at:\n%s\n" % os.path.abspath(file_path), in_args.quiet)

    def _raise_error(_err, tool, check_string=None):
        """
        Custom error handling for Seqbuddy UI
        :param _err: python exception object
        :param tool: name of tool that the error is being raised from
        :param check_string: expected error message
        :return:
        """
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
        _exit(tool)

    def _exit(tool, skip=skip_exit):
        if skip:
            return
        usage = br.Usage()
        memory_footprint = 0 if type(seqbuddy) != SeqBuddy else seqbuddy.memory_footprint
        usage.increment("SeqBuddy", VERSION.short(), tool, memory_footprint)
        usage.save()
        sys.exit()

    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Amend metadata
    if in_args.amend_metadata:
        args = in_args.amend_metadata[0]
        if len(args) == 1:
            attr, sub_value, regex = args[0], "", ".+"
        elif len(args) == 2:
            attr, sub_value, regex = args[0], args[1], ".+"
        else:
            attr, sub_value, regex = args[:3]

        try:
            _print_recs(amend_metadata(seqbuddy, attr.lower(), sub_value, regex))
        except ValueError as err:
            _raise_error(err, "amend_metadata", "Topology values are limited to")
        _exit("amend_metadata")

    # Annotate feature
    if in_args.annotate:
        # _type, location, strand=None, qualifiers=None, pattern=None
        genbank_features = ['assembly_gap', 'attenuator', 'C_region', 'CAAT_signal', 'CDS', 'centromere', 'D-loop',
                            'D_segment', 'enhancer', 'exon', 'gap', 'GC_signal', 'gene', 'iDNA', 'intron', 'J_segment',
                            'LTR', 'mat_peptide', 'misc_binding', 'misc_difference', 'misc_feature', 'misc_recomb',
                            'misc_RNA', 'misc_signal', 'misc_structure', 'mobile_element', 'modified_base', 'mRNA',
                            'ncRNA', 'N_region', 'old_sequence', 'operon', 'oriT', 'polyA_signal', 'polyA_site',
                            'precursor_RNA', 'prim_transcript', 'primer_bind', 'promoter', 'protein_bind', 'RBS',
                            'regulatory', 'repeat_region', 'rep_origin', 'rRNA', 'S_region', 'sig_peptide', 'source',
                            'stem_loop', 'STS', 'TATA_signal', 'telomere', 'terminator', 'tmRNA', 'transit_peptide',
                            'tRNA', 'unsure', 'V_region', 'V_segment', 'variation', "3'UTR", "5'UTR", '-10_signal',
                            '-35_signal']

        strand_types = ['+', 'plus', 'sense', 'pos', '1', 1, '-', 'minus',
                        'anti', 'anti-sense', 'neg', '-1', -1, '0', 0]
        feature_attrs = {"strand": None, "qualifiers": [], "pattern": []}

        def duck_type(_arg, **kwargs):  # Feed feature_attrs into kwargs
            if _arg in strand_types:
                kwargs["strand"] = _arg
            elif re.search(".+[=:].+", _arg):
                kwargs["qualifiers"].append(_arg)
            else:
                kwargs["pattern"].append(_arg)
            return kwargs

        if len(in_args.annotate) < 2:
            _raise_error(AttributeError("Too few parameters provided. Please provide at least a feature type "
                                        "and a location.\n"), "annotate")
        ftype = in_args.annotate[0]
        if ftype not in genbank_features:
            br._stderr("Warning: The provided annotation type is not part of the GenBank format standard\n\n",
                       in_args.quiet)

        if len(ftype) > 16:
            br._stderr("Warning: Feature type is longer than 16 characters and "
                       "will be truncated if printed to GenBank/EMBL format\n\n", in_args.quiet)

        flocation = in_args.annotate[1]

        if len(in_args.annotate) >= 3:
            for _next_arg in in_args.annotate[2:]:
                feature_attrs = duck_type(_next_arg, **feature_attrs)

            if not feature_attrs["qualifiers"]:
                feature_attrs["qualifiers"] = None
            feature_attrs["pattern"] = br.clean_regex(feature_attrs["pattern"], in_args.quiet)
            if not feature_attrs["pattern"]:
                feature_attrs["pattern"] = None
        try:
            seqbuddy = annotate(seqbuddy, ftype, flocation, **feature_attrs)
            if not in_args.out_format:
                seqbuddy.out_format = "gb"
            _print_recs(seqbuddy)
        except AttributeError as e:
            _raise_error(e, "annotate")
        _exit("annotate")

    # Average length of sequences
    if in_args.ave_seq_length:
        clean = False if not in_args.ave_seq_length[0] or in_args.ave_seq_length[0] != "clean" else True
        br._stdout("%s\n" % round(ave_seq_length(seqbuddy, clean), 2))
        _exit("ave_seq_length")

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

        if seqbuddy.alpha != "protein":
            _raise_error(TypeError("The input sequence needs to be protein, not nucleotide"), "back_translate")

        _print_recs(back_translate(seqbuddy, mode=mode, species=species, r_seed=in_args.random_seed))
        _exit("back_translate")

    # BL2SEQ
    if in_args.bl2seq:
        if len(find_repeats(seqbuddy).repeat_ids):
            br._stderr("Warning: There are records with duplicate ids which will be renamed.\n", quiet=in_args.quiet)
        try:
            output_dict = bl2seq(seqbuddy)
            br._stdout("#query\tsubject\t%_ident\tlength\tevalue\tbit_score\n")
            ids_already_seen = []
            for query_id, query_values in output_dict.items():
                ids_already_seen.append(query_id)
                query_values = [(key, value) for key, value in query_values.items()]
                query_values = sorted(query_values, key=lambda l: l[0])
                for subj_id, subj_values in query_values:
                    if subj_id in ids_already_seen:
                        continue

                    ident, length, evalue, bit_score = subj_values
                    br._stdout("%s\t%s\t%s\t%s\t%s\t%s\n" % (query_id, subj_id, ident, length, evalue, bit_score))
        except RuntimeError as e:
            _raise_error(e, "bl2seq", r"not present in \$PATH or working directory")
        _exit("bl2seq")

    # BLAST
    if in_args.blast:
        args = in_args.blast[0]
        params = re.sub(r"\[(.*)\]", "\1", args[1]) if len(args) > 1 else None
        blast_res = None
        try:
            try:
                blast_query = SeqBuddy(args[0])
                blast_res = blast(seqbuddy, blast_query, quiet=in_args.quiet, blast_args=params)
            except br.GuessError:
                blast_res = blast(seqbuddy, args[0], quiet=in_args.quiet, blast_args=params)
            except ValueError as e:
                _raise_error(e, "blast", ["num_threads expects an integer.", "evalue expects a number",
                                          "Trying to compare protein to nucleotide"])

            if len(blast_res) > 0:
                _print_recs(blast_res)
            else:
                br._stdout("No significant matches found\n")

        except (RuntimeError, SystemError) as e:
            _raise_error(e, "blast")
        _exit("blast")

    # Clean Seq
    if in_args.clean_seq:
        args = in_args.clean_seq[0]
        ambig = True
        rep_char = "N"
        lower_args = [str(x).lower() for x in args]
        if "strict" in lower_args:
            ambig = False
            del args[lower_args.index("strict")]
        if args and args[0]:
            rep_char = args[0][0]
        _print_recs(clean_seq(seqbuddy, ambiguous=ambig, rep_char=rep_char))
        _exit("clean_seq")

    # Complement
    if in_args.complement:
        try:
            _print_recs(complement(seqbuddy))
        except TypeError as e:
            _raise_error(e, "complement", "Nucleic acid sequence required, not protein.")
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
            if in_args.count_codons[0] and str(in_args.count_codons[0].lower()) in "concatenate":
                seqbuddy = concat_seqs(seqbuddy)
            codon_table = count_codons(seqbuddy)[1]
            for sequence_id in codon_table:
                br._stdout('#### {0} ####\n'.format(sequence_id))
                br._stdout('Codon\tAA\tNum\tPercent\n')
                for codon in codon_table[sequence_id]:
                    data = codon_table[sequence_id][codon]
                    br._stdout('{0}\t{1}\t{2}\t{3}\n'.format(codon, data[0], data[1], data[2]))
                br._stdout('\n')
        except TypeError as e:
            _raise_error(e, "count_codons", "Nucleic acid sequence required, not protein or other.")
        _exit("count_codons")

    # Count residues
    if in_args.count_residues:
        seqbuddy = replace_subsequence(seqbuddy, "[-.]", "")
        if in_args.count_residues[0] and str(in_args.count_residues[0].lower()) in "concatenate":
            seqbuddy = concat_seqs(seqbuddy)
        count_residues(seqbuddy)
        for rec in seqbuddy.records:
            br._stdout("%s\n" % str(rec.id))
            for residue, counts in rec.buddy_data["res_count"].items():
                try:
                    br._stdout("{0}:\t{1}\t{2} %\n".format(residue, counts[0], round(counts[1] * 100, 2)))
                except TypeError:
                    br._stdout("{0}:\t{1}\n".format(residue, counts))
            br._stdout("\n")
        _exit("count_residues")

    # Delete features
    if in_args.delete_features:
        patterns = br.clean_regex(in_args.delete_features, in_args.quiet)
        for next_pattern in patterns:
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
        try:  # Check to see if the last argument is an integer, which will set number of columns
            if len(in_args.delete_records) == 1:
                columns = 1
            else:
                columns = int(in_args.delete_records[-1])
                del in_args.delete_records[-1]
        except ValueError:
            columns = 1

        # Toggle on/off full search in record (including metadata)
        if "full" in in_args.delete_records:
            description = True
            del in_args.delete_records[in_args.delete_records.index("full")]
        else:
            description = False
        search_terms = []
        for arg in in_args.delete_records:
            if os.path.isfile(arg):
                with open(arg, "r", encoding="utf-8") as ifile:
                    for line in ifile:
                        search_terms.append(line.strip())
            else:
                search_terms.append(arg)

        search_terms = br.clean_regex(search_terms, in_args.quiet)
        if not search_terms:  # If all regular expression are malformed, exit out gracefully
            _print_recs(seqbuddy)
            _exit("delete_records")

        deleted_seqs = []
        for next_pattern in search_terms:
            deleted_seqs += pull_recs(make_copy(seqbuddy), next_pattern, description).records

        seqbuddy = delete_records(seqbuddy, search_terms, description)

        if len(deleted_seqs) > 0 and not in_args.quiet:
            counter = 1
            br._stderr("# ####################### Deleted records ######################## #\n", in_args.quiet)
            for seq in deleted_seqs:
                br._stderr(seq.id, in_args.quiet)
                if counter % columns == 0 or counter == len(deleted_seqs):
                    br._stderr("\n", in_args.quiet)
                else:
                    br._stderr("\t", in_args.quiet)
                counter += 1
            br._stderr("# ################################################################ #\n", in_args.quiet)

        if len(deleted_seqs) == 0:
            br._stderr("# ################################################################ #\n", in_args.quiet)
            br._stderr("# No sequence identifiers match %s\n" % ", ".join(search_terms), in_args.quiet)
            br._stderr("# ################################################################ #\n", in_args.quiet)
        _print_recs(seqbuddy)
        _exit("delete_records")

    # Delete records with feature
    if in_args.delete_recs_with_feature:
        search_terms = []
        for arg in in_args.delete_recs_with_feature:
            if os.path.isfile(arg):
                with open(arg, "r", encoding="utf-8") as ifile:
                    for line in ifile:
                        search_terms.append(line.strip())
            else:
                search_terms.append(arg)

        search_terms = br.clean_regex(search_terms, in_args.quiet)
        if search_terms:
            seqbuddy = delete_recs_with_feature(seqbuddy, search_terms)
        _print_recs(seqbuddy)
        _exit("delete_recs_with_feature")

    # Delete repeats
    if in_args.delete_repeats:
        dlt_repeats = in_args.delete_repeats[0]
        columns = 1
        scope = "all"
        if dlt_repeats:
            for arg in dlt_repeats:
                try:
                    columns = int(arg)
                except ValueError:
                    for scope_option in ["all", "ids", "seqs"]:
                        scope = scope_option if scope_option.startswith(arg) else scope

        find_repeats(seqbuddy)

        if len(seqbuddy.repeat_ids) > 0:
            br._stderr("# ################################################################ #\n", in_args.quiet)

        if len(seqbuddy.repeat_ids) > 0 and scope in ["all", "ids"]:
            br._stderr("# Records with duplicate ids deleted\n", in_args.quiet)
            counter = 1
            for seq in seqbuddy.repeat_ids:
                br._stderr(seq, in_args.quiet)
                if counter % columns == 0 or counter == len(seqbuddy.repeat_ids):
                    br._stderr("\n", in_args.quiet)
                else:
                    br._stderr("\t", in_args.quiet)
                counter += 1
            br._stderr("\n", in_args.quiet)
            seqbuddy = delete_repeats(seqbuddy, 'ids')
            find_repeats(seqbuddy)

        rep_seq_ids = []
        for seq in seqbuddy.repeat_seqs:
            rep_seq_ids.append([])
            for rep_seq_id in seqbuddy.repeat_seqs[seq]:
                rep_seq_ids[-1].append(rep_seq_id)

        if len(rep_seq_ids) > 0 and scope in ["all", "seqs"]:
            br._stderr("# Records with duplicate sequence deleted\n", in_args.quiet)
            counter = 1
            for repeat_seqs in rep_seq_ids:
                br._stderr("[", in_args.quiet)
                for _indx, rep_seq in enumerate(repeat_seqs):
                    br._stderr(rep_seq, in_args.quiet)
                    if _indx + 1 != len(repeat_seqs):
                        br._stderr(", ", in_args.quiet)
                br._stderr("]", in_args.quiet)
                if counter % columns == 0 or counter == len(rep_seq_ids):
                    br._stderr("\n", in_args.quiet)
                else:
                    br._stderr(", ", in_args.quiet)
                counter += 1

        if len(rep_seq_ids) > 0:
            br._stderr("# ################################################################ #\n\n", in_args.quiet)
        else:
            br._stderr("No duplicate records found\n", in_args.quiet)
        _print_recs(delete_repeats(seqbuddy, 'seqs'))
        _exit("delete_repeats")

    # Delete sequences below threshold
    if in_args.delete_small:
        _print_recs(delete_small(seqbuddy, in_args.delete_small))
        _exit("delete_small")

    # Delete taxa
    if in_args.delete_taxa:
        seqbuddy = delete_taxa(seqbuddy, in_args.delete_taxa[0])
        _print_recs(seqbuddy)
        _exit("delete_taxa")

    # degenerate_sequence
    if in_args.degenerate_sequence:
        # if no argument provided will use table 1 first reading frame as default(set above)
        degen_args = -1 if not in_args.degenerate_sequence[0] else in_args.degenerate_sequence[0]
        try:
            degenerate_sequence(seqbuddy, degen_args)
        except KeyError as e:
            br._stderr("Valid dictionaries:\n"
                       "Arg\tDictionary\n"
                       "1\tStandard Genetic Code\n"
                       "2\tVertebrate Mitochondrial Code\n"
                       "3\tYeast Mitochondrial Code\n"
                       "4\tMold / Protozoan / Coelenterate Mitochondrial Code & Mycoplasma / Spiroplasma Code\n"
                       "5\tInvertebrate Mitochondrial Code\n"
                       "6\tCiliate / Dasycladacean / Hexamita Nuclear Code\n"
                       "7\tEchinoderm and Flatworm Mitochondrial Code\n"
                       "8\tEuplotid Nuclear Code\n"
                       "9\tBacterial / Archaeal / Plant Plastid Code\n"
                       "10\tAlternative Yeast Nuclear Code\n"
                       "11\tAscidian Mitochondrial Code\n"
                       "12\tAlternative Flatworm Mitochondrial Code\n\n", in_args.quiet)
            _raise_error(e, "degenerate_sequence", "Could not locate codon dictionary")
        except TypeError as e:
            _raise_error(e, "degenerate_sequence", "Nucleic acid sequence required, not protein")
        _print_recs(seqbuddy)
        _exit('degenerate_sequence')

    # Extact features
    if in_args.extract_feature_sequences:
        patterns = br.clean_regex(in_args.extract_feature_sequences[0], in_args.quiet)
        if patterns:
            seqbuddy = extract_feature_sequences(seqbuddy, patterns)
        _print_recs(seqbuddy)
        _exit("extract_feature_sequences")

    # Extract regions
    if in_args.extract_regions:
        try:
            args = ",".join(in_args.extract_regions[0])
            seqbuddy = extract_regions(seqbuddy, args)
            _print_recs(seqbuddy)
        except ValueError as e:
            br._stderr("""
Extraction regions can be specified in the following ways:
    - Singlets: "2,5,9,-5"
    - Ranges: "40:75,89:100,432:-45"
    - mth of nth: "1/5,3/5"

For examples, please see:
https://github.com/biologyguy/BuddySuite/wiki/SB-Extract-regions

""")
            _raise_error(e, "extract_regions", "Unable to decode the positions string")
        _exit("extract_regions")

    # Find CpG
    if in_args.find_CpG:
        try:
            find_cpg(seqbuddy)
            islands = False
            for rec in seqbuddy.records:
                if rec.buddy_data["cpgs"]:
                    islands = True
                    break

            if islands:
                br._stderr('########### Islands identified ###########\n', in_args.quiet)
                for _indx, rec in enumerate(seqbuddy.records):
                    if rec.buddy_data["cpgs"]:
                        value = ["%s-%s" % (x[0], x[1]) for x in rec.buddy_data["cpgs"]]
                        br._stderr("{0}: {1}".format(rec.id, ", ".join(value)), in_args.quiet)
                        if _indx + 1 != len(seqbuddy.records):
                            br._stderr("\n", in_args.quiet)

                br._stderr('\n##########################################\n\n', in_args.quiet)
            else:
                br._stderr("# No Islands identified\n\n", in_args.quiet)
            _print_recs(seqbuddy)
            _exit("find_CpG")

        except TypeError as e:
            _raise_error(e, "find_CpG", "DNA sequence required, not protein or RNA.")

    # Find orfs
    if in_args.find_orfs:
        min_size = None
        rev_comp = True
        for arg in in_args.find_orfs[0]:
            if arg.lower() == "false":
                rev_comp = False
                continue
            try:
                min_size = int(arg)
            except ValueError:
                pass

        try:
            find_orfs(seqbuddy, min_size=min_size, rev_comp=rev_comp)
            for rec in seqbuddy.records:
                pos_indices = rec.buddy_data['find_orfs']['+']
                neg_indices = rec.buddy_data['find_orfs']['-']
                br._stderr("# {0}\n".format(rec.id), in_args.quiet)
                if len(pos_indices) <= 0:
                    br._stderr("(+) ORFs: None\n", in_args.quiet)
                else:
                    orfs = ", ".join(["%s:%s" % (x[0], x[1]) for x in pos_indices if len(x) > 0])
                    br._stderr("(+) ORFs: %s\n" % orfs, in_args.quiet)
                if rev_comp:
                    if len(neg_indices) <= 0:
                        br._stderr("(-) ORFs: None\n", in_args.quiet)
                    else:
                        orfs = ", ".join(["%s:%s" % (x[1], x[0]) for x in neg_indices if len(x) > 0])
                        br._stderr("(-) ORFs: %s\n" % orfs, in_args.quiet)
            br._stderr("\n", in_args.quiet)
            _print_recs(seqbuddy)
        except TypeError as err:
            _raise_error(err, "find_orfs", "Nucleic acid sequence required, not protein.")
        except ValueError as err:
            _raise_error(err, "find_orfs", "Open reading frames cannot be smaller than 6 residues.")
        _exit("find_orfs")

    # Find pattern
    if in_args.find_pattern:
        ambig = True if 'ambig' in in_args.find_pattern else False
        if ambig:
            del in_args.find_pattern[in_args.find_pattern.index("ambig")]

        patterns = br.clean_regex(in_args.find_pattern, in_args.quiet)
        if patterns:
            find_pattern(seqbuddy, *patterns, ambig=ambig)
        for pattern in patterns:
            num_matches = 0
            for rec in seqbuddy.records:
                indices = rec.buddy_data['find_patterns'][pattern]
                num_matches += len(indices)
            br._stderr("#### {0} matches found across {1} sequences for "
                       "pattern '{2}' ####\n".format(num_matches, len(seqbuddy), pattern), in_args.quiet)
            for rec in seqbuddy.records:
                indices = rec.buddy_data['find_patterns'][pattern]
                if not len(indices):
                    br._stderr("{0}: None\n".format(rec.id), in_args.quiet)
                else:
                    matches = ", ".join(["%s-%s" % (start + 1, end) for start, end in indices])
                    br._stderr("{0}: {1}\n".format(rec.id, matches), in_args.quiet)
                    num_matches += len(indices)

            br._stderr("\n", in_args.quiet)
        _print_recs(seqbuddy)
        _exit("find_pattern")

    # Find repeat sequences or ids
    if in_args.find_repeats:
        columns = 1 if not in_args.find_repeats[0] else in_args.find_repeats[0]
        find_repeats(seqbuddy)

        if len(seqbuddy.repeat_ids) > 0:
            br._stdout("#### Records with duplicate IDs: ####\n")
            counter = 1
            for next_id in seqbuddy.repeat_ids:
                if counter % columns == 0 or counter == len(seqbuddy.repeat_ids):
                    br._stdout("%s\n" % next_id)
                else:
                    br._stdout("%s\t" % next_id)
                counter += 1

            br._stdout("\n")

        else:
            br._stdout("#### No records with duplicate IDs ####\n\n")

        if len(seqbuddy.repeat_seqs) > 0:
            br._stdout("#### Records with duplicate sequences: ####\n")
            counter = 1
            for next_id in seqbuddy.repeat_seqs:
                br._stdout("[")
                for _indx, seq_id in enumerate(seqbuddy.repeat_seqs[next_id]):
                    br._stdout("%s" % seq_id)
                    if _indx + 1 != len(seqbuddy.repeat_seqs[next_id]):
                        br._stdout(", ")
                br._stdout("]")

                if counter % columns != 0 and counter != len(seqbuddy.repeat_seqs):
                    br._stdout(", ")
                else:
                    br._stdout("\n")

                counter += 1
        else:
            br._stdout("#### No records with duplicate sequences ####\n")
        _exit("find_repeats")

    # Find restriction sites
    if in_args.find_restriction_sites:
        if not in_args.out_format:
            seqbuddy.out_format = "gb"

        _enzymes, order, min_cuts, max_cuts, topology = _prepare_restriction_sites(in_args.find_restriction_sites[0])

        clean_seq(seqbuddy)
        try:
            find_restriction_sites(seqbuddy, tuple(_enzymes), min_cuts, max_cuts, topology, quiet=in_args.quiet)
        except TypeError as e:
            _raise_error(e, "find_restriction_sites", check_string='Unable to identify restriction')

        br._stderr('# ### Restriction Sites (indexed at cut-site) ### #\n', in_args.quiet)
        for tup in seqbuddy.restriction_sites:
            br._stderr("{0}\n".format(tup[0]), in_args.quiet)
            restriction_list = tup[1]
            restriction_list = [[key, value] for key, value in restriction_list.items()]
            restriction_list = sorted(restriction_list, key=lambda l: str(l[0])) if order == 'alpha' else \
                sorted(restriction_list, key=lambda l: l[1])

            for _enzyme in restriction_list:
                cut_sites = [str(x) for x in _enzyme[1]]
                br._stderr("{0}\t{1}\n".format(_enzyme[0], ", ".join(cut_sites)), in_args.quiet)
            if tup != seqbuddy.restriction_sites[-1]:
                br._stderr("\n", in_args.quiet)
        br._stderr("# ############################################### #\n\n", in_args.quiet)
        _print_recs(seqbuddy)
        _exit("find_restriction_sites")

    # Group sequences by prefix. I might want to delete this in favour of group_by_regex... Keep them both for now.
    if in_args.group_by_prefix:
        in_args.in_place = True
        check_quiet = in_args.quiet
        in_args.quiet = True  # toggle 'quiet' on so in_place _print_recs() doesn't spam with print messages

        args = in_args.group_by_prefix[0]
        out_dir = os.getcwd()
        num_chars = 0
        split_patterns = []
        for arg in args:
            try:
                num_chars = int(arg)
            except ValueError:
                if os.path.isdir(arg):
                    out_dir = os.path.abspath(arg)
                else:
                    if arg == '':
                        split_patterns.append('-')
                    else:
                        split_patterns.append(arg)

        sp = ["-"] if not split_patterns and not num_chars else split_patterns
        sp = br.clean_regex(sp, check_quiet)
        if not sp:
            in_args.quiet = check_quiet
            _raise_error(ValueError("Split pattern(s) malformed. No files created."), "group_by_prefix")

        taxa_groups = make_groups(seqbuddy, split_patterns=sp, num_chars=num_chars)
        if "".join(split_patterns) != "" and len(taxa_groups) == len(seqbuddy):
            taxa_groups = make_groups(seqbuddy, num_chars=5)

        for next_seqbuddy in taxa_groups:
            in_args.sequence[0] = "%s%s%s.%s" % (out_dir, os.path.sep, next_seqbuddy.identifier,
                                                 br.format_to_extension[next_seqbuddy.out_format])
            br._stderr("New file: %s\n" % in_args.sequence[0], check_quiet)
            open(in_args.sequence[0], "w", encoding="utf-8").close()
            _print_recs(next_seqbuddy)
        _exit("group_by_prefix")

    # Group sequences by regex. This is really flexible.
    if in_args.group_by_regex:
        in_args.in_place = True
        check_quiet = in_args.quiet
        in_args.quiet = True  # toggle 'quiet' on so in_place _print_recs() doesn't spam with print messages

        args = in_args.group_by_regex[0]
        out_dir = os.getcwd()
        regexes = []
        for arg in args:
            if os.path.isdir(arg):
                out_dir = os.path.abspath(arg)
            else:
                regexes.append(arg)

        regexes = br.clean_regex(regexes, check_quiet)
        if not regexes:
            in_args.quiet = False
            _raise_error(ValueError("You must provide at least one valid regular expression."), "group_by_regex")

        taxa_groups = make_groups(seqbuddy, regex=regexes)

        for next_seqbuddy in taxa_groups:
            in_args.sequence[0] = "%s%s%s.%s" % (out_dir, os.path.sep, next_seqbuddy.identifier,
                                                 br.format_to_extension[next_seqbuddy.out_format])
            br._stderr("New file: %s\n" % in_args.sequence[0], check_quiet)
            open(in_args.sequence[0], "w", encoding="utf-8").close()
            _print_recs(next_seqbuddy)
        _exit("group_by_regex")

    # Guess alphabet
    if in_args.guess_alphabet:
        for seq_set in in_args.sequence:
            try:
                seqbuddy = SeqBuddy(seq_set)
            except Exception:  # This should NOT be made more specific. If it throws errors, it's unknown.
                seqbuddy = SeqBuddy("", in_format="raw")

            if str(type(seq_set)) != "<class '_io.TextIOWrapper'>":
                path, seq_set = os.path.split(seq_set)
                br._stdout("%s\t-->\t" % seq_set)
            else:
                br._stdout("PIPE\t-->\t")

            if seqbuddy.alpha == "protein":
                br._stdout("prot\n")
            elif seqbuddy.alpha == "DNA":
                br._stdout("dna\n")
            elif seqbuddy.alpha == "RNA":
                br._stdout("rna\n")
            else:
                br._stdout("Undetermined\n")
        _exit("guess_alphabet")

    # Guess format
    if in_args.guess_format:
        for seq_set in in_args.sequence:
            if str(type(seq_set)) == "<class '_io.TextIOWrapper'>":
                _format = br.guess_format(seq_set)
                if _format:
                    br._stdout("PIPE\t-->\t%s\n" % _format)
                else:
                    br._stdout("PIPE\t-->\tUnknown\n")
            else:
                try:
                    file_format = br.guess_format(seq_set)
                except Exception:  # This should NOT be made more specific. If it throws errors, it's unknown.
                    file_format = None

                path, seq_set = os.path.split(seq_set)
                if not file_format:
                    br._stdout("%s\t-->\tUnknown\n" % seq_set)
                else:
                    br._stdout("%s\t-->\t%s\n" % (seq_set, file_format))
        _exit("guess_format")

    # Hash sequence ids
    if in_args.hash_ids:
        if in_args.hash_ids[0] == 0:
            hash_length = 0
        elif not in_args.hash_ids[0]:
            hash_length = 10
        else:
            hash_length = in_args.hash_ids[0]

        if hash_length < 1:
            br._stderr("Warning: The hash_length parameter was passed in with the value %s. This is not a positive "
                       "integer, so the hash length as been set to 10.\n\n" % hash_length, quiet=in_args.quiet)
            hash_length = 10

        if 32 ** hash_length <= len(seqbuddy) * 2:
            holder = ceil(log(len(seqbuddy) * 2, 32))
            br._stderr("Warning: The hash_length parameter was passed in with the value %s. "
                       "This is too small to properly cover all sequences, "
                       "so it has been increased to %s.\n\n" % (hash_length, holder), in_args.quiet)
            hash_length = holder

        hash_ids(seqbuddy, hash_length=hash_length, r_seed=in_args.random_seed)

        hash_table = "# Hash table\n"
        hash_table += seqbuddy.print_hashmap()
        hash_table += "\n"

        br._stderr(hash_table, in_args.quiet)
        _print_recs(seqbuddy)
        _exit("hash_ids")

    # Head
    if in_args.head:
        head_size = 1 if in_args.head[0] is None else in_args.head[0]
        head(seqbuddy, head_size)
        _print_recs(seqbuddy)
        _exit("head")

    # Insert Seq
    if in_args.insert_seq:
        args = in_args.insert_seq[0]
        if len(args) < 2:
            _raise_error(AttributeError("The insert_seq tool requires at least two arguments (sequence and position)"),
                         "insert_seq")

        sequence = ""
        location = 0
        try:
            int(args[0])
            _raise_error(AttributeError("The first argment must be your insert sequence, not location."), "insert_seq")
        except ValueError:
            sequence = args[0]

        try:
            location = int(args[1])
        except ValueError:
            _raise_error(AttributeError("The second argment must be location, not insert sequence or regex."),
                         "insert_seq")

        regex = [] if len(args) < 2 else args[2:]

        _print_recs(insert_sequence(seqbuddy, sequence, location, regex))
        _exit("insert_seq")

    # In silico digest
    if in_args.in_silico_digest:
        _enzymes, order, min_cuts, max_cuts, topology = _prepare_restriction_sites(in_args.in_silico_digest[0])

        if not in_args.in_silico_digest[0] or _enzymes == ["commercial"]:
            br._stderr("Error: Please provide a list of enzymes you wish to cut your sequences with.\n")
            _exit("in_silico_digest")

        clean_seq(seqbuddy)
        try:
            seqbuddy = in_silico_digest(seqbuddy, tuple(_enzymes), topology=topology)
        except TypeError as e:
            _raise_error(e, "in_silico_digest", check_string='Unable to identify restriction')

        _print_recs(seqbuddy)
        _exit("in_silico_digest")

    # Isoelectric Point
    if in_args.isoelectric_point:
        if seqbuddy.alpha != "protein":
            br._stderr("Nucleic acid sequences detected, converting to protein.\n\n")
            seqbuddy = translate_cds(seqbuddy, quiet=True)

        isoelectric_point(seqbuddy)
        br._stderr("ID\tpI\n")
        for rec in seqbuddy.records:
            br._stdout("{0}\t{1}\n".format(rec.id, round(rec.features[-1].qualifiers["value"], 3)))
        _exit("isoelectric_point")

    # keep taxon
    if in_args.keep_taxa:
        seqbuddy = keep_taxa(seqbuddy, in_args.keep_taxa[0])
        _print_recs(seqbuddy)
        _exit("keep_taxa")

    # List features
    if in_args.list_features:
        for rec in seqbuddy.records:
            br._stdout('#### {0} ####\n'.format(rec.id))
            if len(rec.features) > 0:
                for feat in rec.features:
                    br._stdout('{0}\n'.format(feat.type))
                    if isinstance(feat.location, CompoundLocation):
                        br._stdout('\tLocations:\n')
                        for part in feat.location.parts:
                            br._stdout('\t\t{0}-{1}\n'.format(part.start, part.end))
                    elif isinstance(feat.location, FeatureLocation):
                        br._stdout('\tLocation:\n')
                        br._stdout('\t\t{0}-{1}\n'.format(feat.location.start, feat.location.end))
                    if feat.strand == 1:
                        br._stdout('\tStrand: Sense(+)\n')
                    elif feat.strand == -1:
                        br._stdout('\tStrand: Antisense(-)\n')
                    if str(feat.id) != '<unknown id>':
                        br._stdout('\tID: {0}\n'.format(feat.id))
                    if len(feat.qualifiers) > 0:
                        br._stdout('\tQualifiers:\n')
                        qualifs = OrderedDict(sorted(feat.qualifiers.items()))
                        for key, qual in qualifs.items():
                            br._stdout('\t\t{0}: {1}\n'.format(key, str(qual).strip("[]'")))
            else:
                br._stdout('None\n')
            br._stdout("\n")
        _exit("list_features")

    # List identifiers
    if in_args.list_ids:
        columns = 1 if not in_args.list_ids[0] else abs(in_args.list_ids[0])
        for indx, rec in enumerate(seqbuddy.records):
            if (indx + 1) % columns == 0 or (indx + 1) == len(seqbuddy.records):
                br._stdout("%s\n" % rec.id, )
            else:
                br._stdout("%s\t" % rec.id)
        _exit("list_ids")

    # Lowercase
    if in_args.lowercase:
        _print_recs(lowercase(seqbuddy))
        _exit("lowercase")

    # Make unique IDs
    if in_args.make_ids_unique:
        args = in_args.make_ids_unique[0]
        padding = 0
        sep = ""
        if args:
            for arg in args:
                try:
                    padding = int(arg)
                except ValueError:
                    sep = arg
        _print_recs(make_ids_unique(seqbuddy, sep=sep.strip(), padding=padding))
        _exit("make_ids_unique")

    # Map features from cDNA over to protein
    if in_args.map_features_nucl2prot:
        if len(in_args.sequence) < 2:
            _raise_error(ValueError("You must provide one DNA file and one protein file"), "map_features_nucl2prot")
        file1, file2 = in_args.sequence[:2]
        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)
        if file1.alpha == file2.alpha or (file1.alpha != "protein" and file2.alpha != "protein"):
            _raise_error(ValueError("You must provide one DNA file and one protein file"), "map_features_nucl2prot")
        if file1.alpha == "protein":
            prot = file1
            dna = file2
        else:
            in_args.sequence[0] = in_args.sequence[1]  # in case the -i flag is thrown
            prot = file2
            dna = file1
        try:
            seqbuddy = map_features_nucl2prot(dna, prot, quiet=in_args.quiet)
        except RuntimeError as e:
            _raise_error(e, "map_features_nucl2prot", "There are repeat IDs in self.records")

        seqbuddy.out_format = in_args.out_format if in_args.out_format else "gb"
        _print_recs(seqbuddy)
        _exit("map_features_nucl2prot")

    # Map features from protein over to cDNA
    if in_args.map_features_prot2nucl:
        if len(in_args.sequence) < 2:
            _raise_error(ValueError("You must provide one DNA file and one protein file"), "map_features_nucl2prot")
        file1, file2 = in_args.sequence[:2]
        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)
        if file1.alpha == file2.alpha or (file1.alpha != "protein" and file2.alpha != "protein"):
            _raise_error(ValueError("You must provide one DNA file and one protein file"), "map_features_nucl2prot")
        if file1.alpha != "protein":
            dna = file1
            prot = file2
        else:
            in_args.sequence[0] = in_args.sequence[1]  # in case the -i flag is thrown
            dna = file2
            prot = file1
        try:
            seqbuddy = map_features_prot2nucl(prot, dna, quiet=in_args.quiet)
        except RuntimeError as e:
            _raise_error(e, "map_features_nucl2prot", "There are repeat IDs in self.records")

        seqbuddy.out_format = in_args.out_format if in_args.out_format else "gb"
        _print_recs(seqbuddy)
        _exit("map_features_prot2nucl")

    # Max record length
    if in_args.max_recs:
        num_returned = 1 if not in_args.max_recs[0] else in_args.max_recs[0]
        _print_recs(max_records(seqbuddy, num_returned))
        _exit("max_recs")

    # Merge together multiple files into a single file
    if in_args.merge:
        seqbuddy_objs = [SeqBuddy(x) for x in in_args.sequence]
        try:
            _print_recs(merge(*seqbuddy_objs))
        except RuntimeError as e:
            _raise_error(e, "merge")
        _exit("merge")

    # Min record length
    if in_args.min_recs:
        num_returned = 1 if not in_args.min_recs[0] else in_args.min_recs[0]
        _print_recs(min_records(seqbuddy, num_returned))
        _exit("min_recs")

    # Molecular Weight
    if in_args.molecular_weight:
        try:
            molecular_weight(seqbuddy)
            mws = seqbuddy.molecular_weights
            if seqbuddy.alpha == "DNA":
                br._stderr("ID\tssDNA\tdsDNA\n")
            elif seqbuddy.alpha == "RNA":
                br._stderr("ID\tssRNA\n")
            else:
                br._stderr("ID\tProtein\n")
            for indx, value in enumerate(mws['ids']):
                if len(mws['masses_ds']) != 0:
                    print("{0}\t{1}\t{2}".format(value, mws['masses_ss'][indx], mws['masses_ds'][indx]))
                else:
                    print("{0}\t{1}".format(value, mws['masses_ss'][indx]))
        except KeyError as e:
            _raise_error(e, "molecular_weight")
        _exit("molecular_weight")

    # Number of sequences
    if in_args.num_seqs:
        br._stdout("%s\n" % num_seqs(seqbuddy))
        _exit("num_seqs")

    # Order sequence features alphabetically
    if in_args.order_features_alphabetically:
        ofa = in_args.order_features_alphabetically
        if ofa[0] and type(ofa[0]) == str:
            reverse = True if "reverse".startswith(ofa[0].lower()) else False
        else:
            reverse = False
        _print_recs(order_features_alphabetically(seqbuddy, reverse))
        _exit("order_features_alphabetically")

    # Order sequence features by their position in the sequence
    if in_args.order_features_by_position:
        ofp = in_args.order_features_by_position
        if ofp[0] and type(ofp[0]) == str:
            reverse = True if "reverse".startswith(ofp[0].lower()) else False
        else:
            reverse = False
        _print_recs(order_features_by_position(seqbuddy, reverse))
        _exit("order_features_by_position")

    # Order ids
    if in_args.order_ids:
        if in_args.order_ids[0] and type(in_args.order_ids[0]) == str:
            reverse = True if "reverse".startswith(in_args.order_ids[0].lower()) else False
        else:
            reverse = False
        _print_recs(order_ids(seqbuddy, reverse=reverse))
        _exit("order_ids")

    # Order ids randomly
    if in_args.order_ids_randomly:
        _print_recs(order_ids_randomly(seqbuddy, r_seed=in_args.random_seed))
        _exit("order_ids_randomly")

    # Order records by length
    if in_args.order_recs_by_len:
        rev = True if in_args.order_recs_by_len[0] and "reverse".startswith(in_args.order_recs_by_len[0].lower()) \
            else False
        _print_recs(order_recs_by_len(seqbuddy, rev))
        _exit("order_seqs_by_size")

    # Prepend organism
    if in_args.prepend_organism:
        try:
            length = in_args.prepend_organism[0]
            seqbuddy = prepend_organism(seqbuddy) if length is None else prepend_organism(seqbuddy, length)
        except ValueError as err:
            _raise_error(err, "prepend_organism", check_string=["Multiple species", "Prefix length must be"])

        prefix_map = sorted(["%s: %s" % (prefix, organism) for prefix, organism in seqbuddy.prefix_map.items()])

        if prefix_map:
            br._stderr("# ######################## Prefix Mapping ######################## #\n", in_args.quiet)
            br._stderr("%s\n" % "\n".join(prefix_map), in_args.quiet)
            br._stderr("# ################################################################ #\n\n", in_args.quiet)

        else:
            br._stderr("# ######################## Prefix Mapping ###################### #\n", in_args.quiet)
            br._stderr("# No organism information was identified in the supplied records #\n", in_args.quiet)
            br._stderr("# ############################################################## #\n\n", in_args.quiet)

        _print_recs(seqbuddy)
        _exit("prepend_organism")

    # Prosite Scan
    if in_args.prosite_scan:
        try:
            ps_scan = PrositeScan(seqbuddy, quiet=in_args.quiet, r_seed=in_args.random_seed)
            seqbuddy = ps_scan.run()
        except urllib.error.URLError as err:
            if "Errno 8" in str(err):
                print("Unable to contact EBI, are you connected to the internet?")
                _raise_error(err, "prosite_scan")
            else:
                raise err
        if not in_args.out_format:
            seqbuddy.out_format = "gb"
        else:
            seqbuddy.out_format = in_args.out_format
        _print_recs(seqbuddy)
        _exit("prosite")

    # Pull random records
    if in_args.pull_random_record:
        count = 1 if not in_args.pull_random_record[0] else in_args.pull_random_record[0]
        _print_recs(pull_random_recs(seqbuddy, count=count, r_seed=in_args.random_seed))
        _exit("pull_random_record")

    # Pull record ends
    if in_args.pull_record_ends:
        _print_recs(pull_record_ends(seqbuddy, in_args.pull_record_ends))
        _exit("pull_record_ends")

    # Pull records
    if in_args.pull_records:
        if "full" in in_args.pull_records:
            description = True
            del in_args.pull_records[in_args.pull_records.index("full")]
        else:
            description = False

        search_terms = []
        regex_file = False
        for arg in in_args.pull_records:
            if os.path.isfile(arg):
                regex_file = True
                with open(arg, "r", encoding="utf-8") as ifile:
                    for line in ifile:
                        line = line.strip()
                        if line:
                            search_terms.append(line.strip())
            else:
                search_terms.append(arg)

        search_terms = br.clean_regex(search_terms, in_args.quiet)
        if search_terms:
            seqbuddy = pull_recs(seqbuddy, search_terms, description)
        elif regex_file:  # If an empty regex file is passed in then return no records
            seqbuddy.records = []
        _print_recs(seqbuddy)
        _exit("pull_records")

    # Pull records with feature
    if in_args.pull_records_with_feature:
        search_terms = []
        for arg in in_args.pull_records_with_feature:
            if os.path.isfile(arg):
                with open(arg, "r", encoding="utf-8") as ifile:
                    for line in ifile:
                        search_terms.append(line.strip())
            else:
                search_terms.append(arg)

        search_terms = br.clean_regex(search_terms, in_args.quiet)
        if search_terms:
            seqbuddy = pull_recs_with_feature(seqbuddy, search_terms)
        _print_recs(seqbuddy)
        _exit("pull_records_with_feature")

    # Purge
    if in_args.purge:
        purge(seqbuddy, in_args.purge)
        br._stderr("### Deleted record mapping ###\n", in_args.quiet)
        for indx1, rec in enumerate(seqbuddy.records):
            br._stderr("%s\n" % rec.id, in_args.quiet)
            if rec.buddy_data["purge_set"]:
                for indx2, del_seq_id in enumerate(rec.buddy_data["purge_set"]):
                    br._stderr(del_seq_id, in_args.quiet)
                    if indx2 + 1 != len(rec.buddy_data["purge_set"]):
                        br._stderr(", ", in_args.quiet)
            if indx1 + 1 != len(seqbuddy.records):
                br._stderr("\n\n", in_args.quiet)
            else:
                br._stderr("\n", in_args.quiet)
        br._stderr("##############################\n\n", in_args.quiet)
        _print_recs(seqbuddy)
        _exit("purge")

    # Renaming
    if in_args.rename_ids:
        args = in_args.rename_ids[0]
        if len(args) < 2:
            _raise_error(AttributeError("Please provide at least a query and a replacement string"), "rename_ids")

        query, replace = args[0:2]
        if not br.clean_regex(query, in_args.quiet):
            _raise_error(ValueError("Malformed regular expression."), "rename_ids")
        num = 0
        store = False

        if len(args) > 2:
            args = args[2:]
            if "store" in args:
                store = True
                del args[args.index("store")]

            try:
                num = num if not len(args) else int(args[0])
            except ValueError:
                _raise_error(ValueError("Max replacements argument must be an integer"), "rename_ids")
        try:
            _print_recs(rename(seqbuddy, query=query, replace=replace, num=num, store_old_id=store))
        except AttributeError as e:
            _raise_error(e, "rename_ids", "There are more replacement")
        _exit("rename_ids")

    # Replace sub-sequences
    if in_args.replace_subseq:
        args = in_args.replace_subseq[0]
        args = args[:2] if len(args) > 1 else args
        if not br.clean_regex(args[0], in_args.quiet):
            _raise_error(ValueError("Max replacements argument must be an integer"), "replace_subseq")
        _print_recs(replace_subsequence(seqbuddy, *args))
        _exit("replace_subseq")

    # Reverse complement
    if in_args.reverse_complement:
        try:
            _print_recs(reverse_complement(seqbuddy))
        except TypeError as e:
            _raise_error(e, "reverse_complement", "Nucleic acid sequences required.")
        _exit("reverse_complement")

    # Reverse Transcribe
    if in_args.reverse_transcribe:
        try:
            _print_recs(rna2dna(seqbuddy))
        except TypeError as e:
            _raise_error(e, "reverse_transcribe", "RNA sequence required, not")
        _exit("reverse_transcribe")

    # Screw formats
    if in_args.screw_formats:
        if in_args.screw_formats.lower() not in OUTPUT_FORMATS:
            _raise_error(IOError("Error: unknown format '%s'\n" % in_args.screw_formats), "screw_formats")

        seqbuddy.out_format = in_args.screw_formats
        if in_args.in_place:  # Need to change the file extension
            _path, ext = os.path.splitext(os.path.abspath(in_args.sequence[0]))
            _path = "%s.%s" % (_path, br.format_to_extension[seqbuddy.out_format])

            os.remove(in_args.sequence[0])
            in_args.sequence[0] = _path
            open(in_args.sequence[0], "w", encoding="utf-8").close()
        _print_recs(seqbuddy)
        _exit("screw_formats")

    # Shift reading frame
    if in_args.select_frame:
        try:
            _print_recs(select_frame(seqbuddy, in_args.select_frame))
        except TypeError as e:
            _raise_error(e, "reverse_complement", "Select frame requires nucleic acid, not protein.")
        _exit("select_frame")

    # Shuffle Seqs
    if in_args.shuffle_seqs:
        _print_recs(shuffle_seqs(seqbuddy, r_seed=in_args.random_seed))
        _exit("shuffle_seqs")

    # Split by file number
    if in_args.split_by_x_files:
        file_num = 0
        output_path = ""
        if len(in_args.split_by_x_files[0]) > 2:
            _raise_error(AttributeError("Please provide one or two arguments"), "split_by_x_files")
        if os.path.isfile(str(in_args.sequence[0])):
            input_file_name = in_args.sequence[0]
        else:
            input_file_name = "split_seq.fa"
        for arg in in_args.split_by_x_files[0]:
            try:
                arg = int(arg)
            except ValueError:
                pass
            if type(arg) == int:
                if file_num == 0:
                    file_num = arg
                else:
                    _raise_error(AttributeError("Please provide only one number of files."),
                                 "split_by_x_files")
            elif os.path.isdir(arg):
                output_path = arg
            elif type(arg) == str:
                _raise_error(AttributeError(str(arg) + " is not an existing directory."), "split_by_x_files")
        input_file_name = os.path.splitext(input_file_name)[0] + "_"
        if output_path:
            out_dir = output_path
        else:
            out_dir = os.getcwd()
        if file_num < 1:
            _raise_error(AttributeError("Please provide a valid number of output files.\n"), "split_by_x_files")
        else:
            for idx, sb in enumerate(split_by_x_files(seqbuddy, file_num)):
                output_file_name = "%s%s%s.%s" % (out_dir, os.path.sep, input_file_name + str(idx),
                                                  br.format_to_extension[sb.out_format])
                br._stderr("New file: %s\n" % output_file_name, in_args.quiet)
                sb.write(output_file_name)
        _exit("split_by_x_files")

    # Split by x seqs
    if in_args.split_by_x_seqs:
        seq_num = 0
        output_path = ""
        if len(in_args.split_by_x_seqs[0]) > 2:
            _raise_error(AttributeError("Please provide one or two arguments"), "split_by_x_seqs")
        if os.path.isfile(str(in_args.sequence[0])):
            input_file_name = in_args.sequence[0]
        else:
            input_file_name = "split_seq.fa"
        for arg in in_args.split_by_x_seqs[0]:
            try:
                arg = int(arg)
            except ValueError:
                pass
            if type(arg) == int:
                if seq_num == 0:
                    seq_num = arg
                else:
                    _raise_error(AttributeError("Please provide only one number of sequences."),
                                 "split_by_x_seqs")
            elif os.path.isdir(arg):
                output_path = arg
            elif type(arg) == str:
                _raise_error(AttributeError(str(arg) + " is not an existing directory."), "split_by_x_seqs")
        input_file_name = os.path.splitext(input_file_name)[0] + "_"
        if output_path:
            out_dir = output_path
        else:
            out_dir = os.getcwd()
        if seq_num < 1:
            _raise_error(AttributeError("Please provide a valid number of sequences."), "split_by_x_seqs")
        else:
            for idx, sb in enumerate(split_by_x_seqs(seqbuddy, seq_num)):
                output_file_name = "%s%s%s.%s" % (out_dir, os.path.sep, input_file_name + str(idx),
                                                  br.format_to_extension[sb.out_format])
                br._stderr("New file: %s\n" % output_file_name, in_args.quiet)
                sb.write(output_file_name)
        _exit("split_by_x_seqs")

    # Tail
    if in_args.tail:
        head_size = 1 if in_args.tail[0] is None else in_args.tail[0]
        tail(seqbuddy, head_size)
        _print_recs(seqbuddy)
        _exit("tail")

    # Taxonomic breakdown
    if in_args.taxonomic_breakdown:
        depth = 5 if in_args.taxonomic_breakdown[0] is None else abs(in_args.taxonomic_breakdown[0])
        print(taxonomic_breakdown(seqbuddy, depth))
        _exit("taxonomic_breakdown")

    # Transcribe
    if in_args.transcribe:
        try:
            _print_recs(dna2rna(seqbuddy))
        except TypeError as e:
            _raise_error(e, "transcribe", "DNA sequence required, not")
        _exit("transcribe")

    # Translate CDS
    if in_args.translate:
        if seqbuddy.alpha == "protein":
            _raise_error(TypeError("Nucleic acid sequence required, not protein."), "translate")
        try:
            _print_recs(translate_cds(seqbuddy, quiet=in_args.quiet))
        except TypeError as e:
            _raise_error(e, "translate", ["Nucleic acid sequence required, not protein.", "Record .+ is protein."])
        _exit("translate")

    # Translate 6 reading frames
    if in_args.translate6frames:
        if seqbuddy.alpha == 'protein':
            _raise_error(TypeError("You need to supply DNA or RNA sequences to translate"), "translate6frames")
        try:
            seqbuddy = translate6frames(seqbuddy)
        except TypeError as e:
            _raise_error(e, "translate6frames", ["Nucleic acid sequence required, not protein.", " is protein."])
        if in_args.out_format:
            seqbuddy.out_format = in_args.out_format
        _print_recs(seqbuddy)
        _exit("translate6frames")

    # Transmembrane domains
    if in_args.transmembrane_domains:
        try:
            if not in_args.transmembrane_domains[0]:
                seqbuddy = transmembrane_domains(seqbuddy, quiet=in_args.quiet, keep_temp=in_args.keep_temp)
            else:
                seqbuddy = transmembrane_domains(seqbuddy, job_ids=in_args.transmembrane_domains[0],
                                                 quiet=in_args.quiet, keep_temp=in_args.keep_temp)

            if not in_args.out_format and seqbuddy.out_format not in ["gb", "genbank", "embl"]:
                seqbuddy.out_format = "gb"

        except ImportError as e:
            _raise_error(e, "transmembrane_domains", "Please install the 'suds' package")
        except ValueError as e:
            _raise_error(e, "transmembrane_domains", "is too large to send to TOPCONS. Max record size is 9Mb")
        except ConnectionError as e:
            _raise_error(e, "transmembrane_domains",
                         ["Failed to submit TOPCONS job.", "Job failed...\nServer message",
                          "The job seems to have been lost by the server."])
        except FileNotFoundError as e:
            _raise_error(e, "transmembrane_domains", ["File lost.", "SeqBuddy does not have the necessary hash-map"])
        _print_recs(seqbuddy)
        _exit("transmembrane_domains")

    # Uppercase
    if in_args.uppercase:
        _print_recs(uppercase(seqbuddy))
        _exit("uppercase")


def main():
    br.preparse_flags()
    initiation = []
    try:
        initiation = argparse_init()  # initiation = [in_agrs, seqbuddy]
        command_line_ui(*initiation)
    except (KeyboardInterrupt, br.GuessError) as _e:
        print(_e)
        return False
    except SystemExit:
        return False
    except Exception as _e:
        func = ""
        for next_arg in vars(initiation[0]):
            if getattr(initiation[0], next_arg) and next_arg in br.sb_flags:
                func = next_arg
                break
        br.send_traceback("SeqBuddy", func, _e, VERSION)
        return False
    return True


if __name__ == '__main__':
    main()
