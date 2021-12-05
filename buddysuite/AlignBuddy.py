#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: AlignBuddy.py
author: Stephen R. Bond
email: biologyguy@gmail.com
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description:
AlignmentBuddy is a general wrapper for popular DNA and protein alignment programs that handles format conversion
and allows maintenance of rich feature annotation following alignment.
"""

# ##################################################### IMPORTS ###################################################### #
from __future__ import print_function

# BuddySuite specific
try:
    import buddy_resources as br
    import SeqBuddy as Sb
except ImportError:
    try:
        import buddysuite.buddy_resources as br
        import buddysuite.SeqBuddy as Sb
    except AttributeError:
        from . import buddy_resources as br
        from . import SeqBuddy as Sb

# Standard library
import sys
import os
from copy import deepcopy
from io import StringIO, TextIOWrapper
import random
import re
from collections import OrderedDict
from shutil import which
from subprocess import Popen, PIPE, CalledProcessError
from math import log, ceil

# Third party
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Data import IUPACData

# ##################################################### WISH LIST #################################################### #
# - Map features from a sequence file over to the alignment
# - Annotate alignment (entire columns or subsets of sequences)
# - Pull random records
# - Back-translate
# - Support for MEGA, NBRF/PIR
# - Pairwise percent identity/similarity (build in support for all substitution matricies)
# - Generate jackknifes

# ################################################ GLOBALS ###################################################### #
GAP_CHARS = ["-", ".", " "]
VERSION = br.Version("AlignBuddy", 1, "4.0", br.contributor_list, {"year": 2021, "month": 12, "day": 5})


# #################################################### ALIGNBUDDY #################################################### #
class AlignBuddy(object):
    """
    Core class.
    Open a file or read a handle and parse, or convert raw into a Seq object
    """
    def __init__(self, _input, in_format=None, out_format=None):
        # ####  IN AND OUT FORMATS  #### #
        # Holders for input type. Used for some error handling below
        in_handle = None
        raw_seq = None
        in_file = None
        self.hash_map = OrderedDict()  # This is only filled if hash_ids() is called.
        self.align_tool = OrderedDict()  # This is only filled if generate_msa() is called

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
            in_handle = _input.read()
            _input.seek(0)

        # Plain text in a specific format
        if type(_input) == str and not os.path.isfile(_input):
            raw_seq = br.utf_encode(_input)
            temp = StringIO(_input)
            _input = temp
            #_input.seek(0)

        # File paths
        try:
            if os.path.isfile(_input):
                in_file = _input
                with open(_input, "r", encoding="utf-8") as ifile:
                    _input = StringIO(ifile.read())
        except TypeError:  # This happens when testing something other than a string.
            pass

        self.in_format = br.parse_format(in_format) if in_format else br.guess_format(_input)

        if self.in_format == "empty file":
            raise br.GuessError("Empty file")

        if not self.in_format:
            if in_file:
                raise br.GuessError("Could not determine format from _input file '%s'.\n"
                                    "Try explicitly setting with -f flag." % in_file)
            elif raw_seq:
                raise br.GuessError("Could not determine format from raw input\n --> %s ...\n"
                                    "Try explicitly setting with -f flag." % str(raw_seq)[:50])
            elif in_handle:
                raise br.GuessError("Could not determine format from input file-like object\n"
                                    "Try explicitly setting with -f flag.")
            else:  # This should be unreachable.
                raise br.GuessError("Unable to determine format or input type. "
                                    "Please check how AlignBuddy is being called.")

        self.out_format = self.in_format if not out_format else br.parse_format(out_format)
        # ####  ALIGNMENTS  #### #
        if _input.__class__.__name__ == 'AlignBuddy':
            alignments = _input.alignments

        elif isinstance(_input, list):
            # make sure that the list is actually MultipleSeqAlignment objects
            sample = _input if len(_input) < 5 else random.sample(_input, 5)
            for _seq in sample:
                if type(_seq) != MultipleSeqAlignment:
                    raise TypeError("Alignlist is not populated with MultipleSeqAlignment.")
            alignments = _input

        elif str(type(_input)) == "<class '_io.TextIOWrapper'>" or isinstance(_input, StringIO):
            if self.in_format == "phylipss":
                alignments = list(br.phylip_sequential_read(_input.read(), relaxed=False))
            elif self.in_format == "phylipsr":
                alignments = list(br.phylip_sequential_read(_input.read()))
            else:
                alignments = list(AlignIO.parse(_input, self.in_format))

        elif os.path.isfile(_input):
            with open(_input, "r", encoding='utf-8') as _input:
                if self.in_format == "phylipss":
                    alignments = list(br.phylip_sequential_read(_input.read(), relaxed=False))
                elif self.in_format == "phylipsr":
                    alignments = list(br.phylip_sequential_read(_input.read()))
                else:
                    alignments = list(AlignIO.parse(_input, self.in_format))

        else:  # May be unreachable
            alignments = None

        self.alpha = guess_alphabet(alignments)
        for alignment in alignments:
            alignment._alphabet = self.alpha
            for rec in alignment:
                rec.annotations["molecule_type"] = self.alpha
        self.alignments = alignments
        self.memory_footprint = sum([len(rec) for rec in self.records()])

    def __str__(self):
        empty_alignments = []
        for indx, alignment in enumerate(self.alignments):
            if not len(alignment):
                empty_alignments.append(indx)
        empty_alignments = sorted(empty_alignments, reverse=True)
        for indx in empty_alignments:
            del self.alignments[indx]

        if len(self.alignments) == 0:
            return "AlignBuddy object contains no alignments.\n"

        # There is a weird bug in genbank write() that concatenates dots to the organism name (if set).
        # The following is a work around...
        if self.out_format in ["gb", "genbank"]:
            for rec in self.records_iter():
                try:
                    if re.search(r"(\. )+", rec.annotations['organism']):
                        rec.annotations['organism'] = "."
                except KeyError:
                    pass

        self.out_format = self.out_format.lower()
        multiple_alignments_unsupported = ["fasta", "gb", "genbank", "nexus"]
        if self.out_format in multiple_alignments_unsupported and len(self.alignments) > 1:
            raise ValueError("%s format does not support multiple alignments in one file.\n" % self.out_format)

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
        if not output:
            tmp_file = br.TempFile()
            with open(tmp_file.path, "w", encoding="utf-8") as ofile:
                try:
                    AlignIO.write(self.alignments, ofile, self.out_format)
                except ValueError as e:
                    if "Sequences must all be the same length" in str(e):
                        br._stderr("Warning: Alignment format detected but sequences are different lengths. "
                                   "Format changed to fasta to accommodate proper printing of records.\n\n")
                        AlignIO.write(self.alignments, ofile, "fasta")
                    elif "Repeated name" in str(e) and self.out_format == "phylip":
                        br._stderr("Warning: Phylip format returned a 'repeat name' error, probably due to truncation. "
                                   "Format changed to phylip-relaxed.\n")
                        AlignIO.write(self.alignments, ofile, "phylip-relaxed")
                    else:
                        raise e

            with open(tmp_file.path, "r", encoding="utf-8") as ifile:
                output = ifile.read()
        if self.out_format == "clustal":
            return "%s\n\n" % output.rstrip()
        else:
            return "%s\n" % output.rstrip()

    def __len__(self):
        return len(self.alignments)

    def set_format(self, in_format):
        self.out_format = br.parse_format(in_format)

    def records_iter(self):
        for alignment in self.alignments:
            for rec in alignment:
                yield rec

    def records(self):
        seq_recs = []
        for alignment in self.alignments:
            for rec in alignment:
                seq_recs.append(rec)
        return seq_recs

    def records_dict(self):  # Note that multiple records can have the same ID. Each item in the dict is a list of recs
        seq_recs = OrderedDict()
        for rec in self.records():
            seq_recs.setdefault(rec.id, [])
            seq_recs[rec.id].append(rec)
        return seq_recs

    def lengths(self):
        lengths = [alignment.get_alignment_length() for alignment in self.alignments]
        return lengths

    def write(self, file_path, out_format=None):
        with open(file_path, "w", encoding="utf-8") as ofile:
            if out_format:
                out_format_save = str(self.out_format)
                self.set_format(out_format)
                ofile.write(str(self))
                self.set_format(out_format_save)
            else:
                ofile.write(str(self))
        return


# ################################################# HELPER FUNCTIONS ################################################# #
def guess_alphabet(alignments):
    """
    :param alignments: Duck typed --> AlignBuddy object, list of alignment objects, or a single alignment object
    :return:
    """
    if alignments.__class__.__name__ == 'AlignBuddy':
        align_list = alignments.alignments
    elif type(alignments) == list:
        align_list = alignments
    else:
        align_list = [alignments]
    seq_list = []
    for alignment in align_list:
        seq_list += [str(x.seq) for x in alignment]

    sequence = "".join(seq_list).upper()
    sequence = re.sub(r"[NX\-?]", "", sequence)

    if len(sequence) == 0:
        return None

    if 'U' in sequence:  # U is unique to RNA
        return "RNA"

    percent_dna = len(re.findall("[ATCG]", sequence)) / float(len(sequence))
    if percent_dna > 0.85:  # odds that a sequence with no Us and such a high ATCG count be anything but DNA is low
        return "DNA"
    else:
        return "protein"


def make_copy(alignbuddy):
    alphabet_list = [rec.annotations["molecule_type"] for rec in alignbuddy.records()]
    _copy = deepcopy(alignbuddy)
    _copy.alpha = alignbuddy.alpha
    for indx, rec in enumerate(_copy.records()):
        rec.annotations["molecule_type"] = alphabet_list[indx]
    return _copy


class FeatureReMapper(object):
    """
    Build a list that maps original alignment columns to new positions if columns have been removed
    This will not work if new columns are being added.
    :usage: Instantiate a new object, and for each column in the original alignment, call the 'extend' method,
            specifying whether that column exists in the new alignment or not. Remap the features on the new alignment
            by calling the remap_features method.
    """
    def __init__(self, max_len):
        self.position_map = [None for _ in range(max_len)]
        self.position_map_indx = 0
        self.starting_position_filled = False

    def extend(self, exists=True):
        """
        Iterates the position map, adding an index and whether the new alignment contains the column.
        :param exists: Specify whether the next column exists or not
        """
        if self.position_map[0] is None:
            if exists:
                self.starting_position_filled = True
            self.position_map[self.position_map_indx] = (0, exists)
        else:
            if exists and not self.starting_position_filled:
                self.position_map[self.position_map_indx] = (0, True)
                self.starting_position_filled = True
            elif exists and self.starting_position_filled:
                self.position_map[self.position_map_indx] = (self.position_map[self.position_map_indx - 1][0] + 1, True)
            else:
                self.position_map[self.position_map_indx] = (self.position_map[self.position_map_indx - 1][0], False)
        self.position_map_indx += 1
        return

    def remap_features(self, old_alignment, new_alignment):
        """
        Add all the features from old_alignment that still exist onto new_alignment
        This is quite efficient, as it doesn't nedd to touch the whole sequence
        :param old_alignment: AlignRecord
        :param new_alignment: AlignRecord
        """
        new_records = list(new_alignment)
        for indx, rec in enumerate(old_alignment):
            new_features = []
            for feature in rec.features:
                feature = self._remap(feature)
                if feature:
                    new_features.append(feature)
            new_records[indx].features = new_features
            new_records[indx].annotations = rec.annotations
            new_records[indx].dbxrefs = rec.dbxrefs
        return

    def _remap(self, feature):
        """
        Deal with the weirdness that is SeqRecord features... Compares old features against self.position_map
        :param feature: A feature from the old alignment
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

    def append_pos_map(self, alignment):
        alignment.position_map = self.position_map[:self.position_map_indx + 1]
        return


# ################################################ MAIN API FUNCTIONS ################################################ #
def alignment_lengths(alignbuddy):
    """
    Returns a list of alignment lengths
    :param alignbuddy: The AlignBuddy object to be analyzed
    :type alignbuddy: AlignBuddy
    :return: A list of alignment lengths
    :rtype: list
    """
    output = []
    for alignment in alignbuddy.alignments:
        output.append(alignment.get_alignment_length())
    return output


def bootstrap(alignbuddy, num_bootstraps=1, r_seed=None):
    """
    Sample len(alignbuddy) columns with replacement, and make new alignment(s)
    :param alignbuddy: The AlignBuddy object to be bootstrapped
    :type alignbuddy: AlignBuddy
    :param num_bootstraps: The number of new alignments to be generated
    :param r_seed: Set a seed value so 'random' numbers are reproducible
    :rtype: AlignBuddy
    """
    rand_gen = random.Random(r_seed)

    new_alignments = []
    for align_indx, alignment in enumerate(alignbuddy.alignments):
        alignment = [[res for res in str(seq.seq)] for seq in alignment]
        align_length = len(alignment[0])
        num_seqs = len(alignment)
        for _ in range(num_bootstraps):
            # Instantiate an empty list of lists
            align_as_lists = [['' for _ in range(align_length)] for _ in range(num_seqs)]
            for res_indx in range(align_length):
                column = rand_gen.randint(0, align_length - 1)
                for seq_indx in range(num_seqs):
                    align_as_lists[seq_indx][res_indx] = alignment[seq_indx][column]

            align_as_lists = ["".join(seq) for seq in align_as_lists]
            alb_copy = make_copy(alignbuddy)
            for copy_indx, seq in enumerate(alb_copy.alignments[align_indx]):
                seq.seq = Seq(align_as_lists[copy_indx])
                seq.features = []
                align_as_lists[copy_indx] = seq

            new_alignment = MultipleSeqAlignment(align_as_lists)
            new_alignments.append(new_alignment)

    alignbuddy = AlignBuddy(new_alignments, out_format=alignbuddy.out_format)
    if alignbuddy.out_format == "nexus":
        alignbuddy.out_format = "phylip-relaxed"
    return alignbuddy


def clean_seq(alignbuddy, ambiguous=True, rep_char="N", skip_list=None):
    """
    Remove all non-sequence charcters from sequence strings (wraps SeqBuddy function)
    :param alignbuddy: AlignBuddy object
    :param rep_char: What character should be used to replace ambiguous characters
    :param ambiguous: Specifies whether ambiguous characters should be kept or not
    :param skip_list: Optional list of characters to be left alone
    :return: The cleaned AlignBuddy object
    :rtype: AlignBuddy
    """
    records = alignbuddy.records()
    # Protect gaps from being cleaned by Sb.clean_seq
    for rec in records:
        if rec.annotations["molecule_type"] == "protein":
            rec.seq = Seq(re.sub(r"\*", "-", str(rec.seq)))
        rec.seq = Seq(re.sub("-", "\a", str(rec.seq)))

    skip_list = "\a" if not skip_list else "\a" + "".join(skip_list)

    seqbuddy = Sb.SeqBuddy(records, alpha=alignbuddy.alpha)
    Sb.clean_seq(seqbuddy, ambiguous, rep_char, skip_list)
    for rec in records:
        rec.seq = Seq(re.sub("\a", "-", str(rec.seq)))
    return alignbuddy


def concat_alignments(alignbuddy, group_pattern=None, align_name_pattern="", suppress_position=False):
    """
    Concatenates two or more alignments together, end-to-end
    :param alignbuddy: AlignBuddy object
    :param group_pattern: Regex that matches some regular part of the sequence IDs, dictating who is bound to who
    :param align_name_pattern: Regex that matches something for the whole alignment
    :param suppress_position: Do not append new features denoting the location of each concatenated sequence
    :return: AlignBuddy object containing a single concatenated alignment
    :rtype: AlignBuddy
    """
    if len(alignbuddy.alignments) < 2:
        raise AttributeError("Please provide at least two alignments.")

    concat_groups = OrderedDict()

    if not group_pattern:
        min_length = 1
        for indx, align in enumerate(alignbuddy.alignments):
            max_rec_id_len = 0
            for rec in align:
                max_rec_id_len = len(rec.id) if len(rec.id) > max_rec_id_len else max_rec_id_len
            breakout = False
            while min_length <= max_rec_id_len and not breakout:
                breakout = True
                test_list = []
                for indx2, rec in enumerate(align):
                    if indx2 != 0 and rec.id[:min_length] in test_list:
                        min_length += 1
                        breakout = False
                        break
                    else:
                        test_list.append(rec.id[:min_length])

        group_pattern = "." * min_length

    for indx, align in enumerate(alignbuddy.alignments):
        for rec in align:
            match = re.search(group_pattern, rec.id)
            if match:
                if match.groups():
                    match = "".join(match.groups())
                else:
                    match = match.group(0)
                concat_groups.setdefault(match, [None for _ in range(len(alignbuddy.alignments))])
            else:
                raise ValueError("No match found for record %s in Alignment #%s" % (rec.id, indx + 1))

    for indx, align in enumerate(alignbuddy.alignments):
        for rec in align:
            match = re.search(group_pattern, rec.id)
            if match:
                if match.groups():
                    match = "".join(match.groups())
                else:
                    match = match.group(0)

                if concat_groups[match][indx]:
                    raise ValueError("Replicate matches '%s' in Alignment #%s" % (match, indx + 1))
                concat_groups[match][indx] = rec

        for group, seqs in concat_groups.items():
            if not seqs[indx]:
                concat_groups[group][indx] = SeqRecord(Seq("-" * align.get_alignment_length()))
                concat_groups[group][indx].annotations["molecule_type"] = alignbuddy.alpha

    new_records = [SeqRecord(Seq(""), id=group, features=[]) for group in concat_groups]
    for rec in new_records:
        rec.annotations["molecule_type"] = alignbuddy.alpha
    indx = 0
    for group, seqs in concat_groups.items():
        new_length = 0
        for rec in seqs:
            new_records[indx].seq = Seq(str(new_records[indx].seq) + str(rec.seq))
            rec.features = br.shift_features(rec.features, new_length, new_length + len(rec.seq))
            new_records[indx].features += rec.features
            new_length += len(rec.seq)
        indx += 1

    group_indx = 0
    for group, seqs in concat_groups.items():
        new_length = 0
        align_features = []
        for rec_indx, rec in enumerate(seqs):
            if not suppress_position:
                location = FeatureLocation(new_length, new_length + len(rec.seq))
                match = re.search(align_name_pattern, rec.id)
                if align_name_pattern != "" and match:
                    if match.groups():
                        match = "".join(match.groups())
                    else:
                        match = match.group(0)
                    feature = SeqFeature(location=location, type=match)
                else:
                    if str(rec.id) != "<unknown id>":
                        feature = SeqFeature(location=location, type=rec.id)
                    else:
                        feature = SeqFeature(location=location, type="Alignment_%s" % (rec_indx + 1))
                align_features.append(feature)
            new_length += len(rec.seq)

        new_records[group_indx].features = align_features + new_records[group_indx].features
        group_indx += 1

    alignbuddy.alignments = [MultipleSeqAlignment(new_records)]
    for alignment in alignbuddy.alignments:
        alignment.annotations["molecule_type"] = alignbuddy.alpha
    return alignbuddy


def consensus_sequence(alignbuddy, mode="weighted"):
    # ToDo: include an ambiguous mode that will pull the degenerate nucleotide alphabet in the case of frequency ties.
    # ToDo: Retain sequence annotations --> options to include ALL or CONSERVED
    """
    Generates a simple majority-rule consensus sequence
    :param alignbuddy: The AlignBuddy object to be processed
    :param mode: Simple majority rule or weighted (Henikoff & Henikoff J. Mol. Biol. 1994) [simple, weighted]
    :return: The modified AlignBuddy object (with a single record in each alignment)
    :rtype: AlignBuddy
    """
    consensus_sequences = []

    if mode == "weighted":
        for alignment in alignbuddy.alignments:
            alpha = guess_alphabet(alignment)
            ambig_char = "X" if alpha == "protein" else "N"
            new_seq = ""
            scores_all_pos = {}
            sequence_weights = {}
            for indx in range(alignment.get_alignment_length()):
                residues = {}
                scores_this_pos = {}
                for rec in alignment:
                    residue = rec.seq[indx]
                    residues.setdefault(residue, 0)
                    residues[residue] += 1
                num_residues = len(residues)
                for aminoacid in residues:
                    s = residues[aminoacid]
                    score = float(1) / (num_residues * s)
                    scores_this_pos[aminoacid] = score
                    scores_all_pos[indx] = scores_this_pos

            for rec in alignment:
                seq_score = 0
                for indx in range(alignment.get_alignment_length()):
                    residue_type = str(rec.seq[indx])
                    residue_score = scores_all_pos[indx][residue_type]
                    seq_score += residue_score
                seq_id = rec.id
                sequence_weights[seq_id] = seq_score

            for indx in range(alignment.get_alignment_length()):
                aa_vote_tracker = {}
                for rec in alignment:
                    seq_id = rec.id
                    residue = rec.seq[indx]
                    residue_vote = sequence_weights[seq_id]
                    aa_vote_tracker.setdefault(residue, 0)
                    aa_vote_tracker[residue] += residue_vote

                aa_vote_tracker = sorted(aa_vote_tracker.items(), key=lambda x: x[1], reverse=True)

                if len(aa_vote_tracker) == 1 or aa_vote_tracker[0][1] != aa_vote_tracker[1][1]:
                    new_seq += aa_vote_tracker[0][0]

                else:
                    new_seq += ambig_char

            new_seq = Seq(new_seq)
            description = "Original sequences: %s" % ", ".join([rec.id for rec in alignment])
            new_seq = SeqRecord(new_seq, id="consensus", name="consensus",
                                description=description)
            new_seq.annotations["molecule_type"] = alpha
            consensus_sequences.append(MultipleSeqAlignment([new_seq]))
            for rec in consensus_sequences:
                rec.annotations["molecule_type"] = alpha
    elif mode == "simple":
        for alignment in alignbuddy.alignments:
            alpha = guess_alphabet(alignment)
            ambig_char = "X" if alpha == "protein" else "N"
            new_seq = ""
            for indx in range(alignment.get_alignment_length()):
                residues = {}
                for rec in alignment:
                    residue = rec.seq[indx]
                    residues.setdefault(residue, 0)
                    residues[residue] += 1
                residues = sorted(residues.items(), key=lambda x: x[1], reverse=True)
                if len(residues) == 1 or residues[0][1] != residues[1][1]:
                    new_seq += residues[0][0]
                else:
                    new_seq += ambig_char
            new_seq = Seq(new_seq)
            description = "Original sequences: %s" % ", ".join([rec.id for rec in alignment])
            new_seq = SeqRecord(new_seq, id="consensus", name="consensus",description=description)
            new_seq.annotations["molecule_type"] = alpha
            consensus_sequences.append(MultipleSeqAlignment([new_seq]))
            for rec in consensus_sequences:
                rec.annotations["molecule_type"] = alpha

    else:
        raise ValueError("No valid consensus mode specified (valid modes are 'simple' and 'weighted')")

    alignbuddy.alignments = consensus_sequences
    alignbuddy = trimal(alignbuddy, 'clean')
    return alignbuddy


def delete_invariant_sites(alignbuddy, consider_ambiguous=True):
    """
    Deletes columns where all residues are the same, either explicitly or ambiguously.
    :param alignbuddy: AlignBuddy object
    :param consider_ambiguous: Specify how ambiguous residues are treated ('True' means they can overlap
    the standard code and may be deleted)
    :return: The modified AlignBuddy object
    :rtype: AlignBuddy
    """
    from Bio.Data import IUPACData
    nucs = {'A', 'U', 'C', 'G'} if alignbuddy.alpha == "RNA" else {'A', 'T', 'C', 'G'}
    ambiguiity_codes = IUPACData.ambiguous_rna_values if alignbuddy.alpha == "RNA" \
        else IUPACData.ambiguous_dna_values

    for alignment_index, alignment in enumerate(alignbuddy.alignments):
        keep_columns = []
        if not alignment:
            continue  # Prevent crash if the alignment doesn't have any records in it

        for col_indx in range(alignment.get_alignment_length()):
            col = set([i.upper() for i in alignment[:, col_indx]])

            # Fully invariant columns discarded
            if len(col) == 1:
                continue

            # Simplest scenario, there are differences in the column and we don't care about ambiguous residues
            elif not consider_ambiguous:
                keep_columns.append(col_indx)
                continue

            # Everything else takes ambiguity into consideration
            # Discard fully ambiguous residues
            for i in ["X", "-", ".", "*"]:
                col.discard(i)
            if alignbuddy.alpha != "protein":
                col.discard("N")

            # Again, fully invariant columns discarded
            if len(col) <= 1:
                pass
            elif alignbuddy.alpha == "protein":
                # No partially ambiguous protein characters considered
                keep_columns.append(col_indx)
            # Now it gets trickier... There are many ambiguities possible in nucleotide seqs
            elif len(nucs - col) <= 2:
                # At least two standard residues are present
                keep_columns.append(col_indx)
            else:
                # If we've gotten this far, then there must be at least one ambiguous character
                # Go through each standard residue to see if all col residues can be matched to it
                invariant = False
                for i in nucs:
                    counter = 0
                    for j in col:
                        if i in ambiguiity_codes[j]:
                            counter += 1
                    if counter == len(col):
                        invariant = True
                        break
                if not invariant:
                    keep_columns.append(col_indx)
        keep_columns = [str(i + 1) for i in keep_columns]
        alignment = extract_regions(AlignBuddy([alignment]), ",".join(keep_columns)).alignments[0]
        alignbuddy.alignments[alignment_index] = alignment
    return alignbuddy


def delete_records(alignbuddy, regex):
    """
    Deletes rows with names/IDs matching a search pattern
    :param alignbuddy: AlignBuddy object
    :param regex: The regex pattern to search with (duck typed for list or string)
    :return: The modified AlignBuddy object
    :rtype: AlignBuddy
    """
    if type(regex) == str:
        regex = [regex]
    for indx, pattern in enumerate(regex):
        regex[indx] = ".*" if pattern == "*" else pattern

    regex = "|".join(regex)
    alignments = []
    for alignment in alignbuddy.alignments:
        matches = []
        for record in alignment:
            if not re.search(regex, record.id):
                matches.append(record)
        alignment._records = matches
        alignments.append(alignment)
    alignbuddy.alignments = alignments
    trimal(alignbuddy, "clean")
    return alignbuddy


def dna2rna(alignbuddy):  # Transcribe
    """
    Convert DNA into RNA
    :param alignbuddy: AlignBuddy object
    :return: Modified AlignBuddy object
    :rtype: AlignBuddy
    """
    records = alignbuddy.records()
    seqbuddy = Sb.SeqBuddy(records)
    Sb.dna2rna(seqbuddy)
    alignbuddy.alpha = "RNA"
    return alignbuddy


def enforce_triplets(alignbuddy):
    """
    Organizes nucleotide alignment into triplets
    :param alignbuddy: AlignBuddy object
    :return: The rearranged AlignBuddy object
    :rtype: AlignBuddy
    """
    if alignbuddy.alpha == "protein":
        raise TypeError("Nucleic acid sequence required, not protein.")

    alignbuddy_copy = make_copy(alignbuddy)
    for rec in alignbuddy.records_iter():
        if rec.annotations["molecule_type"] == "protein":
            raise TypeError("Record '%s' is protein. Nucleic acid sequence required." % rec.name)

        seq_list = list(str(rec.seq))
        output = [seq_list[0]] + ([None] * (len(seq_list) - 1))
        output_indx = 0
        held_residues = []
        position = 2
        for residue in seq_list[1:]:
            if position == 1:
                while len(held_residues) >= 3:
                    for res in held_residues[:3]:
                        output_indx += 1
                        output[output_indx] = res
                    held_residues = held_residues[3:]
                position = len(held_residues) + 1
                for res in held_residues:
                    output_indx += 1
                    output[output_indx] = res
                held_residues = []

            if position == 1 \
                or (output[output_indx] not in GAP_CHARS and residue not in GAP_CHARS) \
                    or (output[output_indx] in GAP_CHARS and residue in GAP_CHARS):
                output_indx += 1
                output[output_indx] = residue

            else:
                held_residues.append(residue)
                continue

            if position != 3:
                position += 1
            else:
                position = 1
        back_indx = -1
        while output[back_indx] in ["-", None]:
            back_indx -= 1
        back_indx += 1
        for res in held_residues:
            output[back_indx] = res
            back_indx += 1

        for i in range(back_indx * -1):
            output[back_indx] = "-"
            back_indx += 1
        rec.seq = Seq("".join(output))
    br.remap_gapped_features(alignbuddy_copy.records(), alignbuddy.records())
    trimal(alignbuddy, "clean")
    return alignbuddy


def extract_feature_sequences(alignbuddy, patterns):
    """
    Pull out specific features from annotated sequences
    :param alignbuddy: SeqBuddy object
    :type alignbuddy: AlignBuddy
    :param patterns: The feature(s) to be extracted
    :type patterns: list str
    :return: Modified AlignBuddy object
    :rtype: AlignBuddy
    """
    if type(patterns) == str:
        patterns = [patterns]

    range_patterns = []
    single_patterns = []
    for pattern in patterns:
        if ":" in pattern:
            range_patterns.append(pattern.split(":"))
        else:
            single_patterns.append(pattern)

    # Note that there isn't currently a way to store multiple annotated alignments, but still treat it like there is a
    # way in case this changes in the future
    new_alignments = []
    for alignment in alignbuddy.alignments:
        keep_ranges = []
        for pat in single_patterns:
            matches = []
            for rec in alignment:
                for feat in rec.features:
                    if re.search(pat, feat.type):
                        matches.append([int(feat.location.start), int(feat.location.end)])
                    else:
                        breakout = False
                        for qual_type, quals in feat.qualifiers.items():
                            for qual in quals:
                                if re.search(pat, qual):
                                    matches.append([int(feat.location.start), int(feat.location.end)])
                                    breakout = True
                                    break
                            if breakout:
                                break
            if matches:
                matches = sorted(matches, key=lambda x: x[0])
                start, end = matches[0]
                for next_start, next_end in matches[1:]:
                    if end < next_start:
                        keep_ranges.append([start, end])
                        start, end = next_start, next_end
                    elif end < next_end:
                        end = next_end
                keep_ranges.append([start, end])

        for pat in range_patterns:
            start, end = len(alignment[0]), 0
            pat1, pat2 = False, False
            for rec in alignment:
                for feat in rec.features:
                    if re.search(pat[0], feat.type):
                        start = int(feat.location.start) if int(feat.location.start) < start else start
                        end = int(feat.location.end) if int(feat.location.end) > end else end
                        pat1 = True
                    if re.search(pat[1], feat.type):
                        start = int(feat.location.start) if int(feat.location.start) < start else start
                        end = int(feat.location.end) if int(feat.location.end) > end else end
                        pat2 = True
            if pat1 and pat2:
                keep_ranges.append([start, end])

        if not keep_ranges:
            for rec in alignment:
                rec.seq = Seq("")
                rec.features = []
            new_alignments.append(alignment)
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
            alignment = AlignBuddy([alignment])
            alignment = extract_regions(alignment, final_positions)
            new_alignments.append(alignment.alignments[0])
    alignbuddy.alignments = new_alignments
    return alignbuddy


def extract_regions(alignbuddy, positions):
    """
    Extracts all columns within a given range
    :param alignbuddy: AlignBuddy object
    :param positions: Position code describing which residues to pull (str)
    :return: The modified AlignBuddy object
    :rtype: AlignBuddy

    Position Code:  - Always a string
                    - Comma-separated
                    - Three types of extraction:
                        - Singlets: "2,5,9,-5"
                        - Ranges: "40:75,89:100,432:-45"
                        - mth of nth: "1/5,3/5"
    """
    for indx, alignment in enumerate(alignbuddy.alignments):
        seqbuddy = Sb.SeqBuddy([rec for rec in alignment])
        seqbuddy = Sb.extract_regions(seqbuddy, positions)
        alignbuddy.alignments[indx] = MultipleSeqAlignment(seqbuddy.records)
        alignbuddy.alignments[indx].annotations["molecule_type"] = alignbuddy.alpha
    return alignbuddy


def faux_alignment(seqbuddy, size=0, r_seed=None):
    """
    Creates a meaningless alignment out of a collection of sequences. Gaps are added into each sequence at random
    to create the appropriate length
    :param seqbuddy:
    :param size: The final length of the alignment (defaults to the length of the longest sequence)
    :type size: int
    :param r_seed: Set the random generator seed value
    :return:
    """
    from scipy.stats import gamma
    rand_gen = random.Random(r_seed)
    gamma_50 = gamma.ppf(0.5, 4)
    gamma_scale = 20 / gamma_50

    def gamma_draw(remaining_gaps):
        draw = gamma.ppf(rand_gen.random(), 4)  # sample the gamma distribution
        draw *= gamma_scale  # Scale to the empirically determined value
        draw /= 100  # Convert to percentage
        draw = 1 if draw > 1 else draw  # Ensure we don't go over 100%
        return ceil(remaining_gaps * draw)

    seqbuddy = Sb.clean_seq(seqbuddy, skip_list="*")
    max_length = Sb.max_records(Sb.make_copy(seqbuddy))
    max_length = len(max_length.records[0])
    max_length = size if size and size > max_length else max_length
    records = []
    for rec in Sb.make_copy(seqbuddy).records:
        gaps = max_length - len(rec)
        seq = list(str(rec.seq))
        while gaps > 0:
            next_gap = gamma_draw(gaps)
            gaps -= next_gap
            r = rand_gen.randint(1, len(seq) - 1)
            seq = seq[:r] + ['-' * next_gap] + seq[r:]
        rec.seq = Seq("".join(seq))
        records.append(rec)
    gapped_seqbuddy = Sb.SeqBuddy(records, out_format='fasta')
    alignbuddy = AlignBuddy(str(gapped_seqbuddy), out_format=seqbuddy.out_format)
    alignbuddy = map_features2alignment(seqbuddy, alignbuddy)
    return alignbuddy


# ToDo: Completely refactor the handling of output formats
def generate_msa(seqbuddy, alias, params=None, keep_temp=None, quiet=False):
    """
    Calls sequence aligning tools to generate multiple sequence alignments
    :param seqbuddy: The SeqBuddy object containing the sequences to be aligned
    :param alias: The alignment tool to be used (pagan/prank/muscle/clustalw2/clustalomega/mafft)
    :param params: Additional parameters to be passed to the alignment tool
    :param keep_temp: Determines if/where the temporary files will be kept
    :param quiet: Suppress stderr output
    :return: An AlignBuddy object containing the alignment produced.
    :rtype: AlignBuddy
    """
    if params is None:
        params = ''

    # Figure out what tool is being used
    tool = br.identify_msa_program(alias)

    if not tool:
        raise AttributeError("{0} is not a recognized alignment tool. "
                             "Please check your spelling (case sensitive)".format(alias))

    if keep_temp and os.path.exists(keep_temp):
        check = br.ask("{0} already exists, so files may be overwritten. Proceed [yes]/no?".format(keep_temp))
        if not check:
            sys.exit()
        keep_temp = os.path.abspath(keep_temp)

    if not which(alias):
        error_msg = '#### Could not find %s on your system. ####\n ' \
                    'Please check that your spelling is correct (case sensitive) or go to %s to install %s.' \
                    % (alias, tool["url"], tool["name"])
        raise SystemError(error_msg)
    else:
        valve = br.SafetyValve(global_reps=10)
        Sb.hash_ids(seqbuddy, 8)
        while True:
            valve.step("Generate alignment is failing to create temporary files. Please report this to "
                       "the BuddySuite developers if recurring.")
            try:
                tmp_dir = br.TempDir()
                tmp_in = "%s%stmp.fa" % (tmp_dir.path, os.sep)

                params = re.split(' ', params)

                copy_outfmt = str(seqbuddy.out_format)
                seqbuddy.out_format = 'fasta'
                seqbuddy.write(tmp_in)
                seqbuddy.out_format = copy_outfmt

                # Catch output parameters if passed into the third party program
                for indx, param in enumerate(params):
                    # PAGAN
                    if param == "-f":
                        try:
                            seqbuddy.out_format = br.parse_format(params[indx + 1])
                        except (TypeError, AttributeError, IndexError):
                            pass
                        del params[indx + 1]
                        del params[indx]
                        break
                    # PRANK
                    elif param.startswith("-f="):
                        param = re.match("-f=(.+)", param)
                        try:
                            seqbuddy.out_format = br.parse_format(param.group(1))
                        except (TypeError, AttributeError):
                            pass
                        del params[indx]
                        break
                    # ClustalOmega
                    elif param.startswith("--outfmt="):
                        param = re.match("--outfmt=(.+)", param)
                        try:
                            seqbuddy.out_format = br.parse_format(param.group(1))
                        except (TypeError, AttributeError):
                            pass
                        del params[indx]
                        break
                    # ClustalW2
                    elif param.startswith("-output="):
                        param = re.match("-output=(.+)", param)
                        try:
                            seqbuddy.out_format = br.parse_format(param.group(1))
                        except (TypeError, AttributeError):
                            pass
                        del params[indx]
                        break

                params = ' '.join(params)

                if tool["name"] == 'clustalo':
                    command = '{0} {1} -i {2} -o {3}{4}result -v'.format(alias, params, tmp_in, tmp_dir.path, os.sep)
                elif tool["name"] == 'clustalw':
                    command = '{0} -infile={1} {2} -outfile={3}{4}result'.format(alias, tmp_in, params,
                                                                                 tmp_dir.path, os.sep)
                elif tool["name"] == 'muscle':
                    command = '{0} -in {1} {2}'.format(alias, tmp_in, params)
                elif tool["name"] == 'prank':
                    command = '{0} -d={1} {2} -o={3}{4}result'.format(alias, tmp_in, params, tmp_dir.path, os.sep)
                elif tool["name"] == 'pagan':
                    command = '{0} -s {1} {2} -o {3}{4}result'.format(alias, tmp_in, params, tmp_dir.path, os.sep)
                    print(command)
                elif tool["name"] == 'mafft':
                    command = '{0} {1} {2}'.format(alias, params, tmp_in)

                try:
                    if tool["name"] in ['prank', 'pagan', 'clustalo']:
                        if tool["name"] == 'pagan' and len(seqbuddy) < 4:
                            raise ValueError("PAGAN cannot run when less than 4 sequences are passed in.")
                        if quiet:
                            output = Popen(command, shell=True, universal_newlines=True,
                                           stdout=PIPE, stderr=PIPE).communicate()
                        else:
                            output = Popen(command, shell=True, universal_newlines=True,
                                           stdout=sys.stderr).communicate()
                    else:
                        if quiet:
                            output = Popen(command, shell=True, stdout=PIPE, stderr=PIPE).communicate()
                        else:
                            output = Popen(command, shell=True, stdout=PIPE).communicate()
                        output = output[0].decode("utf-8")
                except CalledProcessError:
                    br._stderr('\n#### {0} threw an error. Scroll up for more info. ####\n\n'.format(tool["name"]),
                               quiet)
                    sys.exit()

                if tool["name"].startswith('clustal'):
                    with open('{0}{1}result'.format(tmp_dir.path, os.path.sep), "r", encoding="utf-8") as result:
                        output = result.read()
                elif tool["name"] == 'prank':
                    possible_files = os.listdir(tmp_dir.path)
                    filename = 'result.best.fas'
                    for _file in possible_files:
                        if 'result.best' in _file and "fas" in _file:
                            filename = _file
                    with open('{0}{1}{2}'.format(tmp_dir.path, os.path.sep, filename), "r", encoding="utf-8") as result:
                        output = result.read()
                elif tool["name"] == 'pagan':
                    with open('{0}{1}result.fas'.format(tmp_dir.path, os.path.sep), "r", encoding="utf-8") as result:
                        output = result.read()
                    if os.path.isfile(".%swarnings" % os.path.sep):  # Pagan spits out this empty file
                        os.remove(".%swarnings" % os.path.sep)

                # Fix broken outputs to play nicely with AlignBuddy parsers
                if (tool["name"] == 'mafft' and '--clustalout' in params) or \
                        (tool["name"] == 'clustalw' and '-output' not in params) or \
                        (tool["name"] == 'clustalo' and ('clustal' in params or '--outfmt clu' in params or
                         '--outfmt=clu' in params)):
                    # Clustal format extra spaces
                    contents = ''
                    prev_line = ''
                    for line in output.splitlines(keepends=True):
                        if line.startswith(' ') and len(line) == len(prev_line) + 1:
                            contents += line[1:]
                        else:
                            contents += line
                        prev_line = line
                    output = contents
                alignbuddy = AlignBuddy(output, out_format=seqbuddy.out_format)

                seqbuddy_recs = []
                for alb_rec in alignbuddy.records():
                    for indx, sb_rec in enumerate(seqbuddy.records):
                        if sb_rec.id == alb_rec.id:
                            seqbuddy_recs.append(sb_rec)
                            del seqbuddy.records[indx]
                            break

                seqbuddy.records = seqbuddy_recs
                # ToDo: Change remap_gapped_features to multicore
                br.remap_gapped_features(seqbuddy_recs, alignbuddy.records())

                for _hash, sb_rec in seqbuddy.hash_map.items():
                    rename(alignbuddy, _hash, sb_rec)

                if keep_temp:
                    # Loop through each saved file and rename any hashes that have been carried over
                    for root, dirs, files in os.walk(tmp_dir.path):
                        for next_file in files:
                            with open("%s%s%s" % (root, os.path.sep, next_file), "r", encoding="utf-8") as ifile:
                                contents = ifile.read()
                            for _hash, sb_rec in seqbuddy.hash_map.items():
                                contents = re.sub(_hash, sb_rec, contents)
                            with open("%s%s%s" % (root, os.path.sep, next_file), "w", encoding="utf-8") as ofile:
                                ofile.write(contents)

                    br.copydir(tmp_dir.path, keep_temp)

                break

            except FileNotFoundError:
                pass
        version = Popen("%s%s" % (alias, tool["ver"]), shell=True, stderr=PIPE, stdout=PIPE).communicate()
        version = version[0].decode() + "\n" + version[1].decode()
        version = re.search(tool['ver_num'], version).group(1)

        alignbuddy.align_tool["tool"] = tool["name"].upper()
        alignbuddy.align_tool["version"] = version

        br._stderr("Returning to AlignBuddy...\n\n", quiet)
        return alignbuddy


def generate_hmm(alignbuddy, alias="hmmbuild"):
    """
    Call hmmerbuild and create an hmm file
    :param alignbuddy:
    :param alias: Specify the binary name of hmmbuild
    :return:
    """
    if not which(alias):
        raise SystemError("Could not find %s on your system. Please check your spelling or install HMMER3." % alias)

    tmp_dir = br.TempDir()
    for align in alignbuddy.alignments:
        align_file = tmp_dir.subfile("align.sto")
        hmm_file = tmp_dir.subfile("align.hmm")
        with open(align_file, "w") as ofile:
            AlignIO.write(align, ofile, "stockholm")
        out, err = Popen("%s %s %s" % (alias, hmm_file, align_file), shell=True, stdout=PIPE, stderr=PIPE).communicate()
        with open(hmm_file, "r") as ifile:
            hmm = ifile.read()
            if not hmm:
                raise SystemError("No output detected after running %s." % alias)
            align.hmm = hmm
            align.hmm_out = out.decode()
    return alignbuddy


def hash_ids(alignbuddy, hash_length=10, r_seed=None):
    """
    Replace all IDs with random hashes
    :param alignbuddy: AlignBuddy object
    :param hash_length: Specifies the length of the new hashed IDs
    :param r_seed: Set the random generator seed value
    :return: The modified AlignBuddy object, with a new attribute `hash_map` added to each alignment object
    :rtype: AlignBuddy
    """
    try:
        hash_length = int(hash_length)
    except ValueError:
        raise TypeError("Hash length argument must be an integer, not %s" % type(hash_length))

    if hash_length < 1:
        raise ValueError("Hash length must be greater than 0")

    if 32 ** hash_length <= len(alignbuddy.records()) * 2:
        raise ValueError("Insufficient number of hashes available to cover all sequences. "
                         "Hash length must be increased.")
    # If a hash_map already exists and fits all the specs, re-apply it.
    if alignbuddy.hash_map:
        alignbuddy_copy = make_copy(alignbuddy)
        re_apply_hash_map = True
        records = alignbuddy_copy.records_dict()
        reverse_hashmap_ids = [rec_id for _hash, rec_id in alignbuddy.hash_map.items()]
        reverse_hashmap_hashes = [_hash for _hash, rec_id in alignbuddy.hash_map.items()]

        for rec_id, rec_list in records.items():
            if rec_id not in reverse_hashmap_ids:
                re_apply_hash_map = False
                break
            for rec in rec_list:
                _hash = reverse_hashmap_hashes[reverse_hashmap_ids.index(rec_id)]
                rec.id = _hash
                rec.name = _hash

        if re_apply_hash_map:
            alignbuddy.alignments = alignbuddy_copy.alignments
            return alignbuddy

    hash_map = OrderedDict()
    for alignment in alignbuddy.alignments:
        temp_seqbuddy = Sb.SeqBuddy(list(alignment))
        Sb.hash_ids(temp_seqbuddy, hash_length, r_seed=r_seed)
        for _hash, _id in temp_seqbuddy.hash_map.items():
            hash_map[_hash] = _id
    alignbuddy.hash_map = hash_map
    return alignbuddy


def lowercase(alignbuddy):
    """
    Converts all sequence residues to lowercase.
    :param alignbuddy: The AlignBuddy object to be modified.
    :return: The modified AlignBuddy object
    :rtype: AlignBuddy
    """
    for rec in alignbuddy.records_iter():
        rec.seq = Seq(str(rec.seq).lower())
    return alignbuddy


def map_features2alignment(seqbuddy, alignbuddy):
    """
    Copy features from an annotated sequence over to its corresponding record in an alignment
    :param seqbuddy: SeqBuddy object
    :param alignbuddy: AlignBuddy object
    :return: The modified AlignBuddy object
    :rtype: AlignBuddy
    """
    def feat_map(feat, _sb_rec, _alb_rec):
        new_location = feat.location
        if type(feat.location) == FeatureLocation:
            chars = 0
            start = None
            end = None
            for indx, residue in enumerate(str(alb_rec.seq)):
                if residue not in GAP_CHARS:
                    chars += 1

                if start is None and chars - 1 == feat.location.start:
                    start = int(indx)

                if chars == feat.location.end:
                    end = int(indx) + 1

                if None not in [start, end]:
                    new_location = FeatureLocation(start, end)
                    break

        else:  # CompoundLocation
            parts = []
            for sub_feature in feat.location.parts:
                parts.append(feat_map(SeqFeature(sub_feature), _sb_rec, _alb_rec).location)
            new_location = CompoundLocation(parts, operator='order')
        feat.location = new_location
        return feat

    alb_recs = alignbuddy.records_dict()
    Sb.clean_seq(seqbuddy)
    for sb_rec in seqbuddy.records:
        sb_rec.features = br.shift_features(sb_rec.features, 0, len(sb_rec.seq))  # Cleans weird start/end positions
        if sb_rec.id in alb_recs:
            for alb_rec in alb_recs[sb_rec.id]:
                for feature in sb_rec.features:
                    alb_rec.features.append(feat_map(feature, sb_rec, alb_rec))
    return alignbuddy


def order_ids(alignbuddy, reverse=False):
    """
    Sorts the alignments by ID, alphabetically
    :param alignbuddy: AlignBuddy object
    :param reverse: Reverses the order
    :return: The modified AlignBuddy object
    :rtype: AlignBuddy
    """
    for indx, alignment in enumerate(alignbuddy.alignments):
        alignment = Sb.SeqBuddy(list(alignment))
        Sb.order_ids(alignment, reverse=reverse)
        alignbuddy.alignments[indx] = MultipleSeqAlignment(alignment.records)
    return alignbuddy


def percent_id(alignbuddy):
    """
    Pairwise comparison of all sequences in the alignment, returning a matrix of percent match (exact).
    Note that the matrix itself is appended to the alignment objects, and does not alter the actual sequences at all
    :param alignbuddy:
    :return:
    """
    for alignment in alignbuddy.alignments:
        matrix = OrderedDict()
        for i, rec1 in enumerate(alignment):
            matrix.setdefault(rec1.id, OrderedDict())
            for rec2 in alignment[i+1:]:
                matrix.setdefault(rec2.id, OrderedDict())
                id_counter = 0
                gap_counter = 0
                for a, b in zip(list(rec1.seq), list(rec2.seq)):
                    if a == b == "-":
                        gap_counter += 1
                        continue
                    id_counter += 1 if a == b else 0
                len_align = len(rec1) - gap_counter
                matrix[rec1.id][rec2.id] = id_counter / len_align
                matrix[rec2.id][rec1.id] = id_counter / len_align
        alignment.percent_ids = matrix
    return alignbuddy


def position_frequency_matrix(alignbuddy):
    """
    Count the frequency of each residue at each position.
    :param alignbuddy:
    :return:
    """
    for alignment, length in zip(alignbuddy.alignments, alignbuddy.lengths()):
        if alignbuddy.alpha == "DNA":
            residues = {res for res in IUPACData.unambiguous_dna_letters}
        elif alignbuddy.alpha == "RNA":
            residues = {res for res in IUPACData.unambiguous_rna_letters}
        else:
            residues = {res for res in IUPACData.protein_letters}

        columns = [{} for _ in range(length)]
        for col_indx in range(length):
            for res in alignment[:, col_indx]:
                residues.add(res)
                columns[col_indx].setdefault(res, 0)
                columns[col_indx][res] += 1

        residues = sorted(residues)
        non_alpha = [char for char in residues if char.isalpha() is False]
        residues = residues[len(non_alpha):] + non_alpha
        pfm = [OrderedDict([(res, 0) for res in residues]) for _ in range(length)]
        num_seqs = len(alignment)
        for indx, col in enumerate(columns):
            for res, count in col.items():
                pfm[indx][res] = round(count / num_seqs, 6)
        setattr(alignment, "pfm", pfm)
    return alignbuddy


def pull_records(alignbuddy, regex, description=False):
    """
    Retrieves rows with names/IDs matching a search pattern
    :param alignbuddy: The AlignBuddy object to be pulled from
    :param regex: List of regex expressions or single regex
    :type regex: str list tuple
    :param description: Allow search in description string
    :return: The modified AlignBuddy object
    :rtype: AlignBuddy
    """
    if type(regex) == str:
        regex = [regex]
    for indx, pattern in enumerate(regex):
        regex[indx] = ".*" if pattern == "*" else pattern

    regex = "|".join(regex)
    alignments = []
    for alignment in alignbuddy.alignments:
        matches = []
        for rec in alignment:
            if re.search(regex, rec.id) or (description and re.search(regex, rec.description)):
                matches.append(rec)
        alignment._records = matches
        alignments.append(alignment)
    alignbuddy.alignments = alignments
    trimal(alignbuddy, "clean")
    return alignbuddy


def rename(alignbuddy, query, replace="", num=0):
    """
    Rename an alignment's sequence IDs
    :param alignbuddy: The AlignBuddy object to be modified
    :param query: The pattern to be searched for
    :param replace: The string to be substituted
    :param num: The maximum number of substitutions to make
    :return: The modified AlignBuddy object
    :rtype: AlignBuddy
    """
    seqbuddy = Sb.SeqBuddy(alignbuddy.records())
    Sb.rename(seqbuddy, query, replace, num)
    return alignbuddy


def rna2dna(alignbuddy):  # Reverse-transcribe
    """
    Convert RNA into DNA.
    :param alignbuddy: AlignBuddy object
    :return: Modified AlignBuddy object
    :rtype: AlignBuddy
    """
    records = alignbuddy.records()
    seqbuddy = Sb.SeqBuddy(records)
    Sb.rna2dna(seqbuddy)
    alignbuddy.alpha = "DNA"
    return alignbuddy


def translate_cds(alignbuddy):
    """
    Translates a nucleotide alignment into a protein alignment.
    :param alignbuddy: The AlignBuddy object to be translated
    :return: The translated AlignBuddy object
    :rtype: AlignBuddy
    """
    if alignbuddy.alpha == "protein":
        raise TypeError("Nucleic acid sequence required, not protein.")

    enforce_triplets(alignbuddy)

    new_aligns = []
    for alignment in alignbuddy.alignments:
        seqbuddy = Sb.SeqBuddy(list(alignment))
        Sb.replace_subsequence(seqbuddy, "\\".join(GAP_CHARS), "-")
        Sb.translate_cds(seqbuddy, alignment=True)
        alignment = MultipleSeqAlignment(seqbuddy.records)
        alignment.annotations = "protein"
        new_aligns.append(alignment)
    alignbuddy.alpha = "protein"
    alignbuddy.alignments = new_aligns
    return alignbuddy


def trimal(alignbuddy, threshold):
    """
    Trims alignment gaps using algorithms from trimAl
    Capella-Gutiérrez, S., Silla-Martínez, J. M., and Gabaldón, T. (2009).
    trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses.
    Bioinformatics 25, 1972–1973. doi:10.1093/bioinformatics/btp348.
    :param alignbuddy: The AlignBuddy object to be trimmed
    :param threshold: The threshold value or trimming algorithm to be used
    :return: The trimmed AlignBuddy object
    :rtype: AlignBuddy
    """
    def gappyout(_gap_distr, _position_map):
        _max_gaps = 0
        # If there are no columns with zero gaps, scan through the distribution to find where the columns start
        for i in _gap_distr:
            if i == 0:
                _max_gaps = i + 1
            else:
                break

        max_slope = -1
        slopes = [-1 for _ in range(len(alignment) + 1)]
        max_iter = len(alignment) + 1
        active_pointer = 0
        while active_pointer < max_iter:
            for i in _gap_distr:
                if i == 0:
                    active_pointer += 1
                else:
                    break

            prev_pointer1 = int(active_pointer)
            if active_pointer + 1 >= max_iter:
                break

            while True:
                active_pointer += 1
                if active_pointer + 1 >= max_iter or _gap_distr[active_pointer] != 0:
                    break

            prev_pointer2 = int(active_pointer)
            if active_pointer + 1 >= max_iter:
                break

            while True:
                active_pointer += 1
                if active_pointer + 1 >= max_iter or _gap_distr[active_pointer] != 0:
                    break

            if active_pointer + 1 >= max_iter:
                break

            slopes[active_pointer] = (active_pointer - prev_pointer2) / len(alignment)
            slopes[active_pointer] /= (_gap_distr[active_pointer] + _gap_distr[prev_pointer2]) / num_columns

            if slopes[prev_pointer1] != -1:
                if slopes[active_pointer] / slopes[prev_pointer1] > max_slope:
                    max_slope = slopes[active_pointer] / slopes[prev_pointer1]
                    _max_gaps = prev_pointer1
            elif slopes[prev_pointer2] != -1:
                if slopes[active_pointer] / slopes[prev_pointer2] > max_slope:
                    max_slope = slopes[active_pointer] / slopes[prev_pointer2]
                    _max_gaps = prev_pointer1

            active_pointer = prev_pointer2

        cuts = []
        for col, gaps in enumerate(each_column):
            if gaps <= _max_gaps:
                cuts.append(col)
                _position_map.extend(True)
            else:
                _position_map.extend(False)

        _new_alignment = alignment[:, 0:0]
        for _col in cuts:
            _new_alignment += alignment[:, _col:_col + 1]

        return _new_alignment

    for alignment_index, alignment in enumerate(alignbuddy.alignments):
        if not alignment:
            continue  # Prevent crash if the alignment doesn't have any records in it
        # gap_distr is the number of columns w/ each possible number of gaps; the index is == to number of gaps
        gap_distr = [0 for _ in range(len(alignment) + 1)]
        num_columns = alignment.get_alignment_length()
        each_column = [0 for _ in range(num_columns)]

        # Each position_map index corresponds to the original column position, values are tuples of the new position
        # and whether the column still exists (True) or has been deleted (False)
        position_map = FeatureReMapper(len(alignment[0]))

        max_gaps = 0
        for indx in range(num_columns):
            num_gaps = len(re.findall("-", str(alignment[:, indx])))
            gap_distr[num_gaps] += 1
            each_column[indx] = num_gaps

        # Remove any columns with any gaps
        if threshold in ["no_gaps", "all"]:
            threshold = 0
            new_alignment = alignment[:, 0:0]
            for next_col, num_gaps in enumerate(each_column):
                if num_gaps <= max_gaps:
                    new_alignment += alignment[:, next_col:next_col + 1]
                    position_map.extend(True)
                else:
                    position_map.extend(False)

        # Remove any columns that contain nothing but gaps
        elif threshold == "clean":
            max_gaps = len(alignment) - 1
            new_alignment = alignment[:, 0:0]
            for next_col, num_gaps in enumerate(each_column):
                if num_gaps <= max_gaps:
                    new_alignment += alignment[:, next_col:next_col + 1]
                    position_map.extend(True)
                else:
                    position_map.extend(False)

        # trimAl algorithm for removing gaps, depending on size of alignment and distribution of seqs
        elif threshold == "gappyout":
            new_alignment = gappyout(gap_distr, position_map)

        elif threshold == "strict":  # ToDo: Implement
            new_alignment = False
        elif threshold == "strictplus":  # ToDo: Implement
            new_alignment = False
        elif type(threshold) in [int, float]:
            if threshold >= 1:
                max_gaps = round(threshold)
            else:
                threshold = 0.0001 if threshold == 0 else threshold
                max_gaps = round(len(alignment) * threshold)

            new_alignment = alignment[:, 0:0]
            for next_col, num_gaps in enumerate(each_column):
                if num_gaps <= max_gaps:
                    new_alignment += alignment[:, next_col:next_col + 1]
                    position_map.extend(True)
                else:
                    position_map.extend(False)
        else:
            raise NotImplementedError("%s not an implemented trimal method" % threshold)

        position_map.remap_features(alignbuddy.alignments[alignment_index], new_alignment)
        position_map.append_pos_map(new_alignment)
        alignbuddy.alignments[alignment_index] = new_alignment

    return alignbuddy


def uppercase(alignbuddy):
    """
    Converts all sequence residues to uppercase.
    :param alignbuddy: The AlignBuddy object to be modified.
    :return: The modified AlignBuddy object
    :rtype: AlignBuddy
    """
    for rec in alignbuddy.records_iter():
        rec.seq = Seq(str(rec.seq).upper())
    return alignbuddy


# ################################################# COMMAND LINE UI ################################################## #
def argparse_init():
    """
    Catching params to prevent weird collisions with alignment program arguments
    :return: in_args from argparse, new AlignBuddy object
    """
    if '--generate_alignment' in sys.argv:
        sys.argv[sys.argv.index('--generate_alignment')] = '-ga'
    if '-ga' in sys.argv:
        ga_indx = sys.argv.index('-ga')
        if len(sys.argv) > ga_indx + 1:
            extra_args = None
            for indx, param in enumerate(sys.argv[ga_indx + 1:]):
                if param in ["-f", "--in_format", "-i", "--in_place", "-k", "--keep_temp", "-o", "--out_format",
                             "-q", "--quiet", "-t", "--test"]:
                    extra_args = ga_indx + 1 + indx
                    break

            if extra_args in [ga_indx + 1, ga_indx + 2]:
                # No conflicts possible, continue on your way
                pass

            elif len(sys.argv) > ga_indx + 2 or extra_args:
                # There must be optional arguments being passed into the alignment tool
                sys.argv[ga_indx + 2] = " %s" % sys.argv[ga_indx + 2].rstrip()

    import argparse

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

    br.flags(parser, ("alignments", "Supply file path(s) or raw alignments. If piping sequences into AlignBuddy this "
                                    "argument must be left blank."),
             br.alb_flags, br.alb_modifiers, VERSION)

    in_args = parser.parse_args()
    br.check_garbage_flags(in_args, "AlignBuddy")

    alignbuddy = []
    align_set = ""

    if in_args.out_format:
        try:
            in_args.out_format = br.parse_format(in_args.out_format)
        except TypeError as e:
            br._stderr("%s\n" % str(e))
            sys.exit()

    try:
        # Some tools do not start with AlignBuddy objs, so skip this for those rare cases
        if not in_args.generate_alignment and not in_args.faux_align:
            for align_set in in_args.alignments:
                if isinstance(align_set, TextIOWrapper) and align_set.buffer.raw.isatty():
                    br._stderr("Warning: No input detected so AlignBuddy is aborting...\n"
                               "For more information, try:\n%s --help\n" % sys.argv[0])
                    sys.exit()
                align_set = AlignBuddy(align_set, in_args.in_format, in_args.out_format)
                alignbuddy += align_set.alignments

            alignbuddy = AlignBuddy(alignbuddy, align_set.in_format, align_set.out_format)

    except br.GuessError as e:
        br._stderr("GuessError: %s\n" % e, in_args.quiet)
        sys.exit()

    except ValueError as e:
        br._stderr("ValueError: %s\n" % e, in_args.quiet)
        sys.exit()

    except br.PhylipError as e:
        br._stderr("PhylipError: %s\n" % e, in_args.quiet)
        sys.exit()

    except TypeError as e:
        if "Format type '%s' is not recognized/supported" % in_args.in_format in str(e):
            br._stderr("TypeError: %s\n" % str(e), in_args.quiet)
            sys.exit()
        else:
            raise e

    in_args.random_seed = None if not in_args.random_seed else in_args.random_seed
    return in_args, alignbuddy


def command_line_ui(in_args, alignbuddy, skip_exit=False, pass_through=False):  # ToDo: Convert to a class
    # ############################################# INTERNAL FUNCTIONS ############################################## #
    def _print_aligments(_alignbuddy):
        try:
            _output = str(_alignbuddy)
        except ValueError as _err:
            br._stderr("ValueError: %s\n" % str(_err))
            return False
        except TypeError as _err:  # Catching biopython/nexus issue
            if " unexpected keyword argument 'interleave'" in str(_err) and _alignbuddy.out_format.startswith("nexus"):
                br._stderr("WARNING: Attempting to write NEXUS format '%s' with an old version of "
                           "BioPython. Switching to default nexus.\n" % _alignbuddy.out_format, in_args.quiet)
                _alignbuddy.out_format = "nexus"
                _output = str(_alignbuddy)
            else:
                raise _err
        if in_args.test:
            br._stderr("*** Test passed ***\n", in_args.quiet)

        elif in_args.in_place:
            _in_place(_output, in_args.alignments[0])

        else:
            br._stdout("%s" % _output)
        return True

    def create_seqbuddy():
        sb = []
        seq_set = None
        for seq_set in in_args.alignments:
            if isinstance(seq_set, TextIOWrapper) and seq_set.buffer.raw.isatty():
                br._stderr("Warning: No input detected so AlignBuddy is aborting...\n"
                           "For more information, try:\n%s --help\n" % sys.argv[0])
                sys.exit()

            seq_set = Sb.SeqBuddy(seq_set, in_args.in_format, in_args.out_format)
            sb += seq_set.records
        if seq_set:
            sb = Sb.SeqBuddy(sb, seq_set.in_format, seq_set.out_format)
        else:
            sb = Sb.SeqBuddy(sb, in_args.in_format, in_args.out_format)
        return sb

    def _in_place(_output, file_path):
        if not os.path.exists(file_path):
            br._stderr("Warning: The -i flag was passed in, but the positional argument doesn't seem to be a "
                       "file. Nothing was written.\n", in_args.quiet)
            br._stderr("%s" % _output, in_args.quiet)
        else:
            with open(os.path.abspath(file_path), "w", encoding="utf-8") as _ofile:
                _ofile.write(_output)
            br._stderr("File overwritten at:\n%s\n" % os.path.abspath(file_path), in_args.quiet)

    def _exit(_tool, skip=skip_exit):
        if skip:
            return
        usage = br.Usage()
        memory_footprint = 0 if type(alignbuddy) != AlignBuddy else alignbuddy.memory_footprint
        usage.increment("AlignBuddy", VERSION.short(), _tool, memory_footprint)
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
    # Alignment lengths
    if in_args.alignment_lengths:
        counts = alignment_lengths(alignbuddy)
        if len(counts) == 1:
            br._stdout("%s\n" % counts[0])
        else:
            for indx, count in enumerate(counts):
                br._stderr("# Alignment %s\n" % (indx + 1), quiet=in_args.quiet)
                br._stdout("%s\n" % count)
        _exit("alignment_lengths")

    # Bootstrap
    if in_args.bootstrap:
        num_bootstraps = in_args.bootstrap[0] if in_args.bootstrap[0] else 1
        _print_aligments(bootstrap(alignbuddy, num_bootstraps=num_bootstraps, r_seed=in_args.random_seed))
        _exit("bootstrap")

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
        _print_aligments(clean_seq(alignbuddy, ambiguous=ambig, rep_char=rep_char))
        _exit("clean_seq")

    # Concatenate Alignments
    if in_args.concat_alignments:
        try:
            args = in_args.concat_alignments[0]
            if args:
                try:
                    group_int = int(args[0])
                    if group_int >= 0:
                        group_pattern = "^%s" % ("." * group_int)
                    else:
                        group_pattern = "%s$" % ("." * abs(group_int))

                except ValueError:
                    group_pattern = br.clean_regex(args[0], in_args.quiet)
                    group_pattern = None if not group_pattern else group_pattern[0]

            else:
                group_pattern = None

            if len(args) > 1:
                try:
                    align_int = int(args[1])
                    if align_int >= 0:
                        align_pattern = "%s$" % ("." * align_int)
                    else:
                        align_pattern = "^%s" % ("." * abs(align_int))

                except ValueError:
                    align_pattern = br.clean_regex(args[1], in_args.quiet)
                    align_pattern = "" if not align_pattern else align_pattern[0]

            else:
                align_pattern = ""
            _print_aligments(concat_alignments(alignbuddy, group_pattern, align_pattern))

        except AttributeError as e:
            _raise_error(e, "concat_alignments", "Please provide at least two alignments.")
        except ValueError as e:
            _raise_error(e, "concat_alignments", ["No match found for record", "Replicate matches"])
        _exit("concat_alignments")

    # Consensus sequence
    if in_args.consensus:
        con = "w" if in_args.consensus[0] is None else in_args.consensus[0]
        con = con.lower()
        con = "simple" if "simple".startswith(con) else con
        con = "weighted" if "weighted".startswith(con) else con
        try:
            alignbuddy = consensus_sequence(alignbuddy, con)
        except ValueError as err:
            _raise_error(err, "consensus", "No valid consensus mode")
        _print_aligments(alignbuddy)
        _exit("consensus")

    # Delete invariant sites
    if in_args.delete_invariant_sites:
        ambig = False if in_args.delete_invariant_sites[0] \
                        and in_args.delete_invariant_sites[0].startswith("amb") else True
        alignbuddy = delete_invariant_sites(alignbuddy, ambig)
        _print_aligments(alignbuddy)
        _exit("delete_invariant_sites")

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

        args = []
        for arg in in_args.delete_records:
            if os.path.isfile(arg):
                with open(arg, "r", encoding="utf-8") as ifile:
                    for line in ifile:
                        args.append(line.strip())
            else:
                args.append(arg)

        args = br.clean_regex(args, in_args.quiet)
        if not args:  # If all regular expression are malformed, exit out gracefully
            _print_aligments(alignbuddy)
            _exit("delete_records")

        pulled = pull_records(make_copy(alignbuddy), args)
        alignbuddy = delete_records(alignbuddy, args)
        deleted_recs = []
        num_deleted = 0
        for alignment in pulled.alignments:
            if alignment:
                deleted_recs.append([])
                for rec in alignment:
                    num_deleted += 1
                    deleted_recs[-1].append(rec.id)
            else:
                deleted_recs.append(None)
        if num_deleted:
            stderr = "# ####################### Deleted records ######################## #\n"
            for indx, alignment in enumerate(deleted_recs):
                if alignment:
                    counter = 1
                    stderr += "# Alignment %s\n" % (indx + 1)
                    for rec_id in alignment:
                        stderr += "%s\t" % rec_id
                        if counter % columns == 0:
                            stderr = "%s\n" % stderr.strip()
                        counter += 1
                    stderr = "%s\n\n" % stderr.strip()
            stderr = "%s\n# ################################################################ #\n" % stderr.strip()
        else:
            stderr = "# ################################################################ #\n"
            stderr += "#     No sequence identifiers match '%s'\n" % "|".join(args)
            stderr += "# ################################################################ #\n"

        br._stderr(stderr, in_args.quiet)
        _print_aligments(alignbuddy)
        _exit("delete_records")

    # Enforce triplets
    if in_args.enforce_triplets:
        try:
            _print_aligments(enforce_triplets(alignbuddy))
        except TypeError as e:
            _raise_error(e, "enforce_triplets", "Nucleic acid sequence required")
        _exit("enforce_triplets")

    # Extact features
    if in_args.extract_feature_sequences:
        patterns = br.clean_regex(in_args.extract_feature_sequences[0], in_args.quiet)
        if patterns:
            alignbuddy = extract_feature_sequences(alignbuddy, patterns)
        _print_aligments(alignbuddy)
        _exit("extract_feature_sequences")

    # Extract regions
    if in_args.extract_regions:
        try:
            args = ",".join(in_args.extract_regions[0])
            alignbuddy = extract_regions(alignbuddy, args)
            _print_aligments(alignbuddy)
        except ValueError as e:
            br._stderr("""
Extraction regions can be specified in the following ways:
    - Singlets: "2,5,9,-5"
    - Ranges: "40:75,89:100,432:-45"
    - mth of nth: "1/5,3/5"

For examples, please see:
https://github.com/biologyguy/BuddySuite/wiki/AB-Extract-regions

""")
            _raise_error(e, "extract_regions", "Unable to decode the positions string")
        _exit("extract_regions")

    # Faux alignment
    if in_args.faux_align:
        seqbuddy = create_seqbuddy()
        in_args.faux_align = 0 if not in_args.faux_align[0] else in_args.faux_align[0]
        _print_aligments(faux_alignment(seqbuddy, in_args.faux_align, r_seed=in_args.random_seed))
        _exit("faux_align")

    # Generate Alignment
    if in_args.generate_alignment:
        # ToDo: The extra arguments parameter probably doesn't need to be dependent on the tool parameter being passed
        args = in_args.generate_alignment[0]
        if not args:
            for tool in ['mafft', 'pagan', 'muscle', 'clustalo', 'clustalomega', 'prank', 'clustalw2', 'clustalw']:
                if which(tool):
                    args = [tool]
                    break
        if not args:
            _raise_error(AttributeError("Unable to identify any supported alignment tools on your system."),
                         "generate_alignment")

        seqbuddy = create_seqbuddy()

        params = re.sub(r'\[(.*)]', r'\1', args[1]) if len(args) > 1 else None

        try:
            alignbuddy = generate_msa(seqbuddy, args[0], params, in_args.keep_temp, in_args.quiet)
            if in_args.out_format:
                alignbuddy.set_format(in_args.out_format)
            _print_aligments(alignbuddy)
        except AttributeError as e:
            _raise_error(e, "generate_alignment", "is not a recognized alignment tool")
        except SystemError as e:
            _raise_error(e, "generate_alignment", "Could not find")
        except ValueError as err:
            _raise_error(err, "generate_alignment", "PAGAN cannot run when less than 4 sequences are passed in")
        _exit("generate_alignment")

    # Generate HMM
    if in_args.generate_hmm:
        alias = "hmmbuild" if not in_args.generate_hmm[0] else in_args.generate_hmm[0]
        try:
            alignbuddy = generate_hmm(alignbuddy, alias)
        except SystemError as err:
            _raise_error(err, "generate_hmm", ["Please check your spelling or install HMMER3",
                                               "No output detected after running"])
        for align in alignbuddy.alignments:
            if not in_args.quiet:
                br._stderr("%s\n" % align.hmm_out)
            br._stdout("%s\n" % align.hmm)
        _exit("generate_hmm")

    # Hash ids
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

        if 32 ** hash_length <= len(alignbuddy.records()) * 2:
            holder = ceil(log(len(alignbuddy.records()) * 2, 32))
            br._stderr("Warning: The hash_length parameter was passed in with the value %s. "
                       "This is too small to properly cover all sequences, so it has been increased to %s.\n\n" %
                       (hash_length, holder), in_args.quiet)
            hash_length = holder
        hash_ids(alignbuddy, hash_length=hash_length, r_seed=in_args.random_seed)

        hash_table = "# Hash table\n"
        for _hash, orig_id in alignbuddy.hash_map.items():
            hash_table += "%s\t%s\n" % (_hash, orig_id)
        hash_table += "\n"

        br._stderr(hash_table, in_args.quiet)
        _print_aligments(alignbuddy)
        _exit("hash_seq_ids")

    # List identifiers
    if in_args.list_ids:
        columns = 1 if not in_args.list_ids[0] else abs(in_args.list_ids[0])
        output = ""
        for indx, alignment in enumerate(alignbuddy.alignments):
            output += "# Alignment %s\n" % str(indx + 1)
            for rec_indx, rec in enumerate(alignment):
                output += "%s\t" % rec.id
                if (rec_indx + 1) % columns == 0:
                    output = "%s\n" % output.strip()
            output = "%s\n\n" % output.strip()
        br._stdout("%s\n" % output.strip())
        _exit("list_ids")

    # Lowercase
    if in_args.lowercase:
        _print_aligments(lowercase(alignbuddy))
        _exit("lowercase")

    # Map features to alignment
    if in_args.mapfeat2align:
        reference_records = []
        for path in in_args.mapfeat2align:
            reference_records += Sb.SeqBuddy(path).records
        seqbuddy = Sb.SeqBuddy(reference_records)
        alignbuddy = map_features2alignment(seqbuddy, alignbuddy)
        in_args.out_format = "genbank" if not in_args.out_format else in_args.out_format
        alignbuddy.set_format(in_args.out_format)
        _print_aligments(alignbuddy)
        _exit("mapfeat2align")

    # Number sequences per alignment
    if in_args.num_seqs:
        if len(alignbuddy.alignments) == 1:
            br._stdout("%s\n" % len(alignbuddy.alignments[0]))
        else:
            output = ""
            for indx, alignment in enumerate(alignbuddy.alignments):
                output += "# Alignment %s\n%s\n\n" % (indx + 1, len(alignment))
            br._stdout("%s\n" % output.strip())
        _exit("num_seqs")

    # Order IDs
    if in_args.order_ids:
        reverse = True if in_args.order_ids[0] else False
        _print_aligments(order_ids(alignbuddy, reverse=reverse))
        _exit("order_ids")

    # Percent IDs
    if in_args.percent_id:
        alignbuddy = percent_id(alignbuddy)
        output = ""
        for indx, alignment in enumerate(alignbuddy.alignments):
            output += "### Alignment %s ###\n\t" % (indx + 1)
            for id1 in alignment.percent_ids:
                output += "%s\t" % id1

            output = output.strip() + "\n"
            for id1 in alignment.percent_ids:
                output += "%s\t" % id1
                for id2 in alignment.percent_ids:
                    output += "\t" if id1 == id2 else "%s\t" % round(alignment.percent_ids[id1][id2], 3)
                output = output.strip() + "\n"
            output = output.strip() + "\n\n"
        output = output.strip() + "\n"
        print(output)
        _exit("percent_ids")

    # Position Frequency Matrix
    if in_args.pos_freq_mat:
        position_frequency_matrix(alignbuddy)
        output = ""
        for indx, alignment in enumerate(alignbuddy.alignments):
            output += "### Alignment %s ###\n" % (indx + 1)
            matrix = [["{0:.3f}".format(pos[res]) for pos in alignment.pfm] for res in alignment.pfm[0]]
            matrix = ["\t".join(res) for res in matrix]
            counter = 0
            for res in alignment.pfm[0]:
                output += "%s\t%s\n" % (res, matrix[counter])
                counter += 1
            output += "\n"
        br._stdout(output)
        _exit("pos_freq_mat")

    # Pull records
    if in_args.pull_records:
        description = False
        args = in_args.pull_records[0]
        for indx, arg in enumerate(args):
            if arg == "full":
                description = True
                del args[indx]
                break
        args = br.clean_regex(args, in_args.quiet)
        if args:
            alignbuddy = pull_records(alignbuddy, args, description)
        _print_aligments(alignbuddy)
        _exit("pull_records")

    # Rename IDs
    if in_args.rename_ids:
        args = in_args.rename_ids[0]
        if not br.clean_regex(args[0], in_args.quiet):
            _print_aligments(alignbuddy)  # Exit gracefully if regex is malformed
            _exit("rename_ids")

        if len(args) not in [2, 3]:
            _raise_error(AttributeError("rename_ids requires two or three argments: "
                                        "query, replacement, [max replacements]"), "rename_ids")
        num = 0
        try:
            num = num if len(args) == 2 else int(args[2])
        except ValueError:
            _raise_error(ValueError("Max replacements argument must be an integer"), "rename_ids")

        try:
            _print_aligments(rename(alignbuddy, args[0], args[1], num))
        except AttributeError as e:
            _raise_error(e, "rename_ids", "There are more replacement")
        _exit("rename_ids")

    # Reverse Transcribe
    if in_args.reverse_transcribe:
        try:
            _print_aligments(rna2dna(alignbuddy))
        except TypeError as e:
            _raise_error(e, "reverse_transcribe", "RNA sequence required, not")
        _exit("reverse_transcribe")

    # Screw formats
    if in_args.screw_formats:
        try:
            alignbuddy.set_format(in_args.screw_formats)
        except TypeError as e:
            _raise_error(e, "screw_formats", "Format type '%s' is not recognized/supported" % in_args.screw_formats)
        if in_args.in_place:  # Need to change the file extension
            _path, ext = os.path.splitext(os.path.abspath(in_args.alignments[0]))
            _path = "%s.%s" % (_path, br.format_to_extension[alignbuddy.out_format])

            os.remove(in_args.alignments[0])
            in_args.alignments[0] = _path
            open(in_args.alignments[0], "w", encoding="utf-8").close()
        _print_aligments(alignbuddy)
        _exit("screw_formats")

    # Split alignments into files
    if in_args.split_to_files:
        if len(alignbuddy.alignments) == 1:
            _raise_error(ValueError("Only one alignment present, nothing written."), "split_to_files")

        args = in_args.split_to_files[0]
        if len(args) > 2:
            br._stderr("Warning: Only one prefix can be accepted, %s where provided. Using the first.\n" % (len(args) - 1))
        in_args.in_place = True

        out_dir = os.path.abspath(args[0])
        prefix = "Alignment_" if len(args) == 1 else args[1]

        os.makedirs(out_dir, exist_ok=True)
        check_quiet = in_args.quiet  # 'quiet' must be toggled to 'on' for _print_recs()
        in_args.quiet = True
        alignments = alignbuddy.alignments
        padding = '{:0>%sd}' % int(ceil(log(len(alignments), 10)))  # Allows sequential numbering, padded with zeros
        for indx, alignment in enumerate(alignments):
            alignbuddy.alignments = [alignment]
            ext = br.format_to_extension[alignbuddy.out_format]
            in_args.alignments[0] = "%s%s%s%s.%s" % (out_dir, os.sep, prefix, padding.format(indx + 1), ext)
            br._stderr("New file: %s\n" % in_args.alignments[0], check_quiet)
            open(in_args.alignments[0], "w", encoding="utf-8").close()
            _print_aligments(alignbuddy)
        _exit("split_to_files")

    # Translate CDS
    if in_args.translate:
        try:
            _print_aligments(translate_cds(alignbuddy))
        except TypeError as e:
            _raise_error(e, "translate", "Nucleic acid sequence required, not protein.")
        _exit("translate")

    # Trimal
    if in_args.trimal:
        args = "gappyout" if not in_args.trimal[0] else in_args.trimal[0]
        try:
            args = abs(float(args))
        except ValueError as e:
            if "could not convert string to float" in str(e):
                pass
            else:
                _raise_error(e, "trimal")
        try:
            _print_aligments(trimal(alignbuddy, args))
        except NotImplementedError as e:
            _raise_error(e, "trimal", "not an implemented trimal method")
        _exit("trimal")

    # Transcribe
    if in_args.transcribe:
        try:
            _print_aligments(dna2rna(alignbuddy))
        except TypeError as e:
            _raise_error(e, "transcribe", "DNA sequence required, not")
        _exit("transcribe")

    # Uppercase
    if in_args.uppercase:
        _print_aligments(uppercase(alignbuddy))
        _exit("uppercase")


def main():
    br.preparse_flags()
    initiation = []
    try:
        initiation = argparse_init()  # initiation = [in_agrs, alignbuddy]
        command_line_ui(*initiation)
    except (KeyboardInterrupt, br.GuessError) as _e:
        print(_e)
        return False
    except SystemExit:
        return False
    except Exception as _e:
        function = ""
        for next_arg in vars(initiation[0]):
            if getattr(initiation[0], next_arg) and next_arg in br.alb_flags:
                function = next_arg
                break
        br.send_traceback("AlignBuddy", function, _e, VERSION)
        return False
    return True


if __name__ == '__main__':
    main()
