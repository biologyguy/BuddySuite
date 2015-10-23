#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: AlignBuddy.py
version: 1, alpha
author: Stephen R. Bond
email: steve.bond@nih.gov
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
# Standard library
import sys
import os
from copy import deepcopy
from io import StringIO, TextIOWrapper
import random
import re
from collections import OrderedDict
from shutil import *
from subprocess import Popen, PIPE, CalledProcessError

# Third party
sys.path.insert(0, "./")  # For stand alone executable, where dependencies are packaged with BuddySuite
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError

# BuddySuite specific
import buddy_resources as br
import MyFuncs

# ##################################################### WISH LIST #################################################### #
# - Map features from a sequence file over to the alignment
# - Back-translate


# #################################################### CHANGE LOG #################################################### #
# ##################################################### GLOBALS ###################################################### #
GAP_CHARS = ["-", ".", " "]
VERSION = br.Version("AlignBuddy", 1, 'alpha', br.contributors)
OUTPUT_FORMATS = ["clustal", "embl", "fasta", "genbank", "gb", "nexus", "stockholm",
                  "phylip", "phylipis", "phylip-strict", "phylip-interleaved-strict",
                  "phylipi", "phylip-relaxed", "phylip-interleaved", "phylipr",
                  "phylips", "phylipsr", "phylip-sequential", "phylip-sequential-relaxed",
                  "phylipss", "phylip-sequential-strict"]


# #################################################### ALIGNBUDDY #################################################### #
class AlignBuddy:  # Open a file or read a handle and parse, or convert raw into a Seq object
    def __init__(self, _input, in_format=None, out_format=None):
        # ####  IN AND OUT FORMATS  #### #
        # Holders for input type. Used for some error handling below
        in_handle = None
        raw_seq = None
        in_file = None

        # Handles
        if str(type(_input)) == "<class '_io.TextIOWrapper'>":
            if not _input.seekable():  # Deal with input streams (e.g., stdout pipes)
                temp = StringIO(_input.read())
                _input = temp
            _input.seek(0)
            in_handle = _input.read()
            _input.seek(0)

        # Plain text in a specific format
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

        self._in_format = parse_format(in_format) if in_format else guess_format(_input)
        if not self._in_format:
            if in_file:
                raise br.GuessError("Could not determine format from _input file '{0}'.\n"
                                    "Try explicitly setting with -f flag.".format(in_file))
            elif raw_seq:
                raise br.GuessError("Could not determine format from raw input\n{0} ..."
                                    "Try explicitly setting with -f flag.".format(raw_seq)[:50])
            elif in_handle:
                raise br.GuessError("Could not determine format from input file-like object\n{0} ..."
                                    "Try explicitly setting with -f flag.".format(in_handle)[:50])
            else:  # This should be unreachable.
                raise br.GuessError("Unable to determine format or input type. "
                                    "Please check how AlignBuddy is being called.")

        self._out_format = self._in_format if not out_format else parse_format(out_format)
        if self._out_format not in OUTPUT_FORMATS:
            raise(TypeError("Output type '%s' is not recognized/supported" % self._out_format))

        # ####  ALIGNMENTS  #### #
        if type(_input) == AlignBuddy:
            alignments = _input.alignments

        elif isinstance(_input, list):
            # make sure that the list is actually MultipleSeqAlignment objects
            sample = _input if len(_input) < 5 else random.sample(_input, 5)
            for _seq in sample:
                if type(_seq) != MultipleSeqAlignment:
                    raise TypeError("Seqlist is not populated with SeqRecords.")
            alignments = _input

        elif str(type(_input)) == "<class '_io.TextIOWrapper'>" or isinstance(_input, StringIO):
            alignments = list(AlignIO.parse(_input, self._in_format))

        elif os.path.isfile(_input):
            with open(_input, "r") as _input:
                alignments = list(AlignIO.parse(_input, self._in_format))
        else:  # May be unreachable
            alignments = None

        self.alpha = guess_alphabet(alignments)
        for alignment in alignments:
            alignment._alphabet = self.alpha
            for rec in alignment:
                rec.seq.alphabet = self.alpha
        self.alignments = alignments

    def set_format(self, in_format):
        self._out_format = parse_format(in_format)

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

    def print(self):
        print(self)
        return

    def __str__(self):
        if len(self.alignments) == 0:
            return "AlignBuddy object contains no alignments.\n"

        self._out_format = self._out_format.lower()
        multiple_alignments_unsupported = ["fasta", "gb", "genbank", "nexus"]
        if self._out_format in multiple_alignments_unsupported and len(self.alignments) > 1:
            raise ValueError("%s format does not support multiple alignments in one file.\n" % self._out_format)

        if self._out_format == "phylipsr":
            output = _phylip_sequential_out(self)

        elif self._out_format == "phylipss":
            output = _phylip_sequential_out(self, relaxed=False)

        else:
            tmp_dir = MyFuncs.TempDir()
            with open("%s/aligns.tmp" % tmp_dir.path, "w") as ofile:
                try:
                    AlignIO.write(self.alignments, ofile, self._out_format)
                except ValueError as e:
                    if "Repeated name " in str(e):
                        raise br.PhylipError(str(e))
            with open("%s/aligns.tmp" % tmp_dir.path, "r") as ifile:
                output = ifile.read()

        return output

    def write(self, file_path):
        with open(file_path, "w") as ofile:
            ofile.write("{0}\n".format(str(self).rstrip()))
        return


# ################################################# HELPER FUNCTIONS ################################################# #
class GuessError(Exception):
    """Raised when input format cannot be guessed"""
    def __init__(self, _value):
        self.value = _value

    def __str__(self):
        return self.value


def guess_alphabet(alignbuddy):
    align_list = alignbuddy if isinstance(alignbuddy, list) else alignbuddy.alignments
    seq_list = []
    for alignment in align_list:
        seq_list += [str(x.seq) for x in alignment]

    sequence = "".join(seq_list).upper()
    sequence = re.sub("[NX\-?]", "", sequence)

    if len(sequence) == 0:
        return None

    if 'U' in sequence:  # U is unique to RNA
        return IUPAC.ambiguous_rna

    percent_dna = len(re.findall("[ATCG]", sequence)) / float(len(sequence))
    if percent_dna > 0.85:  # odds that a sequence with no Us and such a high ATCG count be anything but DNA is low
        return IUPAC.ambiguous_dna
    else:
        return IUPAC.protein


def guess_format(_input):  # _input can be list, SeqBuddy object, file handle, or file path.
    # If input is just a list, there is no BioPython in-format. Default to stockholm.
    if isinstance(_input, list):
        return "stockholm"

    # Pull value directly from object if appropriate
    if type(_input) == AlignBuddy:
        return _input._in_format

    # If input is a handle or path, try to read the file in each format, and assume success if not error and # seqs > 0
    if os.path.isfile(str(_input)):
        _input = open(_input, "r")

    if str(type(_input)) == "<class '_io.TextIOWrapper'>" or isinstance(_input, StringIO):
        # Die if file is empty
        if _input.read() == "":
            sys.exit("Input file is empty.")
        _input.seek(0)

        possible_formats = ["gb", "phylip-relaxed", "phylips", "stockholm", "fasta", "nexus", "clustal"]
        for _format in possible_formats:
            try:
                _input.seek(0)
                if _format == "phylips":
                    if _phylip_sequential_read(_input.read()):
                        _input.seek(0)
                        return parse_format(_format, "in")
                    else:
                        continue
                if list(AlignIO.parse(_input, _format)):
                    _input.seek(0)
                    return parse_format(_format, "in")
                else:
                    continue
            except StopIteration:
                continue
            except br.PhylipError:
                continue
            except ValueError as e:
                if "Found a record of length" in str(e):
                    raise ValueError("Malformed Phylip --> %s" % e)
                continue

        return None  # Unable to determine format from file handle

    else:
        raise GuessError("Unsupported _input argument in guess_format(). %s" % _input)


def parse_format(_format, mode="out"):
    _format = _format.lower()
    if _format in ["phylip", "phylipis", "phylip-strict", "phylip-interleaved-strict"]:
        return "phylip"

    if _format in ["phylipi", "phylip-relaxed", "phylip-interleaved", "phylipr"]:
        return "phylip-relaxed"

    if _format in ["phylips", "phylipsr", "phylip-sequential", "phylip-sequential-relaxed"]:
        return "phylipsr" if mode == "out" else "phylip-sequential"

    if _format in ["phylipss", "phylip-sequential-strict"]:
        return "phylipss" if mode == "out" else "phylip-sequential"

    return _format


def make_copy(alignbuddy):
    alphabet_list = [rec.seq.alphabet for rec in alignbuddy.records()]
    _copy = deepcopy(alignbuddy)
    _copy.alpha = alignbuddy.alpha
    for indx, rec in enumerate(_copy.records()):
        rec.seq.alphabet = alphabet_list[indx]
    return _copy


def _phylip_sequential_out(alignbuddy, relaxed=True):
    output = ""
    ids = []
    for alignment in alignbuddy.alignments:
        id_check = []
        for rec in alignment:
            if rec.id in id_check:
                raise br.PhylipError("Malformed Phylip --> Repeat id '%s'" % rec.id)
            id_check.append(rec.id)

        max_seq_length = 0
        for rec in alignment:
            max_seq_length = len(rec.seq) if len(rec.seq) > max_seq_length else max_seq_length

        output += " %s %s" % (len(alignment), max_seq_length)
        for rec in alignment:
            if relaxed:
                seq_id = re.sub(' \t', '_', rec.id)
                output += "\n%s \n%s" % (seq_id, rec.seq)
            else:
                seq_id = rec.id[:10].ljust(11)
                output += "\n%s%s" % (seq_id, rec.seq)

            if seq_id in ids:
                raise br.PhylipError("Malformed Phylip --> Repeat id '%s' after strict truncation. "
                                     "Try a relaxed Phylip format (phylipr or phylipsr)." % seq_id)
            ids.append(seq_id)

        output += "\n"
    return output


def _phylip_sequential_read(sequence):
    sequence = "\n%s" % sequence
    alignments = re.split("\n ([0-9]+) ([0-9]+)\n", sequence)[1:]
    align_dict = OrderedDict()
    for indx in range(int(len(alignments) / 3)):
        align_dict[(alignments[indx * 3], alignments[indx * 3 + 1])] = alignments[indx * 3 + 2]

    temp_file = MyFuncs.TempFile()
    aligns = []
    for key, seqs in align_dict.items():
        seqs = seqs.strip().split("\n")
        if int(key[0]) != int(len(seqs) / 2):
            raise br.PhylipError("Malformed Phylip --> %s sequences expected, %s found." % (key[0], int(len(seqs) / 2)))

        output = ""
        for indx in range(int(len(seqs) / 2)):
            if int(key[1]) != int(len(seqs[indx * 2 + 1])):
                raise br.PhylipError("Malformed Phylip --> Sequence %s has %s columns, %s expected." %
                                     (seqs[indx * 2], len(seqs[indx * 2 + 1]), key[1]))
            output += ">%s\n%s\n" % (seqs[indx * 2], seqs[indx * 2 + 1])
        temp_file.write(output, "w")
        aligns.append(AlignIO.read(temp_file.get_handle("r"), "fasta"))
        temp_file.close()
    return aligns


def _stderr(message, quiet=False):
    if not quiet:
        sys.stderr.write(message)
    return


def _stdout(message, quiet=False):
    if not quiet:
        sys.stdout.write(message)
    return


# ################################################ MAIN API FUNCTIONS ################################################ #
def alignment_lengths(alignbuddy):
    """
    Returns a list of alignment lengths
    :param alignbuddy: The AlignBuddy object to be analyzed
    :return: A list of alignment lengths
    """
    output = []
    for alignment in alignbuddy.alignments:
        output.append(alignment.get_alignment_length())
    return output


def clean_seq(alignbuddy, skip_list=None):
    """
    Remove all non-sequence charcters from sequence strings
    :param alignbuddy: The AlignBuddy object to be cleaned
    :param skip_list: A list of characters to be left alone
    :return: The cleaned AlignBuddy object
    """
    skip_list = "" if not skip_list else "".join(skip_list)
    for rec in alignbuddy.records_iter():
        if alignbuddy.alpha == IUPAC.protein:
            full_skip = "ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy%s" % skip_list
            rec.seq = Seq(re.sub("[^%s]" % full_skip, "", str(rec.seq)),
                          alphabet=alignbuddy.alpha)
        else:
            full_skip = "ATGCUatgcu%s" % skip_list
            rec.seq = Seq(re.sub("[^%s]" % full_skip, "", str(rec.seq)), alphabet=alignbuddy.alpha)

    return alignbuddy


def codon_alignment(alignbuddy):
    """
    Organizes nucleotide alignment into codon triplets
    :param alignbuddy: The AlignBuddy object to be rearranged
    :return: The rearranged AlignBuddy object
    """
    if alignbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")

    for rec in alignbuddy.records_iter():
        if rec.seq.alphabet == IUPAC.protein:
            raise TypeError("Error: Record %s is protein. Nucleic acid sequence required." % rec.name)

        seq_string = str(rec.seq)
        output = seq_string[0]
        held_residues = ""
        position = 2
        for residue in seq_string[1:]:
            if position == 1:
                while len(held_residues) >= 3:
                    output += held_residues[:3]
                    held_residues = held_residues[3:]
                position = len(held_residues) + 1
                output += held_residues
                held_residues = ""

            if position == 1:
                output += residue

            elif output[-1] not in GAP_CHARS and residue not in GAP_CHARS:
                output += residue

            elif output[-1] in GAP_CHARS and residue in GAP_CHARS:
                output += residue

            else:
                held_residues += residue
                continue

            if position != 3:
                position += 1
            else:
                position = 1

        output += held_residues
        rec.seq = Seq(output, alphabet=rec.seq.alphabet)

    return alignbuddy


def concat_alignments(alignbuddy, pattern):
    """
    Concatenates two or more alignments together, end-to-end
    :param alignbuddy: The AlignBuddy object whose alignments will be concatenated
    :param pattern: The pattern to split the sequence names with
    :return: The AlignBuddy object containing a concatenated alignment
    """
    # collapsed multiple genes from single taxa down to one consensus seq
    # detected mixed sequence types
    if len(alignbuddy.alignments) < 2:
        raise AttributeError("Please provide at least two alignments.")

    def organism_list():
        orgnsms = set()
        for align in alignbuddy.alignments:
            for _record in align:
                orgnsms.add(_record.id)
        return list(orgnsms)

    def add_blanks(_record, _num):
        for _ in range(_num):
            _record.seq += '-'

    missing_organisms = organism_list()

    # dict[alignment_index][organism] -> gene_name
    sequence_names = OrderedDict()
    for al_indx, alignment in enumerate(alignbuddy.alignments):
        sequence_names[al_indx] = OrderedDict()
        for record in alignment:
            organism = re.split(pattern, record.id)[0]
            gene = ''.join(re.split(pattern, record.id)[1:])
            if organism in sequence_names[al_indx].keys():
                sequence_names[al_indx][organism].append(gene)
            else:
                sequence_names[al_indx][organism] = [gene]
            record.id = organism
    for al_indx, alignment in enumerate(alignbuddy.alignments):
        duplicate_table = OrderedDict()
        for record in alignment:
            if record.id in duplicate_table.keys():
                duplicate_table[record.id].append(record)
            else:
                duplicate_table[record.id] = [record]
        temp = []
        for record in alignment:
            if len(duplicate_table[record.id]) == 1:
                temp.append(record)
                duplicate_table.pop(record.id)
        alignbuddy.alignments[al_indx] = MultipleSeqAlignment(temp, alphabet=alignbuddy.alpha)
        for gene in duplicate_table:
            consensus = SummaryInfo(MultipleSeqAlignment(duplicate_table[gene]))
            consensus = consensus.gap_consensus(consensus_alpha=alignbuddy.alpha)
            consensus = SeqRecord(seq=consensus, id=gene)
            alignbuddy.alignments[al_indx].append(consensus)

    for al_indx in range(len(alignbuddy.alignments)):
        for _id in missing_organisms:
            organism = re.split(pattern, _id)[0]
            if organism not in sequence_names[al_indx].keys():
                sequence_names[al_indx][organism] = ['missing']

    for x in range(len(missing_organisms)):
        missing_organisms[x] = re.split(pattern, missing_organisms[x])[0]

    base_alignment = deepcopy(alignbuddy.alignments[0])
    for record in base_alignment:
        while record.id in missing_organisms:
            missing_organisms.remove(record.id)
    for organism in missing_organisms:
        new_record = SeqRecord(Seq('', alphabet=alignbuddy.alpha), id=re.split(pattern, organism)[0])
        add_blanks(new_record, base_alignment.get_alignment_length())
        base_alignment.append(new_record)

    for record in base_alignment:
        record.description = 'concatenated_alignments'
        record.name = 'multiple'

    curr_length = 0
    for al_indx, alignment in enumerate(alignbuddy.alignments):
        for base_indx, base_rec in enumerate(base_alignment):
            added = False
            for record in alignment:
                if base_rec.id == record.id:
                    if al_indx > 0:
                        base_rec.seq += record.seq
                        added = True
            if not added and al_indx > 0:
                add_blanks(base_rec, alignment.get_alignment_length())
            feature_location = FeatureLocation(start=curr_length,
                                               end=curr_length + alignment.get_alignment_length())
            feature = SeqFeature(location=feature_location, type='alignment' + str(al_indx + 1),
                                 qualifiers={'name': sequence_names[al_indx][base_rec.id]})
            base_alignment[base_indx].features.append(feature)
        curr_length += alignment.get_alignment_length()

    alignbuddy.alignments = [base_alignment]
    return alignbuddy


def consensus_sequence(_alignbuddy):
    """
    Generates a rough consensus sequence from an alignment
    :param _alignbuddy: The AlignBuddy object to be processed
    :return: A list of SeqRecords containing the consensus sequences
    """
    output = []
    for alignment in _alignbuddy.alignments:
        aln_id = alignment[0].id
        aln = SummaryInfo(alignment)
        print(aln.gap_consensus())
        output.append(SeqRecord(seq=aln.gap_consensus(), id=aln_id))
    return output


def delete_rows(alignbuddy, search):
    """
    Deletes rows with names/IDs matching a search pattern
    :param alignbuddy: The AlignBuddy object to be modified
    :param search: The regex pattern to search with
    :return: The modified AlignBuddy object
    """
    alignments = []
    for alignment in alignbuddy.alignments:
        matches = []
        for record in alignment:
            if not re.search(search, record.id) and not re.search(search, record.description) \
                    and not re.search(search, record.name):
                matches.append(record)
        alignments.append(MultipleSeqAlignment(matches))
    alignbuddy.alignments = alignments
    trimal(alignbuddy, "clean")
    return alignbuddy


def dna2rna(alignbuddy):
    """
    Back-transcribes DNA into RNA sequences
    :param alignbuddy: The AlignBuddy object to be back-transcribed
    :return: The back-transcribed AlignBuddy object
    """
    if alignbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for alignment in alignbuddy.alignments:
        for rec in alignment:
            rec.seq = Seq(str(rec.seq.transcribe()), alphabet=IUPAC.ambiguous_rna)
    alignbuddy.alpha = IUPAC.ambiguous_rna
    return alignbuddy


def extract_range(alignbuddy, start, end):
    """
    Extracts all columns within a given range
    :param alignbuddy: The AlignBuddy object to be extracted from
    :param start: The starting residue (indexed from 1)
    :param end: The end residue (inclusive)
    :return: The extracted AlignBuddy object
    """
    start = 1 if int(start) < 1 else start
    # Don't use the standard index-starts-at-0... _end must be left for the range to be inclusive
    start, end = int(start) - 1, int(end)
    if end < start:
        raise ValueError("Error at extract range: The value given for end of range is smaller than for the start "
                         "of range.")
    for alignment in alignbuddy.alignments:
        for rec in alignment:
            rec.seq = Seq(str(rec.seq)[start:end], alphabet=rec.seq.alphabet)
            rec.description += " Sub-sequence extraction, from residue %s to %s" % (start + 1, end)
            features = []
            for feature in rec.features:
                if feature.location.end < start:
                    continue
                if feature.location.start > end:
                    continue

                feat_start = feature.location.start - start
                if feat_start < 0:
                    feat_start = 0

                feat_end = feature.location.end - start
                if feat_end > len(str(rec.seq)):
                    feat_end = len(str(rec.seq))

                new_location = FeatureLocation(feat_start, feat_end)
                feature.location = new_location
                features.append(feature)
            rec.features = features
    return alignbuddy


# ToDo: clustalomega may be clustalo
# ToDo: Completely refactor the handling of output formats
def generate_msa(seqbuddy, tool, params=None, keep_temp=None, quiet=False):
    """
    Calls sequence aligning tools to generate multiple sequence alignments
    :param seqbuddy: The SeqBuddy object containing the sequences to be aligned
    :param tool: The alignment tool to be used (pagan/prank/muscle/clustalw2/clustalomega/mafft)
    :param params: Additional parameters to be passed to the alignment tool
    :param keep_temp: Determines if/where the temporary files will be kept
    :return: An AlignBuddy object containing the alignment produced.
    """
    if params is None:
        params = ''
    tool = tool.lower()

    if keep_temp:
        if os.path.exists(keep_temp):
            _stderr("Warning: {0} already exists. Please specify a different path.\n".format(keep_temp), quiet)
            sys.exit()

    tool_urls = {'mafft': 'http://mafft.cbrc.jp/alignment/software/',
                 'prank': 'http://wasabiapp.org/software/prank/prank_installation/',
                 'pagan': 'http://wasabiapp.org/software/pagan/pagan_installation/',
                 'muscle': 'http://www.drive5.com/muscle/downloads.htm',
                 'clustalw': 'http://www.clustal.org/clustal2/#Download',
                 'clustalw2': 'http://www.clustal.org/clustal2/#Download',
                 'clustalomega': 'http://www.clustal.org/omega/#Download'}

    if tool not in tool_urls:
        raise AttributeError("{0} is not a supported alignment tool.".format(tool))
    if which(tool) is None:
        _stderr('#### Could not find {0} in $PATH. ####\n'.format(tool), quiet)
        _stderr('Please go to {0} to install {1}.\n'.format(tool_urls[tool], tool))
        sys.exit()
    else:
        tmp_dir = MyFuncs.TempDir()
        tmp_in = "{0}/tmp.fa".format(tmp_dir.path)

        try:
            from SeqBuddy import hash_sequence_ids
        except ImportError:
            _stderr("SeqBuddy is needed to use generate_msa(). Please install it and try again.")
            sys.exit()

        params = re.split(' ', params, )
        for indx, token in enumerate(params):
            if os.path.exists(token):
                params[indx] = os.path.abspath(token)
        params = ' '.join(params)

        hash_sequence_ids(seqbuddy, 8)

        output = ''

        seqbuddy.out_format = 'fasta'
        with open("{0}/tmp.fa".format(tmp_dir.path), 'w') as out_file:
            out_file.write(str(seqbuddy))
        if tool == 'clustalomega':
            command = '{0} {1} -i {2}'.format(tool, params, tmp_in)
        elif tool == 'clustalw2':
            command = '{0} -infile={1} {2} -outfile={3}/result'.format(tool, tmp_in, params, tmp_dir.path)
        elif tool == 'muscle':
            command = '{0} -in {1} {2}'.format(tool, tmp_in, params)
        elif tool == 'prank':
            command = '{0} -d={1} {2} -o={3}/result'.format(tool, tmp_in, params, tmp_dir.path)
        elif tool == 'pagan':
            command = '{0} -s {1} {2} -o {3}/result'.format(tool, tmp_in, params, tmp_dir.path)
        else:
            command = '{0} {1} {2}'.format(tool, params, tmp_in)

        try:
            if tool in ['prank', 'pagan']:
                Popen(command, shell=True, universal_newlines=True, stdout=sys.stderr).wait()
            else:
                output = Popen(command, shell=True, stdout=PIPE).communicate()
                output = output[0].decode()
        except CalledProcessError:
            _stderr('\n#### {0} threw an error. Scroll up for more info. ####\n\n'.format(tool), quiet)
            sys.exit()

        if tool.startswith('clustalw'):
            with open('{0}/result'.format(tmp_dir.path)) as result:
                output = result.read()
        elif tool == 'prank':
            extension = 'fas'
            if '-f=nexus' in params or '-f=nexus' in params:
                extension = 'nex'
            elif '-f=phylipi' in params or '-f=phylips' in params:
                extension = 'phy'
            possible_files = os.listdir(tmp_dir.path)
            filename = 'result.best.{0}'.format(extension)
            for _file in possible_files:
                if 'result.best' in _file and extension in _file:
                    filename = _file
            with open('{0}/{1}'.format(tmp_dir.path, filename)) as result:
                output = result.read()
        elif tool == 'pagan':
            extension = 'fas'
            if '-f nexus' in params:
                extension = 'nex'
            elif '-f phylipi' in params or '-f phylips' in params:
                extension = 'phy'
            with open('{0}/result.{1}'.format(tmp_dir.path, extension)) as result:
                output = result.read()

        # Fix broken outputs to play nicely with AlignBuddy parsers
        if (tool == 'mafft' and '--clustalout' in params) or \
                (tool.startswith('clustalw2') and '-output' not in params) or \
                (tool == 'clustalomega' and ('clustal' in params or '--outfmt clu' in params or
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
        if tool.startswith('clustalw') and 'nexus' in params or '=NEXUS' in params:  # if taxa has gap chars
            contents = ''
            for indx, line in enumerate(output.splitlines(keepends=True)):
                if indx > 6 and not line.startswith('\n') and ';' not in line:
                    split = re.split(' ', line, maxsplit=1)
                    contents += "'" + split[0] + "' " + split[1]
                else:
                    contents += line
            output = contents
        alignbuddy = AlignBuddy(output)

        for _hash, sb_rec in seqbuddy.hash_map.items():
            rename(alignbuddy, _hash, sb_rec)
            for alignment in alignbuddy.alignments:
                for rec in alignment:
                    if _hash in rec.annotations:
                        rec.annotations.pop(_hash)
                    to_pop = []
                    for key in rec.annotations:
                        if _hash in rec.annotations[key]:
                            to_pop.append(key)
                    for key in to_pop:
                        rec.annotations.pop(key)
                    rec.name = re.sub(_hash, '', rec.name)
                    rec.description = re.sub(_hash, '', rec.description)
        _stderr('\n')

        if keep_temp:
            try:
                copytree(tmp_dir.path, keep_temp)
            except FileExistsError:
                # Should never get here
                pass

        _stderr("Returning to AlignBuddy...\n\n", quiet)

        return alignbuddy


def list_ids(alignbuddy):
    """
    Returns a list of lists of sequence IDs
    :param alignbuddy: The AlignBuddy object to be analyzed
    :return: A list of lists of sequence IDs
    """
    output = []
    for al_indx, _alignment in enumerate(alignbuddy.alignments):
        output.append([])
        for _rec in _alignment:
            output[al_indx].append(_rec.id)
    return output


def lowercase(alignbuddy):
    """
    Converts all sequence residues to lowercase.
    :param alignbuddy: The AlignBuddy object to be modified.
    :return: The modified AlignBuddy object
    """
    for rec in alignbuddy.records_iter():
        rec.seq = Seq(str(rec.seq).lower(), alphabet=rec.seq.alphabet)
    return alignbuddy


def num_seqs(alignbuddy):
    """
    Returns a list of alignment lengths
    :param alignbuddy:
    :return: A list of alignment lengths
    """
    return [len(_alignment) for _alignment in alignbuddy.alignments]


def order_ids(alignbuddy, reverse=False):
    """
    Sorts the alignments by ID, alphabetically
    :param alignbuddy: The AlignBuddy object to be sorted
    :param reverse: Reverses the order
    :return: A sorted AlignBuddy object
    """
    for al_indx, alignment in enumerate(alignbuddy.alignments):
        output = [(rec.id, rec) for rec in alignment]
        output = sorted(output, key=lambda x: x[0], reverse=reverse)
        output = [_rec[1] for _rec in output]
        alignbuddy.alignments[al_indx] = MultipleSeqAlignment(output)
    return alignbuddy


def pull_rows(alignbuddy, search):
    """
    Retrieves rows with names/IDs matching a search pattern
    :param alignbuddy: The AlignBuddy object to be pulled from
    :param search: The regex pattern to search with
    :return: The modified AlignBuddy object
    """
    alignments = []
    for alignment in alignbuddy.alignments:
        matches = []
        for record in alignment:
            if re.search(search, record.id) or re.search(search, record.description) \
                    or re.search(search, record.name):
                matches.append(record)
        alignments.append(MultipleSeqAlignment(matches))
    alignbuddy.alignments = alignments
    trimal(alignbuddy, "clean")
    return alignbuddy


def rename(alignbuddy, query, replace="", num=0):  # TODO Allow a replacement pattern increment (like numbers)
    """
    Rename an alignment's sequence IDs
    :param alignbuddy: The AlignBuddy object to be modified
    :param query: The pattern to be searched for
    :param replace: The string to be substituted
    :param num: The maximum number of substitutions to make
    :return: The modified AlignBuddy object
    """
    for alignment in alignbuddy.alignments:
        for rec in alignment:
            new_name = re.sub(query, replace, rec.id, num)
            rec.id = new_name
            rec.name = new_name
    return alignbuddy


def rna2dna(alignbuddy):
    """
    Transcribes RNA into DNA sequences.
    :param alignbuddy: The AlignBuddy object to be transcribed
    :return: The transcribed AlignBuddy object
    """
    if alignbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for alignment in alignbuddy.alignments:
        for rec in alignment:
            rec.seq = Seq(str(rec.seq.back_transcribe()), alphabet=IUPAC.ambiguous_dna)
    alignbuddy.alpha = IUPAC.ambiguous_dna
    return alignbuddy


def split_alignbuddy(alignbuddy):
    """
    Splits each alignment into its own AlignBuddy object
    :param alignbuddy: The AlignBuddy object to be split
    :return: A list of AlignBuddy objects
    """
    ab_objs_list = []
    for alignment in alignbuddy.alignments:
        ab = AlignBuddy([alignment])
        ab._in_format = alignbuddy._in_format
        ab.set_format(alignbuddy._out_format)
        ab_objs_list.append(ab)
    return ab_objs_list


def translate_cds(alignbuddy, quiet=False):  # adding 'quiet' will suppress the errors thrown by translate(cds=True)
    """
    Translates a nucleotide alignment into a protein alignment.
    :param alignbuddy: The AlignBuddy object to be translated
    :param quiet: Suppress error messages and warnings
    :return: The translated AlignBuddy object
    """
    if alignbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")

    def trans(in_seq):
        try:
            in_seq.seq = in_seq.seq.translate(cds=True, to_stop=True)
            return in_seq

        except TranslationError as e:
            if not quiet:
                sys.stderr.write("Warning: %s in %s\n" % (e, in_seq.id))
            return e

    def map_gaps(nucl, pep):
        nucl = str(nucl.seq)
        pep_string = str(pep.seq)
        new_seq = ""
        for _codon in [nucl[j:j + 3] for j in range(0, len(nucl), 3)]:
            if _codon[0] not in GAP_CHARS:
                new_seq += pep_string[0]
                pep_string = pep_string[1:]

            else:
                new_seq += "-"
        pep.seq = Seq(new_seq, alphabet=IUPAC.protein)
        return pep

    codon_alignment(alignbuddy)
    copy_alignbuddy = make_copy(alignbuddy)
    clean_seq(copy_alignbuddy, skip_list="RYWSMKHBVDNXrywsmkhbvdnx")
    for align_indx, alignment in enumerate(copy_alignbuddy.alignments):
        for rec_indx, rec in enumerate(alignment):
            rec.features = []
            while True:
                test_trans = trans(deepcopy(rec))
                # success
                if str(type(test_trans)) == "<class 'Bio.SeqRecord.SeqRecord'>":
                    break

                # not standard length
                if re.search("Sequence length [0-9]+ is not a multiple of three", str(test_trans)):
                    orig_rec = alignbuddy.alignments[align_indx][rec_indx]
                    orig_rec_seq = str(orig_rec.seq)
                    for _ in range(len(str(rec.seq)) % 3):
                        orig_rec_seq = re.sub(".([%s]+)$" % "".join(GAP_CHARS), r"\1-", orig_rec_seq)

                    orig_rec.seq = Seq(orig_rec_seq, alphabet=orig_rec.seq.alphabet)

                    rec.seq = Seq(str(rec.seq)[:(len(str(rec.seq)) - len(str(rec.seq)) % 3)],
                                  alphabet=rec.seq.alphabet)
                    continue

                # not a start codon
                if re.search("First codon '[A-Za-z]{3}' is not a start codon", str(test_trans)):
                    rec.seq = Seq("ATG" + str(rec.seq)[3:], alphabet=rec.seq.alphabet)
                    continue

                # not a stop codon
                if re.search("Final codon '[A-Za-z]{3}' is not a stop codon", str(test_trans)):
                    rec.seq = Seq(str(rec.seq) + "TGA", alphabet=rec.seq.alphabet)
                    continue

                # non-standard characters
                if re.search("Codon '[A-Za-z]{3}' is invalid", str(test_trans)):
                    regex = re.findall("Codon '([A-Za-z]{3})' is invalid", str(test_trans))
                    regex = "(?i)%s" % regex[0]
                    rec.seq = Seq(re.sub(regex, "NNN", str(rec.seq), count=1), alphabet=rec.seq.alphabet)
                    continue

                # internal stop codon(s) found
                if re.search("Extra in frame stop codon found", str(test_trans)):
                    for i in range(round(len(str(rec.seq)) / 3) - 1):
                        codon = str(rec.seq)[(i * 3):(i * 3 + 3)]
                        if codon.upper() in ["TGA", "TAG", "TAA"]:
                            stop_removed = str(rec.seq)[:(i * 3)] + "NNN" + str(rec.seq)[(i * 3 + 3):]
                            rec.seq = Seq(stop_removed, alphabet=rec.seq.alphabet)
                    continue

                break  # Safety valve, should be unreachable

            try:
                rec.seq = rec.seq.translate()
                rec.seq.alphabet = IUPAC.protein
                rec = map_gaps(alignbuddy.alignments[align_indx][rec_indx], rec)
                alignbuddy.alignments[align_indx][rec_indx].seq = Seq(str(rec.seq), alphabet=rec.seq.alphabet)

            except TranslationError as e1:  # Should be unreachable
                raise TranslationError("%s failed to translate  --> %s\n" % (rec.id, e1))

    alignbuddy.alpha = IUPAC.protein
    return alignbuddy


# http://trimal.cgenomics.org/_media/manual.b.pdf
# ftp://trimal.cgenomics.org/trimal/
def trimal(alignbuddy, threshold, window_size=1):  # ToDo: This might be broken, test it.
    """
    Trims alignment gaps using algorithms from trimal
    :param alignbuddy: The AlignBuddy object to be trimmed
    :param threshold: The threshold value or trimming algorithm to be used
    :param window_size: The window size
    :return: The trimmed AlignBuddy object
    """
    for alignment_index, alignment in enumerate(alignbuddy.alignments):

        # gap_distr is the number of columns w/ each possible number of gaps; the index is == to number of gaps
        gap_distr = [0 for _ in range(len(alignment) + 1)]
        num_columns = alignment.get_alignment_length()
        each_column = [0 for _ in range(num_columns)]

        max_gaps = 0

        for indx in range(num_columns):
                num_gaps = len(re.findall("-", str(alignment[:, indx])))
                gap_distr[num_gaps] += 1
                each_column[indx] = num_gaps

        def remove_cols(cuts):
            _new_alignment = alignment[:, 0:0]
            for _col in cuts:
                _new_alignment += alignment[:, _col:_col + 1]
            return _new_alignment

        def gappyout():
            for i in gap_distr:
                if i == 0:
                    _max_gaps = i + 1
                else:
                    break

            max_slope = -1
            slopes = [-1 for _ in range(len(alignment) + 1)]
            max_iter = len(alignment) + 1
            active_pointer = 0
            while active_pointer < max_iter:
                for i in gap_distr:
                    if i == 0:
                        active_pointer += 1
                    else:
                        break

                prev_pointer1 = int(active_pointer)
                if active_pointer + 1 >= max_iter:
                    break

                while True:
                    active_pointer += 1
                    if active_pointer + 1 >= max_iter or gap_distr[active_pointer] != 0:
                        break

                prev_pointer2 = int(active_pointer)
                if active_pointer + 1 >= max_iter:
                    break

                while True:
                    active_pointer += 1
                    if active_pointer + 1 >= max_iter or gap_distr[active_pointer] != 0:
                        break

                if active_pointer + 1 >= max_iter:
                    break

                slopes[active_pointer] = (active_pointer - prev_pointer2) / len(alignment)
                slopes[active_pointer] /= (gap_distr[active_pointer] + gap_distr[prev_pointer2]) / num_columns

                if slopes[prev_pointer1] != -1:
                    if slopes[active_pointer] / slopes[prev_pointer1] > max_slope:
                        max_slope = slopes[active_pointer] / slopes[prev_pointer1]
                        _max_gaps = prev_pointer1
                elif slopes[prev_pointer2] != -1:
                    if slopes[active_pointer] / slopes[prev_pointer2] > max_slope:
                        max_slope = slopes[active_pointer] / slopes[prev_pointer2]
                        _max_gaps = prev_pointer1

                active_pointer = prev_pointer2

            def mark_cuts():
                cuts = []
                for col, gaps in enumerate(each_column):
                    if gaps <= _max_gaps:
                        cuts.append(col)
                return cuts

            return remove_cols(mark_cuts())

        if threshold in ["no_gaps", "all"]:
            threshold = 0
            new_alignment = alignment[:, 0:0]
            for next_col, num_gaps in enumerate(each_column):
                if num_gaps <= max_gaps:
                    new_alignment += alignment[:, next_col:next_col + 1]

        if threshold == "clean":
            max_gaps = len(alignment) - 1
            new_alignment = alignment[:, 0:0]
            for next_col, num_gaps in enumerate(each_column):
                if num_gaps <= max_gaps:
                    new_alignment += alignment[:, next_col:next_col + 1]

        elif threshold == "gappyout":
            new_alignment = gappyout()

        elif threshold == "strict":  # ToDo: Implement
            new_alignment = False
            pass
        elif threshold == "strictplus":  # ToDo: Implement
            new_alignment = False
            pass
        else:
            try:
                if threshold == 1:
                    _stderr("Warning: Ambiguous threshold of '1'. Assuming 100%, use 0.01 for 1%")

                threshold = float(threshold)
                threshold = 0.0001 if threshold == 0 else threshold
                threshold = (threshold / 100) if threshold > 1 else threshold  # Allows percent or fraction

                max_gaps = round(len(alignment) * threshold)

                new_alignment = alignment[:, 0:0]
                for next_col, num_gaps in enumerate(each_column):
                    if num_gaps <= max_gaps:
                        new_alignment += alignment[:, next_col:next_col + 1]

            except ValueError:
                raise ValueError("Unable to understand the threshold parameter provided -> %s)" % threshold)

        if new_alignment:
            alignbuddy.alignments[alignment_index] = new_alignment
        else:
            raise NotImplementedError("%s has not been implemented" % threshold)
    return alignbuddy


def uppercase(alignbuddy):
    """
    Converts all sequence residues to uppercase.
    :param alignbuddy: The AlignBuddy object to be modified.
    :return: The modified AlignBuddy object
    """
    for rec in alignbuddy.records_iter():
        rec.seq = Seq(str(rec.seq).upper(), alphabet=rec.seq.alphabet)
    return alignbuddy


# ################################################# COMMAND LINE UI ################################################## #
def argparse_init():
    # Catching params to prevent weird collisions with alignment program arguments
    if '--params' in sys.argv:
        sys.argv[sys.argv.index('--params')] = '-p'
    if '-p' in sys.argv:
        arg_index = sys.argv.index('-p')
        params = '-p' + sys.argv.pop(arg_index + 1)
        sys.argv[arg_index] = params

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

    br.flags(parser, ("alignments", "The file(s) you want to start working on"),
             br.alb_flags, br.alb_modifiers, VERSION)

    in_args = parser.parse_args()

    alignbuddy = []
    align_set = ""

    if in_args.out_format and in_args.out_format.lower() not in OUTPUT_FORMATS:
        _stderr("Error: Output type %s is not recognized/supported\n" % in_args.out_format)
        sys.exit()

    try:
        if not in_args.generate_alignment:  # If passing in sequences to do alignment, don't make AlignBuddy obj
            for align_set in in_args.alignments:
                if isinstance(align_set, TextIOWrapper) and align_set.buffer.raw.isatty():
                    sys.exit("Warning: No input detected. Process will be aborted.")
                align_set = AlignBuddy(align_set, in_args.in_format, in_args.out_format)
                alignbuddy += align_set.alignments

            alignbuddy = AlignBuddy(alignbuddy, align_set._in_format, align_set._out_format)

    except GuessError as e:
        _stderr("GuessError: %s\n" % e, in_args.quiet)
        sys.exit()

    except ValueError as e:
        _stderr("ValueError: %s\n" % e, in_args.quiet)
        sys.exit()

    return in_args, alignbuddy


def command_line_ui(in_args, alignbuddy, skip_exit=False):
    # ############################################# INTERNAL FUNCTIONS ############################################## #
    def _print_aligments(_alignbuddy):
        try:
            _output = str(_alignbuddy)
        except ValueError as err:
            _stderr("ValueError: %s\n" % str(err))
            return False
        except TypeError as err:
            _stderr("TypeError: %s\n" % str(err))
            return False
        except br.PhylipError as err:
            _stderr("PhylipError: %s\n" % str(err))
            return False

        if in_args.test:
            _stderr("*** Test passed ***\n", in_args.quiet)
            pass

        elif in_args.in_place:
            _in_place(_output, in_args.alignment[0])

        else:
            _stdout("{0}\n".format(_output.rstrip()))
        return True

    def _in_place(_output, _path):
        if not os.path.exists(_path):
            _stderr("Warning: The -i flag was passed in, but the positional argument doesn't seem to be a "
                    "file. Nothing was written.\n", in_args.quiet)
            _stderr("%s\n" % _output.rstrip(), in_args.quiet)
        else:
            with open(os.path.abspath(_path), "w") as _ofile:
                _ofile.write(_output)
            _stderr("File over-written at:\n%s\n" % os.path.abspath(_path), in_args.quiet)

    def _exit(tool, skip=skip_exit):
        if skip:
            return
        usage = br.Usage()
        usage.increment("AlignBuddy", VERSION.short(), tool)
        usage.save()
        sys.exit()

    def _raise_error(_err, tool, check_string=None):
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
        _exit(tool)

    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Alignment lengths
    if in_args.alignment_lengths:
        counts = alignment_lengths(alignbuddy)
        output = ""
        for indx, count in enumerate(counts):
            output += "# Alignment %s\n%s\n\n" % (indx + 1, count) if len(counts) > 1 else "%s\n" % count

        _stdout("%s\n" % output.strip())
        _exit("alignment_lengths")

    # Back Transcribe
    if in_args.back_transcribe:
        if alignbuddy.alpha != IUPAC.ambiguous_rna:
            _raise_error(ValueError("You need to provide an RNA sequence."), "back_transcribe")
        try:
            _print_aligments(rna2dna(alignbuddy))
        except TypeError as e:
            _raise_error(e, "back_transcribe", "Nucleic acid sequence required, not protein.")
        _exit("back_transcribe")

    # Clean Seq
    if in_args.clean_seq:
        if in_args.clean_seq[0] == "strict":
            _print_aligments(clean_seq(alignbuddy))
        else:
            _print_aligments(clean_seq(alignbuddy, skip_list="RYWSMKHBVDNXrywsmkhbvdnx"))
        _exit("clean_seq")

    # Codon alignment
    if in_args.codon_alignment:
        try:
            _print_aligments(codon_alignment(alignbuddy))
        except TypeError as e:
            _raise_error(e, "codon_alignment", "Nucleic acid sequence required.")
        _exit("codon_alignment")

    # Concatenate Alignments
    if in_args.concat_alignments:
        try:
            _print_aligments(concat_alignments(alignbuddy, in_args.concat_alignments))
        except AttributeError as e:
            _raise_error(e, "concat_alignments", "Please provide at least two alignments.")
        _exit("concat_alignments")

    # Delete rows
    if in_args.delete_rows:
        _print_aligments(delete_rows(alignbuddy, in_args.delete_rows))
        _exit("delete_rows")

    # Extract range
    if in_args.extract_range:
        try:
            _print_aligments(extract_range(alignbuddy, *in_args.extract_range))
        except ValueError as e:
            _raise_error(e, "extract_range", "Error at extract range: The value given for end of range is smaller")
        _exit("extract_range")

    # Generate Alignment
    if in_args.generate_alignment:
        seqbuddy = []
        seq_set = None
        try:
            import SeqBuddy as Sb
        except ImportError:
            _raise_error(ImportError("SeqBuddy is needed to use generate_msa(). Please install it and try again."),
                         "generate_alignment")

        for seq_set in in_args.alignments:
            if isinstance(seq_set, TextIOWrapper) and seq_set.buffer.raw.isatty():
                _raise_error(BufferError("No input detected. Process will be aborted."), "generate_alignment")

            seq_set = Sb.SeqBuddy(seq_set, in_args.in_format, in_args.out_format)
            seqbuddy += seq_set.records
        if seq_set:
            seqbuddy = Sb.SeqBuddy(seqbuddy, seq_set.in_format, seq_set.out_format)
        else:
            seqbuddy = Sb.SeqBuddy(seqbuddy, in_args.in_format, in_args.out_format)
        params = in_args.params if in_args.params is None else in_args.params[0]
        generated_msas = generate_msa(seqbuddy, in_args.generate_alignment[0], params, in_args.keep_temp, in_args.quiet)
        if in_args.out_format:
            generated_msas.set_format(in_args.out_format)
        try:
            _stdout(str(generated_msas))
        except AttributeError as e:
            _raise_error(e, "generate_msa", "is not a valid alignment tool")
        _exit("generate_alignment")

    # List identifiers
    if in_args.list_ids:
        columns = 1 if not in_args.list_ids[0] or in_args.list_ids == 0 else abs(in_args.list_ids[0])
        listed_ids = list_ids(alignbuddy)
        output = ""
        for indx, alignment in enumerate(listed_ids):
            count = 1
            output += "# Alignment %s\n" % str(indx + 1)
            for identifier in alignment:
                if count < columns:
                    output += "%s\t" % identifier
                    count += 1
                else:
                    output += "%s\n" % identifier
                    count = 1
            output += "\n"
        _stdout(output)
        _exit("list_ids")

    # Lowercase
    if in_args.lowercase:
        _print_aligments(lowercase(alignbuddy))
        _exit("lowercase")

    # Number sequences per alignment
    if in_args.num_seqs:
        counts = num_seqs(alignbuddy)
        output = ""
        for indx, count in enumerate(counts):
            output += "# Alignment %s\n%s\n\n" % (indx + 1, count) if len(counts) > 1 else "%s\n" % count
        _stdout("%s\n" % output.strip())
        _exit("num_seqs")

    # Order IDs
    if in_args.order_ids:
        reverse = True if in_args.order_ids[0] and in_args.order_ids[0] == "rev" else False
        _print_aligments(order_ids(alignbuddy, reverse=reverse))
        _exit("order_ids")

    # Pull rows
    if in_args.pull_rows:
        _print_aligments(pull_rows(alignbuddy, in_args.pull_rows))
        _exit("pull_rows")

    # Rename IDs
    if in_args.rename_ids:
        num = 0 if not in_args.params else int(in_args.params[0])
        _print_aligments(rename(alignbuddy, in_args.rename_ids[0], in_args.rename_ids[1], num))
        _exit("rename_ids")

    # Screw formats
    if in_args.screw_formats:
        alignbuddy.set_format(in_args.screw_formats)
        _print_aligments(alignbuddy)
        _exit("screw_formats")

    # Split alignments into files
    if in_args.split_to_files:
        in_args.in_place = True
        out_dir = os.path.abspath(in_args.split_to_files[0])
        filename = in_args.split_to_files[1]
        os.makedirs(out_dir, exist_ok=True)
        check_quiet = in_args.quiet  # 'quiet' must be toggled to 'on' _print_recs() here.
        in_args.quiet = True
        for indx, buddy in enumerate(split_alignbuddy(alignbuddy)):
            alignbuddy.alignments = buddy.alignments
            ext = br.format_to_extension[alignbuddy._out_format]
            in_args.alignment[0] = "%s/%s_%s.%s" % (out_dir, filename, '{:0>4d}'.format(indx + 1), ext)
            _stderr("New file: %s\n" % in_args.alignment[0], check_quiet)
            open(in_args.alignment[0], "w").close()
            _print_aligments(alignbuddy)
        _exit("split_to_files")

    # Translate CDS
    if in_args.translate:
        try:
            _print_aligments(translate_cds(alignbuddy, quiet=in_args.quiet))
        except TypeError as e:
            _raise_error(e, "translate", "Nucleic acid sequence required, not protein.")
        _exit("translate")

    # Trimal
    if in_args.trimal:
        in_args.trimal = 1.0 if not in_args.trimal[0] else in_args.trimal[0]
        _print_aligments(trimal(alignbuddy, in_args.trimal))
        _exit("trimal")

    # Transcribe
    if in_args.transcribe:
        if alignbuddy.alpha != IUPAC.ambiguous_dna:
            _raise_error(ValueError("You need to provide a DNA sequence."), "transcribe")
        try:
            _print_aligments(dna2rna(alignbuddy))
        except TypeError as e:
            _raise_error(e, "transcribe", "Nucleic acid sequence required, not protein.")
        _exit("transcribe")

    # Uppercase
    if in_args.uppercase:
        _print_aligments(uppercase(alignbuddy))
        _exit("uppercase")

if __name__ == '__main__':
    try:
        command_line_ui(*argparse_init())
    except (KeyboardInterrupt, GuessError) as _e:
        print(_e)
    except SystemExit:
        pass
    except Exception as _e:
        br.send_traceback("AlignBuddy", _e)
