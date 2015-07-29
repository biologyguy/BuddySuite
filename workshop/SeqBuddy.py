#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 44 tools and counting

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
from copy import copy, deepcopy
from random import sample, choice, randint, random
from math import floor, ceil, log
from tempfile import TemporaryDirectory
from subprocess import Popen, PIPE
from shutil import which
from multiprocessing import Manager
import ctypes
from hashlib import md5
from io import StringIO
from collections import OrderedDict

# Third party package imports
sys.path.insert(0, "./")  # For stand alone executable, where dependencies are packaged with BuddySuite
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
from Bio import Restriction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError
from Bio.SeqUtils import GC

# My functions
from MyFuncs import run_multicore_function


# ##################################################### WISH LIST #################################################### #
def sim_ident(matrix):  # Return the pairwise similarity and identity scores among sequences
    x = matrix
    return x


def codon_counter():
    # generate frequency statistics for codon composition
    return


def predict_orfs():
    return


def shuffle_seqs():
    # randomly reorder the residues in each sequence
    return


def delete_pattern():
    # remove residues that match a given pattern from all records
    return


def insert_sequence():
    # Add a specific sequence at a defined location in all records. E.g., adding a barcode
    return


def mutate():
    # Apply some model of evolution to generate new sequences from input
    return


def random_aa(_length, number, matrix):
    # create random protein sequences
    x = [_length, number, matrix]
    return x


def random_dna(_length, number, matrix):
    # create random DNA sequences
    x = [_length, number, matrix]
    return x


def divergence_value():
    # http://bioinformatics.org/sms/uneven.html
    return

# - Allow batch calls. E.g., if 6 files are fed in as input, run the SeqBuddy command provided independently on each
# - Add FASTQ support... More generally, support letter annotation mods
# - Add Clustal support
# - Make an 'installer' that puts SeqBuddy into path and adds `sb` sym link
# - Get BuddySuite into PyPi
# - Check on memory requirements before execution
# - Execution timer, for long running jobs
# - Handle all stderr output from private function (to allow quiet execution)
# - Sort out a good way to manage 'lazy' imports

# ################################################# CHANGE LOG for V2 ################################################ #
# - New flag -t/--test, which runs a function but suppresses all stdout (stderr still returned)
# - New function split_by_taxa(). Writes individual files for groups of sequences based on an identifier in their ids
# - Standard-in is handled as input now, allowing SeqBuddy to be daisy chained with pipes (|)
# - Remove the the -p flag dependencies for -prr, -li, -btr, -asl, -cts, -hsi, -frp, -drp, -ofa, -ofp, and -oi
# - Add print method to SeqBuddy class that outputs all sequences to string
# - Add molecular_weight() method that calculates molecular weight
# - Add isoelectric_point() method that calculates isoelectric point
# - Unit tests


# ################################################# HELPER FUNCTIONS ################################################# #
class GuessError(Exception):
    """Raised when input format cannot be guessed"""
    def __init__(self, _value):
        self.value = _value

    def __str__(self):
        return self.value


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


def _stderr(message, quiet=False):
    if not quiet:
        sys.stderr.write(message)
    return


def _stdout(message, quiet=False):
    if not quiet:
        sys.stdout.write(message)
    return


def _format_to_extension(_format):
    format_to_extension = {'fasta': 'fa', 'fa': 'fa', 'genbank': 'gb', 'gb': 'gb', 'nexus': 'nex',
                           'nex': 'nex', 'phylip': 'phy', 'phy': 'phy', 'phylip-relaxed': 'phyr', 'phyr': 'phyr',
                           'stockholm': 'stklm', 'stklm': 'stklm'}
    return format_to_extension[_format]
# ##################################################### SEQ BUDDY #################################################### #


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
                temp = StringIO(_input.read())
                _input = temp
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
            self.alpha = guess_alphabet(_sequences)
        elif self.alpha in ['protein', 'prot', 'p', 'pep', IUPAC.protein]:
            self.alpha = IUPAC.protein
        elif self.alpha in ['dna', 'd', 'cds', IUPAC.ambiguous_dna]:
            self.alpha = IUPAC.ambiguous_dna
        elif self.alpha in ['rna', 'r', IUPAC.ambiguous_rna]:
            self.alpha = IUPAC.ambiguous_rna
        else:
            _stderr("WARNING: Alphabet not recognized. Correct alphabet will be guessed.\n")
            self.alpha = guess_alphabet(_sequences)

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
            _output = phylipi(self)

        elif self.out_format == "phylipis":
            _output = phylipi(self, "strict")

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


# Does not attempt to explicitly deal with weird cases (e.g., ambigous residues).
# The user will need to specify an alphabet with the -a flag if using many non-standard characters in their sequences.
def guess_alphabet(_seqbuddy):
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


def guess_format(_input):  # _input can be list, SeqBuddy object, file handle, or file path.
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

        possible_formats = ["phylip-relaxed", "stockholm", "fasta", "gb", "fastq", "nexus"]  # ToDo: Glean CLUSTAL
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
        return None  # Unable to determine format from file handle

    else:
        raise GuessError("Unsupported _input argument in guess_format(). %s" % _input)


def phylipi(_seqbuddy, _format="relaxed"):  # _format in ["strict", "relaxed"]
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
# #################################################################################################################### #


def blast(_seqbuddy, blast_db, blast_path=None, blastdbcmd=None):  # ToDo: Allow weird binary names to work
    if not blast_path:
        blast_path = which("blastp") if _seqbuddy.alpha == IUPAC.protein else which("blastn")

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
            raise FileNotFoundError("blastp binary not found")

        if not os.path.isfile("%s.pin" % blast_db) or not os.path.isfile("%s.phr" % blast_db) \
                or not os.path.isfile("%s.psq" % blast_db):
            raise RuntimeError("Blastp database not found at '%s'" % blast_db)
    elif blast_check == "blastn":
        if not which(blast_path):
            raise FileNotFoundError("blastn binary not found")

        if not os.path.isfile("%s.nin" % blast_db) or not os.path.isfile("%s.nhr" % blast_db) \
                or not os.path.isfile("%s.nsq" % blast_db):
            raise RuntimeError("Blastn database not found at '%s'" % blast_db)
    else:
        raise RuntimeError("Blast binary doesn't seem to work, at %s" % blast_path)

    if not blastdbcmd:
        blastdbcmd = "blastdbcmd"

    if not which(blastdbcmd):
        raise FileNotFoundError("blastdbcmd")

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


def shuffle(_seqbuddy):
    _output = []
    for _ in range(len(_seqbuddy.records)):
        random_index = randint(1, len(_seqbuddy.records)) - 1
        _output.append(_seqbuddy.records.pop(random_index))
    _seqbuddy.records = _output
    return _seqbuddy


def order_ids(_seqbuddy, _reverse=False):
    _output = [(_rec.id, _rec) for _rec in _seqbuddy.records]
    _output = sorted(_output, key=lambda x: x[0], reverse=_reverse)
    _output = [_rec[1] for _rec in _output]
    _seqbuddy.records = _output
    return _seqbuddy


def rna2dna(_seqbuddy):
    if _seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq.back_transcribe()), alphabet=IUPAC.ambiguous_dna)
    _seqbuddy.alpha = IUPAC.ambiguous_dna
    return _seqbuddy


def dna2rna(_seqbuddy):
    if _seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq.transcribe()), alphabet=IUPAC.ambiguous_rna)
    _seqbuddy.alpha = IUPAC.ambiguous_rna
    return _seqbuddy


def complement(_seqbuddy):
    if _seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for _rec in _seqbuddy.records:
        _rec.seq = _rec.seq.complement()
    return _seqbuddy


def reverse_complement(_seqbuddy):
    if _seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for _rec in _seqbuddy.records:
        _rec.seq = _rec.seq.reverse_complement()
        seq_len = len(_rec.seq)
        shifted_features = [_feature_rc(_feature, seq_len) for _feature in _rec.features]
        _rec.features = shifted_features
    return _seqbuddy


# ToDo: Deal with alignments...
def translate_cds(_seqbuddy, quiet=False):  # adding 'quiet' will suppress the errors thrown by translate(cds=True)
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
                    codon = str(temp_seq.seq)[(_i * 3):(_i * 3 + 3)]
                    if codon.upper() in ["TGA", "TAG", "TAA"]:
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
    return _output


def select_frame(_seqbuddy, frame):  # ToDo: record the deleted residues so the earlier frame can be returned to.
    if _seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Select frame requires nucleic acid, not protein.")
    for _rec in _seqbuddy.records:
        _rec.features = _shift_features(_rec.features, (frame - 1) * -1, len(_rec.seq))
        _rec.seq = Seq(str(_rec.seq)[frame - 1:], alphabet=_rec.seq.alphabet)
    return _seqbuddy


def translate6frames(_seqbuddy):
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


def back_translate(_seqbuddy, _mode='random', _species=None):
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


def ave_seq_length(_seqbuddy, _clean=False):
    if _clean:  # Strip out all gaps and stuff before counting
        clean_seq(_seqbuddy)

    sum_length = 0.
    for _rec in _seqbuddy.records:
        sum_length += len(_rec.seq)
    return sum_length / len(_seqbuddy.records)


def concat_seqs(_seqbuddy, _clean=False):
    if _clean:
        clean_seq(_seqbuddy)

    _new_seq = ""
    concat_ids = []
    features = []
    for _rec in _seqbuddy.records:
        _shift = len(_new_seq)
        full_seq_len = len(_new_seq) + len(str(_rec.seq))
        _rec.features = _shift_features(_rec.features, _shift, full_seq_len)

        location = FeatureLocation(len(_new_seq), len(_new_seq) + len(str(_rec.seq)))
        feature = SeqFeature(location=location, id=_rec.id, type=_rec.id[:15])
        features.append(feature)
        features += _rec.features
        concat_ids.append(_rec.id)
        _new_seq += str(_rec.seq)

    _new_seq = [SeqRecord(Seq(_new_seq, alphabet=_seqbuddy.alpha),
                          description="", id="concatination", features=features)]
    _seqbuddy = SeqBuddy(_new_seq)
    _seqbuddy.out_format = "gb"
    return _seqbuddy


def clean_seq(_seqbuddy, skip_list=None, ambiguous=True):
    """remove all non-sequence charcters from sequence strings"""
    skip_list = "" if not skip_list else "".join(skip_list)
    for _rec in _seqbuddy.records:
        if _seqbuddy.alpha == IUPAC.protein:
            _rec.seq = Seq(re.sub("[^ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy%s]" % skip_list, "", str(_rec.seq)),
                           alphabet=_seqbuddy.alpha)
        else:
            if ambiguous:
                _rec.seq = Seq(re.sub("[^ATGCURYWSMKHBVDNXatgcurywsmkhbvdnx%s]" % skip_list, "", str(_rec.seq)), alphabet=_seqbuddy.alpha)
            else:
                _rec.seq = Seq(re.sub("[^ATGCUatgcu%s]" % skip_list, "", str(_rec.seq)), alphabet=_seqbuddy.alpha)

    return _seqbuddy


def delete_metadata(_seqbuddy):
    _new_seqs = []
    for _rec in _seqbuddy.records:
        _new_seqs.append(SeqRecord(Seq(str(_rec.seq), alphabet=_seqbuddy.alpha), id=_rec.id, name='', description=''))
    _seqbuddy.records = _new_seqs
    return _seqbuddy


# Apply DNA features to protein sequences
def map_features_dna2prot(dna_seqbuddy, prot_seqbuddy):
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
            location = FeatureLocation(floor(_start), floor(_end))
            _feature.location = location

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


# Apply Protein features to DNA sequences
def map_features_prot2dna(prot_seqbuddy, dna_seqbuddy):
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
            location = FeatureLocation(_start, _end)
            _feature.location = location

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


# Merge feature lists
def combine_features(_seqbuddy1, _seqbuddy2):  # ToDo: rewrite this to accept any number of input files.
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


def order_features_by_position(_seqbuddy, _reverse=False):
    for _rec in _seqbuddy.records:
        new_feature_list = [(int(_feature.location.start), _feature) for _feature in _rec.features]
        new_feature_list = sorted(new_feature_list, key=lambda x: x[0], reverse=_reverse)
        new_feature_list = [_feature[1] for _feature in new_feature_list]
        _rec.features = new_feature_list
    return _seqbuddy


def order_features_alphabetically(_seqbuddy, _reverse=False):
    for _rec in _seqbuddy.records:
        new_feature_list = [(_feature.type, _feature) for _feature in _rec.features]
        new_feature_list = sorted(new_feature_list, key=lambda x: x[0], reverse=_reverse)
        new_feature_list = [_feature[1] for _feature in new_feature_list]
        _rec.features = new_feature_list
    return _seqbuddy


def hash_sequence_ids(_seqbuddy, _hash_length=10):
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

    _hash_map = []
    for i in range(len(hash_list)):
        _hash_map.append((hash_list[i], seq_ids[i]))

    _hash_table = "# Hash table\n"
    for _seq in _hash_map:
        _hash_table += "%s,%s\n" % (_seq[0], _seq[1])

    return [_seqbuddy, _hash_map, _hash_table]


def pull_recs(_seqbuddy, _search):
    matched_records = []
    for _rec in _seqbuddy.records:
        if re.search(_search, _rec.description) or re.search(_search, _rec.id) or re.search(_search, _rec.name):
            matched_records.append(_rec)
    _seqbuddy.records = matched_records
    return _seqbuddy


def pull_random_recs(_seqbuddy, _count=1):  # Return a random set of sequences (without replacement)
    random_recs = []
    _count = abs(_count) if abs(_count) <= len(_seqbuddy.records) else len(_seqbuddy.records)
    for i in range(_count):
        rand_index = randint(0, len(_seqbuddy.records) - 1)
        random_recs.append(_seqbuddy.records.pop(rand_index))

    _seqbuddy.records = random_recs
    return _seqbuddy


def pull_record_ends(_seqbuddy, _amount, _which_end):
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
            raise AttributeError("You much pick 'front' or 'rear' for the '_which_end' argument in pull_record_ends.")

        seq_ends.append(_rec)

    _seqbuddy.records = seq_ends
    return _seqbuddy


def extract_range(_seqbuddy, _start, _end):
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


def find_repeats(_seqbuddy, _columns=1):
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


def delete_records(_seqbuddy, search_str):
    retained_records = []
    _deleted = pull_recs(copy(_seqbuddy), search_str).records
    for _rec in _seqbuddy.records:
        if _rec in _deleted:
            continue
        else:
            retained_records.append(_rec)
    _seqbuddy.records = retained_records
    return _seqbuddy


def delete_large(_seqbuddy, max_value):
    retained_records = []
    for _rec in _seqbuddy.records:
        if len(str(_rec.seq)) <= max_value:
            retained_records.append(_rec)
    _seqbuddy.records = retained_records
    return _seqbuddy


def delete_small(_seqbuddy, min_value):
    retained_records = []
    for _rec in _seqbuddy.records:
        if len(str(_rec.seq)) >= min_value:
            retained_records.append(_rec)
    _seqbuddy.records = retained_records
    return _seqbuddy


def delete_features(_seqbuddy, _pattern):
    for _rec in _seqbuddy.records:
        retained_features = []
        for _feature in _rec.features:
            if not re.search(_pattern, _feature.type):
                retained_features.append(_feature)
        _rec.features = retained_features
    return _seqbuddy


def delete_repeats(_seqbuddy, scope='all'):  # scope in ['all', 'ids', 'seqs']
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


def rename(_seqbuddy, query, replace="", _num=0):  # TODO Allow a replacement pattern increment (like numbers)
    for _rec in _seqbuddy.records:
        new_name = re.sub(query, replace, _rec.id, _num)
        _rec.id = new_name
        _rec.name = new_name
    return _seqbuddy


def purge(_seqbuddy, threshold):  # ToDo: Implement a way to return a certain # of seqs (i.e. auto-determine threshold)
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


def bl2seq(_seqbuddy):  # Does an all-by-all analysis, and does not return sequences
    """
    Note on blast2seq: Expect (E) values are calculated on an assumed database size of (the rather large) nr, so the
    threshold may need to be increased quite a bit to return short alignments
    """
    def mc_blast(_query, args):
        _subject_file = args[0]

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
        raise RuntimeError("%s not present in $PATH." % blast_bin)  # ToDo: Implement -p flag

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
            output_dict[query] = {subj: [float(_ident), int(_length), float(_evalue), int(_bit_score)]}
        else:
            output_dict[query][subj] = [float(_ident), int(_length), float(_evalue), int(_bit_score)]

        if subj not in output_dict:
            output_dict[subj] = {query: [float(_ident), int(_length), float(_evalue), int(_bit_score)]}
        else:
            output_dict[subj][query] = [float(_ident), int(_length), float(_evalue), int(_bit_score)]

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


def uppercase(_seqbuddy):
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq).upper(), alphabet=_rec.seq.alphabet)
    return _seqbuddy


def lowercase(_seqbuddy):
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq).lower(), alphabet=_rec.seq.alphabet)
    return _seqbuddy


def split_by_taxa(_seqbuddy, split_pattern):
    recs_by_taxa = {}
    for _rec in _seqbuddy.records:
        split = re.split(split_pattern, _rec.id)
        recs_by_taxa.setdefault(split[1], []).append(_rec)
    return recs_by_taxa


def molecular_weight(_seqbuddy):
    # get the mass of each sequence in daltons

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
    elif _seqbuddy.alpha == (IUPAC.ambiguous_dna or IUPAC.unambiguous_dna):
        _dict = deoxynucleotide_weights
        _dna = True
    elif _seqbuddy.alpha == (IUPAC.ambiguous_rna or IUPAC.unambiguous_rna):
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
        if _dna:
            _output['masses_ds'].append(round(_rec.mass_ds, 3))
        _output['ids'].append(_rec.id)
    return _output


def isoelectric_point(_seqbuddy):
    if _seqbuddy.alpha is not IUPAC.protein:
        raise TypeError("Protein sequence required, not nucleic acid.")
    _isoelectric_points = []
    for _rec in _seqbuddy.records:
        _pI = ProteinAnalysis(str(_rec.seq))
        _pI = round(_pI.isoelectric_point(), 10)
        _isoelectric_points.append((_rec.id, _pI))
    return _isoelectric_points


def count_residues(_seqbuddy):
    # TODO add % Neutral
    # TODO distinguish between ambiguous and unambiguous DNA
    _output = []
    for _rec in _seqbuddy.records:
        if _seqbuddy.alpha is IUPAC.protein:
            resid_count = OrderedDict.fromkeys(('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
                                                'P', 'S', 'T', 'W', 'Y', 'V', 'X', '% Ambiguous', '% Positive', 
                                                '% Negative'), 0)
            
        else:
            resid_count = OrderedDict.fromkeys(('A', 'G', 'C', 'T', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B',
                                                'X', 'N', 'U'), 0)
        num_acids = 0
        for char in str(_rec.seq).upper():
            if char in resid_count:
                resid_count[char] += 1
                num_acids += 1
        if _seqbuddy.alpha is IUPAC.protein:
            resid_count['% Ambiguous'] = round(100 * (resid_count['X']) / num_acids, 2)
            resid_count['% Positive'] = round(100 * (resid_count['H'] + resid_count['K'] +
                                                     resid_count['R']) / num_acids, 2)
            resid_count['% Negative'] = round(100 * (resid_count['D'] + resid_count['E'] +
                                                     resid_count['C'] + resid_count['Y']) / num_acids, 2)
        else:
            resid_count['% Ambiguous'] = round(100 * ((num_acids - (resid_count['A'] + resid_count['T'] +
                                                                    resid_count['C'] + resid_count['G'] +
                                                                    resid_count['U'])) / num_acids), 2)
        _output.append((_rec.id, resid_count))
    return sorted(_output, key=lambda x: x[0])


def raw_seq(_seqbuddy):
    _seqbuddy = clean_seq(_seqbuddy)
    _output = ""
    for _rec in _seqbuddy.records:
        _output += "%s\n\n" % _rec.seq

    return "%s\n" % _output.strip()


def list_ids(_seqbuddy, _columns=1):
    _columns = 1 if _columns == 0 else abs(_columns)
    _output = ""
    _counter = 1
    for rec in _seqbuddy.records:
        _output += "%s\t" % rec.id
        if _counter % _columns == 0:
            _output = "%s\n" % _output.strip()
        _counter += 1
    return "%s\n" % _output.strip()


def num_seqs(_seqbuddy):
    return len(_seqbuddy.records)


def merge(_seqbuddy_list):
    _output = _seqbuddy_list[0]
    for _seqbuddy in _seqbuddy_list[1:]:
        _output.records += _seqbuddy.records
    return _output


def split_file(_seqbuddy):
    # Split the records in a SeqBuddy object up into a collection of new SeqBuddy objects
    sb_objs_list = []
    for _rec in _seqbuddy.records:
        _sb = SeqBuddy([_rec])
        _sb.in_format = _seqbuddy.in_format
        _sb.out_format = _seqbuddy.out_format
        sb_objs_list.append(_sb)
    return sb_objs_list


def find_restriction_sites(_seqbuddy, _commercial=True, _single_cut=True):
    sites = []
    blacklist = [Restriction.AbaSI, Restriction.FspEI, Restriction.MspJI, Restriction.SgeI, Restriction.AspBHI,
                 Restriction.SgrTI, Restriction.YkrI, Restriction.BmeDI]
    for rec in _seqbuddy.records:
        if _commercial:
            analysis = Restriction.Restriction.Analysis(Restriction.CommOnly, rec.seq)
        else:
            analysis = Restriction.Restriction.Analysis(Restriction.AllEnzymes, rec.seq)
        if _single_cut:
            cut_keys = set(analysis.blunt().keys())
        else:
            cut_keys = set(analysis.overhang5()).union(set(analysis.overhang3()))
        comm_keys = set(analysis.with_sites().keys())
        keys = cut_keys & comm_keys
        for _key in blacklist:
            keys.discard(_key)
        result = {_key: analysis.with_sites()[_key] for _key in keys}
        sites.append((rec.id, result))
    print_output = ''
    for tup in sites:
        print_output += "{0}\n".format(tup[0])
        for _key in sorted(tup[1]):
            print_output += "{0}:\t{1}\n".format(_key, str(tup[1][_key]))
        print_output += "\n"
    return [sites, print_output]


def find_pattern(_seqbuddy, _pattern):  # TODO ambiguous letters mode
    # search through sequences for regex matches. For example, to find micro-RNAs
    _pattern = _pattern.upper()
    _output = OrderedDict()
    for _rec in _seqbuddy.records:
        indices = []
        matches = re.finditer(_pattern, str(_rec.seq).upper())
        for match in matches:
            indices.append(match.start())
        _output[_rec.id] = indices
    return _output


def find_CpG(_seqbuddy):
    # Predicts locations of CpG islands in DNA sequences
    _seqbuddy = clean_seq(_seqbuddy)
    if _seqbuddy.alpha not in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna]:
        raise TypeError("DNA sequence required, not protein or RNA.")

    _output = OrderedDict()
    records = []

    def cpg_calc(_seq):  # Returns observed/expected value of a sequence
        _seq = _seq.upper()
        observed_cpg = len(re.findall("CG", _seq)) * len(_seq)
        expected = (len(re.findall("[CG]", _seq)) / 2) ** 2
        return observed_cpg / expected

    def cg_percent(_seq):  # Returns the CG % of a sequence
        _seq = _seq.upper()
        return len(re.findall("[CG]", _seq)) / len(_seq)

    def find_islands(cg_percents, oe_values):  # Returns a list of tuples containing the start and end of an island
        _indx = 0
        _output = []
        while _indx < len(cg_percents):
            if cg_percents[_indx] > .5 and oe_values[_indx] > .6:
                start = _indx
                while _indx < len(cg_percents) and cg_percents[_indx] > .5 and oe_values[_indx] > .6:
                    _indx += 1
                end = _indx
                _output.append((start, end))
            _indx += 1
        return _output

    def map_cpg(_seq, island_ranges):  # Maps CpG islands onto a sequence as capital letters
        cpg_seq = _seq.lower()
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

# ################################################# COMMAND LINE UI ################################################## #
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="SeqBuddy.py", formatter_class=argparse.RawTextHelpFormatter,
                                     description="Commandline wrapper for all the fun functions in this file. "
                                                 "Play with your sequences!")

    parser.add_argument("sequence", help="Supply a file path or a raw sequence", nargs="*", default=[sys.stdin])

    parser.add_argument('-v', '--version', action='version',
                        version='''\
SeqBuddy 2.alpha (2015)

Gnu General Public License, Version 2.0 (http://www.gnu.org/licenses/gpl.html)
This is free software; see the source for detailed copying conditions.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
Questions/comments/concerns can be directed to Steve Bond, steve.bond@nih.gov''')
    # TODO Fix --clean_seq for .phy, .phyr, .stklm, .nex
    parser.add_argument('-cs', '--clean_seq', action='append', nargs="?",
                        help="Strip out non-sequence characters, such as stops (*) and gaps (-). Pass in the word "
                             "'strict' to remove all characters except the unambiguous letter codes.")
    parser.add_argument('-uc', '--uppercase', action='store_true',
                        help='Convert all sequences to uppercase')  # TODO Fix for genbank
    parser.add_argument('-lc', '--lowercase', action='store_true',
                        help='Convert all sequences to lowercase')
    parser.add_argument('-dm', '--delete_metadata', action='store_true',
                        help="Remove meta-data from file (only id is retained)")
    parser.add_argument('-rs', '--raw_seq', action='store_true',
                        help="Return line break separated sequences")
    parser.add_argument('-tr', '--translate', action='store_true',
                        help="Convert coding sequences into amino acid sequences")
    parser.add_argument('-sfr', '--select_frame', action='store', metavar='<frame (int)>', type=int, choices=[1, 2, 3],
                        help="Change the reading frame of sequences by deleting characters off of the front")
    parser.add_argument('-tr6', '--translate6frames', action='store_true',
                        help="Translate nucleotide sequences into all six reading frames")
    parser.add_argument('-btr', '--back_translate', action='append', nargs="*",
                        help="Convert amino acid sequences into codons. Optionally, select mode by passing in "
                             "['random', 'r', 'optimized', 'o'] "
                             "['human', 'h', 'mouse', 'm', 'yeast', 'y', 'ecoli', 'e']")
    parser.add_argument('-d2r', '--transcribe', action='store_true',
                        help="Convert DNA sequences to RNA")
    parser.add_argument('-r2d', '--back_transcribe', action='store_true',
                        help="Convert RNA sequences to DNA")
    parser.add_argument('-cmp', '--complement', action='store_true',
                        help="Return complement of nucleotide sequence")
    parser.add_argument('-rc', '--reverse_complement', action='store_true',
                        help="Return reverse complement of nucleotide sequence")
    parser.add_argument('-li', '--list_ids', nargs='?', action='append', type=int, metavar='int (optional)',
                        help="Output all the sequence identifiers in a file. Optionally, pass in an integer to "
                             "specify the # of columns to write")
    parser.add_argument('-ns', '--num_seqs', action='store_true',
                        help="Counts how many sequences are present in an input file")
    parser.add_argument('-asl', '--ave_seq_length', action='append', nargs="?",
                        help="Return the average length of all sequences. Pass in the word 'clean' to remove gaps etc "
                             "from the sequences before counting.")
    parser.add_argument('-cts', '--concat_seqs', action='append', nargs="?",
                        help="Concatenate a bunch of sequences into a single solid string. Pass in the word 'clean' to "
                             "remove stops, gaps, etc., from the sequences before concatenating.")
    parser.add_argument('-fd2p', '--map_features_dna2prot', action='store_true',
                        help="Take the features annotated onto nucleotide sequences and map to protein sequences. "
                             "Both a protein and cDNA file must be passed in.")
    parser.add_argument('-fp2d', '--map_features_prot2dna', action='store_true',
                        help="Arguments: one cDNA file and one protein file")
    parser.add_argument('-ri', '--rename_ids', action='store', metavar=('<pattern>', '<substitution>'), nargs=2,
                        help="Replace some pattern in ids with something else. Limit number of replacements with -p.")
    parser.add_argument('-cf', '--combine_features', action='store_true',
                        help="Takes the features in two files and combines them for each sequence")
    parser.add_argument('-sh', '--shuffle', action='store_true',
                        help="Randomly reorder the position of records in the file.")
    parser.add_argument('-oi', '--order_ids', action='append', nargs="?",
                        help="Sort all sequences by id in alpha-numeric order. Pass in the word 'rev' to reverse order")
    parser.add_argument('-ofp', '--order_features_by_position', action='append', nargs="?",
                        help="Change the output order of sequence features, based on sequence position. Pass in 'rev' "
                             "to reverse order.")
    parser.add_argument('-ofa', '--order_features_alphabetically', action='append', nargs="?",
                        help="Change the output order of sequence features, based on sequence position. Pass in 'rev' "
                             "to reverse order.")
    parser.add_argument('-sf', '--screw_formats', action='store', metavar="<out_format>",
                        help="Change the file format to something else.")
    parser.add_argument('-hsi', '--hash_seq_ids', action='append', nargs="?", type=int,
                        help="Rename all the identifiers in a sequence list to a fixed length hash. Default 10; "
                             "override by passing in an integer.")
    parser.add_argument('-pr', '--pull_records', action='store', metavar="<regex pattern>",
                        help="Get all the records with ids containing a given string")
    parser.add_argument('-prr', '--pull_random_record', nargs='?', action='append', type=int, metavar='int (optional)',
                        help="Extract random sequences. Optionally, pass in an integer to increase the number of "
                             "sequences returned")
    parser.add_argument('-pre', '--pull_record_ends', nargs=2, metavar=('<amount (int)>', "<front|rear>"),
                        action='store', help="Get the ends (front or rear) of all sequences in a file.")
    parser.add_argument('-er', '--extract_region', action='store', nargs=2, metavar=("<start (int)>", "<end (int)>"),
                        type=int, help="Pull out sub-sequences.")
    parser.add_argument('-dr', '--delete_records', action='store', nargs="+", metavar="<regex pattern>",  # ToDo: update wiki to show multiple regex inputs
                        help="Remove records from a file. The deleted IDs are sent to stderr.")
    parser.add_argument('-dsm', '--delete_small', help='Delete sequences with length below threshold', type=int,
                        metavar='<threshold (int)>', action='store')
    parser.add_argument('-dlg', '--delete_large', help='Delete sequences with length above threshold', type=int,
                        metavar='<threshold (int)>', action='store')
    parser.add_argument('-df', '--delete_features', action='store', nargs="+", metavar="<regex pattern>",  # ToDo: update wiki to show multiple regex inputs
                        help="Remove specified features from all records.")
    parser.add_argument('-drp', '--delete_repeats', action='append', nargs="?", type=int,
                        help="Strip repeat records (ids and/or identical sequences. Optional, pass in an integer to "
                             "specify # columns for deleted ids")
    parser.add_argument('-frp', '--find_repeats', action='append', nargs="?", type=int,
                        help="Identify whether a file contains repeat sequences and/or sequence ids. The number of "
                             "output columns can be modified by passing in an integer.")
    parser.add_argument("-mg", "--merge", action="store_true",
                        help="Group a bunch of seq files together",)
    parser.add_argument("-bl", "--blast", metavar="<BLAST database>", action="store",
                        help="BLAST your sequence file using common settings, return the hits from blastdb")
    parser.add_argument("-bl2s", "--bl2seq", action="store_true",
                        help="All-by-all blast among sequences using bl2seq. Only Returns top hit from each search")
    parser.add_argument("-prg", "--purge", metavar="Max BLAST score", type=int, action="store",
                        help="Delete sequences with high similarity")
    parser.add_argument("-sbt", "--split_by_taxa", action='store', nargs=2, metavar=("<Split Pattern>", "<out dir>"),
                        help="")
    parser.add_argument("-frs", "--find_restriction_sites", action='store', nargs=2,
                        metavar=("<commercial>", "<num_cuts>"), help="")
    parser.add_argument("-fp", "--find_pattern", action='store', nargs="+", metavar="<regex pattern(s)>",
                        help="Search for subsequences, returning the start positions of all matches.")
    parser.add_argument("-mw", "--molecular_weight", action='store_true', help="")
    parser.add_argument("-ip", "--isoelectric_point", action='store_true', help="")
    parser.add_argument("-fcpg", "--find_CpG", action='store_true', help="")
    parser.add_argument("-cr", "--count_residues", action='store_true', help="")
    parser.add_argument("-spf", "--split_file", action='store', help="")
    parser.add_argument('-ga', '--guess_alphabet', action='store_true')
    parser.add_argument('-gf', '--guess_format', action='store_true')
    parser.add_argument('-ho', '--hash_output', help="For ")
    parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')
    parser.add_argument('-p', '--params', help="Free form arguments for some functions", nargs="+", action='store')
    parser.add_argument('-q', '--quiet', help="Suppress stderr messages", action='store_true')
    parser.add_argument('-t', '--test', action='store_true',
                        help="Run the function and return any stderr/stdout other than sequences.")
    parser.add_argument('-o', '--out_format', help="If you want a specific format output", action='store')
    parser.add_argument('-f', '--in_format', help="If SeqBuddy can't guess the file format, just specify it directly.",
                        action='store')
    parser.add_argument('-a', '--alpha', help="If you want the file read with a specific alphabet", action='store')

    in_args = parser.parse_args()

    seqbuddy = []
    seq_set = ""

    try:
        for seq_set in in_args.sequence:
            seq_set = SeqBuddy(seq_set, in_args.in_format, in_args.out_format, in_args.alpha)
            seqbuddy += seq_set.records

        seqbuddy = SeqBuddy(seqbuddy, seq_set.in_format, seq_set.out_format, seq_set.alpha)
    except GuessError:
        sys.exit("Error: SeqBuddy could not understand your input. "
                 "Check the file path or try specifying an input type with -f")

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
            for param in in_args.params:
                binary = Popen("%s -version" % param, stdout=PIPE, shell=True).communicate()
                binary = re.search("([a-z])*[^:]", binary[0].decode("utf-8"))
                binary = binary.group(0)
                if binary == "blastp":
                    blastp = param
                elif binary == "blastn":
                    blastn = param
                elif binary == "blastdbcmd":
                    blastdbcmd = param

        blastp = blastp if blastp else which("blastp")
        blastn = blastn if blastn else which("blastn")
        blastdbcmd = blastdbcmd if blastdbcmd else which("blastdbcmd")

        return {"blastdbcmd": blastdbcmd, "blastp": blastp, "blastn": blastn}


    def _raise_error(_err):
        _stderr("{0}: {1}\n".format(_err.__class__.__name__, str(_err)))


    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Purge
    if in_args.purge:
        purged_seqs, deleted, record_map = purge(seqbuddy, in_args.purge)
        _stderr(record_map, in_args.quiet)
        _print_recs(purged_seqs)

    # BL2SEQ
    if in_args.bl2seq:
        output = bl2seq(seqbuddy)
        _stdout(output[1])

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

    # Shuffle
    if in_args.shuffle:
        _print_recs(shuffle(seqbuddy))

    # Order ids
    if in_args.order_ids:
        reverse = True if in_args.order_ids[0] and in_args.order_ids[0] == "rev" else False
        _print_recs(order_ids(seqbuddy, _reverse=reverse))

    # Find repeat sequences or ids
    if in_args.find_repeats:
        columns = 1 if not in_args.find_repeats[0] else in_args.find_repeats[0]
        unique, rep_ids, rep_seqs, out_string = find_repeats(seqbuddy, columns)
        _stdout(out_string)

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

    # Delete sequences above threshold
    if in_args.delete_large:
        _print_recs(delete_large(seqbuddy, in_args.delete_large))

    # Delete sequences below threshold
    if in_args.delete_small:
        _print_recs(delete_small(seqbuddy, in_args.delete_small))

    # Delete features
    if in_args.delete_features:
        for next_pattern in in_args.delete_features:
            delete_features(seqbuddy, next_pattern)
        _print_recs(seqbuddy)

    # Merge
    if in_args.merge:
        _print_recs(seqbuddy)

    # Screw formats
    if in_args.screw_formats:
        seqbuddy.out_format = in_args.screw_formats
        if in_args.in_place:  # Need to change the file extension
            os.remove(in_args.sequence[0])
            in_args.sequence[0] = ".".join(os.path.abspath(in_args.sequence[0]).split(".")[:-1]) + \
                                  "." + seqbuddy.out_format
            open(in_args.sequence[0], "w").close()
        _print_recs(seqbuddy)

    # Renaming
    if in_args.rename_ids:
        num = 0 if not in_args.params else int(in_args.params[0])
        _print_recs(rename(seqbuddy, in_args.rename_ids[0], in_args.rename_ids[1], num))

    # Uppercase
    if in_args.uppercase:
        _print_recs(uppercase(seqbuddy))

    # Lowercase
    if in_args.lowercase:
        _print_recs(lowercase(seqbuddy))

    # Transcribe
    if in_args.transcribe:
        if seqbuddy.alpha != IUPAC.ambiguous_dna:
            raise ValueError("You need to provide a DNA sequence.")
        _print_recs(dna2rna(seqbuddy))

    # Back Transcribe
    if in_args.back_transcribe:
        if seqbuddy.alpha != IUPAC.ambiguous_rna:
            raise ValueError("You need to provide an RNA sequence.")
        _print_recs(rna2dna(seqbuddy))

    # Complement
    if in_args.complement:
        _print_recs(complement(seqbuddy))

    # Reverse complement
    if in_args.reverse_complement:
        _print_recs(reverse_complement(seqbuddy))

    # List identifiers
    if in_args.list_ids:
        columns = 1 if not in_args.list_ids[0] else in_args.list_ids[0]
        _stdout(list_ids(seqbuddy, columns))

    # Translate CDS
    if in_args.translate:
        if seqbuddy.alpha == IUPAC.protein:
            raise ValueError("You need to supply DNA or RNA sequences to translate")

        _print_recs(translate_cds(seqbuddy, quiet=in_args.quiet))

    # Shift reading frame
    if in_args.select_frame:
        _print_recs(select_frame(seqbuddy, in_args.select_frame))

    # Translate 6 reading frames
    if in_args.translate6frames:
        if seqbuddy.alpha == IUPAC.protein:
            raise ValueError("You need to supply DNA or RNA sequences to translate")

        seqbuddy = translate6frames(seqbuddy)
        if in_args.out_format:
            seqbuddy.out_format = in_args.out_format

        _print_recs(seqbuddy)

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

    # Concatenate sequences
    if in_args.concat_seqs:
        clean = False if not in_args.concat_seqs[0] or in_args.concat_seqs[0] != "clean" else True
        seqbuddy = concat_seqs(seqbuddy, clean)
        if in_args.out_format:
            seqbuddy.out_format = in_args.out_format
        _print_recs(seqbuddy)

    # Count number of sequences in a file
    if in_args.num_seqs:
        _stdout("%s\n" % num_seqs(seqbuddy))

    # Average length of sequences
    if in_args.ave_seq_length:
        clean = False if not in_args.ave_seq_length[0] or in_args.ave_seq_length[0] != "clean" else True
        _stdout("%s\n" % round(ave_seq_length(seqbuddy, clean), 2))

    # Pull sequence ends
    if in_args.pull_record_ends:
        _print_recs(pull_record_ends(seqbuddy, *in_args.pull_record_ends))

    # Extract regions
    if in_args.extract_region:
        _print_recs(extract_range(seqbuddy, *in_args.extract_region))

    # Pull records
    if in_args.pull_records:
        _print_recs(pull_recs(seqbuddy, in_args.pull_records))

    # Pull random records
    if in_args.pull_random_record:
        count = 1 if not in_args.pull_random_record[0] else in_args.pull_random_record[0]
        _print_recs(pull_random_recs(seqbuddy, count))

    # Hash sequence ids
    if in_args.hash_seq_ids:
        hash_length = in_args.hash_seq_ids[0] if in_args.hash_seq_ids[0] else 10
        seqbuddy, hash_map, hash_table = hash_sequence_ids(seqbuddy, hash_length)
        _stderr(hash_table, in_args.quiet)
        _print_recs(seqbuddy)

    # Delete metadata
    if in_args.delete_metadata:
        _print_recs(delete_metadata(seqbuddy))

    # Raw Seq
    if in_args.raw_seq:
        output = raw_seq(seqbuddy)
        if in_args.in_place:
            _in_place(output, in_args.sequence[0])
        else:
            _stdout(output)

    # Clean Seq
    if in_args.clean_seq:
        if in_args.clean_seq[0] == "strict":
            _print_recs(clean_seq(seqbuddy, ambiguous=False))
        else:
            _print_recs(clean_seq(seqbuddy, ambiguous=True))

    # Guess format
    if in_args.guess_format:
        if str(type(in_args.sequence[0])) == "<class '_io.TextIOWrapper'>":
            _stdout("{0}\n".format(seqbuddy.in_format))
        else:
            for seq_set in in_args.sequence:
                _stdout("%s\t-->\t%s\n" % (seq_set, SeqBuddy(seq_set).in_format))

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

    # Combine feature sets from two files into one
    if in_args.combine_features:
        file1, file2 = in_args.sequence[:2]
        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)
        _print_recs(combine_features(file1, file2))

    # Order sequence features by their position in the sequence
    if in_args.order_features_by_position:
        reverse = True if in_args.order_features_by_position[0] \
            and in_args.order_features_by_position[0] == "rev" else False

        _print_recs(order_features_by_position(seqbuddy, reverse))

    # Order sequence features alphabetically
    if in_args.order_features_alphabetically:
        reverse = True if in_args.order_features_alphabetically[0] \
            and in_args.order_features_alphabetically[0] == "rev" else False
        _print_recs(order_features_alphabetically(seqbuddy, reverse))

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

    # Calculate Molecular Weight
    if in_args.molecular_weight:
        lists = molecular_weight(seqbuddy)
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

    # Calculate Isoelectric Point
    if in_args.isoelectric_point:
        try:
            isoelectric_points = isoelectric_point(seqbuddy)
            _stderr("ID\tpI\n")
            for pI in isoelectric_points:
                print("{0}\t{1}".format(pI[0], pI[1]))
        except ValueError as e:
            _raise_error(e)

    # Count residues
    if in_args.count_residues:
        output = count_residues(seqbuddy)
        for sequence in output:
            print(sequence[0])
            if seqbuddy.alpha in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna, IUPAC.ambiguous_rna,
                                  IUPAC.unambiguous_rna]:
                for residue in ['% Ambiguous', 'A', 'T', 'C', 'G', 'U']:
                    print("{0}: {1}".format(residue, sequence[1].pop(residue, None)))
            for residue in (sequence[1]):
                print("{0}: {1}".format(residue, sequence[1][residue]))
            print()

    # Split file
    if in_args.split_file:
        in_args.in_place = True
        out_dir = os.path.abspath(in_args.split_file)
        os.makedirs(out_dir, exist_ok=True)
        taxa_groups = split_by_taxa(seqbuddy, in_args.split_file)
        check_quiet = in_args.quiet  # 'quiet' must be toggled to 'on' _print_recs() here.
        in_args.quiet = True
        for buddy in split_file(seqbuddy):
            seqbuddy.records = buddy.records
            ext = _format_to_extension(seqbuddy.out_format)
            in_args.sequence[0] = "%s/%s.%s" % (out_dir, seqbuddy.records[0].id, ext)
            _stderr("New file: %s\n" % in_args.sequence[0], check_quiet)
            open(in_args.sequence[0], "w").close()
            _print_recs(seqbuddy)

    # Find restriction sites
    if in_args.find_restriction_sites:
        if in_args.find_restriction_sites[0] in ['commercial', 'com', 'c']:
            commercial = True
        elif in_args.find_restriction_sites[0] in ['all', 'a', 'noncommercial', 'nc', 'noncom']:
            commercial = False
        else:
            commercial = True
        if in_args.find_restriction_sites[1] in ['1', 'single', 's']:
            single_cut = True
        elif in_args.find_restriction_sites[1] in ['2', 'double', 'd']:
            single_cut = False
        else:
            single_cut = True
        sb = clean_seq(seqbuddy)
        output = find_restriction_sites(sb, commercial, single_cut)
        _stdout(output[1])

    # Find pattern
    if in_args.find_pattern:
        for pattern in in_args.find_pattern:
            output = find_pattern(seqbuddy, pattern)
            out_string = ""
            num_matches = 0
            for key in output:
                out_string += "{0}: {1}\n".format(key, ", ".join([str(x) for x in output[key]]))
                num_matches += len(output[key])
            _stderr("#### {0} matches found across {1} sequences for "
                    "pattern '{2}' ####\n".format(num_matches, len(output), pattern))
            _stdout("%s\n" % out_string)

    # Find CpG
    if in_args.find_CpG:
        output = find_CpG(seqbuddy)
        out_string = ""
        for key in output[1]:
            out_string += "{0}: {1}\n".format(key, str(output[1][key]))
        print(output[0])
        _stderr('\n'+out_string, in_args.quiet)

