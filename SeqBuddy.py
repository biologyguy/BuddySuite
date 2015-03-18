#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Nov 20 2014 
# 43 tools and counting

"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, version 2 of the License (GPLv2).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details at http://www.gnu.org/licenses/.

name: SeqBuddy.py
date: Nov-20-2014
version: 0, unstable
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

# Third party package imports
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError

# My functions
from MyFuncs import run_multicore_function


# ##################################################### WISH LIST #################################################### #
def sim_ident():  # Return the pairwise similarity and identity scores among sequences
    x = 1
    return x


# - Add FASTQ support... More generally, support letter annotation mods
# - Check on memory requirements before execution
# - Execution timer, for long running jobs
# - Implement daisy chaining
# - Handle all stderr output from private function (to allow quiet execution)
# - Unit Tests
# - Allow batch calls. E.g., if 6 files are fed in as input, run the SeqBuddy command provided independently on each
# - Flag that will suppress stdout, while still outputting stderr
# ################################################# HELPER FUNCTIONS ################################################# #
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
# #################################################################################################################### #


class SeqBuddy():  # Open a file or read a handle and parse, or convert raw into a Seq object
    def __init__(self, _input, _in_format=None, _out_format=None):
        if not _in_format:
            self.in_format = guess_format(_input)
            self.out_format = str(self.in_format) if not _out_format else _out_format

        else:
            self.in_format = _in_format

        if not self.in_format:
            raise AttributeError("Could not determine sequence format from _input. "
                                 "Try explicitly setting with -f flag.")

        self.out_format = self.in_format if not _out_format else _out_format

        if str(type(_input)) == "<class '__main__.SeqBuddy'>":
            _sequences = _input.records

        elif isinstance(_input, list):
            # make sure that the list is actually SeqIO records (just test a few...)
            for _seq in _input[:3]:
                if type(_seq) != SeqRecord:
                    raise TypeError("Seqlist is not populated with SeqRecords.")
            _sequences = _input

        elif str(type(_input)) == "<class '_io.TextIOWrapper'>":
            _sequences = list(SeqIO.parse(_input, self.in_format))

        elif os.path.isfile(_input):
            _input = open(_input, "r")
            _sequences = list(SeqIO.parse(_input, self.in_format))
            _input.close()
        else:
            _sequences = [SeqRecord(Seq(_input))]

        self.alpha = guess_alphabet(_sequences)

        for _i in range(len(_sequences)):
            _sequences[_i].seq.alphabet = self.alpha

        self.records = _sequences

    def to_dict(self):
        _unique, _rep_ids, _rep_seqs = find_repeats(self)
        if len(_rep_ids) > 0:
            raise RuntimeError("There are repeat IDs in self.records\n%s" % _rep_ids)

        records_dict = {}
        for _rec in self.records:
            records_dict[_rec.id] = _rec
        return records_dict


def guess_alphabet(_seqbuddy):  # Does not handle ambiguous dna
    _seq_list = _seqbuddy if isinstance(_seqbuddy, list) else _seqbuddy.records
    _sequence = ""
    for next_seq in _seq_list:
        if len(_sequence) > 1000:
            break
        _sequence += re.sub("[NX\-?]", "", str(next_seq.seq))
        _sequence = _sequence.upper()

    if len(_sequence) == 0:
        return None
    percent_dna = float(_sequence.count("A") + _sequence.count("G") + _sequence.count("T") +
                        _sequence.count("C") + _sequence.count("U")) / float(len(_sequence))
    if percent_dna > 0.95:
        nuc = IUPAC.ambiguous_rna if float(_sequence.count("U")) / float(len(_sequence)) > 0.05 else IUPAC.ambiguous_dna
        return nuc
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

    if str(type(_input)) == "<class '_io.TextIOWrapper'>":
        # Die if file is empty
        if _input.read() == "":
            sys.exit("_input file is empty.")
        _input.seek(0)

        possible_formats = ["phylip-relaxed", "stockholm", "fasta", "gb", "fastq", "nexus"]
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
        raise AttributeError("Unsupported _input argument in guess_format(). %s" % _input)


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
    with open("%s/tmp.fa" % tmp_dir.name, "w") as ofile:
        SeqIO.write(_seqbuddy.records, ofile, "fasta")

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

    ofile = open("%s/seqs.fa" % tmp_dir.name, "w")
    for hit_id in hit_ids:
        hit = Popen("blastdbcmd -db %s -entry 'lcl|%s'" % (blast_db, hit_id), stdout=PIPE, shell=True).communicate()
        hit = hit[0].decode("utf-8")
        hit = re.sub("lcl\|", "", hit)
        ofile.write("%s\n" % hit)

    ofile.close()

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
    _output = []
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq.back_transcribe()), alphabet=IUPAC.ambiguous_dna)
        _output.append(_rec)
    _seqbuddy.records = _output
    return _seqbuddy


def dna2rna(_seqbuddy):
    _output = []
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq.transcribe()), alphabet=IUPAC.ambiguous_rna)
        _output.append(_rec)
    _seqbuddy.records = _output
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


def translate_cds(_seqbuddy, quiet=False):  # adding 'quiet' will suppress the errors thrown by translate(cds=True)
    def trans(in_seq):
        try:
            in_seq.seq = in_seq.seq.translate(cds=True, to_stop=True)
            return in_seq

        except TranslationError as _e1:
            if not quiet:
                sys.stderr.write("Warning: %s in %s\n" % (_e1, in_seq.id))
            return _e1

    _translation = deepcopy(_seqbuddy)
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
                temp_seq.seq = Seq(re.sub(regex[0], "NNN", str(temp_seq.seq), count=1), alphabet=temp_seq.seq.alphabet)
                _rec.seq = Seq(re.sub(regex[0], "NNN", str(_rec.seq), count=1), alphabet=_rec.seq.alphabet)
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
        raise AttributeError("The input sequence needs to be protein, not %s" % _seqbuddy.alpha)

    if not _species:
        lookup_table = rand_table
    elif _species.upper() in ["HUMAN", "H"]:
        lookup_table = h_sapi
    elif _species.upper() in ["MOUSE", "M"]:
        lookup_table = m_muscul
        print("mouse")
    elif _species.upper() in ["ECOLI", "E"]:
        lookup_table = e_coli
        print("ecoli")
    elif _species.upper() in ["YEAST", "Y"]:
        lookup_table = s_cerev
        print(_species)
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
        for aa in _rec.seq:
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


def clean_seq(_seqbuddy):
    """remove all non-sequence charcters from sequence strings"""
    _output = []
    for _rec in _seqbuddy.records:
        _rec.seq = str(_rec.seq).upper()
        if _seqbuddy.alpha == IUPAC.protein:
            _rec.seq = Seq(re.sub("[^ACDEFGHIKLMNPQRSTVWXY]", "", str(_rec.seq)), alphabet=_seqbuddy.alpha)
        else:
            _rec.seq = Seq(re.sub("[^ATGCXNU]", "", str(_rec.seq)), alphabet=_seqbuddy.alpha)

        _output.append(_rec)

    _seqbuddy.records = _output
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

    prot_dict = SeqIO.to_dict(prot_seqbuddy.records)
    dna_dict = SeqIO.to_dict(dna_seqbuddy.records)
    _new_seqs = {}
    stderr_written = False
    for _seq_id in dna_dict:
        if _seq_id not in prot_dict:
            stderr_written = True
            sys.stderr.write("Warning: %s is in the cDNA file, but not in the protein file\n" % _seq_id)
            continue

        _new_seqs[_seq_id] = prot_dict[_seq_id]
        for feature in dna_dict[_seq_id].features:
            prot_dict[_seq_id].features.append(_feature_map(feature))

    for _seq_id in prot_dict:
        if _seq_id not in dna_dict:
            stderr_written = True
            sys.stderr.write("Warning: %s is in the protein file, but not in the cDNA file\n" % _seq_id)

    if stderr_written:
        sys.stderr.write("\n")

    _seqs_list = [_new_seqs[_seq_id] for _seq_id in _new_seqs]
    _seqbuddy = SeqBuddy(_seqs_list)
    _seqbuddy.out_format = "gb"
    return _seqbuddy


# Apply DNA features to protein sequences
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

    prot_dict = SeqIO.to_dict(prot_seqbuddy.records)
    dna_dict = SeqIO.to_dict(dna_seqbuddy.records)
    _new_seqs = {}
    stderr_written = False
    for _seq_id in prot_dict:
        if _seq_id not in dna_dict:
            stderr_written = True
            sys.stderr.write("Warning: %s is in the protein file, but not in the cDNA file\n" % _seq_id)
            continue

        _new_seqs[_seq_id] = dna_dict[_seq_id]
        for feature in prot_dict[_seq_id].features:
            dna_dict[_seq_id].features.append(_feature_map(feature))

    for _seq_id in dna_dict:
        if _seq_id not in prot_dict:
            stderr_written = True
            sys.stderr.write("Warning: %s is in the cDNA file, but not in the protein file\n" % _seq_id)

    if stderr_written:
        sys.stderr.write("\n")

    _seqs_list = [_new_seqs[_seq_id] for _seq_id in _new_seqs]
    _seqbuddy = SeqBuddy(_seqs_list)
    _seqbuddy.out_format = "gb"
    return _seqbuddy


# Merge feature lists
def combine_features(_seqbuddy1, _seqbuddy2):  # ToDo: rewrite this to accept any number of input files.
    # make sure there are no repeat ids
    _unique, _rep_ids, _rep_seqs = find_repeats(_seqbuddy1)
    if len(_rep_ids) > 0:
        raise RuntimeError("There are repeat IDs in the first file provided\n%s" % _rep_ids)

    _unique, _rep_ids, _rep_seqs = find_repeats(_seqbuddy2)
    if len(_rep_ids) > 0:
        raise RuntimeError("There are repeat IDs in the second file provided\n%s" % _rep_ids)

    seq_dict1 = {}
    seq_dict2 = {}

    for _rec in _seqbuddy1.records:
        seq_dict1[_rec.id] = _rec

    for _rec in _seqbuddy2.records:
        seq_dict2[_rec.id] = _rec

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

    _seqbuddy = SeqBuddy([_new_seqs[_seq_id] for _seq_id in _new_seqs], _out_format=_seqbuddy1.in_format)
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


def hash_seqeunce_ids(_seqbuddy, _hash_length=10):
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

    hash_map = []
    for i in range(len(hash_list)):
        hash_map.append((hash_list[i], seq_ids[i]))

    return [hash_map, _seqbuddy]


def pull_recs(_seqbuddy, _search):
    _output = []
    for _rec in _seqbuddy.records:
        if re.search(_search, _rec.description) or re.search(_search, _rec.id) or re.search(_search, _rec.name):
            _output.append(_rec)
    _seqbuddy.records = _output
    return _seqbuddy


def pull_random_recs(_seqbuddy, _count=1):  # Return a random set of sequences (without replacement)
    _output = []
    _count = abs(_count) if abs(_count) <= len(_seqbuddy.records) else len(_seqbuddy.records)
    for i in range(_count):
        rand_index = randint(0, len(_seqbuddy.records) - 1)
        _output.append(_seqbuddy.records.pop(rand_index))

    _seqbuddy.records = _output
    return _seqbuddy


def pull_record_ends(_seqbuddy, _amount, _which_end):
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
        raise AttributeError("Error at extract range: The value given for end of range is smaller than for the start "
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


def find_repeats(_seqbuddy):
    unique_seqs = {}
    repeat_ids = {}
    repeat_seqs = {}

    # First find replicate IDs
    for _rec in _seqbuddy.records:
        if _rec.id in repeat_ids:
            repeat_ids[_rec.id].append(_rec)
        elif _rec.id in unique_seqs:
            repeat_ids[_rec.id] = [_rec]
            repeat_ids[_rec.id].append(unique_seqs[_rec.id])
            del(unique_seqs[_rec.id])
        else:
            unique_seqs[_rec.id] = _rec

    # Then look for replicate sequences ToDo: use MD5 hashes instead of whole sequences as keys
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
            del(unique_seqs[key])

    for key, value in repeat_ids.items():  # find duplicates in the repeat ID list
        for _rep_seq in value:
            _rep_seq = str(_rep_seq.seq)
            if _rep_seq not in flip_uniqe:
                flip_uniqe[_rep_seq] = [key]
            else:
                if _rep_seq not in repeat_seqs:
                    repeat_seqs[_rep_seq] = [key]
                    repeat_seqs[_rep_seq] += flip_uniqe[_rep_seq]

                else:
                    repeat_seqs[_rep_seq].append(key)
    return [unique_seqs, repeat_ids, repeat_seqs]


def delete_records(_seqbuddy, search_str):
    _output = []
    _deleted = pull_recs(copy(_seqbuddy), search_str).records
    for _rec in _seqbuddy.records:
        if _rec in _deleted:
            continue
        else:
            _output.append(_rec)
    _seqbuddy.records = _output
    return _seqbuddy


def delete_large(_seqbuddy, max_value):
    _output = []
    for _rec in _seqbuddy.records:
        if len(str(_rec.seq)) <= max_value:
            _output.append(_rec)
    _seqbuddy.records = _output
    return _seqbuddy


def delete_small(_seqbuddy, min_value):
    _output = []
    for _rec in _seqbuddy.records:
        if len(str(_rec.seq)) >= min_value:
            _output.append(_rec)
    _seqbuddy.records = _output
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
        _unique, _rep_ids, _rep_seqs = find_repeats(_seqbuddy)
        if len(_rep_ids) > 0:
            for _rep_id in _rep_ids:
                store_one_copy = pull_recs(copy(_seqbuddy), "^%s$" % _rep_id).records[0]
                delete_records(_seqbuddy, "^%s$" % _rep_id)
                _seqbuddy.records.append(store_one_copy)

    # Then remove duplicate sequences
    if scope in ['all', 'seqs']:
        _unique, _rep_ids, _rep_seqs = find_repeats(_seqbuddy)
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


def purge(_seqbuddy, threshold):  # ToDo: Implement a way to return a certain # of sequences (i.e. auto-determine threshold)
    keep_set = {}
    purged = []
    _blast_res = bl2seq(_seqbuddy)
    for _query_id in _blast_res:
        if _query_id in purged:
            continue
        else:
            keep_set[_query_id] = []
            for _subj_id in _blast_res[_query_id]:
                _ident, _length, _evalue, _bit_score = _blast_res[_query_id][_subj_id]

                if _bit_score >= threshold:
                    purged.append(_subj_id)
                    keep_set[_query_id].append(_subj_id)

    _output = []

    for _rec in _seqbuddy.records:
        if _rec.id in keep_set:
            _output.append(_rec)

    _seqbuddy.records = _output
    return [_seqbuddy, keep_set]


def bl2seq(_seqbuddy, cores=4):  # Does an all-by-all analysis, and does not return sequences
    """
    Note on blast2seq: Expect (E) values are calculated on an assumed database size of (the rather large) nr, so the
    threshold may need to be increased quite a bit to return short alignments
    """
    def mc_blast(_query, args):
        _values, _subject = args

        if subject.id == _query.id:
            return

        _blast_res = Popen("echo '%s' | %s -subject %s -outfmt 6" %
                           (_query.format("fasta"), blast_bin, subject_file), stdout=PIPE, shell=True).communicate()
        _blast_res = _blast_res[0].decode().split("\n")[0].split("\t")

        while True:
            indx = randint(0, len(_values) - 1)
            try:
                if len(_blast_res) == 1:
                    _values[indx].value += "%s\t%s\t0\t0\t0\t0\n" % (subject.id, _query.id)
                else:
                    # values are: query, subject, %_ident, length, evalue, bit_score
                    if _blast_res[10] == '0.0':
                        _blast_res[10] = '1e-180'
                    _values[indx].value += "%s\t%s\t%s\t%s\t%s\t%s\n" % (_blast_res[0], _blast_res[1], _blast_res[2],
                                                                         _blast_res[3], _blast_res[10],
                                                                         _blast_res[11].strip())
                break
            except ConnectionRefusedError:
                continue
        return

    blast_bin = "blastp" if _seqbuddy.alpha == IUPAC.protein else "blastn"
    if not which(blast_bin):
        raise RuntimeError("%s not present in $PATH." % blast_bin)  # ToDo: Implement -p flag

    tmp_dir = TemporaryDirectory()

    # Prepare to store multicore output in list of Value objects
    values = []
    manager = Manager()
    for _ in range(cores * 2):
        values.append(manager.Value(ctypes.c_char_p, ''))

    # Copy the seqbuddy records into new list, so they can be iteratively deleted below
    _seqs_copy = _seqbuddy.records[:]
    subject_file = "%s/subject.fa" % tmp_dir.name
    _output = ''
    for subject in _seqbuddy.records:
        with open(subject_file, "w") as ifile:
            SeqIO.write(subject, ifile, "fasta")

        run_multicore_function(_seqs_copy, mc_blast, [values, subject_file], out_type=sys.stderr, quiet=True)

        for i in range(len(values)):
            _output += values[i].value
            values[i].value = ''

        _seqs_copy = _seqs_copy[1:]

    # Push output into a dictionary of dictionaries, for more flexible use outside of this function
    output_list = _output.strip().split("\n")
    output_list = [x.split("\t") for x in output_list]
    output_dir = {}
    for match in output_list:
        query, subj, _ident, _length, _evalue, _bit_score = match
        if query not in output_dir:
            output_dir[query] = {subj: [float(_ident), int(_length), float(_evalue), int(_bit_score)]}
        else:
            output_dir[query][subj] = [float(_ident), int(_length), float(_evalue), int(_bit_score)]

        if subj not in output_dir:
            output_dir[subj] = {query: [float(_ident), int(_length), float(_evalue), int(_bit_score)]}
        else:
            output_dir[subj][query] = [float(_ident), int(_length), float(_evalue), int(_bit_score)]

    return output_dir


def uppercase(_seqbuddy):
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq).upper(), alphabet=_rec.seq.alphabet)
    return _seqbuddy


def lowercase(_seqbuddy):
    for _rec in _seqbuddy.records:
        _rec.seq = Seq(str(_rec.seq).lower(), alphabet=_rec.seq.alphabet)
    return _seqbuddy


# ################################################# COMMAND LINE UI ################################################## #
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="SeqBuddy.py", description="Commandline wrapper for all the fun functions in "
                                                                     "this file. Play with your sequences!")

    parser.add_argument("sequence", help="Supply a file path or a raw sequence", nargs="+", default=sys.stdin)

    parser.add_argument('-cs', '--clean_seq', action='store_true',
                        help="Strip out non-sequence characters, such as stops (*) and gaps (-)")
    parser.add_argument('-uc', '--uppercase', action='store_true',
                        help='Convert all sequences to uppercase')
    parser.add_argument('-lc', '--lowercase', action='store_true',
                        help='Convert all sequences to lowercase')
    parser.add_argument('-dm', '--delete_metadata', action='store_true',
                        help="Remove meta-data from file (only id is retained)")
    parser.add_argument('-rs', '--raw_seq', action='store_true',
                        help="Return line break separated sequences")
    parser.add_argument('-tr', '--translate', action='store_true',
                        help="Convert coding sequences into amino acid sequences")
    parser.add_argument('-sfr', '--select_frame', action='store', metavar='<frame (int)>', type=int, choices=[1, 2, 3],
                        help="Change the reading from of sequences by deleting characters off of the front")
    parser.add_argument('-tr6', '--translate6frames', action='store_true',
                        help="Translate nucleotide sequences into all six reading frames")
    parser.add_argument('-btr', '--back_translate', action='store_true',
                        help="Convert amino acid sequences into codons. Select mode and species with -p flag "
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
    parser.add_argument('-li', '--list_ids', action='store_true',
                        help="Output all the sequence identifiers in a file. Use -p to specify # columns to write")
    parser.add_argument('-ns', '--num_seqs', action='store_true',
                        help="Counts how many sequences are present in an input file")
    parser.add_argument('-asl', '--ave_seq_length', action='store_true',
                        help="Return the average length of all sequences. Use '-p clean' to remove gaps etc from the "
                             "sequences before counting.")
    parser.add_argument('-cts', '--concat_seqs', action='store_true',
                        help="Concatenate a bunch of sequences into a single solid string. Use '-p clean' to remove "
                             "stops, gaps, etc., from the sequences before concatenating.")
    parser.add_argument('-fd2p', '--map_features_dna2prot', action='store_true',
                        help="Take the features annotated onto nucleotide sequences and map to protein sequences. "
                             "Both a protein and cDNA file must be passed in.")
    parser.add_argument('-fp2d', '--map_features_prot2dna', action='store_true',
                        help="Arguments: one cDNA file and one protein file")
    parser.add_argument('-ri', '--rename_ids', action='store', metavar=('<pattern>', '<substitution>'), nargs=2,
                        help="Replace some pattern in ids with something else.")
    parser.add_argument('-cf', '--combine_features', action='store_true',
                        help="Takes the features in two files and combines them for each sequence")
    parser.add_argument('-sh', '--shuffle', action='store_true',
                        help="Randomly reorder the position of records in the file.")
    parser.add_argument('-oi', '--order_ids', action='store_true',
                        help="Sort all sequences by id in alpha-numeric order. Use -p 'rev' for reverse order")
    parser.add_argument('-ofp', '--order_features_by_position', action='store_true',
                        help="Change the output order of sequence features, based on sequence position")
    parser.add_argument('-ofa', '--order_features_alphabetically', action='store_true',
                        help="Change the output order of sequence features, based on sequence position")
    parser.add_argument('-sf', '--screw_formats', action='store', metavar="<out_format>",
                        help="Change the file format to something else.")
    parser.add_argument('-hsi', '--hash_seq_ids', action='store_true',
                        help="Rename all the identifiers in a sequence list to a 10 character hash.")
    parser.add_argument('-pr', '--pull_records', action='store', metavar="<regex pattern>",
                        help="Get all the records with ids containing a given string")
    parser.add_argument('-prr', '--pull_random_record', action='store_true',
                        help="Extract random sequences. Use the -p flag to increase the number of sequences returned")
    parser.add_argument('-pre', '--pull_record_ends', action='store', nargs=2, metavar="<amount (int)> <front|rear>",
                        help="Get the ends (front or rear) of all sequences in a file.")
    parser.add_argument('-er', '--extract_region', action='store', nargs=2, metavar="<start (int)> <end (int)>",
                        type=int, help="Pull out sub-sequences.")
    parser.add_argument('-dr', '--delete_records', action='store', nargs="+", metavar="<regex pattern>",
                        help="Remove records from a file. The deleted IDs are sent to stderr.")
    parser.add_argument('-dsm', '--delete_small', help='Delete sequences with length below threshold', type=int,
                        action='store')
    parser.add_argument('-dlg', '--delete_large', help='Delete sequences with length above threshold', type=int,
                        action='store')
    parser.add_argument('-df', '--delete_features', action='store', nargs="+", metavar="<regex pattern>",
                        help="Remove specified features from all records.")
    parser.add_argument('-drp', '--delete_repeats', action='store_true',
                        help="Strip repeat records (ids and/or identical sequences")
    parser.add_argument('-frp', '--find_repeats', action='store_true',
                        help="Identify whether a file contains repeat sequences and/or sequence ids")
    parser.add_argument("-mg", "--merge", action="store_true",
                        help="Group a bunch of seq files together",)
    parser.add_argument("-bl", "--blast", metavar="<BLAST database>", action="store",
                        help="BLAST your sequence file using common settings, return the hits from blastdb")
    parser.add_argument("-bl2s", "--bl2seq", action="store_true",
                        help="All-by-all blast among sequences using bl2seq. Only Returns top hit from each search")
    parser.add_argument("-prg", "--purge", metavar="Max BLAST score", type=int, action="store",
                        help="Delete sequences with high similarity")
    parser.add_argument('-ga', '--guess_alphabet', action='store_true')
    parser.add_argument('-gf', '--guess_format', action='store_true')

    parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')
    parser.add_argument('-p', '--params', help="Free form arguments for some functions", nargs="+", action='store')
    parser.add_argument('-q', '--quiet', help="Suppress stderr messages", action='store_true')  # ToDo: implement this everywhere
    parser.add_argument('-o', '--out_format', help="If you want a specific format output", action='store')
    parser.add_argument('-f', '--in_format', help="If SeqBuddy can't guess the file format, just specify it directly.",
                        action='store')
    
    in_args = parser.parse_args()

    in_place_allowed = True  # This might be deprecated, because no tools that call _print_recs() are False

    seqbuddy = []
    seq_set = ""
    for seq_set in in_args.sequence:
        seq_set = SeqBuddy(seq_set, in_args.in_format)
        seqbuddy += seq_set.records

    seqbuddy = SeqBuddy(seqbuddy)

    seqbuddy.out_format = in_args.out_format if in_args.out_format else seq_set.out_format

    # ############################################# INTERNAL FUNCTION ################################################ #
    def _print_recs(_seqbuddy):
        if len(_seqbuddy.records) == 0:
            sys.stderr.write("Nothing returned.\n")
            return False

        if _seqbuddy.out_format == "phylipi":
            _output = phylipi(_seqbuddy)

        elif _seqbuddy.out_format == "phylipis":
            _output = phylipi(_seqbuddy, "strict")

        else:
            tmp_dir = TemporaryDirectory()
            with open("%s/seqs.tmp" % tmp_dir.name, "w") as ofile:
                SeqIO.write(_seqbuddy.records, ofile, _seqbuddy.out_format)

            with open("%s/seqs.tmp" % tmp_dir.name, "r") as ifile:
                _output = ifile.read()

        if in_args.in_place and in_place_allowed:
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

    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Purge
    if in_args.purge:
        purged_seqs, deleted = purge(seqbuddy, in_args.purge)

        if not in_args.quiet:
            stderr_output = "### Deleted record mapping ###\n"
            for seq_id in deleted:
                stderr_output += "%s\n" % seq_id
                for del_seq_id in deleted[seq_id]:
                    stderr_output += "%s, " % del_seq_id
                stderr_output = stderr_output.strip(", ") + "\n\n"

            sys.stderr.write(stderr_output.strip())
            sys.stderr.write("\n##############################\n\n")

        _print_recs(purged_seqs)

    # BL2SEQ
    if in_args.bl2seq:
        output = bl2seq(seqbuddy)
        sys.stdout.write("#query\tsubject\t%_ident\tlength\tevalue\tbit_score\n")
        ids_already_seen = []
        for query_id in output:
            ids_already_seen.append(query_id)
            for subj_id in output[query_id]:
                if subj_id in ids_already_seen:
                    continue

                ident, length, evalue, bit_score = output[query_id][subj_id]
                sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (query_id, subj_id, ident, length, evalue, bit_score))

    # BLAST
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
        reverse = True if in_args.params and in_args.params[0] == "rev" else False
        _print_recs(order_ids(seqbuddy, _reverse=reverse))

    # Find repeat sequences or ids
    if in_args.find_repeats:
        if in_args.params:
            columns = int(in_args.params[0])
        else:
            columns = 1

        unique, rep_ids, rep_seqs = find_repeats(seqbuddy)
        output = ""
        if len(rep_ids) > 0:
            output += "#### Records with duplicate IDs: ####\n"
            counter = 1
            for next_id in rep_ids:
                output += "%s\t" % next_id
                if counter % columns == 0:
                    output = "%s\n" % output.strip()
                counter += 1

            output = "%s\n\n" % output.strip()

        else:
            output += "#### No records with duplicate IDs ####\n\n"

        if len(rep_seqs) > 0:
            output += "#### Records with duplicate sequences: ####\n"
            counter = 1
            for next_id in rep_seqs:
                output += "["
                for seq_id in rep_seqs[next_id]:
                    output += "%s, " % seq_id
                output = "%s], " % output.strip(", ")

                if counter % columns == 0:
                    output = "%s\n" % output.strip(", ")

                counter += 1

            output = "%s\n\n" % output.strip(", ")
        else:
            output += "#### No records with duplicate sequences ####\n\n"

        sys.stdout.write("%s\n" % output)

    # Delete repeats
    if in_args.delete_repeats:
        if in_args.params:
            columns = int(in_args.params[0])
        else:
            columns = 1

        unique, rep_ids, rep_seqs = find_repeats(seqbuddy)
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
            unique, rep_ids, rep_seqs = find_repeats(seqbuddy)

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
            sys.stderr.write("# ################################################################ #\n")
            sys.stderr.write("%s\n" % stderr_output.strip())
            sys.stderr.write("# ################################################################ #\n\n")

            _print_recs(delete_repeats(seqbuddy, 'seqs'))

        else:
            sys.stderr.write("No duplicate records found\n")

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
                sys.stderr.write(output)

        if len(deleted_seqs) == 0:
            sys.stderr.write("# ################################################################ #\n")
            sys.stderr.write("# No sequence identifiers match %s\n" % ", ".join(in_args.delete_records))
            sys.stderr.write("# ################################################################ #\n")

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
        new_list = SeqBuddy([])
        for infile in in_args.sequence:
            new_list.records += SeqBuddy(infile).records

        new_list.out_format = in_args.out_format if in_args.out_format else seqbuddy.out_format
        _print_recs(new_list)

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
        seqbuddy = rename(seqbuddy, in_args.rename_ids[0], in_args.rename_ids[1], num)
        _print_recs(seqbuddy)

    # Uppercase
    if in_args.uppercase:
        _print_recs(uppercase(seqbuddy))

    # Lowercase
    if in_args.lowercase:
        _print_recs(lowercase(seqbuddy))

    # Transcribe
    if in_args.transcribe:
        if seqbuddy.alpha != IUPAC.ambiguous_dna:
            raise AttributeError("You need to provide a DNA sequence.")
        seqbuddy = dna2rna(seqbuddy)
        _print_recs(seqbuddy)

    # Back Transcribe
    if in_args.back_transcribe:
        if seqbuddy.alpha != IUPAC.ambiguous_rna:
            raise AttributeError("You need to provide an RNA sequence.")
        seqbuddy = rna2dna(seqbuddy)
        _print_recs(seqbuddy)

    # Complement
    if in_args.complement:
        _print_recs(complement(seqbuddy))

    # Reverse complement
    if in_args.reverse_complement:
        _print_recs(reverse_complement(seqbuddy))

    # List identifiers
    if in_args.list_ids:
        if in_args.params:
            columns = int(in_args.params[0])
        else:
            columns = 1
        output = ""
        counter = 1
        for rec in seqbuddy.records:
            output += "%s\t" % rec.id
            if counter % columns == 0:
                output = "%s\n" % output.strip()
            counter += 1
        sys.stdout.write("%s\n" % output.strip())

    # Translate CDS
    if in_args.translate:
        if seqbuddy.alpha == IUPAC.protein:
            raise AttributeError("You need to supply DNA or RNA sequences to translate")

        if in_args.quiet:
            _print_recs(translate_cds(seqbuddy, quiet=True))
        else:
            _print_recs(translate_cds(seqbuddy))

    # Shift reading frame
    if in_args.select_frame:
        _print_recs(select_frame(seqbuddy, in_args.select_frame))

    # Translate 6 reading frames
    if in_args.translate6frames:
        if seqbuddy.alpha == IUPAC.protein:
            raise AttributeError("You need to supply DNA or RNA sequences to translate")

        seqbuddy = translate6frames(seqbuddy)
        if in_args.out_format:
            seqbuddy.out_format = in_args.out_format

        _print_recs(seqbuddy)

    # Back translate CDS
    if in_args.back_translate:
        if in_args.params:
            in_args.params = [i.upper() for i in in_args.params]
            mode = [i for i in in_args.params if i in ['RANDOM', 'R', "OPTIMIZED", "O"]]
            mode = "RANDOM" if len(mode) == 0 else mode[0]
            species = [i for i in in_args.params if i in ['HUMAN', 'H', "MOUSE", "M",
                                                          "YEAST", "Y", "ECOLI", "E"]]
            species = None if len(species) == 0 else species[0]
        else:
            mode = "RANDOM"
            species = None

        _print_recs(back_translate(seqbuddy, mode, species))

    # Concatenate sequences
    if in_args.concat_seqs:
        clean = False if not in_args.params or in_args.params[0] != "clean" else True
        seqbuddy = concat_seqs(seqbuddy, clean)
        if in_args.out_format:
            seqbuddy.out_format = in_args.out_format
        _print_recs(seqbuddy)

    # Count number of sequences in a file
    if in_args.num_seqs:
        sys.stdout.write("%s\n" % len(seqbuddy.records))

    # Average length of sequences
    if in_args.ave_seq_length:
        clean = False if not in_args.params or in_args.params[0] != "clean" else True
        sys.stdout.write("%s\n" % round(ave_seq_length(seqbuddy, clean), 2))

    # Pull sequence ends
    if in_args.pull_record_ends:
        amount, which_end = in_args.pull_record_ends
        amount = int(amount)
        new_seqs = pull_record_ends(seqbuddy, amount, which_end)
        _print_recs(new_seqs)

    # Extract regions
    if in_args.extract_region:
        start, end = in_args.extract_region
        seqbuddy = extract_range(seqbuddy, start, end)
        _print_recs(seqbuddy)

    # Pull records
    if in_args.pull_records:
        search = in_args.pull_records
        records = pull_recs(seqbuddy, search)
        _print_recs(records)

    # Pull random records
    if in_args.pull_random_record:
        count = 1 if not in_args.params else in_args.params[0]
        try:
            count = int(count)
        except ValueError:
            sys.exit("Error: When passing in the -p flag to --pull_random_recs, the value must be convertable into an "
                     "integer. You passed in '%s'..." % count)
        _print_recs(pull_random_recs(seqbuddy, count))

    # Hash sequence ids
    if in_args.hash_seq_ids:
        hash_length = in_args.params[0] if in_args.params else 10
        hashed = hash_seqeunce_ids(seqbuddy, hash_length)
        hash_table = "# Hash table\n"
        for seq in hashed[0]:
            hash_table += "%s,%s\n" % (seq[0], seq[1])
        _stderr("%s\n" % hash_table, in_args.quiet)
        _print_recs(hashed[1])

    # Delete metadata
    if in_args.delete_metadata:
        _print_recs(delete_metadata(seqbuddy))

    # Raw Seq
    if in_args.raw_seq:
        seqbuddy = clean_seq(seqbuddy)
        output = ""
        for rec in seqbuddy.records:
            output += "%s\n\n" % rec.seq
        sys.stdout.write("%s\n" % output.strip())

    # Clean Seq
    if in_args.clean_seq:
        _print_recs(clean_seq(seqbuddy))

    # Guess format
    if in_args.guess_format:
        for seq_set in in_args.sequence:
            sys.stdout.write("%s\t-->\t%s\n" % (seq_set, SeqBuddy(seq_set).in_format))

    # Guess alphabet
    if in_args.guess_alphabet:
        for seq_set in in_args.sequence:
            seqbuddy = SeqBuddy(seq_set)
            sys.stdout.write("%s\t-->\t" % seq_set)
            if seqbuddy.alpha == IUPAC.protein:
                sys.stdout.write("prot\n")
            elif seqbuddy.alpha == IUPAC.ambiguous_dna:
                sys.stdout.write("dna\n")
            elif seqbuddy.alpha == IUPAC.ambiguous_rna:
                sys.stdout.write("rna\n")
            else:
                sys.stdout.write("Undetermined\n")

    # Map features from cDNA over to protein
    if in_args.map_features_dna2prot:
        file1, file2 = in_args.sequence[:2]

        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)

        if file1.alpha == file2.alpha:
            raise AttributeError("You must provide one DNA file and one protein file")

        if file1.alpha == IUPAC.protein:
            prot = file1
            dna = file2
        else:
            in_args.sequence[0] = in_args.sequence[1]  # in case the -i flag is thrown
            prot = file2
            dna = file1

        new_seqs = map_features_dna2prot(dna, prot)
        _print_recs(new_seqs)

    # Map features from protein over to cDNA
    if in_args.map_features_prot2dna:
        file1, file2 = in_args.sequence[:2]

        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)

        if file1.alpha == file2.alpha:
            raise AttributeError("You must provide one DNA file and one protein file")

        if file1.alpha != IUPAC.protein:
            dna = file1
            prot = file2
        else:
            in_args.sequence[0] = in_args.sequence[1]  # in case the -i flag is thrown
            dna = file2
            prot = file1

        new_seqs = map_features_prot2dna(prot, dna)
        _print_recs(new_seqs)

    # Combine feature sets from two files into one
    if in_args.combine_features:
        file1, file2 = in_args.sequence[:2]
        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)
        new_seqs = combine_features(file1, file2)
        _print_recs(new_seqs)

    # Order sequence features by their position in the sequence
    if in_args.order_features_by_position:
        reverse = True if in_args.params and in_args.params[0] == "rev" else False
        _print_recs(order_features_by_position(seqbuddy, reverse))

    # Order sequence features alphabetically
    if in_args.order_features_alphabetically:
        reverse = True if in_args.params and in_args.params[0] == "rev" else False
        _print_recs(order_features_alphabetically(seqbuddy, reverse))