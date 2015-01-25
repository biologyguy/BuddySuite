#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Nov 20 2014 
# 38 tools and counting

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
import pdb
import sys
import os
import re
import string
from copy import copy
from random import sample, choice, randint
from math import ceil
from tempfile import TemporaryDirectory
from subprocess import Popen, PIPE
from shutil import which

# Third party package imports
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError


# ##################################################### WISH LIST #################################################### #


# ################################################# HELPER FUNCTIONS ################################################# #


class SeqBuddy():  # Open a file or read a handle and parse, or convert raw into a Seq object
    def __init__(self, _input, _in_format=None, _out_format=None):
        if not _in_format:
            self.in_format = guess_format(_input)
            self.out_format = str(self.in_format) if not _out_format else _out_format
        if not self.in_format:
            sys.exit("Error: could not determine the seq format in SeqBuddy(). "
                     "Try explicitly setting with -f flag.")

        self.out_format = self.in_format if not _out_format else _out_format

        if str(type(_input)) == "<class '__main__.SeqBuddy'>":
            _sequences = _input.seqs

        elif isinstance(_input, list):
            # make sure that the list is actually SeqIO records (just test a few...)
            for _seq in _input[:3]:
                if type(_seq) != SeqRecord:
                    sys.exit("Error: Seqlist is not populated with SeqRecords.")
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

        self.seqs = _sequences


def guess_alphabet(_seqs):  # Does not handle ambiguous dna
    _seqs = _seqs if isinstance(_seqs, list) else _seqs.seqs
    _sequence = ""
    for next_seq in _seqs:
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
    if str(type(_input)) == "<class '__main__.SeqBuddy'>":
        return _input.in_format

    # If input is a handle or path, try to read the file in each format, and assume success if not error and # seqs > 0
    if os.path.isfile(str(_input)):
        _input = open(_input, "r")

    if str(type(_input)) == "<class '_io.TextIOWrapper'>":
        possible_formats = ["phylip-relaxed", "stockholm", "fasta", "gb", "nexus"]
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
        sys.exit("Error: Unsupported _input argument in guess_format(). %s" % _input)


def phylipi(_input, _format="relaxed"):  # _format in ["strict", "relaxed"]
    max_id_length = 0
    max_seq_length = 0
    for _seq in _input.seqs:
        max_id_length = len(_seq.id) if len(_seq.id) > max_id_length else max_id_length
        max_seq_length = len(_seq.seq) if len(_seq.seq) > max_seq_length else max_seq_length

    _output = " %s %s\n" % (len(_input.seqs), max_seq_length)
    for _seq in _input.seqs:
        _seq_id = _seq.id.ljust(max_id_length) if _format == "relaxed" else _seq.id[:10].ljust(10)
        _output += "%s  %s\n" % (_seq_id, _seq.seq)

    return _output
# #################################################################################################################### #


def blast(_seqs, blast_db, blast_path=None, blastdbcmd=None):  # ToDo: Allow weird binary names to work
    if not blast_path:
        blast_path = which("blastp") if _seqs.alpha == IUPAC.protein else which("blastn")

    blast_check = Popen("%s -version" % blast_path, stdout=PIPE, shell=True).communicate()
    blast_check = re.search("([a-z])*[^:]", blast_check[0].decode("utf-8"))
    if blast_check:
        blast_check = blast_check.group(0)

    # ToDo Check NCBI++ tools are a conducive version (2.2.29 and above, I think [maybe .28])

    # Check to make sure blast is in path and ensure that the blast_db is present
    blast_db = os.path.abspath(blast_db)
    if blast_check == "blastp":
        if not which(blast_path):
            raise FileNotFoundError("blastp")

        if not os.path.isfile("%s.pin" % blast_db) or not os.path.isfile("%s.phr" % blast_db) \
                or not os.path.isfile("%s.psq" % blast_db):
            sys.exit("Error:\tBlastp database not found at '%s'" % blast_db)
    elif blast_check == "blastn":
        if not which(blast_path):
            raise FileNotFoundError("blastn")

        if not os.path.isfile("%s.nin" % blast_db) or not os.path.isfile("%s.nhr" % blast_db) \
                or not os.path.isfile("%s.nsq" % blast_db):
            sys.exit("Error:\tBlastn database not found at '%s'" % blast_db)
    else:
        sys.exit("Blast binary doesn't seem to work, at %s" % blast_path)

    if not blastdbcmd:
        blastdbcmd = "blastdbcmd"

    if not which(blastdbcmd):
        raise FileNotFoundError("blastdbcmd")

    # Check that blastdb was made with the -parse_seqids flag
    extensions = ["pog", "psd", "psi"] if blast_check == "blastp" else ["nog", "nsd", "nsi"]
    if not os.path.isfile("%s.%s" % (blast_db, extensions[0])) or not \
            os.path.isfile("%s.%s" % (blast_db, extensions[1])) or not \
            os.path.isfile("%s.%s" % (blast_db, extensions[2])):
        sys.exit("Error: Incorrect blastdb. When making the blast database, please use the -parse_seqids flag.")

    tmp_dir = TemporaryDirectory()
    with open("%s/tmp.fa" % tmp_dir.name, "w") as ofile:
        SeqIO.write(_seqs.seqs, ofile, "fasta")

    Popen("%s -db %s -query %s/tmp.fa -out %s/out.txt -num_threads 20 -evalue 0.01 -outfmt 6" %
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


def shuffle(_seqs):
    _output = []
    for _ in range(len(_seqs.seqs)):
        random_index = randint(1, len(_seqs.seqs)) - 1
        _output.append(_seqs.seqs.pop(random_index))
    _seqs.seqs = _output
    return _seqs


def rna2dna(_seqs):
    _output = []
    for _seq in _seqs.seqs:
        _seq.seq = Seq(str(_seq.seq.back_transcribe()), alphabet=IUPAC.ambiguous_dna)
        _output.append(_seq)
    _seqs.seqs = _output
    return _seqs


def dna2rna(_seqs):
    _output = []
    for _seq in _seqs.seqs:
        _seq.seq = Seq(str(_seq.seq.transcribe()), alphabet=IUPAC.ambiguous_rna)
        _output.append(_seq)
    _seqs.seqs = _output
    return _seqs


def complement(_seqs):
    if _seqs.alpha == IUPAC.protein:
        sys.exit("Error: The complement function requires a nucleic acid sequence, not protein.")
    for _seq in _seqs.seqs:
        _seq.seq = _seq.seq.complement()
    return _seqs


def reverse_complement(_seqs):
    if _seqs.alpha == IUPAC.protein:
        sys.exit("Error: The complement function requires a nucleic acid sequence, not protein.")
    for _seq in _seqs.seqs:
        _seq.seq = _seq.seq.reverse_complement()
    return _seqs


def translate_cds(_seqs):
    _output = []
    for _seq in _seqs.seqs:
        try:
            _seq.seq = _seq.seq.translate(cds=True, to_stop=True)
        except TranslationError as e1:
            _seq.seq = Seq(str(_seq.seq)[:(len(str(_seq.seq)) - len(str(_seq.seq)) % 3)])
            try:
                _seq.seq = _seq.seq.translate()
                sys.stderr.write("Warning: %s is not a standard CDS\t-->\t%s\n" % (_seq.id, e1))
            except TranslationError as e2:
                sys.stderr.write("Error: %s failed to translate\t-->\t%s\n" % (_seq.id, e2))

        _seq.seq.alphabet = IUPAC.protein
        _output.append(_seq)
    _seqs.seqs = _output
    return _seqs


def translate6frames(_seqs):
    _output = []
    for _seq in _seqs.seqs:
        for i in range(3):
            x = len(str(_seq.seq)) if len(str(_seq.seq)) % 3 * (-1) == 0 else len(str(_seq.seq)) % 3 * (-1)
            temp_seq = str(_seq.seq)[:x]
            temp_seq = Seq(temp_seq, alphabet=_seq.seq.alphabet)
            temp_seq = Seq(str(temp_seq.translate()), alphabet=IUPAC.protein)
            _output.append(SeqRecord(temp_seq, description="", id="%s_%s" % (_seq.id, i + 1),
                                     name="%s" % _seq.name))
            _seq.seq = Seq(str(_seq.seq)[1:], alphabet=_seq.seq.alphabet)

        revcomp = _seq.seq.reverse_complement()
        for i in range(3):
            x = len(str(revcomp)) if len(str(revcomp)) % 3 * (-1) == 0 else len(str(revcomp)) % 3 * (-1)
            temp_seq = str(revcomp)[:x]
            temp_seq = Seq(temp_seq, alphabet=revcomp.alphabet)
            temp_seq = Seq(str(temp_seq.translate()), alphabet=IUPAC.protein)
            _output.append(SeqRecord(temp_seq, description="", id="%s_revcomp_%s" % (_seq.id, i + 1),
                                     name="%s" % _seq.name))
            revcomp = Seq(str(revcomp)[1:], alphabet=revcomp.alphabet)

    _seqs = SeqBuddy(_output, _out_format="fasta")
    return _seqs


def back_translate(_seqs, _mode='random'):  # available modes --> random ToDo: Implement other modes
    codon_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
                   'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
                   'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
                   'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
                   'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                   'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
                   'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
                   'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
                   'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                   'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
                   'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
                   'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
                   'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*'}

    # codon preferences
    # H. sapiense
    h_sapiense = {'H': ([0.58, 0.42], ['CAC', 'CAT']), 'T': ([0.25, 0.11, 0.28, 0.36], ['ACT', 'ACG', 'ACA', 'ACC']),
                  'A': ([0.23, 0.27, 0.40, 0.11], ['GCA', 'GCT', 'GCC', 'GCG']), 'I': ([], ['ATC', 'ATT', 'ATA']),
                  'V': ([], ['GTG', 'GTT', 'GTA', 'GTC']), 'L': ([], ['CTC', 'CTG', 'TTG', 'CTA', 'CTT', 'TTA']),
                  'K': ([], ['AAA', 'AAG']), 'F': ([], ['TTT', 'TTC']), 'D': ([], ['GAC', 'GAT']), 'Y': ([], ['TAC', 'TAT']),
                  'Q': ([], ['CAA', 'CAG']), 'S': ([], ['TCT', 'TCC', 'TCG', 'TCA', 'AGT', 'AGC']),
                  'C': ([], ['TGC', 'TGT']), 'P': ([], ['CCC', 'CCA', 'CCT', 'CCG']), '*': ([], ['TAG', 'TGA', 'TAA']),
                  'R': ([], ['CGT', 'AGA', 'AGG', 'CGC', 'CGG', 'CGA']), 'E': ([], ['GAG', 'GAA']),
                  'N': ([], ['AAT', 'AAC']), 'G': ([], ['GGC', 'GGT', 'GGA', 'GGG']), 'W': ([], ['TGG']), 'M': ([], ['ATG'])}

    UUU F 0.46 17.6 (714298)  UCU S 0.19 15.2 (618711)  UAU Y 0.44 12.2 (495699)  UGU C 0.46 10.6 (430311)
    UUC F 0.54 20.3 (824692)  UCC S 0.22 17.7 (718892)  UAC Y 0.56 15.3 (622407)  UGC C 0.54 12.6 (513028)
    UUA L 0.08  7.7 (311881)  UCA S 0.15 12.2 (496448)  UAA * 0.30  1.0 ( 40285)  UGA * 0.47  1.6 ( 63237)
    UUG L 0.13 12.9 (525688)  UCG S 0.05  4.4 (179419)  UAG * 0.24  0.8 ( 32109)  UGG W 1.00 13.2 (535595)

    CUU L 0.13 13.2 (536515)  CCU P 0.29 17.5 (713233)  CGU R 0.08  4.5 (184609)
    CUC L 0.20 19.6 (796638)  CCC P 0.32 19.8 (804620)  CGC R 0.18 10.4 (423516)
    CUA L 0.07  7.2 (290751)  CCA P 0.28 16.9 (688038)  CAA Q 0.27 12.3 (501911)  CGA R 0.11  6.2 (250760)
    CUG L 0.40 39.6 (1611801)  CCG P 0.11  6.9 (281570)  CAG Q 0.73 34.2 (1391973)  CGG R 0.20 11.4 (464485)

    AUU I 0.36 16.0 (650473)  AAU N 0.47 17.0 (689701)  AGU S 0.15 12.1 (493429)
    AUC I 0.47 20.8 (846466)  AAC N 0.53 19.1 (776603)  AGC S 0.24 19.5 (791383)
    AUA I 0.17  7.5 (304565)  AAA K 0.43 24.4 (993621)  AGA R 0.21 12.2 (494682)
    AUG M 1.00 22.0 (896005)  AAG K 0.57 31.9 (1295568)  AGG R 0.21 12.0 (486463)

    GUU V 0.18 11.0 (448607)  GAU D 0.46 21.8 (885429)  GGU G 0.16 10.8 (437126)
    GUC V 0.24 14.5 (588138)  GAC D 0.54 25.1 (1020595)  GGC G 0.34 22.2 (903565)
    GUA V 0.12  7.1 (287712)  GAA E 0.42 29.0 (1177632)  GGA G 0.25 16.5 (669873)
    GUG V 0.46 28.1 (1143534)  GAG E 0.58 39.6 (1609975)  GGG G 0.25 16.5 (669768)

    # E. coli
    e_coli = {'H': ([], ['CAC', 'CAT']), 'T': ([], ['ACT', 'ACG', 'ACA', 'ACC']),
              'A': ([], ['GCA', 'GCT', 'GCC', 'GCG']), 'I': ([], ['ATC', 'ATT', 'ATA']),
              'V': ([], ['GTG', 'GTT', 'GTA', 'GTC']), 'L': ([], ['CTC', 'CTG', 'TTG', 'CTA', 'CTT', 'TTA']),
              'K': ([], ['AAA', 'AAG']), 'F': ([], ['TTT', 'TTC']), 'D': ([], ['GAC', 'GAT']), 'Y': ([], ['TAC', 'TAT']),
              'Q': ([], ['CAA', 'CAG']), 'S': ([], ['TCT', 'TCC', 'TCG', 'TCA', 'AGT', 'AGC']),
              'C': ([], ['TGC', 'TGT']), 'P': ([], ['CCC', 'CCA', 'CCT', 'CCG']), '*': ([], ['TAG', 'TGA', 'TAA']),
              'R': ([], ['CGT', 'AGA', 'AGG', 'CGC', 'CGG', 'CGA']), 'E': ([], ['GAG', 'GAA']),
              'N': ([], ['AAT', 'AAC']), 'G': ([], ['GGC', 'GGT', 'GGA', 'GGG']), 'W': ([], ['TGG']), 'M': ([], ['ATG'])}

    UUU F 0.57 25.5 (     4)  UCU S 0.29 25.5 (     4)  UAU Y 0.78 44.6 (     7)  UGU C 1.00 19.1 (     3)
    UUC F 0.43 19.1 (     3)  UCC S 0.07  6.4 (     1)  UAC Y 0.22 12.7 (     2)  UGC C 0.00  0.0 (     0)
    UUA L 0.42 51.0 (     8)  UCA S 0.29 25.5 (     4)  UAA * 1.00  6.4 (     1)  UGA * 0.00  0.0 (     0)
    UUG L 0.21 25.5 (     4)  UCG S 0.21 19.1 (     3)  UAG * 0.00  0.0 (     0)  UGG W 1.00  6.4 (     1)

    CUU L 0.16 19.1 (     3)  CCU P 0.40 12.7 (     2)  CAU H 1.00 19.1 (     3)  CGU R 0.20  6.4 (     1)
    CUC L 0.05  6.4 (     1)  CCC P 0.00  0.0 (     0)  CAC H 0.00  0.0 (     0)  CGC R 0.00  0.0 (     0)
    CUA L 0.00  0.0 (     0)  CCA P 0.20  6.4 (     1)  CAA Q 0.50  6.4 (     1)  CGA R 0.00  0.0 (     0)
    CUG L 0.16 19.1 (     3)  CCG P 0.40 12.7 (     2)  CAG Q 0.50  6.4 (     1)  CGG R 0.00  0.0 (     0)

    AUU I 0.73 51.0 (     8)  ACU T 0.00  0.0 (     0)  AAU N 0.69 57.3 (     9)  AGU S 0.00  0.0 (     0)
    AUC I 0.00  0.0 (     0)  ACC T 0.33  6.4 (     1)  AAC N 0.31 25.5 (     4)  AGC S 0.14 12.7 (     2)
    AUA I 0.27 19.1 (     3)  ACA T 0.67 12.7 (     2)  AAA K 0.83 31.8 (     5)  AGA R 0.60 19.1 (     3)
    AUG M 1.00 31.8 (     5)  ACG T 0.00  0.0 (     0)  AAG K 0.17  6.4 (     1)  AGG R 0.20  6.4 (     1)

    GUU V 0.33 12.7 (     2)  GCU A 0.30 19.1 (     3)  GAU D 0.88 44.6 (     7)  GGU G 0.22 12.7 (     2)
    GUC V 0.50 19.1 (     3)  GCC A 0.20 12.7 (     2)  GAC D 0.12  6.4 (     1)  GGC G 0.11  6.4 (     1)
    GUA V 0.17  6.4 (     1)  GCA A 0.30 19.1 (     3)  GAA E 0.71 76.4 (    12)  GGA G 0.56 31.8 (     5)
    GUG V 0.00  0.0 (     0)  GCG A 0.20 12.7 (     2)  GAG E 0.29 31.8 (     5)  GGG G 0.11  6.4 (     1)

    # S. cerevisiae
    s_cerevisiae = {'H': ([], ['CAC', 'CAT']), 'T': ([], ['ACT', 'ACG', 'ACA', 'ACC']),
                    'A': ([], ['GCA', 'GCT', 'GCC', 'GCG']), 'I': ([], ['ATC', 'ATT', 'ATA']),
                    'V': ([], ['GTG', 'GTT', 'GTA', 'GTC']), 'L': ([], ['CTC', 'CTG', 'TTG', 'CTA', 'CTT', 'TTA']),
                    'K': ([], ['AAA', 'AAG']), 'F': ([], ['TTT', 'TTC']), 'D': ([], ['GAC', 'GAT']), 'Y': ([], ['TAC', 'TAT']),
                    'Q': ([], ['CAA', 'CAG']), 'S': ([], ['TCT', 'TCC', 'TCG', 'TCA', 'AGT', 'AGC']),
                    'C': ([], ['TGC', 'TGT']), 'P': ([], ['CCC', 'CCA', 'CCT', 'CCG']), '*': ([], ['TAG', 'TGA', 'TAA']),
                    'R': ([], ['CGT', 'AGA', 'AGG', 'CGC', 'CGG', 'CGA']), 'E': ([], ['GAG', 'GAA']),
                    'N': ([], ['AAT', 'AAC']), 'G': ([], ['GGC', 'GGT', 'GGA', 'GGG']), 'W': ([], ['TGG']), 'M': ([], ['ATG'])}
    UUU F 0.59 26.1 (170666)  UCU S 0.26 23.5 (153557)  UAU Y 0.56 18.8 (122728)  UGU C 0.63  8.1 ( 52903)
    UUC F 0.41 18.4 (120510)  UCC S 0.16 14.2 ( 92923)  UAC Y 0.44 14.8 ( 96596)  UGC C 0.37  4.8 ( 31095)
    UUA L 0.28 26.2 (170884)  UCA S 0.21 18.7 (122028)  UAA * 0.47  1.1 (  6913)  UGA * 0.30  0.7 (  4447)
    UUG L 0.29 27.2 (177573)  UCG S 0.10  8.6 ( 55951)  UAG * 0.23  0.5 (  3312)  UGG W 1.00 10.4 ( 67789)

    CUU L 0.13 12.3 ( 80076)  CCU P 0.31 13.5 ( 88263)  CAU H 0.64 13.6 ( 89007)  CGU R 0.14  6.4 ( 41791)
    CUC L 0.06  5.4 ( 35545)  CCC P 0.15  6.8 ( 44309)  CAC H 0.36  7.8 ( 50785)  CGC R 0.06  2.6 ( 16993)
    CUA L 0.14 13.4 ( 87619)  CCA P 0.42 18.3 (119641)  CAA Q 0.69 27.3 (178251)  CGA R 0.07  3.0 ( 19562)
    CUG L 0.11 10.5 ( 68494)  CCG P 0.12  5.3 ( 34597)  CAG Q 0.31 12.1 ( 79121)  CGG R 0.04  1.7 ( 11351)

    AUU I 0.46 30.1 (196893)  ACU T 0.35 20.3 (132522)  AAU N 0.59 35.7 (233124)  AGU S 0.16 14.2 ( 92466)
    AUC I 0.26 17.2 (112176)  ACC T 0.22 12.7 ( 83207)  AAC N 0.41 24.8 (162199)  AGC S 0.11  9.8 ( 63726)
    AUA I 0.27 17.8 (116254)  ACA T 0.30 17.8 (116084)  AAA K 0.58 41.9 (273618)  AGA R 0.48 21.3 (139081)
    AUG M 1.00 20.9 (136805)  ACG T 0.14  8.0 ( 52045)  AAG K 0.42 30.8 (201361)  AGG R 0.21  9.2 ( 60289)

    GUU V 0.39 22.1 (144243)  GCU A 0.38 21.2 (138358)  GAU D 0.65 37.6 (245641)  GGU G 0.47 23.9 (156109)
    GUC V 0.21 11.8 ( 76947)  GCC A 0.22 12.6 ( 82357)  GAC D 0.35 20.2 (132048)  GGC G 0.19  9.8 ( 63903)
    GUA V 0.21 11.8 ( 76927)  GCA A 0.29 16.2 (105910)  GAA E 0.70 45.6 (297944)  GGA G 0.22 10.9 ( 71216)
    GUG V 0.19 10.8 ( 70337)  GCG A 0.11  6.2 ( 40358)  GAG E 0.30 19.2 (125717)  GGG G 0.12  6.0 ( 39359)


    aa_lookup_table = {'H': ['CAC', 'CAT'], 'T': ['ACT', 'ACG', 'ACA', 'ACC'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'],
                       'I': ['ATC', 'ATT', 'ATA'], 'V': ['GTG', 'GTT', 'GTA', 'GTC'],
                       'L': ['CTC', 'CTG', 'TTG', 'CTA', 'CTT', 'TTA'], 'K': ['AAA', 'AAG'], 'F': ['TTT', 'TTC'],
                       'D': ['GAC', 'GAT'], 'Y': ['TAC', 'TAT'], 'Q': ['CAA', 'CAG'],
                       'S': ['TCT', 'TCC', 'TCG', 'TCA', 'AGT', 'AGC'], 'C': ['TGC', 'TGT'],
                       'P': ['CCC', 'CCA', 'CCT', 'CCG'], '*': ['TAG', 'TGA', 'TAA'],
                       'R': ['CGT', 'AGA', 'AGG', 'CGC', 'CGG', 'CGA'], 'E': ['GAG', 'GAA'],
                       'N': ['AAT', 'AAC'], 'G': ['GGC', 'GGT', 'GGA', 'GGG'], 'W': ['TGG'], 'M': ['ATG']}

    if _seqs.alpha != IUPAC.protein:
        sys.exit("Error: The input sequence needs to be protein, not %s" % _seqs.alpha)

    for _seq in _seqs.seqs:
        dna_seq = ""
        if _mode == 'random':
            for aa in _seq.seq:
                dna_seq += choice(aa_lookup_table[aa])
            _seq.seq = Seq(dna_seq, alphabet=IUPAC.ambiguous_dna)

        else:
            sys.exit("Error: Mode '%s' not implemented. Valid choices are random, blahh, blahh, or blahh" % _mode)

    return _seqs


def ave_seq_length(_seqs):
    sum_length = 0.
    for _seq in _seqs.seqs:
        sum_length += len(_seq.seq)
    return sum_length / len(_seqs.seqs)


def concat_seqs(_seqs):
    _new_seq = ""
    concat_ids = []
    features = []
    for _seq in _seqs.seqs:
        location = FeatureLocation(len(_new_seq), len(_new_seq) + len(str(_seq.seq)))
        feature = SeqFeature(location=location, id=_seq.id, type=_seq.id[:15])
        features.append(feature)
        concat_ids.append(_seq.id)
        _new_seq += str(_seq.seq)

    concat_ids = "|".join(concat_ids)
    _new_seq = [SeqRecord(Seq(_new_seq, alphabet=_seqs.alpha),
                          description=concat_ids, id="concatination", features=features)]
    _seqs = SeqBuddy(_new_seq)
    _seqs.out_format = "gb"
    return _seqs


def clean_seq(_seqs):
    """remove all non-sequence charcters from sequence strings"""
    _output = []
    for _seq in _seqs.seqs:
        _seq.seq = str(_seq.seq).upper()
        if _seqs.alpha == IUPAC.protein:
            _seq.seq = Seq(re.sub("[^ACDEFGHIKLMNPQRSTVWXY]", "", str(_seq.seq)), alphabet=_seqs.alpha)
        else:
            _seq.seq = Seq(re.sub("[^ATGCXNU]", "", str(_seq.seq)), alphabet=_seqs.alpha)

        _output.append(_seq)

    _seqs.seqs = _output
    return _seqs


def delete_metadata(_seqs):
    _new_seqs = []
    for _seq in _seqs.seqs:
        _new_seqs.append(SeqRecord(Seq(str(_seq.seq), alphabet=_seqs.alpha), id=_seq.id, name='', description=''))
    _seqs.seqs = _new_seqs
    return _seqs


# Apply DNA features to protein sequences
def map_features_dna2prot(dna_seqs, prot_seqs):
    prot_dict = SeqIO.to_dict(prot_seqs.seqs)
    dna_dict = SeqIO.to_dict(dna_seqs.seqs)
    _new_seqs = {}
    for _seq_id in dna_dict:
        if _seq_id not in prot_dict:
            sys.stderr.write("Warning: %s is in protein file, but not cDNA file\n" % _seq_id)
            continue

        _new_seqs[_seq_id] = prot_dict[_seq_id]

        for feature in dna_dict[_seq_id].features:
            _start = feature.location.start / 3
            _end = feature.location.end / 3
            location = FeatureLocation(ceil(_start), ceil(_end))
            feature.location = location
            prot_dict[_seq_id].features.append(feature)

    for _seq_id in prot_dict:
        if _seq_id not in dna_dict:
            sys.stderr.write("Warning: %s is in cDNA file, but not protein file\n" % _seq_id)

    _seqs_list = [_new_seqs[_seq_id] for _seq_id in _new_seqs]
    _seqs = SeqBuddy(_seqs_list)
    _seqs.out_format = "gb"
    return _seqs


# Apply DNA features to protein sequences
def map_features_prot2dna(prot_seqs, dna_seqs):
    prot_dict = SeqIO.to_dict(prot_seqs.seqs)
    dna_dict = SeqIO.to_dict(dna_seqs.seqs)
    _new_seqs = {}
    for _seq_id in prot_dict:
        if _seq_id not in dna_dict:
            sys.stderr.write("Warning: %s is in protein file, but not cDNA file\n" % _seq_id)
            continue

        _new_seqs[_seq_id] = dna_dict[_seq_id]

        for feature in prot_dict[_seq_id].features:
            _start = feature.location.start * 3
            _end = feature.location.end * 3
            location = FeatureLocation(_start, _end)
            feature.location = location
            dna_dict[_seq_id].features.append(feature)

    for _seq_id in dna_dict:
        if _seq_id not in prot_dict:
            sys.stderr.write("Warning: %s is in cDNA file, but not protein file\n" % _seq_id)

    _seqs_list = [_new_seqs[_seq_id] for _seq_id in _new_seqs]
    _seqs = SeqBuddy(_seqs_list)
    _seqs.out_format = "gb"
    return _seqs


# Merge feature lists
def combine_features(seqs1, seqs2):
    # make sure there are no repeat ids
    _unique, _rep_ids, _rep_seqs = find_repeats(seqs1)
    if len(_rep_ids) > 0:
        sys.exit("Error: There are repeat IDs in the first file provided\n%s" % _rep_ids)

    _unique, _rep_ids, _rep_seqs = find_repeats(seqs2)
    if len(_rep_ids) > 0:
        sys.exit("Error: There are repeat IDs in the second file provided\n%s" % _rep_ids)

    seq_dict1 = {}
    seq_dict2 = {}

    for _seq in seqs1.seqs:
        seq_dict1[_seq.id] = _seq

    for _seq in seqs2.seqs:
        seq_dict2[_seq.id] = _seq

    # make sure that we're comparing apples to apples across all sequences (i.e., same alphabet)
    reference_alphabet = sample(seq_dict1.items(), 1)[0][1].seq.alphabet
    for _seq_id in seq_dict1:
        if type(seq_dict1[_seq_id].seq.alphabet) != type(reference_alphabet):
            error_mes = "You have mixed multiple alphabets into your sequences. Make sure everything is the same.\n" \
                        "\t%s in first set\n\tOffending alphabet: %s\n\tReference alphabet: %s" \
                        % (_seq_id, seq_dict1[_seq_id].seq.alphabet, reference_alphabet)
            sys.exit(error_mes)

    for _seq_id in seq_dict2:
        if type(seq_dict2[_seq_id].seq.alphabet) != type(reference_alphabet):
            error_mes = "You have mixed multiple alphabets into your sequences. Make sure everything is the same.\n" \
                        "\t%s in first set\n\tOffending alphabet: %s\n\tReference alphabet: %s" \
                        % (_seq_id, seq_dict2[_seq_id].seq.alphabet, reference_alphabet)
            sys.exit(error_mes)

    _new_seqs = {}
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
            sys.stderr.write("Warning: %s is only in the first set of sequences\n" % _seq_id)

        _new_seqs[_seq_id] = seq_dict1[_seq_id]

    for _seq_id in seq_dict2:
        if _seq_id not in seq_dict1:
            sys.stderr.write("Warning: %s is only in the first set of sequences\n" % _seq_id)
            _new_seqs[_seq_id] = seq_dict2[_seq_id]

    _new_seqs = SeqBuddy([_new_seqs[_seq_id] for _seq_id in _new_seqs], _out_format=seqs1.in_format)
    return _new_seqs


def order_features_by_position(_seqs):
    for _seq in _seqs.seqs:
        new_feature_list = [(int(_feature.location.start), _feature) for _feature in _seq.features]
        new_feature_list = sorted(new_feature_list, key=lambda x: x[0])
        new_feature_list = [_feature[1] for _feature in new_feature_list]
        _seq.features = new_feature_list
    return _seqs


def order_features_alphabetically(_seqs):
    for _seq in _seqs.seqs:
        new_feature_list = [(_feature.type, _feature) for _feature in _seq.features]
        new_feature_list = sorted(new_feature_list, key=lambda x: x[0])
        new_feature_list = [_feature[1] for _feature in new_feature_list]
        _seq.features = new_feature_list
    return _seqs


def hash_seqeunce_ids(_seqs):
    hash_list = []
    seq_ids = []
    for i in range(len(_seqs.seqs)):
        new_hash = ""
        seq_ids.append(_seqs.seqs[i].id)
        while True:
            new_hash = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])
            if new_hash in hash_list:
                continue
            else:
                hash_list.append(new_hash)
                break
        _seqs.seqs[i].id = new_hash
        _seqs.seqs[i].name = new_hash

    hash_map = []
    for i in range(len(hash_list)):
        hash_map.append((hash_list[i], seq_ids[i]))

    return [hash_map, _seqs]


def pull_recs(_seqs, _search):
    _output = []
    for _seq in _seqs.seqs:
        if re.search(_search, _seq.description) or re.search(_search, _seq.id) or re.search(_search, _seq.name):
            _output.append(_seq)
    _seqs.seqs = _output
    return _seqs


def pull_seq_ends(_seqs, _amount, _which_end):
    seq_ends = []
    for _seq in _seqs.seqs:
        if _which_end == 'front':
            _seq.seq = _seq.seq[:_amount]

        elif _which_end == "rear":
            _seq.seq = _seq.seq[-1 * _amount:]

        else:
            sys.exit("Error: you much pick 'front' or 'rear' as the third argument in pull_seq_ends.")
        seq_ends.append(_seq)
    _seqs.seqs = seq_ends
    return _seqs


def extract_range(_seqs, _start, _end):
    _start = 1 if int(_start) < 1 else _start
    # Don't use the standard index-starts-at-0... _end must be left for the range to be inclusive
    _start, _end = int(_start) - 1, int(_end)
    if _end < _start:
        sys.exit("Error at extract range: The value given for end of range is smaller than for the start of range.")

    for _seq in _seqs.seqs:
        _seq.seq = Seq(str(_seq.seq)[_start:_end], alphabet=_seq.seq.alphabet)
        _seq.description += " Sub-sequence extraction, from residue %s to %s" % (_start, _end)
        _features = []
        for _feature in _seq.features:
            if _feature.location.end < _start:
                continue
            if _feature.location.start > _end:
                continue

            feat_start = _feature.location.start - _start
            if feat_start < 0:
                feat_start = 0

            feat_end = _feature.location.end - _start
            if feat_end > len(str(_seq.seq)):
                feat_end = len(str(_seq.seq))

            new_location = FeatureLocation(feat_start, feat_end)
            _feature.location = new_location
            _features.append(_feature)
        _seq.features = _features
    return _seqs


def find_repeats(_seqs):
    unique_seqs = {}
    repeat_ids = {}
    repeat_seqs = {}

    # First find replicate IDs
    for _seq in _seqs.seqs:
        if _seq.id in repeat_ids:
            repeat_ids[_seq.id].append(_seq)
        elif _seq.id in unique_seqs:
            repeat_ids[_seq.id] = [_seq]
            repeat_ids[_seq.id].append(unique_seqs[_seq.id])
            del(unique_seqs[_seq.id])
        else:
            unique_seqs[_seq.id] = _seq

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


def delete_records(_seqs, search_str):
    _output = []
    _deleted = pull_recs(copy(_seqs), search_str).seqs
    for _seq in _seqs.seqs:
        if _seq in _deleted:
            continue
        else:
            _output.append(_seq)
    _seqs.seqs = _output
    return _seqs


def delete_features(_seqs, _pattern):
    for _seq in _seqs.seqs:
        retained_features = []
        for _feature in _seq.features:
            if not re.search(_pattern, _feature.type):
                retained_features.append(_feature)
        _seq.features = retained_features
    return _seqs


def delete_repeats(_seqs, scope='all'):  # scope in ['all', 'ids', 'seqs']
    # First, remove duplicate IDs
    if scope in ['all', 'ids']:
        _unique, _rep_ids, _rep_seqs = find_repeats(_seqs)
        if len(_rep_ids) > 0:
            for _rep_id in _rep_ids:
                store_one_copy = pull_recs(copy(_seqs), _rep_id).seqs[0]
                delete_records(_seqs, _rep_id)
                _seqs.seqs.append(store_one_copy)

    # Then remove duplicate sequences
    if scope in ['all', 'seqs']:
        _unique, _rep_ids, _rep_seqs = find_repeats(_seqs)
        if len(_rep_seqs) > 0:
            _rep_seq_ids = []
            for _seq in _rep_seqs:
                _rep_seq_ids.append([])
                for _rep_seq_id in _rep_seqs[_seq]:
                    _rep_seq_ids[-1].append(_rep_seq_id)

            for _rep_seqs in _rep_seq_ids:
                for _rep_seq in _rep_seqs[1:]:
                    delete_records(_seqs, _rep_seq)

    return _seqs


def rename(_seqs, query, replace=""):  # TODO Allow a replacement pattern increment (like numbers)
    for _seq in _seqs.seqs:
        new_name = re.sub(query, replace, _seq.id)
        _seq.id = new_name
        _seq.name = new_name
    return _seqs


def purge(_seqs, threshold):  # ToDo: Implement a way to return a certain # of sequences (i.e. auto-determine threshold)
    purge_set = {}
    for _seq in _seqs.seqs:
        _unique = True
        for _seq_id in purge_set:
            purge_seq = purge_set[_seq_id]["seq"]
            blast_seqs = SeqBuddy([_seq, purge_seq])
            _blast_res = bl2seq(blast_seqs)
            bit_score = float(_blast_res.split("\t")[5])
            if bit_score >= threshold:
                _unique = False
                purge_set[_seq_id]["sim_seq"].append(_seq.id)

        if _unique:
            purge_set[_seq.id] = {"seq": _seq, "sim_seq": []}

    _output = {"seqs": [], "deleted": {}}
    for _seq_id in purge_set:
        _output["seqs"].append(purge_set[_seq_id]["seq"])
        _output["deleted"][_seq_id] = purge_set[_seq_id]["sim_seq"]

    _seqs.seqs = _output["seqs"]
    return [_seqs, _output["deleted"]]


def bl2seq(_seqs):  # Does an all-by-all analysis, and does not return sequences
    blast_bin = "blastp" if _seqs.alpha == IUPAC.protein else "blastn"
    if not which(blast_bin):
        sys.exit("Error: %s not present in $PATH.")  # ToDo: Implement -p flag

    tmp_dir = TemporaryDirectory()
    _seqs_copy = _seqs.seqs[1:]
    subject_file = "%s/subject.fa" % tmp_dir.name
    query_file = "%s/query.fa" % tmp_dir.name
    _output = ""
    for subject in _seqs.seqs:
        with open(subject_file, "w") as ifile:
            SeqIO.write(subject, ifile, "fasta")

        for query in _seqs_copy:
            with open(query_file, "w") as ifile:
                SeqIO.write(query, ifile, "fasta")

            _blast_res = Popen("%s -subject %s -query %s -outfmt 6" %
                               (blast_bin, subject_file, query_file), stdout=PIPE, shell=True).communicate()
            _blast_res = _blast_res[0].decode().split("\n")[0].split("\t")
            if len(_blast_res) == 1:
                _output += "%s\t%s\t0\t0\t0\t0\n" % (subject.id, query.id)
            else:
                # values are: query, subject, %_ident, length, evalue, bit_score
                _output += "%s\t%s\t%s\t%s\t%s\t%s\n" % (_blast_res[0], _blast_res[1], _blast_res[2],
                                                         _blast_res[3], _blast_res[10], _blast_res[11].strip())

        _seqs_copy = _seqs_copy[1:]
    return _output.strip()


def uppercase(_seqs):
    for _seq in _seqs.seqs:
        _seq.seq = Seq(str(_seq.seq).upper(), alphabet=_seq.seq.alphabet)
    return _seqs


def lowercase(_seqs):
    for _seq in _seqs.seqs:
        _seq.seq = Seq(str(_seq.seq).lower(), alphabet=_seq.seq.alphabet)
    return _seqs


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
    parser.add_argument('-tr6', '--translate6frames', action='store_true',
                        help="Translate nucleotide sequences into all six reading frames")
    parser.add_argument('-btr', '--back_translate', action='store_true',
                        help="Convert amino acid sequences into codons. Select mode with -p flag ['random', <others>]")
    parser.add_argument('-d2r', '--transcribe', action='store_true',
                        help="Convert DNA sequences to RNA")
    parser.add_argument('-r2d', '--back_transcribe', action='store_true',
                        help="Convert RNA sequences to DNA")
    parser.add_argument('-cmp', '--complement', action='store_true',
                        help="Return complement of nucleotide sequence")
    parser.add_argument('-rcmp', '--reverse_complement', action='store_true',
                        help="Return reverse complement of nucleotide sequence")
    parser.add_argument('-li', '--list_ids', action='store_true',
                        help="Output all the sequence identifiers in a file. Use -p to specify # columns to write")
    parser.add_argument('-ns', '--num_seqs', action='store_true',
                        help="Counts how many sequences are present in an input file")
    parser.add_argument('-asl', '--ave_seq_length', action='store_true',
                        help="Return the average length of all sequences")
    parser.add_argument('-cts', '--concat_seqs', action='store_true',
                        help="Concatenate a bunch of sequences into a single solid string.")
    parser.add_argument('-fd2p', '--map_features_dna2prot', action='store_true',
                        help="Take the features annotated onto nucleotide sequences and map to protein sequences. "
                             "Both a protein and cDNA file must be passed in.")
    parser.add_argument('-fp2d', '--map_features_prot2dna', action='store_true',
                        help="Arguments: one cDNA file and one protein file")
    parser.add_argument('-ri', '--rename_ids', action='store', metavar=('<pattern>', '<substitution>'), nargs=2,
                        help="Replace some pattern in ids with something else.")
    parser.add_argument('-cf', '--combine_features', action='store_true',
                        help="Takes the features in two files and combines them for each sequence")
    parser.add_argument('-ofp', '--order_features_by_position', action='store_true',
                        help="Change the output order of sequence features, based on sequence position")
    parser.add_argument('-ofa', '--order_features_alphabetically', action='store_true',
                        help="Change the output order of sequence features, based on sequence position")
    parser.add_argument('-sf', '--screw_formats', action='store', metavar="<out_format>",
                        help="Change the file format to something else.")
    parser.add_argument('-sh', '--shuffle', action='store_true',
                        help="Randomly reorder the position of records in the file.")
    parser.add_argument('-hsi', '--hash_seq_ids', action='store_true',
                        help="Rename all the identifiers in a sequence list to a 10 character hash.")
    parser.add_argument('-pr', '--pull_records', action='store', metavar="<regex pattern>",
                        help="Get all the records with ids containing a given string")
    parser.add_argument('-pre', '--pull_record_ends', action='store', nargs=2, metavar="<amount (int)> <front|rear>",
                        help="Get the ends (front or rear) of all sequences in a file.")
    parser.add_argument('-er', '--extract_region', action='store', nargs=2, metavar="<start (int)> <end (int)>",
                        type=int, help="Pull out sub-sequences.")
    parser.add_argument('-dr', '--delete_records', action='store', nargs="+", metavar="<regex pattern>",
                        help="Remove records from a file. The deleted IDs are sent to stderr.")
    parser.add_argument('-df', '--delete_features', action='store', nargs="+", metavar="<regex pattern>",
                        help="Remove specified features from all records.")
    parser.add_argument('-drp', '--delete_repeats', action='store_true',
                        help="Strip repeat records (ids and/or identical sequences")
    parser.add_argument('-fr', '--find_repeats', action='store_true',
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
    parser.add_argument('-o', '--out_format', help="If you want a specific format output", action='store')
    parser.add_argument('-f', '--in_format', help="If SeqBuddy can't guess the file format, just specify it directly.",
                        action='store')
    
    in_args = parser.parse_args()

    in_place_allowed = False

    seqs = []
    seq_set = ""
    for seq_set in in_args.sequence:
        seq_set = SeqBuddy(seq_set, in_args.in_format)
        seqs += seq_set.seqs

    seqs = SeqBuddy(seqs)

    seqs.out_format = in_args.out_format if in_args.out_format else seq_set.out_format

    # ############################################# INTERNAL FUNCTION ################################################ #
    def _print_recs(_seqs):
        if len(_seqs.seqs) == 0:
            sys.stderr.write("Nothing returned.\n")
            return False

        if _seqs.out_format == "phylipi":
            _output = phylipi(_seqs)

        elif _seqs.out_format == "phylipis":
            _output = phylipi(_seqs, "strict")

        else:
            tmp_dir = TemporaryDirectory()
            with open("%s/seqs.tmp" % tmp_dir.name, "w") as ofile:
                SeqIO.write(_seqs.seqs, ofile, _seqs.out_format)

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
        in_place_allowed = True
        purged_seqs, deleted = purge(seqs, in_args.purge)

        stderr_output = "# Deleted record mapping #\n"
        for seq_id in deleted:
            stderr_output += "%s\n%s\n\n" % (seq_id, deleted[seq_id])

        sys.stderr.write(stderr_output)
        _print_recs(purged_seqs)

    # BL2SEQ
    if in_args.bl2seq:
        sys.stderr.write("#query\tsubject\t%_ident\tlength\tevalue\tbit_score\n")
        sys.stdout.write("%s\n" % bl2seq(seqs))

    # BLAST
    if in_args.blast:
        blast_binaries = _get_blast_binaries()
        blast_binary_path = blast_binaries["blastp"] if seqs.alpha == IUPAC.protein else blast_binaries["blastn"]
        try:
            blast_res = blast(seqs, in_args.blast, blast_path=blast_binary_path,
                              blastdbcmd=blast_binaries["blastdbcmd"])

        except FileNotFoundError as e:
            sys.exit("%s binary not found, explicitly set with the -p flag.\n"
                     "To pass in the path to both blast(p/n) and blastdbcmd, separate them with a space." % e)
        _print_recs(blast_res)

    # Shuffle
    if in_args.shuffle:
        in_place_allowed = True
        _print_recs(shuffle(seqs))

    # Delete repeats
    if in_args.delete_repeats:
        in_place_allowed = True
        if in_args.params:
            columns = int(in_args.params[0])
        else:
            columns = 1

        unique, rep_ids, rep_seqs = find_repeats(seqs)
        stderr_output = ""
        if len(rep_ids) > 0:
            stderr_output += "# Records with duplicate ids deleted (first instance retained)\n"
            counter = 1
            for seq in rep_ids:
                stderr_output += "%s\t" % seq
                if counter % columns == 0:
                    stderr_output = "%s\n" % stderr_output.strip()
                counter += 1
            stderr_output += "\n\n"

        rep_seq_ids = []
        for seq in rep_seqs:
            rep_seq_ids.append([])
            for rep_seq_id in rep_seqs[seq]:
                rep_seq_ids[-1].append(rep_seq_id)

        if len(rep_seq_ids) > 0:
            stderr_output += "# Records with duplicate sequence deleted (first instance retained)\n"
            counter = 1
            for rep_seqs in rep_seq_ids:
                for rep_seq in rep_seqs[1:]:
                    stderr_output += "%s\t" % rep_seq
                    if counter % columns == 0:
                        stderr_output = "%s\n" % stderr_output.strip()
                    counter += 1
            stderr_output += "\n"

        if stderr_output != "":
            sys.stderr.write("# ################################################################ #\n")
            sys.stderr.write("%s\n" % stderr_output.strip())
            sys.stderr.write("# ################################################################ #\n")

            _print_recs(delete_repeats(seqs))

        else:
            sys.stderr.write("No duplicate records found\n")

    # Delete records
    if in_args.delete_records:
        in_place_allowed = True
        if in_args.params:
            columns = int(in_args.params[0])
        else:
            columns = 1

        new_list = SeqBuddy(list(seqs.seqs))
        deleted_seqs = []
        for next_pattern in in_args.delete_records:
            deleted_seqs += pull_recs(copy(new_list), next_pattern).seqs
            delete_records(new_list, next_pattern)

        if len(deleted_seqs) > 0:
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

        new_list.out_format = in_args.out_format if in_args.out_format else seqs.out_format
        _print_recs(new_list)

    # Delete features
    if in_args.delete_features:
        in_place_allowed = True
        for next_pattern in in_args.delete_features:
            delete_features(seqs, next_pattern)
        _print_recs(seqs)

    # Merge
    if in_args.merge:
        new_list = SeqBuddy([])
        for infile in in_args.sequence:
            new_list.seqs += SeqBuddy(infile).seqs

        new_list.out_format = in_args.out_format if in_args.out_format else seqs.out_format
        _print_recs(new_list)

    # Screw formats
    if in_args.screw_formats:
        in_place_allowed = True
        seqs.out_format = in_args.screw_formats
        if in_args.in_place:  # Need to change the file extension
            os.remove(in_args.sequence[0])
            in_args.sequence[0] = ".".join(os.path.abspath(in_args.sequence[0]).split(".")[:-1]) + "." + seqs.out_format
            open(in_args.sequence[0], "w").close()

        _print_recs(seqs)

    # Renaming
    if in_args.rename_ids:
        in_place_allowed = True
        seqs = rename(seqs, in_args.rename_ids[0], in_args.rename_ids[1])
        _print_recs(seqs)

    # Uppercase
    if in_args.uppercase:
        in_place_allowed = True
        _print_recs(uppercase(seqs))

    # Lowercase
    if in_args.lowercase:
        in_place_allowed = True
        _print_recs(lowercase(seqs))

    # Transcribe
    if in_args.transcribe:
        in_place_allowed = True
        if seqs.alpha != IUPAC.ambiguous_dna:
            sys.exit("Error: You need to provide a DNA sequence.")
        seqs = dna2rna(seqs)
        _print_recs(seqs)

    # Back Transcribe
    if in_args.back_transcribe:
        in_place_allowed = True
        if seqs.alpha != IUPAC.ambiguous_rna:
            sys.exit("Error: You need to provide an RNA sequence.")
        seqs = rna2dna(seqs)
        _print_recs(seqs)

    # Complement
    if in_args.complement:
        in_place_allowed = True
        _print_recs(complement(seqs))

    # Reverse complement
    if in_args.reverse_complement:
        in_place_allowed = True
        _print_recs(reverse_complement(seqs))

    # List identifiers
    if in_args.list_ids:
        if in_args.params:
            columns = int(in_args.params[0])
        else:
            columns = 1
        output = ""
        counter = 1
        for seq in seqs.seqs:
            output += "%s\t" % seq.id
            if counter % columns == 0:
                output = "%s\n" % output.strip()
            counter += 1
        sys.stdout.write("%s\n" % output.strip())

    # Translate CDS
    if in_args.translate:
        in_place_allowed = True
        if seqs.alpha == IUPAC.protein:
            sys.exit("Error: you need to supply DNA or RNA sequences to translate")
        _print_recs(translate_cds(seqs))

    # Translate 6 reading frames
    if in_args.translate6frames:
        in_place_allowed = True
        if seqs.alpha == IUPAC.protein:
            sys.exit("Error: you need to supply DNA or RNA sequences to translate")

        seqs = translate6frames(seqs)
        if in_args.out_format:
            seqs.out_format = in_args.out_format

        _print_recs(seqs)

    # Back translate CDS
    if in_args.back_translate:
        in_place_allowed = True
        mode = in_args.params if in_args.params else 'random'
        _print_recs(back_translate(seqs, mode))

    # Concatenate sequences
    if in_args.concat_seqs:
        seqs = concat_seqs(seqs)
        if in_args.out_format:
            seqs.out_format = in_args.out_format
        _print_recs(seqs)

    # Count number of sequences in a file
    if in_args.num_seqs:
        sys.stdout.write("%s\n" % len(seqs.seqs))

    # Average length of sequences
    if in_args.ave_seq_length:
        sys.stdout.write("%s\n" % round(ave_seq_length(seqs), 2))

    # Find repeat sequences or ids
    if in_args.find_repeats:
        unique, rep_ids, rep_seqs = find_repeats(seqs)
        output = ""
        if len(rep_ids) > 0:
            output += "Records with duplicate IDs:\n"
            for next_id in rep_ids:
                output += "%s, " % next_id
            output = "%s\n\n" % output.strip(", ")

        else:
            output += "No records with duplicate IDs\n\n"
        if len(rep_seqs) > 0:
            output += "Records with duplicate sequences:\n"
            for next_id in rep_seqs:
                if len(rep_seqs[next_id]) > 1:
                    output += "("
                    for seq_id in rep_seqs[next_id]:
                        output += "%s, " % seq_id
                    output = "%s), " % output.strip(", ")
                else:
                    output += "%s, " % next_id
            output = "%s\n\n" % output.strip(", ")
        else:
            output += "No records with duplicate sequences\n\n"
        if len(unique) > 0:
            output += "Unique records:\n"
            for next_id in unique:
                output += "%s, " % next_id
            output = "%s" % output.strip(", ")
        else:
            output += "No unique records"
        sys.stdout.write("%s\n" % output)

    # Pull sequence ends
    if in_args.pull_record_ends:
        amount, which_end = in_args.pull_record_ends
        amount = int(amount)
        new_seqs = pull_seq_ends(seqs, amount, which_end)
        _print_recs(new_seqs)

    # Extract regions
    if in_args.extract_region:
        start, end = in_args.extract_region
        seqs = extract_range(seqs, start, end)
        _print_recs(seqs)

    # Pull records
    if in_args.pull_records:
        search = in_args.pull_records
        records = pull_recs(seqs, search)
        _print_recs(records)

    # Hash sequence ids
    if in_args.hash_seq_ids:
        in_place_allowed = True
        hashed = hash_seqeunce_ids(seqs)
        hash_table = "# Hash table\n"
        for seq in hashed[0]:
            hash_table += "%s,%s\n" % (seq[0], seq[1])
        sys.stderr.write("%s\n" % hash_table)
        _print_recs(hashed[1])

    # Guess alphabet
    if in_args.guess_alphabet:
        if seqs.alpha == IUPAC.protein:
            sys.stdout.write("prot\n")
        elif seqs.alpha == IUPAC.ambiguous_dna:
            sys.stdout.write("dna\n")
        elif seqs.alpha == IUPAC.ambiguous_rna:
            sys.stdout.write("rna\n")
        else:
            sys.stdout.write("Undetermined\n")

    # Delete metadata
    if in_args.delete_metadata:
        in_place_allowed = True
        _print_recs(delete_metadata(seqs))

    # Raw Seq
    if in_args.raw_seq:
        seqs = clean_seq(seqs)
        output = ""
        for seq in seqs.seqs:
            output += "%s\n\n" % seq.seq
        sys.stdout.write("%s\n" % output.strip())

    # Clean Seq
    if in_args.clean_seq:
        in_place_allowed = True
        _print_recs(clean_seq(seqs))

    # Guess format
    if in_args.guess_format:
        for seq_set in in_args.sequence:
            sys.stdout.write("%s\t-->\t%s\n" % (seq_set, SeqBuddy(seq_set).in_format))

    # Map features from cDNA over to protein
    if in_args.map_features_dna2prot:
        in_place_allowed = True
        file1, file2 = in_args.sequence[:2]

        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)

        if file1.alpha == file2.alpha:
            sys.exit("Error: You must provide one DNA file and one protein file")

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
        in_place_allowed = True
        file1, file2 = in_args.sequence[:2]

        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)

        if file1.alpha == file2.alpha:
            sys.exit("Error: You must provide one DNA file and one protein file")

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
        in_place_allowed = True
        _print_recs(order_features_by_position(seqs))

    # Order sequence features alphabetically
    if in_args.order_features_alphabetically:
        in_place_allowed = True
        _print_recs(order_features_alphabetically(seqs))