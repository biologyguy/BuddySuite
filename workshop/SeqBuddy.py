#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: SeqBuddy.py
version: 1, beta
author: Stephen R. Bond
email: steve.bond@nih.gov
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

# Third party
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
from Bio.Nexus.Trees import TreeError

# BuddySuite specific
import MyFuncs
import buddy_resources as br


# ##################################################### WISH LIST #################################################### #
def sim_ident(matrix):  # Return the pairwise similarity and identity scores among sequences
    x = matrix
    return x


def predict_orfs():
    # Add all predicted open reading frames to seqrecord features list
    # http://www.ncbi.nlm.nih.gov/gorf/gorf.html
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


def incremental_rename(query, replace):
    # Append a number to the end of each replacement to ensure unique ids
    x = (query, replace)
    return x


# - Allow batch calls. E.g., if 6 files are fed in as input, run the SeqBuddy command independently on each
# - Add support for selecting individual sequences to modify (as a global ability for any tool)
# - Add FASTQ support... More generally, support letter annotation mods
# - Add Clustal support
# - Get BuddySuite into PyPi
# - Check on memory requirements before execution
# - Execution timer, for long running jobs
# - Sort out a good way to manage 'lazy' imports (might not be that important)

# ##################################################### CHANGE LOG ################################################### #
# ###################################################### GLOBALS ##################################################### #
VERSION = br.Version("SeqBuddy", 1, 'beta', br.contributors)
OUTPUT_FORMATS = ["ids", "accessions", "summary", "full-summary", "clustal", "embl", "fasta", "fastq", "fastq-sanger",
                  "fastq-solexa", "fastq-illumina", "genbank", "gb", "imgt", "nexus", "phd", "phylip", "phylipi",
                  "phylipis", "raw", "seqxml", "sff", "stockholm", "tab", "qual"]


# ##################################################### SEQBUDDY ##################################################### #
class SeqBuddy:  # Open a file or read a handle and parse, or convert raw into a Seq object
    def __init__(self, sb_input, in_format=None, out_format=None, alpha=None):
        # ####  IN AND OUT FORMATS  #### #
        # Holders for input type. Used for some error handling below
        in_handle = None
        raw_sequence = None
        in_file = None
        self.alpha = alpha

        # Handles
        if str(type(sb_input)) == "<class '_io.TextIOWrapper'>":
            if not sb_input.seekable():  # Deal with input streams (e.g., stdout pipes)
                temp = StringIO(sb_input.read())
                sb_input = temp
            sb_input.seek(0)
            in_handle = sb_input.read()
            sb_input.seek(0)

        # Raw sequences
        if in_format == "raw":
            in_format = "fasta"
            out_format = "fasta"
            if type(sb_input) == str:
                sb_input = [SeqRecord(Seq(sb_input), id="raw_input", description="")]
            else:
                sb_input = [SeqRecord(Seq(sb_input.read()), id="raw_input", description="")]

        # Plain text in a specific format
        if type(sb_input) == str and not os.path.isfile(sb_input):
            raw_sequence = sb_input
            temp = StringIO(sb_input)
            sb_input = temp
            sb_input.seek(0)

        # File paths
        try:
            if os.path.isfile(sb_input):
                in_file = sb_input
        except TypeError:  # This happens when testing something other than a string.
            pass

        if not in_format:
            self.in_format = _guess_format(sb_input)
            if self.in_format == "empty file":
                self.in_format = "fasta"
            self.out_format = str(self.in_format) if not out_format else out_format

        else:
            self.in_format = in_format

        if not self.in_format:
            if in_file:
                raise GuessError("Could not determine format from sb_input file '{0}'.\n"
                                 "Try explicitly setting with -f flag.".format(in_file))
            elif raw_sequence:
                raise GuessError("File not found, or could not determine format from raw input\n{0} ..."
                                 "Try explicitly setting with -f flag.".format(raw_sequence)[:60])
            elif in_handle:
                raise GuessError("Could not determine format from input file-like object\n{0} ..."
                                 "Try explicitly setting with -f flag.".format(in_handle)[:50])
            else:
                raise GuessError("Unable to determine format or input type. Please check how SeqBuddy is being called.")

        self.out_format = self.in_format if not out_format else out_format

        # ####  RECORDS  #### #
        if type(sb_input) == SeqBuddy:
            sequences = sb_input.records

        elif isinstance(sb_input, list):
            # make sure that the list is actually SeqIO records (just test a few...)
            rand_sample = sb_input if len(sb_input) < 5 else sample(sb_input, 5)
            for seq in rand_sample:
                if type(seq) != SeqRecord:
                    raise TypeError("Seqlist is not populated with SeqRecords.")
            sequences = sb_input

        elif str(type(sb_input)) == "<class '_io.TextIOWrapper'>" or isinstance(sb_input, StringIO):
            sequences = list(SeqIO.parse(sb_input, self.in_format))

        elif os.path.isfile(sb_input):
            with open(sb_input, "r") as sb_input:
                sequences = list(SeqIO.parse(sb_input, self.in_format))
        else:
            sequences = [SeqRecord(Seq(sb_input))]  # may be unreachable?

        if self.alpha is None:
            self.alpha = _guess_alphabet(sequences)
        elif self.alpha in ['protein', 'prot', 'p', 'pep', IUPAC.protein]:
            self.alpha = IUPAC.protein
        elif self.alpha in ['dna', 'd', 'cds', IUPAC.ambiguous_dna]:
            self.alpha = IUPAC.ambiguous_dna
        elif self.alpha in ['rna', 'r', IUPAC.ambiguous_rna]:
            self.alpha = IUPAC.ambiguous_rna
        else:
            _stderr("WARNING: Alphabet not recognized. Correct alphabet will be guessed.\n")
            self.alpha = _guess_alphabet(sequences)

        for seq in sequences:
            seq.seq.alphabet = self.alpha

        # The NEXUS parser adds '.copy' to any repeat taxa, strip that off...
        if self.in_format == "nexus":
            for rec in sequences:
                rec.id = re.sub("\.copy[0-9]*$", "", rec.id)

        self.records = sequences

    def to_dict(self):
        sb_copy = find_repeats(_make_copy(self))
        if len(sb_copy.repeat_ids) > 0:
            raise RuntimeError("There are repeat IDs in self.records\n%s" %
                               ", ".join([key for key, recs in sb_copy.repeat_ids.items()]))

        records_dict = OrderedDict()
        for rec in self.records:
            records_dict[rec.id] = rec
        return records_dict

    def print(self):
        print(str(self).strip())
        return

    def __str__(self):
        if len(self.records) == 0:
            return "Error: No sequences in object.\n"

        # There is a weird bug in genbank write() that concatenates dots to the organism name (if set).
        # The following is a work around...
        if self.out_format in ["gb", "genbank"]:
            for rec in self.records:
                try:
                    if re.search("(\. )+", rec.annotations['organism']):
                        rec.annotations['organism'] = "."
                except KeyError:
                    pass

        if self.out_format == "phylipi":
            output = _phylipi(self)

        elif self.out_format == "phylipis":
            output = _phylipi(self, "strict")

        elif self.out_format == "raw":
            output = "\n\n".join([str(rec.seq) for rec in self.records])
        else:
            tmp_dir = MyFuncs.TempDir()
            with open("%s/seqs.tmp" % tmp_dir.path, "w") as _ofile:
                try:
                    SeqIO.write(self.records, _ofile, self.out_format)
                except ValueError as e:
                    if "Sequences must all be the same length" in str(e):
                        _stderr("Warning: Alignment format detected but sequences are different lengths. "
                                "Format changed to fasta to accommodate proper printing of records.\n")
                        SeqIO.write(self.records, _ofile, "fasta")
                    elif "Repeated name" in str(e) and self.out_format == "phylip":
                        _stderr("Warning: Phylip format returned a 'repeat name' error, probably due to truncation. "
                                "Format changed to phylip-relaxed.\n")
                        SeqIO.write(self.records, _ofile, "phylip-relaxed")
                    else:
                        raise e

            with open("%s/seqs.tmp" % tmp_dir.path, "r") as ifile:
                output = ifile.read()

        return output

    def write(self, file_path):
        with open(file_path, "w") as ofile:
            ofile.write(str(self))
        return


# ################################################# HELPER FUNCTIONS ################################################# #
class GuessError(Exception):
    """Raised when input format cannot be guessed"""

    def __init__(self, _value):
        self.value = _value

    def __str__(self):
        return self.value


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

    current_dir = os.getcwd()
    script_location = os.path.realpath(__file__)
    script_location = re.sub('SeqBuddy\.py', '', script_location)

    prompt = MyFuncs.ask("%s binary not found. Try to download it? [yes]/no: " % blast_bin)
    if prompt:
        os.chdir(script_location)
        if _download_blast_binaries(**{blast_bin: True}):
            _stderr("%s downloaded.\n" % blast_bin)
        else:
            _stderr("Failed to download %s.\n" % blast_bin)
        os.chdir(current_dir)

    return False if not which(blast_bin) else True


def _download_blast_binaries(blastn=True, blastp=True, blastdcmd=True, **kwargs):
    """
    :param blastn: Install blastn? (bool)
    :param blastp: Install blastp? (bool)
    :param blastdcmd: Install blastdbcmd? (bool)
    :param kwargs: Just used for unit test options at the moment
                   ignore_pre_install: bool
                   system: {'darwin', 'linux', 'win'}
    :return: True on success, False on failure
    """

    if os.path.isdir(os.path.abspath('~/.buddysuite')) and "ignore_pre_install" not in kwargs:
        current_path = os.path.abspath('~/.buddysuite')
    else:
        current_path = os.getcwd()

    binary_source = 'https://raw.github.com/biologyguy/BuddySuite/master/workshop/build_dir/blast_binaries/'
    bins_to_dl = []
    current_os = sys.platform if "system" not in kwargs else kwargs["system"]
    if current_os.startswith('darwin'):
        if blastdcmd:
            bins_to_dl.append('Darwin_blastdbcmd.zip')
        if blastn:
            bins_to_dl.append('Darwin_blastn.zip')
        if blastp:
            bins_to_dl.append('Darwin_blastp.zip')
    elif current_os.startswith('linux'):
        system_bits = "32" if sys.maxsize < 3000000000 else "64"
        if blastdcmd:
            bins_to_dl.append('Linux_blastdbcmd%s.zip' % system_bits)
        if blastn:
            bins_to_dl.append('Linux_blastn%s.zip' % system_bits)
        if blastp:
            bins_to_dl.append('Linux_blastp%s.zip' % system_bits)
    elif current_os.startswith('win'):
        if blastdcmd:
            bins_to_dl.append('Win32_blastdbcmd.zip')
        if blastn:
            bins_to_dl.append('Win32_blastn.zip')
        if blastp:
            bins_to_dl.append('Win32_blastp.zip')
    else:
        return False

    file_to_name = {'Darwin_blastdbcmd.zip': 'blastdbcmd', 'Darwin_blastn.zip': 'blastn',
                    'Darwin_blastp.zip': 'blastp', 'Linux_blastdbcmd32.zip': 'blastdbcmd',
                    'Linux_blastn32.zip': 'blastn', 'Linux_blastp32.zip': 'blastp',
                    'Linux_blastdbcmd64.zip': 'blastdbcmd', 'Linux_blastn64.zip': 'blastn',
                    'Linux_blastp64.zip': 'blastp', 'Win32_blastdbcmd.zip': 'blastdbcmd.exe',
                    'Win32_blastn.zip': 'blastn.exe', 'Win32_blastp.zip': 'blastp.exe'}

    os.makedirs("{0}/__tempdir__".format(current_path), exist_ok=True)
    try:
        for blast_bin in bins_to_dl:
            with request.urlopen('{0}{1}'.format(binary_source, blast_bin)) as reader, \
                    open("{0}/__tempdir__/{1}".format(current_path, blast_bin), mode='wb') as writer:
                shutil.copyfileobj(reader, writer)
            zip_file = zipfile.ZipFile("{0}/__tempdir__/{1}".format(current_path, blast_bin))
            zip_file.extractall(path=current_path)
            if current_os.startswith('win'):
                os.rename('{0}/{1}'.format(current_path, re.sub('\.zip', '.exe', blast_bin)),
                          '{0}/{1}'.format(current_path, file_to_name[blast_bin]))
            else:
                os.rename('{0}/{1}'.format(current_path, re.sub('\.zip', '', blast_bin)),
                          '{0}/{1}'.format(current_path, file_to_name[blast_bin]))
            os.chmod('{0}/{1}'.format(current_path, file_to_name[blast_bin]), 0o755)
            print("File added: {0}/{1}".format(current_path, file_to_name[blast_bin]))
        shutil.rmtree("{0}/__tempdir__".format(current_path))
    except error.URLError:
        return False

    return True


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
        feature = _shift_features(feature, shift, seq_len)[0]
        feature.strand *= -1
    else:
        raise TypeError("_feature_rc requires a feature with either FeatureLocation or CompoundLocation, "
                        "not %s" % type(feature.location))
    return feature


def _guess_alphabet(seqbuddy):
    """
    Looks through the characters in the SeqBuddy records to determine the most likely alphabet
    Does not attempt to explicitly deal with weird cases (e.g., ambiguous residues).
    The user will need to specify an alphabet with the -a flag if using many non-standard characters in their sequences.
    :param seqbuddy: SeqBuddy object
    :return: IUPAC alphebet object
    """
    seq_list = seqbuddy if isinstance(seqbuddy, list) else seqbuddy.records
    seq_list = [str(x.seq) for x in seq_list]
    sequence = "".join(seq_list).upper()
    sequence = re.sub("[NX\-?]", "", sequence)

    if len(sequence) == 0:
        return None

    if 'U' in sequence:  # U is unique to RNA
        return IUPAC.ambiguous_rna

    percent_dna = len(re.findall("[ATCG]", sequence)) / float(len(sequence))
    percent_protein = len(re.findall("[ACDEFGHIKLMNPQRSTVWXY]", sequence)) / float(len(sequence))
    if percent_dna > 0.85:  # odds that a sequence with no Us and such a high ATCG count be anything but DNA is low
        return IUPAC.ambiguous_dna
    elif percent_protein > 0.85:
        return IUPAC.protein
    else:
        return None


def _guess_format(_input):
    """
    Loop through many possible formats that BioPython has a parser for, and return the format that is
    actually able to return records.
    :param _input: Duck-typed; can be list, SeqBuddy object, file handle, or file path.
    :return: str or None
    """
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
        if _input.read() == "":
            return "empty file"
        _input.seek(0)

        # ToDo: Glean CLUSTAL
        possible_formats = ["phylip-relaxed", "stockholm", "fasta", "gb", "fastq", "nexus", "embl", "seqxml"]
        for next_format in possible_formats:
            try:
                _input.seek(0)
                seqs = SeqIO.parse(_input, next_format)
                if next(seqs):
                    _input.seek(0)
                    return next_format
                else:
                    continue
            except StopIteration:  # ToDo check that other types of error are not possible
                continue
            except ValueError:
                continue
            except SAXParseException:  # Thrown by seqxml parser
                continue
            except TreeError:  # Thrown by NEXUS tree files
                continue
        return None  # Unable to determine format from file handle

    else:
        raise GuessError("Unsupported _input argument in guess_format(). %s" % _input)


def _make_copy(seqbuddy):
    """
    Deepcopy a SeqBuddy object. The alphabet objects are not handled properly when deepcopy is called,
    so need to wrap it.
    :param seqbuddy: SeqBuddy object
    :return: SeqBuddy object
    """
    alphabet_list = [rec.seq.alphabet for rec in seqbuddy.records]
    _copy = deepcopy(seqbuddy)
    _copy.alpha = seqbuddy.alpha
    for indx, rec in enumerate(_copy.records):
        rec.seq.alphabet = alphabet_list[indx]
    return _copy


def _phylipi(seqbuddy, _format="relaxed"):
    """
    Convert sequences to inline Phylip
    :param seqbuddy: SeqBuddy object
    :param _format: {"strict", "relaxed"}
    :return: The sequences formated in phylip format
    :ToDo: Throw errors when sequences are not the same length or truncation of ids leads to repeat taxa. Also
    put up a warning when ids are truncated.
    """
    max_id_length = 0
    max_seq_length = 0
    for rec in seqbuddy.records:
        max_id_length = len(rec.id) if len(rec.id) > max_id_length else max_id_length
        max_seq_length = len(rec.seq) if len(rec.seq) > max_seq_length else max_seq_length

    output = " %s %s\n" % (len(seqbuddy.records), max_seq_length)
    for rec in seqbuddy.records:
        seq_id = rec.id.ljust(max_id_length) if _format == "relaxed" else rec.id[:10].ljust(10)
        output += "%s  %s\n" % (seq_id, rec.seq)

    return output


def _shift_features(features, shift, full_seq_len):
    """
    Adjust the location of features
    :param features: Either a single SeqFeature object, or a list of them
    :param shift: int, how far the new feature should move from 0
    :param full_seq_len: The full length of the original sequence
    :return: List of SeqFeatures
    """
    if type(features) != list:  # Duck type for single feature input
        features = [features]

    shifted_features = []
    for feature in features:
        if type(feature.location) == CompoundLocation:  # Recursively call _shift_features() for compound locations
            new_compound_location = []
            for sub_feature in feature.location.parts:
                sub_feature = _shift_features(SeqFeature(sub_feature), shift, full_seq_len)
                if not sub_feature:
                    continue
                new_compound_location.append(sub_feature[0].location)

            if not new_compound_location:
                continue

            elif len(new_compound_location) == 1:
                feature.location = new_compound_location[0]

            else:
                feature.location = CompoundLocation(new_compound_location, feature.location.operator)

        elif type(feature.location) == FeatureLocation:
            start = feature.location.start + shift
            end = feature.location.end + shift
            if start > full_seq_len or end < 0:
                continue

            start = start if start >= 0 else 0
            end = end if end <= full_seq_len else full_seq_len

            feature.location = FeatureLocation(start, end, feature.strand)

        else:
            raise TypeError("_shift_feature requires a feature with either FeatureLocation or CompoundLocation, "
                            "not %s" % type(feature.location))
        shifted_features.append(feature)

    return shifted_features


def _stderr(message, quiet=False):
    """
    Send text to stderr
    :param message: Text to write
    :param quiet: Suppress message with True
    :return: None
    """
    if not quiet:
        sys.stderr.write(message)
        sys.stderr.flush()
    return


def _stdout(message, quiet=False):
    """
    Send text to stdout
    :param message: Text to write
    :param quiet: Suppress message with True
    :return: None
    """
    if not quiet:
        sys.stdout.write(message)
        sys.stdout.flush()
    return


# ################################################ MAIN API FUNCTIONS ################################################ #
def annotate(seqbuddy, _type, location, strand=None, qualifiers=None, pattern=None):
    """
    Adds a feature annotation to sequences in the SeqBuddy object
    :param seqbuddy: SeqBuddy object
    :param _type: Name/designation of the annotation
    :param location: The location of the new annotation in the sequence.
    If a single SeqFeature, use a tuple (start, end) or FeatureLocation object
    If a CompoundFeature, us a list of tuples [(start1, end1), (start2, end2)] or CompoundFeature object
    NOTE!!! If feeding in tuples, the 'start' index begins at 1, while Feature objects start at 0.
    :param strand: The feature's orientation (+/-/None)
    :param qualifiers: Further information to append to the new feature
    The argument can be a dictionary or a list ["foo: bar", "fizz: buzz"]
    :param pattern: List of regex patterns to specify which sequences to add the feature to
    :return: The updated SeqBuddy object
    """
    # http://www.insdc.org/files/feature_table.html
    old = _make_copy(seqbuddy)
    if pattern:
        recs = pull_recs(seqbuddy, pattern).records
    else:
        recs = seqbuddy.records

    for rec1 in recs:
        for rec2 in old.records:
            if rec1.id == rec2.id:
                old.records.remove(rec2)

    if isinstance(location, list) or isinstance(location, tuple):
        locations = []
        if isinstance(location[0], int):
            locations.append(FeatureLocation(start=location[0] - 1, end=location[1]))
        elif isinstance(location[0], tuple) or isinstance(location[0], list):
            for _tup in location:
                locations.append(FeatureLocation(start=_tup[0] - 1, end=_tup[1]))
        elif isinstance(location[0], str):
            for substr in location:
                substr = re.sub('[ ()]', '', substr)
                substr = re.sub('-|\.\.', ',', substr)
                locations.append(FeatureLocation(start=int(re.split(',', substr)[0]) - 1,
                                                 end=int(re.split(',', substr)[1])))
        location = CompoundLocation(sorted(locations, key=lambda x: x.start), operator='order') \
            if len(locations) > 1 else locations[0]
    elif isinstance(location, str):
        location = re.sub('[ ()]', '', location)
        location = re.split(',', location)
        locations = []
        for substr in location:
            locations.append(FeatureLocation(start=int(substr.split('-')[0]) - 1, end=int(substr.split('-')[1])))
        location = CompoundLocation(sorted(locations, key=lambda x: x.start), operator='order') \
            if len(locations) > 1 else locations[0]
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

    if seqbuddy.alpha == IUPAC.ambiguous_dna:
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
            _stderr("Warning: strand input not recognized. Value set to None.")
    else:
        strand = None
    qualifiers = [qualifiers] if isinstance(qualifiers, str) else qualifiers

    if isinstance(qualifiers, list):
        qual_dict = {}
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
    return sum_length / len(seqbuddy.records)


def back_translate(seqbuddy, mode='random', species=None):
    """
    Back-translates protein sequences into DNA sequences
    :param seqbuddy: SeqBuddy object
    :param mode: The codon selection mode (random/optimized)
    :param species: The model to use for optimized codon selection (human/mouse/yeast/ecoli)
    codon preference tables derived from the data at http://www.kazusa.or.jp
    :return: Modified SeqBuddy object
    """
    # Homo sapiens, species=9606
    if mode.upper() not in ['RANDOM', 'R', 'OPTIMIZED', 'O']:
        raise AttributeError("Back_translate modes accepted are 'random' or 'r' and 'optimized' or 'o'. "
                             "You entered '%s'" % mode)

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

    if seqbuddy.alpha != IUPAC.protein:
        raise TypeError("The input sequence needs to be IUPAC.protein'>, not %s" %
                        str(type(seqbuddy.alpha)))

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
    originals = _make_copy(seqbuddy)
    for rec in seqbuddy.records:
        rec.features = []
        dna_seq = ""
        for aa in rec.seq.upper():
            rand_num = random()
            sum_probs = 0.
            for i in range(len(lookup_table[aa][1])):
                sum_probs += lookup_table[aa][1][i]
                if sum_probs >= rand_num:
                    dna_seq += lookup_table[aa][0][i]
                    break
            rec.seq = Seq(dna_seq, alphabet=IUPAC.ambiguous_dna)

    mapped_featuresseqbuddy = map_features_prot2nucl(originals, seqbuddy, mode="list")
    mapped_featuresseqbuddy.out_format = seqbuddy.out_format
    return mapped_featuresseqbuddy


def bl2seq(seqbuddy):
    """
    Does an all-by-all analysis of the sequences
    :param seqbuddy: SeqBuddy object
    :return: OrderedDict of results dict[key][matches]
    """
    # Note on blast2seq: Expect (E) values are calculated on an assumed database size of (the rather large) nr, so the
    # threshold may need to be increased quite a bit to return short alignments

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

    if seqbuddy.alpha == IUPAC.protein and not _check_for_blast_bin("blastp"):
        raise RuntimeError("Blastp not present in $PATH or working directory.")

    elif seqbuddy.alpha in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna, IUPAC.ambiguous_rna, IUPAC.unambiguous_rna] \
            and not _check_for_blast_bin("blastn"):
        raise RuntimeError("Blastn not present in $PATH or working directory.")

    blast_bin = "blastp" if seqbuddy.alpha == IUPAC.protein else "blastn"

    from multiprocessing import Lock
    lock = Lock()
    tmp_dir = TemporaryDirectory()

    # Copy the seqbuddy records into new list, so they can be iteratively deleted below
    make_ids_unique(seqbuddy)
    seqs_copy = seqbuddy.records[:]
    subject_file = "%s/subject.fa" % tmp_dir.name
    for subject in seqbuddy.records:
        with open(subject_file, "w") as ifile:
            SeqIO.write(subject, ifile, "fasta")

        MyFuncs.run_multicore_function(seqs_copy, mc_blast, [subject_file], out_type=sys.stderr, quiet=True)
        seqs_copy = seqs_copy[1:]

    with open("%s/blast_results.txt" % tmp_dir.name, "r") as _ifile:
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


def blast(seqbuddy, blast_db):
    """
    ToDo: - Implement makeblastdb
          - Allow extra blast parameters

    Runs a BLAST search against a specified database, returning all significant matches.
    :param seqbuddy: SeqBuddy object
    :param blast_db: The location of the BLAST database to run the sequences against
    :return: A SeqBuddy object containing all of the BLAST database matches
    """
    blast_bin = "blastp" if seqbuddy.alpha == IUPAC.protein else "blastn"
    if not _check_for_blast_bin(blast_bin):
        raise SystemError("%s not found in system path." % blast_bin)

    if not _check_for_blast_bin("blastdbcmd"):
        raise SystemError("blastdbcmd not found in system path.")

    extensions = {"blastp": ["phr", "pin", "pog", "psd", "psi", "psq"],
                  "blastn": ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]}

    # Try to catch the common variations of the database names that might be given as input
    if blast_db[-2:] in [".p", ".n"]:
        blast_db = blast_db[:-2]

    if blast_db[-3:] in extensions[blast_bin]:
        blast_db = blast_db[:-4]

    blast_db = os.path.abspath(blast_db)

    # ToDo Check NCBI++ tools are a conducive version (2.2.29 and above, I think [maybe .28])
    # Check that complete blastdb is present and was made with the -parse_seqids flag
    for extension in extensions[blast_bin]:
        if not os.path.isfile("%s.%s" % (blast_db, extension)):
            raise RuntimeError("The .%s file of your blast database was not found. Ensure the -parse_seqids flag was "
                               "used with makeblastdb." % extension)

    seqbuddy = clean_seq(seqbuddy)  # in case there are gaps or something in the sequences

    tmp_dir = MyFuncs.TempDir()
    with open("%s/tmp.fa" % tmp_dir.path, "w") as ofile:
        SeqIO.write(seqbuddy.records, ofile, "fasta")

    Popen("%s -db %s -query %s/tmp.fa -out %s/out.txt -num_threads 4 -evalue 0.01 -outfmt 6" %
          (blast_bin, blast_db, tmp_dir.path, tmp_dir.path), shell=True).wait()

    with open("%s/out.txt" % tmp_dir.path, "r") as ifile:
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

    with open("%s/seqs.fa" % tmp_dir.path, "w") as ofile:
        for hit_id in hit_ids:
            hit = Popen("blastdbcmd -db %s -entry 'lcl|%s'" % (blast_db, hit_id), stdout=PIPE, shell=True).communicate()
            hit = hit[0].decode("utf-8")
            hit = re.sub("lcl\|", "", hit)
            ofile.write("%s\n" % hit)

    new_seqs = SeqBuddy("%s/seqs.fa" % tmp_dir.path)

    return new_seqs


def clean_seq(seqbuddy, ambiguous=True, rep_char="N", skip_list=None):
    """
    Removes all non-sequence characters, and converts ambiguous characters to 'X' if ambiguous=False
    :param seqbuddy: SeqBuddy object
    :param skip_list: Optional list of characters to be left alone
    :param ambiguous: Specifies whether ambiguous characters should be kept or not
    :return: The cleaned SeqBuddy object
    """
    skip_list = "" if not skip_list else "".join(skip_list)
    for rec in seqbuddy.records:
        if rec.seq.alphabet == IUPAC.protein:
            full_skip = "ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy%s" % skip_list
            rec.seq = Seq(re.sub("[^%s]" % full_skip, "", str(rec.seq)),
                          alphabet=rec.seq.alphabet)
        else:
            full_skip = "ATGCURYWSMKHBVDNXatgcurywsmkhbvdnx%s" % skip_list
            rec.seq = Seq(re.sub("[^%s]" % full_skip, "", str(rec.seq)),
                          alphabet=rec.seq.alphabet)
            if not ambiguous:
                full_skip = "ATGCUatgcu%s" % skip_list
                rec.seq = Seq(re.sub("[^%s]" % full_skip, rep_char, str(rec.seq)), alphabet=rec.seq.alphabet)
    return seqbuddy


def complement(seqbuddy):
    """
    Converts DNA/RNA sequences to their complementary sequence
    :param seqbuddy: SeqBuddy object
    :return: The modified SeqBuddy object
    """
    if seqbuddy.alpha == IUPAC.protein:
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
    for rec in seqbuddy.records:
        shift = len(new_seq)
        full_seq_len = len(new_seq) + len(str(rec.seq))
        rec.features = _shift_features(rec.features, shift, full_seq_len)

        location = FeatureLocation(len(new_seq), len(new_seq) + len(str(rec.seq)))
        feature = SeqFeature(location=location, id=rec.id, type=rec.id[:15])
        features.append(feature)
        features += rec.features
        concat_ids.append(rec.id)
        new_seq += str(rec.seq)

    new_seq = [SeqRecord(Seq(new_seq, alphabet=seqbuddy.alpha),
                         description="", id="concatination", features=features)]
    seqbuddy = SeqBuddy(new_seq)
    seqbuddy.out_format = "gb"
    return seqbuddy


def count_codons(seqbuddy):
    """
    Generate frequency statistics for codon composition
    :param seqbuddy: SeqBuddy object
    :return: A tuple containing the original SeqBuddy object and a dictionary - dict[id][codon] = (Amino acid, num, %)
    """
    if seqbuddy.alpha not in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna, IUPAC.ambiguous_rna, IUPAC.unambiguous_rna]:
        raise TypeError("Nucleic acid sequence required, not protein or other.")
    if seqbuddy.alpha in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna]:
        codontable = CodonTable.ambiguous_dna_by_name['Standard'].forward_table
    else:
        codontable = CodonTable.ambiguous_rna_by_name['Standard'].forward_table
    output = OrderedDict()
    for rec in seqbuddy.records:
        sequence = rec.seq
        if len(sequence) % 3 != 0:
            while len(sequence) % 3 != 0:
                sequence = sequence[:-1]
        data_table = {}
        num_codons = len(sequence) / 3
        while len(sequence) > 0:
            codon = str(sequence[:3]).upper()
            if codon in data_table:
                data_table[codon][1] += 1
            else:
                if codon.upper() in ['ATG', 'AUG']:
                    data_table[codon] = ['M', 1, 0.0]
                elif codon.upper() == 'NNN':
                    data_table[codon] = ['X', 1, 0.0]
                elif codon.upper() in ['TAA', 'TAG', 'TGA', 'UAA', 'UAG', 'UGA']:
                    data_table[codon] = ['*', 1, 0.0]
                else:
                    try:
                        data_table[codon] = [codontable[codon.upper()], 1, 0.0]
                    except KeyError:
                        _stderr("Warning: Codon '{0}' is invalid. Codon will be skipped.\n".format(codon))
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

        if seqbuddy.alpha is IUPAC.protein:
            ambig = len(re.findall("X", seq))
            if ambig > 0:
                resid_count['% Ambiguous'] = round(100 * ambig / seq_len, 2)

            pos = len(re.findall("[HKR]", seq))
            resid_count['% Positive'] = round(100 * pos / seq_len, 2)

            neg = len(re.findall("[DEC]", seq))
            resid_count['% Negative'] = round(100 * neg / seq_len, 2)

            neut = len(re.findall("[GAVLIPFYWSTNQM]", seq))
            resid_count['% Uncharged'] = round(100 * neut / seq_len, 2)

            hyrdophobic = len(re.findall("[AVLIPYFWMC]", seq))
            resid_count['% Hyrdophobic'] = round(100 * hyrdophobic / seq_len, 2)

            hyrdophilic = len(re.findall("[NQSTKRHDE]", seq))
            resid_count['% Hyrdophilic'] = round(100 * hyrdophilic / seq_len, 2)

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


def delete_features(seqbuddy, pattern):
    """
    Deletes features with IDs matching a regex pattern
    :param seqbuddy: SeqBuddy object
    :param pattern: The regex pattern to search with
    :return: The modified SeqBuddy object
    """
    for rec in seqbuddy.records:
        retained_features = []
        for _feature in rec.features:
            if not re.search(pattern, _feature.type):
                retained_features.append(_feature)
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
    temp_seqbuddy = SeqBuddy(">seq1\nATGCGCTAGCATGTCA", in_format="fasta")
    for rec in seqbuddy.records:
        rec.name = ''
        rec.description = ''
        rec.annotations = temp_seqbuddy.records[0].annotations
        rec.features = temp_seqbuddy.records[0].features
    return seqbuddy


def delete_records(seqbuddy, patterns):
    """
    Deletes records with IDs matching a regex pattern
    :param seqbuddy: SeqBuddy object
    :param patterns: A single regex pattern, or list of patterns, to search with
    :return: The modified SeqBuddy object
    """
    if type(patterns) == str:
        patterns = [patterns]
    if type(patterns) != list:
        raise ValueError("'patterns' must be a list or a string.")

    retained_records = []
    for pattern in patterns:
        deleted = [rec.id for rec in pull_recs(_make_copy(seqbuddy), pattern).records]
        for rec in seqbuddy.records:
            if rec.id in deleted:
                continue
            else:
                retained_records.append(rec)
        seqbuddy.records = retained_records
        retained_records = []
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
            for rep_id in seqbuddy.repeat_ids:
                store_one_copy = pull_recs(_make_copy(seqbuddy), "^%s$" % rep_id).records[0]
                delete_records(seqbuddy, "^%s$" % rep_id)
                seqbuddy.records.append(store_one_copy)

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
                    rep_seq = re.sub("([|.*?^\[\]()])", r"\\\1", rep_seq)
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


def dna2rna(seqbuddy):
    """
    Transcribes DNA into RNA
    :param seqbuddy: SeqBuddy object
    :return: Modified SeqBuddy object
    """
    if seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for rec in seqbuddy.records:
        rec.seq = Seq(str(rec.seq.transcribe()), alphabet=IUPAC.ambiguous_rna)
    seqbuddy.alpha = IUPAC.ambiguous_rna
    return seqbuddy


def extract_range(seqbuddy, start, end):
    """
    Retrieves subsequences in a specified range
    :param seqbuddy: SeqBuddy object
    :param start: The starting point
    :param end: The end point
    :return: The modified SeqBuddy object
    """
    start = 1 if int(start) < 1 else start
    # Don't use the standard index-starts-at-0... end must be left for the range to be inclusive
    start, end = int(start) - 1, int(end)
    if end < start:
        raise ValueError("The value given for end of range is smaller than for the start of range.")

    for rec in seqbuddy.records:
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
    return seqbuddy


def find_cpg(seqbuddy):
    """
    Predicts locations of CpG islands in DNA sequences
    :param seqbuddy: SeqBuddy object
    :return: Modified SeqBuddy object (buddy_data["cpgs"] appended to all records)
    """
    seqbuddy = clean_seq(seqbuddy)
    if seqbuddy.alpha not in [IUPAC.ambiguous_dna, IUPAC.unambiguous_dna]:
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


def find_pattern(seqbuddy, *patterns):  # TODO ambiguous letters mode
    """
    Finds ﻿occurrences of a sequence pattern
    :param seqbuddy: SeqBuddy object
    :param pattern: regex patterns
    :return: Annotated SeqBuddy object. The match indices are also stored in rec.buddy_data["find_patters"].
    """
    # search through sequences for regex matches. For example, to find micro-RNAs
    lowercase(seqbuddy)
    for pattern in patterns:
        for rec in seqbuddy.records:
            _add_buddy_data(rec, 'find_patterns')
            indices = []
            matches = re.finditer(pattern, str(rec.seq), flags=re.IGNORECASE)
            new_seq = ""
            last_match = 0
            for match in matches:
                indices.append(match.start())
                rec.features.append(SeqFeature(location=FeatureLocation(start=match.start(), end=match.end()),
                                               type='match', qualifiers={'regex': pattern, 'added_by': 'SeqBuddy'}))
                if match.start() > 0:
                    new_seq += str(rec.seq[last_match:match.start() - 1])
                new_seq += str(rec.seq[match.start():match.end()]).upper()
                last_match = match.end() + 1
            new_seq += str(rec.seq[last_match:])
            rec.seq = Seq(new_seq, alphabet=rec.seq.alphabet)

            if not rec.buddy_data['find_patterns']:
                rec.buddy_data['find_patterns'] = OrderedDict({pattern: indices})
            else:
                rec.buddy_data['find_patterns'][pattern] = indices
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
    seqbuddy_copy = _make_copy(seqbuddy)
    for rec in seqbuddy_copy.records:
        seq = str(rec.seq).encode()
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


# ToDo: Make sure cut sites are not already in the features list
def find_restriction_sites(seqbuddy, enzyme_group=(), min_cuts=1, max_cuts=None):
    """
    Finds the restriction sites in the sequences in the SeqBuddy object
    :param seqbuddy: SeqBuddy object
    :param enzyme_group: "commercial", "all", or a list of specific enzyme names
    :param min_cuts: The minimum cut threshold
    :param max_cuts: The maximum cut threshold
    :return: annotated SeqBuddy object, and a dictionary of restriction sites added as the `restriction_sites` attribute
    """
    if seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Unable to identify restriction sites in protein sequences.")
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
            for res in AllEnzymes:
                if str(res) not in blacklist:
                    batch.add(res)

        else:
            try:
                batch.add(enzyme)
            except ValueError:
                _stderr("Warning: %s not a known enzyme\n" % enzyme)

    sites = []
    for rec in seqbuddy.records:
        rec.res_sites = {}
        analysis = Analysis(batch, rec.seq)
        result = analysis.with_sites()
        for key, value in result.items():
            if key.cut_twice():
                _stderr("Warning: Double-cutters not supported.\n")
                pass
            elif min_cuts <= len(value) <= max_cuts:
                try:
                    for zyme in value:
                        cut_start = zyme + key.fst3 - 1
                        cut_end = zyme + key.fst5 + abs(key.ovhg) - 1
                        rec.features.append(SeqFeature(FeatureLocation(start=cut_start, end=cut_end), type=str(key)))
                except TypeError:
                    _stderr("Warning: No-cutters not supported.\n")
                    pass
                rec.res_sites[key] = value
        rec.res_sites = OrderedDict(sorted(rec.res_sites.items(), key=lambda x: x[0]))
        sites.append((rec.id, rec.res_sites))
    order_features_alphabetically(seqbuddy)
    seqbuddy.restriction_sites = sites
    return seqbuddy


def hash_sequence_ids(seqbuddy, hash_length=10):
    """
    Replaces the sequence IDs with random hashes
    :param seqbuddy: SeqBuddy object
    :param hash_length: Specifies the length of the new hashed IDs
    :return: The modified SeqBuddy object, with a new attribute `hash_map` added
    """
    hash_list = []
    seq_ids = []

    try:
        hash_length = int(hash_length)
    except ValueError:
        raise TypeError("Hash length argument must be an integer, not %s" % type(hash_length))

    if hash_length < 1:
        raise ValueError("Hash length must be greater than 0")

    if 32 ** hash_length <= len(seqbuddy.records) * 2:
        raise ValueError("Insufficient number of hashes available to cover all sequences. "
                         "Hash length must be increased.")

    for i in range(len(seqbuddy.records)):
        new_hash = ""
        seq_ids.append(seqbuddy.records[i].id)
        while True:
            new_hash = "".join([choice(string.ascii_letters + string.digits) for _ in range(hash_length)])
            if new_hash in hash_list:
                continue
            else:
                hash_list.append(new_hash)
                break
        seqbuddy.records[i].id = new_hash
        seqbuddy.records[i].name = new_hash

    hash_map = OrderedDict()
    for i in range(len(hash_list)):
        hash_map[hash_list[i]] = seq_ids[i]

    seqbuddy.hash_map = hash_map
    return seqbuddy


def insert_sequence(seqbuddy, sequence, location=0, regexes=None):
    """
    Add a specific sequence at a defined location in all records. E.g., adding a barcode (zero-indexed)
    :param seqbuddy: SeqBuddy object
    :param sequence: The sequence to be inserted
    :param location: The location to insert the sequences at
    :return: The modified SeqBuddy object
    """
    if regexes:
        recs_to_update = pull_recs(_make_copy(seqbuddy), regexes).to_dict()
    else:
        recs_to_update = seqbuddy.to_dict()

    for seq_id, rec in recs_to_update.items():
        if location >= 0:
            new_seq = str(rec.seq)[:location] + sequence + str(rec.seq)[location:]
        elif location == -1:
            new_seq = str(rec.seq) + sequence
        else:
            new_seq = str(rec.seq)[:location + 1] + sequence + str(rec.seq)[location + 1:]

        rec.seq = Seq(new_seq, rec.seq.alphabet)

    for indx, rec in enumerate(seqbuddy.records):
        if rec.id in recs_to_update:
            seqbuddy.records[indx] = recs_to_update[rec.id]

    return seqbuddy


def isoelectric_point(seqbuddy):
    """
    Calculate the isoelectric points
    :param seqbuddy: SeqBuddy object
    :return: SeqBuddy object with isoelectric point appended to each record as the last feature in the feature list
    """
    if seqbuddy.alpha is not IUPAC.protein:
        raise TypeError("Protein sequence required, not nucleic acid.")
    isoelectric_points = OrderedDict()
    for rec in seqbuddy.records:
        iso_point = ProteinAnalysis(str(rec.seq))
        iso_point = round(iso_point.isoelectric_point(), 10)
        isoelectric_points[rec.id] = iso_point
        rec.features.append(SeqFeature(location=FeatureLocation(start=1, end=len(rec.seq)), type='pI',
                                       qualifiers={'value': iso_point}))
    return seqbuddy


def lowercase(seqbuddy):
    """
    Converts all sequence characters to lowercase.
    :param seqbuddy: SeqBuddy object
    :return: The modified SeqBuddy object
    """
    for rec in seqbuddy.records:
        rec.seq = Seq(str(rec.seq).lower(), alphabet=rec.seq.alphabet)
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
                if re.search("\([^()]*\)", regex):
                    split = "".join(list(split.groups()))
                else:
                    split = split.group(0)
        if split != "":
            recs_by_identifier.setdefault(split, []).append(rec)
        else:
            recs_by_identifier["Unknown"].append(rec)

    if not recs_by_identifier["Unknown"]:
        del recs_by_identifier["Unknown"]

    new_seqbuddies = [(identifier, _make_copy(seqbuddy)) for identifier in recs_by_identifier]
    for identifier, sb in new_seqbuddies:
        sb.records = recs_by_identifier[identifier]
        sb.identifier = identifier
    new_seqbuddies = [sb for identifier, sb in new_seqbuddies]
    return new_seqbuddies


def make_ids_unique(seqbuddy):
    """
    Rename all repeat IDs
    Note: the edge case where a new ID creates a new conflict is not handled
    :param seqbuddy: SeqBuddy object
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
                rec.id = "%s-%s" % (rec.id, indx + 1)
        records += recs
    seqbuddy.records = records
    return seqbuddy


def map_features_nucl2prot(dnaseqbuddy, protseqbuddy, mode="key", quiet=False):
    """
    Applies cDNA/mRNA features to protein sequences
    :param dnaseqbuddy: A DNA SeqBuddy object with features to map
    :param protseqbuddy: A protein SeqBuddy
    :param mode: Specify how sequences should be matched up {list, key}
            - list = records are mapped in index order
            - key = records are converted to a dict and matched by key
    :param quiet: Suppress _stderr messages
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
                            "not %s" % type(feature.location))
        return feature

    protseqbuddy = clean_seq(protseqbuddy, skip_list="*")
    dnaseqbuddy = clean_seq(dnaseqbuddy)
    stderr_written = False
    if mode == "list":
        if len(protseqbuddy.records) != len(dnaseqbuddy.records):
            raise ValueError("The two input files do not contain the same number of sequences")

        record_map = list(zip(dnaseqbuddy.records, protseqbuddy.records))

    elif mode == "key":
        prot_dict = protseqbuddy.to_dict()
        nucl_dict = dnaseqbuddy.to_dict()

        record_map = []
        for seq_id, nucl_rec in nucl_dict.items():
            if seq_id not in prot_dict:
                stderr_written = True
                _stderr("Warning: %s is in the cDNA file, but not in the protein file\n" % seq_id, quiet)
                continue
            else:
                record_map.append((nucl_rec, prot_dict[seq_id]))

        for seq_id, prot_rec in prot_dict.items():
            if seq_id not in nucl_dict:
                stderr_written = True
                _stderr("Warning: %s is in the protein file, but not in the cDNA file\n" % seq_id, quiet)
                record_map.append((None, prot_rec))
    else:
        raise ValueError("'mode' must be either 'key' or 'position'.")

    for nucl_rec, prot_rec in record_map:
        # len(cds) or len(cds minus stop)
        if not nucl_rec:
            continue

        if len(prot_rec.seq) * 3 not in [len(nucl_rec.seq), len(nucl_rec.seq) - 3]:
            _stderr("Warning: size mismatch between aa and nucl seqs for %s --> %s, %s\n" %
                    (nucl_rec.id, len(nucl_rec.seq), len(prot_rec.seq)), quiet)

        prot_feature_hashes = []
        for feat in prot_rec.features:
            prot_feature_hashes.append(md5(str(feat).encode()).hexdigest())

        for feat in nucl_rec.features:
            feat = _feature_map(feat)
            if md5(str(feat).encode()).hexdigest() not in prot_feature_hashes:
                prot_rec.features.append(feat)

    if stderr_written:
        _stderr("\n", quiet)

    seqs_list = [prot for nucl, prot in record_map]
    seqbuddy = SeqBuddy(seqs_list)
    seqbuddy.out_format = dnaseqbuddy.out_format
    seqbuddy.in_format = dnaseqbuddy.in_format
    return seqbuddy


def map_features_prot2nucl(protseqbuddy, dnaseqbuddy, mode="key", quiet=False):
    """
    Applies protein features to cDNA/mRNA sequences
    :param protseqbuddy: A protein SeqBuddy object with features to map
    :param dnaseqbuddy: A DNA SeqBuddy
    :param mode: Specify how sequences should be matched up {list, key}
            - list = records are mapped in index order
            - key = records are converted to a dict and matched by key
    :param quiet: Suppress _stderr messages
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
                            "not %s" % type(feature.location))
        return feature

    protseqbuddy = clean_seq(protseqbuddy, skip_list="*")
    dnaseqbuddy = clean_seq(dnaseqbuddy)
    stderr_written = False
    if mode == "list":
        if len(protseqbuddy.records) != len(dnaseqbuddy.records):
            raise ValueError("The two input files do not contain the same number of sequences")

        record_map = list(zip(protseqbuddy.records, dnaseqbuddy.records))

    elif mode == "key":
        prot_dict = protseqbuddy.to_dict()
        nucl_dict = dnaseqbuddy.to_dict()

        record_map = []
        for seq_id, prot_rec in prot_dict.items():
            if seq_id not in nucl_dict:
                stderr_written = True
                _stderr("Warning: %s is in the protein file, but not in the cDNA file\n" % seq_id, quiet)
                continue
            else:
                record_map.append((prot_rec, nucl_dict[seq_id]))

        for seq_id, dna_rec in nucl_dict.items():
            if seq_id not in prot_dict:
                stderr_written = True
                _stderr("Warning: %s is in the cDNA file, but not in the protein file\n" % seq_id, quiet)
                record_map.append((None, dna_rec))

    else:
        raise ValueError("'mode' must be either 'key' or 'position'.")

    for prot_rec, nucl_rec in record_map:
        # len(cds) or len(cds minus stop)
        if not prot_rec:
            continue

        if len(prot_rec.seq) * 3 not in [len(nucl_rec.seq), len(nucl_rec.seq) - 3]:
            _stderr("Warning: size mismatch between aa and nucl seqs for %s --> %s, %s\n" %
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

    if stderr_written:
        _stderr("\n", quiet)

    seqs_list = [nucl for prot, nucl in record_map]
    seqbuddy = SeqBuddy(seqs_list)
    seqbuddy.out_format = protseqbuddy.out_format
    seqbuddy.in_format = protseqbuddy.in_format
    return seqbuddy


def merge(*seqbuddy):
    """
    Merges the feature lists of SeqBuddy objects
    :param *seqbuddy: SeqBuddy objects to be combined
    :return: A new SeqBuddy object
    """

    def merge_records(rec1, rec2):
        if str(rec1.seq) != str(rec2.seq):
            raise RuntimeError("Record mismatch: ID %s" % rec1.id)
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

    seqbuddy = SeqBuddy([rec for _id, rec in seq_dict.items()])
    seqbuddy = order_ids(seqbuddy)
    seqbuddy = order_features_by_position(seqbuddy)
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
    output = {'masses_ss': [], 'masses_ds': [], 'ids': []}
    aa_dict = amino_acid_weights
    if seqbuddy.alpha == IUPAC.protein:
        aa_dict = amino_acid_weights
    elif seqbuddy.alpha in [IUPAC.ambiguous_dna or IUPAC.unambiguous_dna]:
        aa_dict = deoxynucleotide_weights
        dna = True
    elif seqbuddy.alpha in [IUPAC.ambiguous_rna or IUPAC.unambiguous_rna]:
        aa_dict = deoxyribonucleotide_weights
    for rec in seqbuddy.records:
        rec.mass_ds = 0
        rec.mass_ss = 0
        if seqbuddy.alpha == IUPAC.protein:
            rec.mass_ss += 18.02  # molecular weight of a water molecule
        else:
            if dna:
                rec.mass_ss += 79.0  # molecular weight of 5' monophosphate in ssDNA
                rec.mass_ds += 157.9  # molecular weight of the 5' triphosphate in dsDNA
            else:
                rec.mass_ss += 159.0  # molecular weight of a 5' triphosphate in ssRNA
        for indx, value in enumerate(str(rec.seq).upper()):
            rec.mass_ss += aa_dict[value]
            if dna:
                rec.mass_ds += aa_dict[value] + deoxynucleotide_weights[deoxynucleotide_compliments[value]]
        output['masses_ss'].append(round(rec.mass_ss, 3))

        qualifiers = {}
        if seqbuddy.alpha == IUPAC.protein:
            qualifiers["peptide_value"] = round(rec.mass_ss, 3)
        elif dna:
            qualifiers["ssDNA_value"] = round(rec.mass_ss, 3)
            qualifiers["dsDNA_value"] = round(rec.mass_ds, 3)
            output['masses_ds'].append(round(rec.mass_ds, 3))
        elif seqbuddy.alpha in [IUPAC.ambiguous_rna or IUPAC.unambiguous_rna]:
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
    return len(seqbuddy.records)


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
    output = [(_rec.id, _rec) for _rec in seqbuddy.records]
    output = sorted(output, key=lambda x: x[0], reverse=reverse)
    output = [rec[1] for rec in output]
    seqbuddy.records = output
    return seqbuddy


def order_ids_randomly(seqbuddy):
    """
    Reorders seqbuddy.records. The order will always be changed if more than 2 recs are fed in.
    :param seqbuddy: SeqBuddy object
    :return: The reordered SeqBuddy object
    """
    if len(seqbuddy.records) < 2:
        return seqbuddy
    elif len(seqbuddy.records) == 2:
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

    output = []
    valve = MyFuncs.SafetyValve(global_reps=1000)
    while valve.step("order_ids_randomly() was unable to reorder your sequences. This shouldn't happen, so please"
                     "contact the developers to let then know about this error."):
        sb_copy = _make_copy(seqbuddy)
        for _ in range(len(sb_copy.records)):
            random_index = randint(1, len(sb_copy.records)) - 1
            output.append(sb_copy.records.pop(random_index))
        if ["%s%s" % (rec.id, rec.seq) for rec in seqbuddy.records] != ["%s%s" % (rec.id, rec.seq) for rec in output]:
            break
        output = []

    seqbuddy.records = output
    return seqbuddy


def pull_random_recs(seqbuddy, count=1):
    """
    Return a random record or subset of records (without replacement)
    :param seqbuddy: SeqBuddy object
    :param count (int): The number of random records to pull
    :return: The original SeqBuddy object with only the selected records remaining
    """
    count = abs(count) if abs(count) <= len(seqbuddy.records) else len(seqbuddy.records)
    random_recs = []
    for i in range(count):
        rand_index = randint(0, len(seqbuddy.records) - 1)
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
            rec.seq = Seq(str(rec.seq)[:amount], alphabet=rec.seq.alphabet)
            rec.features = _shift_features(rec.features, 0, len(str(rec.seq)))

        else:
            shift = -1 * (len(str(rec.seq)) + amount) if abs(amount) <= len(str(rec.seq)) else 0
            rec.features = _shift_features(rec.features, shift, len(str(rec.seq)))
            rec.seq = rec.seq[amount:]

        seq_ends.append(rec)

    seqbuddy.records = seq_ends
    return seqbuddy


def pull_recs(seqbuddy, search):
    """
    Retrieves sequences with names/IDs matching a search pattern
    :param seqbuddy: SeqBuddy object
    :param search: List of regex expressions or single regex
    :return: The modified SeqBuddy object
    """
    search = "|".join(search) if type(search) == list else search
    matched_records = []
    for rec in seqbuddy.records:
        if re.search(search, rec.description) or re.search(search, rec.id) or re.search(search, rec.name):
            matched_records.append(rec)
    seqbuddy.records = matched_records
    return seqbuddy


def purge(seqbuddy, threshold):
    """
    Deletes highly similar sequences
    ToDo: Implement a way to return a certain # of seqs (i.e. auto-determine threshold)
        - This would probably be a different flag in the UI
    :param seqbuddy: SeqBuddy object
    :param threshold: Sets the similarity threshold
    :return: The purged SeqBuddy object
    """
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


def rename(seqbuddy, query, replace="", num=0):
    """
    Rename sequence IDs
    :param seqbuddy: SeqBuddy object
    :param query: The pattern to be searched for
    :param replace: The string to be substituted
    :param num: The maximum number of substitutions to make
    :return: The modified SeqBuddy object
    """
    replace = re.sub("\s+", "_", replace)
    check_parentheses = re.findall("\([^()]*\)", query)
    check_replacement = re.findall(r"\\[0-9]+", replace)
    check_replacement = sorted([int(match[1:]) for match in check_replacement])
    if check_replacement and check_replacement[-1] > len(check_parentheses):
        raise AttributeError("There are more replacement match values specified than query parenthesized groups")

    for rec in seqbuddy.records:
        if num < 0:
            if check_replacement:
                for indx in sorted(range(check_replacement[-1]), reverse=True):
                    indx += 1
                    replace = re.sub(r"\\%s" % indx, r"\\%s" % (indx + 1), replace)
                right_replace = "\\%s" % (len(check_replacement) + 2)
            else:
                right_replace = "\\2"
            leftmost = str(rec.id)
            rightmost = ""
            hash_to_split_on = "UPNFSZ7FQ6RBhfFzwt0Cku4Yr1n2VvwVUG7x97G7"
            for _ in range(abs(num)):
                if leftmost == "":
                    break
                new_name = re.sub(r"(.*)%s(.*)" % query,
                                  r"\1%s%s%s" % (hash_to_split_on, replace, right_replace), leftmost, 1)
                new_name = new_name.split(hash_to_split_on)
                if len(new_name) == 2:
                    leftmost = new_name[0]
                    rightmost = new_name[1] + rightmost
                    new_name = leftmost + rightmost
                else:
                    new_name = leftmost + rightmost
                    break

        else:
            new_name = re.sub(query, replace, rec.id, num)
        rec.description = rec.description[len(rec.id) + 1:]
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
        rec.seq = Seq(new_seq, alphabet=rec.seq.alphabet)
    return seqbuddy


def reverse_complement(seqbuddy):
    """
    Converts DNA/RNA sequences to their reverse complementary sequence
    :param seqbuddy: SeqBuddy object
    :return: The modified SeqBuddy object
    """
    if seqbuddy.alpha == IUPAC.protein:
        raise TypeError("SeqBuddy object is protein. Nucleic acid sequences required.")
    for rec in seqbuddy.records:
        try:
            rec.seq = rec.seq.reverse_complement()
        except ValueError as e:
            if "Proteins do not have complements!" in str(e):
                raise TypeError("Record '%s' is protein. Nucleic acid sequences required." % rec.id)
            else:
                raise e  # Hopefully never gets here
        seq_len = len(rec.seq)
        shifted_features = [_feature_rc(feature, seq_len) for feature in rec.features]
        rec.features = shifted_features
    return seqbuddy


def rna2dna(seqbuddy):
    """
    Back-transcribes RNA into cDNA
    :param seqbuddy: SeqBuddy object
    :return: Modified SeqBuddy object
    """
    if seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Nucleic acid sequence required, not protein.")
    for rec in seqbuddy.records:
        rec.seq = Seq(str(rec.seq.back_transcribe()), alphabet=IUPAC.ambiguous_dna)
    seqbuddy.alpha = IUPAC.ambiguous_dna
    return seqbuddy


def select_frame(seqbuddy, frame, add_metadata=True):
    """
    Changes the reading frame of the sequences
    :param seqbuddy: SeqBuddy object
    :param frame: The reading frame to shift to
    :return: The shifted SeqBuddy object
    """
    def reset_frame(_rec, _residues):
        _rec.seq = Seq("%s%s" % (_residues, str(_rec.seq)), alphabet=_rec.seq.alphabet)
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

        _rec.features = _shift_features(_rec.features, len(_residues), len(_rec.seq))
        _rec.description = re.sub("\(frame[23][A-Za-z]{1,2}\)", "", _rec.description).strip()
        return _rec

    if seqbuddy.alpha == IUPAC.protein:
        raise TypeError("Select frame requires nucleic acid, not protein.")

    for rec in seqbuddy.records:
        for indx, feature in enumerate(rec.features):
            if feature.type == "frame_shift" and "residues" in feature.qualifiers:
                rec = reset_frame(rec, feature.qualifiers["residues"][0])
                del rec.features[indx]
                continue

            if feature.location.start + 1 < frame and add_metadata:
                feature.qualifiers["shift"] = [str(feature.location.start + 1 - frame)]

        check_description = re.search("\(frame[23]([A-Za-z]{1,2})\)", rec.description)
        if check_description:
            rec = reset_frame(rec, check_description.group(1))

        rec.features = _shift_features(rec.features, (frame - 1) * -1, len(rec.seq))
        if frame in [2, 3] and add_metadata:
            location = FeatureLocation(-1, 0)
            residues = str(rec.seq)[:frame - 1]
            rec.features.append(SeqFeature(location=location, type="frame_shift", strand=1,
                                           qualifiers={"residues": [residues]}))
            rec.description += " (frame%s%s)" % (frame, residues)
        rec.seq = Seq(str(rec.seq)[frame - 1:], alphabet=rec.seq.alphabet)
    return seqbuddy


def shuffle_seqs(seqbuddy):
    """
    Randomly reorder the residues in each sequence
    :param seqbuddy: SeqBuddy object
    :return: The shuffled SeqBuddy object
    """
    for rec in seqbuddy.records:
        tokens = []
        for letter in rec.seq:
            tokens.append(letter)
        new_seq = ''
        while len(tokens) > 0:
            rand_indx = randint(0, len(tokens) - 1)
            new_seq += tokens.pop(rand_indx)
        rec.seq = Seq(data=new_seq, alphabet=seqbuddy.alpha)
    return seqbuddy


def split_file(seqbuddy):
    """
    Split the records in a SeqBuddy object up into a collection of new SeqBuddy objects
    :param seqbuddy: SeqBuddy object
    :return: A list of SeqBuddy objects
    """
    sb_objs_list = []
    for rec in seqbuddy.records:
        sb = SeqBuddy([rec])
        sb.in_format = seqbuddy.in_format
        sb.out_format = seqbuddy.out_format
        sb_objs_list.append(sb)
    return sb_objs_list


def translate6frames(seqbuddy):
    """
    Translates a nucleotide sequence into a protein sequence across all six reading frames.
    :param seqbuddy: SeqBuddy object
    :return: The translated SeqBuddy object
    """
    frame1, frame2, frame3 = _make_copy(seqbuddy), _make_copy(seqbuddy), _make_copy(seqbuddy)
    seqbuddy = reverse_complement(seqbuddy)

    rframe1, rframe2, rframe3 = _make_copy(seqbuddy), _make_copy(seqbuddy), _make_copy(seqbuddy)

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
    for i in range(len(frame1.records)):
        frame1.records[i].id = "%s_f1" % frame1.records[i].id
        frame2.records[i].id = "%s_f2" % frame2.records[i].id
        frame3.records[i].id = "%s_f3" % frame3.records[i].id
        rframe1.records[i].id = "%s_rf1" % rframe1.records[i].id
        rframe2.records[i].id = "%s_rf2" % rframe2.records[i].id
        rframe3.records[i].id = "%s_rf3" % rframe3.records[i].id

        output += [frame1.records[i], frame2.records[i], frame3.records[i],
                   rframe1.records[i], rframe2.records[i], rframe3.records[i]]

    seqbuddy = SeqBuddy(output, out_format=seqbuddy.out_format)
    return seqbuddy


def translate_cds(seqbuddy, quiet=False):
    """
    Translates a nucleotide sequence into a protein sequence.
    :param seqbuddy: SeqBuddy object
    :param quiet: Suppress the errors thrown by translate(cds=True)
    :return: The translated SeqBuddy object
    """
    def trans(in_seq):
        try:
            in_seq.seq = in_seq.seq.translate(cds=True, to_stop=True)
            return in_seq

        except TranslationError as e1:
            _stderr("Warning: %s in %s\n" % (e1, in_seq.id), quiet)
            return e1

        except ValueError as e1:
            if "Proteins cannot be translated" in str(e1):
                raise TypeError("Record %s is protein." % in_seq.id)
            else:
                raise e1  # Hopefully never get here.

    clean_seq(seqbuddy, ambiguous=False, rep_char="N")
    replace_subsequence(seqbuddy, "x", "N")
    for rec in seqbuddy.records:  # Removes any frame shift annotations before translating
        for indx, feature in enumerate(rec.features):
            if feature.type == "frame_shift" and "residues" in feature.qualifiers:
                del rec.features[indx]
                break
        rec.description = re.sub("\(frame[23]([A-Za-z]{1,2})\)", "", rec.description).strip()

    translated_sb = _make_copy(seqbuddy)
    for rec in translated_sb.records:
        rec.features = []
        # Modify a copy of the sequence as needed to complete the cds translation
        temp_seq = deepcopy(rec)
        temp_seq.seq.alphabet = rec.seq.alphabet
        valve = MyFuncs.SafetyValve(global_reps=1000)
        while valve.test(str(temp_seq.seq), str(temp_seq)):
            test_trans = trans(temp_seq)
            # success
            if str(type(test_trans)) == "<class 'Bio.SeqRecord.SeqRecord'>":
                break

            # not standard length
            if re.search("Sequence length [0-9]+ is not a multiple of three", str(test_trans)):
                temp_seq.seq = Seq(str(temp_seq.seq)[:(len(str(temp_seq.seq)) - len(str(temp_seq.seq)) % 3)],
                                   alphabet=temp_seq.seq.alphabet)
                rec.seq = Seq(str(rec.seq)[:(len(str(rec.seq)) - len(str(rec.seq)) % 3)],
                              alphabet=rec.seq.alphabet)
                continue

            # not a start codon
            if re.search("First codon '[A-Za-z]{3}' is not a start codon", str(test_trans)):
                temp_seq.seq = Seq("ATG" + str(temp_seq.seq)[3:], alphabet=temp_seq.seq.alphabet)
                continue

            # not a stop codon
            if re.search("Final codon '[A-Za-z]{3}' is not a stop codon", str(test_trans)):
                temp_seq.seq = Seq(str(temp_seq.seq) + "TGA", alphabet=temp_seq.seq.alphabet)
                continue

            # non-standard characters, this should be unreachable.
            if re.search("Codon '[A-Za-z]{3}' is invalid", str(test_trans)):
                regex = re.findall("Codon '([A-Za-z]{3})' is invalid", str(test_trans))
                regex = "(?i)%s" % regex[0]
                temp_seq.seq = Seq(re.sub(regex, "NNN", str(temp_seq.seq), count=1), alphabet=temp_seq.seq.alphabet)
                rec.seq = Seq(re.sub(regex, "NNN", str(rec.seq), count=1), alphabet=rec.seq.alphabet)
                continue

            # internal stop codon(s) found
            if re.search("Extra in frame stop codon found", str(test_trans)):
                for i in range(round(len(str(temp_seq.seq)) / 3) - 1):
                    codon = str(temp_seq.seq)[(i * 3):(i * 3 + 3)]
                    if codon.upper() in ["TGA", "TAG", "TAA", "UGA", "UAG", "UAA"]:
                        new_seq = str(temp_seq.seq)[:(i * 3)] + "NNN" + str(temp_seq.seq)[(i * 3 + 3):]
                        temp_seq.seq = Seq(new_seq, alphabet=temp_seq.seq.alphabet)
                continue

            # Should be unreachable
            raise RuntimeError("There were uncaught translation exceptions\n%s" % test_trans)

        rec.seq = rec.seq.translate()
        rec.seq.alphabet = IUPAC.protein

    output = map_features_nucl2prot(seqbuddy, translated_sb, mode="list", quiet=quiet)
    output.out_format = seqbuddy.out_format
    seqbuddy = output
    return seqbuddy


def uppercase(seqbuddy):
    """
    Converts all sequence characters to uppercase.
    :param seqbuddy: SeqBuddy object
    :return: The modified SeqBuddy object
    """
    for rec in seqbuddy.records:
        rec.seq = Seq(str(rec.seq).upper(), alphabet=rec.seq.alphabet)
    return seqbuddy


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
''',)

    br.flags(parser, ("sequence", "Supply file path(s) or raw sequence. If piping sequences "
                                  "into SeqBuddy this argument can be left blank."),
             br.sb_flags, br.sb_modifiers, VERSION)

    in_args = parser.parse_args()

    seqbuddy = []
    seq_set = ""

    try:
        for seq_set in in_args.sequence:
            if isinstance(seq_set, TextIOWrapper) and seq_set.buffer.raw.isatty():
                _stderr("Warning: No input detected. Process will be aborted.")
                sys.exit()
            seq_set = SeqBuddy(seq_set, in_args.in_format, in_args.out_format, in_args.alpha)
            seqbuddy += seq_set.records

        seqbuddy = SeqBuddy(seqbuddy, seq_set.in_format, seq_set.out_format, seq_set.alpha)
    except GuessError as e:
        _stderr("GuessError: %s\n" % e, in_args.quiet)
        sys.exit()

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

    def _exit(tool, skip=skip_exit):
        if skip:
            return
        usage = br.Usage()
        usage.increment("SeqBuddy", VERSION.short(), tool)
        usage.save()
        sys.exit()

    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Add feature
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
            _stderr("Warning: The provided annotation type is not part of the GenBank format standard\n\n",
                    in_args.quiet)

        if len(ftype) > 16:
            _stderr("Warning: Feature type is longer than 16 characters and "
                    "will be truncated if printed to GenBank/EMBL format\n\n", in_args.quiet)

        flocation = in_args.annotate[1]

        if len(in_args.annotate) >= 3:
            for next_arg in in_args.annotate[2:]:
                feature_attrs = duck_type(next_arg, **feature_attrs)

            if not feature_attrs["qualifiers"]:
                feature_attrs["qualifiers"] = None
            if not feature_attrs["pattern"]:
                feature_attrs["pattern"] = None

        seqbuddy = annotate(seqbuddy, ftype, flocation, **feature_attrs)
        if in_args.out_format:
            seqbuddy.out_format = in_args.out_format
        _print_recs(seqbuddy)
        _exit("annotate")

    # Average length of sequences
    if in_args.ave_seq_length:
        clean = False if not in_args.ave_seq_length[0] or in_args.ave_seq_length[0] != "clean" else True
        _stdout("%s\n" % round(ave_seq_length(seqbuddy, clean), 2))
        _exit("ave_seq_length")

    # Back Transcribe
    if in_args.back_transcribe:
        if seqbuddy.alpha != IUPAC.ambiguous_rna:
            _raise_error(ValueError("You need to provide an RNA sequence."), "back_transcribe")
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

        if seqbuddy.alpha != IUPAC.protein:
            _raise_error(TypeError("The input sequence needs to be protein, not nucleotide"), "back_translate")

        _print_recs(back_translate(seqbuddy, mode, species))
        _exit("back_translate")

    # BL2SEQ
    if in_args.bl2seq:
        if len(find_repeats(seqbuddy).repeat_ids):
            _stderr("Warning: There are records with duplicate ids which will be renamed.\n", quiet=in_args.quiet)
        try:
            output_dict = bl2seq(seqbuddy)
            output_str = "#query\tsubject\t%_ident\tlength\tevalue\tbit_score\n"
            ids_already_seen = []
            for query_id, query_values in output_dict.items():
                ids_already_seen.append(query_id)
                query_values = [(key, value) for key, value in query_values.items()]
                query_values = sorted(query_values, key=lambda l: l[0])
                for subj_id, subj_values in query_values:
                    if subj_id in ids_already_seen:
                        continue

                    ident, length, evalue, bit_score = subj_values
                    output_str += "%s\t%s\t%s\t%s\t%s\t%s\n" % (query_id, subj_id, ident, length, evalue, bit_score)
            _stdout(output_str)
        except RuntimeError as e:
            _raise_error(e, "bl2seq", "not present in $PATH or working directory")
        _exit("bl2seq")

    # BLAST
    if in_args.blast:
        try:
            blast_res = blast(seqbuddy, in_args.blast)
            if len(blast_res.records) > 0:
                _print_recs(blast_res)
            else:
                _stdout("No significant matches found\n")
            _exit("blast")

        except (RuntimeError, SystemError) as e:
            _raise_error(e, "blast")

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
            output = ""
            for sequence_id in codon_table:
                output += '#### {0} ####\n'.format(sequence_id)
                output += 'Codon\tAA\tNum\tPercent\n'
                for codon in codon_table[sequence_id]:
                    data = codon_table[sequence_id][codon]
                    output += '{0}\t{1}\t{2}\t{3}\n'.format(codon, data[0], data[1], data[2])
                output += '\n'
            _stdout(output)
        except TypeError as e:
            _raise_error(e, "count_codons", "Nucleic acid sequence required, not protein or other.")
        _exit("count_codons")

    # Count residues
    if in_args.count_residues:
        if in_args.count_residues[0] and str(in_args.count_residues[0].lower()) in "concatenate":
            seqbuddy = concat_seqs(seqbuddy)

        count_residues(seqbuddy)
        output = ""
        for rec in seqbuddy.records:
            output += "%s\n" % str(rec.id)
            for residue, counts in rec.buddy_data["res_count"].items():
                try:
                    output += "{0}:\t{1}\t{2} %\n".format(residue, counts[0], round(counts[1] * 100, 2))
                except TypeError:
                    output += "{0}:\t{1}\n".format(residue, counts)
            output += "\n"
        _stdout(output)
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
        try:
            if len(in_args.delete_records) == 1:
                columns = 1
            else:
                columns = int(in_args.delete_records[-1])
                del in_args.delete_records[-1]
        except ValueError:
            columns = 1

        deleted_seqs = []
        for next_pattern in in_args.delete_records:
            deleted_seqs += pull_recs(_make_copy(seqbuddy), next_pattern).records

        seqbuddy = delete_records(seqbuddy, in_args.delete_records)

        if len(deleted_seqs) > 0 and not in_args.quiet:
            counter = 1
            raw = "# ####################### Deleted records ######################## #\n"
            for seq in deleted_seqs:
                raw += "%s\t" % seq.id
                if counter % columns == 0:
                    raw = "%s\n" % raw.strip()
                counter += 1
            raw = "%s\n# ################################################################ #\n" % raw.strip()
            _stderr(raw, in_args.quiet)

        if len(deleted_seqs) == 0:
            stderr_out = "# ################################################################ #\n"
            stderr_out += "# No sequence identifiers match %s\n" % ", ".join(in_args.delete_records)
            stderr_out += "# ################################################################ #\n"
            _stderr(stderr_out, in_args.quiet)
        _print_recs(seqbuddy)
        _exit("delete_records")

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
        stderr_output = ""
        if len(seqbuddy.repeat_ids) > 0 and scope in ["all", "ids"]:
            stderr_output += "# Records with duplicate ids deleted\n"
            counter = 1
            for seq in seqbuddy.repeat_ids:
                stderr_output += "%s\t" % seq
                if counter % columns == 0:
                    stderr_output = "%s\n" % stderr_output.strip()
                counter += 1
            stderr_output = "%s\n\n" % stderr_output.strip()
            seqbuddy = delete_repeats(seqbuddy, 'ids')
            find_repeats(seqbuddy)

        rep_seq_ids = []
        for seq in seqbuddy.repeat_seqs:
            rep_seq_ids.append([])
            for rep_seq_id in seqbuddy.repeat_seqs[seq]:
                rep_seq_ids[-1].append(rep_seq_id)

        if len(rep_seq_ids) > 0 and scope in ["all", "seqs"]:
            stderr_output += "# Records with duplicate sequence deleted\n"
            counter = 1
            for repeat_seqs in rep_seq_ids:
                stderr_output += "["
                for rep_seq in repeat_seqs:
                    stderr_output += "%s, " % rep_seq
                stderr_output = "%s], " % stderr_output.strip(", ")
                if counter % columns == 0:
                    stderr_output = "%s\n" % stderr_output.strip(", ")
                counter += 1
            stderr_output = "%s\n" % stderr_output.strip(", ")

        if stderr_output != "":
            stderr_output = "# ################################################################ #\n%s\n" \
                            "# ################################################################ #\n\n" \
                            % stderr_output.strip()
            _stderr(stderr_output, in_args.quiet)

        else:
            _stderr("No duplicate records found\n", in_args.quiet)

        _print_recs(delete_repeats(seqbuddy, 'seqs'))
        _exit("delete_repeats")

    # Delete sequences below threshold
    if in_args.delete_small:
        _print_recs(delete_small(seqbuddy, in_args.delete_small))
        _exit("delete_small")

    # Extract regions
    if in_args.extract_region:
        try:
            extract_range(seqbuddy, *in_args.extract_region)
            _print_recs(seqbuddy)
        except ValueError as e:
            _raise_error(e, "extract_region", "The value given for end of range is smaller than for the start")

        _exit("extract_region")

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
                output = ""
                for rec in seqbuddy.records:
                    if rec.buddy_data["cpgs"]:
                        value = ["%s-%s" % (x[0], x[1]) for x in rec.buddy_data["cpgs"]]
                        output += "{0}: {1}\n".format(rec.id, ", ".join(value))

                _stderr('########### Islands identified ###########\n%s\n'
                        '##########################################\n\n' %
                        output.strip(), in_args.quiet)
            else:
                _stderr("# No Islands identified\n\n", in_args.quiet)
            _print_recs(seqbuddy)
            _exit("find_CpG")

        except TypeError as e:
            _raise_error(e, "find_CpG", "DNA sequence required, not protein or RNA.")

    # Find pattern
    if in_args.find_pattern:
        find_pattern(seqbuddy, *in_args.find_pattern)
        for pattern in in_args.find_pattern:
            output = ""
            num_matches = 0
            for rec in seqbuddy.records:
                indices = rec.buddy_data['find_patterns'][pattern]
                if not len(indices):
                    output += "{0}: None\n".format(rec.id)
                else:
                    output += "{0}: {1}\n".format(rec.id, ", ".join([str(x) for x in indices]))
                    num_matches += len(indices)

            _stderr("#### {0} matches found across {1} sequences for "
                    "pattern '{2}' ####\n".format(num_matches, len(seqbuddy.records), pattern), in_args.quiet)
            _stderr("%s\n" % output, in_args.quiet)
        _print_recs(seqbuddy)
        _exit("find_pattern")

    # Find repeat sequences or ids
    if in_args.find_repeats:
        columns = 1 if not in_args.find_repeats[0] else in_args.find_repeats[0]
        find_repeats(seqbuddy)

        output_str = ""
        if len(seqbuddy.repeat_ids) > 0:
            output_str += "#### Records with duplicate IDs: ####\n"
            counter = 1
            for next_id in seqbuddy.repeat_ids:
                output_str += "%s\t" % next_id
                if counter % columns == 0:
                    output_str = "%s\n" % output_str.strip()
                counter += 1

            output_str = "%s\n\n" % output_str.strip()

        else:
            output_str += "#### No records with duplicate IDs ####\n\n"

        if len(seqbuddy.repeat_seqs) > 0:
            output_str += "#### Records with duplicate sequences: ####\n"
            counter = 1
            for next_id in seqbuddy.repeat_seqs:
                output_str += "["
                for seq_id in seqbuddy.repeat_seqs[next_id]:
                    output_str += "%s, " % seq_id
                output_str = "%s], " % output_str.strip(", ")

                if counter % columns == 0:
                    output_str = "%s\n" % output_str.strip(", ")

                counter += 1

            output_str = "%s\n\n" % output_str.strip(", ")
        else:
            output_str += "#### No records with duplicate sequences ####\n\n"

        output_str = "{0}\n".format(output_str.strip())

        _stdout(output_str)
        _exit("find_repeats")

    # Find restriction sites
    if in_args.find_restriction_sites:
        min_cuts, max_cuts, _enzymes, order = None, None, [], 'position'
        if not in_args.out_format:
            seqbuddy.out_format = "gb"

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
                    max_cuts = param if param >= min_cuts else min_cuts
                    min_cuts = param if param <= min_cuts else min_cuts
                elif param > max_cuts:
                    max_cuts = param
                elif param < min_cuts:
                    min_cuts = param

            elif param in ['alpha', 'position']:
                order = param
            else:
                _enzymes.append(param)

        _enzymes = ["commercial"] if len(_enzymes) == 0 else _enzymes
        max_cuts = min_cuts if min_cuts and not max_cuts else max_cuts
        min_cuts = 1 if not min_cuts else min_cuts

        clean_seq(seqbuddy)
        try:
            find_restriction_sites(seqbuddy, tuple(_enzymes), min_cuts, max_cuts)
        except TypeError as e:
            _raise_error(e, "find_restriction_sites")

        output = '# ### Restriction Sites (indexed at cut-site) ### #\n'
        for tup in seqbuddy.restriction_sites:
            output += "{0}\n".format(tup[0])
            restriction_list = tup[1]
            restriction_list = [[key, value] for key, value in restriction_list.items()]
            restriction_list = sorted(restriction_list, key=lambda l: str(l[0])) if order == 'alpha' else \
                sorted(restriction_list, key=lambda l: l[1])

            for _enzyme in restriction_list:
                cut_sites = [str(x) for x in _enzyme[1]]
                output += "{0}\t{1}\n".format(_enzyme[0], ", ".join(cut_sites))
            output += "\n"
        output = "%s\n# ############################################### #\n\n" % output.strip()
        _stderr(output, quiet=in_args.quiet)
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
                    split_patterns.append(arg)

        sp = ["-"] if not split_patterns and not num_chars else split_patterns

        taxa_groups = make_groups(seqbuddy, split_patterns=sp, num_chars=num_chars)
        if "".join(split_patterns) != "" and len(taxa_groups) == len(seqbuddy.records):
            taxa_groups = make_groups(seqbuddy, num_chars=5)

        for next_seqbuddy in taxa_groups:
            in_args.sequence[0] = "%s/%s.%s" % (out_dir, next_seqbuddy.identifier,
                                                br.format_to_extension[next_seqbuddy.out_format])
            _stderr("New file: %s\n" % in_args.sequence[0], check_quiet)
            open(in_args.sequence[0], "w").close()
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

        if not regexes:
            in_args.quiet = False
            _raise_error(ValueError("You must provide at least one regular expression."), "group_by_regex")

        taxa_groups = make_groups(seqbuddy, regex=regexes)

        for next_seqbuddy in taxa_groups:
            in_args.sequence[0] = "%s/%s.%s" % (out_dir, next_seqbuddy.identifier,
                                                br.format_to_extension[next_seqbuddy.out_format])
            _stderr("New file: %s\n" % in_args.sequence[0], check_quiet)
            open(in_args.sequence[0], "w").close()
            _print_recs(next_seqbuddy)

        _exit("group_by_regex")

    # Guess alphabet
    if in_args.guess_alphabet:
        output = ""
        for seq_set in in_args.sequence:
            if str(type(seq_set)) != "<class '_io.TextIOWrapper'>":
                try:
                    seqbuddy = SeqBuddy(seq_set)
                except Exception:  # This should NOT be made more specific. If it throws errors, it's unknown.
                    seqbuddy = SeqBuddy("", in_format="raw")
                seq_set = seq_set.split("/")[-1]
                output += "%s\t-->\t" % seq_set
            else:
                output += "PIPE\t-->\t"

            if seqbuddy.alpha == IUPAC.protein:
                output += "prot\n"
            elif seqbuddy.alpha == IUPAC.ambiguous_dna:
                output += "dna\n"
            elif seqbuddy.alpha == IUPAC.ambiguous_rna:
                output += "rna\n"
            else:
                output += "Undetermined\n"
        _stdout(output)
        _exit("guess_alphabet")

    # Guess format
    if in_args.guess_format:
        output = ""
        for seq_set in in_args.sequence:
            if str(type(seq_set)) == "<class '_io.TextIOWrapper'>":
                if seqbuddy.in_format:
                    output += "PIPE\t-->\t%s\n" % seqbuddy.in_format
                else:
                    output += "PIPE\t-->\tUnknown\n"
            else:
                try:
                    file_format = _guess_format(seq_set)
                except Exception:  # This should NOT be made more specific. If it throws errors, it's unknown.
                    file_format = None

                seq_set = seq_set.split("/")[-1]
                if not file_format:
                    output += "%s\t-->\tUnknown\n" % seq_set
                else:
                    output += "%s\t-->\t%s\n" % (seq_set, file_format)
        _stdout(output)
        _exit("guess_format")

    # Hash sequence ids
    if in_args.hash_seq_ids:
        if in_args.hash_seq_ids[0] == 0:
            hash_length = 0
        elif not in_args.hash_seq_ids[0]:
            hash_length = 10
        else:
            hash_length = in_args.hash_seq_ids[0]

        if hash_length < 1:
            _stderr("Warning: The hash_length parameter was passed in with the value %s. This is not a positive "
                    "integer, so the hash length as been set to 10.\n\n" % hash_length, quiet=in_args.quiet)
            hash_length = 10

        if 32 ** hash_length <= len(seqbuddy.records) * 2:
            holder = ceil(log(len(seqbuddy.records) * 2, 32))
            _stderr("Warning: The hash_length parameter was passed in with the value %s. This is too small to properly "
                    "cover all sequences, so it has been increased to %s.\n\n" % (hash_length, holder), in_args.quiet)
            hash_length = holder

        hash_sequence_ids(seqbuddy, hash_length)

        hash_table = "# Hash table\n"
        for _hash, orig_id in seqbuddy.hash_map.items():
            hash_table += "%s,%s\n" % (_hash, orig_id)
        hash_table += "\n"

        _stderr(hash_table, in_args.quiet)
        _print_recs(seqbuddy)
        _exit("hash_seq_ids")

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

    # Isoelectric Point
    if in_args.isoelectric_point:
        if seqbuddy.alpha != IUPAC.protein:
            _stderr("Nucleic acid sequences detected, converting to protein.\n\n")
            seqbuddy = translate_cds(seqbuddy, quiet=True)

        isoelectric_point(seqbuddy)
        _stderr("ID\tpI\n")
        output = ""
        for rec in seqbuddy.records:
            output += "{0}\t{1}\n".format(rec.id, round(rec.features[-1].qualifiers["value"], 3))
        _stdout(output)
        _exit("isoelectric_point")

    # List features
    if in_args.list_features:
        output = ""
        for rec in seqbuddy.records:
            output += '#### {0} ####\n'.format(rec.id)
            if len(rec.features) > 0:
                for feat in rec.features:
                    output += '{0}\n'.format(feat.type)
                    if isinstance(feat.location, CompoundLocation):
                        output += '\tLocations:\n'
                        for part in feat.location.parts:
                            output += '\t\t{0}-{1}\n'.format(part.start, part.end)
                    elif isinstance(feat.location, FeatureLocation):
                        output += '\tLocation:\n'
                        output += '\t\t{0}-{1}\n'.format(feat.location.start, feat.location.end)
                    if feat.strand == 1:
                        output += '\tStrand: Sense(+)\n'
                    elif feat.strand == -1:
                        output += '\tStrand: Antisense(-)\n'
                    if str(feat.id) != '<unknown id>':
                        output += '\tID: {0}\n'.format(feat.id)
                    if len(feat.qualifiers) > 0:
                        output += '\tQualifiers:\n'
                        qualifs = OrderedDict(sorted(feat.qualifiers.items()))
                        for key, qual in qualifs.items():
                            output += '\t\t{0}: {1}\n'.format(key, str(qual).strip("[]'"))
                    if feat.ref is not None:
                        output += '\ref: {0}'.format(feat.ref)
            else:
                output += 'None\n'

            output = "%s\n\n" % output.strip()
        _stdout(output)
        _exit("list_features")

    # List identifiers
    if in_args.list_ids:
        columns = 1 if not in_args.list_ids[0] else abs(in_args.list_ids[0])
        output = ""
        for indx, rec in enumerate(seqbuddy.records):
            output += "%s\t" % rec.id
            if (indx + 1) % columns == 0:
                output = "%s\n" % output.strip()
        _stdout("%s\n" % output.strip())
        _exit("list_ids")

    # Lowercase
    if in_args.lowercase:
        _print_recs(lowercase(seqbuddy))
        _exit("lowercase")

    # Make unique IDs
    if in_args.make_ids_unique:
        _print_recs(make_ids_unique(seqbuddy))
        _exit("make_ids_unique")

    # Map features from cDNA over to protein
    if in_args.map_features_nucl2prot:
        if len(in_args.sequence) < 2:
            _raise_error(ValueError("You must provide one DNA file and one protein file"), "map_features_nucl2prot")
        file1, file2 = in_args.sequence[:2]
        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)
        if file1.alpha == file2.alpha or (file1.alpha != IUPAC.protein and file2.alpha != IUPAC.protein):
            _raise_error(ValueError("You must provide one DNA file and one protein file"), "map_features_nucl2prot")
        if file1.alpha == IUPAC.protein:
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

        if in_args.out_format:
            seqbuddy.out_format = in_args.out_format
        _print_recs(seqbuddy)
        _exit("map_features_nucl2prot")

    # Map features from protein over to cDNA
    if in_args.map_features_prot2nucl:
        if len(in_args.sequence) < 2:
            _raise_error(ValueError("You must provide one DNA file and one protein file"), "map_features_nucl2prot")
        file1, file2 = in_args.sequence[:2]
        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)
        if file1.alpha == file2.alpha or (file1.alpha != IUPAC.protein and file2.alpha != IUPAC.protein):
            _raise_error(ValueError("You must provide one DNA file and one protein file"), "map_features_nucl2prot")
        if file1.alpha != IUPAC.protein:
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

        if in_args.out_format:
            seqbuddy.out_format = in_args.out_format
        _print_recs(seqbuddy)
        _exit("map_features_prot2nucl")

    # Merge together multiple files into a single file
    if in_args.merge:
        seqbuddy_objs = [SeqBuddy(x) for x in in_args.sequence]
        try:
            _print_recs(merge(*seqbuddy_objs))
        except RuntimeError as e:
            _raise_error(e, "merge")
        _exit("merge")

    # Molecular Weight
    if in_args.molecular_weight:
        molecular_weight(seqbuddy)
        mws = seqbuddy.molecular_weights
        if seqbuddy.alpha == (IUPAC.ambiguous_dna or IUPAC.unambiguous_dna):
            _stderr("ID\tssDNA\tdsDNA\n")
        elif seqbuddy.alpha == (IUPAC.ambiguous_rna or IUPAC.unambiguous_rna):
            _stderr("ID\tssRNA\n")
        else:
            _stderr("ID\tProtein\n")
        for indx, value in enumerate(mws['ids']):
            if len(mws['masses_ds']) != 0:
                print("{0}\t{1}\t{2}".format(value, mws['masses_ss'][indx], mws['masses_ds'][indx]))
            else:
                print("{0}\t{1}".format(value, mws['masses_ss'][indx]))
        _exit("molecular_weight")

    # Number of sequences
    if in_args.num_seqs:
        _stdout("%s\n" % num_seqs(seqbuddy))
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
        _print_recs(order_ids_randomly(seqbuddy))
        _exit("order_ids_randomly")

    # Pull random records
    if in_args.pull_random_record:
        count = 1 if not in_args.pull_random_record[0] else in_args.pull_random_record[0]
        _print_recs(pull_random_recs(seqbuddy, count))
        _exit("pull_random_record")

    # Pull record ends
    if in_args.pull_record_ends:
        _print_recs(pull_record_ends(seqbuddy, in_args.pull_record_ends))
        _exit("pull_record_ends")

    # Pull records
    if in_args.pull_records:
        _print_recs(pull_recs(seqbuddy, in_args.pull_records))
        _exit("pull_records")

    # Purge
    if in_args.purge:
        purge(seqbuddy, in_args.purge)
        record_map = "### Deleted record mapping ###\n"
        for rec in seqbuddy.records:
            record_map += "%s\n" % rec.id
            if rec.buddy_data["purge_set"]:
                for del_seq_id in rec.buddy_data["purge_set"]:
                    record_map += "%s, " % del_seq_id
            record_map = record_map.strip(", ") + "\n\n"
        record_map = record_map.strip() + "\n##############################\n\n"

        _stderr(record_map, in_args.quiet)
        _print_recs(seqbuddy)
        _exit("purge")

    # Renaming
    if in_args.rename_ids:
        args = in_args.rename_ids[0]
        if 2 > len(args) or len(args) > 3:
            _raise_error(AttributeError("rename_ids requires two or three argments: "
                                        "query, replacement, [max replacements]"), "rename_ids")
        num = 0
        try:
            num = num if len(args) == 2 else int(args[2])
        except ValueError:
            _raise_error(ValueError("Max replacements argument must be an integer"), "rename_ids")

        try:
            _print_recs(rename(seqbuddy, args[0], args[1], num))
        except AttributeError as e:
            _raise_error(e, "rename_ids", "There are more replacement")
        _exit("rename_ids")

    # Replace sub-sequences
    if in_args.replace_subseq:
        args = in_args.replace_subseq[0]
        args = args[:2] if len(args) > 1 else args
        _print_recs(replace_subsequence(seqbuddy, *args))

    # Reverse complement
    if in_args.reverse_complement:
        try:
            _print_recs(reverse_complement(seqbuddy))
        except TypeError as e:
            _raise_error(e, "reverse_complement", "Nucleic acid sequences required.")
        _exit("reverse_complement")

    # Screw formats
    if in_args.screw_formats:
        if in_args.screw_formats.lower() not in OUTPUT_FORMATS:
            _raise_error(IOError("Error: unknown format '%s'\n" % in_args.screw_formats), "screw_formats")

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
        try:
            _print_recs(select_frame(seqbuddy, in_args.select_frame))
        except TypeError as e:
            _raise_error(e, "reverse_complement", "Select frame requires nucleic acid, not protein.")

        _exit("select_frame")

    # Shuffle Seqs
    if in_args.shuffle_seqs:
        _print_recs(shuffle_seqs(seqbuddy))
        _exit("shuffle_seqs")

    # Transcribe
    if in_args.transcribe:
        if seqbuddy.alpha != IUPAC.ambiguous_dna:
            _raise_error(ValueError("You need to provide a DNA sequence."), "transcribe")
        _print_recs(dna2rna(seqbuddy))
        _exit("transcribe")

    # Translate CDS
    if in_args.translate:
        if seqbuddy.alpha == IUPAC.protein:
            _raise_error(TypeError("Nucleic acid sequence required, not protein."), "translate")
        try:
            _print_recs(translate_cds(seqbuddy, quiet=in_args.quiet))
        except TypeError as e:
            _raise_error(e, "translate", ["Nucleic acid sequence required, not protein.", " is protein."])
        _exit("translate")

    # Translate 6 reading frames
    if in_args.translate6frames:
        if seqbuddy.alpha == IUPAC.protein:
            _raise_error(TypeError("You need to supply DNA or RNA sequences to translate"), "translate6frames")
        try:
            seqbuddy = translate6frames(seqbuddy)
        except TypeError as e:
            _raise_error(e, "translate6frames", ["Nucleic acid sequence required, not protein.", " is protein."])
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
