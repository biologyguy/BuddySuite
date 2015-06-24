#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Created on: Nov 18 2014 

"""
DESCRIPTION OF PROGRAM
AlignmentBuddy is a general wrapper for popular DNA and protein alignment programs that handles format conversion
and allows maintencance of rich feature annotation following alignment.
"""

# Standard library imports
import sys
import os
import argparse
import shutil
from io import StringIO
from random import sample

# Third party package imports
from Bio import SeqIO, AlignIO


# ##################################################### WISH LIST #################################################### #
# 'Clean' an alignment, as implemented in phyutility


# ################################################# HELPER FUNCTIONS ################################################# #
def _stderr(message, quiet=False):
    if not quiet:
        sys.stderr.write(message)
    return

class GuessError(Exception):
    """Raised when input format cannot be guessed"""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

# #################################################### ALIGN BUDDY ################################################### #
class AlignBuddy:  # Open a file or read a handle and parse, or convert raw into a Seq object
    def __init__(self, _input, _in_format=None, _out_format=None):
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

        # Raw sequences
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

        if not _in_format:
            self.in_format = guess_format(_input)
            self.out_format = str(self.in_format) if not _out_format else _out_format

        else:
            self.in_format = _in_format

        if not self.in_format:
            if in_file:
                raise GuessError("Could not determine format from _input file '{0}'.\n"
                                 "Try explicitly setting with -f flag.".format(in_file))
            elif raw_seq:
                raise GuessError("Could not determine format from raw input\n{0} ..."
                                 "Try explicitly setting with -f flag.".format(raw_seq)[:50])
            elif in_handle:
                raise GuessError("Could not determine format from input file-like object\n{0} ..."
                                 "Try explicitly setting with -f flag.".format(in_handle)[:50])
            else:
                raise GuessError("Unable to determine format or input type. Please check how SeqBuddy is being called.")

        self.out_format = self.in_format if not _out_format else _out_format

        # ####  RECORDS  #### #
        if str(type(_input)) == "<class '__main__.SeqBuddy'>":
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

    def print(self):
        _output = ""
        for _rec in self.records:
            _output += _rec.format(self.out_format)
        return _output

    def write(self, _file_path):
        with open(_file_path, "w") as _ofile:
            SeqIO.write(self.records, _ofile, self.out_format)
        return


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

def screw_formats_align(_alignments, _out_format):
    _output = ""
    if _out_format == "phylipi":
        if len(_alignments) > 1:
            print("Warning: the input file contained more than one alignment, but phylip can only handle one. "
                  "The topmost alignment is shown here.", file=sys.stderr)
        _seqs = list(_alignments[0])
        _output += " %s %s\n" % (len(_seqs), len(_seqs[0].seq))
        max_id_length = 0
        for _seq in _seqs:
            max_id_length = len(_seq.id) if len(_seq.id) > max_id_length else max_id_length

        for _seq in _seqs:
            _seq_id = _seq.id.ljust(max_id_length)
            _output += "%s  %s\n" % (_seq_id, _seq.seq)
    else:
        for alignment in _alignments:
            _output += alignment.format(_out_format)

    return _output

'''
# ################################################ INTERNAL FUNCTIONS ################################################ #
def _set_alphabet(_sequences, alpha=None):  # update sequence alphabet in place
    if not alpha:
        alpha = guess_alphabet(_sequences)
    if alpha == "nucl":
        alpha = IUPAC.ambiguous_dna
    elif alpha == "prot":
        alpha = IUPAC.protein
    else:
        sys.exit("Error: Can't deterimine alphabet in _set_alphabet")
    for i in range(len(_sequences)):
        _sequences[i].seq.alphabet = alpha
    return _sequences


def sequence_list(sequence, _seq_format=None):  # Open a file and parse, or convert raw into a Seq object
    if isinstance(sequence, list):
        _sequences = sequence
    elif os.path.isfile(sequence):
        if not _seq_format:
            _seq_format = guess_format(sequence)
        if not _seq_format:
            sys.exit("Error: could not determine the format of your input sequence file. Explicitly set with -r flag.")
        with open(sequence, "r") as _infile:
            _sequences = list(SeqIO.parse(_infile, _seq_format))
    else:
        # dna_or_prot = IUPAC.protein if guess_alphabet(sequence) == "prot" else IUPAC.ambiguous_dna
        _sequences = [SeqRecord(Seq(sequence))]

    return _sequences


def _print_recs(_sequences, in_place=False):
    if len(_sequences) == 0:
        print("Nothing returned.", file=sys.stderr)
        return False
    _sequences = _set_alphabet(_sequences)
    _output = ""
    for _rec in _sequences:
        try:
            _output += _rec.format(out_format) + "\n"
        except ValueError as e:
            print("Error: %s\n" % e, file=sys.stderr)

    if in_args.in_place and in_place_allowed:  # TODO This is broken! Don't call in_args from here
        if not os.path.exists(in_args.sequence[0]):
            print("Warning: The -i flag was passed in, but the positional argument doesn't seem to be a file. Nothing "
                  "was written.",
                  file=sys.stderr)
            print(_output.strip())
        else:
            with open(os.path.abspath(in_args.sequence[0]), "w") as ofile:
                ofile.write(_output)
            print("File over-written at:\n%s" % os.path.abspath(in_args.sequence[0]), file=sys.stderr)
    else:
        print(_output.strip())

# #################################################################################################################### #
'''

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="alignBuddy", description="Sequience alignment with a splash of Kava")

    parser.add_argument("in_file", help="The file(s) you want to start working on", action="store", nargs="+")
    parser.add_argument("-al", "--align",
                        help="Pick your aligner. If the package is not in your $PATH, specify the location with -ab."
                             "Set any options (all in quotes) with the -p flag",
                        choices=["mafft", "prank", "pagan", "muscle", "clustalw"])
    parser.add_argument("-ab", "--align_binary", help="Specify the path to your alignment package if not in $PATH",
                        action="store")
    parser.add_argument('-sf', '--screw_formats', action='store', help="Arguments: <out_format>")

    parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')
    parser.add_argument('-p', '--params', help="Free form arguments for some functions", nargs="+", action='store')
    parser.add_argument('-o', '--out_format', help="Some functions use this flag for output format", action='store')
    parser.add_argument('-f', '--in_format', help="If the file extension isn't sane, specify the format", action='store')

    in_args = parser.parse_args()

    with open(in_args.in_file[0], "r") as ifile:
        alignment = list(AlignIO.parse(ifile, "nexus"))

    if in_args.screw_formats:
        new_align = screw_formats_align(alignment, in_args.screw_formats)
        print(new_align)


    if in_args.align_binary:
        if not shutil.which(in_args.alignment_package):
            sys.exit("Error: Unable to locate %s in your $PATH. Please specify a path to the binary with the -a flag."
                     % in_args.alignment_package)

    elif False:
        align_binary = os.path.abspath(in_args.align_binary)
        if not os.path.isfile(align_binary):
            sys.exit("Error: Unable to resolve the provided path to %s.\nUser input: %s" %
                     (in_args.alignment_package, align_binary))
