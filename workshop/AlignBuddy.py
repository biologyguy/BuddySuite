#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Created on: Nov 18 2014 

"""
DESCRIPTION OF PROGRAM
AlignmentBuddy is a general wrapper for popular DNA and protein alignment programs that handles format conversion
and allows maintencance of rich feature annotation following alignment.
"""

import argparse
from Bio import SeqIO, AlignIO
import os
import sys
import shutil


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
