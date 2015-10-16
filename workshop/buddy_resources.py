#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, version 2 of the License (GPLv2).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details at http://www.gnu.org/licenses/.

name: shared_resources.py
date: Aug-21-2015
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: Gnu General Public License, Version 2.0 (http://www.gnu.org/licenses/gpl.html)
derivative work: No

Description:
Collection of resources used by all BuddySuite tools
Including dictionaries of the commands available for each Buddy Tool
"""
import sys
import argparse
import datetime
from collections import OrderedDict
import os
from configparser import ConfigParser
import json
import traceback
import re

sys.path.insert(0, "./")
from MyFuncs import TempFile

if __name__ == '__main__':
    sys.exit(datetime.datetime.strptime(str(datetime.date.today()), '%Y-%m-%d'))


class Contributor:
    def __init__(self, first, last, commits=None, github=None):
        self.first = first
        self.last = last
        self.commits = commits
        self.github = github

    def __str__(self):
        if not self.commits:
            commits = "0 commits"
        elif self.commits == 1:
            commits = "1 commit"
        else:
            commits = "%s commits" % self.commits

        _output = "%s %s, %s" % (self.first, self.last, commits)
        return _output if not self.github else "%s, %s" % (_output, self.github)


# Credit to rr- (http://stackoverflow.com/users/2016221/rr)
# http://stackoverflow.com/questions/18275023/dont-show-long-options-twice-in-print-help-from-argparse
class CustomHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string


class Usage:
    def __init__(self):
        self.config = config_values()
        if self.config["install_path"]:
            self.usage_file_path = "%s/usage.json" % self.config["install_path"]
        else:
            self.usage_file_path = "/tmp/usage.json"

        if not os.path.isfile(self.usage_file_path):
            with open(self.usage_file_path, "w") as ofile:
                ofile.write("{}")
        try:
            with open(self.usage_file_path) as ifile:
                self.stats = json.load(ifile)
                if not self.stats:  # Empty file needs to be populated with a little info
                    self.clear_stats()

        except ValueError:  # If the json file can't be read for whatever reason, start from scratch
            print("Error reading usage json file. Starting from scratch.")
            self.clear_stats()

    def clear_stats(self):
        self.stats = {"user_hash": self.config["user_hash"]}

    def increment(self, buddy, version, tool):
        self.stats.setdefault(buddy, {})
        self.stats[buddy].setdefault(version, {})
        self.stats[buddy][version].setdefault(tool, 0)
        self.stats[buddy][version][tool] += 1
        return

    def save(self, send_report=True):
        if self.config["diagnostics"] == "True" and send_report:
            self.stats.setdefault("last_upload", datetime.date.today().isoformat())
            if (datetime.datetime.today() - datetime.datetime.strptime(self.stats["last_upload"],
                                                                       '%Y-%m-%d')).days >= 1:
                self.send_report()
                return
        with open(self.usage_file_path, "w") as ofile:
            json.dump(self.stats, ofile)
        return

    def send_report(self):
        self.stats["date"] = str(datetime.date.today())
        from ftplib import FTP, all_errors
        temp_file = TempFile()
        json.dump(self.stats, temp_file.get_handle())
        try:
            ftp = FTP("rf-cloning.org", user="buddysuite", passwd="seqbuddy")
            ftp.storlines("STOR usage_%s" % temp_file.name, temp_file.get_handle("rb"))
            self.clear_stats()
            self.stats["last_upload"] = str(datetime.date.today())
        except all_errors as e:
            print("FTP Error: %s" % e)

        self.save(send_report=False)
        return


class Version:
    def __init__(self, name, major, minor, _contributors, release_date=None):
        self.name = name
        self.major = major
        self.minor = minor
        self.contributors = _contributors  # This needs to be a list of Contributor objects
        if not release_date:
            self.release_date = datetime.date.today()
        else:
            # e.g., release_date = {"year": 2015, "month": 3, "day": 21}
            self.release_date = datetime.datetime(**release_date)

    def contributors_string(self):
        _contributors = sorted(self.contributors, key=lambda x: x.commits, reverse=True)
        _output = ""
        for contributor in _contributors:
            _output += "%s\n" % contributor
        return _output.strip()

    def short(self):
        return "%s.%s" % (self.major, self.minor)

    def __str__(self):
        _output = '''\
%s %s.%s (%s)

Gnu General Public License, Version 2.0 (http://www.gnu.org/licenses/gpl.html)
This is free software; see the source for detailed copying conditions.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
Questions/comments/concerns can be directed to Steve Bond, steve.bond@nih.gov

Contributors:
%s
''' % (self.name, self.major, self.minor, self.release_date, self.contributors_string())
        return _output


def config_values():
    config_file = "%s/.buddysuite/config.ini" % os.path.expanduser('~')
    if os.path.isfile(config_file):
        config = ConfigParser()
        config.read(config_file)
        options = {"install_path": config.get('Install_path', 'path'),
                   "email": config.get('other', 'email'),
                   "diagnostics": config.get('other', 'diagnostics'),
                   "user_hash": config.get('other', 'user_hash')}
    else:
        options = {"install_path": False,
                   "email": "buddysuite@nih.gov",
                   "diagnostics": False,
                   "user_hash": "hashless"}
    return options


# Might want to include date in error file name
def error_report(error_msg):
    from ftplib import FTP, all_errors
    temp_file = TempFile()
    temp_file.write(error_msg)
    try:
        ftp = FTP("rf-cloning.org", user="buddysuite", passwd="seqbuddy")
        ftp.storlines("STOR error_%s" % temp_file.name, open(temp_file.path, "rb"))
    except all_errors as e:
        print("FTP Error: %s" % e)


def flags(parser, _positional=None, _flags=None, _modifiers=None, version=None):
    """
    :param parser: argparse.ArgumentParser object
    :param _positional: tuple e.g., ("user_input", "Specify accession numbers or search terms...")
    :param _flags: dict e.g., db_flags
    :param version: Version object
    :return:
    """
    if _positional:
        positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")
        positional.add_argument(_positional[0], help=_positional[1], nargs="*", default=[sys.stdin])

    if _flags:
        _flags = OrderedDict(sorted(_flags.items(), key=lambda x: x[0]))
        parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")
        for func, in_args in _flags.items():
            args = ("-%s" % in_args["flag"], "--%s" % func)
            kwargs = {}
            for cmd, val in in_args.items():
                if cmd == 'flag':
                    continue
                kwargs[cmd] = val
            parser_flags.add_argument(*args, **kwargs)

    if _modifiers:
        _modifiers = OrderedDict(sorted(_modifiers.items(), key=lambda x: x[0]))
        parser_modifiers = parser.add_argument_group(title="\033[1mModifying options\033[m")
        for func, in_args in _modifiers.items():
            args = ("-%s" % in_args["flag"], "--%s" % func)
            kwargs = {}
            for cmd, val in in_args.items():
                if cmd == 'flag':
                    continue
                kwargs[cmd] = val
            parser_modifiers.add_argument(*args, **kwargs)

    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-h', '--help', action="help", help="show this help message and exit")
    if version:
        misc.add_argument('-v', '--version', action='version', version=str(version))


def send_traceback(tool, e):
    config = config_values()
    tb = "%s\n" % config["user_hash"]
    for _line in traceback.format_tb(sys.exc_info()[2]):
        _line = re.sub('"/.*/(.*)?"', r'"\1"', _line)
        tb += _line
    tb = "%s: %s\n\n%s" % (type(e).__name__, e, tb)
    print("\033[m%s has crashed with the following traceback:\033[91m\n\n%s\n\n\033[m" % (tool, tb))

    send_diagnostic = True if config["diagnostics"] == "True" else False
    if not send_diagnostic:
        prompt = input("\033[1mWould you like to send a crash report with the above "
                       "traceback to the developers ([y]/n)?\033[m")
        if prompt.lower() in ["y", "yes", ""]:
            send_diagnostic = True

    else:
        print("An error report with the above traceback is being sent to the BuddySuite developers because "
              "you have elected to participate in the Software Improvement Program. To opt-out of this "
              "program in the future, re-run the BuddySuite installer and un-check the box on the "
              "'Diagnostics' screen.\n")

    if send_diagnostic:
        print("Preparing error report for FTP upload...\nSending...\n")
        error_report(tb)
        print("Success, thank you.\n")


contributors = [Contributor("Stephen", "Bond", 291, "https://github.com/biologyguy"),
                Contributor("Karl", "Keat", 265, "https://github.com/KarlKeat")]

# VERSIONS
VERSIONS = {"SeqBuddy": Version("SeqBuddy", 2, 'alpha', contributors),
            "DatabaseBuddy": Version("DatabaseBuddy", 1, 'alpha', contributors),
            "AlignBuddy": Version("AlignBuddy", 1, 'alpha', contributors),
            "PhyloBuddy": Version("PhyloBuddy", 1, 'alpha', contributors)}

# flag, action, nargs, metavar, help, choices, type
# #################################################### INSTALLER ##################################################### #
bsi_flags = {"cmd_line": {"flag": "cmd",
                          "action": "store_true",
                          "help": "Command line version of the installer (for non-graphical systems)."}}

bsi_modifiers = {}
# ##################################################### SEQBUDDY ##################################################### #
sb_flags = {"annotate": {"flag": "ano",
                         "nargs": "*",
                         "metavar": "args",
                         "help": "Add a feature (annotation) to selected sequences. "
                                 "Args: <name>, <location (start1-end1,start2-end2...)>, [strand (+|-)], "
                                 "[qualifier (foo=bar) [qualifier]], [regex_pattern [regex_[pattern]]}"},
            "ave_seq_length": {"flag": "asl",
                               "action": "append",
                               "nargs": "?",
                               "metavar": "'clean'",
                               "help": "Calculate average sequence length. Specify 'clean' to remove gaps etc first"},
            "back_transcribe": {"flag": "r2d",
                                "action": "store_true",
                                "help": "Convert RNA sequences to DNA"},
            "back_translate": {"flag": "btr",
                               "action": "append",
                               "nargs": "*",
                               "metavar": 'arg',
                               "help": "Convert amino acid sequences into codons. Optionally, "
                                       "select mode by passing in [{random, r, optimized, o}] "
                                       "[{human, h, mouse, m, yeast, y, ecoli, e}]"},
            "bl2seq": {"flag": "bl2s",
                       "action": "store_true",
                       "help": "All-by-all blast among sequences using bl2seq. "
                               "Only Returns the top hit from each search"},
            "blast": {"flag": "bl",
                      "action": "store",
                      "metavar": "<path/to/blast_db>",
                      "help": "Search a BLAST database with your sequence file using common settings, "
                              "returning the hits"},
            "clean_seq": {"flag": "cs",
                          "action": "append",
                          "nargs": "?",
                          "metavar": "'strict'",
                          "help": "Strip out non-sequence characters, such as stops (*) and gaps (-). "
                                  "Pass in the word 'strict' to remove all characters except the "
                                  "unambiguous letter codes"},
            "complement": {"flag": "cmp",
                           "action": "store_true",
                           "help": "Return complement of nucleotide sequence"},
            "concat_seqs": {"flag": "cts",
                            "action": "append",
                            "nargs": "?",
                            "metavar": "'clean'",
                            "help": "Concatenate a bunch of sequences into a single solid string. Pass in "
                                    "the word 'clean' to remove stops, gaps, etc., from the sequences "
                                    "before concatenating."},
            "count_codons": {"flag": "cc",
                             "action": "append",
                             "nargs": "?",
                             "metavar": "'concatenate'",
                             "help": "Return codon frequency statistics"},
            "count_residues": {"flag": "cr",
                               "action": "append",
                               "nargs": "?",
                               "metavar": "'concatenate'",
                               "help": "Generate a table of sequence compositions"},
            "delete_features": {"flag": "df",
                                "action": "store",
                                "nargs": "+",
                                "metavar": "<regex (str)>",
                                "help": "Remove specified features from all records"},
            "delete_large": {"flag": "dlg",
                             "action": "store",
                             "metavar": "<threshold (int)>",
                             "type": int,
                             "help": "Delete sequences with length above threshold"},
            "delete_metadata": {"flag": "dm",
                                "action": "store_true",
                                "help": "Remove meta-data from file (only id is retained)"},
            "delete_records": {"flag": "dr",
                               "action": "store",
                               "nargs": "+",
                               "metavar": "args",
                               "help": "Remove records from a file (deleted IDs are sent to stderr). "
                                       "Regular expressions are understood, and an int as the final argument will"
                                       "specify number of columns for deleted IDS"},
            "delete_repeats": {"flag": "drp",
                               "action": "append",
                               "nargs": "*",
                               "metavar": "columns (int)",
                               "help": "Strip repeat records (ids and/or identical sequences. "
                                       "Optional, pass in an integer to specify number of columns for deleted IDs"},
            "delete_small": {"flag": "dsm",
                             "action": "store",
                             "metavar": "<threshold (int)>",
                             "type": int,
                             "help": "Delete sequences with length below threshold"},
            "extract_region": {"flag": "er",
                               "action": "store",
                               "nargs": 2,
                               "metavar": ("<start (int)>", "<end (int)>"),
                               "type": int,
                               "help": "Pull out sub-sequences"},
            "find_CpG": {"flag": "fcpg",
                         "action": "store_true",
                         "help": "Predict regions under strong purifying selection based on high CpG content"},
            "find_pattern": {"flag": "fp",
                             "action": "store",
                             "nargs": "+",
                             "metavar": "<regex>",
                             "help": "Search for subsequences, returning the start positions of all matches"},
            "find_repeats": {"flag": "frp",
                             "action": "append",
                             "nargs": "?",
                             "type": int,
                             "metavar": "columns (int)",
                             "help": "Identify whether a file contains repeat sequences "
                                     "and/or sequence ids. The number of output columns "
                                     "can be modified by passing in an integer."},
            "find_restriction_sites": {"flag": "frs",
                                       "action": "append",
                                       "nargs": "*",
                                       "metavar": "",
                                       "help": "Identify restriction sites. Args: [enzymes "
                                               "{specific enzymes, commercial, all}], [Num cuts (int) [num cuts]], "
                                               "[order {alpha, position}]"},
            "guess_alphabet": {"flag": "ga",
                               "action": "store_true",
                               "help": "Glean the alphabet type of input file(s)"},
            "guess_format": {"flag": "gf",
                             "action": "store_true",
                             "help": "Glean the flat file format of input file(s)"},
            "hash_seq_ids": {"flag": "hsi",
                             "action": "append",
                             "nargs": "?",
                             "type": int,
                             "metavar": "hash length (int)",
                             "help": "Rename all sequence IDs to fixed length hashes. "
                                     "Default length is 10."},
            "insert_seq": {"flag": "is",
                           "action": "append",
                           "nargs": "*",
                           "metavar": ("<sequence>", "<front|rear|index(int)>"),
                           "help": "Insert a sequence at the desired location"},
            "isoelectric_point": {"flag": "ip",
                                  "action": "store_true",
                                  "help": "Calculate isoelectric points"},
            "list_features": {"flag": "lf",
                              "action": "store_true",
                              "help": "Print out all sequence annotations"},
            "list_ids": {"flag": "li",
                         "action": "append",
                         "nargs": "?",
                         "type": int,
                         "metavar": "Columns (int)",
                         "help": "Output the sequence identifiers. Optionally, pass in an integer to "
                                 "specify the # of columns to write"},
            "lowercase": {"flag": "lc",
                          "action": "store_true",
                          "help": "Convert all sequences to lowercase"},
            "make_unique_ids": {"flag": "mui",
                                "action": "store_true",
                                "help": "Add a number at the end of replicate ids to make them unique"},
            "map_features_nucl2prot": {"flag": "fn2p",
                                       "action": "store_true",
                                       "help": "Take the features annotated onto nucleotide sequences "
                                               "and map to protein sequences. Both a protein and "
                                               "cDNA file must be passed in"},
            "map_features_prot2nucl": {"flag": "fp2n",
                                       "action": "store_true",
                                       "help": "Take the features annotated onto protein sequences "
                                               "and map to cDNA sequences. Both a protein and "
                                               "cDNA file must be passed in"},
            "merge": {"flag": "mrg",
                      "action": "store_true",
                      "help": "Merge multiple copies of sequence records together, "
                              "combining their feature lists"},
            "molecular_weight": {"flag": "mw",
                                 "action": "store_true",
                                 "help": "Compute the molecular weight of sequences"},
            "num_seqs": {"flag": "ns",
                         "action": "store_true",
                         "help": "Counts how many sequences are present in an input file"},
            "order_features_alphabetically": {"flag": "ofa",
                                              "action": "append",
                                              "nargs": "?",
                                              "metavar": "'rev'",
                                              "help": "Change the output order of sequence features, based "
                                                      "on feature name. Pass in 'rev' to reverse order"},
            "order_features_by_position": {"flag": "ofp",
                                           "action": "append",
                                           "nargs": "?",
                                           "metavar": "'rev'",
                                           "help": "Change the output order of sequence features, based on "
                                                   "sequence position. Pass in 'rev' to reverse order"},
            "order_ids": {"flag": "oi",
                          "action": "append",
                          "nargs": "?",
                          "metavar": "'rev'",
                          "help": "Sort sequences by id alpha-numerically. Pass in the word 'rev' to reverse order"},
            "order_ids_randomly": {"flag": "oir",
                                   "action": "store_true",
                                   "help": "Randomly reorder the position of each record"},
            "pull_random_record": {"flag": "prr",
                                   "action": "append",
                                   "nargs": "?",
                                   "type": int,
                                   "metavar": "Num records (int)",
                                   "help": "Extract random sequences. Optionally, pass in an integer to "
                                           "increase the number of sequences returned"},
            "pull_record_ends": {"flag": "pre",
                                 "action": "store",
                                 "type": int,
                                 "metavar": "<amount (int)>",
                                 "help": "Get the ends of all sequences in a file (use negative numbers to get rear)"},
            "pull_records": {"flag": "pr",
                             "action": "store",
                             "nargs": "+",
                             "metavar": "<regex>",
                             "help": "Get all the records with ids containing a given string"},
            "purge": {"flag": "prg",
                      "action": "store",
                      "metavar": "<Max BLAST score (int)>",
                      "type": int,
                      "help": "Delete sequences with high similarity"},
            "raw_seq": {"flag": "rs",
                        "action": "store_true",
                        "help": "Return line break separated sequences"},
            "rename_ids": {"flag": "ri",
                           "action": "store",
                           "metavar": ("<pattern>", "<substitution>"),
                           "nargs": 2,
                           "help": "Replace some pattern in ids with something else. "
                                   "Limit number of replacements with -p"},
            "reverse_complement": {"flag": "rc",
                                   "action": "store_true",
                                   "help": "Return reverse complement of nucleotide sequence"},
            "screw_formats": {"flag": "sf",
                              "action": "store",
                              "metavar": "<out_format>",
                              "help": "Change the file format to something else"},
            "select_frame": {"flag": "sfr",
                             "action": "store",
                             "metavar": "<frame (int)>",
                             "type": int,
                             "choices": [1, 2, 3],
                             "help": "Change the reading from of sequences by deleting characters off of the front"},
            "shuffle_seqs": {"flag": "ss",
                             "action": "store_true",
                             "help": "Shuffles the letters in all the sequences"},
            "split_by_taxa": {"flag": "sbt",
                              "action": "store",
                              "nargs": 2,
                              "metavar": ("<Split Pattern>", "<out dir>"),
                              "help": "Sort sequences into separate files based on taxa"},
            "split_to_files": {"flag": "stf",
                               "action": "store",
                               "metavar": "<out dir>",
                               "help": "Write individual files for each sequence"},
            "transcribe": {"flag": "d2r",
                           "action": "store_true",
                           "help": "Convert DNA sequences to RNA"},
            "translate": {"flag": "tr",
                          "action": "store_true",
                          "help": "Convert coding sequences into amino acid sequences"},
            "translate6frames": {"flag": "tr6",
                                 "action": "store_true",
                                 "help": "Translate nucleotide sequences into all six reading frames"},
            "uppercase": {"flag": "uc",
                          "action": "store_true",
                          "help": "Convert all sequences to uppercase"}}

sb_modifiers = {"alpha": {"flag": "a",
                          "action": "store",
                          "help": "If you want the file read with a specific alphabet"},
                "in_format": {"flag": "f",
                              "action": "store",
                              "help": "If SeqBuddy can't guess the file format, just specify it directly"},
                "in_place": {"flag": "i",
                             "action": "store_true",
                             "help": "Rewrite the input file in-place. Be careful!"},
                "out_format": {"flag": "o",
                               "action": "store",
                               "help": "If you want a specific format output"},
                "params": {"flag": "p",
                           "action": "store",
                           "nargs": "+",
                           "help": "Free form arguments for some functions"},
                "quiet": {"flag": "q",
                          "action": "store_true",
                          "help": "Suppress stderr messages"},
                "test": {"flag": "t",
                         "action": "store_true",
                         "help": "Run the function and return any stderr/stdout other than sequences"}}

# #################################################### ALIGNBUDDY #################################################### #
alb_flags = {"alignment_lengths": {"flag": "al",
                                   "action": "store_true",
                                   "help": "Returns a list of alignment lengths"},
             "back_transcribe": {"flag": "r2d",
                                 "action": "store_true",
                                 "help": "Convert RNA alignments to DNA"},
             "clean_seq": {"flag": "cs",
                           "action": "append",
                           "nargs": "?",
                           "help": "Strip out non-sequence characters, such as stops (*) "
                                   "and gaps (-). Pass in the word 'strict' to remove all "
                                   "characters except the unambiguous letter codes"},
             "codon_alignment": {"flag": "ca",
                                 "action": "store_true",
                                 "help": "Shift all gaps so the sequence is in triplets"},
             "concat_alignments": {"flag": "cta",
                                   "action": "store",
                                   "help": "Concatenates two or more alignments by splitting and matching the "
                                           "sequence identifiers. Arguments: <split_pattern>"},
             "delete_rows": {"flag": "dr",
                             "action": "store",
                             "help": "Remove selected rows from alignments. Arguments: <search_pattern>"},
             "extract_range": {"flag": "er",
                               "action": "store",
                               "nargs": 2,
                               "metavar": ("<start (int)>", "<end (int)>"),
                               "type": int,
                               "help": "Pull out sub-alignments in a given range"},
             "generate_alignment": {"flag": "ga",
                                    "action": "append",
                                    "help": ""},
             "list_ids": {"flag": "li",
                          "action": "append",
                          "nargs": "?",
                          "type": int,
                          "metavar": "int (optional)",
                          "help": "Output all the sequence identifiers in a file. Optionally, pass in an integer to "
                                  "specify the # of columns to write"},
             "lowercase": {"flag": "lc",
                           "action": "store_true",
                           "help": "Convert all sequences to lowercase"},
             "num_seqs": {"flag": "ns",
                          "action": "store_true",
                          "help": "Counts how many sequences are present in each alignment"},
             "order_ids": {"flag": "oi",
                           "action": "append",
                           "nargs": "?",
                           "help": "Sort all sequences in an alignment by id in alpha-numeric order. "
                                   "Pass in the word 'rev' to reverse order"},
             "pull_rows": {"flag": "pr",
                           "action": "store",
                           "help": "Keep selected rows from alignements. Arguments: <search_pattern>"},
             "rename_ids": {"flag": "ri",
                            "action": "store",
                            "nargs": 2,
                            "metavar": ("<pattern>", "<substitution>"),
                            "help": "Replace some pattern in ids with something else. "
                                    "Limit number of replacements with -p"},
             "screw_formats": {"flag": "sf",
                               "action": "store",
                               "help": "Arguments: <out_format>"},
             "split_to_files": {"flag": "stf",
                                "action": "store",
                                "nargs": 2,
                                "metavar": ("<out dir>", "<out file>"),
                                "help": "Write individual files for each alignment"},
             "translate": {"flag": "tr",
                           "action": "store_true",
                           "help": "Convert coding sequences into amino acid sequences"},
             "trimal": {"flag": "trm",
                        "action": "append",
                        "nargs": "?",
                        "help": "Delete columns with a certain percentage of gaps. Or auto-detect with 'gappyout'"},
             "transcribe": {"flag": "d2r",
                            "action": "store_true",
                            "help": "Convert DNA alignments to RNA"},
             "uppercase": {"flag": "uc",
                           "action": "store_true",
                           "help": "Convert all sequences to uppercase"},
             }

alb_modifiers = {"in_format": {"flag": "f",
                               "action": "store",
                               "help": "If AlignBuddy can't guess the file format, just specify it directly"},
                 "in_place": {"flag": "i",
                              "action": "store_true",
                              "help": "Rewrite the input file in-place. Be careful!"},
                 "keep_temp": {"flag": "k",
                               "action": "store",
                               "help": "Save temporary files created by generate_tree in current working directory"},
                 "out_format": {"flag": "o",
                                "action": "store",
                                "help": "If you want a specific format output"},
                 "params": {"flag": "p",
                            "action": "store",
                            "nargs": "+",
                            "help": "Free form arguments for some functions"},
                 "quiet": {"flag": "q",
                           "action": "store_true",
                           "help": "Suppress stderr messages"},
                 "test": {"flag": "t",
                          "action": "store_true",
                          "help": "Run the function and return any stderr/stdout other than sequences"}}

# #################################################### PHYLOBUDDY #################################################### #

pb_flags = {"consensus_tree": {"flag": "ct",
                               "action": "append",
                               "nargs": "?",
                               "type": float,
                               "metavar": "min frequency (default 0.5)",
                               "help": "Generate a consensus tree"},
            "display_trees": {"flag": "dt",
                              "action": "store_true",
                              "help": "Visualize trees graphically"},
            "distance": {"flag": "dis",
                         "nargs": "?",
                         "action": "append",
                         "choices": ["wrf", "uwrf", "ed"],
                         "help": "Calculate similarity metrics between pairs of trees"},
            "generate_tree": {"flag": "gt",
                              "nargs": "*",
                              "action": "append",
                              "type": str,
                              "metavar": ("{'raxml', 'phyml', 'fasttree'}", "'program specific arguments'"),
                              "help": "Accept alignment file as input, and perform "
                                      "phylogenetic inference with a third party program"},
            "list_ids": {"flag": "li",
                         "action": "append",
                         "nargs": "?",
                         "type": int,
                         "metavar": "Num columns",
                         "help": "Display all taxa ids"},
            "print_trees": {"flag": "ptr",
                            "action": "store_true",
                            "help": "Output trees to the terminal (ascii)"},
            "prune_taxa": {"flag": "pt",
                           "action": "append",
                           "nargs": "+",
                           "metavar": "Regex",
                           "help": "Remove taxa with matching labels/IDs"},
            "rename_ids": {"flag": "ri",
                           "action": "store",
                           "nargs": 2,
                           "metavar": ("<pattern>", "<substitution>"),
                           "help": "Replace some pattern in ids with something else"},
            "root": {"flag": "rt",
                     "action": "append",
                     "nargs": "*",
                     "metavar": "Rooting taxa",
                     "help": "(Re)root a tree at its midpoint or on specified taxa"},
            "screw_formats": {"flag": "sf",
                              "action": "store",
                              "metavar": "<out format>",
                              "help": "Change the format of a tree file"},
            "show_unique": {"flag": "su",
                            "action": "store_true",
                            "help": "Color leaf branches based on whether the taxa is present in both trees"},
            "split_polytomies": {"flag": "sp",
                                 "action": "store_true",
                                 "help": "Create a binary tree by splitting polytomies randomly"},
            "unroot": {"flag": "ur",
                               "action": "store_true",
                               "help": "Remove any roots"}
            }

pb_modifiers = {"in_format": {"flag": "f",
                              "action": "store",
                              "metavar": "<format>",
                              "help": "If PhyloBuddy can't guess the file format, try specifying it directly"},
                "in_place": {"flag": "i",
                             "action": "store_true",
                             "help": "Rewrite the input file in-place. Be careful!"},
                "keep_temp": {"flag": "k",
                              "nargs": "?",
                              "action": "store",
                              "metavar": "path",
                              "help": "Save temporary files, if any; default to current working directory"},
                "out_format": {"flag": "o",
                               "metavar": "<format>",
                               "action": "store",
                               "help": "Choose a specific output format"},
                "quiet": {"flag": "q",
                          "action": "store_true",
                          "help": "Suppress stderr messages"},
                "test": {"flag": "t",
                         "action": "store_true",
                         "help": "Run the function and return any stderr/stdout other than trees"}}

# ################################################## DATABASEBUDDY ################################################### #
db_flags = {"guess_database": {"flag": "gd",
                               "action": "store_true",
                               "help": "List the database that each provided accession belongs to."},
            "live_shell": {"flag": "ls",
                           "action": "store_true",
                           "help": "Interactive database searching. The best tool for sequence discovery."},
            # "retrieve_accessions": {"flag": "ra",
            #                        "action": "store_true",
            #                        "help": "Use search terms to find a list of sequence accession numbers"},
            # "retrieve_sequences": {"flag": "rs",
            #                       "action": "store_true",
            #                       "help": "Get sequences for every included accession"}
            }

db_modifiers = {"database": {"flag": "d",
                             "action": "store",
                             "choices": [],  # This needs to be set to DATABASES in the main program
                             "help": "Specify a specific database or database class to search"},
                "out_format": {"flag": "o",
                               "action": "store",
                               "help": "If you want a specific format output"},
                # "quiet": {"flag": "q",
                #          "action": "store_true",
                #          "help": "Suppress stderr messages"},
                # "test": {"flag": "t",
                #         "action": "store_true",
                #         "help": "Run the function and return any stderr/stdout other than sequences"}
                }
