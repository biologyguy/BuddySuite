#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, version 2 of the License (GPLv2).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details at http://www.gnu.org/licenses/.

name: DbBuddy.py
date: July-16-2015
version: 1, unstable
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: Gnu General Public License, Version 2.0 (http://www.gnu.org/licenses/gpl.html)
derivative work: No

Description:
Collection of functions that interact with public sequence databases. Pull them into a script or run from command line.
"""

import sys


# ##################################################### WISH LIST #################################################### #
def get_genbank_file():
    x = 1
    return x


# ################################################# HELPER FUNCTIONS ################################################# #
class GuessError(Exception):
    """Raised when input format cannot be guessed"""
    def __init__(self, _value):
        self.value = _value

    def __str__(self):
        return self.value


def _stderr(message, quiet=False):
    if not quiet:
        sys.stderr.write(message)
    return


def _stdout(message, quiet=False):
    if not quiet:
        sys.stdout.write(message)
    return


# ##################################################### DB BUDDY ##################################################### #
class DbBuddy:  # Open a file or read a handle and parse, or convert raw into a Seq object
    def __init__(self, _input, _in_format=None, _out_format=None):
        # ####  IN AND OUT FORMATS  #### #
        # Holders for input type. Used for some error handling below
        in_handle = None
        _raw_search = None
        in_file = None


def guess_database(_input):  # _input can be list, SeqBuddy object, file handle, or file path.
    # If input is just a list, there is no BioPython in-format. Default to gb.
    if isinstance(_input, list):
        return "gb"
# #################################################################################################################### #


def search_ensembl(_dbbuddy):
    return _dbbuddy


def download_everything(_dbbuddy):
    return _dbbuddy

# ################################################# COMMAND LINE UI ################################################## #
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="DbBuddy.py", formatter_class=argparse.RawTextHelpFormatter,
                                     description="Commandline wrapper for all the fun functions in this file. "
                                                 "Access those databases!")

    parser.add_argument("user_input", help="Specify accession numbers or search terms, either in a file or as space "
                                           "separated list", nargs="*", default=[sys.stdin])

    parser.add_argument('-v', '--version', action='version',
                        version='''\
DbBuddy 1.alpha (2015)

Gnu General Public License, Version 2.0 (http://www.gnu.org/licenses/gpl.html)
This is free software; see the source for detailed copying conditions.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
Questions/comments/concerns can be directed to Steve Bond, steve.bond@nih.gov''')

    parser.add_argument('-de', '--download_everything',
                        help="Retrieve full records for all search terms and accessions")
    parser.add_argument('-q', '--quiet', help="Suppress stderr messages", action='store_true')
    parser.add_argument('-t', '--test', action='store_true',
                        help="Run the function and return any stderr/stdout other than sequences.")
    parser.add_argument('-o', '--out_format', help="If you want a specific format output", action='store')
    parser.add_argument('-f', '--in_format', help="If SeqBuddy can't guess the file format, just specify it directly.",
                        action='store')

    in_args = parser.parse_args()

    dbbuddy = []
    search_set = ""

    try:
        for search_set in in_args.sequence:
            search_set = DbBuddy(search_set, in_args.in_format, in_args.out_format)
            dbbuddy += search_set.records

        seqbuddy = DbBuddy(dbbuddy, search_set.in_format, search_set.out_format)
    except GuessError:
        sys.exit("Error: SeqBuddy could not understand your input. "
                 "Check the file path or try specifying an input type with -f")

    # ############################################# INTERNAL FUNCTION ################################################ #
    def _print_recs(_seqbuddy):
        if in_args.test:
            _stderr("*** Test passed ***\n", in_args.quiet)
            pass

        else:
            _stdout("{0}\n".format(str(_seqbuddy).rstrip()))

    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Download everything
    if in_args.download_everything:
        _print_recs(download_everything(dbbuddy))
