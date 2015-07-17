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

# Standard library imports
# import pdb
# import timeit
import sys
import requests
from io import StringIO
import re

# Third party package imports
sys.path.insert(0, "./")  # For stand alone executable, where dependencies are packaged with BuddySuite
from bs4 import BeautifulSoup

# My functions
from MyFuncs import *


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
    def __init__(self, _input, _database=None, _out_format="summary"):
        self.search_terms = []
        self.accessions = {}
        self.records = {}
        self.out_format = _out_format.lower()

        # DbBuddy objects
        if type(_input) == list:
            for _dbbuddy in _input:
                if type(_dbbuddy) != DbBuddy:
                    raise TypeError("List of non-DbBuddy objects passed into DbBuddy as _input. %s" % _dbbuddy)

                self.search_terms += _dbbuddy.search_terms
                # ToDo: update() will overwrite any common records between the two dicts, should check whether values are already set first
                self.accessions.update(_dbbuddy.accessions)
                self.records.update(_dbbuddy.records)

            _input = None

        # Handles
        elif str(type(_input)) == "<class '_io.TextIOWrapper'>":
            # This will also deal with input streams (e.g., stdout pipes)
            _input = _input.read().strip()

        # Plain text
        elif type(_input) == str and not os.path.isfile(_input):
            _input = _input.strip()

        # File paths
        elif os.path.isfile(_input):
            with open(_input, "r") as _ifile:
                _input = _ifile.read().strip()

        else:
            raise GuessError("Could not determine the input type.")

        if _input:
            # try to glean accessions first
            accessions_check = re.sub("[\n\r, ]+", "\t", _input)
            accessions_check = accessions_check.split("\t")
            for _accession in accessions_check:
                db = guess_database(_accession)
                if db:
                    self.accessions[_accession] = db if not _database else _database

            # If accessions not identified, assume search terms
            if len(self.accessions) != len(accessions_check):
                search_term_check = re.sub("[\n\r,]+", "\t", _input)
                search_term_check =  search_term_check.split("\t")
                if search_term_check not in self.accessions:
                    self.search_terms.append(search_term_check)

    def __str__(self):
        _output = ""
        if self.out_format == "summary":
            for _accession in self.accessions:
                _output += "%s\t" % _accession
                try:
                    _output += "%s\n" % self.type
                except IndexError:
                    _output += "Unknown\n"
            _output = "%s/n" % _output.strip()
        return _output


class record:
    def __init__(self, _accession, _database=None, _seq=None, _meta=None, _type=None):
        self.accession = _accession
        self.database = _database
        self.seq = _seq
        self.metadata = _meta
        self.type = _type


def guess_database(accession):
    # RefSeq
    if re.match("^[NX][MR]_[0-9]+", accession):
        return "refseq_nucl"
    if re.match("^[ANYXZ]P_[0-9]+", accession):
        return "refseq_prot"

    # UniProt
    if re.match("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", accession):
        return "uniprot"

    # Ensembl stable ids (http://www.ensembl.org/info/genome/stable_ids/index.html)
    if re.match("^(ENS|FB)[A-Z]*[0-9]+", accession):
        return "ensembl"

    return None
    #server = "http://rest.ensembl.org"
    #ext = "/archive/id/ENSG00000157764?"

    #r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    #if not r.ok:
    #  r.raise_for_status()
    #  sys.exit()

    #decoded = r.json()
    #print(repr(decoded))
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

    parser.add_argument('-de', '--download_everything', action="store_true",
                        help="Retrieve full records for all search terms and accessions")
    parser.add_argument('-gd', '--guess_database', action="store_true",
                        help="List the database that each provided accession belongs to.")

    parser.add_argument('-q', '--quiet', help="Suppress stderr messages", action='store_true')
    parser.add_argument('-t', '--test', action='store_true',
                        help="Run the function and return any stderr/stdout other than sequences.")
    parser.add_argument('-o', '--out_format', default="Summary",
                        help="If you want a specific format output", action='store')
    parser.add_argument('-f', '--in_format', help="If DbBuddy can't guess the file format, just specify it directly.",
                        action='store')
    parser.add_argument('-d', '--database', choices=["all", "ensembl", "genbank", "uniprot", "dna", "protein"],
                        help='Specify a specific database or database class to search', action='store')
    in_args = parser.parse_args()

    dbbuddy = []
    search_set = ""

    try:
        if len(in_args.user_input) > 1:
            for search_set in in_args.user_input:
                dbbuddy.append(DbBuddy(search_set, in_args.database, in_args.out_format))

            dbbuddy = DbBuddy(dbbuddy, in_args.database, in_args.out_format)
        else:
            dbbuddy = DbBuddy(in_args.user_input[0], in_args.database, in_args.out_format)

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

    # Guess database
    if in_args.guess_database:
        for accession, database in dbbuddy.accessions.items():
            _stdout("%s:\t%s\n" % (accession, database))