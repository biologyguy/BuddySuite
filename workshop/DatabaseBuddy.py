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
from urllib.parse import urlencode
from urllib.error import HTTPError
from urllib.request import Request, urlopen
from time import time, sleep
import json

# Third party package imports
sys.path.insert(0, "./")  # For stand alone executable, where dependencies are packaged with BuddySuite
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# My functions
from MyFuncs import *


# ##################################################### WISH LIST #################################################### #


# ################################################# HELPER FUNCTIONS ################################################# #
class GuessError(Exception):
    """Raised when input format cannot be guessed"""
    def __init__(self, _value):
        self.value = _value

    def __str__(self):
        return self.value


class DatabaseError(Exception):
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
        if _database and _database.lower() not in ["gb", "genbank", "uniprot", "ensembl"]:
            raise DatabaseError("%s is not currently supported." % _database)

        if _database and _database.lower() == "genbank":
            _database = "gb"

        self.search_terms = []
        self.records = {}
        self.out_format = _out_format.lower()
        self.failures = []

        # DbBuddy objects
        if type(_input) == list:
            for _dbbuddy in _input:
                if type(_dbbuddy) != DbBuddy:
                    raise TypeError("List of non-DbBuddy objects passed into DbBuddy as _input. %s" % _dbbuddy)

                self.search_terms += _dbbuddy.search_terms
                # ToDo: update() will overwrite any common records between the two dicts, should check whether values are already set first
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
            raise GuessError("DbBuddy could not determine the input type.")

        if _input:
            # try to glean accessions first
            accessions_check = re.sub("[\n\r, ]+", "\t", _input)
            accessions_check = accessions_check.split("\t")
            for _accession in accessions_check:
                _record = Record(_accession)
                _record.guess_database()
                if _record.database:
                    if _database:
                        _record.database = _database.lower()
                    self.records[_accession] = _record

            # If accessions not identified, assume search terms
            if len(self.records) != len(accessions_check):
                search_term_check = re.sub("[\n\r,]+", "\t", _input)
                search_term_check = [x.strip() for x in search_term_check.split("\t")]
                for search_term in search_term_check:
                    if search_term not in self.records:
                        self.search_terms.append(search_term)

    def __str__(self):
        _output = ""
        if self.out_format == "summary":
            for _accession, _rec in self.records.items():
                _output += "%s\t%s\n" % (_accession, _rec.database)
            _output = "%s\n" % _output.strip()

        else:
            nuc_recs = [_rec.record for _accession, _rec in self.records.items() if _rec.type == "nucleotide" and _rec.record]
            prot_recs = [_rec.record for _accession, _rec in self.records.items() if _rec.type == "protein" and _rec.record]
            tmp_dir = TemporaryDirectory()

            if len(nuc_recs) > 0:
                with open("%s/seqs.tmp" % tmp_dir.name, "w") as _ofile:
                    SeqIO.write(nuc_recs, _ofile, self.out_format)

                with open("%s/seqs.tmp" % tmp_dir.name, "r") as ifile:
                    _output = "%s\n" % ifile.read()

            if len(prot_recs) > 0:
                with open("%s/seqs.tmp" % tmp_dir.name, "w") as _ofile:
                    SeqIO.write(prot_recs, _ofile, self.out_format)

                with open("%s/seqs.tmp" % tmp_dir.name, "r") as ifile:
                    _output += "%s\n" % ifile.read()

        if _output == "":
            _output = "No records returned\n"

        return _output


class Record:
    def __init__(self, _accession, _database=None, _record=None, _type=None):
        self.accession = _accession
        self.database = _database  # genbank, ensembl, uniprot
        self.record = _record  # SeqIO record
        self.type = _type  # protein, nucleotide

    def guess_database(self):
        # RefSeq
        if re.match("^[NX][MR]_[0-9]+", self.accession):
            self.database = "gb"
            self.type = "nucleotide"

        if re.match("^[ANYXZ]P_[0-9]+", self.accession):
            self.database = "gb"
            self.type = "protein"

        # UniProt
        if re.match("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", self.accession):
            self.database = "uniprot"
            self.type = "protein"

        # Ensembl stable ids (http://www.ensembl.org/info/genome/stable_ids/index.html)
        if re.match("^(ENS|FB)[A-Z]*[0-9]+", self.accession):
            self.database = "ensembl"
            self.type = "nucleotide"

        return

    def __str__(self):
        return "Accession:\t{0}\nDatabase:\t{1}\nRecord:\t{2}\nType:\t{3}\n".format(self.accession, self.database,
                                                                                    self.record, self.type)


# ################################################# Database Clients ################################################# #
class UniProtRestClient:
    def __init__(self, _dbbuddy, server='http://www.uniprot.org/uniprot'):
        self.dbbuddy = _dbbuddy
        self.server = server

    def fetch_proteins(self):
        _records = [_rec for _accession, _rec in self.dbbuddy.records.items() if _rec.database == "uniprot"]
        _ids = ",".join([_rec.accession for _rec in _records]).strip(",")
        request = Request("%s?query=%s&format=txt" % (self.server, _ids))
        response = urlopen(request)
        data = SeqIO.to_dict(SeqIO.parse(response, "swiss"))
        for _rec in _records:
            if _rec.accession not in data:
                self.dbbuddy.failures.append(_rec.accession)
            else:
                self.dbbuddy.records[_rec.accession].record = data[_rec.accession]


class NCBIClient:
    def __init__(self, _dbbuddy):
        Entrez.email = ""
        self.dbbuddy = _dbbuddy

    def fetch_nucliotides(self):
        for _accession, _rec in self.dbbuddy.records.items():
            if _rec.database == "gb" and _rec.type == "nucleotide":
                try:
                    handle = Entrez.efetch(db="nucleotide", id=_rec.accession, rettype="gb", retmode="text")
                    self.dbbuddy.records[_accession].record = SeqIO.read(handle, "gb")
                except HTTPError as e:
                    if e.getcode() == 400:
                        self.dbbuddy.failures.append(_rec.accession)

    def fetch_proteins(self):
        for _accession, _rec in self.dbbuddy.records.items():
            if _rec.database == "gb" and _rec.type == "protein":
                try:
                    handle = Entrez.efetch(db="protein", id=_rec.accession, rettype="gb", retmode="text")
                    self.dbbuddy.records[_accession].record = SeqIO.read(handle, "gb")
                except HTTPError as e:
                    if e.getcode() == 400:
                        self.dbbuddy.failures.append(_rec.accession)


class EnsemblRestClient:
    def __init__(self, _dbbuddy, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.dbbuddy = _dbbuddy
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, _id, headers=None, params=None):
        endpoint = "/%s/%s" % (endpoint.strip("/"), _id)
        if not headers:
            headers = {}

        if 'Content-Type' not in headers:
            headers['Content-Type'] = 'text/x-seqxml+xml'

        if params:
            endpoint += '?' + urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time() - self.last_req
            if delta < 1:
                sleep(1 - delta)
            self.last_req = time()
            self.req_count = 0

        try:
            request = Request(self.server + endpoint, headers=headers)
            response = urlopen(request)
            if headers["Content-Type"] == "application/json":
                content = response.read().decode()
                data = json.loads(content)
            else:
                data = SeqIO.read(response, "seqxml")

            self.req_count += 1

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.getcode() == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    sleep(float(retry))
                    self.perform_rest_action(endpoint, headers, params)
            else:
                self.dbbuddy.failures.append(_id)
                #sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    def fetch_nucleotides(self):
        for _accession, _rec in self.dbbuddy.records.items():
            if _rec.database == "ensembl" and _rec.type == "nucleotide":
                _rec.record = self.perform_rest_action("/sequence/id", _rec.accession)

                if _rec.record:
                    rec_info = self.perform_rest_action("/lookup/id", _rec.accession, headers={"Content-Type": "application/json"})

                if _rec.record and rec_info:
                    _rec.record.name = rec_info["display_name"]
                    _rec.record.description = rec_info["description"]
                    _rec.record.annotations["keywords"] = ['Ensembl']
                    _rec.record.annotations["source"] = rec_info["source"]
                    _rec.record.annotations["organism"] = rec_info["species"]
                    _rec.record.annotations["accessions"] = [rec_info["id"]]
                    _rec.record.annotations["sequence_version"] = int(rec_info["version"])


# #################################################################################################################### #

def download_everything(_dbbuddy):
    # Get sequences from UniProt
    uniprot = UniProtRestClient(_dbbuddy)
    uniprot.fetch_proteins()

    # Get sequences from Ensembl
    ensembl = EnsemblRestClient(_dbbuddy)
    ensembl.fetch_nucleotides()

    # Get sequences from genbank
    refseq = NCBIClient(_dbbuddy)
    refseq.fetch_nucliotides()
    refseq.fetch_proteins()

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
    parser.add_argument('-o', '--out_format', help="If you want a specific format output", action='store')
    parser.add_argument('-f', '--in_format', action='store',
                        help="If DbBuddy can't guess the accession format, try specifying it directly.")
    parser.add_argument('-d', '--database', choices=["all", "ensembl", "genbank", "uniprot", "dna", "protein"],
                        help='Specify a specific database or database class to search', action='store')
    in_args = parser.parse_args()

    dbbuddy = []
    out_format = "summary" if not in_args.out_format else in_args.out_format
    search_set = ""

    try:
        if len(in_args.user_input) > 1:
            for search_set in in_args.user_input:
                dbbuddy.append(DbBuddy(search_set, in_args.database, out_format))

            dbbuddy = DbBuddy(dbbuddy, in_args.database, in_args.out_format)
        else:
            dbbuddy = DbBuddy(in_args.user_input[0], in_args.database, out_format)

    except GuessError:
        sys.exit("Error: SeqBuddy could not understand your input. "
                 "Check the file path or try specifying an input type with -f")

    # ############################################# INTERNAL FUNCTION ################################################ #
    def _print_recs(_dbbuddy):
        if in_args.test:
            _stderr("*** Test passed ***\n", in_args.quiet)
            pass

        else:
            if _dbbuddy.out_format != "summary":
                accession_only = [_accession for _accession, _rec in _dbbuddy.records.items() if not _rec.record]
                if len(accession_only) > 0:
                    _output = "# ################## Accessions without Records ################## #\n"
                    _counter = 1
                    for _next_acc in accession_only:
                        _output += "%s\t" % _next_acc
                        if _counter % 4 == 0:
                            _output = "%s\n" % _output.strip()
                        _counter += 1
                    _output = "%s\n# ################################################################ #\n\n" % _output.strip()
                    _stderr(_output, in_args.quiet)

            _stdout("{0}\n".format(str(_dbbuddy).rstrip()))

    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Download everything
    if in_args.download_everything:
        dbbuddy.out_format = "gb" if not in_args.out_format else in_args.out_format
        download_everything(dbbuddy)

        if len(dbbuddy.failures) > 0:
            output = "# ###################### Accession failures ###################### #\n"
            counter = 1
            for next_acc in dbbuddy.failures:
                output += "%s\t" % next_acc
                if counter % 4 == 0:
                    output = "%s\n" % output.strip()
                counter += 1
            _stderr("%s\n# ################################################################ #\n\n" % output.strip())

        _print_recs(dbbuddy)

    # Guess database  ToDo: Sort by database
    if in_args.guess_database:
        output = ""
        if len(dbbuddy.records) > 0:
            output += "# Accession\tDatabase\n"
            for accession, record in dbbuddy.records.items():
                output += "%s\t%s\n" % (accession, record.database)
            output += "\n"

        if len(dbbuddy.search_terms) > 0:
            output += "# Search terms\n"
            for term in dbbuddy.search_terms:
                output += "%s\n" % term

        if len(dbbuddy.records) == 0 and len(dbbuddy.search_terms) == 0:
            output += "Nothing to return\n"

        _stdout(output)
