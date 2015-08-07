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
from io import StringIO
import re
from urllib.parse import urlencode
from urllib.error import HTTPError
from urllib.request import Request, urlopen
from time import time, sleep
import json
from multiprocessing import Lock
from collections import OrderedDict
from hashlib import md5
import cmd

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
        sys.stderr.flush()
    return


def _stdout(message, quiet=False, color=None):
    if not quiet:
        if color and re.search("\\033\[9[0-6]m", color):
            sys.stdout.write(color)
        sys.stdout.write("%s\033[m" % message)
        sys.stdout.flush()
    return


def terminal_colors():
    colors = ['\033[95m', '\033[94m', '\033[92m', '\033[91m', '\033[93m', '\033[90m']
    _counter = 0
    while True:
        try:
            yield colors[_counter]
        except IndexError:
            _counter = 0
            yield colors[_counter]
        _counter += 1


# ##################################################### DB BUDDY ##################################################### #
class DbBuddy:  # Open a file or read a handle and parse, or convert raw into a Seq object
    def __init__(self, _input, _database=None, _out_format="summary"):
        if _database and _database.lower() not in ["gb", "genbank", "uniprot", "ensembl"]:
            raise DatabaseError("%s is not currently supported." % _database)

        if _database and _database.lower() == "genbank":
            _database = "gb"

        self.search_terms = []
        self.records = OrderedDict()
        self.recycle_bin = {}  # If records are filtered out, send them here instead of deleting them
        self.out_format = _out_format.lower()
        self.failures = {}
        self.databases = [_database] if _database else []

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
                    if _record.database not in self.databases:
                        self.databases.append(_record.database)

            # If accessions not identified, assume search terms
            if len(self.records) != len(accessions_check):
                search_term_check = re.sub("[\n\r,]+", "\t", _input)
                search_term_check = [x.strip() for x in search_term_check.split("\t")]
                for search_term in search_term_check:
                    if search_term not in self.records:
                        self.search_terms.append(search_term)

    def record_breakdown(self):
        _output = {x: [] for x in ["full", "partial", "accession"]}
        _output["full"] = [_accession for _accession, _rec in self.records.items() if _rec.record]
        _output["partial"] = [_accession for _accession, _rec in self.records.items() if not _rec.record and _rec.summary]
        _output["accession"] = [_accession for _accession, _rec in self.records.items() if not _rec.record and not _rec.summary]
        return _output

    def _print_recs(self, _num=0):
        _num = _num if _num > 0 else len(self.records)
        if in_args.test:
            _stderr("*** Test passed ***\n", in_args.quiet)
            pass

        else:
            # First deal with anything that broke or wasn't downloaded
            errors_etc = ""
            if len(self.failures):
                errors_etc += "# ########################## Failures ########################### #\n"
                for failure in self.failures:
                    errors_etc += str(failure)

            if len(self.record_breakdown()["accession"]) > 0:
                errors_etc = "# ################## Accessions without Records ################## #\n"
                _counter = 1
                for _next_acc in self.record_breakdown()["accession"]:
                    errors_etc += "%s\t" % _next_acc
                    if _counter % 4 == 0:
                        errors_etc = "%s\n" % errors_etc.strip()
                    _counter += 1
                errors_etc += "\n"

            if errors_etc != "":
                errors_etc = "%s\n# ################################################################ #\n\n" % errors_etc.strip()
                _stderr(errors_etc, in_args.quiet)

            _output = ""
            # Full records
            if self.out_format not in ["summary", "ids", "accessions"]:
                nuc_recs = [_rec.record for _accession, _rec in self.records.items() if _rec.type == "nucleotide" and _rec.record]
                prot_recs = [_rec.record for _accession, _rec in self.records.items() if _rec.type == "protein" and _rec.record]
                tmp_dir = TemporaryDirectory()
                if len(nuc_recs) > 0:
                    with open("%s/seqs.tmp" % tmp_dir.name, "w") as _ofile:
                        SeqIO.write(nuc_recs[:_num], _ofile, self.out_format)

                    with open("%s/seqs.tmp" % tmp_dir.name, "r") as ifile:
                        _output += "%s\n" % ifile.read()

                if len(prot_recs) > 0:
                    with open("%s/seqs.tmp" % tmp_dir.name, "w") as _ofile:
                        SeqIO.write(prot_recs[:_num], _ofile, self.out_format)

                    with open("%s/seqs.tmp" % tmp_dir.name, "r") as ifile:
                        _output += "%s\n" % ifile.read()
            # Summary outputs
            else:
                saved_headings = []
                for _accession, _rec in list(self.records.items())[:_num]:
                    colors = terminal_colors()
                    if self.out_format in ["ids", "accessions"]:
                        _output += "%s\n" % _accession

                    elif self.out_format == "summary":
                        headings = [heading for heading, _value in _rec.summary.items()]
                        if saved_headings != headings:
                            _output += "%sACCN\t%sDB" % (next(colors), next(colors))
                            for heading in headings:
                                _output += "\t%s%s" % (next(colors), heading)
                            _output += "\n"
                            saved_headings = list(headings)
                            colors = terminal_colors()

                        _output += "%s%s" % (next(colors), _accession)
                        if _rec.database:
                            _output += "\t%s%s" % (next(colors), _rec.database)
                        else:
                            next(colors)

                        for attrib, _value in _rec.summary.items():
                            if len(str(_value)) > 50:
                                _output += "\t%s%s..." % (next(colors), str(_value[:47]))
                            else:
                                _output += "\t%s%s" % (next(colors), _value)
                        _output += "\033[m\n"

            _stdout("{0}\n".format(_output.rstrip()))

    def __str__(self):
        _output = "### DatabaseBuddy object ###\n"
        _output += "DBs searched: %s\n" % " ,".join(self.databases)
        _output += "Out format: %s\n" % self.out_format
        _output += "Search term(s): "
        _output += "None\n" if not self.search_terms else "%s\n" % " ,".join(self.search_terms)

        breakdown = self.record_breakdown()
        _output += "Full Recs: %s\n" % len(breakdown["full"])
        _output += "Partial Recs: %s\n" % len(breakdown["partial"])
        _output += "ACCN only: %s\n" % len(breakdown["accession"])
        _output += "Filtered out: %s\n" % len(self.recycle_bin)
        _output += "Failures: %s\n" % len(self.failures)

        return _output


class Record:
    def __init__(self, _accession, _record=None, summary=None, _database=None, _type=None, _search_term=None):
        self.accession = _accession
        self.record = _record  # SeqIO record
        self.summary = summary if summary else OrderedDict()  # Dictionary of attributes
        self.database = _database  # refseq, genbank, ensembl, uniprot
        self.type = _type  # protein, nucleotide, gi_num
        self.search_term = _search_term  # In case the record was the result of a particular search

    def guess_database(self):
        # RefSeq
        if re.match("^[NX][MR]_[0-9]+", self.accession):
            self.database = "refseq"
            self.type = "nucleotide"

        elif re.match("^[ANYXZ]P_[0-9]+", self.accession):
            self.database = "refseq"
            self.type = "protein"

        # GenBank
        elif re.match("^[A-Z][0-9]{5}$|^[A-Z]{2}[0-9]{6}$", self.accession):  # Nucleotide
            self.database = "gb"
            self.type = "nucleotide"

        elif re.match("^[A-Z]{3}[0-9]{5}$", self.accession):  # Protein
            self.database = "gb"
            self.type = "protein"

        elif re.match("^[A-Z]{4}[0-9]{8,10}$", self.accession):  # Whole Genome
            self.database = "gb"
            self.type = "nucleotide"

        elif re.match("^[A-Z]{5}[0-9]{7}$", self.accession):  # MGA (Mass sequence for Genome Annotation)
            self.database = "gb"
            self.type = "protein"

        elif re.match("^[0-9]+$", self.accession):  # GI number
            self.database = "gb"
            self.type = "gi_num"  # Need to check genbank accession number to figure out that this is

        # UniProt/SwissProt
        elif re.match("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", self.accession):
            self.database = "uniprot"
            self.type = "protein"

        # Ensembl stable ids (http://www.ensembl.org/info/genome/stable_ids/index.html)
        elif re.match("^(ENS|FB)[A-Z]*[0-9]+", self.accession):
            self.database = "ensembl"
            self.type = "nucleotide"

        return

    def search(self, regex):
        for param in [self.accession, self.database, self.type, self.search_term]:
            if re.search(regex, param):
                return True
        for _key, _value in self.summary.items():
            if re.search(regex, _key) or re.search(regex, _value):
                return True
        if self.record:
            if re.search(regex, self.record.format("gb")):
                return True
        # If nothing hits, default to False
        return False

    def __str__(self):
        return "Accession:\t{0}\nDatabase:\t{1}\nRecord:\t{2}\nType:\t{3}\n".format(self.accession, self.database,
                                                                                    self.record, self.type)


class Failure:
    def __init__(self, query, error_message):
        self.query = query
        self.error_msg = error_message

        # Create a unique identifier for identification purposes
        self.hash = "%s%s" % (query, error_message)
        self.hash = md5(self.hash.encode()).hexdigest()

    def __str__(self):
        _output = "%s\n" % self.query
        _output += "%s\n" % self.error_msg
        return _output


# ################################################# Database Clients ################################################# #
class UniProtRestClient:
    # http://www.uniprot.org/help/uniprotkb_column_names
    # Limit URLs to 2,083 characters
    def __init__(self, _dbbuddy, server='http://www.uniprot.org/uniprot'):
        self.dbbuddy = _dbbuddy
        self.server = server
        self.lock = Lock()
        self.temp_dir = TempDir()
        self.http_errors_file = "%s/errors.txt" % self.temp_dir.path
        open(self.http_errors_file, "w").close()
        self.results_file = "%s/results.txt" % self.temp_dir.path
        open(self.results_file, "w").close()

    def mc_rest_query(self, _term, args):  # Multicore ready
        http_errors_file, results_file, request_params = args
        _term = re.sub(" ", "+", _term)
        request_string = ""
        for _param, _value in request_params.items():
            _value = re.sub(" ", "+", _value)
            request_string += "&{0}={1}".format(_param, _value)

        try:
            request = Request("{0}?query={1}{2}".format(self.server, _term, request_string))
            response = urlopen(request)
            response = response.read().decode()
            response = re.sub("^Entry.*\n", "", response, count=1)
            with self.lock:
                with open(results_file, "a") as ofile:

                    ofile.write("%s\n%s//\n" % (_term, response))
            return

        except HTTPError as e:
            with self.lock:
                with open(http_errors_file, "a") as ofile:
                    ofile.write("%s\n%s//\n" % (_term, e))
            return

    def filter_records(self, regex):
        for _id, _rec in self.dbbuddy.records.items():
            if not _rec.search(regex):
                self.dbbuddy.recycle_bin[_id] = _rec
        for _id in self.dbbuddy.records:
            if _id in self.dbbuddy.records:
                del self.dbbuddy.records[_id]
        return

    def _parse_error_file(self):
        with open(self.http_errors_file, "r") as ifile:
            http_errors_file = ifile.read().strip("//\n")
            if http_errors_file != "":
                http_errors_file = http_errors_file.split("//")
                for error in http_errors_file:
                    error = error.split("\n")
                    error = (error[0], "\n".join(error[1:]) if len(error) > 2 else (error[0], error[1]))
                    error = Failure(*error)
                    if error.hash not in self.dbbuddy.failures:
                        self.dbbuddy.failures[error.hash] = error
                open(self.http_errors_file, "w").close()
        return

    def count_hits(self):
        search_terms = "("
        search_terms += "+)OR(+".join(self.dbbuddy.search_terms)  # ToDo: make sure there isn't a size limit
        search_terms += ")"
        search_terms = re.sub(" ", "+", search_terms)
        self.mc_rest_query(search_terms, [self.http_errors_file, self.results_file, {"format": "list"}])
        with open(self.results_file, "r") as ifile:
            _count = len(ifile.read().strip().split("\n")[1:-1])  # The range clips off the search term and trailing //
        open(self.results_file, "w").close()

        self._parse_error_file()
        return _count

    def search_proteins(self, force=False):
        # start by determining how many results we would get from all searches.
        _stderr("Retrieving UniProt hit count... ")
        _count = self.count_hits()

        if _count == 0:
            _stderr("No hits\n")
            return

        else:
            _stderr("%s\n" % _count)

        # download the tab info on all or subset
        params = {"format": "tab", "columns": "id,entry name,length,organism-id,organism,protein names,comments"}
        if len(self.dbbuddy.search_terms) > 1:
            _stderr("Querying UniProt with %s search terms...\n" % len(self.dbbuddy.search_terms))
            run_multicore_function(self.dbbuddy.search_terms, self.mc_rest_query, max_processes=10,
                                   func_args=[self.http_errors_file, self.results_file, params])
        else:
            _stderr("Querying UniProt with the search term '%s'...\n" % self.dbbuddy.search_terms[0])
            self.mc_rest_query(self.dbbuddy.search_terms[0], [self.http_errors_file, self.results_file, params])

        with open(self.results_file, "r") as ifile:
            results = ifile.read().strip("//\n").split("//")

        for result in results:
            result = result.strip().split("\n")
            for hit in result[1:]:
                hit = hit.split("\t")
                if len(hit) == 6:  # In case 'comments' isn't returned
                    raw = OrderedDict([("entry_name", hit[1]), ("length", int(hit[2])), ("organism-id", hit[3]),
                                       ("organism", hit[4]), ("protein_names", hit[5]), ("comments", "")])
                else:
                    raw = OrderedDict([("entry_name", hit[1]), ("length", int(hit[2])), ("organism-id", hit[3]),
                                       ("organism", hit[4]), ("protein_names", hit[5]), ("comments", hit[6])])

                self.dbbuddy.records[hit[0]] = Record(hit[0], _database="uniprot", _type="protein",
                                                      _search_term=result[0], summary=raw)

    def fetch_proteins(self):
        _records = [_rec for _accession, _rec in self.dbbuddy.records.items() if _rec.database == "uniprot"]
        if len(_records) > 0:
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

    def gi2acc(self):
        gi_records = [_rec for _accession, _rec in self.dbbuddy.records.items() if _rec.type == "gi_num"]
        gi_accessions = [_rec.accession for _indx, _rec in enumerate(gi_records)]
        # record = Entrez.read(Entrez.fetch())
        # handle = Entrez.efetch(db="nucleotide", id=gi_accessions[:2], rettype="acc", retmode="text")
        print("Hello")
        handle = Entrez.esummary(db="nucleotide", id=["28864546", "28800981"])
        sys.exit(handle.read())

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

    def search_nucliotides(self):
        for _term in self.dbbuddy.search_terms:
            try:
                count = Entrez.read(Entrez.esearch(db='nucleotide', term=_term, rettype="count"))["Count"]
                handle = Entrez.esearch(db='nucleotide', term=_term, retmax=count)
                result = Entrez.read(handle)
                for _id in result["IdList"]:
                    if _id not in self.dbbuddy.records:
                        self.dbbuddy.records[_id] = Record(_id, _database="genbank", _type="gi_num")
                self.gi2acc()
                print(self.dbbuddy)
                sys.exit()
            except HTTPError as e:
                if e.getcode() == 400:
                    self.dbbuddy.failures.append(_term)
                elif e.getcode() == 503:
                    print("503 'Service unavailable': NCBI is either blocking you or they are "
                          "experiencing some technical issues.")
                else:
                    sys.exit("Error connecting with NCBI. Returned code %s" % e.getcode())


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
                # sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

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
class LiveSearch(cmd.Cmd):
    intro = "Welcome to the DatabaseBuddy live shell. Type 'help' at any time for a list of available commands " \
            "or 'help <command>' for more detailed explaination of a given command.\n"
    prompt = 'DbBuddy> '
    file = None
    dbbuddy = None

    def do_status(self, line):
        _stdout(str(self.dbbuddy))

    def do_show(self, line):
        if line:
            try:
                line = int(line)
                self.dbbuddy._print_recs(_num=line)
            except ValueError:
                _stdout("'%s' is not an integer, nothing displayed\n" % line)
        else:
            if len(self.dbbuddy.records) > 100:
                confirm = input("%s records currently in buffer, show them all (y/[n])?")
                if confirm.lower() in ["yes", "y"]:
                    self.dbbuddy._print_recs()
                else:
                    _stdout("Include an integer value with 'show' to return a specific number of records.\n")
            else:
                self.dbbuddy._print_recs()

    @staticmethod
    def help_status():
        _stdout("Display the current state of your buffer, including how many accessions and full records have been "
                "downloaded.\n\n", color="\033[92m")

    def do_format(self, line):
        self.dbbuddy.out_format = line

    def help_show(self):
        _stdout('''\
Output the records currently held in your buffer (out_format currently set to '\033[94m%s\033[92m')
Optionally include an integer value to limit how many will be shown.''' % self.dbbuddy.out_format, color="\033[92m")

    def help_format(self):
        _stdout('''\
Set the output format:
    ids or accessions ->  Simple list of all accessions in the buffer
    summary           ->  Extended information about each record
    <SeqIO format>    ->  Any sequence file format supported by BioPython (e.g. gb, fasta, clustal)
                          See http://biopython.org/wiki/SeqIO for details
''', color="\033[92m")

    def do_exit(self, line):
        return True

    def help_exit(self):
        _stdout("End the live session.\n\n")

    def do_quit(self, line):
        return True

    def help_quit(self):
        _stdout("End the live session.\n\n")

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    def do_search(self):
        pass

    def help_search(self):
        pass

    def do_filter(self):
        pass

    def help_filter(self):
        pass

    def do_download(self):
        pass

    def help_download(self):
        pass

    def do_write(self):
        pass

    def help_write(self):
        pass

    def do_database(self):
        pass

    def help_database(self):
        pass

    def do_type(self):
        pass

    def help_type(self):
        pass

    def do_reset(self):
        pass

    def help_reset(self):
        pass


def download_everything(_dbbuddy):
    # Get sequences from UniProt
    uniprot = UniProtRestClient(_dbbuddy)
    uniprot.fetch_proteins()

    # Get sequences from Ensembl
    ensembl = EnsemblRestClient(_dbbuddy)
    ensembl.fetch_nucleotides()

    # Get sequences from genbank
    refseq = NCBIClient(_dbbuddy)
    refseq.gi2acc()
    refseq.fetch_nucliotides()
    refseq.fetch_proteins()

    return _dbbuddy


def retrieve_accessions(_dbbuddy):
    check_all = True if len(_dbbuddy.databases) == 0 else False

    if "gb" in _dbbuddy.databases or check_all:
        refseq = NCBIClient(_dbbuddy)
        refseq.gi2acc()
        # refseq.search_nucliotides()

    elif "uniprot" in _dbbuddy.databases or check_all:
        uniprot = UniProtRestClient(_dbbuddy)
        uniprot.search_proteins()

    return _dbbuddy


def retrieve_summary(_dbbuddy):
    check_all = True if len(_dbbuddy.databases) == 0 else False

    if "gb" in _dbbuddy.databases or check_all:
        refseq = NCBIClient(_dbbuddy)
        refseq.gi2acc()
        # refseq.search_nucliotides()

    elif "uniprot" in _dbbuddy.databases or check_all:
        uniprot = UniProtRestClient(_dbbuddy)
        uniprot.search_proteins()

    return _dbbuddy


def retrieve_sequences(_dbbuddy):
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
    parser.add_argument('-ls', '--live_search', action="store_true",
                        help="Interactive database searching. The best tool for sequence discovery.")
    parser.add_argument('-ra', '--retrieve_accessions', action="store_true",
                        help="Use search terms to find a list of sequence accession numbers")
    parser.add_argument('-rs', '--retrieve_sequences', action="store_true",
                        help="Get sequences for every included accession")
    parser.add_argument('-gd', '--guess_database', action="store_true",
                        help="List the database that each provided accession belongs to.")

    parser.add_argument('-q', '--quiet', help="Suppress stderr messages", action='store_true')
    parser.add_argument('-t', '--test', action='store_true',
                        help="Run the function and return any stderr/stdout other than sequences.")
    parser.add_argument('-o', '--out_format', help="If you want a specific format output", action='store')  # accessions/ids, summary, SeqIO formats
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

    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Live Search
    if in_args.live_search:
        if dbbuddy.records:
            retrieve_sequences(dbbuddy)

        if dbbuddy.search_terms:
            retrieve_accessions(dbbuddy)

        live_search = LiveSearch()
        live_search.dbbuddy = dbbuddy
        live_search.cmdloop()
        sys.exit()

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

        dbbuddy._print_recs()
        sys.exit()

    # Retrieve Accessions
    if in_args.retrieve_accessions:
        if not in_args.out_format:
            dbbuddy.out_format = "ids"
        retrieve_accessions(dbbuddy)
        dbbuddy._print_recs()
        sys.exit()

    if in_args.retrieve_sequences:
        sys.exit()

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
        sys.exit()
