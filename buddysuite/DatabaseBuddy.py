#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: DbBuddy.py
version: 1.0
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description:
Collection of functions that interact with public sequence databases. Pull them into a script or run from command line.
"""

# ##################################################### IMPORTS ###################################################### #
# Standard library
# import pdb
# import timeit
from __future__ import print_function

# BuddySuite specific
import buddy_resources as br
import MyFuncs

# Standard library
import sys
import os
import re
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
from time import sleep
import json
from multiprocessing import Lock
from collections import OrderedDict
from hashlib import md5
import cmd
from subprocess import Popen, PIPE
from io import TextIOWrapper, StringIO
import warnings
import readline
import pickle

# Third party
# sys.path.insert(0, "./")  # For stand alone executable, where dependencies are packaged with BuddySuite
from Bio import Entrez
from Bio import SeqIO
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


# ##################################################### WISH LIST #################################################### #
# - Catch memory limits before they are overblown by big fetches. Python std-lib doesn't seem to have a tool for this.
# - Push some functionality off to separate threads, like session dump and stats usage
# - BLAST functionality
# - Load sequence files up into the live session (very useful for BLAST, when implemented)

# #################################################### CHANGE LOG #################################################### #
# ###################################################### GLOBALS ##################################################### #
TRASH_SYNOS = ["t", "tb", "t_bin", "tbin", "trash", "trashbin", "trash-bin", "trash_bin"]
RECORD_SYNOS = ["r", "rec", "recs", "records", "main", "filtered"]
SEARCH_SYNOS = ["st", "search", "search-terms", "search_terms", "terms"]
DATABASES = ["ncbi_nuc", "ncbi_prot", "uniprot", "ensembl"]
RETRIEVAL_TYPES = ["protein", "nucleotide", "gi_num"]
FORMATS = ["ids", "accessions", "summary", "full-summary", "clustal", "embl", "fasta", "fastq", "fastq-sanger",
           "fastq-solexa", "fastq-illumina", "genbank", "gb", "imgt", "nexus", "phd", "phylip", "seqxml", "sff",
           "stockholm", "tab", "qual"]
CONFIG = br.config_values()
VERSION = br.Version("DatabaseBuddy", 1, 0, br.contributors)

GREY = "\033[90m"
RED = "\033[91m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
BLUE = "\033[94m"
MAGENTA = "\033[95m"
CYAN = "\033[96m"
WHITE = "\033[97m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
NO_UNDERLINE = "\033[24m"
DEF_FONT = "\033[39m"


# ##################################################### DB BUDDY ##################################################### #
class DbBuddy(object):  # Open a file or read a handle and parse, or convert raw into a Seq object
    def __init__(self, _input=None, _databases=None, _out_format="summary"):
        self.search_terms = []
        self.records = OrderedDict()  # Record objects
        self.trash_bin = {}  # If records are filtered out, send them here instead of deleting them
        self.out_format = _out_format.lower()
        self.failures = {}  # The key for these is a hash of the Failure, and the values are actual Failure objects
        self.databases = check_database(_databases)
        _databases = self.databases[0] if len(self.databases) == 1 else None  # This is to check if a specific db is set
        self.server_clients = {"ncbi": False, "ensembl": False, "uniprot": False}

        # Empty DbBuddy object
        if not _input:
            pass

        # DbBuddy objects
        elif type(_input) == list:
            for _dbbuddy in _input:
                if type(_dbbuddy) != DbBuddy:
                    raise TypeError("List of non-DbBuddy objects passed into DbBuddy as _input. %s" % _dbbuddy)

                self.search_terms += _dbbuddy.search_terms
                # ToDo: update() will overwrite any common records between the two dicts,
                # should check whether values are already set first
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
            with open(_input, "r", encoding="utf-8") as _ifile:
                _input = _ifile.read().strip()

        else:
            raise br.GuessError("DbBuddy could not determine the input type.")

        if _input:
            # try to glean accessions first
            accessions_check = re.sub("[\n\r, ]+", "\t", _input)
            accessions_check = accessions_check.split("\t")
            for _accession in accessions_check:
                _record = Record(_accession)
                _record.guess_database()
                if _record.database:
                    if _databases:
                        _record.database = _databases
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

    def __hash__(self):
        _records = tuple([(_key, _value) for _key, _value in self.records.items()])
        return hash(_records) ^ hash(self.out_format)  # The ^ is bitwise XOR, returning a string of bits

    def __eq__(self, other):
        return isinstance(other, type(self)) and ((self.records, self.out_format) == (other.records, other.out_format))

    def __str__(self):
        _output = "############################\n"
        _output += "### DatabaseBuddy object ###\n"
        _output += "Databases:    %s\n" % ", ".join(self.databases)
        _output += "Out format:   %s\n" % self.out_format
        _output += "Searches:     "
        _output += "None\n" if not self.search_terms else "%s\n" % ", ".join(self.search_terms)

        breakdown = self.record_breakdown()
        _output += "Full Recs:    %s\n" % len(breakdown["full"])
        _output += "Summary Recs: %s\n" % len(breakdown["summary"])
        _output += "ACCN only:    %s\n" % len(breakdown["accession"])
        _output += "Trash bin:  %s\n" % len(self.trash_bin)
        _output += "Failures:     %s\n" % len(self.failures)
        _output += "############################\n"
        return _output

    def filter_records(self, regex, mode):
        if mode not in ["keep", "remove", "restore"]:
            raise ValueError("The 'mode' argument in filter() must be 'keep', 'remove', or 'restore', not %s." % mode)

        column_errors = {"KeyError": [], "ValueError": []}
        for _id, _rec in self.trash_bin.items() if mode == 'restore' else self.records.items():
            try:
                if mode == "keep" and not _rec.search(regex):
                    self.trash_bin[_id] = _rec
                elif mode == "remove" and _rec.search(regex):
                    self.trash_bin[_id] = _rec
                elif mode == "restore" and _rec.search(regex):
                    self.records[_id] = _rec
            except KeyError as _e:
                if str(_e) not in column_errors["KeyError"]:
                    column_errors["KeyError"].append(str(_e))
            except ValueError as _e:
                if str(_e) not in column_errors["ValueError"]:
                    column_errors["ValueError"].append(str(_e))

        if mode == "restore":
            for _id in self.records:
                if _id in self.trash_bin:
                    del self.trash_bin[_id]
        else:
            for _id in self.trash_bin:
                if _id in self.records:
                    del self.records[_id]

        return column_errors

    def record_breakdown(self):
        _output = {x: [] for x in ["full", "summary", "accession"]}
        _output["full"] = [_accession for _accession, _rec in self.records.items() if _rec.record]
        _output["summary"] = [_accession for _accession, _rec in self.records.items()
                              if not _rec.record and _rec.summary]
        _output["accession"] = [_accession for _accession, _rec in self.records.items()
                                if not _rec.record and not _rec.summary]
        return _output

    def server(self, _server):  # _server in ["uniprot", "ncbi", "ensembl"]
        if self.server_clients[_server]:
            return self.server_clients[_server]
        if _server == "uniprot":
            client = UniProtRestClient(self)
        elif _server == "ncbi":
            client = NCBIClient(self)
        elif _server == "ensembl":
            client = EnsemblRestClient(self)
        else:
            raise ValueError('"uniprot", "ncbi", and "ensembl" are the only valid options, not %s' % _server)
        self.server_clients[_server] = client
        return client

    def trash_breakdown(self):
        _output = {x: [] for x in ["full", "summary", "accession"]}
        _output["full"] = [_accession for _accession, _rec in self.trash_bin.items() if _rec.record]
        _output["summary"] = [_accession for _accession, _rec in self.trash_bin.items()
                              if not _rec.record and _rec.summary]
        _output["accession"] = [_accession for _accession, _rec in self.trash_bin.items()
                                if not _rec.record and not _rec.summary]
        return _output

    def print(self, _num=0, quiet=False, columns=None, destination=None, group="records"):
        """
        :param _num: Limit the number of rows (records) returned, otherwise everything is output
        :param quiet: suppress stderr
        :param columns: Variable, list of column names to include in summary output
        :param destination: a file path or handle to write to
        :param group: Either 'records' or 'trash_bin'
        :return: Nothing.
        """
        group = self.trash_bin if group == "trash_bin" else self.records

        _num = _num if _num > 0 else len(group)

        # First deal with anything that broke or wasn't downloaded
        errors_etc = ""
        if len(self.failures) > 0:
            errors_etc += "# ########################## Failures ########################### #\n"
            for _hash, failure in self.failures.items():
                errors_etc += str(failure)

        if len(self.record_breakdown()["accession"]) > 0 and _num == len(group):
            errors_etc += "# ################## Accessions without Records ################## #\n"
            _counter = 1
            for _next_acc in self.record_breakdown()["accession"]:
                errors_etc += "%s\t" % _next_acc
                if _counter % 4 == 0:
                    errors_etc = "%s\n" % errors_etc.strip()
                _counter += 1
            errors_etc += "\n"

        if errors_etc != "":
            errors_etc = "%s\n# ################################################################ #\n\n" \
                         % errors_etc.strip()
            _stderr(errors_etc, quiet)

        _output = ""
        # Summary outputs
        if self.out_format in ["summary", "full-summary", "ids", "accessions"]:
            def pad_columns(line_group, col_widths, all_lines):
                for indx_x, next_line in enumerate(line_group):
                    for indx_y, _col in enumerate(next_line):
                        line_group[indx_x][indx_y] = str(_col).ljust(col_widths[indx_y] + 2)
                all_lines += line_group
                all_lines.append([])
                return all_lines

            lines = []
            current_group = []
            saved_headings = []
            column_widths = []
            for _accession, _rec in list(group.items())[:_num]:
                if self.out_format in ["ids", "accessions"]:
                    lines.append([_accession])

                elif self.out_format in ["summary", "full-summary"]:
                    headings = ["ACCN", "DB"]
                    headings += [heading for heading, _value in _rec.summary.items()]
                    headings += ["record"]

                    if columns:
                        headings = [heading for heading in headings if heading in columns]

                    if saved_headings != headings:
                        if saved_headings:
                            lines = pad_columns(current_group, column_widths, lines)
                            column_widths = []
                            current_group = []

                        for heading in headings:
                            column_widths.append(len(str(heading)))
                        current_group.append(headings)
                        saved_headings = list(headings)

                    current_group.append([])
                    attrib_counter = 0
                    if "ACCN" in headings:
                        current_group[-1].append(_accession)
                        if len(str(_accession)) > column_widths[attrib_counter]:
                            column_widths[attrib_counter] = len(str(_accession))
                        attrib_counter += 1

                    if "DB" in headings:
                        if _rec.database:
                            current_group[-1].append(_rec.database)
                            if len(_rec.database) > column_widths[attrib_counter]:
                                column_widths[attrib_counter] = len(_rec.database)
                        else:
                            current_group[-1].append("")
                        attrib_counter += 1

                    for attrib, _value in _rec.summary.items():
                        if attrib in headings:
                            if len(str(_value)) > 50 and self.out_format != "full-summary":
                                current_group[-1].append("%s..." % _value[:47])
                                column_widths[attrib_counter] = 50
                            else:
                                current_group[-1].append(_value)
                                if len(str(_value)) > column_widths[attrib_counter]:
                                    column_widths[attrib_counter] = len(str(_value))
                            attrib_counter += 1

                    if self.out_format not in ["ids", "accessions"] and "record" in headings:
                        if _rec.record:
                            current_group[-1].append("full")
                        else:
                            current_group[-1].append("summary")

            if self.out_format in ["summary", "full-summary"] and current_group:
                lines = pad_columns(current_group, column_widths, lines)

            _output = "\033[m\033[40m\033[97m"
            for line in lines:
                colors = terminal_colors()
                _output += "%s\n" % "".join(["%s%s" % (next(colors), _col) for _col in line])

        # Full records
        else:
            nuc_recs = [_rec.record for _accession, _rec in group.items() if
                        _rec.type == "nucleotide" and _rec.record]
            prot_recs = [_rec.record for _accession, _rec in group.items() if
                         _rec.type == "protein" and _rec.record]
            tmp_dir = MyFuncs.TemporaryDirectory()
            if len(nuc_recs) > 0:
                with open("%s/seqs.tmp" % tmp_dir.name, "w", encoding="utf-8") as _ofile:
                    SeqIO.write(nuc_recs[:_num], _ofile, self.out_format)

                with open("%s/seqs.tmp" % tmp_dir.name, "r", encoding="utf-8") as ifile:
                    _output += "%s\n" % ifile.read()

            if len(prot_recs) > 0:
                with open("%s/seqs.tmp" % tmp_dir.name, "w", encoding="utf-8") as _ofile:
                    SeqIO.write(prot_recs[:_num], _ofile, self.out_format)

                with open("%s/seqs.tmp" % tmp_dir.name, "r", encoding="utf-8") as ifile:
                    _output += "%s\n" % ifile.read()

        if not destination:
            _stdout("{0}\n".format(_output.rstrip()))
        else:
            # remove any escape characters and convert space padding to tabs if writing the file
            _output = re.sub("\\033\[[0-9]*m", "", _output)
            _output = re.sub("  +", "\t", _output)
            destination.write(_output)


# ################################################# SUPPORT CLASSES ################################################## #
class Record(object):
    def __init__(self, _accession, gi=None, _version=None, _record=None, summary=None, _size=None,
                 _database=None, _type=None, _search_term=None):
        self.accession = _accession
        self.gi = gi  # This is for NCBI records
        self.version = _version
        self.record = _record  # SeqIO record
        self.summary = summary if summary else OrderedDict()  # Dictionary of attributes
        self.size = _size if _size in [None, ''] else int(_size)
        self.database = _database
        self.type = check_type(_type)
        self.search_term = _search_term  # In case the record was the result of a particular search

    def ncbi_accn(self):
        if not self.version:
            return self.accession
        else:
            return "%s.%s" % (self.accession, self.version)

    def guess_database(self):
        # RefSeq
        if re.match("^[NX][MR]_[0-9]+", self.accession):
            self.database = "ncbi_nuc"
            self.type = "nucleotide"

        if re.match("^[NX][C]_[0-9]+", self.accession):  # Chromosome
            self.database = "ncbi_nuc"
            self.type = "nucleotide"

        elif re.match("^[ANYXZ]P_[0-9]+", self.accession):
            self.database = "ncbi_prot"
            self.type = "protein"

        # UniProt/SwissProt
        elif re.match("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", self.accession):
            self.database = "uniprot"
            self.type = "protein"

        # Ensembl stable ids (http://www.ensembl.org/info/genome/stable_ids/index.html)
        elif re.match("^(ENS|FB)[A-Z]*[0-9]+", self.accession):
            self.database = "ensembl"
            self.type = "nucleotide"

        # GenBank
        elif re.match("^[A-Z][0-9]{5}$|^[A-Z]{2}[0-9]{6}$", self.accession):  # Nucleotide
            self.database = "ncbi_nuc"
            self.type = "nucleotide"

        elif re.match("^[A-Z]{3}[0-9]{5}$", self.accession):  # Protein
            self.database = "ncbi_prot"
            self.type = "protein"

        elif re.match("[0-9][A-Z0-9]{3}(_[A-Z0-9])?$", self.accession):  # PDB
            self.database = "ncbi_prot"
            self.type = "protein"

        elif re.match("^[A-Z]{4}[0-9]{8,10}$", self.accession):  # Whole Genome
            self.database = "ncbi_nuc"
            self.type = "nucleotide"

        elif re.match("^[A-Z]{5}[0-9]{7}$", self.accession):  # MGA (Mass sequence for Genome Annotation)
            self.database = "ncbi_prot"
            self.type = "protein"

        elif re.match("^[0-9]+$", self.accession):  # GI number
            self.database = "ncbi_nuc"
            self.type = "gi_num"  # Need to check genbank accession number to figure out what this is
            self.gi = str(self.accession)

        if self.database in ["ncbi_nuc", "ncbi_prot"]:
            # Catch accn/version
            with_version = re.search("^(.*?)\.([0-9]+)$", self.accession)
            if with_version:
                accn, ver = with_version.group(1), with_version.group(2)
                self.accession = accn
                self.version = ver

        # ToDo: This is for testing, needs to be removed for production
        # else:
        #    raise TypeError("Unable to guess database for accession '%s'" % self.accession)

        return

    def search(self, regex):
        regex = "." if regex == "*" else regex
        column = re.match("\((.*)\)", regex)
        if column:
            column = column.group(1)
            # Special case, if user is searching sequence length
            if re.match("length.+", column):
                if not re.match("^length[ =<>]*[0-9]*$", column):
                    raise ValueError("Invalid syntax for seaching 'length': %s" % column)

                limit = re.search("length[ =<>]*([0-9]*)", column).group(1)
                try:
                    limit = int(limit)
                except ValueError:
                    limit = re.search("length[ =<>]*(.*)", column).group(1)[:-1]
                    raise ValueError("Unable to recast limit: %s" % limit)

                operator = re.search("length *([=<>]*) *[0-9]", column).group(1)
                if operator not in ["=", ">", ">=", "<", "<="]:
                    raise ValueError("Invalid operator: %s" % operator)

                if operator == "=" and int(self.summary["length"]) == limit:
                    return True
                elif operator == ">" and int(self.summary["length"]) > limit:
                    return True
                elif operator == ">=" and int(self.summary["length"]) >= limit:
                    return True
                elif operator == "<" and int(self.summary["length"]) < limit:
                    return True
                elif operator == "<=" and int(self.summary["length"]) <= limit:
                    return True
                else:
                    return False

            column = "?i" if column == "i?" else column  # Catch an easy syntax error. Maybe a bad idea?
            regex = regex[len(column) + 2:] if column != "?i" else regex
            try:
                if re.search(str(regex), str(self.summary[column])):
                    return True
                else:
                    return False
            except KeyError:
                if column != "?i":
                    raise KeyError(column)

        for param in [self.accession, self.database, self.type, self.search_term]:
            if re.search(regex, str(param)):
                return True

        for _key, _value in self.summary.items():
            if re.search(regex, _key) or re.search(regex, str(_value)):
                return True

        if self.record:
            if re.search(regex, self.record.format("gb")):
                return True
        # If nothing hits, default to False
        return False

    def update(self, new_rec):
        # ToDo: automate this by looping through Record object dir()
        self.accession = new_rec.accession if new_rec.accession else self.accession
        self.gi = new_rec.gi if new_rec.gi else self.gi
        self.version = new_rec.version if new_rec.version else self.version
        self.record = new_rec.record if new_rec.record else self.record
        self.summary = new_rec.summary if new_rec.summary else self.summary
        self.size = new_rec.size if new_rec.size else self.size
        self.database = new_rec.database if new_rec.database else self.database
        self.type = new_rec.type if new_rec.type else self.type
        self.search_term = new_rec.search_term if new_rec.search_term else self.search_term

    def __str__(self):
        return "Accession:\t{0}\nDatabase:\t{1}\nRecord:\t{2}\nType:\t{3}\n".format(self.accession, self.database,
                                                                                    self.record, self.type)


class Failure(object):
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


# ################################################# HELPER FUNCTIONS ################################################# #
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


def _stdout(message, quiet=False, format_in=None, format_out=None):
    if format_in and type(format_in) == list:
        format_in = "".join(format_in)
    if format_out and type(format_out) == list:
        format_out = "".join(format_out)
    if not quiet:
        if format_in and re.search("\\033\[[0-9]*m", format_in):
            sys.stdout.write(format_in)
        if format_out and re.search("\\033\[[0-9]*m", format_out):
            sys.stdout.write("%s%s" % (message, format_out))
        else:
            sys.stdout.write("%s\033[m" % message)
        sys.stdout.flush()
    return


def terminal_colors():
    colors = [MAGENTA, CYAN, GREEN, RED, YELLOW, GREY]
    _counter = 0
    while True:
        try:
            yield colors[_counter]
        except IndexError:
            _counter = 0
            yield colors[_counter]
        _counter += 1


def check_database(_database):
    _output = []
    if type(_database) == list:
        for _db in [x.lower() for x in _database]:
            if _db == "all":
                _output = DATABASES
                break
            elif _db in DATABASES:
                _output.append(_db)
            else:
                _stderr("Warning: '%s' is not a valid database choice, omitted.\n" % _db)

        if not _output:
            _stderr("Warning: No valid database choice provided. Setting to default 'all'.\n")
            _output = DATABASES

    elif not _database:
        _output = DATABASES
    else:
        _database = _database.lower()
        if _database == "all":
            _output = DATABASES
        elif _database in DATABASES:
            _output = [_database]
        else:
            _stderr("Warning: '%s' is not a valid database choice. Setting to default 'all'.\n")
            _output = DATABASES
    return _output


def check_type(_type):
    _type = None if not _type else _type.lower()
    if _type in ["p", "pr", "prt", "prtn", "prn", "prot", "protn", "protien", "protein"]:
        _type = "protein"
    elif _type in ["n", "ncl", "nuc", "dna", "nt", "gene", "transcript", "nucleotide"]:
        _type = "nucleotide"
    elif _type in ["g", "gi", "gn", "gin", "gi_num", "ginum", "gi_number"]:
        _type = "gi_num"

    if _type and _type not in RETRIEVAL_TYPES:
        _stderr("Warning: '%s' is not a valid choice for '_type'. Setting to default 'protein'.\n" % _type)
        _type = "protein"
    return _type


# ################################################# Database Clients ################################################# #
class UniProtRestClient(object):
    # http://www.uniprot.org/help/uniprotkb_column_names
    def __init__(self, _dbbuddy, server='http://www.uniprot.org/uniprot'):
        self.dbbuddy = _dbbuddy
        self.server = server
        self.temp_dir = MyFuncs.TempDir()
        self.http_errors_file = "%s/errors.txt" % self.temp_dir.path
        open(self.http_errors_file, "w", encoding="utf-8").close()
        self.results_file = "%s/results.txt" % self.temp_dir.path
        open(self.results_file, "w", encoding="utf-8").close()
        self.max_url = 1000

    def query_uniprot(self, _term, args):  # Multicore ready
        http_errors_file, results_file, request_params, lock = args
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
            with lock:
                with open(results_file, "a", encoding="utf-8") as ofile:

                    ofile.write("# Search: %s\n%s//\n" % (_term, response))
            return

        except HTTPError as _e:
            with lock:
                with open(http_errors_file, "a", encoding="utf-8") as ofile:
                    ofile.write("%s\n%s//\n" % (_term, _e))

        except URLError as _e:
            with lock:
                with open(http_errors_file, "a", encoding="utf-8") as ofile:
                    ofile.write("%s\n%s//\n" % (_term, _e))

        except KeyboardInterrupt:
            _stderr("\n\tUniProt query interrupted by user\n")

    def _parse_error_file(self):
        with open(self.http_errors_file, "r", encoding="utf-8") as ifile:
            http_errors_file = ifile.read().strip("//\n")
        if http_errors_file != "":
            _output = ""
            http_errors_file = http_errors_file.split("//")
            for error in http_errors_file:
                error = error.split("\n")
                error = (error[0], "\n".join(error[1:])) if len(error) > 2 else (error[0], error[1])
                error = Failure(*error)
                if error.hash not in self.dbbuddy.failures:
                    self.dbbuddy.failures[error.hash] = error
                    _output += "%s\n" % error
            open(self.http_errors_file, "w", encoding="utf-8").close()
            return _output  # Errors found
        else:
            return False  # No errors to report

    def count_hits(self):
        # Limit URLs to 2,083 characters
        search_terms = []
        for _term in self.dbbuddy.search_terms:
            if len(_term) > self.max_url:
                raise ValueError("Search term exceeds size limit of %s characters." % self.max_url)

            _term = "(%s)" % _term  # Parentheses to keep search terms together
            _term = re.sub(" ", "+", _term)
            if not search_terms:
                search_terms.append(_term)

            elif (len(search_terms[-1]) + len(_term) + 4) <= self.max_url:
                search_terms[-1] += "+OR+%s" % _term

            else:
                search_terms.append(_term)

        search_terms = search_terms[0] if len(search_terms) == 1 else search_terms
        if not search_terms:
            return 0
        self.query_uniprot(search_terms, [self.http_errors_file, self.results_file, {"format": "list"}, Lock()])
        with open(self.results_file, "r", encoding="utf-8") as ifile:
            _count = len(ifile.read().strip().split("\n")[1:-1])  # The range clips off the search term and trailing //
        open(self.results_file, "w", encoding="utf-8").close()

        errors = self._parse_error_file()
        if errors:
            _stderr("{0}{1}The following errors were encountered while querying UniProt with "
                    "count_hits():{2}\n\n{3}{4}".format(RED, UNDERLINE, NO_UNDERLINE, errors, DEF_FONT))
        return _count

    def search_proteins(self):
        # start by determining how many results we would get from all searches.
        open(self.results_file, "w", encoding="utf-8").close()
        _count = self.count_hits()

        if _count == 0:
            _stderr("Uniprot returned no results\n\n")
            return

        else:
            _stderr("Retrieving summary data for %s records from UniProt\n" % _count)

        # download the tab info on all or subset
        params = {"format": "tab", "columns": "id,entry name,length,organism-id,organism,protein names,comments"}
        runtime = MyFuncs.RunTime(prefix="\t")
        if len(self.dbbuddy.search_terms) > 1:
            _stderr("Querying UniProt with %s search terms (Ctrl+c to abort)\n" % len(self.dbbuddy.search_terms))
            runtime.start()
            MyFuncs.run_multicore_function(self.dbbuddy.search_terms, self.query_uniprot, max_processes=10, quiet=True,
                                           func_args=[self.http_errors_file, self.results_file, params, Lock()])
        else:
            _stderr("Querying UniProt with the search term '%s'...\n" % self.dbbuddy.search_terms[0])
            runtime.start()
            self.query_uniprot(self.dbbuddy.search_terms[0], [self.http_errors_file, self.results_file, params, Lock()])
        runtime.end()
        errors = self._parse_error_file()
        if errors:
            _stderr("{0}{1}The following errors were encountered while querying UniProt with "
                    "search_proteins():{2}\n\n{3}{4}".format(RED, UNDERLINE, NO_UNDERLINE, errors, DEF_FONT))

        with open(self.results_file, "r", encoding="utf-8") as ifile:
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
                                                      _search_term=result[0], summary=raw, _size=int(hit[2]))
        _stderr("\n")

    def fetch_proteins(self):
        open(self.results_file, "w", encoding="utf-8").close()
        _records = [_rec for _accession, _rec in self.dbbuddy.records.items() if
                    _rec.database == "uniprot" and _rec.database == "uniprot" and not _rec.record]

        if len(_records) > 0:
            _stderr("Retrieving %s full records from UniProt...\n" % len(_records))
            accessions = [_records[0].accession]
            for _rec in _records[1:]:
                if len(accessions[-1]) + len(_rec.accession) + 1 > self.max_url:
                    accessions.append(_rec.accession)
                else:
                    accessions[-1] += ",%s" % _rec.accession

            runtime = MyFuncs.RunTime(prefix="\t")
            runtime.start()
            params = {"format": "txt"}
            MyFuncs.run_multicore_function(accessions, self.query_uniprot, max_processes=10, quiet=True,
                                           func_args=[self.http_errors_file, self.results_file, params, Lock()])
            runtime.end()
            errors = self._parse_error_file()
            if errors:
                _stderr("{0}{1}The following errors were encountered while querying UniProt with "
                        "fetch_proteins():{2}\n{3}{4}".format(RED, UNDERLINE, NO_UNDERLINE, errors, DEF_FONT))

            with open(self.results_file, "r", encoding="utf-8") as ifile:
                data = ifile.read().strip().split("//\n//")

            if data[0] == "":
                _stderr("No sequences returned\n\n")
                return

            clean_recs = []
            for _rec in data:
                if _rec:
                    # Strip the first line from multi-core searches
                    clean_recs.append(re.sub("# Search.*\n", "", _rec.strip()))

            with open(self.results_file, "w", encoding="utf-8") as ifile:
                ifile.write("//\n".join(clean_recs))
                ifile.write("\n//")

            with open(self.results_file, "r", encoding="utf-8") as ifile:
                _records = SeqIO.parse(ifile, "swiss")
                for _rec in _records:
                    if _rec.id not in self.dbbuddy.records:
                        # ToDo: fix failures
                        print(_rec.id)
                        self.dbbuddy.failures.setdefault("# Uniprot fetch: Ids not"
                                                         " in dbbuddy.records", []).append(_rec.id)
                    else:
                        self.dbbuddy.records[_rec.id].record = _rec


class NCBIClient(object):
    def __init__(self, _dbbuddy):
        Entrez.email = CONFIG["email"]
        Entrez.tool = "buddysuite"
        self.dbbuddy = _dbbuddy
        self.temp_dir = MyFuncs.TempDir()
        self.http_errors_file = "%s/errors.txt" % self.temp_dir.path
        open(self.http_errors_file, "w", encoding="utf-8").close()
        self.results_file = "%s/results.txt" % self.temp_dir.path
        open(self.results_file, "w", encoding="utf-8").close()
        self.max_url = 1000
        self.max_attempts = 5  # NCBI throws a lot of 503 errors, so keep trying until we get through...

    def _clear_files(self):
        open(self.http_errors_file, "w", encoding="utf-8").close()
        open(self.results_file, "w", encoding="utf-8").close()
        return

    def _parse_error_file(self):
        with open(self.http_errors_file, "r", encoding="utf-8") as ifile:
            http_errors_file = ifile.read().strip("//\n")
        if http_errors_file != "":
            _output = ""
            http_errors_file = http_errors_file.split("//")
            for error in http_errors_file:
                error = error.split("\n")
                error = (error[0], "\n".join(error[1:])) if len(error) > 2 else (error[0], error[1])
                error = Failure(*error)
                if error.hash not in self.dbbuddy.failures:
                    self.dbbuddy.failures[error.hash] = error
                    _output += "%s\n" % error
            open(self.http_errors_file, "w", encoding="utf-8").close()
            return _output  # Errors found
        else:
            return False  # No errors to report

    def _split_for_url(self, accessions):
        _groups = [""]
        for accn in accessions:
            accn = str(accn)
            if _groups[-1] == "":
                _groups[-1] = accn
            elif len(_groups[-1] + accn) + 1 <= self.max_url:
                _groups[-1] += ",%s" % accn
            else:
                _groups.append(accn)
        return _groups

    def _mc_taxa(self, _taxa_ids, args):
        lock = args[0]
        error = False
        handle = False
        timer = MyFuncs.time()
        for i in range(self.max_attempts):
            try:
                handle = Entrez.esummary(db="taxonomy", id=_taxa_ids, retmax=10000)
                '''
                Example output: esummary.fcgi?db=taxonomy&id=649
                    <eSummaryResult>
                        <DocSum>
                            <Id>649</Id>
                            <Item Name="Status" Type="String">active</Item>
                            <Item Name="Rank" Type="String">species</Item>
                            <Item Name="Division" Type="String">g-proteobacteria</Item>
                            <Item Name="ScientificName" Type="String">Aeromonas eucrenophila</Item>
                            <Item Name="CommonName" Type="String"/>
                            <Item Name="TaxId" Type="Integer">649</Item>
                            <Item Name="AkaTaxId" Type="Integer">0</Item>
                            <Item Name="Genus" Type="String">Aeromonas</Item>
                            <Item Name="Species" Type="String">eucrenophila</Item>
                            <Item Name="Subsp" Type="String"/>
                            <Item Name="ModificationDate" Type="Date">2014/12/30 00:00</Item>
                        </DocSum>
                    </eSummaryResult>
                '''
                timer = MyFuncs.time() - timer
                if timer < 1:
                    sleep(1 - timer)
                break
            except HTTPError as _e:
                if i == self.max_attempts - 1:
                    error = _e
                sleep(1)

        with lock:
            if error:
                with open(self.http_errors_file, "r", encoding="utf-8") as ifile:
                    ifile.write("%s\n%s//\n" % (_taxa_ids, error))
            else:
                with open(self.results_file, "a", encoding="utf-8") as ofile:
                    ofile.write("%s### END ###\n" % handle.read())
        return

    def _get_taxa(self, _taxa_ids):
        self._clear_files()
        _taxa_ids = self._split_for_url(_taxa_ids)
        MyFuncs.run_multicore_function(_taxa_ids, self._mc_taxa, [Lock()], max_processes=3, quiet=True)
        with open(self.results_file, "r", encoding="utf-8") as ifile:
            results = ifile.read().split("\n### END ###\n")
            results = [x for x in results if x]

        _output = {}
        for result in results:
            for summary in Entrez.parse(StringIO(result)):
                _output[summary["TaxId"]] = summary["ScientificName"]
        return _output

    def _mc_accn2gi(self, accns, args):
        lock = args[0]
        error = False
        handle = False
        timer = MyFuncs.time()
        for i in range(self.max_attempts):
            try:
                handle = Entrez.efetch(db="nucleotide", id=accns, rettype="gi", retmax=10000)
                '''
                Example output: efetch.fcgi?db=nucleotide&id=XP_010103297.1,XP_010103298.1,XP_010103299.1&rettype=gi
                    703125407
                    703125412
                    703125416
                '''
                timer = MyFuncs.time() - timer
                if timer < 1:
                    sleep(1 - timer)
                break
            except HTTPError as _e:
                if i == self.max_attempts - 1:
                    error = _e
                sleep(1)

        with lock:
            if error:
                with open(self.http_errors_file, "r", encoding="utf-8") as ifile:
                    ifile.write("%s\n%s//\n" % (accns, error))
            else:
                with open(self.results_file, "a", encoding="utf-8") as ofile:
                    ofile.write("%s### END ###\n" % handle.read())
        return

    def _get_gis(self, accns):  # These accns should include version numbers
        self._clear_files()
        accns = self._split_for_url(accns)
        runtime = MyFuncs.RunTime(prefix="\t")
        _stderr("Converting NCBI accessions to gi numbers...\n")
        runtime.start()
        MyFuncs.run_multicore_function(accns, self._mc_accn2gi, [Lock()], max_processes=3, quiet=True)
        runtime.end()
        with open(self.results_file, "r", encoding="utf-8") as ifile:
            results = ifile.read().split("\n### END ###\n")
            results = [x.split("\n") for x in results]
            results = [x for sublist in results for x in sublist if x]
        _stderr("\tDone\n")
        return results

    def _mc_summaries(self, gi_nums, args):
        lock = args[0]
        error = False
        handle = False
        timer = MyFuncs.time()
        for i in range(self.max_attempts):
            try:
                # db needs to be set to something, but if using gi nums it doesn't matter if protein or nucleotide.
                handle = Entrez.esummary(db="nucleotide", id=gi_nums, retmax=10000)
                '''
                Example output: esummary.fcgi?db=nucleotide&id=728840875
                    <eSummaryResult>
                        <DocSum>
                            <Id>728840875</Id>
                            <Item Name="Caption" Type="String">KHG20318</Item>
                            <Item Name="Title" Type="String">
                                Proline-rich receptor-like protein kinase PERK1 [Gossypium arboreum]
                            </Item>
                            <Item Name="Extra" Type="String">
                                gi|728840875|gb|KHG20318.1||gnl|WGS:JRRC|F383_09126[728840875]
                            </Item>
                            <Item Name="Gi" Type="Integer">728840875</Item>
                            <Item Name="CreateDate" Type="String">2014/12/04</Item>
                            <Item Name="UpdateDate" Type="String">2014/12/04</Item>
                            <Item Name="Flags" Type="Integer">0</Item>
                            <Item Name="TaxId" Type="Integer">29729</Item>
                            <Item Name="Length" Type="Integer">649</Item>
                            <Item Name="Status" Type="String">live</Item>
                            <Item Name="ReplacedBy" Type="String"/>
                            <Item Name="Comment" Type="String">
                                <![CDATA[ ]]>
                            </Item>
                        </DocSum>
                    </eSummaryResult>
                '''
                timer = MyFuncs.time() - timer
                if timer < 1:
                    sleep(1 - timer)
                break
            except [HTTPError, RuntimeError] as _e:
                if i == self.max_attempts - 1:
                    error = _e
                sleep(1)

        with lock:
            if error:
                with open(self.http_errors_file, "r", encoding="utf-8") as ifile:
                    ifile.write("%s\n%s//\n" % (gi_nums, error))
            else:
                with open(self.results_file, "a", encoding="utf-8") as ofile:
                    ofile.write("%s### END ###\n" % handle.read())
        return

    def _fetch_summaries(self, gi_nums):
        self._clear_files()
        gi_nums = self._split_for_url(gi_nums)
        runtime = MyFuncs.RunTime(prefix="\t")
        _stderr("Retrieving record summaries from NCBI...\n")
        runtime.start()
        MyFuncs.run_multicore_function(gi_nums, self._mc_summaries, [Lock()], max_processes=3, quiet=True)
        runtime.end()
        with open(self.results_file, "r", encoding="utf-8") as _ifile:
            results = _ifile.read().split("\n### END ###\n")
            results = [x for x in results if x != ""]

        _output = {}
        taxa = []
        for result in results:
            for summary in Entrez.parse(StringIO(result)):
                _rec = OrderedDict()
                _rec["gi_num"] = str(summary["Gi"])
                # status can be 'live', 'dead', 'withdrawn', 'replaced'
                _rec["status"] = summary["Status"] if summary["ReplacedBy"] == '' else \
                    "%s->%s" % (summary["Status"], summary["ReplacedBy"])
                _rec["TaxId"] = summary["TaxId"]
                if summary["TaxId"] not in taxa:
                    taxa.append(summary["TaxId"])
                _rec["organism"] = ""
                _rec["length"] = summary["Length"]
                _rec["comments"] = summary["Title"]

                _output[summary["Caption"]] = Record(summary["Caption"], gi=str(summary["Gi"]),
                                                     summary=_rec, _size=_rec["length"])
                _output[summary["Caption"]].guess_database()
        taxa = self._get_taxa(taxa)
        for accn, rec in _output.items():
            rec.summary["organism"] = taxa[rec.summary["TaxId"]]
        _stderr("\t%s records received.\n" % len(_output))
        return _output

    def fetch_summary(self):
        # EUtils esummary will only take gi numbers
        self._clear_files()
        accns = [rec.ncbi_accn() for accn, rec in self.dbbuddy.records.items() if
                 rec.database in ["ncbi_nuc", "ncbi_prot"] and not rec.gi]

        if accns:
            gi_nums = self._get_gis(accns)
            summaries = self._fetch_summaries(gi_nums)
            for accn in accns:
                if re.search("\.[0-9]+$", accn):  # Add record the version number if present
                    wrong_accn = re.search("^(.*?)\.", accn).group(1)
                    summaries[wrong_accn].accession = accn
                    summaries[accn] = summaries[wrong_accn]
                    del summaries[wrong_accn]

                if accn not in summaries:
                    failure = Failure("ACCN: %s" % accn, "Unable to fetch summary from NCBI")
                    if failure.hash not in self.dbbuddy.failures:
                        self.dbbuddy.failures[failure.hash] = failure

            for accn, rec in summaries.items():
                if accn in self.dbbuddy.records:
                    self.dbbuddy.records[accn].update(rec)
                else:
                    self.dbbuddy.records[accn] = rec

            open(self.results_file, "w", encoding="utf-8").close()

        gi_nums = [accn for accn, rec in self.dbbuddy.records.items() if rec.type == "gi_num" and not rec.summary]

        if gi_nums:
            summaries = self._fetch_summaries(gi_nums)
            summary_gis = [_rec.gi for accn, _rec in summaries.items()]
            for gi_num in gi_nums:
                if gi_num not in summary_gis:
                    failure = Failure("gi: %s" % gi_num, "gi_nums: Unable to fetch summary from NCBI")
                    self.dbbuddy.failures[failure.hash] = failure

            for accn, rec in summaries.items():
                if accn in self.dbbuddy.records:
                    self.dbbuddy.records[accn].update(rec)
                elif rec.gi in self.dbbuddy.records:
                    self.dbbuddy.records[rec.gi].update(rec)
                    self.dbbuddy.records[accn] = self.dbbuddy.records[rec.gi]
                    del self.dbbuddy.records[rec.gi]
                else:
                    self.dbbuddy.records[accn] = rec

    def search_ncbi(self, database):  # database in ["nucleotide", "protein"]
        for _term in self.dbbuddy.search_terms:
            try:
                count = Entrez.read(Entrez.esearch(db=database, term=_term, rettype="count"))["Count"]
                handle = Entrez.esearch(db=database, term=_term, retmax=count)
                '''
                Example output: esearch.fcgi?db=nucleotide&term=perk1&retmax=5
                <eSearchResult>
                    <Count>456</Count>
                    <RetMax>5</RetMax>
                    <RetStart>0</RetStart>
                    <IdList>
                        <Id>909549231</Id>
                        <Id>909549227</Id>
                        <Id>909549224</Id>
                        <Id>909546647</Id>
                        <Id>306819620</Id>
                    </IdList>
                    <TranslationSet/>
                    <TranslationStack>
                        <TermSet>
                            <Term>perk1[All Fields]</Term>
                            <Field>All Fields</Field>
                            <Count>456</Count>
                            <Explode>N</Explode>
                        </TermSet>
                        <OP>GROUP</OP>
                    </TranslationStack>
                    <QueryTranslation>perk1[All Fields]</QueryTranslation>
                </eSearchResult>
                '''
                result = Entrez.read(handle)
                for _id in result["IdList"]:
                    if _id not in self.dbbuddy.records:
                        self.dbbuddy.records[_id] = Record(_id)
                        self.dbbuddy.records[_id].guess_database()

                self.fetch_summary()

            except [HTTPError, RuntimeError] as _e:
                if _e.getcode() == 503:
                    failure = Failure(_term, "503 'Service unavailable': NCBI is either blocking you or they are "
                                             "experiencing some technical issues.")
                    if failure.hash not in self.dbbuddy.failures:
                        self.dbbuddy.failures[failure.hash] = failure
                else:
                    failure = Failure(_term, str(_e))
                    if failure.hash not in self.dbbuddy.failures:
                        self.dbbuddy.failures[failure.hash] = failure

            except KeyboardInterrupt:
                _stderr("\n\tNCBI query interrupted by user\n")

    def _mc_seq(self, accns, args):
        database, lock = args
        error = False
        handle = False
        timer = MyFuncs.time()
        for i in range(self.max_attempts):
            try:
                handle = Entrez.efetch(db=database, id=accns, rettype="gb", retmode="text", retmax=10000)
                '''
                Example output: efetch.fcgi?db=protein&id=920714169&rettype=gb&retmode=text
                LOCUS       KOM54257                 441 aa            linear   PLN 21-AUG-2015
                DEFINITION  hypothetical protein LR48_Vigan10g014900 [Vigna angularis].
                ACCESSION   KOM54257
                VERSION     KOM54257.1  GI:920714169
                DBLINK      BioProject: PRJNA261643
                            BioSample: SAMN03074979
                DBSOURCE    accession CM003380.1
                KEYWORDS    .
                SOURCE      Vigna angularis (adzuki bean)
                  ORGANISM  Vigna angularis
                            Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                            Spermatophyta; Magnoliophyta; eudicotyledons; Gunneridae;
                            Pentapetalae; rosids; fabids; Fabales; Fabaceae; Papilionoideae;
                            Phaseoleae; Vigna.
                REFERENCE   1  (residues 1 to 441)
                  AUTHORS   Wan,P.
                  TITLE     The draft genome sequence of high starch accumulation legume
                            species adzuki bean (Vigna angularis)
                  JOURNAL   Unpublished
                REFERENCE   2  (residues 1 to 441)
                  AUTHORS   Wan,P.
                  TITLE     Direct Submission
                  JOURNAL   Submitted (27-FEB-2015) College of Plant Science and Technology,
                            Beijing University of Agriculture, Huilongguan Beinonglu 7
                            Changping District, Beijing 102206, China
                COMMENT     Method: conceptual translation.
                FEATURES             Location/Qualifiers
                     source          1..441
                                     /organism="Vigna angularis"
                                     /cultivar="Jingnong 6"
                                     /db_xref="taxon:3914"
                                     /chromosome="10"
                                     /tissue_type="seedling"
                                     /country="China: Beijing"
                     Protein         1..441
                                     /product="hypothetical protein"
                     CDS             1..441
                                     /locus_tag="LR48_Vigan10g014900"
                                     /coded_by="complement(join(CM003380.1:1168320..1168634,
                                     CM003380.1:1168784..1168945,CM003380.1:1169144..1169291,
                                     CM003380.1:1169589..1169665,CM003380.1:1169837..1169907,
                                     CM003380.1:1169999..1170085,CM003380.1:1170242..1170707))"
                                     /note="GO_function: GO:0004672 - protein kinase activity
                                     [Evidence IEA];
                                     GO_function: GO:0004674 - protein serine/threonine kinase
                                     activity [Evidence IEA];
                                     GO_function: GO:0004713 - protein tyrosine kinase activity
                                     [Evidence IEA];
                                     GO_function: GO:0005524 - ATP binding [Evidence IEA];
                                     GO_process: GO:0006468 - protein phosphorylation [Evidence
                                     IEA]"
                                     /db_xref="InterPro:IPR000719"
                                     /db_xref="InterPro:IPR001245"
                                     /db_xref="InterPro:IPR002290"
                                     /db_xref="InterPro:IPR008271"
                                     /db_xref="InterPro:IPR017441"
                                     /db_xref="InterPro:IPR020635"
                                     /db_xref="UniProtKB/Swiss-Prot:PERK1"
                ORIGIN
                        1 mppkpspppa payaaqpppp pppfiissgg sgsnysggep lpppspgisl gfskstftye
                       61 elaratdgfs danllgqggf gyvhrgilpn gkevavkqlk agsgqgeref qaeveiisrv
                      121 hhkhlvslvg ycitgsqrll vyefvpnntm efhlhgrgrp tmdwptrlri algsakglay
                      181 lhedchpkii hrdiksanil ldfkfeakva dfglakfssd vnthvstrvm gtfgylapey
                      241 assgkltdks dvfsygvmll elitgrrpvd ktqtfmedsl vdwarplltr aleeddfdsi
                      301 idprlqndyd pnemarmvac aaactrhsak rrprmsqvvr alegdvslad lnegikpghs
                      361 tmysshessd ydtvqyredm kkfrkmalgt qeygasseys aatseyglnp sgssseaqsr
                      421 qttrememrk mknsqgfsgs s
                //
                '''
                timer = MyFuncs.time() - timer
                if timer < 1:
                    sleep(1 - timer)
                break
            except [HTTPError, RuntimeError] as _e:
                if i == self.max_attempts - 1:
                    error = _e
                sleep(1)

        with lock:
            if error:
                with open(self.http_errors_file, "r", encoding="utf-8") as ifile:
                    ifile.write("%s\n%s//\n" % (accns, error))
            else:
                with open(self.results_file, "a", encoding="utf-8") as ofile:
                    ofile.write(handle.read())
        return

    def _get_seq(self, gi_nums, database):
        self._clear_files()
        gi_nums = self._split_for_url(gi_nums)
        runtime = MyFuncs.RunTime(prefix="\t")
        _stderr("Fetching full sequence records from NCBI...\n")
        runtime.start()
        MyFuncs.run_multicore_function(gi_nums, self._mc_seq, [database, Lock()], max_processes=3, quiet=True)
        runtime.end()
        with open(self.results_file, "r", encoding="utf-8") as ifile:
            results = SeqIO.to_dict(SeqIO.parse(ifile, "gb"))
        _stderr("\tDone\n")
        return results

    def fetch_sequence(self, database):  # database in ["nucleotide", "protein"]
        db = "ncbi_nuc" if database == "nucleotide" else "ncbi_prot"
        gi_nums = [_rec.gi for accn, _rec in self.dbbuddy.records.items() if _rec.database == db]
        if len(gi_nums) > 0:
            try:
                records = self._get_seq(gi_nums, database)
                for accn, rec in records.items():
                    # Catch accn/version
                    with_version = re.search("^(.*?)\.([0-9]+)$", accn)
                    if with_version:
                        accn, ver = with_version.group(1), with_version.group(2)
                        self.dbbuddy.records[accn].record = rec
                        self.dbbuddy.records[accn].version = ver
                    else:
                        try:
                            self.dbbuddy.records[accn].record = rec
                        except KeyError:
                            self.dbbuddy.failures.setdefault("# NCBI fetch: Ids not in dbbuddy.records",
                                                             []).append(rec.id)
            except KeyboardInterrupt:
                _stderr("\n\tNCBI query interrupted by user\n")


class EnsemblRestClient(object):
    def __init__(self, _dbbuddy, server='http://rest.ensembl.org/'):
        self.dbbuddy = _dbbuddy
        self.temp_dir = MyFuncs.TempDir()
        self.http_errors_file = "%s/errors.txt" % self.temp_dir.path
        open(self.http_errors_file, "w", encoding="utf-8").close()
        self.results_file = "%s/results.txt" % self.temp_dir.path
        open(self.results_file, "w", encoding="utf-8").close()
        self.server = server
        self.species = self.perform_rest_action("info/species", headers={"Content-type": "application/json",
                                                                         "Accept": "application/json"})["species"]
        self.species = {x["display_name"]: x for x in self.species if x["display_name"]}

    def perform_rest_action(self, endpoint, **kwargs):
        """
        :param endpoint:
        :param kwargs: requires 'headers' {'Content-type': [text/x-seqxml+xml, application/json],
                                           "Accept": "application/json"} and can also take 'data'
        :return:
        """
        endpoint = endpoint.strip("/")
        kwargs_backup = dict(kwargs)
        try:
            if "data" in kwargs:
                data = '{'
                for key, value in kwargs["data"].items():
                    data += '"%s": %s, ' % (key, value)
                data = "%s}" % data.strip(", ")
                data = re.sub("'", '"', data)
                kwargs["data"] = data.encode('utf-8')

            request = Request(self.server + endpoint, **kwargs)
            response = urlopen(request)

            if request.get_header("Content-type") == "application/json":
                content = response.read().decode()
                data = json.loads(content)
            elif request.get_header("Content-type") == "text/x-seqxml+xml":
                data = SeqIO.parse(response, "seqxml")
            else:
                raise ValueError(request.headers)

            return data

        except HTTPError as _e:
            # check if we are being rate limited by the server
            if _e.getcode() == 429:
                if 'Retry-After' in _e.headers:
                    retry = _e.headers['Retry-After']
                    sleep(float(retry) + 1)
                    self.perform_rest_action(endpoint, **kwargs_backup)
            else:
                failure = Failure("%s" % self.server + endpoint, "Ensemble request failed. %s" % _e)
                self.dbbuddy.failures[failure.hash] = failure

    def fetch_nucleotide(self):
        accns = [accn for accn, rec in self.dbbuddy.records.items() if rec.database == "ensembl"]
        if len(accns) > 0:
            _stderr("Fetching sequence from Ensembl...\n")
            runtime = MyFuncs.RunTime(prefix="\t")
            runtime.start()
            data = self.perform_rest_action("sequence/id", data={"ids": accns},
                                            headers={"Content-type": "text/x-seqxml+xml"})
            runtime.end()
            for rec in data:
                summary = self.dbbuddy.records[rec.id].summary
                rec.description = summary['comments']
                rec.accession = rec.id
                rec.name = rec.id
                species = re.search("([a-z])[a-z]*_([a-z]{1,3})", summary['organism'])
                species = "%s%s" % (species.group(1).upper(), species.group(2))
                new_id = "%s-%s" % (species, summary['name'])
                self.dbbuddy.records[rec.id].record = rec
                self.dbbuddy.records[rec.id].record.id = new_id

    def _mc_search(self, species, args):
        identifier, lock = args
        self.dbbuddy.failures = {}
        data = self.perform_rest_action("lookup/symbol/%s/%s" % (species, identifier),
                                        headers={"Content-type": "application/json", "Accept": "application/json"})
        with lock:
            with open(self.results_file, "a", encoding="utf-8") as ofile:
                ofile.write("%s\n### END ###\n" % data)
            if self.dbbuddy.failures:
                with open(self.http_errors_file, "a", encoding="utf-8") as ofile:
                    for _hash, failure in self.dbbuddy.failures.items():
                        ofile.write("%s\n" % failure)

    def _parse_summary(self, summary):
        accn = summary['id']
        size = abs(summary["start"] - summary["end"])
        _version = summary['version']

        required_keys = ['display_name', 'species', 'biotype', 'object_type',
                         'strand', 'assembly_name', 'description', 'version']

        for key in required_keys:
            if key not in summary:
                summary[key] = ''

        summary = OrderedDict([('name', summary['display_name']), ('length', size),
                               ('organism', summary['species']),
                               ('organism-id', self.species[summary['species']]['taxon_id']),
                               ('biotype', summary['biotype']), ('object_type', summary['object_type']),
                               ('strand', summary['strand']), ('assembly_name', summary['assembly_name']),
                               ('comments', summary['description'])])

        rec = Record(accn, summary=summary, _version=_version,
                     _size=size, _database="ensembl", _type="nucleotide")
        return rec

    def search_ensembl(self):
        open(self.http_errors_file, "w", encoding="utf-8").close()
        open(self.results_file, "w", encoding="utf-8").close()
        species = [_name for _name, _info in self.species.items()]
        for search_term in self.dbbuddy.search_terms:
            _stderr("Searching Ensembl for %s...\n" % search_term)
            runtime = MyFuncs.RunTime(prefix="\t")
            runtime.start()
            MyFuncs.run_multicore_function(species, self._mc_search, [search_term, Lock()], quiet=True)
            runtime.end()
            with open(self.results_file, "r", encoding="utf-8") as ifile:
                results = ifile.read().split("\n### END ###")

            counter = 0
            for rec in results:
                rec = rec.strip()
                if not rec or rec in ["None", ""]:
                    continue
                counter += 1
                rec = re.sub("'", '"', rec)
                rec = self._parse_summary(json.loads(rec))
                if rec.accession in self.dbbuddy.records:
                    self.dbbuddy.records[rec.accession].update(rec)
                else:
                    self.dbbuddy.records[rec.accession] = rec
            _stderr("Returned %s results\n\n" % counter)

    def fetch_summary(self):
        accns = [accn for accn, rec in self.dbbuddy.records.items() if rec.database == "ensembl"]
        data = self.perform_rest_action("lookup/id", data={"ids": accns},
                                        headers={"Content-type": "application/json", "Accept": "application/json"})

        for accn, results in data.items():
            size = abs(results["start"] - results["end"])
            required_keys = ['display_name', 'species', 'biotype', 'object_type',
                             'strand', 'assembly_name', 'description']
            for key in required_keys:
                if key not in results:
                    results[key] = ''

            summary = OrderedDict([('name', results['display_name']), ('length', size),
                                   ('organism', results['species']), ('biotype', results['biotype']),
                                   ('object_type', results['object_type']), ('strand', results['strand']),
                                   ('assembly_name', results['assembly_name']), ('comments', results['description'])])

            rec = Record(accn, summary=summary, _version=results['version'],
                         _size=size, _database="ensembl", _type="nucleotide")
            self.dbbuddy.records[accn].update(rec)
        return


# ################################################ MAIN API FUNCTIONS ################################################ #
class LiveSearch(cmd.Cmd):
    def __init__(self, _dbbuddy, crash_file):
        """
        :param _dbbuddy: pre-instantiated DbBuddy object
        :param crash_file: MyFuncs.TempFile object instantiated in binary mode
        """
        self.terminal_default = "\033[m\033[40m%s" % WHITE
        cmd.Cmd.__init__(self)
        hash_heading = ""
        colors = terminal_colors()
        for _ in range(23):
            hash_heading += "%s#" % next(colors)
        _stdout('''{1}

{0} {1}{3}{2}Welcome to the DatabaseBuddy live shell{1} {0}{1}

{2}Type 'help' for a list of available commands or 'help <command>' for further details.
                  To end the session, use the 'quit' command.{1}

'''.format(hash_heading, self.terminal_default, BOLD, UNDERLINE))
        self.prompt = '{0}{1}DbBuddy>{0} '.format(self.terminal_default, BOLD)

        self.doc_leader = '''\

{0}{1}      {2}{3}DatabaseBuddy Help{1}{4}      {0}{1}

A general workflow: 1) {5}search{1} databases with search terms or accession numbers
                    2) {5}show{1} summary information (no sequence has been downloaded yet)
                    3) Create filtered set of results with {5}keep{1} and {5}remove{1}
                    4) {5}restore{1} some records from the trash bin to the filtered set
                    5) {5}fetch{1} full sequence records for the filtered set
                    6) Switch to a {5}format{1} that includes sequences, like fasta or genbank
                    7) {5}write{1} sequences to a file or {5}save{1} the live session

Further details about each command can be accessed by typing 'help <command>'
'''.format("".join(["%s-" % next(colors) for _ in range(29)]), self.terminal_default,
           UNDERLINE, BOLD, NO_UNDERLINE, GREEN)
        self.doc_header = "Available commands:                                                         "
        self.dbbuddy = _dbbuddy
        self.crash_file = crash_file
        self.dump_session()

        if CONFIG["install_path"]:
            self.history_path = "%s/.cmd_history" % CONFIG["install_path"]

        else:
            self.history_path = "/tmp/.db_cmd_history"

        if not os.path.isfile(self.history_path):
            open(self.history_path, "w", encoding="utf-8").close()
        readline.read_history_file(self.history_path)

        # As implemented, one UnDo is possible (reload the most recent dump). Set self.undo to true every time a dump
        # occurs, and back to False if undo is used.
        self.undo = False

        _stderr(self.terminal_default)  # This needs to be called here if stderr is going to format correctly
        if self.dbbuddy.records or self.dbbuddy.search_terms:
            retrieve_summary(_dbbuddy)
        else:
            _stdout("Your session is currently unpopulated. Use 'search' to retrieve records.\n",
                    format_out=self.terminal_default)
        self.hash = None
        self.shell_execs = []  # Only populate this if "bash" is called by the user
        self.usage = br.Usage()
        self.usage.increment("DatabaseBuddy", VERSION.short(), "LiveSearch")
        self.cmdloop()

    # @staticmethod
    # def do_crash(line=None):
    #    open("a file that doesn't exist")

    def precmd(self, line):
        readline.write_history_file(self.history_path)
        return line

    def postcmd(self, stop, line):
        command = line.split(" ")[0]
        self.usage.increment("DatabaseBuddy", VERSION.short(), command)
        return stop

    def dump_session(self):
        import pickle
        self.crash_file.save("%s_undo" % self.crash_file.path)
        self.crash_file.open()
        pickle.dump(self.dbbuddy, self.crash_file.handle, protocol=-1)
        self.crash_file.close()
        self.undo = True

    def default(self, line):
        _stdout('*** Unknown syntax: %s\n\n' % line, format_in=RED, format_out=self.terminal_default)

    @staticmethod
    def _append_slash_if_dir(p):  # Used for expanding file patsh
            if p and os.path.isdir(p) and p[-1] != os.sep:
                return p + os.sep
            else:
                return p

    def get_headings(self):  # ToDo: This does not get all possible headings. Fix that...
        headings = []
        if len(self.dbbuddy.records) > 0:
            _rec = []
            for _accn, _rec in self.dbbuddy.records.items():
                break
            headings = ["ACCN", "DB"] + [heading for heading, _value in _rec.summary.items()]
        return headings

    def filter(self, line, mode="keep"):
        if mode not in ["keep", "remove", "restore"]:
            raise ValueError("The 'mode' argument in filter() must be "
                             "'keep', 'remove', or 'restore', not %s." % mode)

        if not line:
            if mode == "keep":
                action = "Specify a search string to be used as a filter (records will be retained): "
            elif mode == "remove":
                action = "Specify a search string to be used as a filter (records will be removed): "
            else:
                action = "Specify a string to search the trash bin with: "
            line = input("%s%s%s" %
                         (RED, action, self.terminal_default))

        # If the user doesn't supply a string, do nothing
        if not line:
            _stdout("Error: you must specify a search string.\n", format_in=RED, format_out=self.terminal_default)
            return

        # Kill the command if the user is mixing quote types to separate search terms
        error_message = "Error: It appears that you are trying to mix quote types (\" and ') while specifying " \
                        "multiple filters. Please pick one or the other.\n\n"
        if line[0] == "'":
            if line.strip()[-1] == '"':
                _stdout(error_message, format_in=RED, format_out=self.terminal_default)
                return
            line = line.strip("'").split("' '")
        else:
            if line.strip()[-1] == "'":
                _stdout(error_message, format_in=RED, format_out=self.terminal_default)
                return
            line = line.strip('"').split('" "')
        max_regex_len = 6  # length of the string 'filter'
        for _filter in line:
            max_regex_len = len(_filter) if len(_filter) > max_regex_len else max_regex_len
        tabbed = "{0: <%s}{1}\n" % (max_regex_len + 2)
        _stdout("\033[4m%s\n" % (" " * (max_regex_len + 16)), format_in=RED, format_out=self.terminal_default)
        heading = "# Recs recovered" if mode == "restore" else "# Recs removed"
        _stdout(tabbed.format("Filter", heading), format_out=self.terminal_default)

        _errors = {"KeyError": [], "ValueError": []}
        current_count = len(self.dbbuddy.records)
        for _filter in line:
            for _key, _value in self.dbbuddy.filter_records(_filter, mode=mode).items():
                _errors[_key] += _value
            _stdout(tabbed.format(_filter, abs(current_count - len(self.dbbuddy.records))),
                    format_out=self.terminal_default)
            current_count = len(self.dbbuddy.records)

        if _errors["KeyError"]:
            _stderr("%s\nThe following column headings were not present in all records (ignored):\n"
                    "%s%s\n" % (RED, ", ".join(_errors["KeyError"]), DEF_FONT))

        if _errors["ValueError"]:
            _stderr("%s\nThe following errors occurred:\n"
                    "%s%s\n" % (RED, ", ".join(_errors["ValueError"]), DEF_FONT))

        output_message = "\n%s records remain.\n\n" % len(self.dbbuddy.records) if mode != "restore" \
            else "\n%s records remain in the trash bin.\n\n" % len(self.dbbuddy.trash_bin)

        _stdout(output_message, format_in=GREEN, format_out=self.terminal_default)
        self.dump_session()

    def do_bash(self, line):
        _stdout("", format_out=CYAN)
        if not line:
            line = input("Bash> ")
        # Need to strip out leading/trailing quotes for this to work
        line = re.sub('^["](.*)["]$', r"\1", line)
        line = re.sub("^['](.*)[']$", r"\1", line)
        if line[:2] == "cd":
            line = line.lstrip("cd ")
            try:
                _path = os.path.abspath(line)
                os.chdir(_path)
                _stdout("%s\n" % _path, format_out=CYAN)
            except FileNotFoundError:
                _stdout("-sh: cd: %s: No such file or directory\n" % line, format_in=RED,
                        format_out=self.terminal_default)
        else:
            Popen(line, shell=True).wait()
        _stdout("\n", format_out=self.terminal_default)

    def do_database(self, line):
        if not line:
            line = input("%sSpecify database:%s " % (RED, self.terminal_default))

        line = line.split(" ")
        new_database_list = []
        for l in line:
            if l not in DATABASES and l != "all":
                _stdout("Error: %s is not a valid database choice.\n"
                        "Please select from %s\n" % (l, ["all"] + DATABASES), format_in=RED,
                        format_out=self.terminal_default)
            else:
                new_database_list.append(l)
        if new_database_list:
            if "all" in new_database_list:
                self.dbbuddy.databases = DATABASES
            else:
                self.dbbuddy.databases = new_database_list

            _stdout("Database search list updated to %s\n\n" % self.dbbuddy.databases, format_in=GREEN,
                    format_out=self.terminal_default)
        else:
            _stdout("Database search list not changed.\n\n", format_in=RED, format_out=self.terminal_default)
        self.dump_session()

    def do_delete(self, line="all"):
        if not self.dbbuddy.trash_bin and not self.dbbuddy.records and not self.dbbuddy.search_terms:
            _stdout("The live session is already empty.\n\n", format_in=RED, format_out=self.terminal_default)
            return

        line = line.lower()
        if line not in ["", "a", "all", "failures", "f"] + TRASH_SYNOS + RECORD_SYNOS + SEARCH_SYNOS:
            _stdout("Sorry, I don't understand what you want to delete.\n Select from: all, main, trash-bin, "
                    "failures, search\n\n", format_in=RED, format_out=self.terminal_default)
            return

        if line in ["failures", "f"]:
            if not self.dbbuddy.failures:
                _stdout("Failures list is already empty.\n\n", format_in=RED, format_out=self.terminal_default)
            else:
                confirm = input("%sAre you sure you want to clear all %s failures (y/[n])?%s " %
                                (RED, len(self.dbbuddy.failures), self.terminal_default))

                if confirm.lower() not in ["yes", "y"]:
                    _stdout("Aborted...\n", format_in=RED, format_out=self.terminal_default)
                else:
                    self.dbbuddy.failures = {}
                    _stdout("List of failures removed.\n\n", format_in=GREEN, format_out=self.terminal_default)

        elif line in SEARCH_SYNOS:
            if not self.dbbuddy.search_terms:
                _stdout("Search terms list is already empty.\n\n", format_in=RED, format_out=self.terminal_default)
            else:
                confirm = input("%sAre you sure you want to delete all %s search terms (y/[n])?%s " %
                                (RED, len(self.dbbuddy.search_terms), self.terminal_default))

                if confirm.lower() not in ["yes", "y"]:
                    _stdout("Aborted...\n", format_in=RED, format_out=self.terminal_default)
                else:
                    self.dbbuddy.search_terms = []
                    _stdout("Search terms removed.\n\n", format_in=GREEN, format_out=self.terminal_default)

        elif line in TRASH_SYNOS:
            if not self.dbbuddy.trash_bin:
                _stdout("Trash bin is already empty.\n", format_in=RED, format_out=self.terminal_default)
            else:
                confirm = input("%sAre you sure you want to delete all %s records from your trash bin (y/[n])?%s " %
                                (RED, len(self.dbbuddy.trash_bin), self.terminal_default))

                if confirm.lower() not in ["yes", "y"]:
                    _stdout("Aborted...\n", format_in=RED, format_out=self.terminal_default)
                else:
                    self.dbbuddy.trash_bin = {}
                    _stdout("Trash bin emptied.\n\n", format_in=GREEN, format_out=self.terminal_default)

        elif line in RECORD_SYNOS:
            if not self.dbbuddy.records:
                _stdout("Records list is already empty.\n", format_in=RED, format_out=self.terminal_default)
            else:
                confirm = input("%sAre you sure you want to delete all %s records from your main "
                                "filtered list (y/[n])?%s " % (RED, len(self.dbbuddy.records), self.terminal_default))
                if confirm.lower() not in ["yes", "y"]:
                    _stdout("Aborted...\n", format_in=RED, format_out=self.terminal_default)
                else:
                    self.dbbuddy.records = {}
                    _stdout("All records removed from main list (trash bin is still intact).\n\n",
                            format_in=GREEN, format_out=self.terminal_default)

        else:
            confirm = input("%sAre you sure you want to completely reset your live session (y/[n])?%s " %
                            (RED, self.terminal_default))

            if confirm.lower() not in ["yes", "y"]:
                _stdout("Aborted...\n", format_in=RED, format_out=self.terminal_default)
            else:
                self.dbbuddy.trash_bin = {}
                self.dbbuddy.records = OrderedDict()
                self.dbbuddy.search_terms = []
                self.dbbuddy.failures = {}
                _stdout("Live session cleared of all data.\n\n", format_in=GREEN, format_out=self.terminal_default)

        self.dump_session()

    def do_failures(self, line=None):
        if line != "":
            _stdout("Note: 'failures' does not take any arguments\n", format_in=RED, format_out=self.terminal_default)

        if not self.dbbuddy.failures:
            _stdout("No failures to report\n\n", format_in=GREEN, format_out=self.terminal_default)
        else:
            _stdout("The following failures have occured\n", format_in=[UNDERLINE, GREEN],
                    format_out=self.terminal_default)
            for _hash, _values in self.dbbuddy.failures.items():
                _stdout("%s\n\n" % _values, format_out=self.terminal_default)

    def do_fetch(self, line=None):
        if line != "":
            _stdout("Note: 'fetch' does not take any arguments\n", format_in=RED, format_out=self.terminal_default)

        accn_only = self.dbbuddy.record_breakdown()["accession"]
        if accn_only:
            search_terms = list(self.dbbuddy.search_terms)
            self.dbbuddy.search_terms = []
            retrieve_summary(self.dbbuddy)
            self.dbbuddy.search_terms = search_terms

        amount_seq_requested = 0
        new_records_fetched = []
        for _accn, _rec in self.dbbuddy.records.items():
            if not _rec.record:  # Not fetching sequence if the full record already exists
                amount_seq_requested += _rec.size
                new_records_fetched.append(_accn)

        if amount_seq_requested > 5000000:
            confirm = input("{0}You are requesting {2}{1}{0} residues of sequence data. "
                            "Continue (y/[n])?{3}".format(GREEN, MyFuncs.pretty_number(amount_seq_requested),
                                                          YELLOW, self.terminal_default))
            if confirm.lower() not in ["yes", "y"]:
                _stdout("Aborted...\n\n", format_in=RED, format_out=self.terminal_default)
                return

        retrieve_sequences(self.dbbuddy)
        seq_retrieved = 0
        for _accn in new_records_fetched:
            if self.dbbuddy.records[_accn].record:
                seq_retrieved += self.dbbuddy.records[_accn].size

        _stdout("Retrieved %s residues of sequence data\n\n" % MyFuncs.pretty_number(seq_retrieved),
                format_out=self.terminal_default)
        self.dump_session()

    def do_format(self, line):
        if not line:
            line = input("%sSpecify format:%s " % (GREEN, self.terminal_default))

        if line not in FORMATS:
            _stdout("Sorry, {1}'{2}'{0} is not a valid format. Please select from the "
                    "following BioPython supported formats:\n\t{3}\n\n".format(RED, YELLOW, line, ", ".join(FORMATS)),
                    format_in=RED, format_out=self.terminal_default)
            return

        self.dbbuddy.out_format = line
        _stdout("Output format changed to %s%s\n\n" % (YELLOW, line), format_in=GREEN,
                format_out=self.terminal_default)
        self.dump_session()

    def do_load(self, line=None, quiet=False):
        if not line:
            line = input("%sWhere is the dump_file?%s " % (RED, self.terminal_default))
        try:
            with open(os.path.abspath(line), "rb") as ifile:
                self.dbbuddy = pickle.load(ifile)
            self.dump_session()
            for _db, client in self.dbbuddy.server_clients.items():
                if client:
                    client.temp_dir = MyFuncs.TempDir()
                    client.http_errors_file = "%s/errors.txt" % client.temp_dir.path
                    open(client.http_errors_file, "w", encoding="utf-8").close()
                    client.results_file = "%s/results.txt" % client.temp_dir.path
                    open(client.results_file, "w", encoding="utf-8").close()

            _stdout("Session loaded from file.\n\n", format_in=GREEN, format_out=self.terminal_default, quiet=quiet)

        except IOError as _e:
            _stderr("%s\n" % _e)
            _prompt = input("Specify a path to read your session from, or 'abort' to cancel. ")
            if _prompt == 'abort':
                _stderr("Aborted...\n")
            else:
                self.do_load(_prompt)
        except EOFError:
            _stdout("Error: Unable to read the provided file. Are you sure it's a saved DbBuddy live session?\n\n",
                    format_in=RED, format_out=self.terminal_default)

    def do_keep(self, line=None):
        self.filter(line, mode="keep")

    def do_quit(self, line=None):
        if line != "":
            _stdout("Note: 'quit' does not take any arguments\n", format_in=RED, format_out=self.terminal_default)

        if (self.dbbuddy.records or self.dbbuddy.trash_bin) and self.hash != hash(self.dbbuddy):
            confirm = input("You have unsaved records, are you sure you want to quit (y/[n])?")
            if confirm.lower() in ["yes", "y"]:
                pass
            else:
                _stdout("Aborted...\n\n", format_in=RED, format_out=self.terminal_default)
                return
        self.usage.save()
        _stdout("Goodbye\033[m\n\n")
        sys.exit()

    def do_trash(self, line=None):
        self.do_show(line, "trash_bin")

    def do_remove(self, line=None):
        self.filter(line, mode="remove")

    def do_restore(self, line):
        self.filter(line, "restore")

    def do_save(self, line=None):
        if not line:
            line = input("%sWhere would you like your session saved?%s " % (RED, self.terminal_default))

        # Ensure the specified directory exists
        line = os.path.abspath(line)
        _dir = "/%s" % "/".join(line.split("/")[:-1])
        if not os.path.isdir(_dir):
            _stdout("Error: The specified directory does not exist. Please create it before continuing "
                    "(you can use the 'bash' command from within the DbBuddy Live Session.\n\n", format_in=RED,
                    format_out=self.terminal_default)
            return

        # Set the .db extension
        if not re.search("\.db$", line):
            line += ".db"

        # Warn if file exists
        if os.path.isfile(line):
            confirm = input("%sFile already exists, overwrite [y]/n?%s " % (RED, self.terminal_default))
            if confirm.lower() in ["n", "no"]:
                _stdout("Abort...\n\n", format_in=RED, format_out=self.terminal_default)
                return
        try:
            open(line, "wb").close()
        except PermissionError:
            _stdout("Error: You do not have write privileges in the specified directory.\n\n",
                    format_in=RED, format_out=self.terminal_default)
            return

        self.crash_file.save(line)
        self.hash = hash(self.dbbuddy)
        _stdout("Live session saved\n\n", format_in=GREEN, format_out=self.terminal_default)

    def do_search(self, line):
        if not line:
            line = input("%sSpecify search string:%s " % (RED, self.terminal_default))

        temp_buddy = DbBuddy(line)
        temp_buddy.databases = self.dbbuddy.databases

        if len(temp_buddy.records):
            retrieve_sequences(temp_buddy)

        if len(temp_buddy.search_terms):
            retrieve_summary(temp_buddy)

        for _term in temp_buddy.search_terms:
            if _term not in self.dbbuddy.search_terms:
                self.dbbuddy.search_terms.append(line)

        for _accn, _rec in temp_buddy.records.items():
            if _accn not in self.dbbuddy.records:
                self.dbbuddy.records[_accn] = _rec

        for _hash, failure in temp_buddy.failures.items():
            if _hash not in self.dbbuddy.failures:
                self.dbbuddy.failures[_hash] = failure
        self.dump_session()

    def do_show(self, line=None, group="records"):
        if line:
            line = line.split(" ")

        breakdown = self.dbbuddy.trash_breakdown() if group == "trash_bin" else self.dbbuddy.record_breakdown()

        num_records = len(self.dbbuddy.trash_bin) if group == "trash_bin" else len(self.dbbuddy.records)
        if not num_records:
            _stdout("Nothing in '%s' to show.\n\n" % re.sub("_", " ", group), format_in=RED,
                    format_out=self.terminal_default)
            return

        if self.dbbuddy.out_format not in ["ids", "accessions", "summary", "full-summary"]:
            if not breakdown["full"]:
                _stdout("Warning: only summary data available; there is nothing to display in %s format. "
                        "Use 'fetch' to retrieve sequences first.\n\n"
                        % self.dbbuddy.out_format, format_in=RED,
                        format_out=self.terminal_default)
                return

            if breakdown["summary"] or breakdown["accession"]:
                _stderr("%sWarning: %s records are only summary data, so will not be displayed in %s format. "
                        "Use 'fetch' to retrieve all sequence data.%s\n"
                        % (RED, len(breakdown["summary"] + breakdown["accession"]),
                           self.dbbuddy.out_format, self.terminal_default))

            num_records = len(breakdown["full"])

        columns = []
        for _next in line:
            try:
                num_records = int(_next)
            except ValueError:
                columns.append(_next)

        columns = None if not columns else columns

        if num_records > 100:
            confirm = input("%sShow all %s records (y/[n])?%s " %
                            (RED, num_records, self.terminal_default))
            if confirm.lower() not in ["yes", "y"]:
                _stdout("Include an integer value with 'show' to return a specific number of records.\n\n",
                        format_out=self.terminal_default)
                return
        try:
            self.dbbuddy.print(_num=num_records, columns=columns, group=group)
        except ValueError as _e:
            if re.search("Sequences must all be the same length", str(_e)):
                _stdout("Error: BioPython will not output sequences of different length in '%s' format." %
                        self.dbbuddy.out_format, format_in=RED, format_out=self.terminal_default)
            elif re.search("No suitable quality scores found in letter_annotations of SeqRecord", str(_e)):
                _stdout("Error: BioPython requires quality scores to output in '%s' format, and this data is not "
                        "currently available to DatabaseBuddy." % self.dbbuddy.out_format,
                        format_in=RED, format_out=self.terminal_default)
            else:
                raise ValueError(_e)

        _stderr("%s\n" % self.terminal_default)

    def do_sort(self, line=None):
        def sub_sort(records, headings, _rev=False):
            heading = headings[0]
            subgroups = {}
            if heading == "ACCN":
                return OrderedDict(sorted(records.items(), key=lambda _x: _x[1].accession, reverse=_rev))

            int_headings = 0
            for accn, _rec in records.items():
                if heading == "DB":
                    subgroups.setdefault(_rec.database, {})
                    subgroups[_rec.database][accn] = _rec
                elif heading == "record":
                    if _rec.record:
                        _value = "full"
                    else:
                        _value = "summary"
                    subgroups.setdefault(_value, {})
                    subgroups[_value][accn] = _rec
                else:
                    if heading not in _rec.summary:
                        subgroups.setdefault("zzzzz", {})
                        subgroups["zzzzz"][accn] = _rec
                    else:
                        subgroups.setdefault(_rec.summary[heading], {})
                        subgroups[_rec.summary[heading]][accn] = _rec
                        if isinstance(_rec.summary[heading], int):
                            if _rec.summary[heading] > int_headings:
                                int_headings = _rec.summary[heading]

            if int_headings and "zzzzz" in subgroups:
                subgroups[int_headings + 1] = subgroups["zzzzz"]
                del subgroups["zzzzz"]

            subgroups = OrderedDict(sorted(subgroups.items(), key=lambda _x: _x[0], reverse=rev))
            final_order = OrderedDict()
            for subgroup, _recs in subgroups.items():
                if len(_recs) == 1:
                    _recs = [_rec for accn, _rec in _recs.items()]
                    final_order[_recs[0].accession] = _recs[0]
                else:
                    if len(headings) > 1:
                        _recs = sub_sort(_recs, headings[1:], rev)

                    for x, y in _recs.items():
                        final_order[x] = y
            return final_order

        sort_columns = line.split(" ")
        lower_cols = [col.lower() for col in sort_columns]
        rev = True if "rev" in lower_cols or "reverse" in lower_cols else False
        if rev:
            if "rev" in lower_cols:
                rev_indx = lower_cols.index("rev")
            else:
                rev_indx = lower_cols.index("reverse")
            del sort_columns[rev_indx]

        if not sort_columns or sort_columns[0] == '':
            sort_columns = ["ACCN"]

        self.dbbuddy.records = sub_sort(self.dbbuddy.records, sort_columns, rev)
        self.dump_session()

    def do_status(self, line=None):
        if line != "":
            _stdout("Note: 'status' does not take any arguments\n\n", format_in=RED, format_out=self.terminal_default)
        _stdout("%s\n" % str(self.dbbuddy), format_out=self.terminal_default)

    def do_write(self, line=None):
        if not line:
            line = input("%sWhere would you like your records written?%s " % (RED, self.terminal_default))

        # Ensure the specified directory exists
        line = os.path.abspath(line)
        _dir = "/%s" % "/".join(line.split("/")[:-1])
        if not os.path.isdir(_dir):
            _stdout("Error: The specified directory does not exist. Please create it before continuing "
                    "(you can use the 'bash' command from within the DbBuddy Live Session.\n\n", format_in=RED,
                    format_out=self.terminal_default)
            return

        # Warn if file exists
        if os.path.isfile(line):
            confirm = input("%sFile already exists, overwrite [y]/n?%s " % (RED, self.terminal_default))
            if confirm.lower() in ["n", "no"]:
                _stdout("Abort...\n\n", format_in=RED, format_out=self.terminal_default)
                return

        try:
            ofile = open(line, "w", encoding="utf-8")
        except PermissionError:
            _stdout("Error: You do not have write privileges in the specified directory.\n\n",
                    format_in=RED, format_out=self.terminal_default)
            return

        self.dbbuddy.print(quiet=True, destination=ofile)
        breakdown = self.dbbuddy.record_breakdown()
        if self.dbbuddy.out_format in ["ids", "accessions"]:
            _stdout("%s accessions " % len(breakdown["accession"]), format_in=GREEN,
                    format_out=self.terminal_default)
        elif self.dbbuddy.out_format in ["summary", "full-summary"]:
            _stdout("%s summary records " % (len(breakdown["full"] + breakdown["summary"])), format_in=GREEN,
                    format_out=self.terminal_default)
        else:
            non_full = len(breakdown["summary"] + breakdown["accession"])
            if non_full > 0:
                _stdout('''\
NOTE: There are %s summary records in the Live Session, and only full records can be written
  in '%s' format. Use the 'download' command to retrieve full records.
''' % (non_full, self.dbbuddy.out_format), format_in=RED, format_out=self.terminal_default)
            _stdout("%s %s records  " % (len(breakdown["full"]), self.dbbuddy.out_format), format_in=GREEN,
                    format_out=self.terminal_default)
        _stdout("written to %s.\n\n" % line, format_in=GREEN,
                format_out=self.terminal_default)
        self.hash = hash(self.dbbuddy)
        ofile.close()
        _stdout("Records written to file\n\n", format_in=GREEN, format_out=self.terminal_default)

    def do_undo(self, line=None):
        if line != "":
            _stdout("Note: 'status' does not take any arguments\n", format_in=RED, format_out=self.terminal_default)
        if not self.undo:
            _stdout("There is currently no undo history (only a single undo is possible).\n\n",
                    format_in=RED, format_out=self.terminal_default)
            return

        self.do_load("%s_undo" % self.crash_file.path, quiet=True)
        self.undo = False
        _stdout("Most recent state reloaded\n\n", format_in=GREEN, format_out=self.terminal_default)

    def complete_bash(self, *args):
        text = args[0]
        if not self.shell_execs:
            path_dirs = Popen("echo $PATH", stdout=PIPE, shell=True).communicate()
            path_dirs = path_dirs[0].decode("utf-8").split(":")
            for _dir in path_dirs:
                if not os.path.isdir(_dir):
                    continue
                root, dirs, files = next(os.walk(_dir))
                for _file in files:
                    if os.access("%s/%s" % (root, _file), os.X_OK):
                        self.shell_execs.append(_file.strip())
        return ["%s " % x for x in self.shell_execs if x.startswith(text)]

    @staticmethod
    def complete_database(*args):
        text = args[0]
        startidx = args[2]
        endidx = args[3]
        if startidx and endidx:
            pass
        return [db for db in DATABASES if db.startswith(text)]

    @staticmethod
    def complete_delete(*args):
        text = args[0]
        return [x for x in ["all", "failures", "search", "trash", "records"] if x.startswith(text)]

    @staticmethod
    def complete_format(*args):
        text = args[0]
        return [x for x in FORMATS if x.startswith(text)]

    def complete_keep(self, *args):
        text = args[0]
        return ["(%s) " % x for x in self.get_headings() if x.lower().startswith(text.lower())]

    def complete_load(self, *args):
        line, startidx, endidx = args[1:]
        # ToDo: pulled code from stack overflow, modify or credit.
        import glob

        before_arg = line.rfind(" ", 0, startidx)
        if before_arg == -1:
            return  # arg not found

        fixed = line[before_arg + 1:startidx]  # fixed portion of the arg
        arg = line[before_arg + 1:endidx]
        pattern = arg + '*'

        completions = []
        for path in glob.glob(pattern):
            path = self._append_slash_if_dir(path)
            completions.append(path.replace(fixed, "", 1))
        return completions

    def complete_trash(self, *args):
        text = args[0]
        return ["%s " % x for x in self.get_headings() if x.lower().startswith(text.lower())]

    def complete_remove(self, *args):
        text = args[0]
        return ["(%s) " % x for x in self.get_headings() if x.lower().startswith(text.lower())]

    def complete_restore(self, *args):
        text = args[0]
        return ["(%s) " % x for x in self.get_headings() if x.lower().startswith(text.lower())]

    def complete_save(self, *args):
        line, startidx, endidx = args[1:]
        # ToDo: pulled code from stack overflow, modify or credit.
        import glob

        before_arg = line.rfind(" ", 0, startidx)
        if before_arg == -1:
            return  # arg not found

        fixed = line[before_arg + 1:startidx]  # fixed portion of the arg
        arg = line[before_arg + 1:endidx]
        pattern = arg + '*'

        completions = []
        for path in glob.glob(pattern):
            path = self._append_slash_if_dir(path)
            completions.append(path.replace(fixed, "", 1))
        return completions

    def complete_show(self, *args):
        text = args[0]
        return ["%s " % x for x in self.get_headings() if x.lower().startswith(text.lower())]

    def complete_sort(self, *args):
        text = args[0]
        return ["%s " % x for x in self.get_headings() + ["reverse"] if x.lower().startswith(text.lower())]

    def complete_write(self, *args):
        line, startidx, endidx = args[1:]
        # ToDo: pulled code from stack overflow, modify or credit.
        import glob

        before_arg = line.rfind(" ", 0, startidx)
        if before_arg == -1:
            return  # arg not found

        fixed = line[before_arg + 1:startidx]  # fixed portion of the arg
        arg = line[before_arg + 1:endidx]
        pattern = arg + '*'

        completions = []
        for path in glob.glob(pattern):
            path = self._append_slash_if_dir(path)
            completions.append(path.replace(fixed, "", 1))
        return completions

    def help_bash(self):
        _stdout('''\
Run bash commands from the DbBuddy Live Session.
Be careful!! This is not sand-boxed in any way, so give the 'bash' command
all the respect you would normally give the terminal window.\n
''', format_in=GREEN, format_out=self.terminal_default)

    def help_database(self):
        _stdout('''\
Reset the database(s) to be searched. Separate multiple databases with spaces.
Currently set to: {0}{1}{2}
Valid choices: {0}{3}\n
'''.format(YELLOW, ", ".join(self.dbbuddy.databases), GREEN, ", ".join(["all"] + DATABASES)),
            format_in=GREEN, format_out=self.terminal_default)

    def help_delete(self):
        _stdout('''\
Remove records completely from the Live Session. Be careful, this is permanent.
Choices are:
    search-terms, st: Delete all search terms from live session
    failures, f:      Clear list of failures
    trash-bin, tb:    Empty the trash bin
    records, recs:    Delete all the main list of records (leaving the trash bin alone)
    all:              Delete everything\n
''', format_in=GREEN, format_out=self.terminal_default)

    def help_failures(self):
        _stdout('''\
Print the status of any failures the Live Session has encountered.\n
''', format_in=GREEN, format_out=self.terminal_default)

    def help_fetch(self):
        _stdout('''\
Retrieve full records for all accessions in the main record list.
If requesting more than 50 Mbp of sequence data, you will be prompted to confirm the command.\n
''', format_in=GREEN, format_out=self.terminal_default)

    def help_format(self):
        _stdout('''\
Set the output format:
    Valid choices            ->  ["ids", "accessions", "summary", "full-summary", <SeqIO formats>]
    ids or accessions        ->  Simple list of all accessions in the buffer
    summary or full-summary  ->  Information about each record
    <SeqIO format>           ->  Full sequence records, in any sequence file format
                                 supported by BioPython (e.g. gb, fasta, clustal)
                                 See http://biopython.org/wiki/SeqIO#File_Formats for details\n
''', format_in=GREEN, format_out=self.terminal_default)

    def help_keep(self):
        _stdout('''\
Further refine your results with search terms:
    - Records that DO NOT MATCH your filters are relegated to the 'trash bin'; return them to the main list
      with 'restore' command.
    - Multiple filters can be included at the same time, each enclosed in quotes and separated by spaces.
    - Records are searched exhaustively by default; to restrict the search to a given column/field, prefix
      the filter with '(<column name>)'. E.g., '(organism)Rattus'.
    - To filter by sequence length, the following operators are recognized: =, >, >=, <, and <=
      Use these operators inside the column prefix. E.g., '(length>300)'
    - Regular expressions are understood (https://docs.python.org/3/library/re.html).
    - Searches are case sensitive. To make them insensitive, prefix the filter with '(?i)'. E.g., '(?i)HuMaN'. \n
''', format_in=GREEN, format_out=self.terminal_default)

    def help_quit(self):
        _stdout("End the live session.\n\n", format_in=GREEN, format_out=self.terminal_default)

    def help_load(self):
        _stdout("Recover the contents of a previous session that crashed.\n\n", format_in=GREEN,
                format_out=self.terminal_default)

    def help_trash(self):
        _stdout('''\
Output the records held in the trash bin (out_format currently set to '{0}{1}{2}')
Optionally include an integer value and/or column name(s) to limit
the number of records and amount of information per record displayed.\n
'''.format(YELLOW, self.dbbuddy.out_format, GREEN), format_in=GREEN, format_out=self.terminal_default)

    def help_remove(self):
        _stdout('''\
Further refine your results with search terms:
    - Records that MATCH your filters are relegated to the 'trash bin'; return them to the main list
      with the 'restore' command
    - Multiple filters can be included at the same time, each enclosed in quotes and separated by spaces.
    - Records are searched exhaustively by default; to restrict the search to a given column/field, prefix
      the filter with '(<column name>)'. E.g., '(organism)Rattus'.
    - To filter by sequence length, the following operators are recognized: =, >, >=, <, and <=
      Use these operators inside the column prefix. E.g., '(length>300)'
    - Regular expressions are understood (https://docs.python.org/3/library/re.html).
    - Searches are case sensitive. To make them insensitive, prefix the filter with '(?i)'. E.g., '(?i)HuMaN'.\n
''', format_in=GREEN, format_out=self.terminal_default)

    def help_restore(self):
        _stdout('''\
Return a subset of filtered records back into the main list (use '%srestore *%s' to recover all records)
    - Multiple filters can be included at the same time, each enclosed in quotes and separated by spaces.
    - Records are searched exhaustively by default; to restrict the search to a given column/field, prefix
      the filter with '(<column name>)'. E.g., '(organism)Rattus'.
    - To filter by sequence length, the following operators are recognized: =, >, >=, <, and <=
      Use these operators inside the column prefix. E.g., '(length>300)'
    - Regular expressions are understood (https://docs.python.org/3/library/re.html).
    - Searches are case sensitive. To make them insensitive, prefix the filter with '(?i)'. E.g., '(?i)HuMaN'.\n
''' % (YELLOW, GREEN), format_in=GREEN, format_out=self.terminal_default)

    def help_save(self):
        _stdout('''\
Save your live session in DB format; the session can be restored later with the '{0}load{1}' command.
Supply the path where the file should be written to.\n
'''.format(YELLOW, GREEN), format_in=GREEN, format_out=self.terminal_default)

    def help_search(self):
        _stdout('''\
Search databases (currently set to {0}{1}{2}). If search terms
are supplied summary info will be downloaded, if accession numbers
are supplied then full sequence records will be downloaded.\n
'''.format(YELLOW, self.dbbuddy.databases, GREEN), format_in=GREEN, format_out=self.terminal_default)

    def help_show(self):
        _stdout('''\
Output the records held in the Live Session (output format currently set to '{0}{1}{2}')
Optionally include an integer value and/or column name(s) to limit
the number of records and amount of information per record displayed.\n
'''.format(YELLOW, self.dbbuddy.out_format, GREEN), format_in=GREEN, format_out=self.terminal_default)

    def help_sort(self):
        _stdout('''\
Alter the order that records are shown in.
By default records will be ordered by accession number, although any number of column
headings (case sensitive) can be included for a more customized sort.
To reverse the sort order include the keyword 'rev' or 'reverse'.\n
''', format_in=GREEN, format_out=self.terminal_default)

    def help_status(self):
        _stdout("Display the current state of your Live Session, including how many accessions and full records "
                "have been downloaded.\n\n", format_in=GREEN, format_out=self.terminal_default)

    def help_undo(self):
        _stdout('''\
Revert the most recent change to your live session.
The commands that can be undone are: - keep
                                     - remove
                                     - restore
                                     - search
                                     - fetch
                                     - format
                                     - database
                                     - load

Caution: Undo can only restore a sinlge previous step. All other history is lost.\n
''', format_in=GREEN, format_out=self.terminal_default)

    def help_write(self):
        _stdout('''\
Send records to a file (format currently set to '{0}{1}{2}').
Supply the file name to be written to.\n
'''.format(YELLOW, self.dbbuddy.out_format, GREEN), format_in=GREEN, format_out=self.terminal_default)

# DL everything
"""
def download_everything(_dbbuddy):
    # Get sequences from UniProt
    uniprot = UniProtRestClient(_dbbuddy)
    uniprot.fetch_proteins()
    return
    # Get sequences from Ensembl
    ensembl = EnsemblRestClient(_dbbuddy)
    ensembl.fetch_nucleotides()

    # Get sequences from genbank
    refseq = NCBIClient(_dbbuddy)
    refseq.gi2acc()
    refseq.fetch_nucliotides()
    refseq.fetch_proteins()

    return _dbbuddy
"""

"""
def retrieve_accessions(_dbbuddy):
    check_all = False if _dbbuddy.databases else True

    if "uniprot" in _dbbuddy.databases or check_all:
        uniprot = UniProtRestClient(_dbbuddy)
        uniprot.search_proteins()

    return _dbbuddy  # TEMPORARY

    if "ncbi_nuc" in _dbbuddy.databases or "ncbi_prot" in _dbbuddy.databases or check_all:
        refseq = NCBIClient(_dbbuddy)
        refseq.gi2acc()
        # refseq.search_nucliotides()

    return _dbbuddy
"""


def retrieve_summary(_dbbuddy):
    check_all = False if _dbbuddy.databases else True
    if "uniprot" in _dbbuddy.databases or check_all:
        uniprot = _dbbuddy.server("uniprot")
        uniprot.search_proteins()

    if "ncbi_nuc" in _dbbuddy.databases or check_all:
        refseq = _dbbuddy.server("ncbi")
        refseq.search_ncbi("nucleotide")
        refseq.fetch_summary()

    if "ncbi_prot" in _dbbuddy.databases or check_all:
        refseq = _dbbuddy.server("ncbi")
        refseq.search_ncbi("protein")
        refseq.fetch_summary()

    if "ensembl" in _dbbuddy.databases or check_all:
        ensembl = _dbbuddy.server("ensembl")
        ensembl.search_ensembl()
        ensembl.fetch_summary()
        # ensembl.fetch_summary()

    return _dbbuddy


def retrieve_sequences(_dbbuddy):
    check_all = False if _dbbuddy.databases else True
    if "uniprot" in _dbbuddy.databases or check_all:
        uniprot = _dbbuddy.server("uniprot")
        uniprot.fetch_proteins()

    if "ncbi_nuc" in _dbbuddy.databases or check_all:
        refseq = _dbbuddy.server("ncbi")
        refseq.fetch_sequence("nucleotide")

    if "ncbi_prot" in _dbbuddy.databases or check_all:
        refseq = _dbbuddy.server("ncbi")
        refseq.fetch_sequence("protein")

    if "ensembl" in _dbbuddy.databases or check_all:
        ensembl = _dbbuddy.server("ensembl")
        ensembl.fetch_nucleotide()

    return _dbbuddy


# ################################################# COMMAND LINE UI ################################################## #
def argparse_init():
    import argparse

    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="DbBuddy.py", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                     description='''
\033[1mDatabaseBuddy\033[m
  Go forth to the servers of sequence, and discover.

\033[1mUsage examples\033[m:
  DbBuddy.py -ls (launch empty live session)
  DbBuddy.py "<accn1,accn2,accn3,...>" -<cmd>
  DbBuddy.py "<search term1, search term2,...>" -<cmd>
  DbBuddy.py "<accn1,search term1>" -<cmd>
  DbBuddy.py "/path/to/file_of_accns" -<cmd>
''')

    br.db_modifiers["database"]["choices"] = DATABASES
    br.flags(parser, ("user_input", "Specify accession numbers or search terms, "
                                    "either in a file or as a comma separated list"),
             br.db_flags, br.db_modifiers, VERSION)

    in_args = parser.parse_args()

    dbbuddy = []
    out_format = "summary" if not in_args.out_format else in_args.out_format

    try:
        if isinstance(in_args.user_input[0], TextIOWrapper) and in_args.user_input[0].buffer.raw.isatty():
                dbbuddy = DbBuddy()
                in_args.live_shell = True
        elif len(in_args.user_input) > 1:
            for search_set in in_args.user_input:
                dbbuddy.append(DbBuddy(search_set, in_args.database, out_format))

            dbbuddy = DbBuddy(dbbuddy, in_args.database, in_args.out_format)
        else:
            dbbuddy = DbBuddy(in_args.user_input[0], in_args.database, out_format)

    except br.GuessError:
        sys.exit("Error: SeqBuddy could not understand your input. "
                 "Check the file path or try specifying an input type with -f")

    return in_args, dbbuddy


def command_line_ui(in_args, dbbuddy, skip_exit=False):
    # ############################################## COMMAND LINE LOGIC ############################################## #
    def _print_recs(_dbbuddy):
        _dbbuddy.print()

    # Live Shell
    def launch_live_shell():
        # Create a temp file for crash handling
        temp_file = MyFuncs.TempFile(byte_mode=True)
        temp_file.open()
        try:  # Catch all exceptions and try to send error report to server
            LiveSearch(dbbuddy, temp_file)
        except SystemExit:
            pass
        except (KeyboardInterrupt, br.GuessError) as _e:
            print(_e)
        except Exception as _e:
            import traceback
            save_file = "./DbSessionDump_%s" % temp_file.name
            temp_file.save(save_file)
            tb = "%s\n" % CONFIG["user_hash"]
            for _line in traceback.format_tb(sys.exc_info()[2]):
                _line = re.sub('"/.*/(.*)?"', r'"\1"', _line)
                tb += _line
            tb = "%s: %s\n\n%s" % (type(_e).__name__, _e, tb)
            _stderr("\033[mThe live session has crashed with the following traceback:%s\n\n%s\n\n\033[mYour work has "
                    "been saved to %s, and can be loaded by launching DatabaseBuddy and using the 'load' "
                    "command.\n" % (RED, tb, save_file))

            send_diagnostic = True if CONFIG["diagnostics"] == "True" else False
            if not send_diagnostic:
                prompt = input("%sWould you like to send a crash report with the above "
                               "traceback to the developers ([y]/n)?\033[m" % BOLD)
                if prompt.lower() in ["y", "yes", ""]:
                    send_diagnostic = True

            else:
                _stderr("An error report with the above traceback is being sent to the BuddySuite developers because "
                        "you have elected to participate in the Software Improvement Program. To opt-out of this "
                        "program, re-run the BuddySuite installer and un-check the box on the 'Diagnostics' screen.\n")

            if send_diagnostic:
                _stderr("Preparing error report for FTP upload...\nSending...\n")
                br.error_report(tb)
                _stderr("Success, thank you.\n")

        temp_file.close()
        sys.exit()

    def _exit(tool, skip=skip_exit):
        if skip:
            return
        usage = br.Usage()
        usage.increment("PhyloBuddy", VERSION.short(), tool)
        usage.save()
        sys.exit()

    if in_args.live_shell:
        launch_live_shell()

    """
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

        dbbuddy.print()
        sys.exit()

    # Retrieve Accessions
    if in_args.retrieve_accessions:
        if not in_args.out_format:
            dbbuddy.out_format = "ids"
        retrieve_accessions(dbbuddy)
        dbbuddy.print()
        sys.exit()

    if in_args.retrieve_sequences:
        sys.exit()
    """
    # Guess database
    if in_args.guess_database:
        output = ""
        if len(dbbuddy.records) > 0:
            records_by_db = {}
            for accn, record in dbbuddy.records.items():
                records_by_db.setdefault(record.database, [])
                records_by_db[record.database].append(accn)

            for db in records_by_db:
                records_by_db[db] = sorted(records_by_db[db])

            db_order = sorted(list(records_by_db))
            ordered_records = []
            accn_len = 11
            for db in db_order:
                for accn in records_by_db[db]:
                    accn_len = accn_len if len(accn) <= accn_len else len(accn)
                    ordered_records.append((accn, db))

            output += "%sDatabase\n" % "# Accession".ljust(accn_len + 3)
            for accn, db in ordered_records:
                output += "%s%s\n" % (str(accn).ljust(accn_len + 3), db)
            output += "\n"

        if len(dbbuddy.search_terms) > 0:
            output += "# Search terms\n"
            for term in dbbuddy.search_terms:
                output += "%s\n" % term

        if len(dbbuddy.records) == 0 and len(dbbuddy.search_terms) == 0:
            output += "Nothing to return\n"

        _stdout(output)
        sys.exit()

    # Default to LiveSearch
    launch_live_shell()


def main():
    initiation = []
    try:
        initiation = argparse_init()
        command_line_ui(*initiation)
    except (KeyboardInterrupt, br.GuessError) as e:
        print(e)
    except SystemExit:
        pass
    except Exception as e:
        function = ""
        for next_arg in vars(initiation[0]):
            if getattr(initiation[0], next_arg) and next_arg in br.db_flags:
                function = next_arg
                break
        br.send_traceback("DbBuddy", function, e)

if __name__ == '__main__':
    main()
