#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: DatabaseBuddy.py
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
try:
    import buddy_resources as br
except ImportError:
    try:
        import buddysuite.buddy_resources as br
    except AttributeError:
        from . import buddy_resources as br

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
import dill
import glob

# Third party
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
RETRIEVAL_TYPES = ["protein", "nucleotide"]
FORMATS = ["ids", "accessions", "summary", "full-summary", "clustal", "embl", "fasta", "fastq", "fastq-sanger",
           "fastq-solexa", "fastq-illumina", "genbank", "gb", "imgt", "nexus", "phd", "phylip", "seqxml",
           "stockholm", "tab", "qual"]
CONFIG = br.config_values()
VERSION = br.Version("DatabaseBuddy", 1, "2.7", br.contributors, {"year": 2017, "month": 6, "day": 16})

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
        self.trash_bin = OrderedDict()  # If records are filtered out, send them here instead of deleting them
        self.out_format = _out_format.lower()
        self.failures = OrderedDict()  # The key for these is a hash of the Failure, and the values are Failure objects
        self.databases = check_database(_databases)
        self.server_clients = {"ncbi": False, "ensembl": False, "uniprot": False}
        self.memory_footprint = 0

        # Empty DbBuddy object
        if not _input:
            return

        # DbBuddy objects
        elif type(_input) == list:
            for _dbbuddy in _input:
                if type(_dbbuddy) != DbBuddy:
                    raise TypeError("List of non-DbBuddy objects passed into DbBuddy as _input. %s" % _dbbuddy)

                self.search_terms += _dbbuddy.search_terms
                # Should probably check whether values are already set first, because update() will overwrite any common
                # records between the two dicts
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
        elif type(_input) == str and os.path.isfile(_input):
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
                    self.records[_record.accession] = _record

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
        if isinstance(other, type(self)) and self.out_format == other.out_format:
            recs1 = "".join(["%s%s" % (key, rec) for key, rec in self.records.items()])
            recs2 = "".join(["%s%s" % (key, rec) for key, rec in other.records.items()])
            print(recs1)
            if recs1 == recs2:
                return True
        return False

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
            # try:
            if mode == "keep" and not _rec.search(regex):
                self.trash_bin[_id] = _rec
            elif mode == "remove" and _rec.search(regex):
                self.trash_bin[_id] = _rec
            elif mode == "restore" and _rec.search(regex):
                self.records[_id] = _rec
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

    def server(self, _server):
        if _server not in ["uniprot", "ncbi", "ensembl"]:
            raise ValueError('"uniprot", "ncbi", and "ensembl" are the only valid options, not %s' % _server)
        if self.server_clients[_server]:
            return self.server_clients[_server]
        if _server == "uniprot":
            client = UniProtRestClient(self)
        elif _server == "ncbi":
            client = NCBIClient(self)
        else:  # _server must be "ensembl"
            client = EnsemblRestClient(self)
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
        ToDo: Allow slices of records to be returned (e.g., [5:-9])
        :param _num: Limit the number of rows (records) returned, otherwise everything is output
        :param quiet: suppress stderr
        :param columns: Variable, list of column names to include in summary output
        :param destination: a file path or handle to write to
        :param group: Either 'records' or 'trash_bin'
        :return: Nothing.
        """
        group = self.trash_bin if group == "trash_bin" else self.records
        _num = len(group) if abs(_num) >= len(group) or not _num else _num

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
            br._stderr(errors_etc, quiet)

        _output = ""
        # If negative number requests, return from back of list
        records = list(group.items())[:_num] if _num > 0 else list(group.items())[len(group) + _num:]

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
            for _accession, _rec in records:
                if self.out_format in ["ids", "accessions"]:
                    lines.append([_accession])

                elif self.out_format in ["summary", "full-summary"]:
                    headings = ["ACCN", "DB", "Type"]
                    if "length" in _rec.summary:
                        headings.append("length")
                    headings += [heading for heading, _value in _rec.summary.items()
                                 if heading not in ["comments", "length"]]
                    if "comments" in _rec.summary:
                        headings.append("comments")
                    headings.append("record")

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

                    if "Type" in headings:
                        if _rec.type:
                            current_group[-1].append(_rec.type[:4])
                        else:
                            current_group[-1].append("")
                        attrib_counter += 1

                    for heading in headings:
                        if heading in ["ACCN", "DB", "Type"]:
                            continue
                        if heading in _rec.summary:
                            _value = _rec.summary[heading]
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
            # Make sure IDs are not too long for GenBank format
            if self.out_format in ["gb", "genbank"]:
                for _accession, _rec in group.items():
                    if len(_accession) > 16:
                        br._stderr("Warning: Genbank format returned an 'ID too long' error. Format changed to EMBL.\n\n")
                        self.out_format = "embl"
                        break

            records = [_rec.record for _accession, _rec in records if _rec.record]
            tmp_file = br.TempFile()
            SeqIO.write(records, tmp_file.get_handle("w"), self.out_format)
            _output += "%s\n" % tmp_file.read()
            tmp_file.clear()

        if not destination:
            _stdout("{0}\n".format(_output.rstrip()))
        else:
            # remove any escape characters and convert space padding to tabs if writing the file
            _output = re.sub("\\033\[[0-9]*m", "", _output)
            _output = re.sub(" +\n", "\n", _output)
            destination.write(_output)


# ################################################# SUPPORT CLASSES ################################################## #
class Record(object):
    def __init__(self, _accession, _version=None, _record=None, summary=None, _size=None,
                 _database=None, _type=None, _search_term=None):
        self.accession = _accession
        self.version = _version
        self.record = _record  # SeqIO record
        self.summary = summary if summary else OrderedDict()  # Dictionary of attributes
        self.size = _size if _size in [None, ''] else int(_size)
        self.database = _database
        self.type = check_type(_type)
        self.search_term = _search_term  # In case the record was the result of a particular search

    def guess_database(self):
        import re  # This is to allow the objects to be loaded again later
        # RefSeq
        # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/
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
        # http://www.uniprot.org/help/accession_numbers
        elif re.match("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", self.accession):
            self.database = "uniprot"
            self.type = "protein"

        # Ensembl stable ids
        # http://www.ensembl.org/info/genome/stable_ids/index.html
        # First match is normal ensembl, the second is flybase
        elif re.match("^ENS[A-Z]*[0-9]+", self.accession) or re.match("^FB[a-z]{2}[0-9]+", self.accession):
            self.database = "ensembl"
            self.type = "nucleotide"

        # GenBank
        elif re.match("^[A-Z][0-9]{5}$|^[A-Z]{2}[0-9]{6}(\.([0-9]+))?$", self.accession):  # Nucleotide
            self.database = "ncbi_nuc"
            self.type = "nucleotide"

        elif re.match("^[A-Z]{3}[0-9]{5}(\.([0-9]+))?$", self.accession):  # Protein
            self.database = "ncbi_prot"
            self.type = "protein"

        elif re.match("[0-9][A-Z0-9]{3}(_[A-Z0-9])?(\.([0-9]+))?$", self.accession) \
                and re.search("[A-Z]", self.accession):  # PDB
            self.database = "ncbi_prot"
            self.type = "protein"

        elif re.match("^[A-Z]{4}[0-9]{8,10}(\.([0-9]+))?$", self.accession):  # Whole Genome
            self.database = "ncbi_nuc"
            self.type = "nucleotide"

        elif re.match("^[A-Z]{5}[0-9]{7}(\.([0-9]+))?$", self.accession):  # MGA (Mass sequence for Genome Annotation)
            self.database = "ncbi_prot"
            self.type = "protein"

        # Catch accn.version
        version = re.search("^(.*?)\.([0-9]+)$", self.accession)
        if version:
            self.version = version.group(2)

        return

    def search(self, regex):
        import re  # This is to allow the objects to be loaded again later
        regex = ".*" if regex == "*" else regex  # This prevents a crash
        # Default is case-senstive, so check if the user desires otherwise
        if regex[:2] in ["i?", "?i"]:
            flags = re.IGNORECASE
            regex = regex[2:]
        else:
            flags = 0

        column = re.match("\((.*?)\)", regex)
        if column:
            column = column.group(1)
            # Special case, if user is searching sequence length
            if re.match("length.+", column.strip(), flags=re.IGNORECASE):
                if not re.match("^length[ =<>]+[0-9]+$", column, flags=re.IGNORECASE):
                    raise ValueError("Invalid syntax for seaching 'length': %s" % column)

                limit = re.search("length *([ =<>]+)([0-9]+)", column, flags=re.IGNORECASE)
                operator = limit.group(1).strip()
                limit = int(limit.group(2))
                if "length" not in self.summary:
                    return False

                length = int(self.summary["length"])
                if operator not in ["=", ">", ">=", "<", "<="]:
                    raise ValueError("Invalid operator: %s" % operator)

                if operator == "=" and length == limit:
                    return True
                elif operator == ">" and length > limit:
                    return True
                elif operator == ">=" and length >= limit:
                    return True
                elif operator == "<" and length < limit:
                    return True
                elif operator == "<=" and length <= limit:
                    return True
                else:
                    return False

            # Strip off column syntax
            regex = re.search("^\(.*?\)(.*)", regex, flags=flags)
            regex = None if not regex else regex.group(1).strip()

            if column.lower() == "accn":
                if re.search(regex, str(self.accession), flags=flags):
                    return True

            if column.lower() == "type":
                if re.search(regex, str(self.type), flags=flags):
                    return True

            if column.lower() == "db":
                if re.search(regex, str(self.database), flags=flags):
                    return True

            if column in self.summary and not regex:  # This will return everything with the given column
                return True

            if column in self.summary and re.search(regex, str(self.summary[column]), flags=flags):
                return True
            else:
                return False

        for param in [self.accession, self.database, self.type, self.search_term]:
            if re.search(regex, str(param), flags=flags):
                return True

        for _key, _value in self.summary.items():
            if re.search(regex, _key, flags=flags) or re.search(regex, str(_value), flags=flags):
                return True

        if self.record:
            if re.search(regex, self.record.format("embl"), flags=flags):
                return True
        # If nothing hits, default to False
        return False

    def update(self, new_rec):
        self.accession = new_rec.accession if new_rec.accession else self.accession
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
        """
        Keep track of failed attempts to get records from the databases
        :param query: What was searched for
        :type query: str
        :param error_message: Either an automatically generated message or something custom that is informative
        :type error_message: str
        """
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


def _stdout(message, quiet=False, format_in=None, format_out=None):
    output = ""
    if format_in:
        format_in = format_in if type(format_in) == list else [format_in]
        for _format in format_in:
            if not re.search("\\033\[[0-9]*m", _format):
                raise AttributeError('Malformed format_in attribute escape code')
        output += "".join(format_in)

    if format_out:
        format_out = format_out if type(format_out) == list else [format_out]
        for _format in format_out:
            if not re.search("\\033\[[0-9]*m", _format):
                raise AttributeError('Malformed format_out attribute escape code')
        output += "%s%s" % (message, "".join(format_out))
    else:
        output += "%s\033[m" % message

    if not quiet:
        sys.stdout.write(output)
        sys.stdout.flush()
    return


def terminal_colors():
    colors = [CYAN, GREEN, RED, YELLOW, GREY, MAGENTA]
    _counter = 0
    while True:
        try:
            yield colors[_counter]
        except IndexError:
            _counter = 0
            yield colors[_counter]
        _counter += 1


def check_database(database=None):
    if not database:
        return DATABASES
    if type(database) != list:
        database = [database]
    database = [x.lower() for x in database]
    if 'all' in database:
        return DATABASES
    _output = []
    for _db in database:
        if _db in DATABASES:
            _output.append(_db)
        else:
            br._stderr("Warning: '%s' is not a valid database choice, omitted.\n" % _db)
    if not _output:
        br._stderr("Warning: No valid database choice provided. Setting to default 'all'.\n")
        _output = DATABASES
    return _output


def check_type(_type):
    _type = None if not _type else _type.lower()
    if _type in ["p", "pr", "prt", "prtn", "prn", "prot", "protn", "protien", "protein"]:
        _type = "protein"
    elif _type in ["n", "ncl", "nuc", "dna", "nt", "gene", "transcript", "nucleotide"]:
        _type = "nucleotide"

    if _type and _type not in RETRIEVAL_TYPES:
        br._stderr("Warning: '%s' is not a valid choice for '_type'. Setting to default 'protein'.\n" % _type)
        _type = "protein"
    return _type


def retrieve_summary(_dbbuddy):
    check_all = False if _dbbuddy.databases else True
    if "uniprot" in _dbbuddy.databases or check_all:
        uniprot = _dbbuddy.server("uniprot")
        uniprot.search_proteins()

    if "ncbi_nuc" in _dbbuddy.databases or check_all:
        refseq = _dbbuddy.server("ncbi")
        refseq.search_ncbi("nucleotide")
        refseq.fetch_summaries("ncbi_nuc")

    if "ncbi_prot" in _dbbuddy.databases or check_all:
        refseq = _dbbuddy.server("ncbi")
        refseq.search_ncbi("protein")
        refseq.fetch_summaries("ncbi_prot")

    if "ensembl" in _dbbuddy.databases or check_all:
        ensembl = _dbbuddy.server("ensembl")
        ensembl.search_ensembl()
        ensembl.fetch_summaries()

    return _dbbuddy


def retrieve_sequences(_dbbuddy):
    check_all = False if _dbbuddy.databases else True
    if "uniprot" in _dbbuddy.databases or check_all:
        uniprot = _dbbuddy.server("uniprot")
        uniprot.fetch_proteins()

    if "ncbi_nuc" in _dbbuddy.databases or check_all:
        refseq = _dbbuddy.server("ncbi")
        refseq.fetch_sequences("nucleotide")

    if "ncbi_prot" in _dbbuddy.databases or check_all:
        refseq = _dbbuddy.server("ncbi")
        refseq.fetch_sequences("protein")

    if "ensembl" in _dbbuddy.databases or check_all:
        ensembl = _dbbuddy.server("ensembl")
        ensembl.fetch_nucleotide()

    return _dbbuddy


# ################################################# Database Clients ################################################# #
class GenericClient(object):
    def __init__(self, _dbbuddy, max_url=1000):
        self.dbbuddy = _dbbuddy
        self.http_errors_file = br.TempFile()
        self.results_file = br.TempFile()
        self.max_url = max_url
        self.lock = Lock()

    def parse_error_file(self):
        http_errors_file = self.http_errors_file.read().strip("//\n")
        if http_errors_file != "":
            _output = ""
            http_errors_file = http_errors_file.split("//")
            for error in http_errors_file:
                error = error.strip().split("\n")
                error = (error[0], "\n".join(error[1:])) if len(error) > 2 else (error[0], error[1])
                error = Failure(*error)
                if error.hash not in self.dbbuddy.failures:
                    self.dbbuddy.failures[error.hash] = error
                    _output += "%s\n" % error
            self.http_errors_file.clear()
            return _output  # Errors found
        else:
            return False  # No errors to report

    def write_error(self, msg, err):
        with self.lock:
            self.http_errors_file.write("%s\n%s\n//\n" % (msg, err))
        return

    def group_terms_for_url(self, terms):
        groups = [""]
        for term in terms:
            term = str(term)
            if len(term) + 1 > self.max_url:
                raise ValueError("The provided accession or search term is too long (>%s).\n"
                                 "The problematic string is:\n%s" % (self.max_url, term))
            if len(groups[-1]) + len(term) + 1 > self.max_url:
                groups[-1] = groups[-1].strip(",")
                groups.append("%s," % term)
            else:
                groups[-1] += "%s," % term
        groups[-1] = groups[-1].strip(",")
        return groups


class UniProtRestClient(GenericClient):
    # http://www.uniprot.org/help/uniprotkb_column_names
    def __init__(self, _dbbuddy, server='http://www.uniprot.org/uniprot'):
        GenericClient.__init__(self, _dbbuddy)
        self.server = server

    def query_uniprot(self, search_term, request_params):  # Multicore ready
        if type(request_params) == list:  # In case it's coming in from multicore run
            request_params = request_params[0]
        search_term = re.sub(" ", "+", search_term)
        request_string = ""
        for _param, _value in request_params.items():
            _value = re.sub(" ", "+", _value)
            request_string += "&{0}={1}".format(_param, _value)

        try:
            request = Request("{0}?query={1}{2}".format(self.server, search_term, request_string))
            response = urlopen(request)
            response = response.read().decode("utf-8")
            response = re.sub("^Entry.*\n", "", response, count=1)
            with self.lock:
                self.results_file.write("# Search: %s\n%s//\n" % (search_term, response))

        except HTTPError as err:
            self.write_error("Uniprot search failed for '%s'" % search_term, err)

        except URLError as err:
            if "Errno 8" in str(err):
                self.write_error("Uniprot request failed, are you connected to the internet?", err)
            else:
                self.write_error("Uniprot request failed", err)
        return

    def count_hits(self):
        # Limit URLs to 2,083 characters
        _count = 0
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

        for search_term in search_terms:
            self.query_uniprot(search_term, {"format": "list"})
            content = re.sub("(#.*?\n|[\n /]+$)", "", self.results_file.read())
            content = content.split("\n")
            _count += len(content) if content[0] != '' else 0
            self.results_file.clear()

        self.parse_error_file()
        return _count

    def search_proteins(self):
        # start by determining how many results we would get from all searches.
        self.results_file.clear()
        _count = self.count_hits()

        if _count == 0:
            br._stderr("Uniprot returned no results\n\n")
            return

        else:
            br._stderr("Retrieving summary data for %s records from UniProt\n" % _count)

        # download the tab info on all or subset
        params = {"format": "tab", "columns": "id,entry name,length,organism-id,organism,protein names,comments"}
        runtime = br.RunTime(prefix="\t")
        runtime.start()
        if len(self.dbbuddy.search_terms) > 1:
            br._stderr("Querying UniProt with %s search terms (Ctrl+c to abort)\n" % len(self.dbbuddy.search_terms))
            br.run_multicore_function(self.dbbuddy.search_terms, self.query_uniprot, max_processes=10, quiet=True,
                                      func_args=[params])
        else:
            br._stderr("Querying UniProt with the search term '%s'...\n" % self.dbbuddy.search_terms[0])
            self.query_uniprot(self.dbbuddy.search_terms[0], params)
        runtime.end()
        self.parse_error_file()

        content = re.sub("(#.*?\n|[\n /]+$)", "", self.results_file.read().strip())
        results = content.split("//")
        result_count = 0
        for result in [str(x) for x in results]:
            result = result.strip().split("\n")
            for hit in result:
                result_count += 1
                hit = hit.split("\t")
                if len(hit) == 6:  # In case 'comments' isn't returned
                    raw = OrderedDict([("entry_name", hit[1]), ("length", int(hit[2])), ("TaxId", hit[3]),
                                       ("organism", hit[4]), ("protein_names", hit[5]), ("comments", "")])
                else:
                    raw = OrderedDict([("entry_name", hit[1]), ("length", int(hit[2])), ("TaxId", hit[3]),
                                       ("organism", hit[4]), ("protein_names", hit[5]), ("comments", hit[6])])

                self.dbbuddy.records[hit[0]] = Record(hit[0], _database="uniprot", _type="protein",
                                                      _search_term=result[0], summary=raw, _size=int(hit[2]))
        br._stderr("\t%s records received.\n" % result_count)

    def fetch_proteins(self):
        self.results_file.clear()
        _records = [_rec for _accession, _rec in self.dbbuddy.records.items() if
                    _rec.database == "uniprot" and not _rec.record]

        if not _records:
            return

        br._stderr("Requesting %s full records from UniProt...\n" % len(_records))
        accessions = self.group_terms_for_url([_rec.accession for _rec in _records])
        runtime = br.RunTime(prefix="\t")
        runtime.start()
        params = {"format": "txt"}
        if len(accessions) > 1:
            br.run_multicore_function(accessions, self.query_uniprot, max_processes=10, quiet=True, func_args=[params])
        else:
            self.query_uniprot(accessions[0], params)

        runtime.end()
        errors = self.parse_error_file()
        if errors:
            br._stderr("{0}{1}The following errors were encountered while querying UniProt with "
                       "fetch_proteins():{2}\n{3}{4}".format(RED, UNDERLINE, NO_UNDERLINE, errors, DEF_FONT))

        data = self.results_file.read().strip()
        data = re.sub("# Search.*?\n", "", data)
        data = re.sub("//(\n//)+", "//\n", data)
        data = re.sub("^//\n*", "", data)
        data = re.sub("//\n\n+", "//\n", data)
        if data in ["", "//\n"]:
            br._stderr("No sequences returned\n\n")
            return

        self.results_file.write(data, "w")

        _records = SeqIO.parse(self.results_file.get_handle("r"), "swiss")
        for _rec in _records:
            self.dbbuddy.records[_rec.id].record = _rec
        return


class NCBIClient(GenericClient):
    def __init__(self, _dbbuddy):
        GenericClient.__init__(self, _dbbuddy)
        self.Entrez = Entrez
        self.Entrez.email = CONFIG["email"]
        self.Entrez.tool = "buddysuite"
        self.max_attempts = 5  # NCBI throws a lot of 503 errors, so keep trying until we get through...
        self.tries = 0

    def _mc_query(self, query, func_args):
        """
        Make a request to Entrez for some data
        :param query: Appropriately sized/formatted request string
        :param func_args: tool = "esummary_taxa", "esummary_seq", or "efetch_seq"
        :return:
        """
        tool, db = func_args
        if db in ["ncbi_nuc", "ncbi_prot"]:
            db = "nucleotide" if db == "ncbi_nuc" else "protein"
        if db and db not in ["nucleotide", "protein"]:
            raise ValueError("Unknown type '%s', choose between 'nucleotide' and 'protein" % db)
        handle = None
        timer = br.time()
        counter = int(self.max_attempts)
        while counter > 0:
            counter -= 1
            try:
                if tool == "esummary_taxa":
                    # Example query of taxa ids: "649,734,1009,2302"
                    handle = Entrez.esummary(db="taxonomy", id=query, retmax=10000)
                elif tool == "esummary_seq":
                    # Example query of ACCNs: "XP_010103297.1,XP_010103298.1,XP_010103299.1"
                    handle = Entrez.esummary(db=db, id=query, retmax=10000)
                elif tool == "efetch_seq":
                    # Example query of ACCNs: "XP_010103297.1,XP_010103298.1,XP_010103299.1"
                    handle = Entrez.efetch(db=db, id=query, rettype="gb", retmode="text", retmax=10000)
                elif tool == "esearch":
                    count = Entrez.read(Entrez.esearch(db=db, term=re.sub('[\'"]', '', query), rettype="count"))["Count"]
                    handle = Entrez.esearch(db=db, term=re.sub('[\'"]', '', query), retmax=count, idtype='acc')
                else:
                    raise ValueError("_mc_query() 'tool' argument must be in 'esummary_taxa', "
                                     "'esummary_seq', or 'efetch_seq'")

                # This is a throttle so the NCBI server isn't spammed too rapidly
                timer = br.time() - timer
                if timer < 1:
                    sleep(1 - timer)
                break
            except HTTPError as err:
                if err.getcode() != 503 or counter == 0:
                    self.write_error("NCBI request failed: %s" % query, err)
                    break
                sleep(1)
            except ConnectionResetError as err:
                if "[Errno 54] Connection reset by peer" not in str(err) or counter == 0:
                    self.write_error("NCBI request failed: %s" % query, err)
                    break
                sleep(1)
            except URLError as err:
                if "Errno 8" in str(err):
                    self.write_error("NCBI request failed, are you connected to the internet?", err)
                else:
                    self.write_error("NCBI request failed", err)
                break

        if handle:
            if tool == "efetch_seq":
                result = "%s\n" % handle.read().strip()
            else:
                result = "%s\n### END ###\n" % handle.read().strip()

            with self.lock:
                self.results_file.write(result)
        return

    def search_ncbi(self, _type):
        """
        Query NCBI with search terms (not accns or GI nums)
        :param _type: "nucleotide" or "protein"
        :return:
        """
        if not self.dbbuddy.search_terms:
            return
        self.results_file.clear()
        if len(self.dbbuddy.search_terms) > 1:
            br.run_multicore_function(self.dbbuddy.search_terms, self._mc_query, func_args=["esearch", _type],
                                      max_processes=3, quiet=True)
        else:
            self._mc_query(self.dbbuddy.search_terms[0], func_args=["esearch", _type])

        self.parse_error_file()

        results = self.results_file.read().split("\n### END ###\n")
        results = [x for x in results if x != ""]
        accns = []
        for result in results:
            result = Entrez.read(StringIO(result))
            accns += result["IdList"]
        if not accns:
            br._stderr("NCBI returned no %s results\n\n" % _type)
            return
        for accn, rec in self.dbbuddy.records.items():
            if rec.accession in accns:
                del accns[accns.index(str(rec.accession))]
        database = 'ncbi_nuc' if _type == 'nucleotide' else 'ncbi_prot'
        for accn in accns:
            self.dbbuddy.records[accn] = Record(accn, _database=database, _type=_type)
        return

    def fetch_summaries(self, database):
        """
        Fetch metadata for all records
        :param database: in "ncbi_prot" and "ncbi_nuc"
        :return:
        """
        _type = "protein" if database == "ncbi_prot" else "nucleotide"
        self.results_file.clear()
        accns = [accn for accn, rec in self.dbbuddy.records.items() if rec.database == database]
        if not accns:
            return
        accn_searches = self.group_terms_for_url(accns)

        # Download all of the summaries
        self.results_file.clear()
        br._stderr("Retrieving %s %s record summaries from NCBI...\n" % (len(accns), _type))
        runtime = br.RunTime(prefix="\t")
        runtime.start()
        if len(accn_searches) > 1:
            br.run_multicore_function(accn_searches, self._mc_query, func_args=["esummary_seq", database],
                                      max_processes=3, quiet=True)
        else:
            self._mc_query(accn_searches[0], func_args=["esummary_seq", database])
        runtime.end()
        results = self.results_file.read().split("\n### END ###\n")
        results = [x for x in results if x != ""]

        # Sift through all the results and grab summary information
        records = {}
        taxa = []
        for result in results:
            try:  # This will catch and retry when the server fails on us
                summaries = [x for x in Entrez.parse(StringIO(result))]
                self.tries = 0
            except RuntimeError:
                if self.tries >= self.max_attempts:
                    self.tries = 0
                    break
                else:
                    br._stderr("Problem talking to NCBI, retrying...\n")
                    self.tries += 1
                    self.fetch_summaries(database)
                break

            for summary in [dict(x) for x in summaries]:
                # status can be 'live', 'dead', 'withdrawn', 'replaced'
                status = summary["Status"] if summary["ReplacedBy"] == '' else \
                    "%s->%s" % (summary["Status"], summary["ReplacedBy"])

                keys = ["TaxId",
                        "organism",
                        "length",
                        "comments",
                        "status"]
                values = [summary["TaxId"],
                          "",
                          summary["Length"],
                          summary["Title"],
                          status]
                rec_summary = OrderedDict([(key, value) for key, value in zip(keys, values)])

                if summary["TaxId"] not in taxa:
                    taxa.append(summary["TaxId"])

                accn = summary["AccessionVersion"]
                version = summary["AccessionVersion"].split(".")[-1]
                records[accn] = Record(accn, _version=version, summary=rec_summary, _type=_type,
                                       _size=rec_summary["length"], _database=database)

        # Get taxa names for all of the records retrieved
        self.results_file.clear()
        _taxa_ids = self.group_terms_for_url(taxa)
        if len(_taxa_ids) > 1:
            br.run_multicore_function(_taxa_ids, self._mc_query, func_args=["esummary_taxa", database],
                                      max_processes=3, quiet=True)
        else:
            self._mc_query(_taxa_ids[0], func_args=["esummary_taxa", database])
        self.parse_error_file()

        results = self.results_file.read().split("\n### END ###\n")
        results = [x for x in results if x]

        taxa = {}
        for result in results:
            for summary in Entrez.parse(StringIO(result)):
                taxa[summary["TaxId"]] = "Unclassified" if "ScientificName" not in summary \
                    else summary["ScientificName"]

        # Apply the taxa names that were downloaded
        for accn, rec in records.items():
            if rec.summary["TaxId"] in taxa:
                rec.summary["organism"] = taxa[rec.summary["TaxId"]]
            else:
                rec.summary["organism"] = "Unclassified"
        br._stderr("\t%s records received.\n" % len(records))

        # Update the dbbuddy object with all the new info
        for accn, rec in records.items():
            if str(accn) in self.dbbuddy.records:  # GI only records
                del self.dbbuddy.records[str(accn)]
            if rec.accession.split(".")[0] in self.dbbuddy.records:  # Un-versioned accns
                del self.dbbuddy.records[rec.accession.split(".")[0]]
            if rec.accession in self.dbbuddy.records:
                self.dbbuddy.records[rec.accession].update(rec)
            else:
                self.dbbuddy.records[rec.accession] = rec
        return

    def fetch_sequences(self, database):  # database in ["nucleotide", "protein"]
        db = "ncbi_nuc" if database == "nucleotide" else "ncbi_prot"
        accns = [accn for accn, _rec in self.dbbuddy.records.items() if _rec.database == db]
        if not accns:
            return
        self.results_file.clear()
        accns = self.group_terms_for_url(accns)
        runtime = br.RunTime(prefix="\t")
        br._stderr("Fetching full %s sequence records from NCBI...\n" % database)
        runtime.start()
        if len(accns) > 1:
            br.run_multicore_function(accns, self._mc_query, func_args=["efetch_seq", db],
                                      max_processes=3, quiet=True)
        else:
            self._mc_query(accns[0], func_args=["efetch_seq", db])
        self.parse_error_file()

        runtime.end()
        records = {}
        for rec in SeqIO.parse(self.results_file.get_handle("r"), "gb"):
            if rec.id not in records:
                records[rec.id] = rec
        br._stderr("\tDone\n")
        for accn, rec in records.items():
            self.dbbuddy.records[accn].record = rec
            version = re.search("^.*?\.([0-9]+)$", accn)
            if version:
                self.dbbuddy.records[accn].version = version.group(1)
        return


class EnsemblRestClient(GenericClient):
    def __init__(self, _dbbuddy, server='http://rest.ensembl.org/'):
        GenericClient.__init__(self, _dbbuddy)
        self.server = server
        self.species = self.perform_rest_action("info/species", headers={"Content-type": "application/json",
                                                                         "Accept": "application/json"})
        self.parse_error_file()
        if self.species:
            self.species = self.species["species"]
            self.species = {x["display_name"]: x for x in self.species if x["display_name"]}
        else:
            self.species = {}
        self.max_attempts = 5

    def _mc_search(self, species, args):
        identifier = args[0]
        self.dbbuddy.failures = {}
        data = self.perform_rest_action("lookup/symbol/%s/%s" % (species, identifier),
                                        headers={"Content-type": "application/json", "Accept": "application/json"})
        with self.lock:
            self.results_file.write("%s\n### END ###\n" % data)

    def perform_rest_action(self, endpoint, **kwargs):
        """
        :param endpoint: Ensembl specific REST commands
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
            self.max_attempts = 5  # Reset max_attempts, which may have been reduced if there were some 429 errors
            if request.get_header("Content-type") == "application/json":
                content = response.read().decode()
                data = json.loads(content)
            elif request.get_header("Content-type") == "text/x-seqxml+xml":
                data = SeqIO.parse(response, "seqxml")
            else:
                raise ValueError("Unknown request headers '%s'" % request.headers)
            return data

        except HTTPError as err:
            # check if we are being rate limited by the server
            err_code = err.getcode()
            if err_code == 429:
                if self.max_attempts == 0:
                    self.max_attempts = 5
                    self.write_error("Server Busy", err)
                    pass
                elif 'Retry-After' in err.headers:
                    self.max_attempts -= 1
                    retry = err.headers['Retry-After']
                    sleep(float(retry) + 1)
                    self.perform_rest_action(endpoint, **kwargs_backup)
            elif err_code == 400:
                pass
            else:
                self.write_error("Ensembl request failed: %s" % endpoint, err)

        except URLError as err:
            if "Errno 8" in str(err):
                self.write_error("Ensembl request failed, are you connected to the internet?", err)
            else:
                self.write_error("Ensembl request failed", err)
        return

    def search_ensembl(self):
        self.results_file.clear()
        species = [name for name, info in self.species.items()]
        for search_term in self.dbbuddy.search_terms:
            br._stderr("Searching Ensembl for %s...\n" % search_term)
            # br.run_multicore_function(species, self._mc_search, [search_term], quiet=True)
            # TODO: fix multicore --> Many REST requests are failing unless a single request is sent at a time
            for i in species:
                self._mc_search(i, [search_term])
            self.parse_error_file()
            results = self.results_file.read().split("\n### END ###")
            counter = 0
            for rec in results:
                rec = rec.strip()
                if not rec or rec in ["None", ""]:
                    continue
                counter += 1
                rec = re.sub("'", '"', rec)
                summary = json.loads(rec)
                accn = summary['id']
                size = abs(summary["start"] - summary["end"])
                _version = None if 'version' not in summary else summary['version']

                required_keys = ['display_name', 'species', 'biotype', 'object_type',
                                 'strand', 'assembly_name', 'description', 'version']

                for key in required_keys:
                    if key not in summary:
                        summary[key] = ''

                summary = OrderedDict([('name', summary['display_name']), ('length', size),
                                       ('organism', summary['species']),
                                       ('TaxId', self.species[summary['species']]['taxon_id']),
                                       ('biotype', summary['biotype']), ('object_type', summary['object_type']),
                                       ('strand', summary['strand']), ('assembly_name', summary['assembly_name']),
                                       ('comments', summary['description'])])

                rec = Record(accn, summary=summary, _version=_version,
                             _size=size, _database="ensembl", _type="nucleotide")

                if rec.accession in self.dbbuddy.records:
                    self.dbbuddy.records[rec.accession].update(rec)
                else:
                    self.dbbuddy.records[rec.accession] = rec

            if counter > 0:
                br._stderr("\t%s records received\n" % counter)
            else:
                br._stderr("Ensembl returned no results\n")

    def fetch_summaries(self):
        accns = [accn for accn, rec in self.dbbuddy.records.items() if rec.database == "ensembl"]
        data = {}
        for group in [accns[i:i+50] for i in range(0, len(accns), 50)]:  # Max 50 accessions per request
            query = self.perform_rest_action("lookup/id",
                                             data={"ids": group},
                                             headers={"Content-type": "application/json",
                                                      "Accept": "application/json"})
            data.update(query)
        self.parse_error_file()
        if not data:
            return

        for accn, results in data.items():
            if not results:
                continue
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

            version = None if "version" not in results else results["version"]
            rec = Record(accn, summary=summary, _version=version,
                         _size=size, _database="ensembl", _type="nucleotide")
            self.dbbuddy.records[accn].update(rec)
        return

    def fetch_nucleotide(self):
        accns = [accn for accn, rec in self.dbbuddy.records.items() if rec.database == "ensembl"]
        if len(accns) > 0:
            br._stderr("Fetching sequence from Ensembl...\n")
            runtime = br.RunTime(prefix="\t")
            runtime.start()
            for group in [accns[i:i+50] for i in range(0, len(accns), 50)]:  # Max 50 accessions per request
                data = self.perform_rest_action("sequence/id",
                                                data={"ids": group},
                                                headers={"Content-type": "text/x-seqxml+xml"})

                def_summary = OrderedDict([(x, None) for x in ['comments', 'organism', 'name']])
                for rec in data:
                    if rec.id not in self.dbbuddy.records:
                        self.dbbuddy.records[rec.id] = Record(rec.id, summary=def_summary)

                    summary = self.dbbuddy.records[rec.id].summary
                    rec.description = summary['comments']
                    rec.accession = rec.id
                    rec.name = rec.id
                    if summary['organism']:
                        species = re.search("([a-z])[a-z]*_([a-z]{1,3})", summary['organism'])
                        species = "%s%s" % (species.group(1).upper(), species.group(2))
                    else:
                        species = "Unknown"
                    new_id = "%s-%s" % (species, summary['name'])
                    self.dbbuddy.records[rec.id].record = rec
                    self.dbbuddy.records[rec.id].record.id = new_id
            self.parse_error_file()
            runtime.end()


# ################################################ MAIN API FUNCTIONS ################################################ #
class LiveShell(cmd.Cmd):
    def __init__(self, _dbbuddy, crash_file):
        """
        :param _dbbuddy: pre-instantiated DbBuddy object
        :param crash_file: br.TempFile object instantiated in binary mode
        """
        self.tmpdir = br.TempDir()
        self.terminal_default = "\033[m\033[40m%s" % WHITE
        cmd.Cmd.__init__(self)
        colors = terminal_colors()
        hash_heading = "%s#" % "#".join([next(colors) for _ in range(23)])
        _stdout('''{1}

{0} {1}{3}{2}Welcome to the DatabaseBuddy live shell{1} {0}{1}

{2}Type 'help' for a list of available commands or 'help <command>' for further details.
                  To end the session, use the 'quit' command.{1}

'''.format(hash_heading, self.terminal_default, BOLD, UNDERLINE))
        self.prompt = '{0}{1}DbBuddy{2}{1}>{2} '.format(MAGENTA, BOLD, self.terminal_default)

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

        if CONFIG["data_dir"]:
            self.history_path = "%s%scmd_history" % (CONFIG["data_dir"], os.sep)
        else:
            self.history_path = "%s%scmd_history" % (self.tmpdir.path, os.sep)
        try:
            if not os.path.isfile(self.history_path):
                open(self.history_path, "w", encoding="utf-8").close()
            else:
                open(self.history_path, "r").close()
        except PermissionError:
            self.history_path = "%s%scmd_history" % (self.tmpdir.path, os.sep)
            open(self.history_path, "w", encoding="utf-8").close()

        readline.read_history_file(self.history_path)
        readline.set_history_length(1000)

        # As implemented, one UnDo is possible (reload the most recent dump). Set self.undo to true every time a dump
        # occurs, and back to False if undo is used.
        self.undo = False

        br._stderr(self.terminal_default)  # This needs to be called here if stderr is going to format correctly
        if self.dbbuddy.records or self.dbbuddy.search_terms:
            retrieve_summary(_dbbuddy)
        else:
            _stdout("Your session is currently unpopulated. Use 'search' to retrieve records.\n",
                    format_out=self.terminal_default)
        self.hash = None
        self.shell_execs = []  # Only populate this if "bash" is called by the user
        self.usage = br.Usage()
        breakout = False
        while not breakout:
            try:
                self.cmdloop()
                breakout = True
            except KeyboardInterrupt:
                _stdout("\n")

    def precmd(self, line):
        # ToDo: Long commands are added to history, they are not output correctly in the terminal. For some reason they
        # are not cleared completely when you move to the next history index, leaving a truncated path at the prompt.
        # Need to track this bug down and squash it, just not sure how (i.e., don't want the 'and len(line) < 40' part)
        if line not in ["y", "n", "yes", "no"] and len(line) < 50:
            readline.write_history_file(self.history_path)
        else:
            readline.read_history_file(self.history_path)
        return line

    def postcmd(self, stop, line):
        command = line.split(" ")[0]
        self.usage.increment("LiveShell", VERSION.short(), command)
        return stop

    def dump_session(self):
        # Need to remove Lock()s to pickle
        for client in [client for db, client in self.dbbuddy.server_clients.items() if client]:
            client.lock = False
        self.crash_file.save("%s_undo" % self.crash_file.path)
        self.crash_file.open()
        dill.dump(self.dbbuddy, self.crash_file.handle, protocol=-1)
        self.crash_file.close()
        self.undo = True
        for client in [client for db, client in self.dbbuddy.server_clients.items() if client]:
            client.lock = Lock()

    def default(self, line):
        if line == "exit":
            self.do_quit()
        else:
            _stdout('*** Unknown syntax: %s\n\n' % line, format_in=RED, format_out=self.terminal_default)

    @staticmethod
    def _append_slash_if_dir(p):  # Used for expanding file path
            if p and os.path.isdir(p) and p[-1] != os.sep:
                return p + os.sep
            else:
                return p

    def get_headings(self):
        headings = ["ACCN", "DB", "Type", "record"]
        if len(self.dbbuddy.records) > 0:
            for _accn, _rec in self.dbbuddy.records.items():
                for heading, _value in _rec.summary.items():
                    if heading not in headings:
                        headings.append(heading)
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
        line = br.clean_regex(line)
        for _filter in line:
            for _key, _value in self.dbbuddy.filter_records(_filter, mode=mode).items():
                _errors[_key] += _value
            _stdout(tabbed.format(_filter, abs(current_count - len(self.dbbuddy.records))),
                    format_out=self.terminal_default)
            current_count = len(self.dbbuddy.records)

        output_message = "\n%s records remain.\n\n" % len(self.dbbuddy.records) if mode != "restore" \
            else "\n%s records remain in the trash bin.\n\n" % len(self.dbbuddy.trash_bin)

        _stdout(output_message, format_in=GREEN, format_out=self.terminal_default)
        self.dump_session()

    def do_EOF(self, line):
        return True

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
        line = re.sub("['\"]", "", line.strip())
        line = re.sub("\t+", " ", line)
        line = re.sub(", *", " ", line)
        line = line.lower()
        line = line.split(" ")
        new_database_list = []
        not_a_database = []
        for l in line:
            if l not in DATABASES and l != "all":
                not_a_database.append(l)
            else:
                new_database_list.append(l)
        if not_a_database:
            _stdout("Invalid database choice(s): %s.\n"
                    "Please select from %s\n" % (", ".join(not_a_database), ["all"] + DATABASES), format_in=RED,
                    format_out=self.terminal_default)
        if new_database_list:
            if "all" in new_database_list:
                self.dbbuddy.databases = DATABASES
            else:
                self.dbbuddy.databases = list(set(new_database_list))

            _stdout("Database search list updated to %s\n\n" % self.dbbuddy.databases, format_in=GREEN,
                    format_out=self.terminal_default)
        else:
            _stdout("Database search list not changed.\n\n", format_in=RED, format_out=self.terminal_default)
        self.dump_session()

    def do_delete(self, line="all"):
        if not self.dbbuddy.trash_bin and not self.dbbuddy.records \
                and not self.dbbuddy.search_terms and not self.dbbuddy.failures:
            _stdout("The live session is already empty.\n\n", format_in=RED, format_out=self.terminal_default)
            return

        line = "" if not line else line.lower()
        if line not in TRASH_SYNOS + RECORD_SYNOS + SEARCH_SYNOS \
                and not "all".startswith(line) and not "failures".startswith(line):
            _stdout("Sorry, I don't understand what you want to delete.\n Select from: all, main, trash-bin, "
                    "failures, search\n\n", format_in=RED, format_out=self.terminal_default)
            return

        if "failures".startswith(line) and line != "":
            if not self.dbbuddy.failures:
                _stdout("Failures list is already empty.\n\n", format_in=RED, format_out=self.terminal_default)
            else:
                confirm = br.ask("%sAre you sure you want to clear all %s failures (y/[n])?%s " %
                                 (RED, len(self.dbbuddy.failures), self.terminal_default), default="no")

                if not confirm:
                    _stdout("Aborted...\n", format_in=RED, format_out=self.terminal_default)
                else:
                    self.dbbuddy.failures = OrderedDict()
                    _stdout("List of failures removed.\n\n", format_in=GREEN, format_out=self.terminal_default)

        elif line in SEARCH_SYNOS:
            if not self.dbbuddy.search_terms:
                _stdout("Search terms list is already empty.\n\n", format_in=RED, format_out=self.terminal_default)
            else:
                confirm = br.ask("%sAre you sure you want to delete all %s search terms (y/[n])?%s " %
                                 (RED, len(self.dbbuddy.search_terms), self.terminal_default), default="no")

                if not confirm:
                    _stdout("Aborted...\n", format_in=RED, format_out=self.terminal_default)
                else:
                    self.dbbuddy.search_terms = []
                    _stdout("Search terms removed.\n\n", format_in=GREEN, format_out=self.terminal_default)

        elif line in TRASH_SYNOS:
            if not self.dbbuddy.trash_bin:
                _stdout("Trash bin is already empty.\n\n", format_in=RED, format_out=self.terminal_default)
            else:
                confirm = br.ask("%sAre you sure you want to delete all %s records from your trash bin (y/[n])?%s " %
                                 (RED, len(self.dbbuddy.trash_bin), self.terminal_default), default="no")

                if not confirm:
                    _stdout("Aborted...\n", format_in=RED, format_out=self.terminal_default)
                else:
                    self.dbbuddy.trash_bin = OrderedDict()
                    _stdout("Trash bin emptied.\n\n", format_in=GREEN, format_out=self.terminal_default)

        elif line in RECORD_SYNOS:
            if not self.dbbuddy.records:
                _stdout("Records list is already empty.\n", format_in=RED, format_out=self.terminal_default)
            else:
                confirm = br.ask("%sAre you sure you want to delete all %s records from your main filtered list "
                                 "(y/[n])?%s " % (RED, len(self.dbbuddy.records), self.terminal_default), default="no")
                if not confirm:
                    _stdout("Aborted...\n", format_in=RED, format_out=self.terminal_default)
                else:
                    self.dbbuddy.records = OrderedDict()
                    _stdout("All records removed from main list (trash bin is still intact).\n\n",
                            format_in=GREEN, format_out=self.terminal_default)

        else:
            confirm = br.ask("%sAre you sure you want to completely reset your live session (y/[n])?%s " %
                             (RED, self.terminal_default), default="no")

            if not confirm:
                _stdout("Aborted...\n", format_in=RED, format_out=self.terminal_default)
            else:
                self.dbbuddy.trash_bin = OrderedDict()
                self.dbbuddy.records = OrderedDict()
                self.dbbuddy.search_terms = []
                self.dbbuddy.failures = OrderedDict()
                _stdout("Live session cleared of all data.\n\n", format_in=GREEN, format_out=self.terminal_default)

        self.dump_session()

    def do_failures(self, *_):
        if not self.dbbuddy.failures:
            _stdout("No failures to report\n\n", format_in=GREEN, format_out=self.terminal_default)
        else:
            _stdout("The following failures have occured\n", format_in=[UNDERLINE, GREEN],
                    format_out=self.terminal_default)
            for _hash, _values in self.dbbuddy.failures.items():
                _stdout("%s\n\n" % _values, format_out=self.terminal_default)

    def do_fetch(self, *_):
        accn_only = self.dbbuddy.record_breakdown()["accession"]
        if accn_only:
            search_terms = list(self.dbbuddy.search_terms)
            self.dbbuddy.search_terms = []
            retrieve_summary(self.dbbuddy)
            self.dbbuddy.search_terms = search_terms

        amount_seq_requested = 0
        new_records_fetched = []
        for _accn, _rec in self.dbbuddy.records.items():
            if not _rec.record and _rec.size:  # Not fetching sequence if the full record already exists
                amount_seq_requested += _rec.size
                new_records_fetched.append(_accn)

        if amount_seq_requested > 5000000:
            confirm = br.ask("{0}You are requesting {2}{1}{0} residues of sequence data. "
                             "Continue (y/[n])?{3}".format(GREEN, br.pretty_number(amount_seq_requested),
                                                           YELLOW, self.terminal_default), default="no")
            if not confirm:
                _stdout("Aborted...\n\n", format_in=RED, format_out=self.terminal_default)
                return

        retrieve_sequences(self.dbbuddy)
        seq_retrieved = 0
        for _accn in new_records_fetched:
            if self.dbbuddy.records[_accn].record:
                seq_retrieved += self.dbbuddy.records[_accn].size

        _stdout("Retrieved %s residues of sequence data\n\n" % br.pretty_number(seq_retrieved),
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
                dbbuddy = dill.load(ifile)
                self.dbbuddy.search_terms = dbbuddy.search_terms
                self.dbbuddy.records = dbbuddy.records
                self.dbbuddy.trash_bin = dbbuddy.trash_bin
                self.dbbuddy.out_format = dbbuddy.out_format
                self.dbbuddy.failures = dbbuddy.failures
                self.dbbuddy.databases = dbbuddy.databases
                self.dbbuddy.memory_footprint = dbbuddy.memory_footprint

            for _db, client in self.dbbuddy.server_clients.items():
                if client:
                    client.http_errors_file = br.TempFile()
                    client.results_file = br.TempFile()

            _stdout("Session loaded from file.\n\n", format_in=GREEN, format_out=self.terminal_default, quiet=quiet)
            self.dump_session()
        except (EOFError, IOError):
            _stdout("Error: Unable to read the provided file. Are you sure it's a saved DbBuddy live session?\n\n",
                    format_in=RED, format_out=self.terminal_default)

    def do_keep(self, line=None):
        self.filter(line, mode="keep")

    def do_quit(self, *_):
        if (self.dbbuddy.records or self.dbbuddy.trash_bin) and self.hash != hash(self.dbbuddy):
            confirm = br.ask("You have unsaved records, are you sure you want to quit (y/[n])?", default="no")
            if not confirm:
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

        # Create directory if necessary
        line = os.path.abspath(line)
        _dir, _file = os.path.split(line)
        try:
            os.makedirs(_dir, exist_ok=True)
        except PermissionError:
            _stdout("Error: You do not have write privileges to create a directory in the specified path.\n\n",
                    format_in=RED, format_out=self.terminal_default)
            return

        # Set the .db extension
        if not re.search("\.db$", line):
            line += ".db"

        # Warn if file exists
        if os.path.isfile(line):
            confirm = br.ask("%sFile already exists, overwrite [y]/n?%s " % (RED, self.terminal_default))
            if not confirm:
                _stdout("Abort...\n\n", format_in=RED, format_out=self.terminal_default)
                return
        try:
            open(line, "wb").close()
        except PermissionError:
            _stdout("Error: You do not have write privileges to create a file in the specified directory.\n\n",
                    format_in=RED, format_out=self.terminal_default)
            return

        self.crash_file.save(line)
        self.hash = hash(self.dbbuddy)
        _stdout("Live session saved\n\n", format_in=GREEN, format_out=self.terminal_default)

    def do_search(self, line):
        if not line:
            line = input("%sSpecify search string(s):%s " % (RED, self.terminal_default))

        # Do this on a temp dbbuddy obj so searches are not repeated
        temp_buddy = DbBuddy(line)
        temp_buddy.databases = self.dbbuddy.databases
        retrieve_summary(temp_buddy)
        for _term in temp_buddy.search_terms:
            if _term not in self.dbbuddy.search_terms:
                self.dbbuddy.search_terms.append(_term)

        for _accn, _rec in temp_buddy.records.items():
            if _accn not in self.dbbuddy.records:
                self.dbbuddy.records[_accn] = _rec

        for _hash, failure in temp_buddy.failures.items():
            _stdout("%s\n" % failure, format_in=RED, format_out=self.terminal_default)
            if _hash not in self.dbbuddy.failures:
                self.dbbuddy.failures[_hash] = failure
        self.dump_session()

    def do_show(self, line=None, group="records"):
        line = [] if not line else line.split(" ")

        # Note that trashbin is only shown with the command 'trash' from the UI
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
                br._stderr("%sWarning: %s records are only summary data, so will not be displayed in %s format. "
                           "Use 'fetch' to retrieve all sequence data.%s\n"
                           % (RED, len(breakdown["summary"] + breakdown["accession"]),
                              self.dbbuddy.out_format, self.terminal_default))

            num_records = len(breakdown["full"])

        columns = []
        force_num_records = 0
        for _next in line:
            try:
                force_num_records = int(_next)
            except ValueError:
                columns.append(_next)

        columns = None if not columns else columns

        if num_records > 100 and not force_num_records:
            confirm = br.ask("%sShow all %s records (y/[n])?%s " % (RED, num_records, self.terminal_default), False)
            if not confirm:
                _stdout("Include an integer value with 'show' to return a specific number of records.\n\n",
                        format_out=self.terminal_default)
                return
        try:
            num_records = num_records if not force_num_records else force_num_records
            self.dbbuddy.print(_num=num_records, columns=columns, group=group)
        except ValueError as _e:
            if "Sequences must all be the same length" in str(_e):
                _stdout("Error: '%s' format does not support sequences of different length." %
                        self.dbbuddy.out_format, format_in=RED, format_out=self.terminal_default)
            elif "No suitable quality scores found in letter_annotations of SeqRecord" in str(_e):
                _stdout("Error: BioPython requires quality scores to output in '%s' format, and this data is not "
                        "currently available to DatabaseBuddy." % self.dbbuddy.out_format,
                        format_in=RED, format_out=self.terminal_default)
            else:
                raise ValueError(_e)

        br._stderr("%s\n" % self.terminal_default)

    def do_sort(self, line=None):
        def sub_sort(records, _sort_columns, _rev=False):
            heading = _sort_columns[0]
            subgroups = OrderedDict()
            for accn, _rec in records.items():
                if heading.lower() == "accn":
                    subgroups.setdefault(_rec.accession, {})
                    subgroups[_rec.accession][accn] = _rec
                elif heading.lower() == "type":
                    subgroups.setdefault(_rec.type, {})
                    subgroups[_rec.type][accn] = _rec
                elif heading.lower() == "db":
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

                try:  # If the column is numbers sort numerically, otherwise alphabetically
                    subgroups = OrderedDict(sorted(subgroups.items(), key=lambda _x: int(_x[0]), reverse=_rev))
                except ValueError:
                    subgroups = OrderedDict(sorted(subgroups.items(), key=lambda _x: _x[0], reverse=_rev))

            final_order = OrderedDict()
            for subgroup, _recs in subgroups.items():
                if len(_recs) == 1:
                    _recs = [_rec for accn, _rec in _recs.items()]
                    final_order[_recs[0].accession] = _recs[0]
                else:
                    if len(_sort_columns) > 1:
                        _recs = sub_sort(_recs, _sort_columns[1:], _rev)

                    for x, y in _recs.items():
                        final_order[x] = y
            return final_order

        line = "ACCN" if not line else line
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

    def do_status(self, *_):
        _stdout("%s\n" % str(self.dbbuddy), format_out=self.terminal_default)

    def do_write(self, line=None):
        if not line:
            line = input("%sWhere would you like your records written?%s " % (RED, self.terminal_default))

        # Ensure the specified directory exists
        line = os.path.abspath(line)
        _dir, _file = os.path.split(line)
        if not os.path.isdir(_dir):
            _stdout("The specified directory does not exist. Please create it before continuing "
                    "(you can use the 'bash' command from within the DbBuddy Live Session).\n\n", format_in=RED,
                    format_out=self.terminal_default)
            return

        # Warn if file exists
        if os.path.isfile(line):
            confirm = br.ask("%sFile already exists, overwrite [y]/n?%s " % (RED, self.terminal_default))
            if not confirm:
                _stdout("Abort...\n\n", format_in=RED, format_out=self.terminal_default)
                return

        breakdown = self.dbbuddy.record_breakdown()
        if self.dbbuddy.out_format in ["ids", "accessions", "summary", "full-summary"]:
            if breakdown["full"]:
                confirm = br.ask("%sYou are about to write to a summary format "
                                 "which does not include sequence. Continue [y]/n?%s" % (RED, self.terminal_default))
                if not confirm:
                    _stdout("Abort...\n", format_in=RED, format_out=self.terminal_default)
                    return
            count = len(breakdown["full"] + breakdown["summary"])
            msg = "accession" if self.dbbuddy.out_format in ["ids", "accessions"] else "summary record"
            msg += "s" if count > 1 else ""
            msg = "%s %s" % (count, msg)
        else:
            non_full = len(breakdown["summary"] + breakdown["accession"])
            if non_full > 0:
                _stdout('''\
NOTE: There are %s summary records in the Live Session, and only full records can be written
  in '%s' format. Use the 'fetch' command to retrieve full records.
''' % (non_full, self.dbbuddy.out_format), format_in=RED, format_out=self.terminal_default)
            msg = "%s %s record" % (len(breakdown["full"]), self.dbbuddy.out_format)
            msg += "s" if len(breakdown["full"]) > 1 else ""

        _stdout("%s written to %s.\n\n" % (msg, line), format_in=GREEN,
                format_out=self.terminal_default)
        self.hash = hash(self.dbbuddy)
        try:
            ofile = open(line, "w", encoding="utf-8")
            self.dbbuddy.print(quiet=True, destination=ofile)
            ofile.close()
        except PermissionError:
            _stdout("Error: You do not have write privileges in the specified directory.\n\n",
                    format_in=RED, format_out=self.terminal_default)
        return

    def do_undo(self, *_):
        if not self.undo:
            _stdout("There is currently no undo history (only a single undo is possible).\n\n",
                    format_in=RED, format_out=self.terminal_default)
            return

        self.do_load("%s_undo" % self.crash_file.path, quiet=True)
        self.undo = False
        _stdout("Most recent state reloaded\n\n", format_in=GREEN, format_out=self.terminal_default)

    def complete_bash(self, *args):
        # ToDo: Windows support
        text = args[0]
        if not self.shell_execs:
            path_dirs = Popen("echo $PATH", stdout=PIPE, shell=True).communicate()
            path_dirs = path_dirs[0].decode("utf-8").split(":")
            for _dir in path_dirs:
                if not os.path.isdir(_dir):
                    continue
                root, dirs, files = next(os.walk(_dir))
                for _file in files:
                    if os.access("%s%s%s" % (root, os.path.sep, _file), os.X_OK):
                        self.shell_execs.append(_file.strip())
        return ["%s " % x for x in self.shell_execs if x.startswith(text)]

    @staticmethod
    def complete_database(*args):
        text = args[0]
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
        return sorted(completions)

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
        return sorted(completions)

    def complete_show(self, *args):
        text = args[0]
        return ["%s " % x for x in self.get_headings() if x.lower().startswith(text.lower())]

    def complete_sort(self, *args):
        text = args[0]
        return ["%s " % x for x in self.get_headings() + ["reverse"] if x.lower().startswith(text.lower())]

    def complete_write(self, *args):
        line, startidx, endidx = args[1:]
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
        return sorted(completions)

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
    - To make searches case-insensitive, prefix the filter with 'i?'. E.g., 'i?(organism)HuMaN' or i?hUmAn.. \n
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
    - To make searches case-insensitive, prefix the filter with 'i?'. E.g., 'i?(organism)HuMaN' or i?hUmAn.\n
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
    - To make searches case-insensitive, prefix the filter with 'i?'. E.g., 'i?(organism)HuMaN' or i?hUmAn.\n
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
the number of records and amount of information per record displayed.
Use a negative integer to return records from the bottom of the list.\n
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
    br.check_garbage_flags(in_args, "DatabaseBuddy")

    dbbuddy = []
    out_format = "summary" if not in_args.out_format else in_args.out_format

    if isinstance(in_args.user_input[0], TextIOWrapper) and in_args.user_input[0].buffer.raw.isatty():
            dbbuddy = DbBuddy()
            in_args.live_shell = True
    elif len(in_args.user_input) > 1:
        for search_set in in_args.user_input:
            dbbuddy.append(DbBuddy(search_set, in_args.database, out_format))

        dbbuddy = DbBuddy(dbbuddy, in_args.database, out_format)
    else:
        dbbuddy = DbBuddy(in_args.user_input[0], in_args.database, out_format)

    return in_args, dbbuddy


def command_line_ui(in_args, dbbuddy, skip_exit=False):
    # ############################################## COMMAND LINE LOGIC ############################################## #
    # Live Shell
    temp_file = br.TempFile(byte_mode=True)

    def launch_live_shell():
        # Create a temp file for crash handling
        temp_file.open()
        try:  # Catch all exceptions and try to send error report to server
            LiveShell(dbbuddy, temp_file)
        except SystemExit:
            pass
        except br.GuessError as err:
            print(err)
        except Exception as err:
            save_file = ".%sDbSessionDump_%s" % (os.sep, temp_file.name)
            temp_file.save(save_file)
            br.send_traceback("DatabaseBuddy", "live_shell", err, VERSION)
            br._stderr("\n%sYour work has been saved to %s, and can be loaded by launching DatabaseBuddy and using "
                       "the 'load' command.%s\n" % (GREEN, save_file, DEF_FONT))
        dbbuddy.memory_footprint = int(os.path.getsize(temp_file.path))
        _exit("LiveShell")

    def _exit(tool, skip=skip_exit):
        if skip:
            return
        usage = br.Usage()
        usage.increment("DatabaseBuddy", VERSION.short(), tool, dbbuddy.memory_footprint)
        usage.save()
        temp_file.close()
        sys.exit()

    if in_args.live_shell:
        launch_live_shell()

    """
    def _print_recs(_dbbuddy):
        _dbbuddy.print()

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
            br._stderr("%s\n# ################################################################ #\n\n" % output.strip())

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

    # Default to LiveShell
    launch_live_shell()


def main():
    br.preparse_flags()
    initiation = []
    try:
        initiation = argparse_init()
        command_line_ui(*initiation)
    except br.GuessError as e:
        print(e)
        return False
    except SystemExit:
        return False
    except Exception as e:
        function = ""
        for next_arg in vars(initiation[0]):
            if getattr(initiation[0], next_arg) and next_arg in br.db_flags:
                function = next_arg
                break
        br.send_traceback("DbBuddy", function, e, VERSION)
        return False
    return True

if __name__ == '__main__':
    main()
