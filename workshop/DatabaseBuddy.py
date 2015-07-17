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
import time

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
        self.records = {}
        self.out_format = _out_format.lower()

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
                db = guess_database(_accession)
                if db:
                    self.records[_accession] = Record(_accession, db) if not _database else Record(_accession, _database)

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
            _output = "%s/n" % _output.strip()
        return _output


class Record:
    def __init__(self, _accession, _database, _seq=None, _meta=None, _type=None):
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


# ################################################# Database Clients ################################################# #
class EnsemblRestClient:
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, headers=None, params=None):
        if headers is None:
            {}

        if 'Content-Type' not in headers:
            headers['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = urllib.Request(self.server + endpoint, headers=hdrs)
            response = urllib.urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except urllib.HTTPError, e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    def get_variants(self, species, symbol):
        genes = self.perform_rest_action(
            '/xrefs/symbol/{0}/{1}'.format(species, symbol),
            params={'object_type': 'gene'}
        )
        if genes:
            stable_id = genes[0]['id']
            variants = self.perform_rest_action(
                '/overlap/id/{0}'.format(stable_id),
                params={'feature': 'variation'}
            )
            return variants
        return None
# #################################################################################################################### #


def pull_ensembl_records(_dbbuddy):
    search_recs = [_rec for _accession, _rec in _dbbuddy.records.items() if _rec.database == "ensembl"]
    ids = '"%s' % '", "'.join([_rec.accession for _rec in search_recs])
    ids += '"'
    headers = {"Content-Type": "application/json"}
    for _rec in search_recs:
        #request = requests.get("http://rest.ensembl.org/archive/id/%s" % _rec.accession, headers=headers)
        request = requests.get("http://rest.ensembl.org/archive/id/%s?type=protein" % _rec.accession, headers=headers)
        #request = requests.get("http://rest.ensembl.org/sequence/id/%s" % _rec.accession, headers=headers)
        #request = requests.get("http://rest.ensembl.org/sequence/id/%s?type=protein" % _rec.accession, headers=headers)
        #request = requests.get("http://rest.ensembl.org/sequence/id/%s?type=protein" % _rec.accession, headers=headers)

        if not request.ok:
            request.raise_for_status()

        decoded = request.json()
        print(decoded)



def search_ensembl(_dbbuddy):
    for search_term in _dbbuddy.search_terms:
        url = "http://metazoa.ensembl.org/Multi/Search/Results?q=%s;species=all;collection=all;site=ensemblunit"\
              % search_term
        content = requests.get(url).text

        soup = BeautifulSoup(content, 'html.parser')
        try:
            paginate = soup.find('div', {"class": 'paginate'}).find_all('a')
            max_page = 1
            for page in paginate:
                try:
                    if int(page.text) > max_page:
                        max_page = int(page.text)
                except ValueError:
                    continue
        except AttributeError:
            max_page = 1

        print("%s pages of results were returned" % max_page)
        ids = {}
        for page_num in range(max_page):
            stdout.write("\rCollecting data from results page %s" % (page_num + 1),)
            stdout.flush()

            url = "http://metazoa.ensembl.org/Multi/Search/Results?page=%s;q=%s;species=all;collection=all;site=ensemblunit"\
                  % (page_num + 1, in_args.search_term)
            content = requests.get(url).text

            soup = BeautifulSoup(content)

            for row in soup.find_all('div', {"class": 'row'}):
                sub_soup = BeautifulSoup(str(row))
                lhs = sub_soup.find('div', {"class": "lhs"}).text
                rhs = sub_soup.find('div', {"class": "rhs"}).text

                if lhs == "Gene ID":
                    gene_id = rhs

                if lhs == "Species":
                    if rhs in ids:
                        ids[rhs].append(gene_id)
                    else:
                        ids[rhs] = [gene_id]

        if len(ids) == 0:
            exit("\rNo records found for query '%s'" % in_args.search_term)

        output = ""
        for species in ids:
            output += "%s\n" % species
            for next_id in ids[species]:
                output += "%s\n" % next_id
            output += "\n"

        if in_args.outfile:
            outfile = os.path.abspath(in_args.outfile)
            with open(outfile, "w") as ofile:
                ofile.write(output)
            print("Output written to %s" % outfile)

        else:
            print(output)
    return _dbbuddy


def download_everything(_dbbuddy):
    pull_ensembl_records(_dbbuddy)
    return

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
        download_everything(dbbuddy)

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
