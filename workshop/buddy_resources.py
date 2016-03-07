#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: buddy_resources.py
version: 1.0
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Collection of resources used by all BuddySuite tools,
             including dictionaries of the commands available for each Buddy tool
"""
from __future__ import print_function
import sys
if sys.version_info[0] < 3:
    print("Error: Attempting to run BuddySuite with Python %s.%s. Python 3 required." %
          (sys.version_info[0], sys.version_info[1]))
    sys.exit()
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
from Bio import AlignIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


# ##################################################### CLASSES ###################################################### #
class GuessError(Exception):
    """Raised when input format cannot be guessed"""
    def __init__(self, _value):
        self.value = _value

    def __str__(self):
        return self.value


class PhylipError(Exception):
    """Raised when phylip format is malformed"""
    def __init__(self, _value):
        self.value = _value

    def __str__(self):
        return self.value


class Contributor(object):
    def __init__(self, first, last, middle="", commits=None, github=None):
        self.first = first.strip()
        self.middle = middle.strip()
        self.last = last.strip()
        self.commits = commits
        self.github = github

    def name(self):
        _name = " ".join([self.first, self.middle, self.last])
        _name = re.sub("  ", " ", _name)  # in case there is no middle name
        return _name

    def __str__(self):
        _output = "%s %s" % (self.first, self.last)
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


class Usage(object):
    # ToDo: Check internet connectivity!!!
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


class Version(object):
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
        name_line_length = 0
        for contributor in _contributors:
            if len(contributor.name()) > name_line_length:
                name_line_length = len(contributor.name())

        _output = ""
        for contributor in _contributors:
            _output += "%s%s\n" % (contributor.name().ljust(name_line_length + 2), contributor.github)
        return _output.strip()

    def short(self):
        return "%s.%s" % (self.major, self.minor)

    def __str__(self):
        _output = '''\
%s %s.%s (%s)

Public Domain Notice
This is free software; see the source for detailed copying conditions.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
Questions/comments/concerns can be directed to Steve Bond, steve.bond@nih.gov

Contributors:
%s
''' % (self.name, self.major, self.minor, self.release_date, self.contributors_string())
        return _output


# #################################################### FUNCTIONS ##################################################### #
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
    :param _modifiers: dict of single letter flags that have more global significance
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


def parse_format(_format):
    available_formats = ["clustal", "embl", "fasta", "genbank", "gb", "nexus", "stockholm",
                         "phylip", "phylipis", "phylip-strict", "phylip-interleaved-strict",
                         "phylipi", "phylip-relaxed", "phylip-interleaved", "phylipr",
                         "phylips", "phylipsr", "phylip-sequential", "phylip-sequential-relaxed",
                         "phylipss", "phylip-sequential-strict", "nexml", "newick"]

    _format = _format.lower()
    if _format in ["phylip", "phylipis", "phylip-strict", "phylip-interleaved-strict"]:
        return "phylip"

    if _format in ["phylipi", "phylip-relaxed", "phylip-interleaved", "phylipr"]:
        return "phylip-relaxed"

    if _format in ["phylips", "phylipsr", "phylip-sequential", "phylip-sequential-relaxed"]:
        return "phylipsr"

    if _format in ["phylipss", "phylip-sequential-strict"]:
        return "phylipss"

    if _format not in available_formats:
        raise TypeError("Format type '%s' is not recognized/supported" % _format)

    return _format


def phylip_sequential_out(_input, relaxed=True, _type="alignbuddy"):
    output = ""
    if _type == "alignbuddy":
        alignments = _input.alignments
    else:
        alignments = [_input.records]

    for alignment in alignments:
        ids = []
        id_check = []
        aln_len = 0
        for rec in alignment:
            if rec.id in id_check:
                raise PhylipError("Malformed Phylip --> Repeat id '%s'" % rec.id)
            id_check.append(rec.id)
            if not aln_len:
                aln_len = len(str(rec.seq))

        max_id_len = 0
        for rec in alignment:
            if len(str(rec.seq)) != aln_len:
                raise PhylipError("Malformed Phylip --> The length of record '%s' is incorrect" % rec.id)
            max_id_len = len(rec.id) if len(rec.id) > max_id_len else max_id_len

        output += " %s %s" % (len(alignment), aln_len)
        for rec in alignment:
            if relaxed:
                seq_id = re.sub('[ \t]+', '_', rec.id)
                output += "\n%s%s" % (seq_id.ljust(max_id_len + 2), rec.seq)
            else:
                seq_id = rec.id[:10].ljust(10)
                output += "\n%s%s" % (seq_id, rec.seq)

            if seq_id in ids:
                raise PhylipError("Malformed Phylip --> Repeat id '%s' after strict truncation. "
                                  "Try a relaxed Phylip format (phylipr or phylipsr)." % seq_id)
            ids.append(seq_id)

        output += "\n\n"
    return output


def phylip_sequential_read(sequence, relaxed=True):
    sequence = "\n %s" % sequence.strip()
    sequence = re.sub("\n+", "\n", sequence)
    alignments = re.split("\n *([0-9]+) ([0-9]+)\n", sequence)[1:]
    align_dict = OrderedDict()
    for indx in range(int(len(alignments) / 3)):
        align_dict[(int(alignments[indx * 3]), int(alignments[indx * 3 + 1]), indx)] = alignments[indx * 3 + 2]

    temp_file = TempFile()
    aligns = []
    for key, seqs in align_dict.items():
        records = []
        seqs = re.sub("[\n\t]", " ", seqs).strip()
        while seqs != "":
            if not relaxed:
                _id = seqs[:10]
                seqs = seqs[10:]
            else:
                _id = re.match("([^ ]+) +", seqs).group(1)
                seqs = re.sub("[^ ]+ +", "", seqs, count=1)
            rec = ""
            while len(rec) < key[1]:
                breakdown = re.match("([^ ]+)", seqs)
                if not breakdown:
                    raise PhylipError("Malformed Phylip --> Less sequence found than expected")
                rec += breakdown.group(0)
                if re.match("[^ ]+$", seqs):
                    seqs = ""
                seqs = re.sub("[^ ]+ +", "", seqs, count=1)

            records.append((_id, rec))

        if len(records) != key[0]:
            raise PhylipError("Malformed Phylip --> %s sequences expected, %s found." % (key[0], len(records)))

        key_list = []
        output = ""
        for seq_id, seq in records:
            if key[1] != len(seq):
                raise PhylipError("Malformed Phylip --> Sequence %s has %s columns, %s expected." %
                                  (seq_id, len(seq), key[1]))
            if seq_id in key_list:
                if relaxed:
                    raise PhylipError("Malformed Phylip --> Repeat ID %s." % seq_id)
                else:
                    raise PhylipError("Malformed Phylip --> Repeat id '%s' after strict truncation. "
                                      "Try a relaxed Phylip format (phylipr or phylipsr)." % seq_id)
            key_list.append(seq_id)
            output += ">%s\n%s\n" % (seq_id, seq)
        temp_file.write(output, "w")
        aligns.append(AlignIO.read(temp_file.get_handle("r"), "fasta"))
        temp_file.close()
    return aligns


def replacements(input_str, query, replace="", num=0):
    """
    This will allow fancy positional regular expression replacements from left-to-right, as well as normal right-to-left
    :param input_str:
    :param query:
    :param replace:
    :param num:
    :return:
    """
    check_parentheses = re.findall("\([^()]*\)", query)
    check_replacement = re.findall(r"\\[0-9]+", replace)
    check_replacement = sorted([int(match[1:]) for match in check_replacement])
    if check_replacement and check_replacement[-1] > len(check_parentheses):
        raise AttributeError("There are more replacement match values specified than query parenthesized groups")

    if num < 0:
        if check_replacement:
            for indx in sorted(range(check_replacement[-1]), reverse=True):
                indx += 1
                replace = re.sub(r"\\%s" % indx, r"\\%s" % (indx + 1), replace)
            right_replace = "\\%s" % (len(check_replacement) + 2)
        else:
            right_replace = "\\2"
        leftmost = str(input_str)
        new_str = str(input_str)
        rightmost = ""
        hash_to_split_on = "UPNFSZ7FQ6RBhfFzwt0Cku4Yr1n2VvwVUG7x97G7"
        for _ in range(abs(num)):
            if leftmost == "":
                break
            new_str = re.sub(r"(.*)%s(.*)" % query,
                             r"\1%s%s%s" % (hash_to_split_on, replace, right_replace), leftmost, 1)
            new_str = new_str.split(hash_to_split_on)
            if len(new_str) == 2:
                leftmost = new_str[0]
                rightmost = new_str[1] + rightmost
                new_str = leftmost + rightmost
            else:
                new_str = leftmost + rightmost
                break
    else:
        new_str = re.sub(query, replace, input_str, num)

    return new_str


def send_traceback(tool, e):
    # ToDo: Explicitly state the tool being called in the ErrorReport. It's not always obvious...
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


def shift_features(features, shift, full_seq_len):
    """
    Adjust the location of features
    :param features: Either a single SeqFeature object, or a list of them
    :param shift: int, how far the new feature should move from 0
    :param full_seq_len: The full length of the original sequence
    :return: List of SeqFeatures
    """
    if type(features) != list:  # Duck type for single feature input
        features = [features]

    shifted_features = []
    for feature in features:
        if type(feature.location) == CompoundLocation:  # Recursively call shift_features() for compound locations
            new_compound_location = []
            for sub_feature in feature.location.parts:
                sub_feature = shift_features(SeqFeature(sub_feature), shift, full_seq_len)
                if not sub_feature:
                    continue
                new_compound_location.append(sub_feature[0].location)

            if not new_compound_location:
                continue

            elif len(new_compound_location) == 1:
                feature.location = new_compound_location[0]

            else:
                feature.location = CompoundLocation(new_compound_location, feature.location.operator)

        elif type(feature.location) == FeatureLocation:
            start = feature.location.start + shift
            end = feature.location.end + shift
            if start > full_seq_len or end < 0:
                continue

            start = start if start >= 0 else 0
            end = end if end <= full_seq_len else full_seq_len

            feature.location = FeatureLocation(start, end, feature.strand)

        else:
            raise TypeError("_shift_feature requires a feature with either FeatureLocation or CompoundLocation, "
                            "not %s" % type(feature.location))
        shifted_features.append(feature)

    return shifted_features


def ungap_feature_ends(feat, rec):
    """
    If a feature begins or ends on a gap, it makes it much harder to track changes, so force the feature onto actual
    residues.
    :param feat: Either a FeatureLocation or CompoundLocation object
    :param rec: The original SeqRecord that the feature is derived from
    :return: The modified feature object
    """
    if feat.location.start < 0:
        feat.location = FeatureLocation(0, feat.location.end, feat.location.strand)

    if feat.location.end < 0:
        feat.location = FeatureLocation(feat.location.start, 0, feat.location.strand)

    if feat.location.start > feat.location.end:
            feat.location = FeatureLocation(feat.location.end, feat.location.start, feat.location.strand)

    if type(feat.location) == CompoundLocation:
        parts = []
        for part in feat.location.parts:
            part = ungap_feature_ends(SeqFeature(part), rec)
            parts.append(part.location)
        feat.location = CompoundLocation(parts, feat.location.operator)

    elif type(feat.location) == FeatureLocation:
        extract = str(feat.extract(rec.seq))
        front_gaps = re.search("^-+", extract)
        if front_gaps:
            if not feat.location.strand or feat.location.strand == 1:
                new_start = feat.location.start + len(front_gaps.group(0))
                feat.location = FeatureLocation(new_start, feat.location.end, 1)
            else:
                new_end = feat.location.end - len(front_gaps.group(0))
                feat.location = FeatureLocation(feat.location.start, new_end, -1)

        rear_gaps = re.search("-+$", extract)
        if rear_gaps:
            if not feat.location.strand or feat.location.strand == 1:
                new_end = feat.location.end - len(rear_gaps.group(0))
                feat.location = FeatureLocation(feat.location.start, new_end, 1)
            else:
                new_start = feat.location.start + len(rear_gaps.group(0))
                feat.location = FeatureLocation(new_start, feat.location.end, -1)
    else:
        raise TypeError("FeatureLocation or CompoundLocation object required.")
    return feat


def _old2new(feat, old_rec, new_rec):
    if feat.location.start == feat.location.end == 0:
        return feat

    if type(feat.location) == CompoundLocation:
        parts = []
        for part in feat.location.parts:
            new_part = _old2new(SeqFeature(part), old_rec, new_rec)
            if new_part:
                parts.append(new_part.location)
        if len(parts) == 1:
            feat.location = parts[0]
        elif len(parts) > 1:
            feat.location = CompoundLocation(parts, feat.location.operator)
        else:
            return None
    elif type(feat.location) == FeatureLocation:
        if feat.location.start > feat.location.end:
            start, end = feat.location.end, feat.location.start
        else:
            start, end = feat.location.start, feat.location.end
        old_seq = str(old_rec.seq).lower()
        new_seq = str(new_rec.seq).lower()
        old_front_seq = old_seq[:start]
        old_front_seq = re.sub("-", "", old_front_seq)
        old_feat_seq = old_seq[start:end]
        old_feat_seq = re.sub("-", "", old_feat_seq)
        start, end = 0, 0
        new_front_seq, new_feat_seq = "", ""
        for indx, residue in enumerate(new_seq):
            if residue == "-":
                continue

            if not start:
                if old_front_seq in ["", new_front_seq]:
                    start = indx + 1
                    new_feat_seq += residue
                    if new_feat_seq == old_feat_seq:
                        end = indx + 1
                        break
                else:
                    new_front_seq += residue
            else:
                new_feat_seq += residue
                if new_feat_seq == old_feat_seq:
                    end = indx + 1
                    break
        start -= 1
        if start == -1:
            return None
        end = end if end != 0 else len(new_seq)
        feat.location = FeatureLocation(start, end, feat.location.strand)
    else:
        raise TypeError("FeatureLocation or CompoundLocation object required.")
    return feat


def remap_gapped_features(old_records, new_records):
    """
    If adding, subtracting, or moving around in a sequence, the features need to be shifted to accomodate.
    This only works if all of the original non-gap residues are present in the new record
    :param old_records: Starting sequence (can be gapped as well)
    :param new_records: New sequence with different gap pattern
    :return:
    """
    # Start by forcing feature start-end positions onto actual residues, in cases were they fall on gaps
    for old_rec, new_rec in zip(old_records, new_records):
        features = []
        for feat in old_rec.features:
            features.append(ungap_feature_ends(feat, old_rec))
        old_rec.features = features
        features = []
        for feat in old_rec.features:
            feat = _old2new(feat, old_rec, new_rec)
            if feat:
                features.append(feat)
        new_rec.features = features
    return new_records


# #################################################### VARIABLES ##################################################### #

contributors = [Contributor("Stephen", "Bond", commits=497, github="https://github.com/biologyguy"),
                Contributor("Karl", "Keat", commits=299, github="https://github.com/KarlKeat"),
                Contributor("Jeremy", "Labarge", commits=25, github="https://github.com/biojerm")]

# NOTE: If this is added to, be sure to update the unit test!
format_to_extension = {'fasta': 'fa', 'fa': 'fa', 'genbank': 'gb', 'gb': 'gb', 'newick': 'nwk', 'nwk': 'nwk',
                       'nexus': 'nex', 'nex': 'nex', 'phylip': 'phy', 'phy': 'phy', 'phylip-relaxed': 'phyr',
                       'phyr': 'phyr', 'phylipss': 'physs', 'physs': 'physs', 'phylipsr': 'physr',
                       'physr': 'physr', 'stockholm': 'stklm', 'stklm': 'stklm', 'clustal': 'clus', 'clus': 'clus'}


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
                      "action": "append",
                      "nargs": "+",
                      "metavar": ("subject", "<blast params>"),
                      "help": "Search a BLAST database or subject sequence file with your query sequence file, "
                              "returning the full hits"},
            "clean_seq": {"flag": "cs",
                          "action": "append",
                          "nargs": "*",
                          "metavar": "args",
                          "help": "Strip out non-sequence characters. Args: ['strict'][replace char]. "
                                  "'strict' removes all ambiguous letters, and optionally choose the replacement "
                                  "character (default=N)"},
            "complement": {"flag": "cmp",
                           "action": "store_true",
                           "help": "Return complement of nucleotide sequence"},
            "concat_seqs": {"flag": "cts",
                            "action": "append",
                            "nargs": "?",
                            "metavar": "'clean'",
                            "help": "Concatenate a bunch of sequences into a single solid string. Pass in "
                                    "the word 'clean' to remove stops, gaps, etc., from the sequences "
                                    "before concatenating"},
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
            "degenerate_sequence": {"flag": "dgn",
                                    "action": "append",
                                    "nargs": '?',
                                    "type": int,
                                    "help": "Convert unambiguous codons to degenerate codons"},
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
                               "metavar": ("[columns (int)]", "[scope (all|ids|seqs)]"),
                               "help": "Strip repeat records (ids and/or identical sequences). "
                                       "Defaults: 1 'all'"},
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
            "group_by_prefix": {"flag": "gbp",
                                "action": "append",
                                "nargs": "*",
                                "metavar": "args",
                                "help": "Sort sequences into separate files based on prefix. "
                                        "args: [Split Pattern [Split pattern ...]] [length (int)] [out dir]"},
            "group_by_regex": {"flag": "gbr",
                               "action": "append",
                               "nargs": "*",
                               "metavar": "args",
                               "help": "Sort sequences into separate files based on regular expression. "
                                       "args: <regex> [regex [...]] [out dir]"},
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
            "make_ids_unique": {"flag": "miu",
                                "action": "append",
                                "nargs": "*",
                                "metavar": ("<separator(string)>", "<padding(int)>"),
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
            "rename_ids": {"flag": "ri",
                           "action": "append",
                           "metavar": "args",
                           "nargs": "*",
                           "help": "Replace some pattern in ids with something else. "
                                   "args: <pattern>, <substitution>, [max replacements (int)], ['store']"},
            "replace_subseq": {"flag": "rs",
                               "action": "append",
                               "metavar": "args",
                               "nargs": "+",
                               "help": "Replace some pattern in sequences with something else. "
                                       "args: <query (regex)> [query ...] [replacement]"},
            "reverse_complement": {"flag": "rc",
                                   "action": "store_true",
                                   "help": "Return reverse complement of nucleotide sequence"},
            "reverse_transcribe": {"flag": "r2d",
                                   "action": "store_true",
                                   "help": "Convert RNA sequences to DNA"},
            "screw_formats": {"flag": "sf",
                              "action": "store",
                              "metavar": "<out_format>",
                              "help": "Change the file format to something else"},
            "select_frame": {"flag": "sfr",
                             "action": "store",
                             "metavar": "<frame (int)>",
                             "type": int,
                             "choices": [1, 2, 3],
                             "help": "Change the reading frame of nucleotide sequences"},
            "shuffle_seqs": {"flag": "ss",
                             "action": "store_true",
                             "help": "Randomly rearrange the residues in each record"},
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
                               "metavar": "",
                               "action": "store",
                               "help": "If you want a specific format output"},
                "params": {"flag": "p",
                           "action": "store",
                           "metavar": "",
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
             "clean_seq": {"flag": "cs",
                           "action": "append",
                           "nargs": "*",
                           "metavar": "args",
                           "help": "Strip out non-alignment characters. Args: ['strict'][replace char]. "
                                   "'strict' removes all ambiguous letters, and optionally choose the replacement "
                                   "character (default=N)"},
             "concat_alignments": {"flag": "cta",
                                   "action": "append",
                                   "nargs": "*",
                                   "metavar": "regex|int",
                                   "help": "Concatenates two or more alignments using a regex pattern or fixed length "
                                           "prefix to group record ids."},
             "consensus": {"flag": "con",
                           "action": "store_true",
                           "help": "Create majority-rule consensus sequences"},
             "delete_records": {"flag": "dr",
                                "nargs": "+",
                                "action": "append",
                                "metavar": "regex",
                                "help": "Remove alignment rows with IDs that contain matches to the provided patterns"},
             "enforce_triplets": {"flag": "et",
                                  "action": "store_true",
                                  "help": "Shift gaps so sequences are organized in triplets"},
             "extract_range": {"flag": "er",
                               "action": "store",
                               "nargs": 2,
                               "metavar": ("<start (int)>", "<end (int)>"),
                               "type": int,
                               "help": "Pull out sub-alignments in a given range"},
             "generate_alignment": {"flag": "ga",
                                    "action": "append",
                                    "nargs": "*",
                                    "metavar": "args",
                                    "help": "Create a new alignment from unaligned sequences. "
                                            "args: [alignment program] [optional params]"},
             "hash_ids": {"flag": "hi",
                          "action": "append",
                          "nargs": "?",
                          "type": int,
                          "metavar": "hash length (int)",
                          "help": "Rename all sequence IDs to fixed length hashes. Default length is 10."},
             "list_ids": {"flag": "li",
                          "action": "append",
                          "nargs": "?",
                          "type": int,
                          "metavar": "int",
                          "help": "Output the sequence identifiers. Optionally, pass in an integer to "
                                  "specify the # of columns to write"},
             "lowercase": {"flag": "lc",
                           "action": "store_true",
                           "help": "Convert all sequences to lowercase"},
             "mapfeat2align": {"flag": "mf2a",
                               "action": "store",
                               "nargs": "+",
                               "metavar": "<unaligned file>",
                               "help": "Transfer features from annotated sequences over to an alignment."},
             "num_seqs": {"flag": "ns",
                          "action": "store_true",
                          "help": "Count how many sequences are present in each alignment"},
             "order_ids": {"flag": "oi",
                           "action": "append",
                           "nargs": "?",
                           "choices": ["rev"],
                           "help": "Sort all sequences in an alignment by id in alpha-numeric order. "
                                   "Pass in the word 'rev' to reverse order"},
             "pull_records": {"flag": "pr",
                              "nargs": "+",
                              "action": "append",
                              "metavar": "regex",
                              "help": "Keep alignment rows with IDs that contains matches to the provided patterns"},
             "rename_ids": {"flag": "ri",
                            "action": "append",
                            "metavar": "args",
                            "nargs": "*",
                            "help": "Replace some pattern in ids with something else. "
                                    "args: <pattern>, <substitution>, [max replacements (int)]"},
             "reverse_transcribe": {"flag": "r2d",
                                    "action": "store_true",
                                    "help": "Convert RNA alignments to DNA"},
             "screw_formats": {"flag": "sf",
                               "action": "store",
                               "metavar": "<out_format>",
                               "help": "Change the file format to something else"},
             "split_to_files": {"flag": "stf",
                                "action": "append",
                                "nargs": "+",
                                "metavar": "args",
                                "help": "Write individual files for each alignment. Args: <out_dir> [prefix]"},
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
            "hash_ids": {"flag": "hi",
                         "action": "append",
                         "nargs": "*",
                         "metavar": "args",
                         "help": "Rename all taxon label IDs (and optionally inner node lables) to fixed length hashes."
                                 " args: [hash length (int)] ['nodes']"},
            "num_tips": {"flag": "nt",
                         "action": "store_true",
                         "help": "Display the number of tips in each tree"},
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


if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(prog='buddy_resources')

    main_parser.add_argument('-v', '--version', help='Show module version #s', action='store_true')
    main_parser.add_argument('-t', '--tools', help="List all tools", action='store_true')

    in_args = main_parser.parse_args()

    if in_args.version:
        import SeqBuddy
        import AlignBuddy
        import PhyloBuddy
        import DatabaseBuddy
        print("SeqBuddy: %s" % SeqBuddy.VERSION.short())
        print("AlignBuddy: %s" % AlignBuddy.VERSION.short())
        print("PhyloBuddy: %s" % PhyloBuddy.VERSION.short())
        print("DatabaseBuddy: %s" % DatabaseBuddy.VERSION.short())
        sys.exit(datetime.datetime.strptime(str(datetime.date.today()), '%Y-%m-%d'))

    if in_args.tools:
        print("### SeqBuddy")
        for key in sb_flags:
            print(key)
        print("\n### AlignBuddy")
        for key in alb_flags:
            print(key)
        print("\n### PhyloBuddy")
        for key in pb_flags:
            print(key)
        print("\n### DatabaseBuddy")
        for key in db_flags:
            print(key)
        sys.exit()
