#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: buddy_resources.py
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
if int(sys.version_info[0]) < 3:
    print("Error: Attempting to run BuddySuite with Python %s. Python 3+ required." % sys.version_info[0])
    sys.exit()
import argparse
import datetime
from collections import OrderedDict
import os
from configparser import ConfigParser, NoOptionError
import json
import traceback
import re
import sre_compile
from ftplib import FTP, all_errors
from hashlib import md5
from urllib import request
from urllib.error import URLError, HTTPError, ContentTooShortError
from multiprocessing import Process, cpu_count
from time import time, sleep
from math import floor
from tempfile import TemporaryDirectory
from shutil import copytree, rmtree, copyfile
import string
from random import choice
import signal
from pkg_resources import Requirement, resource_filename, DistributionNotFound

from Bio import AlignIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC


# ################################################## MYFUNCS ################################################### #
class Timer(object):
    def __init__(self):
        self.current_time = round(time())

    def start(self):
        self.current_time = round(time())
        return

    def end(self):
        return pretty_time(round(time()) - self.current_time)


class RunTime(object):
    def __init__(self, prefix=None, postfix=None, out_type="stdout", _sleep=0, final_clear=False):
        """
        Sets up a dynamic counter that lets the user know how long a job has been going for
        :param prefix: Some arbitrary text to go in front of the counter message
        :param postfix:  Some arbitrary text to go behind the counter message
        :param out_type: Either "stdout" or "stderr". Cannot be sys.stdout/err directly, for pickling reasons
        :param _sleep: Pause the loop for however many seconds
        :param final_clear: If set to True, the counter message will be deleted before moving on
        """
        if out_type not in ["stdout", "stderr"]:
            raise ValueError("The 'out_type' parameter must be either 'stdout' or 'stderr', not %s." % out_type)
        self.out_type = out_type
        self.prefix = prefix if prefix else ""
        self.postfix = postfix if postfix else ""
        self.running_process = None
        self.sleep = _sleep
        self.final_clear = final_clear

    def _run(self, check_file_path):
        out_type = sys.stdout if self.out_type == "stdout" else sys.stderr
        d_print = DynamicPrint(out_type)
        start_time = round(time())
        elapsed = 0
        while True:
            prefix = self.prefix if not hasattr(self.prefix, '__call__') else self.prefix()
            postfix = self.postfix if not hasattr(self.postfix, '__call__') else self.postfix()
            with open("%s" % check_file_path, "r") as ifile:
                if round(time()) - start_time == elapsed:
                    continue
                elif ifile.read() == "Running":
                    d_print.write("%s%s%s" % (prefix, pretty_time(elapsed), postfix))
                    elapsed = round(time()) - start_time
                else:
                    if not self.final_clear:
                        d_print.write("%s%s%s\n" % (prefix, pretty_time(elapsed), postfix))
                    else:
                        d_print.clear()
                    break
            sleep(self.sleep)
        return

    def start(self):
        if self.running_process:
            self.end()
        tmp_file = TempFile()
        tmp_file.write("Running")
        p = Process(target=self._run, args=(tmp_file.path,))
        p.daemon = 1
        p.start()
        self.running_process = [tmp_file, p]
        return

    def end(self):
        if not self.running_process:
            return
        self.running_process[0].clear()
        while self.running_process[1].is_alive():
            continue
        self.running_process = None
        return


# maybe use curses library in the future to extend this for multi-line printing
class DynamicPrint(object):
    def __init__(self, out_type="stdout", quiet=False):
        """
        :param out_type: 'stdout', 'stderr', sys.stdout, or sys.stderr
        :param quiet:
        """
        self._last_print = ""
        self._next_print = ""
        self._writer = self._write()

        out_type = sys.stdout if out_type == "stdout" else out_type
        out_type = sys.stderr if out_type == "stderr" else out_type
        self.out_type = out_type
        self.quiet = quiet
        return

    def _write(self):
        try:
            while True:
                self.out_type.write("\r%s\r%s" % (" " * len(self._last_print), self._next_print),)
                self.out_type.flush()
                self._last_print = self._next_print
                yield
        finally:
            pass

    def write(self, content):
        if not self.quiet:
            content = re.sub("\t", "    ", content)
            self._next_print = content
            next(self._writer)
        return

    def new_line(self, number=1):
        if not self.quiet:
            self.out_type.write("\n" * number)
            self.out_type.flush()
            self._last_print = ""
        return

    def clear(self):
        self.write("")
        return


def pretty_time(seconds):
    if seconds < 60:
        output = "%i sec" % seconds
    elif seconds < 3600:
        minutes = floor(seconds / 60)
        seconds -= minutes * 60
        output = "%i min, %i sec" % (minutes, seconds)
    elif seconds < 86400:
        hours = floor((seconds / 60) / 60)
        seconds -= hours * 60 * 60
        minutes = floor(seconds / 60)
        seconds -= minutes * 60
        output = "%i hrs, %i min, %i sec" % (hours, minutes, seconds)
    else:
        days = floor(((seconds / 60) / 60) / 24)
        seconds -= (days * 60 * 60 * 24)
        hours = floor((seconds / 60) / 60)
        seconds -= (hours * 60 * 60)
        minutes = floor(seconds / 60)
        seconds -= (minutes * 60)
        output = "%i days, %i hrs, %i min, %i sec" % (days, hours, minutes, seconds)

    return output


def pretty_number(num, mode='short', precision=2):  # mode in ['short', 'medium', 'long']
    magnitude = 0
    while abs(num) >= 1000 and magnitude < 8:
        magnitude += 1
        num = round(num / 1000.0, precision)
    if mode == 'short':
        return ('%s %s' % (num, ['', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'][magnitude])).strip()
    elif mode == 'medium':
        return ('%s %s' % (num, ['', 'Kilo', 'Mega', 'Giga', 'Tera', 'Peta', 'Exa', 'Zetta',
                                 'Yotta'][magnitude])).strip()
    elif mode == 'long':
        return ('%s %s' % (num, ['', 'Thousand', 'Million', 'Billion', 'Trillion', 'Quadrillion', 'Quintillion',
                                 'Sextillion', 'Septillion'][magnitude])).strip()
    else:
        raise ValueError("Valid 'mode' values are 'short', 'medium', and 'long'")


def usable_cpu_count():
    cpus = cpu_count()
    if cpus > 7:
        max_processes = cpus - 3
    elif cpus > 3:
        max_processes = cpus - 2
    elif cpus > 1:
        max_processes = cpus - 1
    else:
        max_processes = 1

    return max_processes


def run_multicore_function(iterable, function, func_args=False, max_processes=0, quiet=False, out_type=sys.stdout):
        # fun little piece of abstraction here... directly pass in a function that is going to be looped over, and
        # fork those loops onto independent processes. Any arguments the function needs must be provided as a list.
        if func_args and not isinstance(func_args, list):
            raise AttributeError("The arguments passed into the multi-thread function must be provided as a list")

        d_print = DynamicPrint(out_type)
        if max_processes == 0:
            max_processes = usable_cpu_count()

        else:
            cpus = cpu_count()
            if max_processes > cpus:
                max_processes = cpus
            elif max_processes < 1:
                max_processes = 1

        max_processes = max_processes if max_processes < len(iterable) else len(iterable)

        running_processes = 0
        child_list = []
        start_time = round(time())
        elapsed = 0
        counter = 0
        if not quiet:
            d_print.write("Running function %s() on %s cores\n" % (function.__name__, max_processes))
        # fire up the multi-core!!
        if not quiet:
            d_print.write("\tJob 0 of %s" % len(iterable))

        for next_iter in iterable:
            if type(iterable) is dict:
                next_iter = iterable[next_iter]
            if os.name == "nt":  # Multicore doesn't work well on Windows, so for now just run serial
                function(next_iter, func_args)
                continue
            while 1:     # Only fork a new process when there is a free processor.
                if running_processes < max_processes:
                    # Start new process
                    if not quiet:
                        d_print.write("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

                    if func_args:
                        p = Process(target=function, args=(next_iter, func_args))
                    else:
                        p = Process(target=function, args=(next_iter,))
                    p.start()
                    child_list.append(p)
                    running_processes += 1
                    counter += 1
                    break
                else:
                    # processor wait loop
                    while 1:
                        for i in range(len(child_list)):
                            if child_list[i].is_alive():
                                continue
                            else:
                                child_list.pop(i)
                                running_processes -= 1
                                break

                        if not quiet:
                            if (start_time + elapsed) < round(time()):
                                elapsed = round(time()) - start_time
                                d_print.write("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

                        if running_processes < max_processes:
                            break

        # wait for remaining processes to complete --> this is the same code as the processor wait loop above
        if not quiet:
            d_print.write("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

        while len(child_list) > 0:
            for i in range(len(child_list)):
                if child_list[i].is_alive():
                    continue
                else:
                    child_list.pop(i)
                    running_processes -= 1
                    break  # need to break out of the for-loop, because the child_list index is changed by pop

            if not quiet:
                if (start_time + elapsed) < round(time()):
                    elapsed = round(time()) - start_time
                    d_print.write("\t%s total jobs (%s, %s jobs remaining)" % (len(iterable), pretty_time(elapsed),
                                                                               len(child_list)))

        if not quiet:
            d_print.write("\tDONE: %s jobs in %s\n" % (len(iterable), pretty_time(elapsed)))
        # func_args = []  # This may be necessary because of weirdness in assignment of incoming arguments
        return


class TempDir(object):
    def __init__(self):
        self.dir = next(self._make_dir())
        self.path = self.dir.name
        self.subdirs = []
        self.subfiles = []

    def _make_dir(self):
        tmp_dir = TemporaryDirectory()
        yield tmp_dir
        rmtree(self.path)

    def subdir(self, dir_name=None):
        if not dir_name:
            dir_name = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])
            while dir_name in self.subdirs:  # Catch the very unlikely case that a duplicate occurs
                dir_name = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])

        subdir_path = "%s%s%s" % (self.path, os.path.sep, dir_name)
        if not os.path.exists(subdir_path):
            os.mkdir(subdir_path)
        if dir_name not in self.subdirs:
            self.subdirs.append(dir_name)
        return subdir_path

    def del_subdir(self, _dir):
        path, _dir = os.path.split(_dir)
        del self.subdirs[self.subdirs.index(_dir)]
        rmtree("%s%s%s" % (self.path, os.path.sep, _dir))
        return

    def subfile(self, file_name=None):
        if not file_name:
            file_name = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])
            root, dirs, files = next(walklevel(self.path))
            while file_name in files:  # Catch the very unlikely case that a duplicate occurs
                file_name = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])

        open("%s%s%s" % (self.path, os.path.sep, file_name), "w", encoding="utf-8").close()
        self.subfiles.append(file_name)
        return "%s%s%s" % (self.path, os.path.sep, file_name)

    def del_subfile(self, _file):
        path, _file = os.path.split(_file)
        del self.subfiles[self.subfiles.index(_file)]
        os.remove("%s%s%s" % (self.path, os.path.sep, _file))
        return

    def save(self, location, keep_hash=False):
        location = location if not keep_hash else "%s%s%s" % (location, os.path.sep, os.path.split(self.path)[-1])
        if os.path.isdir(location):
            print("Save Error: Indicated output folder already exists in TempDir.save(%s)" % location, file=sys.stderr)
            return False
        else:
            copytree(self.dir.name, location)
            return True


class TempFile(object):
    # I really don't like the behavior of tempfile.[Named]TemporaryFile(), so hack TemporaryDirectory() via TempDir()
    def __init__(self, mode="w", byte_mode=False):
        self._tmp_dir = TempDir()  # This needs to be a persistent (ie self.) variable, or the directory will be deleted
        path, dir_hash = os.path.split(self._tmp_dir.path)
        self.name = dir_hash
        self.path = "%s%s%s" % (self._tmp_dir.path, os.path.sep, dir_hash)
        open(self.path, "w", encoding="utf-8").close()
        self.handle = None
        self.bm = "b" if byte_mode else ""
        self.mode = mode

    def open(self, mode=None):
        mode = "%s%s" % (self.mode, self.bm) if not mode else "%s%s" % (mode, self.bm)
        if self.handle:
            self.close()
        self.handle = open(self.path, mode)

    def close(self):
        if self.handle:
            self.handle.close()
            self.handle = None

    def get_handle(self, mode=None):
        self.open(self.mode) if not mode else self.open(mode)
        return self.handle

    def write(self, content, mode="a"):
        mode = "%s%s" % (mode, self.bm)
        if mode not in ["w", "wb", "a", "ab"]:
            print("Write Error: mode must be 'w' or 'a' in TempFile.write()", file=sys.stderr)
            return False
        already_open = True if self.handle else False
        if not already_open:
            self.open(mode[0])
        if mode in ["a", "ab"]:
            self.handle.write(content)
        else:
            self.handle.truncate(0)
            self.handle.write(content)
        if not already_open:
            self.close()
        return True

    def read(self):
        already_open = True if self.handle else False
        position = 0
        if already_open:
            position = self.handle.tell()
            self.close()
        with open(self.path, "r%s" % self.bm) as ifile:
            content = ifile.read()
        if already_open:
            self.open(mode="a")
            self.handle.seek(position)
        return content

    def clear(self):
        self.close()
        content = "" if self.bm == "" else b""
        self.write(content, mode="w")
        return

    def save(self, location):
        with open(location, "w%s" % self.bm) as ofile:
            ofile.write(self.read())
        return


class SafetyValve(object):  # Use this class if you're afraid of an infinite loop
    def __init__(self, global_reps=1000, state_reps=10, counter=0):
        self.counter = counter

        self._start_global_reps = global_reps
        self.global_reps = global_reps

        self._start_state_reps = state_reps
        self.state_reps = state_reps

        self.state = ""

    def step(self, message=""):  # step() is general, and doesn't care whether useful computation is on going
        self.global_reps -= 1
        self.counter += 1
        if self.global_reps == 0:
            raise RuntimeError("You just popped your global_reps safety valve. %s" % message)
        return True

    def test(self, state, message=""):  # test() keeps track of some variable 'state' to see if its value keeps changing
        if self.state == str(state):
            self.state_reps -= 1
        else:
            self.state_reps = self._start_state_reps
            self.state = str(state)

        if self.state_reps == 0:
            raise RuntimeError("You just popped your state_reps safety valve. %s" % message)
        return True


# Pulled this function off of Stack Overflow -- posted by nosklo
# Iterates over directories only to a specified depth (useful in a for loop)
# Note that this is a generator, so need to use next() or `with` to get a result
def walklevel(some_dir, level=1):
    some_dir = some_dir.rstrip(os.path.sep)
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]


def copydir(source, dest):
    def search(_dir):
        contents = os.listdir(_dir)
        files = []
        dirs = []
        for thing in contents:
            thing_path = "%s%s%s" % (_dir, os.path.sep, thing)
            if os.path.isdir(thing_path):
                dirs.append(thing_path)
            else:
                files.append(thing_path)
        if len(dirs) != 0:
            for _dir in dirs:
                files += search(_dir)
        return files

    file_paths = search(source)
    if not os.path.exists(dest):
        os.makedirs(dest)
    for path in file_paths:
        _path, _file = os.path.split(path)
        copyfile(path, "%s%s%s" % (dest, os.path.sep, _file))


def ask(input_prompt, default="yes", timeout=0):
    if default == "yes":
        yes_list = ["yes", "y", '']
        no_list = ["no", "n", "abort"]
    else:
        yes_list = ["yes", "y"]
        no_list = ["no", "n", "abort", '']

    def kill(*args):
        raise TimeoutError(args)

    try:
        if os.name == "nt":
            import msvcrt
            timeout = timeout if timeout > 0 else 3600
            start_time = time()
            sys.stdout.write(input_prompt)
            sys.stdout.flush()
            _response = ''
            while True:
                if msvcrt.kbhit():
                    _chr = msvcrt.getche()
                    if ord(_chr) == 13:  # enter_key
                        break
                    elif ord(_chr) >= 32:  # space_char
                        _response += _chr.decode()
                if len(_response) == 0 and (time() - start_time) > timeout:
                    _response = default
                    break

            print('')  # needed to move to next line
            if _response.lower() in yes_list:
                return True
            else:
                return False

        else:
            signal.signal(signal.SIGALRM, kill)
            signal.alarm(timeout)
            _response = input(input_prompt)
            signal.alarm(0)
            while True:
                if _response.lower() in yes_list:
                    return True
                elif _response.lower() in no_list:
                    return False
                else:
                    print("Response not understood. Valid options are 'yes' and 'no'.")
                    signal.alarm(timeout)
                    _response = input(input_prompt)
                    signal.alarm(0)

    except TimeoutError:
        return False


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
    def __init__(self):
        self.tmpfile = TempFile()
        self.config = config_values()
        usage_file = None
        if self.config["diagnostics"] and self.config["data_dir"]:
            usage_file = "%s%sbuddysuite_usage.json" % (self.config["data_dir"], os.path.sep)
            try:
                if not os.path.isfile(usage_file):
                    with open(usage_file, "w", encoding="utf-8") as ofile:
                        ofile.write("{}")

                with open(usage_file, "r", encoding="utf-8") as ifile:
                    self.stats = json.load(ifile)
                    if not self.stats:  # Empty file needs to be populated with a little info
                        self.clear_stats()

            except (PermissionError, ValueError):
                usage_file = None

        if not usage_file:
            self.tmpfile.write("{}")
            self.clear_stats()
            usage_file = self.tmpfile.path

        self.usage_file_path = usage_file

    def clear_stats(self):
        self.stats = {"user_hash": self.config["user_hash"]}

    def increment(self, buddy, version, tool, obj_size=None):
        self.stats.setdefault(buddy, {})
        self.stats[buddy].setdefault(version, {})
        self.stats[buddy][version].setdefault(tool, 0)
        self.stats[buddy][version][tool] += 1
        if obj_size:
            self.stats[buddy][version].setdefault("sizes", [])
            self.stats[buddy][version]["sizes"].append(obj_size)
        return

    def save(self, send_report=True):
        if self.config["diagnostics"] and send_report:
            self.stats.setdefault("last_upload", datetime.date.today().isoformat())
            if (datetime.datetime.today() - datetime.datetime.strptime(self.stats["last_upload"],
                                                                       '%Y-%m-%d')).days >= 7:
                self.send_report()
                return
        try:
            with open(self.usage_file_path, "w", encoding="utf-8") as ofile:
                json.dump(self.stats, ofile)
        except PermissionError:
            pass
        return

    def send_report(self):
        self.stats["date"] = str(datetime.date.today())
        temp_file = TempFile()
        json.dump(self.stats, temp_file.get_handle())
        try:
            ftp = FTP("rf-cloning.org", user="buddysuite", passwd="seqbuddy", timeout=1)
            ftp.storlines("STOR usage_%s" % temp_file.name, temp_file.get_handle("rb"))
            self.clear_stats()
            self.stats["last_upload"] = str(datetime.date.today())
        except all_errors as e:
            if "timed out" in str(e):
                return
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
            self.release_date = datetime.date(**release_date)

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
    options = {"email": "buddysuite@nih.gov",
               "diagnostics": False,
               "user_hash": "hashless",
               "shortcuts": ""}
    try:
        config_file = resource_filename(Requirement.parse("buddysuite"),
                                        "buddysuite{0}buddy_data{0}config.ini".format(os.path.sep))
        config = ConfigParser()
        config.read(config_file)
        for _key, value in options.items():
            try:
                if _key in ['diagnostics']:
                    options[_key] = config.getboolean('DEFAULT', _key)
                else:
                    options[_key] = config.get('DEFAULT', _key)
            except KeyError:
                options[_key] = value
        options["shortcuts"] = options["shortcuts"].split(",")
        options["data_dir"] = resource_filename(Requirement.parse("buddysuite"), "buddysuite%sbuddy_data" % os.path.sep)
        if not os.path.isdir(options["data_dir"]):
            options["data_dir"] = False
    except (DistributionNotFound, KeyError, NoOptionError):  # This occurs when buddysuite isn't installed
        options["data_dir"] = False
    return options


def check_garbage_flags(in_args, tool):
    """
    If an unknown flag is thrown immediately after the sequence/alignment/tree positional argument it is not treated
    as a flag, leading to a GuessError. Catch it here and exit more gracefully.
    :param in_args: parser.parse_args object
    :param tool: Which BuddySuite tool is calling?
    :return: None
    """
    flag_list = []
    if tool == "AlignBuddy":
        flag_list = in_args.alignments
    if tool == "DatabaseBuddy":
        flag_list = in_args.user_input
    if tool == "PhyloBuddy":
        flag_list = in_args.trees
    if tool == "SeqBuddy":
        flag_list = in_args.sequence

    for flag in flag_list:
        if flag and re.match(" -", str(flag)):
            _stderr("%s.py: error: unrecognized arguments: %s\n" % (tool, flag))
            sys.exit()
    return True


def error_report(trace_back, permission=False):
    message = ""
    error_hash = re.sub("^#.*?\n{2}", "", trace_back, flags=re.DOTALL)  # Remove error header information before hashing
    error_hash = md5(error_hash.encode("utf-8")).hexdigest()  # Hash the error
    try:  # Check online to see if error has been reported before
        raw_error_data = request.urlopen("https://raw.githubusercontent.com/biologyguy/BuddySuite/error_codes/"
                                         "diagnostics/error_codes", timeout=2)
        error_string = raw_error_data.read().decode("utf-8")  # Read downloaded file
        error_string = re.sub("#.*\n", "", error_string)
        error_json = json.loads(error_string)  # Convert JSON into a data table

        version_str = re.search("# [A-Z]?[a-z]+Buddy: (.*)", trace_back).group(1)

        if error_hash in error_json.keys():  # Check if error is known (if it's in the data table)
            if error_json[error_hash][1] == "None" or error_json[error_hash][1] == version_str:  # If error not resolved
                message += "This is a known bug since version %s, " \
                           "but it has not been resolved yet.\n" % error_json[error_hash][0]

            else:  # If error has been resolved
                print("This bug was resolved in version %s. We recommend you upgrade to the latest version (if you "
                      "downloaded BuddySuite using pip, use the command pip install "
                      "buddysuite --upgrade).\n" % error_json[error_hash][1])
                return

        else:  # If error is unknown
            message += "Uh oh, you've found a new bug! This issue is not currently in our bug tracker.\n"

    except (URLError, HTTPError, ContentTooShortError) as err:  # If there is an error, just blow through
        message += "Failed to locate known error codes:\n%s\n" % str(err)

    if permission:
        message += "An error report with the above traceback is being sent to the BuddySuite developers because " \
                   "you have elected to participate in the Software Improvement Program. You may opt-out of this " \
                   "program at any time by re-running the BuddySuite installer.\n"
        print(message)
    else:
        permission = ask("%s\nAn error report with the above traceback has been prepared and is ready to send to the "
                         "BuddySuite developers.\nWould you like to upload the report? [y]/n " % message, timeout=15)
    try:
        if permission:
            print("\nPreparing error report for FTP upload...")
            temp_file = TempFile()
            temp_file.write(trace_back)
            print("Connecting to FTP server...")
            ftp = FTP("rf-cloning.org", user="buddysuite", passwd="seqbuddy", timeout=5)
            print("Sending...")
            ftp.storlines("STOR error_%s" % temp_file.name, open(temp_file.path, "rb"))  # Upload error to FTP
            print("Success! Thank you.")
    except all_errors as e:
            print("Well... We tried. Seems there was a problem with the FTP upload\n%s" % e)
    return


def preparse_flags():
    func_lists = [sb_flags, alb_flags, pb_flags, db_flags, sb_modifiers, pb_modifiers, alb_modifiers, db_modifiers]
    for indx, arg in enumerate(sys.argv):
        if indx > 0 and re.match('^-', arg):
            change = True
            for func_list in func_lists:
                for func in func_list:
                    if arg == "-" + func_list[func]["flag"] or arg == "--" + func:
                        change = False
            if arg in ["-v", "-t", "-c", "-h", "--version", "--test", "--help"]:
                change = False
            if change:
                sys.argv[indx] = " " + sys.argv[indx]


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
        for func, _in_args in _flags.items():
            args = ("-%s" % _in_args["flag"], "--%s" % func)
            kwargs = {}
            for cmd, val in _in_args.items():
                if cmd == 'flag':
                    continue
                kwargs[cmd] = val
            parser_flags.add_argument(*args, **kwargs)

    if _modifiers:
        _modifiers = OrderedDict(sorted(_modifiers.items(), key=lambda x: x[0]))
        parser_modifiers = parser.add_argument_group(title="\033[1mModifying options\033[m")
        for func, _in_args in _modifiers.items():
            args = ("-%s" % _in_args["flag"], "--%s" % func)
            kwargs = {}
            for cmd, val in _in_args.items():
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
    # If your file is not phylip-relaxed, leaving relaxed as True WILL break your code.
    # If your file is strict you must set relaxed to False.
    # (Strict forces 10 character taxa names, relaxed requires whitespace between name and sequence)
    sequence = "\n %s" % sequence.strip()
    sequence = re.sub("\n+", "\n", sequence)
    alignments = re.split("\n *([0-9]+) ([0-9]+)\n", sequence)[1:]
    align_dict = OrderedDict()
    for indx in range(int(len(alignments) / 3)):
        align_dict[(int(alignments[indx * 3]), int(alignments[indx * 3 + 1]), indx)] = alignments[indx * 3 + 2]

    temp_file = TempFile()
    aligns = []
    for _key, seqs in align_dict.items():
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
            while len(rec) < _key[1]:
                breakdown = re.match("([^ ]+)", seqs)
                if not breakdown:
                    raise PhylipError("Malformed Phylip --> Less sequence found than expected")
                rec += breakdown.group(0)
                if re.match("[^ ]+$", seqs):
                    seqs = ""
                seqs = re.sub("[^ ]+ +", "", seqs, count=1)

            records.append((_id, rec))

        if len(records) != _key[0]:
            raise PhylipError("Malformed Phylip --> %s sequences expected, %s found." % (_key[0], len(records)))

        key_list = []
        output = ""
        for seq_id, seq in records:
            if _key[1] != len(seq):
                raise PhylipError("Malformed Phylip --> Sequence %s has %s columns, %s expected." %
                                  (seq_id, len(seq), _key[1]))
            if seq_id in key_list:
                if relaxed:
                    raise PhylipError("Malformed Phylip --> Repeat ID %s." % seq_id)
                else:
                    raise PhylipError("Malformed Phylip --> Repeat id '%s' after strict truncation. "
                                      "Try a relaxed Phylip format (phylipr or phylipsr)." % seq_id)
            key_list.append(seq_id)
            output += ">%s\n%s\n" % (seq_id, seq)
        with open(temp_file.path, "w", encoding="utf-8") as ofile:
            ofile.write(output)
        with open(temp_file.path, "r", encoding="utf-8") as ifile:
            aligns.append(AlignIO.read(ifile, "fasta"))
    return aligns


def phylip_guess(next_format, _input):
    if next_format == "phylip":
        sequence = "\n %s" % _input.read().strip()
        _input.seek(0)
        alignments = re.split("\n ([0-9]+) ([0-9]+)\n", sequence)[1:]
        align_sizes = []
        for indx in range(int(len(alignments) / 3)):
            align_sizes.append((int(alignments[indx * 3]), int(alignments[indx * 3 + 1])))

        phy = list(AlignIO.parse(_input, "phylip"))
        _input.seek(0)

        indx = 0
        phy_ids = []
        for key in align_sizes:
            phy_ids.append([])
            for rec in phy[indx]:
                assert len(rec.seq) == key[1]
                phy_ids[-1].append(rec.id)
            indx += 1

        phy_rel = list(AlignIO.parse(_input, "phylip-relaxed"))
        _input.seek(0)
        if phy_rel:
            for indx, aln in enumerate(phy_rel):
                for rec in aln:
                    if len(rec.seq) != align_sizes[indx][1]:
                        return parse_format("phylip")
                    if rec.id in phy_ids[indx]:
                        return
                    else:
                        return parse_format("phylip-relaxed")
        return parse_format("phylip")

    if next_format == "phylipss":
        if phylip_sequential_read(_input.read(), relaxed=False):
            _input.seek(0)
            return parse_format(next_format)
        else:
            return
    if next_format == "phylipsr":
        if phylip_sequential_read(_input.read()):
            _input.seek(0)
            return parse_format(next_format)
        else:
            return


def replacements(input_str, query, replace="", num=0):
    """
    This will allow fancy positional regular expression replacements from right-to-left, as well as normal left-to-right
    :param input_str: The string that replacements will be working on
    :param query: Regular expression
    :param replace: What to substitute the matches with
    :param num: Number of matches to replace. 0 = All, +ve nums from left-to-right, -ve nums from right-to-left
    :return: Modified string
    """
    # First make sure the user isn't trying to access more replacement groups than specified with parentheses
    check_parentheses = re.findall("\([^()]*\)", query)
    check_keep_group = re.findall(r"\\[0-9]+", replace)
    check_keep_group = sorted([int(match[1:]) for match in check_keep_group])
    if check_keep_group and check_keep_group[-1] > len(check_parentheses):
        raise AttributeError("There are more replacement match values specified than query parenthesized groups")

    if num < 0:  # Make replacements from right-to-left
        matches = re.finditer(query, input_str)
        matches = list(matches)
        matches.reverse()
        num = abs(num) if abs(num) <= len(matches) else len(matches)
        for match_indx in range(num):
            groups = matches[match_indx].groups()
            next_replace = str(replace)
            for group_indx, group in enumerate(groups):
                group_indx += 1
                next_replace = re.sub(r'\\%s' % group_indx, group, next_replace)
            start = matches[match_indx].start()
            end = matches[match_indx].end()
            input_str = input_str[:start] + next_replace + input_str[end:]
    else:
        input_str = re.sub(query, replace, input_str, num)

    return input_str


def send_traceback(tool, function, e, version):
    now = datetime.datetime.now()
    config = config_values()
    tb = ""
    for _line in traceback.format_tb(sys.exc_info()[2]):
        if os.name == "nt":
            _line = re.sub('"(?:[A-Za-z]:)*\{0}.*\{0}(.*)?"'.format(os.sep), r'"\1"', _line)
        else:
            _line = re.sub('"{0}.*{0}(.*)?"'.format(os.sep), r'"\1"', _line)
        tb += _line
    bs_version = "# %s: %s\n" % (tool, version.short())
    func = "# Function: %s\n" % function
    platform = "# Platform: %s\n" % sys.platform
    python = "# Python: %s\n" % re.sub("[\n\r]", "", sys.version)
    user = "# User: %s\n" % config['user_hash']
    date = "# Date: %s\n\n" % now.strftime('%Y-%m-%d')
    error = "%s: %s\n\n" % (type(e).__name__, e)

    tb = "".join([bs_version, func, python, platform, user, date, error, tb])
    print("\033[m%s::%s has crashed with the following traceback:\033[91m\n\n%s\n\n\033[m" % (tool, function, tb))
    error_report(tb, config["diagnostics"])
    return


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
        if feat.strand == -1 and rec.seq.alphabet == IUPAC.protein:
            feat.strand = None
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
        if feat.location.start == feat.location.end == 0:
            return feat
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
            return None  # I don't think this should ever be hit
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
        new_rec.annotations = old_rec.annotations
        new_rec.dbxrefs = old_rec.dbxrefs
    return new_records


def _stderr(message, quiet=False):
    """
    Send text to stderr
    :param message: Text to write
    :param quiet: Suppress message with True
    :return: None
    """
    if not quiet:
        try:
            sys.stderr.write(message)
            sys.stderr.flush()
        except UnicodeEncodeError:  # Mainly a work around for Windows
            message = message.encode("utf-8")
            sys.stderr.buffer.write(message)
            sys.stderr.flush()
    return


def _stdout(message, quiet=False):
    """
    Send text to stdout
    :param message: Text to write
    :param quiet: Suppress message with True
    :return: None
    """
    if not quiet:
        try:
            sys.stdout.write(message)
            sys.stdout.flush()
        except UnicodeEncodeError:  # Mainly a work around for Windows
            message = message.encode("utf-8")
            sys.stdout.buffer.write(message)
            sys.stdout.flush()
    return


def utf_encode(_input):
    tmp_file = TempFile()
    with open(tmp_file.path, "w", encoding="utf-8") as ofile:
        ofile.write(_input)
    import codecs
    with codecs.open(tmp_file.path, "r", "utf-8", errors="replace") as ifile:
        _input = ifile.read()
    _input = re.sub("\r", "", _input)
    return _input


def isfile_override(path):
    # This is a hack for Windows, which throws errors if the 'path' is too long
    import stat
    try:
        st = os.stat(path)
    except OSError:
        return False
    except ValueError as err:
        if "path too long for Windows" in str(err):
            return False
        else:
            raise err
    return stat.S_ISREG(st.st_mode)

if os.name == "nt":
    os.path.isfile = isfile_override


def clean_regex(patterns, quiet=False):
    """
    Ensure that user provided regular expression are valid
    :param patterns: either a single str regex or list of regexes
    :param quiet: Suppress stderr
    :return:
    """
    patterns = [patterns] if type(patterns) == str else patterns
    failures = []
    failure_str = ""
    for indx, regex_test in enumerate(patterns):
        try:
            re.compile(regex_test)
        except sre_compile.error as err:
            if not failures:
                failure_str += "##### Regular expression failures #####\n"
            failure_str += "%s --> %s\n" % (regex_test, str(err))
            failures.append(indx)
    if failures:
        failure_str += "#######################################\n\n"
        failures = sorted(failures, reverse=True)
        for indx in failures:
            del patterns[indx]
    _stderr(failure_str, quiet)
    return patterns

# #################################################### VARIABLES ##################################################### #

contributors = [Contributor("Stephen", "Bond", commits=1055, github="https://github.com/biologyguy"),
                Contributor("Karl", "Keat", commits=392, github="https://github.com/KarlKeat"),
                Contributor("Jeremy", "Labarge", commits=26, github="https://github.com/biojerm"),
                Contributor("Dustin", "Mitchell", commits=12, github="https://github.com/djmitche"),
                Contributor("Jason", "Bowen", commits=6, github="https://github.com/jwbowen"),
                Contributor("Todd", "Smith", commits=5, github="https://github.com/etiology"),
                Contributor("Sofia", "Barreira", commits=2, github="https://github.com/alicarea"),
                Contributor("Alexander", "Jones", commits=2, github="https://github.com/alexanjm"),
                Contributor("Adam", "Palmer", commits=2, github="https://github.com/apalm112"),
                Contributor("Helena", "Mendes-Soares", commits=1, github="https://github.com/mendessoares")]

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
                            "help": "Concatenate multiple sequences into a single solid string. Pass in "
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
                                       "specify number of columns for deleted IDs"},
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
            "extract_feature_sequences": {"flag": "efs",
                                          "action": "append",
                                          "nargs": "+",
                                          "metavar": "feature",
                                          "help": "Pull specific regions out of sequences based on the presence of"
                                          " specific features in the annotations"},
            "extract_regions": {"flag": "er",
                                "action": "append",
                                "nargs": "+",
                                "metavar": "positions",
                                "help": "Pull out specific residues"},
            "find_CpG": {"flag": "fcpg",
                         "action": "store_true",
                         "help": "Predict regions under strong purifying selection based on high CpG content"},
            "find_orfs": {"flag": "orf",
                          "action": "store_true",
                          "help": "Finds all the open reading frames in the sequences and their reverse complements."},
            "find_pattern": {"flag": "fp",
                             "action": "store",
                             "nargs": "+",
                             "metavar": "<regex>",
                             "help": "Search for subsequences, returning the start positions of all matches. Include "
                                     "the word 'ambig' to search with ambiguous character codes."},
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
            "prosite_scan": {"flag": "psc",
                             "action": "append",
                             "nargs": "?",
                             "help": "Annotate sequence features with Prosite Scan. Include the word 'strict' to "
                                     "exclude common matches. Internet connection required."},
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
            "pull_records_with_feature": {"flag": "prf",
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
            "transmembrane_domains": {"flag": "tmd",
                                      "nargs": "*",
                                      "action": "append",
                                      "metavar": "Job ID",
                                      "help": "Annotate transmembrane domains using TOPCONS2 (internet required)"},
            "uppercase": {"flag": "uc",
                          "action": "store_true",
                          "help": "Convert all sequences to uppercase"}}

sb_modifiers = {"alpha": {"flag": "a",
                          "action": "store",
                          "help": "If you want the file read with a specific alphabet"},
                "in_format": {"flag": "f",
                              "action": "store",
                              "help": "If SeqBuddy can't guess the file format, try specifying it directly"},
                "in_place": {"flag": "i",
                             "action": "store_true",
                             "help": "Rewrite the input file in-place. Be careful!"},
                "keep_temp": {"flag": "k",
                              "action": "store",
                              "help": "Save temporary files created by generate_tree in current working directory"},
                "out_format": {"flag": "o",
                               "metavar": "",
                               "action": "store",
                               "help": "If you want a specific format output"},
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
             "bootstrap": {"flag": "bts",
                           "action": "append",
                           "nargs": "?",
                           "type": int,
                           "metavar": "# of bootstraps (default=1)",
                           "help": "Generate bootstrap alignment(s) from the input alignment(s)."},
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
                                "action": "store",
                                "metavar": "args",
                                "help": "Remove alignment rows with IDs that contain matches to the provided patterns"},
             "enforce_triplets": {"flag": "et",
                                  "action": "store_true",
                                  "help": "Shift gaps so sequences are organized in triplets"},
             "extract_feature_sequences": {"flag": "efs",
                                           "action": "append",
                                           "nargs": "+",
                                           "metavar": "feature",
                                           "help": "Pull specific regions out of alignments based on the presence of"
                                           " specific features in the annotations"},
             "extract_regions": {"flag": "er",
                                 "action": "append",
                                 "nargs": "+",
                                 "metavar": "positions",
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
                               "help": "If AlignBuddy can't guess the file format, try specifying it directly"},
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

pb_flags = {"collapse_polytomies": {"flag": "cpt",
                                    "action": "append",
                                    "nargs": "*",
                                    "metavar": ("threshold", "{'bootstrap', 'length'}"),
                                    "help": "Create a polytomy from any nodes with less support that threshold"},
            "consensus_tree": {"flag": "ct",
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
