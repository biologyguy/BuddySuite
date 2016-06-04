#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: MyFuncs.py
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Collection of useful classes and functions
"""

from multiprocessing import Process, cpu_count
import sys
from time import time
from math import floor, ceil
import os
from tempfile import TemporaryDirectory
from shutil import copytree, rmtree, copyfile
import re
import string
from random import choice


class Timer(object):
    def __init__(self):
        self.current_time = round(time())

    def start(self):
        self.current_time = round(time())
        return

    def end(self):
        return pretty_time(round(time()) - self.current_time)


class RunTime(object):
    def __init__(self, prefix="", postfix="", out_type=sys.stdout):
        self.out_type = out_type
        self.prefix = prefix
        self.postfix = postfix
        self.running_process = None

    def _run(self, check_file_path):
        d_print = DynamicPrint(self.out_type)
        start_time = round(time())
        elapsed = 0
        while True:
            with open("%s" % check_file_path, "r") as ifile:
                if ifile.read() == "Running":
                    d_print.write("%s%s%s" % (self.prefix, pretty_time(elapsed), self.postfix))
                    elapsed = round(time()) - start_time
                else:
                    d_print.write("%s%s%s\n" % (self.prefix, pretty_time(elapsed), self.postfix))
                    break
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
        self._last_print = ""
        self._next_print = ""
        self._writer = self._write()

        out_type = sys.stdout if out_type == "stdout" else out_type
        out_type = sys.stderr if out_type == "stderr" else out_type
        self.out_type = out_type
        self.quiet = quiet

    def _write(self):
        try:
            while True:
                self.out_type.write("\r%s\r%s" % (" " * len(self._last_print), self._next_print),)
                self.out_type.flush()
                self._last_print = self._next_print
                yield
        finally:
            self.out_type.write("")

    def write(self, content):
        if not self.quiet:
            content = re.sub("\t", "    ", content)
            self._next_print = content
            next(self._writer)
        return

    def new_line(self):
        if not self.quiet:
            self.out_type.write("\n")
            self.out_type.flush()
            self._last_print = ""
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
        return '%s %s' % (num, ['', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'][magnitude])
    elif mode == 'medium':
        return '%s %s' % (num, ['', 'Kilo', 'Mega', 'Giga', 'Tera', 'Peta', 'Exa', 'Zetta', 'Yotta'][magnitude])
    elif mode == 'long':
        return '%s %s' % (num, ['', 'Thousand', 'Million', 'Billion', 'Trillion', 'Quadrillion', 'Quintillion',
                                'Sextillion', 'Septillion'][magnitude])
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
            while 1:     # Only fork a new process when there is a free processor.
                if running_processes < max_processes:
                    # Start new process
                    if not quiet:
                        d_print.write("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

                    if func_args:
                        if not isinstance(func_args, list):
                            exit("Error in run_multicore_function(): The arguments passed into the multi-thread "
                                 "function must be provided as a list")
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

    def _make_dir(self):
        tmp_dir = TemporaryDirectory()
        yield tmp_dir
        rmtree(self.path)

    def subdir(self, dir_name=None):
        if not dir_name:
            dir_name = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])
            while dir_name in self.subdirs:  # Catch the very unlikely case that a duplicate occurs
                dir_name = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])

        subdir_path = "%s/%s" % (self.path, dir_name)
        os.mkdir(subdir_path)
        self.subdirs.append(dir_name)
        return subdir_path

    def del_subdir(self, _dir):
        _dir = _dir.split("/")[-1]
        del self.subdirs[self.subdirs.index(_dir)]
        rmtree("%s/%s" % (self.path, _dir))
        return

    def save(self, location, keep_hash=False):
        location = location if not keep_hash else "%s/%s" % (location, self.path.split("/")[-1])
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
        dir_hash = self._tmp_dir.path.split("/")[-1]
        self.name = dir_hash
        self.path = "%s/%s" % (self._tmp_dir.path, dir_hash)
        open(self.path, "w").close()
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
        self.write("", mode="w")
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
    for root, dirs, files in os.walk(source):
        if not os.path.isdir(dest):
            os.makedirs(dest)
        for each_file in files:
            rel_path = re.sub(source, '', root)
            rel_path = rel_path.lstrip(os.sep)
            dest_path = os.path.join(dest, rel_path, each_file)
            copyfile(os.path.join(root, each_file), dest_path)


def normalize(data, trim_ends=1.0):
    if 0. > trim_ends > 1.0:
        raise ValueError("normalize() trim_ends parameter should be between 0.5 and 1.0")

    if trim_ends > 0.5:
        trim_ends = 1 - trim_ends

    max_limit = ceil(len(data) * (1 - trim_ends)) - 1
    min_limit = -1 * (max_limit + 1)

    if type(data) == dict:
        sorted_data = sorted([data[key] for key in data])
        _max = sorted_data[max_limit]
        _min = sorted_data[min_limit]
        data_range = _max - _min
        for key in data:
            data[key] = (data[key] - _min) / data_range
            data[key] = 1. if data[key] > 1. else data[key]
            data[key] = 0. if data[key] < 0. else data[key]

    else:
        sorted_data = sorted(data)
        _max = sorted_data[max_limit]
        _min = sorted_data[min_limit]
        data_range = _max - _min
        for i in range(len(data)):
            data[i] = (data[i] - _min) / data_range
            data[i] = 1. if data[i] > 1. else data[i]
            data[i] = 0. if data[i] < 0. else data[i]

    return data


# This will only work if SMTP is running locally
def sendmail(sender, recipient, subject, message):
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart

    msg = MIMEMultipart()
    msg.preamble = subject
    msg.add_header("From", sender)
    msg.add_header("Subject", subject)
    msg.add_header("To", recipient)

    msg.attach(MIMEText(message))

    smtp = smtplib.SMTP('localhost')
    smtp.starttls()
    smtp.sendmail(sender, recipient, msg.as_string())
    smtp.quit()
    return


def ask(input_prompt, default="yes"):
    if default == "yes":
        yes_list = ["yes", "y", '']
        no_list = ["no", "n", "abort"]
    else:
        yes_list = ["yes", "y"]
        no_list = ["no", "n", "abort", '']

    _response = input(input_prompt)
    while True:
        if _response.lower() in yes_list:
            return True
        elif _response.lower() in no_list:
            return False
        else:
            print("Response not understood. Valid options are 'yes' and 'no'.")
            _response = input(input_prompt)
