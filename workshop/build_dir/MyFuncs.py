#!/usr/bin/env python3

from multiprocessing import Process, cpu_count
from sys import stdout, exit, stderr
from time import time
from math import floor
import os
from tempfile import TemporaryDirectory
from shutil import copytree, rmtree
from re import sub


# maybe use curses library in the future to extend this for multi-line printing
class DynamicPrint():
    def __init__(self, out_type=stdout):
        self._last_print = ""
        self._next_print = ""
        self._writer = self._write()
        self.out_type = out_type

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
        content = sub("\t", "    ", content)
        self._next_print = content
        next(self._writer)


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


def run_multicore_function(iterable, function, func_args=False, max_processes=0, quiet=False, out_type=stdout):
        # fun little piece of abstraction here... directly pass in a function that is going to be looped over, and
        # fork those loops onto independent processes. Any arguments the function needs must be provided as a list.
        d_print = DynamicPrint(out_type)
        cpus = cpu_count()
        if max_processes == 0:
            if cpus > 7:
                max_processes = cpus - 3
            elif cpus > 3:
                max_processes = cpus - 2
            elif cpus > 1:
                max_processes = cpus - 1
            else:
                max_processes = 1

        if max_processes > cpus:
            max_processes = cpus

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


class TempDir():
    def __init__(self):
        self.dir = next(self._make_dir())
        self.path = self.dir.name

    def _make_dir(self):
        tmp_dir = TemporaryDirectory()
        yield tmp_dir
        rmtree(self.path)

    def save(self, location):
        if os.path.isdir(location):
            print("Save Error: Indicated output folder already exists in TempDir.save(%s)" % location, file=stderr)
            return False
        else:
            copytree(self.dir.name, location)
            return True


class TempFile():
    # I really don't like the behavior of tempfile.[Named]TemporaryFile(), so hack TemporaryDirectory() via TempDir()
    def __init__(self):
        self._tmp_dir = TempDir()  # This needs to be a persistent (ie self.) variable, or the directory will be deleted
        dir_hash = self._tmp_dir.path.split("/")[-1]
        self.path = "%s/%s" % (self._tmp_dir.path, dir_hash)
        self.handle = None

    def open(self, mode="w"):
        if not self.handle:
            self.handle = open(self.path, mode)

    def close(self):
        if self.handle:
            self.handle.close()
            self.handle = None

    def write(self, content, mode="a"):
        if mode not in ["w", "a"]:
            print("Write Error: mode must be 'w' or 'a' in TempFile.write()", file=stderr)
            return False
        already_open = True if self.handle else False
        if not already_open:
            self.open(mode)
        if mode == "a":
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
        with open(self.path, "r") as ifile:
            content = ifile.read()
        if already_open:
            self.open(mode="a")
            self.handle.seek(position)
        return content

    def save(self, location):
        with open(location, "w") as ofile:
            ofile.write(self.read())
        return


class SafetyValve():  # Use this class if you're afraid of an infinit loop
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
            exit("Error: You just popped your global_reps safety valve. %s" % message)
    
    def test(self, state, message=""):  # test() keeps track of some variable 'state' to see if its value keeps changing
        if self.state == str(state):
            self.state_reps -= 1
        else:
            self.state_reps = self._start_state_reps
            self.state = str(state)
            
        if self.state_reps == 0:
            exit("Error: You just popped your state_reps safety valve. %s" % message)


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
