#!/usr/bin/env python3

import argparse
import timeit
from buddysuite import SeqBuddy
from buddysuite import AlignBuddy
from buddysuite import PhyloBuddy
from subprocess import Popen


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="performanceScanner", description="Check function time",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("regex", help="pattern to search", action="store")
    # parser.add_argument("replace", help="replacement", action="store")
    # parser.add_argument("file", help="input file", action="store")
    # parser.add_argument("-m", "--module", help="specify a single module", action='store')
    parser.add_argument("-f", "--file", help="specify the input file", action='store')
    parser.add_argument("-c", "--command", help="specify a single module", action='store')
    # use seq buddy for now, then ad the --module option
    in_args = parser.parse_args()

    # create wrapper to use functions with arguments inside timeit:
    def wrapper(func, *args):
        def wrapped():
            return func(*args)
        return wrapped

    # create a dictionary with the default parameters for each function (-c, --command)
    opts_sb = {
        'list_ids': '3',
        'lowercase': '',
    }

    # if in_args.module == "all":
    if in_args.command == "list_ids":
        short_list = [in_args.file, "--" + in_args.command, opts_sb[in_args.command]]
        wrapped = wrapper(SeqBuddy, short_list)
        timeit.timeit(wrapped, number=1000)
        # print(timeit.timeit("sb in_args.file ", number=1000))
        # Popen("echo Hello", shell=True).wait()
