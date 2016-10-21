#!/usr/bin/env python3

import sys
import argparse
import timeit
from buddysuite import SeqBuddy
#from subprocess import Popen
from buddysuite import AlignBuddy
from buddysuite import PhyloBuddy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="performanceScanner", description="Check function time", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("regex", help="pattern to search", action="store")
    # parser.add_argument("replace", help="replacement", action="store")
    # parser.add_argument("file", help="input file", action="store")
    # parser.add_argument("-m", "--module", help="specify a single module", action='store')
    parser.add_argument("-c", "--command", help="specify a single module", action='store')
    # use seq buddy for now, then ad the --module option
    in_args = parser.parse_args()

    # create a dictionary with the default parameters for each function (-c, --command)
    opts_sb = {
        'list_ids': '3',
    }

    # if in_args.module == "all":
    if in_args.command == "list_ids":
        timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa -li opts_sb[in_args.command] -q", shell=True, stderr=PIPE, stdout=PIPE).wait()', number=10)
        newtimer = timer / 10
        print(newtimer)
    elif in_args.command == "num_seqs":
        timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa -ns -q", shell=True, stderr=PIPE, stdout=PIPE).wait()', number=10)
        newtimer = timer / 10
        print(newtimer)
    elif in_args.command == "lowercase":
        timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa -lc", shell=True, stderr=PIPE, stdout=PIPE).wait()', number=10)
        newtimer = timer / 10
        print(newtimer)
