#!/usr/bin/env python3

import sys
import argparse
import timeit

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
        "num_seqs": "",
        "lowercase": ""
    }

    for flag, args in opts_sb.items():
        sys.stdout.write("%s: " % flag)
        sys.stdout.flush()
        timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=10)
        sys.stdout.write("%s\n" % round(timer / 10, 3))