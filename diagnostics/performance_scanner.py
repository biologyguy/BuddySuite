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
        "lowercase": "",
        "ave_seq_length": "",
        "back_translate": "",
        "bl2seq": "",
        "blast": "mnemiopsis_leidyi",
        "clean_seq": "",
        "complement": "",
        "concat_seqs": "",
        "count_codons": "",
    }

    if in_args.command == "all":
        for flag, args in opts_sb.items():
            sys.stdout.write("%s: " % flag)
            sys.stdout.flush()
            if flag == "bl2seq":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy seqs_25.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=10)
                sys.stdout.write("%s\n" % round(timer / 10, 3))
            elif flag == "concat_seqs":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=10)
                sys.stdout.write("%s\n" % round(timer / 10, 3))
            else:
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=10)
                sys.stdout.write("%s\n" % round(timer / 10, 3))
    elif in_args.command != "all" and in_args.command in opts_sb:
        sys.stdout.write("%s: " % in_args.command)
        sys.stdout.flush()
        if in_args.command == "bl2seq":
            timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy seqs_25.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=10)
            sys.stdout.write("%s\n" % round(timer / 10, 3))
        elif in_args.command == "concat_seqs":
            timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=10)
            sys.stdout.write("%s\n" % round(timer / 10, 3))
        else:
            timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=10)
            sys.stdout.write("%s\n" % round(timer / 10, 3))
    else:
        sys.stdout.write("Invalid command\n")
