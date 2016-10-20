#!/usr/bin/env python3

import argparse
from buddysuite import SeqBuddy
from buddysuite import AlignBuddy
from buddysuite import PhyloBuddy


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="performanceScanner", description="Check function time",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    #parser.add_argument("regex", help="pattern to search", action="store")
    #parser.add_argument("replace", help="replacement", action="store")
    #parser.add_argument("file", help="input file", action="store")
    parser.add_argument("-m", "--module", help="specify a single module", action='store')
    

    in_args = parser.parse_args()

    print("Hello World\n")
