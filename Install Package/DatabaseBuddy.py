#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Jan 15 2015 

"""
DESCRIPTION OF PROGRAM
"""

import sys
import os
import re
import shutil

import argparse


# ##################################################### WISH LIST #################################################### #
def get_genbank_file():
    x = 1
    return x


def run_prosite():
    x = 1
    return x


parser = argparse.ArgumentParser(prog="DatabaseBuddy", description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("positional_arg1", help="", action="store")
parser.add_argument("-t", "--true", help="", action="store_true", default=False)
parser.add_argument("-c", "--choice", help="", type=str, choices=["", ""], default=False)
parser.add_argument("-m", "--multi_arg", nargs="+", help="", default=[])

in_args = parser.parse_args()


class NewClass():
    """DESCRIPTION OF CLASS"""
    def __init__(self):
        self.x = 1

    def class_def(self):
        self.x = 1
        return self.x


def def1():
    """DESCRIPTION OF FUNC"""
    x = 1
    return x


def def2():
    """DESCRIPTION OF FUNC"""
    x = 1
    return x


if __name__ == '__main__':
    print('Hello')