#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: alignbuddy_tests.py
version: 1, alpha
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Collection of PyTest unit tests for the buddy_resources.py module
"""

import pytest
from subprocess import Popen, PIPE

try:
    import workshop.buddy_resources as br
except ImportError:
    import buddy_resources as br


def test_versions():
    sb_ver = Popen("sb -v", stdout=PIPE, shell=True).communicate()
    assert sb_ver[0].decode() == str(br.VERSIONS["SeqBuddy"])

    alb_ver = Popen("alb -v", stdout=PIPE, shell=True).communicate()
    assert alb_ver[0].decode() == str(br.VERSIONS["AlignBuddy"])

    pb_ver = Popen("pb -v", stdout=PIPE, shell=True).communicate()
    assert pb_ver[0].decode() == str(br.VERSIONS["PhyloBuddy"])

    db_ver = Popen("db -v", stdout=PIPE, shell=True).communicate()
    assert db_ver[0].decode() == str(br.VERSIONS["DatabaseBuddy"])


# ######################  '_format_to_extension' ###################### #
def test_format_to_extension():
    ext_dict = {'fasta': 'fa', 'fa': 'fa', 'genbank': 'gb', 'gb': 'gb', 'nexus': 'nex', 'nex': 'nex', 'phylip': 'phy',
                'phy': 'phy', 'phylip-relaxed': 'phyr', 'phyr': 'phyr', 'stockholm': 'stklm', 'stklm': 'stklm'}
    for i, j in ext_dict.items():
        assert j == br.format_to_extension(i)