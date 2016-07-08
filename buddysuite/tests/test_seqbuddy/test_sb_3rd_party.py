#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests 3rd party software that is wrapped by SeqBuddy
"""

import pytest
import os
from unittest import mock
from subprocess import Popen, PIPE
import re

from ... import SeqBuddy as Sb
from ... import buddy_resources as br

blast_version = Popen("blastn -version", shell=True, stdout=PIPE).communicate()[0].decode()
blast_version = re.search("[0-9]+\.[0-9]+\.[0-9]+", blast_version).group(0)
if blast_version not in ["2.2.28", "2.2.29", "2.2.30", "2.2.31", "2.3.0", "2.4.0"]:
    raise ValueError("Untested Blast version (%s). Please update the tests as necessary "
                     "(each version of blast seems to do something a little different...)" % blast_version)


# ######################  '-bl2s', '--bl2seq' ###################### #
def test_bl2seq(sb_resources, sb_helpers):
    result = Sb.bl2seq(sb_resources.get_one("d f"))
    assert sb_helpers.string2hash(str(result)) in ['87dbd3baeb59285ad25e6473c87bb5bb',
                                                   '8280eb4010208db891020a96ad783edb']

    result = Sb.bl2seq(sb_resources.get_one("p f"))
    assert sb_helpers.string2hash(str(result)) in ['248d4c53d7947c4c8dfd7c415bfbfbf2',
                                                   '33b393de45d0d628a217bf9107ec9719',
                                                   'ca7105bf6646c1ab3f07efeea57a69df']


# ######################  '-bl', '--blast' ###################### #
def test_blastn(sb_resources, sb_odd_resources, sb_helpers):
    tester = Sb.pull_recs(sb_resources.get_one("d f"), '8', True)
    tester = Sb.blast(tester, sb_odd_resources["blastn"])
    assert sb_helpers.seqs2hash(tester) == "95c417b6c2846d1b7a1a07f50c62ff8a"

    with pytest.raises(RuntimeError) as e:
        tester = sb_resources.get_one("d f")
        Sb.blast(tester, "Mnemiopsis_cds.nhr")
    assert "The .nhr file of your blast database was not found" in str(e.value)

    try:
        with mock.patch("buddysuite.SeqBuddy._check_for_blast_bin", return_value=False):
            with pytest.raises(SystemError) as e:
                Sb.blast(tester, sb_odd_resources["blastn"])
            assert 'blastn not found in system path' in str(e.value)
    except ImportError:
        with mock.patch("SeqBuddy._check_for_blast_bin", return_value=False):
            with pytest.raises(SystemError) as e:
                Sb.blast(tester, sb_odd_resources["blastn"])
            assert 'blastn not found in system path' in str(e.value)
    tester = Sb.SeqBuddy(">Seq1\nATGCGCGCTACGCTAGCTAGCTAGCTCGCATGCAT")
    tester = Sb.blast(tester, sb_odd_resources["blastn"])
    assert len(tester.records) == 0


def test_blastp(sb_resources, sb_odd_resources, sb_helpers):
    seqbuddy = Sb.pull_recs(sb_resources.get_one('p f'), '8', True)
    tester = Sb.blast(seqbuddy, sb_odd_resources["blastp"])
    assert sb_helpers.seqs2hash(tester) in ["4237c79672c1cf1d4a9bdb160a53a4b9",
                                            "118d4f412e2a362b9d16130abbf395c5"]

    with pytest.raises(RuntimeError) as e:
        tester = sb_resources.get_one("p f")
        Sb.blast(tester, "Mnemiopsis_pep.phr")
    assert "The .phr file of your blast database was not found" in str(e.value)

    try:
        with mock.patch("buddysuite.SeqBuddy._check_for_blast_bin", return_value=False):
            with pytest.raises(SystemError) as e:
                Sb.blast(tester, sb_odd_resources["blastp"])
            assert 'blastp not found in system path' in str(e.value)
    except ImportError:
        with mock.patch("SeqBuddy._check_for_blast_bin", return_value=False):
            with pytest.raises(SystemError) as e:
                Sb.blast(tester, sb_odd_resources["blastp"])
            assert 'blastp not found in system path' in str(e.value)
