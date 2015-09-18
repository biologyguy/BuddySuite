#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# NOTE: BioPython 16.6+ required.

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