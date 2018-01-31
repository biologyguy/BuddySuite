#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests 3rd party software that is wrapped by SeqBuddy
"""

import pytest
from unittest import mock
from subprocess import Popen, PIPE
import re
import SeqBuddy as Sb
import buddy_resources as br

blast_version = Popen("blastn -version", shell=True, stdout=PIPE).communicate()[0].decode()
blast_version = re.search("[0-9]+\.[0-9]+\.[0-9]+", blast_version).group(0)
if blast_version not in ["2.2.28", "2.2.29", "2.2.30", "2.2.31", "2.3.0", "2.4.0", "2.5.0", "2.6.0", "2.7.1"]:
    raise ValueError("Untested Blast version (%s). Please update the tests as necessary "
                     "(each version of blast seems to do something a little different...)" % blast_version)


# ######################  '-bl2s', '--bl2seq' ###################### #
def test_bl2seq(sb_resources, hf):
    result = Sb.bl2seq(sb_resources.get_one("d f"))
    assert hf.string2hash(str(result)) in ['87dbd3baeb59285ad25e6473c87bb5bb', '8280eb4010208db891020a96ad783edb']

    result = Sb.bl2seq(sb_resources.get_one("p f"))
    assert hf.string2hash(str(result)) in ['248d4c53d7947c4c8dfd7c415bfbfbf2', '33b393de45d0d628a217bf9107ec9719',
                                           'ca7105bf6646c1ab3f07efeea57a69df']


# ######################  '-bl', '--blast' ###################### #
def test_blastn(sb_resources, sb_odd_resources, hf, monkeypatch):
    tester = Sb.pull_recs(sb_resources.get_one("d f"), '8', True)
    tester = Sb.blast(tester, sb_odd_resources["blastn"])
    assert hf.buddy2hash(tester) == "95c417b6c2846d1b7a1a07f50c62ff8a"

    tester = Sb.SeqBuddy(">Seq1\nATGCGCGCTACGCTAGCTAGCTAGCTCGCATGCAT")
    tester = Sb.blast(tester, sb_odd_resources["blastn"])
    assert len(tester.records) == 0

    with pytest.raises(RuntimeError) as e:
        tester = sb_resources.get_one("d f")
        Sb.blast(tester, "Mnemiopsis_cds.nhr")
    assert "The .nhr file of your BLASTN database was not found" in str(e.value)

    monkeypatch.setattr(Sb, "_check_for_blast_bin", lambda *_: False)
    with pytest.raises(SystemError) as e:
        Sb.blast(tester, sb_odd_resources["blastn"])
    assert 'blastn not found in system path' in str(e.value)


def test_blastp(sb_resources, sb_odd_resources, hf, monkeypatch):
    seqbuddy = Sb.pull_recs(sb_resources.get_one('p f'), '8', True)
    tester = Sb.blast(seqbuddy, sb_odd_resources["blastp"])
    assert hf.buddy2hash(tester) in ["4237c79672c1cf1d4a9bdb160a53a4b9", "118d4f412e2a362b9d16130abbf395c5"]

    with pytest.raises(RuntimeError) as e:
        tester = sb_resources.get_one("p f")
        Sb.blast(tester, "Mnemiopsis_pep.phr")
    assert "The .phr file of your BLASTP database was not found" in str(e.value)

    monkeypatch.setattr(Sb, "_check_for_blast_bin", lambda *_: False)
    with pytest.raises(SystemError) as e:
        Sb.blast(tester, sb_odd_resources["blastp"])
    assert 'blastp not found in system path' in str(e.value)


def test_makeblastdb(monkeypatch, sb_resources, hf):
    def mock_check_blast_bin(binary):
        if binary == "makeblastdb":
            return False
        else:
            return True

    subject = Sb.pull_recs(sb_resources.get_one('p f'), '8', True)
    query = Sb.pull_recs(sb_resources.get_one('p f'), 'α[^8]', True)
    output = Sb.blast(subject, query)
    assert hf.buddy2hash(output) in ["4639da7978256eb8dae0e9e7a1ad3d01", "638938107e33174206e9fce9b789fe64"]

    monkeypatch.setattr(Sb, "_check_for_blast_bin", mock_check_blast_bin)
    with pytest.raises(SystemError) as err:
        Sb.blast(subject, query)
    assert "makeblastdb not found in system path." in str(err)


# #####################  '-psc', '--prosite_scan' ###################### ##
"""  Need to figure out a way of applying a timeout to these...
def test_prosite_scan(sb_resources, hf):
    seqbuddy = sb_resources.get_one("d f")
    ps_scan = Sb.PrositeScan(seqbuddy)
    ps_scan.run()
    assert hf.buddy2hash(ps_scan.seqbuddy) == "e9090efdd362d527a115049dfced42cd"


# ######################  '-tmd', '--transmembrane_domains' ###################### #
def test_transmembrane_domains_pep(sb_resources, hf):
    tester = sb_resources.get_one("p f")
    Sb.pull_recs(tester, "Panxα[234]")
    tester = Sb.transmembrane_domains(tester, quiet=True)
    tester.out_format = "gb"
    assert hf.buddy2hash(tester) == "3d62e2548b9181574d51907d2205b36c"


def test_transmembrane_domains_cds(sb_resources, hf):
    tmp_dir = br.TempDir()
    tmp_dir.subdir("topcons")
    tester = sb_resources.get_one("d f")
    Sb.pull_recs(tester, "Panxα[234]")
    tester = Sb.transmembrane_domains(tester, quiet=True, keep_temp="%s%stopcons" % (tmp_dir.path, os.sep))
    tester.out_format = "gb"
    assert hf.buddy2hash(tester) == "479eb1c8728c959b813c97962cac545a"
    _root, dirs, files = next(br.walklevel("%s%stopcons" % (tmp_dir.path, os.sep)))
    _root, dirs, files = next(br.walklevel("{0}{1}topcons{1}{2}".format(tmp_dir.path, os.sep, dirs[0])))
    assert files
"""