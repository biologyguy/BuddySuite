#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests 3rd party alignment software that is wrapped by AlignBuddy
Note that the Linux and Mac versions of ClustalW and PAGAN return different results, so
it is necessary to check for two different hashes (maybe three, if Windows also acts differently)
"""

import pytest
import os
from unittest import mock
from shutil import which
from subprocess import Popen, PIPE
from ... import AlignBuddy as Alb
from ... import SeqBuddy as Sb
from ... import buddy_resources as br


# ###########################################  'ga', '--generate_alignment' ########################################## #

# ##########   PAGAN   ########## #
# This is tested for PAGAN version 0.61. NOTE: Do not split these up. Only one instance of Pagan can run at a time
def test_pagan(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan')
    assert alb_helpers.align2hash(tester) in ['da1c6bb365e2da8cb4e7fad32d7dafdb',
                                              '1219647676b359a5ad0be6d9dda81c73']
    # NEXUS
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f nexus')
    assert alb_helpers.align2hash(tester) in ['f93607e234441a2577fa7d8a387ef7ec',
                                              '42bfddd38fa4ed75a99841abf2112e54']
    # PHYLIPI
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f phylipi')
    assert alb_helpers.align2hash(tester) in ['09dd492fde598670d7cfee61d4e2eab8',
                                              '438e1551b3f1c8526fc8a44eaf2a3dc1']
    # PHYLIPS
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f phylips')
    assert alb_helpers.align2hash(tester) in ['249c88cb64d41c47388514c65bf8fff1',
                                              '6366e50da5a6b33d2d281d6ea13df0b7']
    # Multi-param
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f nexus --translate')
    assert alb_helpers.align2hash(tester) == 'dd140ec4eb895ce75d574498a58aa28a'

    # A few edge cases
    tester = sb_resources.get_one("d f")
    tester = Sb.pull_recs(tester, "α[2345]")
    Alb.generate_msa(tester, "pagan", "-f foo", quiet=True)

    tester = sb_resources.get_one("d f")
    tester = Sb.pull_recs(tester, "α[2345]")
    Alb.generate_msa(tester, "pagan", "-f nexus", quiet=True)

    tester = sb_resources.get_one("d f")
    tester = Sb.pull_recs(tester, "α[2345]")
    Alb.generate_msa(tester, "pagan", "-f phylipi", quiet=True)


# ##########   PRANK   ########## #
# PRANK is not deterministic, so just test that something reasonable is returned
def test_prank_inputs(sb_resources):
    # FASTA
    tester = Sb.pull_recs(sb_resources.get_one("d f"), 'α1')
    tester = Alb.generate_msa(tester, 'prank', '-once')
    assert tester.out_format == 'fasta'


def test_prank_outputs1(sb_resources):
    # NEXUS
    tester = Sb.pull_recs(sb_resources.get_one("d f"), 'α1')
    tester = Alb.generate_msa(tester, 'prank', '-f=nexus -once')
    assert tester.out_format == 'nexus'


def test_prank_outputs2(sb_resources):
    # PHYLIPI
    tester = Sb.pull_recs(sb_resources.get_one("d f"), 'α1')
    tester = Alb.generate_msa(tester, 'prank', params='-f=phylipi -once')
    assert tester.out_format == 'phylip-relaxed'


def test_prank_outputs3(sb_resources):
    # PHYLIPS
    tester = Sb.pull_recs(sb_resources.get_one("d f"), 'α1')
    tester = Alb.generate_msa(tester, 'prank', params='-f=phylips -once')
    assert tester.out_format == 'phylipsr'


# ##########   MUSCLE   ########## #
def test_muscle_inputs(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'muscle')
    assert alb_helpers.align2hash(tester) == '5ec18f3e0c9f5cf96944a1abb130232f'


def test_muscle_outputs(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'muscle', '-clw')
    assert alb_helpers.align2hash(tester) == '91542667cef761ccaf39d8cb4e877944'


def test_muscle_multi_param(sb_resources, alb_helpers):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'muscle', '-clw -diags')
    assert alb_helpers.align2hash(tester) == '91542667cef761ccaf39d8cb4e877944'

# ##########   CLUSTAL   ########## #
clustalw_bin = 'clustalw' if which('clustalw') else 'clustalw2'


def test_clustalw_inputs(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalw_bin)
    assert alb_helpers.align2hash(tester) in ['955440b5139c8e6d7d3843b7acab8446',
                                              'efc9b04f73c72036aa230a8d72da228b']


def test_clustalw_outputs1(sb_resources, alb_helpers):
    # NEXUS
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalw_bin, '-output=nexus')
    assert alb_helpers.align2hash(tester) in ['f4a61a8c2d08a1d84a736231a4035e2e',
                                              '7ca92360f0787164664843c895dd98f2']


def test_clustalw_outputs2(sb_resources, alb_helpers):
    # PHYLIP
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalw_bin, '-output=phylip')
    assert alb_helpers.align2hash(tester) in ['a9490f124039c6a2a6193d27d3d01205',
                                              'd3cc272a45fbde4b759460faa8e63ebc']


def test_clustalw_outputs3(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalw_bin, '-output=fasta')
    assert alb_helpers.align2hash(tester) in ['955440b5139c8e6d7d3843b7acab8446',
                                              'efc9b04f73c72036aa230a8d72da228b']


def test_clustalw_multi_param(sb_resources, alb_helpers):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalw_bin, '-output=phylip -noweights')
    assert alb_helpers.align2hash(tester) == 'ae9126eb8c482a82d4060d175803c478'


# ##########   CLUSTAL Omega   ########## #
clustalo_bin = 'clustalo' if which('clustalo') else 'clustalomega'
clustalo_version = Popen("{0} --version".format(clustalo_bin), shell=True, stdout=PIPE).communicate()[0].decode()
if clustalo_version not in ["1.2.2", "1.0.3"]:
    raise ValueError("Untested CLustalO version (%s). Please update the tests as necessary." % clustalo_version)

def test_clustalomega_inputs1(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalo_bin)
    assert alb_helpers.align2hash(tester) in ['f5afdc7c76ab822bdc95230329766aba', 'e2f0efc90372b23a1c753629d43a79c4']


def test_clustalomega_inputs2(sb_resources, alb_helpers):
    # PHYLIP
    tester = sb_resources.get_one("d py")
    tester = Alb.generate_msa(tester, clustalo_bin)
    assert alb_helpers.align2hash(tester) in ['8299780bf9485b89a2f3462ead666142', '5d808493da6c0f43a572b2a9257dce4f']


def test_clustalomega_inputs3(sb_resources, alb_helpers):
    # STOCKHOLM
    tester = sb_resources.get_one("d s")
    tester = Alb.generate_msa(tester, clustalo_bin)
    assert alb_helpers.align2hash(tester) in ['d6654e3db3818cc3427cb9241113fdfa', 'aeb2c5926843402cf620299802946224']


def test_clustalomega_outputs1(sb_resources, alb_helpers):
    # CLUSTAL
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalo_bin, '--outfmt=clustal')
    assert alb_helpers.align2hash(tester) in ['970f6e4389f77a30563763937a3d32bc', '50c48b43089f528714c7933cd0f3f91c']


def test_clustalomega_outputs2(sb_resources, alb_helpers):
    # PHYLIP
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalo_bin, '--outfmt=phylip')
    assert alb_helpers.align2hash(tester) in ['692c6af848bd90966f15908903894dbd', '75ec68313368dac249b40fe63b26777e']


def test_clustalomega_outputs3(sb_resources, alb_helpers):
    # STOCKHOLM
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalo_bin, '--outfmt=stockholm')
    assert alb_helpers.align2hash(tester) in ['4c24975c033abcf15911a61cb9663a97', '75ec68313368dac249b40fe63b26777e']


def test_clustalomega_multi_param(sb_resources, alb_helpers):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalo_bin, '--outfmt=clustal --iter=1')
    assert alb_helpers.align2hash(tester) in ['25480f7a9340ff643bb7eeb326e8f981', '9c55bb8be0cc89c5346d8f699e97cc87']


# ##########   MAFFT   ########## #
def test_mafft_inputs(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'mafft')
    assert alb_helpers.align2hash(tester) == 'f94e0fd591dad83bd94201f0af038904'


def test_mafft_outputs(sb_resources, alb_helpers):
    # CLUSTAL
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'mafft', '--clustalout')
    assert alb_helpers.align2hash(tester) == 'd6046c77e2bdb5683188e5de653affe5'


def test_mafft_multi_param(sb_resources, alb_helpers):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'mafft', '--clustalout --noscore')
    assert alb_helpers.align2hash(tester) == 'd6046c77e2bdb5683188e5de653affe5'


def test_generate_alignment_keep_temp(monkeypatch, sb_resources):
    tester = sb_resources.get_one("d f")
    temp_dir = br.TempDir()
    temp_dir.subdir("ga_temp_files")

    def ask_false(*ask_args):
        if ask_args:
            pass
        return False

    def ask_true(*ask_args):
        if ask_args:
            pass
        return True

    try:
        monkeypatch.setattr("buddysuite.buddy_resources.ask", ask_false)
    except ImportError:
        monkeypatch.setattr("buddy_resources.ask", ask_false)
    with pytest.raises(SystemExit):
        Alb.generate_msa(tester, clustalo_bin, keep_temp="%s/ga_temp_files" % temp_dir.path)

    try:
        monkeypatch.setattr("buddysuite.buddy_resources.ask", ask_true)
    except ImportError:
        monkeypatch.setattr("buddy_resources.ask", ask_true)
    Alb.generate_msa(tester, clustalo_bin, keep_temp="%s/ga_temp_files" % temp_dir.path)
    assert os.path.isfile("%s/ga_temp_files/result" % temp_dir.path)
    assert os.path.isfile("%s/ga_temp_files/tmp.fa" % temp_dir.path)


def test_generate_alignments_genbank(sb_resources, alb_helpers):
    tester = sb_resources.get_one("p g")
    tester = Alb.generate_msa(tester, "mafft")
    assert alb_helpers.align2hash(tester) == "f894ff6060ec5c2904f48ba0c5cdc8fd"


def test_generate_alignments_edges1(capsys, sb_resources):
    tester = sb_resources.get_one("d f")

    with pytest.raises(AttributeError) as e:
        Alb.generate_msa(tester, "foo")
    assert "foo is not a supported alignment tool." in str(e)

    # noinspection PyUnresolvedReferences
    with mock.patch.dict('os.environ'):
        del os.environ['PATH']
        with pytest.raises(SystemError) as err:
            Alb.generate_msa(tester, "mafft")
        assert "#### Could not find mafft in $PATH. ####" in str(err)


args = [("prank", "-f=phylipi"), (clustalo_bin, "--outfmt=foo"), (clustalw_bin, "-output=foo"),
        ("prank", "-f=nexus"), ("prank", "-f=foo")]


@pytest.mark.parametrize("tool,params", args)
def test_generate_alignments_edges2(tool, params, sb_resources):
    tester = sb_resources.get_one("d f")
    tester = Sb.pull_recs(tester, "α[2345]")
    Alb.generate_msa(tester, tool, params, quiet=True)
