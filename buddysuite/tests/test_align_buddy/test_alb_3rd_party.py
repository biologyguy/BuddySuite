#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" tests basic functionality of AlignBuddy class """
import pytest
import os
from unittest import mock
try:
    from buddysuite import AlignBuddy as Alb
    from buddysuite import SeqBuddy as Sb
    from buddysuite import MyFuncs
except ImportError:
    import AlignBuddy as Alb
    import SeqBuddy as Sb
    import MyFuncs


# ###########################################  'ga', '--generate_alignment' ########################################## #
# This is tested for PAGAN version 0.61. NOTE: Do not split these up. Only one instance of Pagan can run at a time

def test_pagan(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan')
    assert alb_helpers.align_to_hash(tester) == 'da1c6bb365e2da8cb4e7fad32d7dafdb'

    # NEXUS
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f nexus')
    assert alb_helpers.align_to_hash(tester) == 'f93607e234441a2577fa7d8a387ef7ec'

    # PHYLIPI
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f phylipi')
    assert alb_helpers.align_to_hash(tester) == '09dd492fde598670d7cfee61d4e2eab8'

    # PHYLIPS
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f phylips')
    assert alb_helpers.align_to_hash(tester) == '249c88cb64d41c47388514c65bf8fff1'

    # Multi-param
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f nexus --translate')
    assert alb_helpers.align_to_hash(tester) == 'dd140ec4eb895ce75d574498a58aa28a'

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


def test_muscle_inputs(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'muscle')
    assert alb_helpers.align_to_hash(tester) == '5ec18f3e0c9f5cf96944a1abb130232f'


def test_muscle_outputs(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'muscle', '-clw')
    assert alb_helpers.align_to_hash(tester) == '91542667cef761ccaf39d8cb4e877944'


def test_muscle_multi_param(sb_resources, alb_helpers):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'muscle', '-clw -diags')
    assert alb_helpers.align_to_hash(tester) == '91542667cef761ccaf39d8cb4e877944'


def test_clustalw2_inputs(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalw2')
    assert alb_helpers.align_to_hash(tester) == '955440b5139c8e6d7d3843b7acab8446'


def test_clustalw2_outputs1(sb_resources, alb_helpers):
    # NEXUS
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalw2', '-output=nexus')
    assert alb_helpers.align_to_hash(tester) == 'f4a61a8c2d08a1d84a736231a4035e2e'


def test_clustalw2_outputs2(sb_resources, alb_helpers):
    # PHYLIP
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalw2', '-output=phylip')
    assert alb_helpers.align_to_hash(tester) == 'a9490f124039c6a2a6193d27d3d01205'


def test_clustalw2_outputs3(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalw2', '-output=fasta')
    assert alb_helpers.align_to_hash(tester) == '955440b5139c8e6d7d3843b7acab8446'


def test_clustalw2_multi_param(sb_resources, alb_helpers):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalw2', '-output=phylip -noweights')
    assert alb_helpers.align_to_hash(tester) == 'ae9126eb8c482a82d4060d175803c478'


def test_clustalomega_inputs1(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalomega')
    assert alb_helpers.align_to_hash(tester) == 'f5afdc7c76ab822bdc95230329766aba'


def test_clustalomega_inputs2(sb_resources, alb_helpers):
    # PHYLIP
    tester = sb_resources.get_one("d py")
    tester = Alb.generate_msa(tester, 'clustalomega')
    assert alb_helpers.align_to_hash(tester) == '8299780bf9485b89a2f3462ead666142'


def test_clustalomega_inputs3(sb_resources, alb_helpers):
    # STOCKHOLM
    tester = sb_resources.get_one("d s")
    tester = Alb.generate_msa(tester, 'clustalomega')
    assert alb_helpers.align_to_hash(tester) == 'd6654e3db3818cc3427cb9241113fdfa'


def test_clustalomega_outputs1(sb_resources, alb_helpers):
    # CLUSTAL
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=clustal')
    assert alb_helpers.align_to_hash(tester) == '970f6e4389f77a30563763937a3d32bc'


def test_clustalomega_outputs2(sb_resources, alb_helpers):
    # PHYLIP
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=phylip')
    assert alb_helpers.align_to_hash(tester) == '692c6af848bd90966f15908903894dbd'


def test_clustalomega_outputs3(sb_resources, alb_helpers):
    # STOCKHOLM
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=stockholm')
    assert alb_helpers.align_to_hash(tester) == '4c24975c033abcf15911a61cb9663a97'


def test_clustalomega_multi_param(sb_resources, alb_helpers):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'clustalomega', '--outfmt=clustal --iter=1')
    assert alb_helpers.align_to_hash(tester) == '25480f7a9340ff643bb7eeb326e8f981'


def test_mafft_inputs(sb_resources, alb_helpers):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'mafft')
    assert alb_helpers.align_to_hash(tester) == 'f94e0fd591dad83bd94201f0af038904'


def test_mafft_outputs(sb_resources, alb_helpers):
    # CLUSTAL
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'mafft', '--clustalout')
    assert alb_helpers.align_to_hash(tester) == 'd6046c77e2bdb5683188e5de653affe5'


def test_mafft_multi_param(sb_resources, alb_helpers):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'mafft', '--clustalout --noscore')
    assert alb_helpers.align_to_hash(tester) == 'd6046c77e2bdb5683188e5de653affe5'


def test_generate_alignment_keep_temp(monkeypatch, sb_resources):
    tester = sb_resources.get_one("d f")
    temp_dir = MyFuncs.TempDir()
    temp_dir.subdir("ga_temp_files")

    def ask_false(*ask_args):
        if ask_args:
            pass
        return False

    def ask_true(*ask_args):
        if ask_args:
            pass
        return True

    monkeypatch.setattr("MyFuncs.ask", ask_false)
    monkeypatch.setattr("buddysuite.MyFuncs.ask", ask_false)
    with pytest.raises(SystemExit):
        Alb.generate_msa(tester, "clustalomega", keep_temp="%s/ga_temp_files" % temp_dir.path)

    monkeypatch.setattr("MyFuncs.ask", ask_true)
    monkeypatch.setattr("buddysuite.MyFuncs.ask", ask_true)
    Alb.generate_msa(tester, "clustalomega", keep_temp="%s/ga_temp_files" % temp_dir.path)
    assert os.path.isfile("%s/ga_temp_files/result" % temp_dir.path)
    assert os.path.isfile("%s/ga_temp_files/tmp.fa" % temp_dir.path)


def test_generate_alignments_genbank(sb_resources, alb_helpers):
    tester = sb_resources.get_one("p g")
    tester = Alb.generate_msa(tester, "mafft")
    assert alb_helpers.align_to_hash(tester) == "f894ff6060ec5c2904f48ba0c5cdc8fd"


def test_generate_alignments_edges1(capsys, sb_resources):
    tester = sb_resources.get_one("d f")

    with pytest.raises(AttributeError) as e:
        Alb.generate_msa(tester, "foo")
    assert "foo is not a supported alignment tool." in str(e)

    # noinspection PyUnresolvedReferences
    with mock.patch.dict('os.environ'):
        del os.environ['PATH']
        with pytest.raises(SystemExit):
            Alb.generate_msa(tester, "mafft")
        out, err = capsys.readouterr()
        assert "#### Could not find mafft in $PATH. ####\n" in err


args = [("prank", "-f=phylipi"), ("clustalomega", "--outfmt=foo"), ("clustalw2", "-output=foo"),
        ("prank", "-f=nexus"), ("prank", "-f=foo")]


@pytest.mark.parametrize("tool,params", args)
def test_generate_alignments_edges2(tool, params, sb_resources):
    tester = sb_resources.get_one("d f")
    tester = Sb.pull_recs(tester, "α[2345]")
    Alb.generate_msa(tester, tool, params, quiet=True)
