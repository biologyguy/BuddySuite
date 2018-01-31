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
import re
import AlignBuddy as Alb
import SeqBuddy as Sb
import buddy_resources as br


class MockPopen(object):
    def __init__(self, *_, **__):
        self.foo = [_, __]

    @staticmethod
    def wait(*_, **__):
        return

    @staticmethod
    def communicate(*_, **__):
        return ["out".encode(), "err".encode()]


# ###########################################  'ga', '--generate_alignment' ########################################## #

# ##########   PAGAN   ########## #
# This is tested for PAGAN version 0.61. NOTE: Do not split these up. Only one instance of Pagan can run at a time
@br.skip_windows
def test_pagan(sb_resources, hf):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan')
    assert hf.buddy2hash(tester) in ['da1c6bb365e2da8cb4e7fad32d7dafdb', '1219647676b359a5ad0be6d9dda81c73']
    # NEXUS
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f nexus')
    assert hf.buddy2hash(tester) in ['f93607e234441a2577fa7d8a387ef7ec', '42bfddd38fa4ed75a99841abf2112e54']
    # PHYLIPI
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f phylipi')
    assert hf.buddy2hash(tester) in ['09dd492fde598670d7cfee61d4e2eab8', '438e1551b3f1c8526fc8a44eaf2a3dc1']
    # PHYLIPS
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f phylips')
    assert hf.buddy2hash(tester) in ['249c88cb64d41c47388514c65bf8fff1', '6366e50da5a6b33d2d281d6ea13df0b7']
    # Multi-param
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'pagan', '-f nexus --translate')
    assert hf.buddy2hash(tester) == 'dd140ec4eb895ce75d574498a58aa28a'

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
@br.skip_windows
def test_prank_inputs(sb_resources):
    # FASTA
    tester = Sb.pull_recs(sb_resources.get_one("d f"), 'α1')
    tester = Alb.generate_msa(tester, 'prank', '-once')
    assert tester.out_format == 'fasta'


@br.skip_windows
def test_prank_outputs1(sb_resources):
    # NEXUS
    tester = Sb.pull_recs(sb_resources.get_one("d f"), 'α1')
    tester = Alb.generate_msa(tester, 'prank', '-f=nexus -once')
    assert tester.out_format == 'nexus'


@br.skip_windows
def test_prank_outputs2(sb_resources):
    # PHYLIPI
    tester = Sb.pull_recs(sb_resources.get_one("d f"), 'α1')
    tester = Alb.generate_msa(tester, 'prank', params='-f=phylipi -once')
    assert tester.out_format == 'phylip-relaxed'


@br.skip_windows
def test_prank_outputs3(sb_resources):
    # PHYLIPS
    tester = Sb.pull_recs(sb_resources.get_one("d f"), 'α1')
    tester = Alb.generate_msa(tester, 'prank', params='-f=phylips -once')
    assert tester.out_format == 'phylipsr'


# ##########   MUSCLE   ########## #
@br.skip_windows
def test_muscle_inputs(sb_resources, hf):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'muscle')
    assert hf.buddy2hash(tester) == '5ec18f3e0c9f5cf96944a1abb130232f'


@br.skip_windows
def test_muscle_outputs(sb_resources, hf):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'muscle', '-clw')
    assert hf.buddy2hash(tester) == '91542667cef761ccaf39d8cb4e877944'


@br.skip_windows
def test_muscle_multi_param(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'muscle', '-clw -diags')
    assert hf.buddy2hash(tester) == '91542667cef761ccaf39d8cb4e877944'


# ##########   CLUSTAL   ########## #
clustalw_bin = 'clustalw' if which('clustalw') else 'clustalw2'


@br.skip_windows
def test_clustalw_inputs(sb_resources, hf):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalw_bin)
    assert hf.buddy2hash(tester) in ['955440b5139c8e6d7d3843b7acab8446', 'efc9b04f73c72036aa230a8d72da228b']


@br.skip_windows
def test_clustalw_outputs1(sb_resources, hf):
    # NEXUS
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalw_bin, '-output=nexus')
    assert hf.buddy2hash(tester) in ['f4a61a8c2d08a1d84a736231a4035e2e', '7ca92360f0787164664843c895dd98f2']


@br.skip_windows
def test_clustalw_outputs2(sb_resources, hf):
    # PHYLIP
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalw_bin, '-output=phylip')
    assert hf.buddy2hash(tester) in ['a9490f124039c6a2a6193d27d3d01205', 'd3cc272a45fbde4b759460faa8e63ebc']


@br.skip_windows
def test_clustalw_outputs3(sb_resources, hf):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalw_bin, '-output=fasta')
    assert hf.buddy2hash(tester) in ['955440b5139c8e6d7d3843b7acab8446', 'efc9b04f73c72036aa230a8d72da228b']


@br.skip_windows
def test_clustalw_multi_param(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalw_bin, '-output=phylip -noweights')
    assert hf.buddy2hash(tester) == 'ae9126eb8c482a82d4060d175803c478'


# ##########   CLUSTAL Omega   ########## #
@br.skip_windows
def get_clustalo_version():
    clustalo_bin = 'clustalo' if which('clustalo') else 'clustalomega'
    clustalo_version = Popen("{0} --version".format(clustalo_bin), shell=True,
                             stdout=PIPE).communicate()[0].decode().strip()
    if clustalo_version not in ["1.2.4", "1.2.3", "1.2.2", "1.2.1", "1.2.0", "1.0.3"]:
        raise ValueError("Untested CLustalO version (%s). Please update the tests as necessary." % clustalo_version)
    return clustalo_bin


clustalo_bin = get_clustalo_version()


@br.skip_windows
def test_clustalomega_inputs1(sb_resources, hf):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalo_bin)
    assert hf.buddy2hash(tester) in ['f5afdc7c76ab822bdc95230329766aba', 'e2f0efc90372b23a1c753629d43a79c4']


@br.skip_windows
def test_clustalomega_inputs2(sb_resources, hf):
    # PHYLIP
    tester = sb_resources.get_one("d py")
    tester = Alb.generate_msa(tester, clustalo_bin)
    assert hf.buddy2hash(tester) in ['8299780bf9485b89a2f3462ead666142', '5d808493da6c0f43a572b2a9257dce4f']


@br.skip_windows
def test_clustalomega_inputs3(sb_resources, hf):
    # STOCKHOLM
    tester = sb_resources.get_one("d s")
    tester = Alb.generate_msa(tester, clustalo_bin)
    assert hf.buddy2hash(tester) in ['8983529dd432dc9bb7f9b1a8acb64b18', 'd6654e3db3818cc3427cb9241113fdfa',
                                     'aeb2c5926843402cf620299802946224']


@br.skip_windows
def test_clustalomega_outputs1(sb_resources, hf):
    # CLUSTAL
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalo_bin, '--outfmt=clustal')
    assert hf.buddy2hash(tester) in ['970f6e4389f77a30563763937a3d32bc', '50c48b43089f528714c7933cd0f3f91c']


@br.skip_windows
def test_clustalomega_outputs2(sb_resources, hf):
    # PHYLIP
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalo_bin, '--outfmt=phylip')
    assert hf.buddy2hash(tester) in ['692c6af848bd90966f15908903894dbd', '75ec68313368dac249b40fe63b26777e']


@br.skip_windows
def test_clustalomega_outputs3(sb_resources, hf):
    # STOCKHOLM
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalo_bin, '--outfmt=stockholm')
    assert hf.buddy2hash(tester) in ['4c24975c033abcf15911a61cb9663a97', 'c115d474a16c23ca4219fa3d2fc9f154']


@br.skip_windows
def test_clustalomega_multi_param(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, clustalo_bin, '--outfmt=clustal --iter=1')
    assert hf.buddy2hash(tester) in ['cac11c9b4e2381f4d62af940028d7fe4', '25480f7a9340ff643bb7eeb326e8f981',
                                     '9c55bb8be0cc89c5346d8f699e97cc87']


# ##########   MAFFT   ########## #
@br.skip_windows
def test_mafft_inputs(sb_resources, hf):
    # FASTA
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'mafft')
    assert hf.buddy2hash(tester) == 'f94e0fd591dad83bd94201f0af038904'


@br.skip_windows
def test_mafft_outputs(sb_resources, hf):
    # CLUSTAL
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'mafft', '--clustalout')
    assert hf.buddy2hash(tester) == 'd6046c77e2bdb5683188e5de653affe5'


@br.skip_windows
def test_mafft_multi_param(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    tester = Alb.generate_msa(tester, 'mafft', '--clustalout --noscore')
    assert hf.buddy2hash(tester) == 'd6046c77e2bdb5683188e5de653affe5'


@br.skip_windows
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

    monkeypatch.setattr(br, "ask", ask_false)
    with pytest.raises(SystemExit):
        Alb.generate_msa(tester, clustalo_bin, keep_temp="%s%sga_temp_files" % (temp_dir.path, os.sep))

    monkeypatch.setattr(br, "ask", ask_true)
    Alb.generate_msa(tester, clustalo_bin, keep_temp="%s%sga_temp_files" % (temp_dir.path, os.sep))
    assert os.path.isfile("{0}{1}ga_temp_files{1}result".format(temp_dir.path, os.sep))
    assert os.path.isfile("{0}{1}ga_temp_files{1}tmp.fa".format(temp_dir.path, os.sep))


@br.skip_windows
def test_generate_alignments_genbank(sb_resources, hf):
    tester = sb_resources.get_one("p g")
    tester = Alb.generate_msa(tester, "mafft")
    assert hf.buddy2hash(tester) == "a4ab6b2a2ddda38a4d04abc18c54d18b"


@br.skip_windows
def test_generate_alignments_edges1(sb_resources):
    tester = sb_resources.get_one("d f")

    with pytest.raises(AttributeError) as e:
        Alb.generate_msa(tester, "foo")
    assert "foo is not a recognized alignment tool. Please check your spelling" in str(e)

    # noinspection PyUnresolvedReferences
    with mock.patch.dict('os.environ'):
        del os.environ['PATH']
        with pytest.raises(SystemError) as err:
            Alb.generate_msa(tester, "mafft")
        assert "#### Could not find mafft on your system. ####" in str(err)


@br.skip_windows
def get_args():
    return [("prank", "-f=phylipi"), (clustalo_bin, "--outfmt=foo"), (clustalw_bin, "-output=foo"),
            ("prank", "-f=nexus"), ("prank", "-f=foo")]


args = get_args()


@br.skip_windows
@pytest.mark.parametrize("tool,params", args)
def test_generate_alignments_edges2(tool, params, sb_resources):
    tester = sb_resources.get_one("d f")
    tester = Sb.pull_recs(tester, "α[2345]")
    Alb.generate_msa(tester, tool, params, quiet=True)


@br.skip_windows
def get_hmmer_version():
    hmmer_version = Popen("hmmbuild -h", shell=True, stdout=PIPE).communicate()[0].decode()
    hmmer_version = re.search("# HMMER (.*?) \(", hmmer_version).group(1)
    if hmmer_version not in ["3.1b2"]:
        raise ValueError("Untested HMMER version (%s). Please update the tests as necessary." % hmmer_version)


get_hmmer_version()


@br.skip_windows
def test_generate_hmm(alb_resources, hf, monkeypatch):
    tester = alb_resources.get_one("m p c")
    tester = Alb.generate_hmm(tester)
    for align in tester.alignments:
        align.hmm = re.search("(HMM +A +C +D.+//)", align.hmm, re.DOTALL).group(1)
    assert "  COMPO   2.68250  3.88919  3.04853  2.78121  3.08118  3.13138  3.72607  2.65113  2.67024  2.34849  " \
           "3.43272  3.05512  3.56890  3.09868  3.03335  2.74953  2.90269  2.58958  4.30351" in tester.alignments[0].hmm
    assert "COMPO   2.61975  3.93095  3.12640  2.80659  3.03969  2.94881  3.78599  2.73397  2.73613  2.36723  " \
           "3.48106  3.11755  3.38828  3.15135  3.06078  2.68581  2.82442  2.59321  4.24683" in tester.alignments[1].hmm

    tester = alb_resources.get_one("m d c")
    tester = Alb.generate_hmm(tester, "hmmbuild")
    for align in tester.alignments:
        align.hmm = re.search("(HMM +A +C +G.+//)", align.hmm, re.DOTALL).group(1)
    assert """\
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   1.37149  1.47979  1.42806  1.27722
          1.38629  1.38629  1.38629  1.38629
          0.10249  4.35641  2.47000  1.46634  0.26236  0.00000        *""" in tester.alignments[0].hmm
    assert """\
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   1.26618  1.65166  1.51488  1.18245
          0.93669  1.43354  1.84356  1.55419
          0.72237  1.59420  1.16691  3.55520  0.02899  0.00000        *""" in tester.alignments[1].hmm

    with pytest.raises(SystemError) as err:
        Alb.generate_hmm(tester, "foo")
    assert "Could not find foo on your system." in str(err)

    monkeypatch.setattr(Alb, "Popen", MockPopen)

    with pytest.raises(SystemError) as err:
        Alb.generate_hmm(tester)
    assert "No output detected after running hmmbuild." in str(err)
