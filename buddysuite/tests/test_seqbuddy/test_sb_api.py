#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" tests basic functionality of AlignBuddy class """
import pytest
from Bio.SeqFeature import FeatureLocation, CompoundLocation, Reference, SeqFeature
from Bio.Seq import Seq
from unittest import mock
import os
from os.path import join
import re
import urllib.request
import suds.client
import shutil
import time
from collections import OrderedDict

import SeqBuddy as Sb
import buddy_resources as br


TEMPDIR = br.TempDir()


# ##################### '-amd', '--amend_metadata' ###################### ##
def test_amend_metadata_desc(sb_resources, hf, monkeypatch):
    tester = sb_resources.get_one("p g")
    monkeypatch.setattr(br, "clean_regex", lambda regex: [regex])
    tester = Sb.amend_metadata(tester, "description", "foo", "ML")
    assert hf.buddy2hash(tester) == "988a005283cb6b9de73e18d52631fb79"


def test_amend_metadata_topo(sb_resources, hf, monkeypatch, capsys):
    monkeypatch.setattr(br, "clean_regex", lambda regex: [regex])
    tester = sb_resources.get_one("p g")
    for value, _hash in [('linear', 'b3feec305beee8243be7f7afb02a909d'),
                         ('circular', '7ff89b41edab0a10b53629a28092e062'),
                         ('', '0a8462e72f64fcd22544bb153b51b2b6')]:
        tester = Sb.amend_metadata(tester, "topology", value, ".*")
        assert hf.buddy2hash(tester) == _hash

    with pytest.raises(ValueError) as err:
        Sb.amend_metadata(tester, "topology", "foo_bar", ".*")

    assert "Topology values are limited to ['', 'linear', 'circular']" in str(err)


def test_amend_metadata_str_attr(sb_resources, hf, monkeypatch):
    monkeypatch.setattr(br, "clean_regex", lambda regex: [regex])

    """
    There are limited number of values that biopython will accept
    ['BCT', 'CON', 'ENV', 'EST', 'FUN', 'GSS', 'HTC', 'HTG', 'HUM', 
     'INV', 'MAM', 'MUS', 'PAT', 'PHG', 'PLN', 'PRI', 'PRO', 'ROD', 
     'STS', 'SYN', 'TGN', 'UNA', 'UNC', 'VRL', 'VRT', 'XXX']
     """
    # data_file_division
    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "data_file_division", "ROD", ".*")
    assert hf.buddy2hash(tester) == "a834ade7ffdc453b9c61817c9138a550"

    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "data_file_division", "FOO", ".*")
    assert hf.buddy2hash(tester) == "6d767ae1bdda9ce23ce99cfae35a1a74"

    # date
    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "date", "21-APR-2000", ".*")
    assert hf.buddy2hash(tester) == "3695432066ca4c61727aa139d40a7b8e"

    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "date", "FOO", ".*")
    assert hf.buddy2hash(tester) == "28b8f6425b92755f2936883837b1c452"

    # source
    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "source", "Mnemiopsis", ".*")
    assert hf.buddy2hash(tester) == "ff0d55dc060e28e47440a07414ea62bf"

    # organism
    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "organism", "Foo", "leidyi")
    assert hf.buddy2hash(tester) == "0494f26e437928e8af96ba6e7cbfb7cf"

    # comment
    tester = sb_resources.get_one("p g")
    structured_comment = OrderedDict()
    structured_comment["Hello there"] = OrderedDict()
    structured_comment["Hello there"]["something inner"] = "My text for testing"
    tester.records[0].annotations["structured_comment"] = structured_comment
    tester = Sb.amend_metadata(tester, "comment", "data", "text")
    assert hf.buddy2hash(tester) == "5f919a25b9735e47ccc0f492328120c0"


def test_amend_metadata_refs(sb_resources, hf, monkeypatch):
    tester = sb_resources.get_one("p g")
    monkeypatch.setattr(br, "clean_regex", lambda regex: [regex])
    reference = Reference()
    reference.authors = "Bond SR, Keat KE, Barreira SN, Baxevanis AD"
    reference.title = "BuddySuite: Command-Line Toolkits for Manipulating Sequences, Alignments, and Phylogenetic Trees"
    reference.journal = "Mol Biol Evol."
    reference.pubmed_id = '28333216'
    reference.comment = 'Hurray for published papers!'
    for rec in tester.records:
        rec.annotations["references"] = [reference]
    tester = Sb.amend_metadata(tester, "references", "https://github.com/biologyguy/BuddySuite", "Hurray.*")
    assert hf.buddy2hash(tester) == "b131fd4a1781403a9b870ac4bde2cf3d"

    tester = Sb.amend_metadata(tester, "references", "", "")
    assert hf.buddy2hash(tester) == "b131fd4a1781403a9b870ac4bde2cf3d"


def test_amend_metadata_list_attr(sb_resources, hf, monkeypatch):
    monkeypatch.setattr(br, "clean_regex", lambda regex: [regex])

    # taxonomy
    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "taxonomy", "Eukaryota Opisthokonta Metazoa Eumetazoa Ctenophora Tentaculata"
                                                   " Lobata Bolinopsidae Mnemiopsis", ".*")
    assert hf.buddy2hash(tester) == "cdc30e7ec65c525bac898bbfaa75a0b7"

    tester = Sb.amend_metadata(tester, "taxonomy", "FooBar", "Eukaryota")
    assert hf.buddy2hash(tester) == "bb80f495ade4e8a8f62253642f350d22"

    # keywords
    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "keywords", "Something Else", ".*")
    assert hf.buddy2hash(tester) == "bda9a1405d954aa2542d9017ba5c0796"

    tester.records[0].annotations["keywords"] = []
    tester = Sb.amend_metadata(tester, "keywords", "FooBar", "Else")
    assert hf.buddy2hash(tester) == "809663ca6a6f92772a02301b0ab901a7"


def test_amend_metadata_dbxrefs(sb_resources, hf, monkeypatch):
    monkeypatch.setattr(br, "clean_regex", lambda regex: [regex])
    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "dbxrefs", "Project:1234 Project:4321", ".*")
    assert hf.buddy2hash(tester) == "342286427f1dc495131610a9c02587cf"

    tester.records[0].dbxrefs = []
    tester = Sb.amend_metadata(tester, "dbxrefs", "Activity", "Project")
    assert hf.buddy2hash(tester) == "65d5213aebe744763d3662eb57bbd514"


def test_amend_metadata_version(sb_resources, hf, monkeypatch):
    monkeypatch.setattr(br, "clean_regex", lambda regex: [regex])
    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "version", 5, ".*")
    assert hf.buddy2hash(tester) == "d83c18c44529700eaa3ab32da3ec8d08"

    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "sequence_version", 5, ".*")
    assert hf.buddy2hash(tester) == "d83c18c44529700eaa3ab32da3ec8d08"

    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "sequence_version", "foo", ".*")
    assert hf.buddy2hash(tester) == "0a8462e72f64fcd22544bb153b51b2b6"


def test_amend_metadata_arb_attr(sb_resources, hf, monkeypatch):
    monkeypatch.setattr(br, "clean_regex", lambda regex: [regex])
    tester = sb_resources.get_one("p g")
    tester = Sb.amend_metadata(tester, "id", "Hech", "Mle")
    assert hf.buddy2hash(tester) == "c0a9461be6cb3b7e6ad76200f009c4ab"

    tester = Sb.amend_metadata(tester, "foobar", "Some stuff", ".*")
    assert tester.records[0].foobar == "Some stuff"


# ##################### '-ano', '--annotate' ###################### ##
def test_annotate_pattern(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), pattern='α4')
    assert hf.buddy2hash(tester) == '4f9b3f82ce3c8d6c953f60a5a0e9574e'


def test_annotate_no_pattern(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', (1, 100))
    assert hf.buddy2hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_compoundlocation(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', [(1, 100), (200, 250)])
    assert hf.buddy2hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_list_str(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', ['1-100', '(200-250)'])
    assert hf.buddy2hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_str(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', '1-100, (200-250)')
    assert hf.buddy2hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_fl_obj(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', FeatureLocation(start=0, end=100))
    assert hf.buddy2hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_cl_obj(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', CompoundLocation([FeatureLocation(start=0, end=100),
                                                                   FeatureLocation(start=199, end=250)],
                                                                  operator='order'))
    assert hf.buddy2hash(tester) == 'd22c44cf1a53624b58a86b0fb98c33a6'


def test_annotate_typerror(sb_resources):
    with pytest.raises(TypeError):
        tester = sb_resources.get_one("d g")
        Sb.annotate(tester, 'misc_feature', 5)


def test_annotate_location_invalid(sb_resources):
    with pytest.raises(AttributeError) as err:
        tester = sb_resources.get_one("d g")
        Sb.annotate(tester, 'misc_feature', ["foo", 100])
    assert "The provided location string is invalid" in str(err)

    with pytest.raises(AttributeError) as err:
        tester = sb_resources.get_one("d g")
        Sb.annotate(tester, 'misc_feature', "foo-100")
    assert "The provided location string is invalid" in str(err)


def test_annotate_pos_strand(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand='+')
    assert hf.buddy2hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_neg_strand(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand='-')
    assert hf.buddy2hash(tester) == '08524707f09d7eb273775f791d92964c'


def test_annotate_no_strand(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), strand=0)
    assert hf.buddy2hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_qualifier_dict(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers=OrderedDict([('foo', 'bar'), ('hello', 'world')]))
    assert hf.buddy2hash(tester) == '34e9dfb9cfe62f0a4657c977eda45688'


def test_annotate_qualifier_list(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers=['foo=bar', 'hello=world'])
    assert hf.buddy2hash(tester) == '34e9dfb9cfe62f0a4657c977eda45688'


def test_annotate_qualifier_error(sb_resources):
    tester = sb_resources.get_one("d g")
    with pytest.raises(TypeError):
        Sb.annotate(tester, 'misc_feature', (1, 100), qualifiers=tuple)


def test_annotate_out_of_range(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', [(-10, 100), (200, 10000)])
    assert hf.buddy2hash(tester) == 'a8f90863b2bbeaa519e8230a187532ca'

    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', [(1, 10000)])
    assert hf.buddy2hash(tester) == 'f5c90e3458fbca9b9565dac7877cc248'

    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', [(1, 10000), (20000, 30000)])
    assert hf.buddy2hash(tester) == 'f5c90e3458fbca9b9565dac7877cc248'

    tester = sb_resources.get_one("d g")
    tester = Sb.annotate(tester, 'misc_feature', FeatureLocation(start=-10, end=100))
    assert hf.buddy2hash(tester) == '83a975003d5bb5987c6f2cbaeaa75cf7'


def test_annotate_protein(sb_resources, hf):
    tester = sb_resources.get_one("p g")
    tester = Sb.annotate(tester, 'misc_feature', (1, 100))
    assert hf.buddy2hash(tester) == 'e350f5d02a970332344715a6822c4ab5'


def test_annotate_unrec_strand(capsys, sb_resources):
    tester = sb_resources.get_one("d g")
    Sb.annotate(tester, 'misc_feature', (1, 100), strand='foo')
    out, err = capsys.readouterr()
    assert err == "Warning: strand input not recognized. Value set to None."
    
    
# ######################  '-asl', '--ave_seq_length' ###################### #
def test_ave_seq_length_dna(sb_resources):
    for seqbuddy in sb_resources.get_list("d f g n pr s"):
        assert round(Sb.ave_seq_length(seqbuddy, clean=True), 2) == 1285.15


def test_ave_seq_length_pep(sb_resources):
    for seqbuddy in sb_resources.get_list("p f g n pr s"):
        assert round(Sb.ave_seq_length(seqbuddy, clean=True), 2) == 427.38


# ######################  '-btr', '--back_translate' ###################### #
# Only fasta and genbank
hashes = [('p f', 'human', '1b14489a78bfe8255c777138877b9648', '5e42effd0bb67445200263f4bf59dbeb'),
          ('p g', 'human', '0899fb80a7c7cdd0abb3c839ff9c41b6', 'c717f0c902793d39490ae9f73e53990b'),
          ('p f', 'yeast', '859ecfb88095f51bfaee6a1d1abeb50f', '114b6c4f0a29db275d7832318870bcae'),
          ('p g', 'yeast', 'a98c276af6ed064feb84f9eee8676449', 'ccab4e7618f9b04704c3fd6e52cca47f'),
          ('p f', 'ecoli', '952a91a4506afb57f27136aa1f2a8af9', '592d09df34a6c34caa2d0ff9dbc54930'),
          ('p g', 'ecoli', '64bb1cfb9be3bd96affa13de525e3867', '5a06b4d87e487bb9c460589d920f9837'),
          ('p f', 'mouse', '3a3ee57f8dcde25c99a655494b218928', 'eb29f87da4d022e276e19a71aa03c469'),
          ('p g', 'mouse', 'ef5420f9f0b459b9cff1bd2ebb1a3c64', '193e2d7dde93c93d8d83ea40cd99fb25')]


@pytest.mark.parametrize("key,organism,opt_hash,rand_hash", hashes)
def test_back_translate(key, organism, opt_hash, rand_hash, sb_resources, hf):
    tester = sb_resources.get_one(key)
    tester = Sb.back_translate(tester, mode='OPTIMIZED', species=organism)
    assert hf.buddy2hash(tester) == opt_hash

    tester = sb_resources.get_one(key)
    tester = Sb.back_translate(tester, species=organism, r_seed=12345)
    assert hf.buddy2hash(tester) == rand_hash


def test_back_translate_nucleotide_exception(sb_resources):
    with pytest.raises(TypeError):
        Sb.back_translate(sb_resources.get_one("d g"))


def test_back_translate_bad_mode(sb_resources):
    with pytest.raises(AttributeError):
        Sb.back_translate(sb_resources.get_one("p f"), 'fgsdjkghjdalgsdf', 'human')


def test_back_translate_bad_organism(sb_resources):
    seqbuddy = sb_resources.get_one("p f")
    with pytest.raises(AttributeError):
        Sb.back_translate(seqbuddy, 'OPTIMIZED', 'fgsdjkghjdalgsdf')


# ######################  '-bl2s', '--bl2seq' ###################### #
# ToDo: Mock output data
def test_bl2seq():
    pass


def test_bl2_no_binary(sb_resources):
    # noinspection PyUnresolvedReferences
    with mock.patch.dict(os.environ, {"PATH": ""}):
        with mock.patch('builtins.input', return_value="n"):
            with pytest.raises(RuntimeError):
                Sb.bl2seq(sb_resources.get_one("d f"))

            with pytest.raises(RuntimeError):
                Sb.bl2seq(sb_resources.get_one("p f"))


# ######################  '-bl', '--blast' ###################### #
class MockPopen(object):
    def __init__(self, command, *_, **__):
        self.command = command

    def communicate(self):
        if self.command[:11] == "makeblastdb":
            if "nucl" in self.command:
                output = ["""\
Building a new DB, current time: 01/25/2017 10:35:34
New DB name:   Mnemiopsis_cds
New DB title:  Mnemiopsis_cds.fa
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 12 sequences in 0.000432968 seconds.""".encode()]
            else:
                output = ["""\
Building a new DB, current time: 01/25/2017 10:35:34
New DB name:   {0}{1}Mnemiopsis_pep
New DB title:  {0}{1}Mnemiopsis_pep.fa
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 12 sequences in 0.000432968 seconds.""".encode()]
        elif "blast_error" in self.command:
            output = ["", "Some sort of Error".encode()]
        elif self.command[:6] == "blastn":
            output = ["".encode(), "".encode()]
            out_dir = re.search("-out (.*out.txt)", self.command)
            with open(out_dir.group(1), "w") as ofile:
                ofile.write("Mle-Panx12\tem1vubDvb9\t100.000\t1209\t0\t0\t1\t1209\t1\t1209\t0.0\t2233\n"
                            "Mle-Panx12\towIvkyuadd\t100.000\t1209\t0\t0\t1\t1209\t1\t1209\t0.0\t2233\n"
                            "Mle-Panx2\towIvkyuadd\t100.000\t1314\t0\t0\t1\t1314\t1\t1314\t0.0\t2427\n")
        elif self.command[:10] == "blastdbcmd":
            if "em1vubDvb9" in self.command:
                output = [">em1vubDvb9 cDNA - ML25997a.\nATGGTTATTGACATCCTCTCCGGTTTTAAGGGGATCACGC".encode("utf-8")]
            else:
                output = [">owIvkyuadd cDNA - ML25998a.\nATGGTATTGGATCTCATTTCTGGAAGCTTGCATCACACGA".encode("utf-8")]
        else:
            output = []

        return output


def test_blast(monkeypatch, capsys, sb_resources, hf):
    tmp_dir = br.TempDir()
    for extension in ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]:
        shutil.copyfile("%sblast%sMnemiopsis_cds.%s" % (hf.resource_path, os.sep, extension),
                        "%s%squery_db.%s" % (tmp_dir.path, os.sep, extension))

    monkeypatch.setattr(Sb, "_check_for_blast_bin", lambda _: True)
    monkeypatch.setattr(br, "TempDir", lambda: tmp_dir)
    monkeypatch.setattr(Sb, "Popen", MockPopen)

    subject = sb_resources.get_one("d f")
    subject = Sb.pull_recs(subject, "2")

    query = sb_resources.get_one("d g")
    blast = Sb.blast(Sb.make_copy(subject), query)
    out, err = capsys.readouterr()

    assert "Building a new DB with makeblastdb" in err
    assert "Sequence type: Nucleotide" in err
    assert "New DB" not in err

    assert "blastn -num_threads 4 -evalue 0.01" in err

    assert """\
# ######################## BLAST results ######################## #
Mle-Panx12	""" in err, print(re.findall("Mle-Panx12.*", str(err)))

    assert """\
# ######################## BLAST results ######################## #
Mle-Panx12	Mle-Panxα12	100.000	1209	0	0	1	1209	1	1209	0.0	2233
Mle-Panx12	Mle-Panxα2	100.000	1209	0	0	1	1209	1	1209	0.0	2233
Mle-Panx2	Mle-Panxα2	100.000	1314	0	0	1	1314	1	1314	0.0	2427
# ############################################################### #""" in err, print(str(err))
    assert str(blast) == """\
>Mle-Panxα12 cDNA - ML25997a.
ATGGTTATTGACATCCTCTCCGGTTTTAAGGGGATCACGC
>Mle-Panxα2 cDNA - ML25998a.
ATGGTATTGGATCTCATTTCTGGAAGCTTGCATCACACGA
"""

    # This is calling the pre-indexed database instead of a sequence file
    expected_output = """\
>em1vubDvb9 cDNA - ML25997a.
ATGGTTATTGACATCCTCTCCGGTTTTAAGGGGATCACGC
>owIvkyuadd cDNA - ML25998a.
ATGGTATTGGATCTCATTTCTGGAAGCTTGCATCACACGA
"""

    query = "%s%squery_db" % (tmp_dir.path, os.sep)
    assert str(Sb.blast(Sb.make_copy(subject), query)) == expected_output

    query = "%s%squery_db.n" % (tmp_dir.path, os.sep)
    assert str(Sb.blast(Sb.make_copy(subject), query)) == expected_output

    query = "%s%squery_db.nhr" % (tmp_dir.path, os.sep)
    assert str(Sb.blast(Sb.make_copy(subject), query)) == expected_output


def test_blast_quiet(monkeypatch, capsys, sb_resources, hf):
    tmp_dir = br.TempDir()
    for extension in ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]:
        shutil.copyfile("%sblast%sMnemiopsis_cds.%s" % (hf.resource_path, os.sep, extension),
                        "%s%squery_db.%s" % (tmp_dir.path, os.sep, extension))

    monkeypatch.setattr(Sb, "_check_for_blast_bin", lambda _: True)
    monkeypatch.setattr(br, "TempDir", lambda: tmp_dir)
    monkeypatch.setattr(Sb, "Popen", MockPopen)

    subject = sb_resources.get_one("d f")
    subject = Sb.pull_recs(subject, "2")
    query = sb_resources.get_one("d g")

    blast = Sb.blast(Sb.make_copy(subject), Sb.make_copy(query), quiet=True)
    out, err = capsys.readouterr()
    assert err == ""
    assert str(blast) == """\
>Mle-Panxα12 cDNA - ML25997a.
ATGGTTATTGACATCCTCTCCGGTTTTAAGGGGATCACGC
>Mle-Panxα2 cDNA - ML25998a.
ATGGTATTGGATCTCATTTCTGGAAGCTTGCATCACACGA
"""


def test_blast_blacklist_args(monkeypatch, capsys, sb_resources, hf):
    tmp_dir = br.TempDir()
    for extension in ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]:
        shutil.copyfile("%sblast%sMnemiopsis_cds.%s" % (hf.resource_path, os.sep, extension),
                        "%s%squery_db.%s" % (tmp_dir.path, os.sep, extension))

    monkeypatch.setattr(Sb, "_check_for_blast_bin", lambda _: True)
    monkeypatch.setattr(br, "TempDir", lambda: tmp_dir)
    monkeypatch.setattr(Sb, "Popen", MockPopen)

    subject = sb_resources.get_one("d f")
    subject = Sb.pull_recs(subject, "2")
    query = sb_resources.get_one("d g")

    Sb.blast(Sb.make_copy(subject), Sb.make_copy(query),
             blast_args='-db foo -query bar -subject baz -out there -outfmt 7 -not_black_listed 54')

    out, err = capsys.readouterr()
    black_list = ["db", "query", "subject", "out", "outfmt"]
    for no_no in black_list:
        assert "Warning: Explicitly setting the blast -%s parameter is not supported in SeqBuddy.\n" % no_no

    assert "blastn -num_threads 4 -evalue 0.01 -not_black_listed 54" in err


def test_blast_num_threads(monkeypatch, capsys, sb_resources, hf):
    tmp_dir = br.TempDir()
    for extension in ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]:
        shutil.copyfile("%sblast%sMnemiopsis_cds.%s" % (hf.resource_path, os.sep, extension),
                        "%s%squery_db.%s" % (tmp_dir.path, os.sep, extension))

    monkeypatch.setattr(Sb, "_check_for_blast_bin", lambda _: True)
    monkeypatch.setattr(br, "TempDir", lambda: tmp_dir)
    monkeypatch.setattr(Sb, "Popen", MockPopen)

    subject = sb_resources.get_one("d f")
    subject = Sb.pull_recs(subject, "2")
    query = sb_resources.get_one("d g")

    Sb.blast(Sb.make_copy(subject), Sb.make_copy(query), blast_args="-num_threads 5")
    out, err = capsys.readouterr()
    assert "blastn -num_threads 5 -evalue 0.01 \n" in err

    with pytest.raises(ValueError) as err:
        Sb.blast(Sb.make_copy(subject), Sb.make_copy(query), blast_args="-num_threads foo")
    assert "-num_threads expects an integer." in str(err)


def test_blast_evalue(monkeypatch, capsys, sb_resources, hf):
    tmp_dir = br.TempDir()
    for extension in ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]:
        shutil.copyfile("%sblast%sMnemiopsis_cds.%s" % (hf.resource_path, os.sep, extension),
                        "%s%squery_db.%s" % (tmp_dir.path, os.sep, extension))

    monkeypatch.setattr(Sb, "_check_for_blast_bin", lambda _: True)
    monkeypatch.setattr(br, "TempDir", lambda: tmp_dir)
    monkeypatch.setattr(Sb, "Popen", MockPopen)

    subject = sb_resources.get_one("d f")
    subject = Sb.pull_recs(subject, "2")
    query = sb_resources.get_one("d g")

    Sb.blast(Sb.make_copy(subject), Sb.make_copy(query), blast_args="-evalue 0.1")
    out, err = capsys.readouterr()
    assert "blastn -num_threads 4 -evalue 0.1 \n" in err

    with pytest.raises(ValueError) as err:
        Sb.blast(Sb.make_copy(subject), Sb.make_copy(query), blast_args="-evalue foo")
    assert "-evalue expects a number." in str(err)


def test_blast_errors(monkeypatch, sb_resources, hf):
    # Note that the order of these tests is important because of all the wonky monkeypatching
    tmp_dir = br.TempDir()
    for extension in ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]:
        shutil.copyfile("%sblast%sMnemiopsis_cds.%s" % (hf.resource_path, os.sep, extension),
                        "%s%squery_db.%s" % (tmp_dir.path, os.sep, extension))

    monkeypatch.setattr(Sb, "_check_for_blast_bin", lambda _: True)
    monkeypatch.setattr(br, "TempDir", lambda: tmp_dir)
    monkeypatch.setattr(Sb, "Popen", MockPopen)

    subject = sb_resources.get_one("d f")
    subject = Sb.pull_recs(subject, "2")
    query = sb_resources.get_one("d g")

    with pytest.raises(RuntimeError) as err:
        Sb.blast(Sb.make_copy(subject), Sb.make_copy(query), blast_args="-blast_error 0.1")
    assert "Some sort of Error" in str(err)

    os.remove("%s%squery_db.nhr" % (tmp_dir.path, os.sep))
    with pytest.raises(RuntimeError) as err:
        Sb.blast(subject, query, blast_args="-delete_nhr 0.1")
    assert "The .nhr file of your BLASTN database was not found." in str(err)

    query = sb_resources.get_one("p g")
    with pytest.raises(ValueError) as err:
        Sb.blast(subject, query)
    assert "Trying to compare protein to nucleotide." in str(err)

    monkeypatch.setattr(Sb, "_check_for_blast_bin", lambda program: False if program == "makeblastdb" else True)
    with pytest.raises(SystemError) as err:
        Sb.blast(subject, query)
    assert "makeblastdb not found in system path." in str(err)

    monkeypatch.setattr(Sb, "_check_for_blast_bin", lambda program: False if program == "blastdbcmd" else True)
    with pytest.raises(SystemError) as err:
        Sb.blast(subject, query)
    assert "blastdbcmd not found in system path." in str(err)

    monkeypatch.setattr(Sb, "_check_for_blast_bin", lambda program: False if program == "blastn" else True)
    with pytest.raises(SystemError) as err:
        Sb.blast(subject, query)
    assert "blastn not found in system path." in str(err)


# ######################  '-cs', '--clean_seq'  ###################### #
def test_clean_seq_prot(sb_resources, hf):
    # Protein
    tester = sb_resources.get_one("p f")
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert hf.buddy2hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert hf.buddy2hash(tester) == "dc53f3be7a7c24425dddeea26ea0ebb5"


def test_clean_seq_dna(sb_odd_resources, hf):
    # DNA
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_dna"])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert hf.buddy2hash(tester) == "71b28ad2730a9849f2ba0f70e9e51a9f"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert hf.buddy2hash(tester) == "1912fadb5ec52a38ec707c58085b86ad"
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_dna"])
    tester = Sb.clean_seq(tester, ambiguous=False, rep_char="X")
    assert hf.buddy2hash(tester) == "4c10ba4474d7484652cb633f03db1be1"


def test_clean_seq_rna(sb_odd_resources, hf):
    # RNA
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_rna"])
    tester = Sb.clean_seq(tester, ambiguous=True)
    assert hf.buddy2hash(tester) == "cdb1b963536d57efc7b7f87d2bf4ad22"
    tester = Sb.clean_seq(tester, ambiguous=False)
    assert hf.buddy2hash(tester) == "e66b785c649ad5086bcefd22e9ef9b41"
    tester = Sb.SeqBuddy(sb_odd_resources["ambiguous_rna"])
    tester = Sb.clean_seq(tester, ambiguous=False, rep_char="X")
    assert hf.buddy2hash(tester) == "8ae19ab51b04076112d2f649353a4a79"


def test_clean_seq_align(sb_resources, capsys):
    # Alignment formats are converted to fasta to prevent errors with sequence lengths
    for tester in sb_resources.get_list("d n py pr s"):
        str(Sb.clean_seq(tester))
        out, err = capsys.readouterr()
        assert err == "Warning: Alignment format detected but sequences are different lengths. " \
                      "Format changed to fasta to accommodate proper printing of records.\n\n", print(err)


# ######################  '-cmp', '--complement' ###################### #
hashes = [('d f', 'e4a358ca57aca0bbd220dc6c04c88795'), ('d g', '3366fcc6ead8f1bba4a3650e21db4ec3'),
          ('d n', '365bf5d08657fc553315aa9a7f764286'), ('d py', '520036b49dd7c70b9dbf4ce4d2c0e1d8'),
          ('d pr', 'dc1f7a3769a1e0b007969db1ab405e89'), ('d s', '5891348e8659290c2355fabd0f3ba4f4')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_complement(key, next_hash, sb_resources, hf):
    tester = Sb.complement(sb_resources.get_one(key))
    assert hf.buddy2hash(tester) == next_hash


def test_complement_pep_exception(sb_resources):  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.complement(sb_resources.get_one("p f"))


# ######################  '-cts', '--concat_seqs' ###################### #
hashes = [('d f', '177c8a4bb83271a4380d98075455a436'), ('d g', '4bba7fbae1fd7a675ef5dda95683fba0'),
          ('d n', '81457d1d0089ec13ce3cbecd08a86489'), ('d py', 'c900e968aa9adae8865be374bd814954'),
          ('d pr', '81457d1d0089ec13ce3cbecd08a86489'), ('d s', '81457d1d0089ec13ce3cbecd08a86489'),
          ('p f', '4e56305350df0775b623f3f88dde71fc'), ('p g', 'fcc7414f710b6321547248668268bec7'),
          ('p n', '58834d7fb61ca6b106f8cc68a6bcfcf0'), ('p py', '7e4819287236bf41b95b7952cbcbbd05'),
          ('p pr', '1312e411c015f1e1e8339b0c6b492676'), ('p s', '58834d7fb61ca6b106f8cc68a6bcfcf0')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_concat_seqs(key, next_hash, sb_resources, hf):
    tester = Sb.concat_seqs(sb_resources.get_one(key))
    assert hf.buddy2hash(tester) == next_hash


hashes = [('d g', '4bba7fbae1fd7a675ef5dda95683fba0'), ('d n', '177c8a4bb83271a4380d98075455a436'),
          ('p g', 'e2fabdbdfac33e4af80ad1831524df59'), ('p py', '2d25d6ce7b4391faff110a74d27f074f')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_concat_seqs_clean(key, next_hash, sb_resources, hf):
    tester = Sb.concat_seqs(sb_resources.get_one(key), clean=True)
    tester.out_format = "gb"
    assert hf.buddy2hash(tester) == next_hash, print(tester)


# ##################### '-cc', 'count_codons' ###################### ##
def test_count_codons_dna(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    tester.records[0].seq = tester.records[0].seq[:-4]
    counter = Sb.count_codons(tester)
    assert hf.buddy2hash(counter[0]) == '4dbd9a5c68d85bb200c75b309fdaeeca'
    assert hf.string2hash(str(counter[1])) == '8d47313f3db02aee48575ff8ff4741b4'


def test_count_codons_rna(sb_resources, hf):
    tester = Sb.count_codons(sb_resources.get_one("r f"))[1]
    assert hf.string2hash(str(tester)) == 'b91daa8905533b5885d2067d9d6ffe36'


def test_count_codons_dna_badchar(sb_resources, hf):
    tester = Sb.count_codons(Sb.insert_sequence(sb_resources.get_one("d f"), 'PPP', -1))[1]
    assert hf.string2hash(str(tester)) == '9aba116675fe0e9eaaf43e5c6e0ba99d'


def test_count_codons_pep_exception(sb_resources):
    tester = sb_resources.get_one("p f")
    with pytest.raises(TypeError):
        Sb.count_codons(tester)


# ######################  '-cr', '--count_residues' ###################### #
def test_count_residues_unambig_dna(sb_resources):
    tester = Sb.SeqBuddy(">seq1\nACGCGAAGCGAACGCGCAGACGACGCGACGACGACGACGCA", in_format="fasta")
    tester = Sb.count_residues(tester).to_dict()
    assert tester["seq1"].buddy_data["res_count"]['A'] == [13, 0.3170731707317073]

    tester = sb_resources.get_one("d f")
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert res_count['G'] == [265, 0.21703521703521703]
    assert "% Ambiguous" not in res_count and "U" not in res_count


def test_count_residues_unambig_rna(sb_resources):
    tester = sb_resources.get_one("r f")
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert res_count['U'] == [356, 0.2915642915642916]
    assert "% Ambiguous" not in res_count and "T" not in res_count


def test_count_residues_ambig_dna(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources['ambiguous_dna'])
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert "U" not in res_count
    assert res_count['Y'] == [1, 0.000819000819000819]
    assert res_count['% Ambiguous'] == 0.98


def test_count_residues_ambig_rna(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources['ambiguous_rna'])
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert "T" not in res_count
    assert res_count['U'] == [353, 0.2891072891072891]
    assert res_count["Y"] == [1, 0.000819000819000819]
    assert res_count['% Ambiguous'] == 0.98


def test_count_residues_prot(sb_resources):
    tester = sb_resources.get_one("f p")
    Sb.count_residues(tester)
    res_count = tester.to_dict()['Mle-Panxα6'].buddy_data["res_count"]
    assert res_count['P'] == [17, 0.04176904176904177]
    assert res_count["G"] == [23, 0.056511056511056514]
    assert "% Ambiguous" not in res_count
    res_count = tester.to_dict()['Mle-Panxα8'].buddy_data["res_count"]
    assert res_count["% Ambiguous"] == 1.2
    assert res_count["% Positive"] == 12.23
    assert res_count["% Negative"] == 12.71
    assert res_count["% Uncharged"] == 73.62
    assert res_count["% Hydrophilic"] == 36.93
    assert res_count["% Hydrophobic"] == 55.4


# ######################  '-dgn' '--degenerate_sequence'################### #
def test_degenerate_sequence_without_arguments(sb_resources, hf):
    tester = sb_resources.get_one("f d")
    tester = Sb.degenerate_sequence(tester)
    assert hf.buddy2hash(tester) == '0638bc6546eebd9d50f771367d6d7855'


def test_degenerate_sequence_with_different_codon_tables(sb_resources, hf):
    for indx, next_hash in enumerate(['0638bc6546eebd9d50f771367d6d7855', '72373f8356051e2c6b67642451379054',
                                      '9172ad5947c0961b54dc5adbd03d4249', 'b45ac94ee6a98e495e115bfeb5bd9bcd',
                                      '76c45b4de8f7527b4139446b4551712b', 'baa5b48938cc5cae953c9083a5b21b12',
                                      '0ca67c4740fefbc7a20d806715c3ca12', 'd43ad8f328ff1d30eb1fb7bcd667a345',
                                      'd9d0f5cd8f0c25a0042527cc1cea802e', '4b9790f3f4eeeae1a9667b62b93bc961',
                                      '7ec4365c3571813d63cee4b70ba5dcf5']):
        tester = Sb.degenerate_sequence(sb_resources.get_one("d f"), table=(indx + 1))
        assert hf.buddy2hash(tester) == next_hash


def test_degenerate_sequence_edges(sb_resources, sb_odd_resources, hf):
    tester = sb_resources.get_one("p f")

    # Bad table reference
    with pytest.raises(KeyError) as e:
        Sb.degenerate_sequence(tester, 100)
    assert "Could not locate codon dictionary" in str(e.value)

    # Protein input
    with pytest.raises(TypeError) as e:
        Sb.degenerate_sequence(tester, 1)
    assert "Nucleic acid sequence required, not protein." in str(e.value)

    # RNA input
    tester = Sb.degenerate_sequence(Sb.SeqBuddy(sb_odd_resources["ambiguous_rna"]))
    assert hf.buddy2hash(tester) == "ff52e05971aeafd24c73a3b543901e4b"


# ######################  '-df', '--delete_features' ###################### #
def test_delete_features(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.delete_features(tester, 'donor')
    assert hf.buddy2hash(tester) == 'f84df6a77063c7def13babfaa0555bbf'


# ######################  '-dl', '--delete_large' ###################### #
def test_delete_large(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    tester = Sb.delete_large(tester, 1285)
    assert hf.buddy2hash(tester) == '25859dc69d46651a1e04a70c07741b35'


# ######################  '-dm', '--delete_metadata' ###################### #
hashes = [('d f', 'aa92396a9bb736ae6a669bdeaee36038'), ('d g', 'ad7ca097144843b8c13856e2a40afe09'),
          ('d n', 'cb1169c2dd357771a97a02ae2160935d'), ('d py', '503e23720beea201f8fadf5dabda75e4'),
          ('d pr', '52c23bd793c9761b7c0f897d3d757c12'), ('d s', 'a50943ccd028b6f5fa658178fa8cf54d'),
          ('p f', 'bac5dc724b1fee092efccd2845ff2513'), ('p g', '945327e8f1e483bb99f78540f507c4e1'),
          ('p n', '17ff1b919cac899c5f918ce8d71904f6'), ('p py', '968ed9fa772e65750f201000d7da670f'),
          ('p pr', 'ce423d5b99d5917fbef6f3b47df40513'), ('p s', 'e224c16f6c27267b5f104c827e78df33')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_delete_metadata(key, next_hash, sb_resources, hf):
    tester = Sb.delete_metadata(sb_resources.get_one(key))
    assert hf.buddy2hash(tester) == next_hash


# ######################  '-dr', '--delete_records' ###################### #
hashes = [('d f', '54bdb42423b1d331acea18218101e5fc'), ('d g', 'e2c03f1fa21fd27b2ff55f7f721a1a99'),
          ('d n', '6bc8a9409b1ef38e4f6f12121368883e'), ('d py', 'bda7be10061b0dcaeb66bebe3d736fee'),
          ('d pr', '6e2fce2239e2669b23f290049f87fbc4'), ('d s', '4c97144c5337f8a40c4fd494e622bf0d')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_delete_records(key, next_hash, sb_resources, hf):
    tester = Sb.delete_records(sb_resources.get_one(key), 'α2')
    assert hf.buddy2hash(tester) == next_hash


def test_delete_records2(sb_resources, hf):
    tester = Sb.delete_records(sb_resources.get_one("d f"), ['α1', 'α2'])
    assert hf.buddy2hash(tester) == "eca4f181dae3d7998464ff71e277128f"

    tester = Sb.delete_records(sb_resources.get_one("d f"), 'α1|α2')
    assert hf.buddy2hash(tester) == "eca4f181dae3d7998464ff71e277128f"

    with pytest.raises(ValueError) as e:
        Sb.delete_records(sb_resources.get_one("d f"), dict)
    assert "'patterns' must be a list or a string." in str(e.value)

    # Doesn't find anything to delete
    tester = Sb.delete_records(sb_resources.get_one("p g"), 'ML2', description=False)
    assert len(tester.records) == len(sb_resources.get_one("p g").records)

    tester = Sb.delete_records(sb_resources.get_one("p g"), 'ML2', description=True)
    assert hf.buddy2hash(tester) == "64049b9afd347f4507e264847e5f0500"


# ######################  '-drf', '--delete_recs_with_feature' ###################### #
hashes = [('p g', 'c486218295ec6d6d1a9c47023d952d40'), ('d g', 'fc91bfaed2df6926983144637cf0ba0f')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_delete_recs_with_feature(key, next_hash, sb_resources, hf):
    tester = Sb.delete_recs_with_feature(sb_resources.get_one(key), 'splice_.')
    assert hf.buddy2hash(tester) == next_hash, print(tester)


# #####################  '-drp', '--delete_repeats' ###################### ##
def test_delete_repeats(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources["duplicate"])
    tester = Sb.delete_repeats(tester)
    tester = Sb.find_repeats(tester)
    assert len(tester.repeat_ids) == 0
    assert len(tester.repeat_seqs) == 0


# ######################  '-ds', '--delete_small' ###################### #
def test_delete_small(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    tester = Sb.delete_small(tester, 1285)
    assert hf.buddy2hash(tester) == '196adf08d4993c51050289e5167dacdf'


# ######################  '-dt', '--delete_taxa' ###################### #
def test_delete_taxa(sb_resources, hf):
    tester = Sb.delete_taxa(sb_resources.get_one("p g"), "Lobata")
    assert hf.buddy2hash(tester) == "129c253374dd6171620884c92bece557"

    tester = Sb.delete_taxa(sb_resources.get_one("p g"), ["Lobata"])
    assert hf.buddy2hash(tester) == "129c253374dd6171620884c92bece557"

    tester = Sb.delete_taxa(sb_resources.get_one("p g"), ["leidyi", "Homo"])
    assert hf.buddy2hash(tester) == "96d74ce4bba524b4847fb2363f51e112"


# ######################  '-d2r', '--transcribe' and 'r2d', '--back_transcribe' ###################### #
hashes = [('d f', 'd2db9b02485e80323c487c1dd6f1425b', 'b831e901d8b6b1ba52bad797bad92d14'),
          ('d g', '360ad6806711e37a0a8aa5208536656b', '2e02a8e079267bd9add3c39f759b252c'),
          ('d n', 'f3bd73151645359af5db50d2bdb6a33d', 'cb1169c2dd357771a97a02ae2160935d'),
          ('d py', 'e55bd18b6d82a7fc3150338173e57e6a', '503e23720beea201f8fadf5dabda75e4'),
          ('d pr', 'a083e03b4e0242fa3c23afa80424d670', '52c23bd793c9761b7c0f897d3d757c12'),
          ('d s', '45b511f34653e3b984e412182edee3ca', '228e36a30e8433e4ee2cd78c3290fa6b')]


@pytest.mark.parametrize("key,d2r_hash,r2d_hash", hashes)
def test_transcribe(key, d2r_hash, r2d_hash, sb_resources, hf):
    tester = Sb.dna2rna(sb_resources.get_one(key))
    assert hf.buddy2hash(tester) == d2r_hash
    tester = Sb.rna2dna(tester)
    assert hf.buddy2hash(tester) == r2d_hash


def test_transcribe_pep_exception(sb_resources):  # Asserts that a ValueError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.dna2rna(sb_resources.get_one("p f"))


def test_back_transcribe_pep_exception(sb_resources):  # Asserts that a TypeError will be thrown if user inputs protein
    with pytest.raises(TypeError):
        Sb.rna2dna(sb_resources.get_one("p f"))


# ######################  '-efs', '--extract_feature_sequences' ###################### #
hashes = [('CDS', '956b6a14e02c9c2a2faa11ffb7e2bbed'), (["TMD"], 'd23b3ecdd5d432518c20572e7af03dc1'),
          (["TMD", "splice_a"], '344ffeb8e86442e0ae7e38d5b49072e1'),
          (["TMD2:TMD3"], 'fb54774a4a7d35dfe43e4ae31de0f44b'),
          (["TMD3:TMD2"], 'fb54774a4a7d35dfe43e4ae31de0f44b'), (["TMD2:foo"], '3cdbd5c8790f12871f8e04e40e315c93'),
          ("foo", '3cdbd5c8790f12871f8e04e40e315c93'), ([], '3cdbd5c8790f12871f8e04e40e315c93')]


@pytest.mark.parametrize("args,next_hash", hashes)
def test_extract_feature_sequences(args, next_hash, sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.extract_feature_sequences(tester, args)
    assert hf.buddy2hash(tester) == next_hash


# ######################  '-er', '--extract_regions' ###################### #
hashes = [('d f', '8c2fac57aedf6b0dab3d0f5bcf88e99f'), ('d g', '4211d7ea855794a657f6c3d73c67cd5a'),
          ('d n', '4063ab66ced2fafb080ceba88965d2bb'), ('d py', '33e6347792aead3c454bac0e05a292c6'),
          ('d pr', '9a5c491aa293c6cedd48c4c249d55aff'), ('d s', 'cd8d857feba9b6e459b8a9d56f11b7f5'),
          ('d q', "c312424429d3204f2e2c4cd42907d3d5"),
          ('p f', '2586d1e3fc283e6f5876251c1c57efce'), ('p g', 'cf5da50e9fdb3690d8d2732492c187e9'),
          ('p n', '6a27222d8f60ee8496cbe0c41648a116'), ('p py', 'c9a1dd913190f95bba5eca6a89685c75'),
          ('p pr', '6f579144a43dace285356ce6eb326d3b'), ('p s', '727099e0abb89482760eeb20f7edd0cd')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_extract_regions_multiformat(key, next_hash, sb_resources, hf):
    tester = Sb.extract_regions(sb_resources.get_one(key), "50:300")
    assert hf.buddy2hash(tester) == next_hash

    tester = Sb.extract_regions(sb_resources.get_one(key), "300:50")
    assert hf.buddy2hash(tester) == next_hash


hashes = [('0', 'a67e8747d5c7fb9e4077b5d9675009b8'), ('1', 'a67e8747d5c7fb9e4077b5d9675009b8'),
          ('-10000000', 'a67e8747d5c7fb9e4077b5d9675009b8'), (',1/', 'a67e8747d5c7fb9e4077b5d9675009b8'),
          ('1000000', '22822e8d5c7b9d087b3c80303dd3bcf0'), ('2,5,9,-5', '989c67ea5e2c1036b36b546004076109')]


@pytest.mark.parametrize("args,next_hash", hashes)
def test_extract_regions_singlets(args, next_hash, sb_resources, hf):
    tester = Sb.extract_regions(sb_resources.get_one("p g"), args)
    assert hf.buddy2hash(tester) == next_hash


hashes = [('0:10', '5873fb4b611edf38f652c756c3861c05'), ('1:10', '5873fb4b611edf38f652c756c3861c05'),
          ('10:1', '5873fb4b611edf38f652c756c3861c05'), (':10', '5873fb4b611edf38f652c756c3861c05'),
          ('-10:', '7869fea073804e8d023201024783bd2b'), ('40:75,89:100,432:-45', '487da42bfd6484620f3574890234b849')]


@pytest.mark.parametrize("args,next_hash", hashes)
def test_extract_regions_ranges(args, next_hash, sb_resources, hf):
    tester = Sb.extract_regions(sb_resources.get_one("p g"), args)
    assert hf.buddy2hash(tester) == next_hash


hashes = [('1/50', '40a5e3b46fbd4467a2ec0deab675292f'), ('-1/50', '1bae9ec3fe61eba3f47343111228c124'),
          ('1/-50', 'fc63bd3abb0a062ba8c2145b9391692b'), ('50/1', '0a8462e72f64fcd22544bb153b51b2b6'),
          ('50/25', 'a829efaf6997f6fa099930968ed543fa'), ('1:5/50', '7abf1717abc47c93dcee30b724ce0714'),
          ('-5:/50', '25d9ae5a167d9b1d828f710c2d4a5257'), (':5/50', '7abf1717abc47c93dcee30b724ce0714'),
          ('1:10,1/50,-1', '371f31b2d61d515f09aafb9f9105652c')]


@pytest.mark.parametrize("args,next_hash", hashes)
def test_extract_regions_mth_of_nth(args, next_hash, sb_resources, hf):
    tester = Sb.extract_regions(sb_resources.get_one("p g"), args)
    assert hf.buddy2hash(tester) == next_hash


hashes = [('2,5,9,-5', '5fded3dfa757a5ed90de0461c5095d7b'), ('10:45,60:75,90:-5', 'e014aa8b4d8cf9447d4a45cc4a1ada9c'),
          ('1:3,1/20,-1', '4258dfc66a07e849ac9c396aa2763c71')]


@pytest.mark.parametrize("args,next_hash", hashes)
def test_extract_regions_fastq(args, next_hash, sb_resources, hf):
    tester = Sb.extract_regions(sb_resources.get_one("d q"), args)
    assert hf.buddy2hash(tester) == next_hash


def test_extract_regions_edges(sb_resources):
    with pytest.raises(ValueError) as err:
        Sb.extract_regions(sb_resources.get_one("p g"), "foo")
    assert "Unable to decode the positions string 'foo'" in str(err)


# #####################  '-fcpg', '--find_CpG' ###################### ##
def test_find_cpg(sb_resources, hf):
    tester = sb_resources.get_one("d g")
    tester = Sb.find_cpg(tester)
    assert hf.buddy2hash(tester) == "9499f524da0c35a60502031e94864928"


# #####################  '-orf', '--find_orf' ###################### ##
def test_find_orf(sb_resources, hf):
    tester = Sb.SeqBuddy("ATGAAATTTCCCGGGTAG", in_format='raw', out_format='gb')
    tester = Sb.find_orfs(tester)
    assert hf.buddy2hash(tester) == "b34b2f3197b81ce6838a70a9784b79c2"

    tester = Sb.find_orfs(sb_resources.get_one("d g"))
    assert hf.buddy2hash(tester) == "214851fb21af6e6710414f132ae379ad"

    tester.out_format = "fasta"
    assert hf.buddy2hash(tester) == "9f912b59d69b82f72295e9e325f6f87a"

    tester = Sb.find_orfs(sb_resources.get_one("d g"), include_feature=False)
    assert hf.buddy2hash(tester) == "908744b00d9f3392a64b4b18f0db9fee"

    tester = Sb.find_orfs(sb_resources.get_one("r f"))
    assert hf.buddy2hash(tester) == "d2db9b02485e80323c487c1dd6f1425b"

    tester.out_format = "gb"
    assert hf.buddy2hash(tester) == "91e5934e1c688a35efaa4c98b1650701"

    tester = Sb.find_orfs(sb_resources.get_one("r f"), include_feature=False)
    tester.out_format = "gb"
    assert hf.buddy2hash(tester) == "456da121b26cb567d363b39765ca0dce"

    tester = Sb.find_orfs(sb_resources.get_one("d g"), min_size=500)
    assert hf.buddy2hash(tester) == "4f8a1825e1a1e2e1f2e18b5ce887c1a8"

    tester = Sb.find_orfs(sb_resources.get_one("d g"), rev_comp=False)
    assert hf.buddy2hash(tester) == "4f8a1825e1a1e2e1f2e18b5ce887c1a8"

    with pytest.raises(ValueError) as err:
        Sb.find_orfs(tester, min_size=2)
    assert "Open reading frames cannot be smaller than 6 residues." in str(err)

    tester = sb_resources.get_one("p g")
    with pytest.raises(TypeError) as err:
        Sb.find_orfs(tester)
    assert "Nucleic acid sequence required, not protein." in str(err)


# #####################  '-fp', '--find_pattern' ###################### ##
def test_find_pattern(sb_resources, hf):
    tester = Sb.find_pattern(sb_resources.get_one("d g"), "ATGGT")
    assert hf.buddy2hash(tester) == "f2a07ced2523081ba54249ae30775a40"
    tester.out_format = "fasta"
    assert hf.buddy2hash(tester) == "9160702ded16a80db79886d9de1abdf2"


hashes = [("d g", ['ATg{2}T'], {}, '57d80590f2e6596d33d5af936a584617'),
          ("d g", ["ATg{2}T", "tga.{1,6}tg"], {}, 'a13217987f5dd23f6fab71eb733271ff'),
          ("d g", ["ATg{2}T", "tga.{1,6}tg"], {'include_feature': False}, '2e02a8e079267bd9add3c39f759b252c'),
          ("p g", ['[bz]{2}x{50,100}[bz]{2}'], {'ambig': True}, '0e8d9d1a4c9a20f140acbad901292246'),
          ("d g", ['ATGGN{6}'], {'ambig': True}, '22b29f5d3aa45d7a2c7c5f3fdff2e210'),
          ("r f", ['AUGGN{6}'], {'ambig': True}, 'b7abcb4334232e38dfbac9f46234501a')]


@pytest.mark.parametrize("key,args,kwargs,next_hash", hashes)
def test_find_pattern2(key, args, kwargs, next_hash, sb_resources, hf):
    tester = Sb.find_pattern(sb_resources.get_one(key), *args, **kwargs)
    assert hf.buddy2hash(tester) == next_hash


# #####################  '-frp', '--find_repeats' ###################### ##
def test_find_repeats(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources["duplicate"])
    tester.unique_seqs, tester.repeat_ids, tester.repeat_seqs = {}, {}, {}
    Sb.find_repeats(tester)
    assert 'Seq1' in tester.unique_seqs
    assert 'Seq12' in tester.repeat_ids

    for key in tester.repeat_seqs:
        assert 'Seq12' in tester.repeat_seqs[key] or 'Seq10A' in tester.repeat_seqs[key]


# ######################  '-frs', '--find_restriction_sites' ###################### #
def test_restriction_sites_no_args(sb_resources):
    # No arguments passed in = commercial REs and any number of cut sites
    tester = Sb.find_restriction_sites(sb_resources.get_one("d g"))
    assert """AciI            232..235
     AciI            317..320
     AciI            462..465
     AciI            511..514
     AciI            641..644
     AciI            758..761
     AciI            788..791
     AciI            1095..1098
     AciI            1112..1115
     AciI            1133..1136
     AciI            1205..1208
     AciI            1210..1213""" in str(tester)

    new_res_dict = {}
    for key, value in tester.restriction_sites[-1][1].items():
        new_res_dict[str(key)] = value
    assert new_res_dict["AciI"] == [233, 318, 463, 512, 642, 759, 789,
                                    1096, 1113, 1134, 1206, 1211]


def test_restriction_sites_all_emzymes(sb_resources):
    # All enzymes
    tester = Sb.find_restriction_sites(sb_resources.get_one("d g"), enzyme_group=["all"])
    assert """AciI            232..235
     AciI            317..320
     AciI            462..465
     AciI            511..514
     AciI            641..644
     AciI            758..761
     AciI            788..791
     AciI            1095..1098
     AciI            1112..1115
     AciI            1133..1136
     AciI            1205..1208
     AciI            1210..1213""" in str(tester)

    new_res_dict = {}
    for key, value in tester.restriction_sites[-1][1].items():
        new_res_dict[str(key)] = value
    assert new_res_dict["AciI"] == [233, 318, 463, 512, 642, 759, 789,
                                    1096, 1113, 1134, 1206, 1211]


def test_restriction_sites_circular(sb_resources, sb_odd_resources, hf):
    # circular
    tester = Sb.find_restriction_sites(sb_resources.get_one("d g"), topology="circular")
    assert """LpnPI           1227..1230
     LpnPI           66..69
     LpnPI           103..106
     LpnPI           146..149
     LpnPI           167..170
     LpnPI           223..226
     LpnPI           281..284
     LpnPI           308..311
     LpnPI           320..323
     LpnPI           397..400""" in str(tester)

    res_sites = [x for x in tester.restriction_sites if x[0] == "Mle-Panxα1"][0][1]
    new_res_dict = {}
    for key, value in res_sites.items():
        new_res_dict[str(key)] = value
    assert new_res_dict["LpnPI"] == [1333, 61, 67, 104, 202, 209, 236, 291, 304, 349, 446, 488, 506, 515, 589, 598,
                                     697, 810, 811, 1001, 1046, 1058, 1072, 1085, 1110, 1205, 1237, 1256, 1296, 1319]

    # circular using genbank annotation
    tester = Sb.find_restriction_sites(Sb.SeqBuddy(sb_odd_resources["circular"]),
                                       enzyme_group=["EcoRI", "KspI", "TasI"])
    assert hf.buddy2hash(tester) == "6f032158c44ea7e4c52ee8843e270dfe"
    assert hf.string2hash(str(tester.restriction_sites)) == "4594f7ffa3b8afc2c070c1f33cbe128c"


def test_restriction_sites_limit_cuts(capsys, sb_resources, hf):
    # Specify a few REs and limit the number of cuts
    tester = sb_resources.get_one("d g")
    tester = Sb.find_restriction_sites(tester, min_cuts=2, max_cuts=4,
                                       enzyme_group=["EcoRI", "KspI", "TasI", "Bme1390I", "FooBR"])
    out, err = capsys.readouterr()
    # The c42b3bf and 0d2e5fdb hashes are for BioPython 1.70
    assert hf.buddy2hash(tester) in ['c42b3bf0367557383000b897432fed2d', '04cd62ab44f1479616370d04800fd54a']
    assert hf.string2hash(str(tester.restriction_sites)) in ["0d2e5fdba6fed434495481397a91e56a",
                                                             "e16a3aabf4681e7a4d186e7c7685f545"]
    assert "Warning: FooBR not a known enzyme" in err

    # Rerun to ensure that enzymes are not listed again if they are already in the feature list
    tester = Sb.find_restriction_sites(tester, min_cuts=2, max_cuts=4,
                                       enzyme_group=["EcoRI", "KspI", "TasI", "Bme1390I", "FooBR"])
    assert hf.buddy2hash(tester) in ['c42b3bf0367557383000b897432fed2d', '04cd62ab44f1479616370d04800fd54a']
    assert hf.string2hash(str(tester.restriction_sites)) in ["0d2e5fdba6fed434495481397a91e56a",
                                                             "e16a3aabf4681e7a4d186e7c7685f545"]
    # RNA
    tester = Sb.find_restriction_sites(sb_resources.get_one("r g"), min_cuts=2, max_cuts=4,
                                       enzyme_group=["EcoRI", "KspI", "TasI"])
    assert hf.buddy2hash(tester) == 'f440f8f7cbe21aad026d8cc7f41f98b6'

    with pytest.raises(TypeError) as e:
        Sb.find_restriction_sites(sb_resources.get_one("p g"))
    assert str(e.value) == "Unable to identify restriction sites in protein sequences."

    with pytest.raises(ValueError) as e:
        Sb.find_restriction_sites(tester, min_cuts=4, max_cuts=2)
    assert str(e.value) == "min_cuts parameter has been set higher than max_cuts."

    # 2-cutters and non-cutters
    Sb.find_restriction_sites(tester, enzyme_group=["AjuI", "AlwFI"])
    out, err = capsys.readouterr()
    assert "Warning: Double-cutters not supported." in err
    assert "Warning: No-cutters not supported." in err


# ######################  '-hsi', '--hash_sequence_ids' ###################### #
def test_hash_seq_ids(sb_resources):
    tester = sb_resources.get_one("d f")
    Sb.hash_ids(tester)
    assert len(tester.records[0].id) == 10

    tester = sb_resources.get_one("d f")
    tester = Sb.hash_ids(tester, 25)
    assert len(tester.records[0].id) == 25
    assert len(tester.hash_map) == 13

    tester_copy = Sb.make_copy(tester)
    tester_copy.reverse_hashmap()
    tester_copy = Sb.hash_ids(tester_copy, 25)
    assert tester_copy.records[0].id == tester.records[0].id

    tester_copy.reverse_hashmap()
    tester_copy.records[0].id = tester_copy.records[0].id[:-1]
    tester_copy = Sb.hash_ids(tester_copy, 25)
    assert tester_copy.records[0].id != tester.records[0].id


def test_hash_seq_ids_errors(sb_resources):
    tester = sb_resources.get_one("d f")
    with pytest.raises(TypeError) as e:
        Sb.hash_ids(tester, "foo")
    assert str(e.value) == "Hash length argument must be an integer, not <class 'str'>"

    with pytest.raises(ValueError) as e:
        Sb.hash_ids(tester, 0)
    assert str(e.value) == "Hash length must be greater than 0"

    tester.records *= 10
    with pytest.raises(ValueError) as e:
        Sb.hash_ids(tester, 1)
    assert "Insufficient number of hashes available to cover all sequences." in str(e.value)


# ##################### '-is', 'insert_seq' ###################### ##
def test_insert_seqs_start(sb_resources, hf):
    insert = 'AACAGGTCGAGCA'
    tester = sb_resources.get_one("d f")
    assert hf.buddy2hash(Sb.insert_sequence(tester, insert)) == 'f65fee08b892af5ef93caa1bf3cb3980'

    tester = sb_resources.get_one("d f")
    assert hf.buddy2hash(Sb.insert_sequence(tester, insert, -9000)) == 'f65fee08b892af5ef93caa1bf3cb3980'

    tester = sb_resources.get_one("d f")
    assert hf.buddy2hash(Sb.insert_sequence(tester, insert, -1)) == '792397e2e32e95b56ddc15b8b2310ec0'

    tester = sb_resources.get_one("d f")
    assert hf.buddy2hash(Sb.insert_sequence(tester, insert, 9000)) == '792397e2e32e95b56ddc15b8b2310ec0'

    tester = sb_resources.get_one("d f")
    Sb.insert_sequence(tester, insert, 100, ["α[23]", "α5"])
    assert hf.buddy2hash(tester) == 'edcd7934eb026ac3ea4b603ac85ca79f'

    tester = sb_resources.get_one("d f")
    assert hf.buddy2hash(Sb.insert_sequence(tester, insert, -25)) == '29cab1e72ba95572c3aec469270071e9'


# ######################  '-isd', '--in_silico_digest' ###################### #
def test_in_silico_digest(capsys, sb_resources, sb_odd_resources, hf):
    tester = Sb.in_silico_digest(sb_resources.get_one("d g"), enzyme_group=["NheI", "XhoI", "TseI", "FooBR"],
                                 topology="linear")
    out, err = capsys.readouterr()
    tester.records = [rec for rec in tester.records if rec.name == "Mle-Panxα10A"]
    assert hf.buddy2hash(tester) == "3466f95a9c7959c3697a60d8992f0277"
    assert "Warning: FooBR not a known enzyme" in err

    with pytest.raises(TypeError) as e:
        Sb.in_silico_digest(sb_resources.get_one("p g"))
    assert str(e.value) == "Unable to identify restriction sites in protein sequences."

    with pytest.raises(ValueError) as e:
        Sb.in_silico_digest(sb_resources.get_one("d g"), enzyme_group=["NheI", "XhoI", "TseI"], topology="cirt")
    assert str(e.value) == "Invalid topology. Accepted values are None, 'circular' and 'linear' "

    # 2-cutters and non-cutters
    tester = sb_resources.get_one("d g")
    Sb.in_silico_digest(tester, enzyme_group=["AjuI", "AlwFI"])
    out, err = capsys.readouterr()
    assert "Warning: Double-cutters not supported." in err
    assert "Warning: No-cutters not supported." in err

    # Circular cuts
    tester = Sb.SeqBuddy(sb_odd_resources["circular_digest"])
    tester = Sb.in_silico_digest(tester, enzyme_group=["EcoRI", "HpaI"], topology="circular")
    assert hf.buddy2hash(tester) == '7b4a311446f845cb7d0b401fec908f03'

    tester = Sb.SeqBuddy(sb_odd_resources["circular_digest"])
    tester = Sb.in_silico_digest(tester, enzyme_group=["EcoRI", "HpaI"])
    assert hf.buddy2hash(tester) == 'e7b43b83d52c2a3a1ab53ef35c7918a4'


# ######################  '-ip', '--isoelectric_point' ###################### #
def test_isoelectric_point(sb_resources, hf):
    for tester in sb_resources.get_list("p f g n pr s"):
        tester = Sb.isoelectric_point(Sb.clean_seq(tester))
        assert tester.to_dict()["Mle-Panxα12"].features[-1].qualifiers["value"] == 6.0117797852
        if tester.out_format == "gb":
            assert hf.buddy2hash(tester) == "39342fedfb3e43d6dea3e455e7e8bbb6", print(tester)

        with pytest.raises(TypeError):
            Sb.isoelectric_point(sb_resources.get_one("d f"))


# ######################  '-kt', '--keep_taxa' ###################### #
def test_keep_taxa(sb_resources):
    tester = Sb.keep_taxa(sb_resources.get_one("p g"), "Lobata")
    assert len(tester) == 2

    tester = Sb.keep_taxa(sb_resources.get_one("p g"), ["leidyi"])
    assert len(tester) == 3

    tester = Sb.keep_taxa(sb_resources.get_one("p g"), ["Homo"])
    assert len(tester) == 0

    tester = Sb.keep_taxa(sb_resources.get_one("p g"), ["leidyi", "Homo"])
    assert len(tester) == 3

    tester = Sb.keep_taxa(sb_resources.get_one("p g"), ["Lobata", "Homo"], match_all=True)
    assert len(tester) == 0

    tester = Sb.keep_taxa(sb_resources.get_one("p g"), ["Lobata", "leidyi"], match_all=True)
    assert len(tester) == 2


# ######################  '-lc', '--lowercase' and 'uc', '--uppercase'  ###################### #
hashes = [('d f', '25073539df4a982b7f99c72dd280bb8f', 'b831e901d8b6b1ba52bad797bad92d14'),
          ('d g', '2e02a8e079267bd9add3c39f759b252c', '2e02a8e079267bd9add3c39f759b252c'),
          ('d n', '52e74a09c305d031fc5263d1751e265d', 'cb1169c2dd357771a97a02ae2160935d'),
          ('d py', 'cfe6cb9c80aebd353cf0378a6d284239', '503e23720beea201f8fadf5dabda75e4'),
          ('d pr', '6e5542f41d17ff33afb530b4d07408a3', '52c23bd793c9761b7c0f897d3d757c12'),
          ('d s', 'b82538a4630810c004dc8a4c2d5165ce', '228e36a30e8433e4ee2cd78c3290fa6b'),
          ('p f', 'c10d136c93f41db280933d5b3468f187', '14227e77440e75dd3fbec477f6fd8bdc'),
          ('p g', '0a8462e72f64fcd22544bb153b51b2b6', '0a8462e72f64fcd22544bb153b51b2b6'),
          ('p n', '8b6737fe33058121fd99d2deee2f9a76', '17ff1b919cac899c5f918ce8d71904f6'),
          ('p py', '968ed9fa772e65750f201000d7da670f', 'aacda2f5d4077f23926400f74afa2f46'),
          ('p pr', 'ce423d5b99d5917fbef6f3b47df40513', 'e3dc2e0347f40fffec45d053f4f34c96'),
          ('p s', 'f35cbc6e929c51481e4ec31e95671638', 'c0dce60745515b31a27de1f919083fe9')]


@pytest.mark.parametrize("key,uc_hash,lc_hash", hashes)
def test_cases(key, uc_hash, lc_hash, sb_resources, hf):  # NOTE: Biopython always writes genbank to lower case
    tester = Sb.uppercase(sb_resources.get_one(key))
    assert hf.buddy2hash(tester) == uc_hash
    tester = Sb.lowercase(tester)
    assert hf.buddy2hash(tester) == lc_hash


# ######################  '-mui', '--make_ids_unique' ###################### #
def test_make_ids_unique(sb_odd_resources, hf):
    tester = Sb.SeqBuddy(sb_odd_resources["duplicate"])
    Sb.make_ids_unique(tester)
    assert hf.buddy2hash(tester) == "363c7ed14be59bcacede092b8f334a52"

    tester = Sb.SeqBuddy(sb_odd_resources["duplicate"])
    Sb.make_ids_unique(tester, sep="-", padding=4)
    assert hf.buddy2hash(tester) == "0054df3003ba16287159147f3b85dc7b"


# ######################  '-fn2p', '--map_features_nucl2prot' ###################### #
# Map the genbank DNA file to all protein files, and the fasta DNA file to fasta protein
hashes = [('p f', 'a31e54081bf7cf594a1a48ddb298d748'), ('p g', '98386831f03b40dbcbe52a6f5685a475'),
          ('p n', 'dc9ba63eabe524e7721c72523e288dab'), ('p py', 'fe0e5174ccfd87f97e90dcff6e7a94e6'),
          ('p pr', '0c40a1e59634b67f0ef13515024b72fb'), ('p s', 'a128d1b61939d2aee7fed1f1225d19e9')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_map_features_nucl2prot(key, next_hash, sb_resources, hf):
    tester = Sb.map_features_nucl2prot(sb_resources.get_one("d g"), sb_resources.get_one(key))
    tester.out_format = "gb"
    assert hf.buddy2hash(tester) == next_hash

    if key == "p f":
        tester = Sb.map_features_nucl2prot(sb_resources.get_one("d f"), sb_resources.get_one(key))
        tester.out_format = "gb"
        assert hf.buddy2hash(tester) == "dffab18027b2c445e442b423d9e999f0"


def test_map_features_nucl2prot_2(capsys, sb_resources, hf):
    tester = sb_resources.get_one("d g")
    Sb.pull_record_ends(tester, 1300)
    tester = Sb.annotate(tester, "foo", [(20, 40), (50, 60)], pattern="α9")
    mapped = Sb.map_features_nucl2prot(Sb.make_copy(tester), sb_resources.get_one("p f"))
    mapped.out_format = "gb"
    assert hf.buddy2hash(mapped) == "35c20722fa38101356a3c11accc02691"
    out, err = capsys.readouterr()
    assert hf.string2hash(err) == "6c840e4acaaf4328672ca164f854000a"

    prot_tester = sb_resources.get_one("p g")
    prot_tester.records = sorted(prot_tester.records, key=lambda x: x.id)
    dna_tester = sb_resources.get_one("d g")
    dna_tester.records = sorted(dna_tester.records, key=lambda x: x.id)
    Sb.rename(dna_tester, "α4", "A4")

    Sb.map_features_nucl2prot(Sb.make_copy(dna_tester), prot_tester, mode="list")
    assert hf.buddy2hash(prot_tester) == "f7beb0fd652a9d1910ead5bf79120173"

    Sb.map_features_nucl2prot(Sb.make_copy(dna_tester), prot_tester, mode="key")
    assert hf.buddy2hash(prot_tester) == "f7beb0fd652a9d1910ead5bf79120173"
    out, err = capsys.readouterr()
    assert hf.string2hash(err) == "6fd2b5f2a7a3995d3f49c4919c3358b0"

    mock_obj = type('test', (object,), {})()
    mock_obj.start = 2
    mock_obj.end = 6
    dna_tester.records[0].features[0].location = mock_obj
    with pytest.raises(TypeError) as e:
        Sb.map_features_nucl2prot(Sb.make_copy(dna_tester), sb_resources.get_one("p f"))
    assert "FeatureLocation or CompoundLocation object required." in str(e.value)

    tester.records[0].features[0].location = tester.records[1].features[0].location
    with pytest.raises(ValueError) as e:
        Sb.map_features_nucl2prot(Sb.make_copy(tester), sb_resources.get_one("p f"), mode="foo")
    assert "'mode' must be either 'key' or 'position'" in str(e.value)

    with pytest.raises(ValueError) as e:
        Sb.pull_recs(tester, "α[1-8]")
        Sb.map_features_nucl2prot(Sb.make_copy(tester), sb_resources.get_one("p f"), mode="list")
    assert "The two input files do not contain the same number of sequences" in str(e.value)


# ######################  '-fp2n', '--map_features_prot2nucl' ###################### #
hashes = [('d f', '47a7b6cf12399a3c58995d53b334a0c4'), ('d g', 'feceaf5e17935afb100b4b6030e27fee'),
          ('d n', '7bdacec654710d4280c8c1d3239b1ada'), ('d py', '91e29792e174e97a9c2d49b296bf33a2'),
          ('d pr', 'cba3651d35e33ad363f0cef1c1596f96'), ('d s', '9c88f299ebf8f4bf50f0f9e2655394fe')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_map_features_prot2nucl(key, next_hash, sb_resources, hf):
    tester = Sb.map_features_prot2nucl(sb_resources.get_one("p g"), sb_resources.get_one(key))
    tester.out_format = "gb"
    assert hf.buddy2hash(tester) == next_hash


def test_map_features_prot2nucl_2(capsys, sb_resources, hf):
    prot_tester = sb_resources.get_one("p g")
    Sb.pull_record_ends(prot_tester, 450)
    prot_tester = Sb.annotate(prot_tester, "foo", [(20, 40), (50, 60)], pattern="α9")
    mapped = Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), sb_resources.get_one("d f"))
    mapped.out_format = "gb"
    assert hf.buddy2hash(mapped) == "12bfe195548fde9e539a3426c9f0dc40", print(mapped)
    out, err = capsys.readouterr()
    assert err == "Warning: size mismatch between aa and nucl seqs for Mle-Panxα7A --> 450, 1875\n"

    prot_tester = sb_resources.get_one("p g")
    prot_tester.records = sorted(prot_tester.records, key=lambda x: x.id)
    dna_tester = sb_resources.get_one("d g")
    dna_tester.records = sorted(dna_tester.records, key=lambda x: x.id)
    Sb.rename(prot_tester, "α4", "A4")
    Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), dna_tester, mode="list")
    assert hf.buddy2hash(dna_tester) == "c32f43cc5205867c0eb1d3873e27319b"

    Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), dna_tester, mode="key")
    assert hf.buddy2hash(dna_tester) == "c32f43cc5205867c0eb1d3873e27319b"
    out, err = capsys.readouterr()
    assert hf.string2hash(err) == "c1f51587b07d2cc6156b8aac07384834"

    mock_obj = type('test', (object,), {})()
    mock_obj.start = 2
    mock_obj.end = 6
    prot_tester.records[0].features[0].location = mock_obj
    with pytest.raises(TypeError) as e:
        Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), sb_resources.get_one("d f"))
    assert "FeatureLocation or CompoundLocation object required." in str(e.value)

    prot_tester.records[0].features[0].location = prot_tester.records[1].features[0].location
    with pytest.raises(ValueError) as e:
        Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), sb_resources.get_one("p f"), mode="foo")
    assert "'mode' must be either 'key' or 'position'" in str(e.value)

    with pytest.raises(ValueError) as e:
        Sb.pull_recs(prot_tester, "α[1-8]")
        Sb.map_features_prot2nucl(Sb.make_copy(prot_tester), sb_resources.get_one("p f"), mode="list")
    assert "The two input files do not contain the same number of sequences" in str(e.value)


# #####################  make_groups' ###################### ##
def test_make_groups(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources["cnidaria_pep"])
    sb_list = Sb.make_groups(tester)
    assert len(sb_list) == 20
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier == seqbuddy.records[0].id

    sb_list = Sb.make_groups(tester, split_patterns=["u", "h"])
    assert len(sb_list) == 4
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Unknown", "Hv", "C", "Pp"]

    sb_list = Sb.make_groups(tester, num_chars=1)
    assert len(sb_list) == 5
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["A", "H", "C", "P", "N"]

    sb_list = Sb.make_groups(tester, split_patterns=["l"], num_chars=3)
    assert len(sb_list) == 3
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Unknown", "Ae", "C"]

    sb_list = Sb.make_groups(tester, regex="([ACH]).*([βγ])")
    assert len(sb_list) == 5
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Unknown", "Aβ", "Cβ", "Hβ", "Cγ"]

    sb_list = Sb.make_groups(tester, regex="Ate")
    assert len(sb_list) == 2
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["Ate", "Unknown"]

    sb_list = Sb.make_groups(tester, regex="Panx.(G*)")
    assert len(sb_list) == 2
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        assert seqbuddy.identifier in ["G", "Unknown"]


# #####################  '-max', '--max_recs' ###################### ##
def test_max_recs(sb_resources):
    tester = Sb.max_records(sb_resources.get_one("p f"))
    assert tester.records[0].id == "Mle-Panxα7A"

    tester = Sb.max_records(sb_resources.get_one("p f"), 3)
    assert len(tester.records) == 3


# #####################  '-mg', '--merge' ###################### ##
def test_merge(sb_resources, sb_odd_resources, hf):
    tester = Sb.SeqBuddy(sb_odd_resources['dummy_feats'])
    tester = Sb.merge(tester, sb_resources.get_one("d g"))
    assert hf.buddy2hash(tester) == "bae5aeb130b3d5319378a122a6f61df5"

    tester.records[0].seq = tester.records[1].seq
    with pytest.raises(RuntimeError) as e:
        Sb.merge(tester, sb_resources.get_one("d g"))
    assert "Sequence mismatch for record 'Mle-Panxα1'" in str(e.value)


# #####################  '-min', '--min_recs' ###################### ##
def test_min_recs(sb_resources):
    tester = Sb.min_records(sb_resources.get_one("p f"))
    assert tester.records[0].id == "Mle-Panxα10B"

    tester = Sb.min_records(sb_resources.get_one("p f"), 3)
    assert len(tester.records) == 3


# ######################  '-mw', '--molecular_weight' ###################### #
def test_molecular_weight(sb_resources, sb_odd_resources, hf):
    # Unambiguous DNA
    tester = Sb.molecular_weight(sb_resources.get_one("d g"))
    assert tester.molecular_weights['masses_ds'][0] == 743477.1
    assert tester.molecular_weights['masses_ss'][0] == 371242.6
    assert hf.buddy2hash(tester) == "08c8ab50e6b66adb8e579df3c923c2bc", tester.write("temp.del")
    # Ambiguous DNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(sb_odd_resources['ambiguous_dna']))
    assert tester.molecular_weights['masses_ds'][0] == 743477.08
    assert tester.molecular_weights['masses_ss'][0] == 371202.59
    # Unambiguous RNA
    tester = Sb.molecular_weight(sb_resources.get_one("r f"))
    assert tester.molecular_weights['masses_ss'][0] == 387372.6
    # Ambiguous RNA
    tester = Sb.molecular_weight(Sb.SeqBuddy(sb_odd_resources['ambiguous_rna']))
    assert tester.molecular_weights['masses_ss'][0] == 387371.6
    # Protein
    tester = Sb.molecular_weight(sb_resources.get_one("p g"))
    assert tester.molecular_weights['masses_ss'][0] == 45692.99
    assert hf.buddy2hash(tester) == "6e3b7b4fda2e126d4c0214133f34209a"

    tester = sb_resources.get_one("d f")
    seq = str(tester.records[0].seq)
    seq = "j" + seq
    tester.records[0].seq = Seq(seq, alphabet=tester.records[0].seq.alphabet)
    with pytest.raises(KeyError) as err:
        Sb.molecular_weight(tester)
    assert "Invalid residue \'J\' in record Mle-Panxα9. \'J\' is not valid a valid character in IUPACAmbiguousDNA()." \
           in str(err)


# ######################  '-ns', '--num_seqs' ###################### #
def test_num_seqs(sb_resources):
    for tester in sb_resources.get_list("d p f g pr n s"):
        assert Sb.num_seqs(tester) == 13
    for tester in sb_resources.get_list("d p py ps"):
        assert Sb.num_seqs(tester) == 8


def test_empty_file(sb_odd_resources):
    tester = Sb.SeqBuddy(sb_odd_resources["blank"])
    assert type(tester) == Sb.SeqBuddy
    assert len(tester.records) == 0


# ######################  '-ofa', '--order_features_alphabetically' ###################### #
hashes = [('d f', 'b831e901d8b6b1ba52bad797bad92d14', 'b831e901d8b6b1ba52bad797bad92d14'),
          ('d g', '21547b4b35e49fa37e5c5b858808befb', '3b718ec3cb794bcb658d900e517110cc'),
          ('d n', 'cb1169c2dd357771a97a02ae2160935d', 'cb1169c2dd357771a97a02ae2160935d'),
          ('d py', '503e23720beea201f8fadf5dabda75e4', '503e23720beea201f8fadf5dabda75e4'),
          ('d pr', '52c23bd793c9761b7c0f897d3d757c12', '52c23bd793c9761b7c0f897d3d757c12'),
          ('d s', '228e36a30e8433e4ee2cd78c3290fa6b', '228e36a30e8433e4ee2cd78c3290fa6b'),
          ('p f', '14227e77440e75dd3fbec477f6fd8bdc', '14227e77440e75dd3fbec477f6fd8bdc'),
          ('p g', 'ae394ce8a346e2294a6981015a372916', 'd96a10c51bbfabb935164ca834997554'),
          ('p n', '17ff1b919cac899c5f918ce8d71904f6', '17ff1b919cac899c5f918ce8d71904f6'),
          ('p py', '968ed9fa772e65750f201000d7da670f', '968ed9fa772e65750f201000d7da670f'),
          ('p pr', 'ce423d5b99d5917fbef6f3b47df40513', 'ce423d5b99d5917fbef6f3b47df40513'),
          ('p s', 'c0dce60745515b31a27de1f919083fe9', 'c0dce60745515b31a27de1f919083fe9')]


@pytest.mark.parametrize("key,fwd_hash,rev_hash", hashes)
def test_order_features_alphabetically(key, fwd_hash, rev_hash, sb_resources, hf):
    tester = Sb.order_features_alphabetically(sb_resources.get_one(key))
    assert hf.buddy2hash(tester) == fwd_hash
    tester = Sb.order_features_alphabetically(sb_resources.get_one(key), reverse=True)
    assert hf.buddy2hash(tester) == rev_hash


# ######################  '-ofp', '--order_features_by_position' ###################### #
hashes = [('d f', 'b831e901d8b6b1ba52bad797bad92d14', 'b831e901d8b6b1ba52bad797bad92d14'),
          ('d g', '2e02a8e079267bd9add3c39f759b252c', '4345a14fe27570b3c837c30a8cb55ea9'),
          ('d n', 'cb1169c2dd357771a97a02ae2160935d', 'cb1169c2dd357771a97a02ae2160935d'),
          ('d py', '503e23720beea201f8fadf5dabda75e4', '503e23720beea201f8fadf5dabda75e4'),
          ('d pr', '52c23bd793c9761b7c0f897d3d757c12', '52c23bd793c9761b7c0f897d3d757c12'),
          ('d s', '228e36a30e8433e4ee2cd78c3290fa6b', '228e36a30e8433e4ee2cd78c3290fa6b'),
          ('p f', '14227e77440e75dd3fbec477f6fd8bdc', '14227e77440e75dd3fbec477f6fd8bdc'),
          ('p g', '0a8462e72f64fcd22544bb153b51b2b6', '121ba21504e09c97cedb4f1fc4cf40bb'),
          ('p n', '17ff1b919cac899c5f918ce8d71904f6', '17ff1b919cac899c5f918ce8d71904f6'),
          ('p py', '968ed9fa772e65750f201000d7da670f', '968ed9fa772e65750f201000d7da670f'),
          ('p pr', 'ce423d5b99d5917fbef6f3b47df40513', 'ce423d5b99d5917fbef6f3b47df40513'),
          ('p s', 'c0dce60745515b31a27de1f919083fe9', 'c0dce60745515b31a27de1f919083fe9')]


@pytest.mark.parametrize("key,fwd_hash,rev_hash", hashes)
def test_order_features_position(key, fwd_hash, rev_hash, sb_resources, hf):
    tester = Sb.order_features_by_position(sb_resources.get_one(key))
    assert hf.buddy2hash(tester) == fwd_hash
    tester = Sb.order_features_by_position(sb_resources.get_one(key), reverse=True)
    assert hf.buddy2hash(tester) == rev_hash


# ######################  '-oi', '--order_ids' ###################### #
hashes = [('d f', 'e9090efdd362d527a115049dfced42cd', 'fcb016fff87d26822fa518d62e355c65'),
          ('d g', 'c0d656543aa5d20a266cffa790c035ce', '2507c667a304fdc003bc68255e094d7b'),
          ('d n', '132757da01b3caf174d024efdb2c3acd', '286bac7a213997924203622c3357457c'),
          ('d py', '3c49bdc1b0fe4e1d6bfc148eb0293e21', 'd6e79a5faeaff396aa7eab0b460c3eb9'),
          ('d pr', '684a99151e8cbf1e5eb96e44b875ba08', 'a5fcbd95e837b4807408124c2396db6e'),
          ('d s', 'c06985b566277a29e598dea0dc41baef', '20b43a724670a1151b1f0418478046ef')]


@pytest.mark.parametrize("key,fwd_hash,rev_hash", hashes)
def test_order_ids1(key, fwd_hash, rev_hash, sb_resources, hf):
    tester = Sb.order_ids(sb_resources.get_one(key))
    assert hf.buddy2hash(tester) == fwd_hash
    tester = Sb.order_ids(sb_resources.get_one(key), reverse=True)
    assert hf.buddy2hash(tester) == rev_hash


def test_order_ids2(sb_resources, hf):
    seqbuddy = sb_resources.get_one("p n")
    Sb.rename(seqbuddy, "Mle-Panxα4", "Mle004-Panxα4")
    Sb.rename(seqbuddy, "Mle-Panxα5", "Mle05-Panxα5")
    Sb.rename(seqbuddy, "Mle-Panxα9", "aMle-PanxαBlahh")
    Sb.order_ids(seqbuddy)
    assert hf.buddy2hash(seqbuddy) == "5c1316e18205432b044101e720646cd5"


# ######################  '-oir', '--order_ids_randomly' ###################### #
hashes = [('d f', '78fa4ce6cf7fa4e8e82f0a7ccee260dd'), ('d g', '220ae6ddfe74d46127b95f2715e28d0a'),
          ('d n', '83cff49333f9c3c46ee1f4cf4f5e963e'), ('p py', 'fd91fb622d9f5dad099c7a566dc7bd5b'),
          ('p pr', '82238767f3d3793e9504eab5b8c286dc'), ('p s', '0c9c4769daed5c765ea9e52a0454ee28')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_order_ids_randomly(key, next_hash, sb_resources, hf):
    tester = sb_resources.get_one(key)
    tester = Sb.order_ids_randomly(tester, r_seed=12345)
    assert hf.buddy2hash(tester) == next_hash


def test_order_ids_randomly2(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    tester = Sb.pull_recs(tester, "α[789]")
    tester = Sb.order_ids_randomly(tester, r_seed=12345)
    assert hf.buddy2hash(tester) == "e2e95a6f695b5f98746664db179a9aea"

    Sb.pull_recs(tester, "α[89]")
    Sb.order_ids_randomly(tester, r_seed=12345)
    assert hf.buddy2hash(tester) == "6fc8502070994a24ea61e752a5bc28cb"

    Sb.pull_recs(tester, "α9")
    Sb.order_ids_randomly(tester, r_seed=12345)
    assert hf.buddy2hash(tester) == "d1d16f79697d7502460ca0d9e65e51b8"

    tester = Sb.SeqBuddy(tester.records * 3, out_format="fasta")
    Sb.order_ids_randomly(tester, r_seed=12345)
    assert hf.buddy2hash(tester) == "bb75e7fc15f131e31271ea5006241615", print(tester)


# ######################  '-obl', '--order_recs_by_len' ###################### #
def test_order_recs_by_len(sb_resources, hf):
    tester = sb_resources.get_one("p f")
    assert hf.buddy2hash(Sb.order_recs_by_len(tester)) == "bb114c02bfda1d1ad90bfb3375dc3a3b"
    assert hf.buddy2hash(Sb.order_recs_by_len(tester, rev=True)) == "e99cf3d600d725e6dbd0cd5a3800face"


# #####################  '-ppo', '--prepend_organism' ###################### ##
def test_prepend_organism(sb_resources, hf):
    tester = sb_resources.get_one("p g")
    tester.records[4].annotations["organism"] = "Testus robustis"
    tester = Sb.prepend_organism(tester)
    tester.out_format = "fasta"
    assert hf.buddy2hash(tester) == "12af6bc1c299f3aa1034825ceacb51a3", print(tester)
    assert len(tester.prefix_map) == 2
    assert "Mlei" in tester.prefix_map
    assert "Trob" in tester.prefix_map

    tester = sb_resources.get_one("p g")
    tester.records[4].annotations["organism"] = "Testus"
    tester = Sb.prepend_organism(tester)
    tester.out_format = "fasta"
    assert hf.buddy2hash(tester) == "c0c13e6224893ba1700ae4811d667b44", print(tester)
    assert len(tester.prefix_map) == 2
    assert "Mlei" in tester.prefix_map
    assert "Test" in tester.prefix_map

    tester = Sb.prepend_organism(tester, 5)
    tester.out_format = "fasta"
    assert hf.buddy2hash(tester) == "54921aef13643851e4b097be235044b1", print(tester)
    assert len(tester.prefix_map) == 2
    assert "Mleid" in tester.prefix_map
    assert "Testu" in tester.prefix_map

    tester = sb_resources.get_one("p g")
    tester.records[4].annotations["organism"] = "Moby leily"
    tester = Sb.prepend_organism(tester)
    assert hf.buddy2hash(tester) == "2fd79883e44b2d3a514eddc9cbef4d54", print(tester)

    with pytest.raises(ValueError) as err:
        Sb.prepend_organism(tester, 0)

    assert "Prefix length must be > 2" in str(err)


# #####################  '-psc', '--prosite_scan' ###################### ##
def test_prosite_scan_init(sb_resources):
    seqbuddy = sb_resources.get_one("d f")
    ps_scan = Sb.PrositeScan(seqbuddy)
    assert hash(ps_scan.seqbuddy) == hash(seqbuddy)
    assert not ps_scan.quiet
    assert ps_scan.r_seed is None
    assert ps_scan.base_url == 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5'
    assert ps_scan.check_interval == 30
    assert len(ps_scan.http_headers) == 1
    assert 'User-Agent' in ps_scan.http_headers, print(ps_scan.http_headers)
    for key in ['data_dir', 'diagnostics', 'email', 'user_hash']:
        assert key in ps_scan.user_deets


def test_prosite_scan_rest_request(sb_resources, monkeypatch):
    def mock_urlopen(req, request_data=None):
        tmp_file = br.TempFile(byte_mode=True)
        file_text = "Hello world\n%s\n%s" % (req.full_url, request_data)
        tmp_file.write(file_text.encode())
        return tmp_file.get_handle("r")

    monkeypatch.setattr(urllib.request, "urlopen", mock_urlopen)
    seqbuddy = sb_resources.get_one("d f")
    ps_scan = Sb.PrositeScan(seqbuddy)
    assert ps_scan._rest_request("http://www.foo.bar", {"test1": 1}) == "Hello world\nhttp://www.foo.bar\n{'test1': 1}"
    assert ps_scan._rest_request("http://www.foo.bar") == "Hello world\nhttp://www.foo.bar\nNone"


def test_prosite_scan_mc_run_prosite(sb_resources, hf, monkeypatch, capsys):

    def status():
        for next_status in ["RUNNING", "PENDING", "FINISHED"]:
            yield next_status
    status_obj = status()

    def mock_rest_request(self, url, *args):
        if "run" in url:
            file_text = "job_1"
        elif "status" in url:
            file_text = next(status_obj)
        elif "result" in url:
            file_text = """\
XP_009950933.1\tdd1675f26044a927878cc3d92ea431e5\t565\tProSiteProfiles\tPS50835\tIg-like domain profile.\t\
21\t111\t6.542\tT\t25-01-2018\tIPR007110\tImmunoglobulin-like domain
XP_009950933.1\tdd1675f26044a927878cc3d92ea431e5\t565\tProSiteProfiles\tPS50835\tIg-like domain profile.\t\
134\t186\t8.157\tT\t25-01-2018\tIPR007110\tImmunoglobulin-like domain
XP_009950933.1\tdd1675f26044a927878cc3d92ea431e5\t565\tProSiteProfiles\tPS50104\tTIR domain profile.\t\
382\t403\t34.751\tT\t25-01-2018\tIPR000157\tToll/interleukin-1 receptor homology (TIR) domain
"""
        else:
            raise RuntimeError("This shouldn't ever happen", self, args)

        return file_text

    monkeypatch.setattr(Sb.PrositeScan, "_rest_request", mock_rest_request)
    monkeypatch.setattr(Sb.time, "sleep", lambda _: True)
    out_file = br.TempFile()
    seqbuddy = sb_resources.get_one("p f")
    seqbuddy.records = [seqbuddy.records[0]]
    ps_scan = Sb.PrositeScan(seqbuddy)
    ps_scan._mc_run_prosite(seqbuddy.records[0], [out_file.path, Sb.Lock()])
    output = out_file.read()
    assert hf.string2hash(output) == "8ad55e8a179dc4de67ed09a738cac814", print(output)

    def raise_attrib_error(_, url, *__):
        if "result" in url:
            raise AttributeError("Just raising something.")
        else:
            return "DONE"

    monkeypatch.setattr(Sb.PrositeScan, "_rest_request", raise_attrib_error)
    with pytest.raises(AttributeError) as err:
        ps_scan._mc_run_prosite(seqbuddy.records[0], [out_file.path, Sb.Lock()])

    assert "Just raising something." in str(err)

    def raise_keyboard_error(_, url, *__):
        if "result" in url:
            raise KeyboardInterrupt("Interrupted!")
        else:
            return "DONE"

    monkeypatch.setattr(Sb.PrositeScan, "_rest_request", raise_keyboard_error)
    ps_scan._mc_run_prosite(seqbuddy.records[0], [out_file.path, Sb.Lock()])
    out, err = capsys.readouterr()
    assert "Mle-Panxα12 killed with KeyboardInterrupt" in err, print(out, err)

    def raise_http_error(_, url, *__):
        if "result" in url:
            raise urllib.error.HTTPError("", 400, "blahh", {}, br.TempFile())
        else:
            return "DONE"

    monkeypatch.setattr(Sb.PrositeScan, "_rest_request", raise_http_error)
    ps_scan._mc_run_prosite(seqbuddy.records[0], [out_file.path, Sb.Lock()])
    out, err = capsys.readouterr()
    assert "Error: Failed to retrieve Mle-Panxα12" in err

    def raise_http_error(_, url, *__):
        if "result" in url:
            raise urllib.error.HTTPError("", 200, "Weird 200 error", {}, br.TempFile())
        else:
            return "DONE"

    monkeypatch.setattr(Sb.PrositeScan, "_rest_request", raise_http_error)
    with pytest.raises(urllib.error.HTTPError) as err:
        ps_scan._mc_run_prosite(seqbuddy.records[0], [out_file.path, Sb.Lock()])
    assert "Weird 200 error" in str(err), print(err)


def test_prosite_scan_run(sb_resources, hf, monkeypatch, capsys):
    def mock_mc_run_prosite(_, _rec, args):
        out_file_path, lock = args
        temp_seq = Sb.SeqBuddy([_rec], out_format="gb")
        feature1 = SeqFeature(FeatureLocation(0, 100), type="Region", qualifiers={"note": "Foo"})
        feature2 = SeqFeature(FeatureLocation(150, 200), type="Region", qualifiers={"note": "Bar"})
        temp_seq.records[0].features = [feature1, feature2]
        with lock:
            with open(out_file_path, "a") as out_file:
                out_file.write("%s\n" % str(temp_seq))

    monkeypatch.setattr(Sb.PrositeScan, "_mc_run_prosite", mock_mc_run_prosite)
    seqbuddy = sb_resources.get_one("d g")
    seqbuddy.records = [seqbuddy.records[0], seqbuddy.records[3]]

    ps_scan = Sb.PrositeScan(seqbuddy)
    seqbuddy = ps_scan.run()
    seqbuddy.records = sorted(seqbuddy.records, key=lambda rec: rec.id)
    assert hf.buddy2hash(seqbuddy) == "722f03e74786908a54c3d80de0653ed6", print(seqbuddy)

    seqbuddy = sb_resources.get_one("p g")
    seqbuddy.records = seqbuddy.records[:2]
    ps_scan = Sb.PrositeScan(seqbuddy)
    seqbuddy = ps_scan.run()
    seqbuddy.records = sorted(seqbuddy.records, key=lambda rec: rec.id)
    assert hf.buddy2hash(seqbuddy) == "848d2308deeaf4c57341d58018fd3905", print(seqbuddy)

    seqbuddy = sb_resources.get_one("p f")
    seqbuddy.records = [seqbuddy.records[0]]
    rec_seq = str(seqbuddy.records[0].seq)
    new_seq = rec_seq[:50] + "*" + rec_seq[50:150] + "*" + rec_seq[150:]
    seqbuddy.records[0].seq = Seq(new_seq, alphabet=seqbuddy.records[0].seq.alphabet)
    ps_scan = Sb.PrositeScan(seqbuddy)
    seqbuddy = ps_scan.run()
    seqbuddy.out_format = "gb"
    assert hf.buddy2hash(seqbuddy) == "1d420ba126c440698e85adb03f0d651d", print(seqbuddy)

    def raise_keyboard_error(*_, **__):
        raise KeyboardInterrupt("Interrupted!")
    monkeypatch.setattr(br, "run_multicore_function", raise_keyboard_error)

    capsys.readouterr()
    seqbuddy = sb_resources.get_one("p g")
    seqbuddy.records = sorted(seqbuddy.records[:2], key=lambda rec: rec.id)
    initial_hash = hf.buddy2hash(seqbuddy)
    ps_scan = Sb.PrositeScan(seqbuddy)
    seqbuddy = ps_scan.run()
    seqbuddy.records = sorted(seqbuddy.records, key=lambda rec: rec.id)
    assert hf.buddy2hash(seqbuddy) == initial_hash, print(seqbuddy)


# #####################  '-prr', '--pull_random_recs' ###################### ##
hashes = [('d f', '8a27843a57fc5fdcfc2e3552565a6c1d'), ('d g', 'd5e75f41571a5123769afc9814a571a7'),
          ('d n', '1f0e0124b2ae310284c773ea01a4f709'), ('p py', '5713b6020069c45e51146b4c16978ded'),
          ('p pr', 'f1e8bf5bba47a6bbeb229cda55e2bec1'), ('p s', '4e09f49bff7ea0a4903e0906e8c502ed')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_pull_random_recs(key, next_hash, sb_resources, hf):
    tester = sb_resources.get_one(key)
    tester = Sb.pull_random_recs(tester, count=3, r_seed=12345)
    assert hf.buddy2hash(tester) == next_hash


# #####################  '-pre', '--pull_record_ends' ###################### ##
def test_pull_record_ends(sb_resources, hf):
    tester = Sb.pull_record_ends(sb_resources.get_one("d g"), 10)
    assert hf.buddy2hash(tester) == 'd46867e4ca7a9f474c45473fc3495413'

    tester = Sb.pull_record_ends(sb_resources.get_one("d g"), 2000)
    assert hf.buddy2hash(tester) == '908744b00d9f3392a64b4b18f0db9fee'

    tester = Sb.pull_record_ends(sb_resources.get_one("d g"), -10)
    assert hf.buddy2hash(tester) == 'd7970570d65872993df8a3e1d80f9ff5'

    tester = Sb.pull_record_ends(sb_resources.get_one("d g"), -2000)
    assert hf.buddy2hash(tester) == '908744b00d9f3392a64b4b18f0db9fee'

    with pytest.raises(ValueError):
        Sb.pull_record_ends(sb_resources.get_one("d f"), 'foo')


# ######################  '-pr', '--pull_records' ###################### #
hashes = [('d f', '5b4154c2662b66d18776cdff5af89fc0'), ('d g', 'e196fdc5765ba2c47f97807bafb6768c'),
          ('d n', 'bc7dbc612bc8139eba58bf896b7eaf2f'), ('d py', '7bb4aac2bf50381ef1d27d82b7dd5a53'),
          ('d pr', '72328eddc751fd79406bb911dafa57a2'), ('d s', 'b006b40ff17ba739929448ae2f9133a6')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_pull_recs(key, next_hash, sb_resources, hf):
    tester = Sb.pull_recs(sb_resources.get_one(key), 'α2')
    assert hf.buddy2hash(tester) == next_hash


def test_pull_recs2(sb_resources, hf):
    tester = Sb.pull_recs(sb_resources.get_one("p g"), 'ML2', description=False)
    assert len(tester.records) == 0

    tester = Sb.pull_recs(sb_resources.get_one("p g"), 'ML2', description=True)
    assert hf.buddy2hash(tester) == "466acff4d79969ea30cfd94e1f996a27"


# ######################  '-prf', '--pull_records_with_feature' ###################### #
hashes = [('p g', '8c41bd906501628f987a055ec829c9b6'), ('d g', '36757409966ede91ab19deb56045d584')]


@pytest.mark.parametrize("key, next_hash", hashes)
def test_pull_records_with_feature(key, next_hash, sb_resources, hf):
    tester = Sb.pull_recs_with_feature(sb_resources.get_one(key), 'splice_acceptor')
    assert hf.buddy2hash(tester) == next_hash


# #####################  '-prg', '--purge' ###################### ##
def test_purge(sb_resources, hf, monkeypatch):
    tester = sb_resources.get_one("p f")
    Sb.pull_recs(tester, "α1[02]")
    bl2seq_output = OrderedDict([('Mle-Panxα10A', OrderedDict([('Mle-Panxα10B', [100.0, 235, 7e-171, 476.0]),
                                                               ('Mle-Panxα12', [56.28, 398, 5e-171, 478.0])])),
                                 ('Mle-Panxα10B', OrderedDict([('Mle-Panxα10A', [100.0, 235, 7e-171, 476.0]),
                                                               ('Mle-Panxα12', [47.51, 381, 2e-128, 366.0])])),
                                 ('Mle-Panxα12', OrderedDict([('Mle-Panxα10A', [56.28, 398, 5e-171, 478.0]),
                                                              ('Mle-Panxα10B', [47.51, 381, 2e-128, 366.0])]))])
    monkeypatch.setattr(Sb, "bl2seq", lambda *_: bl2seq_output)
    Sb.purge(tester, 200)
    assert hf.buddy2hash(tester) == '256681ed87c67f8f3a8c5771572767f1'


# ######################  '-ri', '--rename_ids' ###################### #
hashes = [('d f', '8b4a9e3d3bb58cf8530ee18b9df67ff1'), ('d g', '78c73f97117bd937fd5cf52f4bd6c26e'),
          ('d n', '243024bfd2f686e6a6e0ef65aa963494'), ('d py', '98bb9b57f97555d863054ddb526055b4'),
          ('d pr', '2443c47a712f19099e94fc015dc980a9'), ('d s', '65196fd4f2a4e339e1545f6ed2a6acc3')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_rename_ids(key, next_hash, sb_resources, hf):
    tester = Sb.rename(sb_resources.get_one(key), 'Panx', 'Test', 0)
    assert hf.buddy2hash(tester) == next_hash


def test_rename_ids2(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "Panxα([1-6])", "testing\\1")
    assert hf.buddy2hash(tester) == "f0d7ec055b3b2d9ec11a86634b32e9ef"

    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "Panxα([1-6])", "testing\\1", -1)
    assert hf.buddy2hash(tester) == "f0d7ec055b3b2d9ec11a86634b32e9ef"

    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "[a-z]", "?", 2)
    assert hf.buddy2hash(tester) == "08aaee2e7b9997b512c7d1b2fe748d40"

    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "[a-z]", "?", -2)
    assert hf.buddy2hash(tester) == "9b3946afde20c991099463d099be22e0"

    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "[A-Z]", "?", -20)
    assert hf.buddy2hash(tester) == "451993d7e816881e2700697263b1d8fa"

    tester = sb_resources.get_one("d f")
    Sb.rename(tester, "[a-z]", "?", 2, store_old_id=True)
    assert hf.buddy2hash(tester) == "959fe04c0366c9a143052a02f090707e"

    with pytest.raises(AttributeError) as e:
        Sb.rename(tester, "[a-z]", "\\1?", -2)
    assert "There are more replacement match" in str(e)


# ##################### '-rs', 'replace_subseq' ###################### ##
def test_replace_subsequence(sb_resources, hf):
    tester = sb_resources.get_one("d f")
    Sb.replace_subsequence(tester, "atg(.{5}).{3}", "FOO\\1BAR")
    assert hf.buddy2hash(tester) == "f12707c2b0ef866f0039bac96abb29e0"


# ######################  '-rc', '--reverse_complement' ###################### #
hashes = [('d f', 'e77be24b8a7067ed54f06e0db893ce27'), ('d g', '47941614adfcc5bd107f71abef8b3e00'),
          ('d n', 'f549c8dc076f6b3b4cf5a1bc47bf269d'), ('d py', '0dd20827c011a0a7a0e78881b38ae06a'),
          ('d pr', '0b954f4a263bf38ddeac61ab54f77dc2'), ('d s', '0d6b7deda824b4fc42b65cb87e1d4d14')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_reverse_complement(key, next_hash, sb_resources, hf):
    tester = Sb.reverse_complement(sb_resources.get_one(key))
    assert hf.buddy2hash(tester) == next_hash


def test_reverse_complement_pep_exception(sb_resources):  # Asserts a TypeError will be thrown if user inputs protein
    tester = sb_resources.get_one('p f')
    with pytest.raises(TypeError) as e:
        Sb.reverse_complement(tester)
    assert str(e.value) == "SeqBuddy object is protein. Nucleic acid sequences required."

    prot_record = tester.records[0]
    tester = sb_resources.get_one('d f')
    tester.records.append(prot_record)
    with pytest.raises(TypeError) as e:
        Sb.reverse_complement(tester)
    assert str(e.value) == "Record 'Mle-Panxα12' is protein. Nucleic acid sequences required."


# ######################  '-sfr', '--select_frame' ###################### #
hashes = [('d f', 1, "b831e901d8b6b1ba52bad797bad92d14"), ('d f', 2, "2de033b2bf2327f2795fe425db0bd78f"),
          ('d f', 3, "1c29898d4964e0d1b03207d7e67e1958"), ('d g', 1, "908744b00d9f3392a64b4b18f0db9fee"),
          ('d g', 2, "49c176dec7cc43890a059e0f0f4a9de4"), ('d g', 3, "826d5ae1d4f0ab295d9e39e33999e35f"),
          ('d n', 1, "cb1169c2dd357771a97a02ae2160935d"), ('d n', 2, "87d784f197b55f812d2fc82774da43d1"),
          ('d n', 3, "5d6d2f337ecdc6f9a85e981c975f3e08")]


@pytest.mark.parametrize("key,shift,next_hash", hashes)
def test_select_frame(key, shift, next_hash, sb_resources, hf):
    tester = Sb.select_frame(sb_resources.get_one(key), shift)
    assert hf.buddy2hash(tester) == next_hash


def test_select_frame_edges(sb_resources, hf):
    tester = Sb.select_frame(sb_resources.get_one("d f"), 2)
    temp_file = br.TempFile()
    tester.write(temp_file.path)
    tester = Sb.select_frame(Sb.SeqBuddy(temp_file.path), 1)
    assert hf.buddy2hash(tester) == "b831e901d8b6b1ba52bad797bad92d14"

    tester = Sb.select_frame(sb_resources.get_one("d g"), 2)
    tester = Sb.select_frame(tester, 1)
    assert hf.buddy2hash(tester) == "908744b00d9f3392a64b4b18f0db9fee"

    with pytest.raises(TypeError) as e:  # If protein is input
        Sb.select_frame(sb_resources.get_one("p f"), 2)
    assert "Select frame requires nucleic acid, not protein." in str(e.value)


# ##################### '-ss', 'shuffle_seqs' ###################### ##
hashes = [('d f', 'a86372ed83afc8ba31001919335017bc'), ('d g', '1222dd65c103fa15a50352670dc95f82'),
          ('d n', 'ba022b722620d8688c1f70d535142a5b'), ('d py', '14068e77f8c78f98ff9462116c6781f6'),
          ('d pr', '175821885e688971a818548ff67d0226'), ('d s', '70ba6a80a50c3a99ee3fbb6f5e47c11f'),
          ('p f', '5194c37d5388792119ec988bd7acfbf5'), ('p g', '73396d34998512d1fef4ae6a27c74e49'),
          ('p n', 'e19f273ea91434427c8c3f585fabf1fc'), ('p py', 'fab43d82a2974e8ed0a8b983278ebdd7'),
          ('p pr', 'b6a632a612b1249f29a88383b86d1c1c'), ('p s', 'd34648a0c505644fcd4d29045c9fa502')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_shuffle_seqs(key, next_hash, sb_resources, hf):
    tester = sb_resources.get_one(key)
    tester = Sb.shuffle_seqs(tester, r_seed=12345)
    assert hf.buddy2hash(tester) == next_hash


# ##################### '-sxf', 'split_by_x_files' ###################### ##
def test_split_by_x_files(sb_resources):
    tester = Sb.SeqBuddy(sb_resources.get_one("d f"))
    # File number % Seq number != 0
    sb_list = Sb.split_by_x_files(tester, file_number=3)
    assert len(sb_list) == 3
    counter = 0
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        for record in seqbuddy.records:
            assert record.id == tester.records[counter].id
            assert record.seq == tester.records[counter].seq
            counter += 1

    # File number % Seq number == 0
    tester = Sb.SeqBuddy(sb_resources.get_one("p g"))
    sb_list = Sb.split_by_x_files(tester, file_number=13)
    assert len(sb_list) == 13
    counter = 0
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        for record in seqbuddy.records:
            assert record.id == tester.records[counter].id
            assert record.seq == tester.records[counter].seq
            counter += 1

    # File number > Seq number
    tester = Sb.SeqBuddy(sb_resources.get_one("d py"))
    sb_list = Sb.split_by_x_files(tester, file_number=13)
    assert len(sb_list) == 8


# ##################### '-sxs', 'split_by_x_seqs' ###################### ##
def test_split_by_x_seqs(sb_resources):
    tester = Sb.SeqBuddy(sb_resources.get_one("d f"))
    sb_list = Sb.split_by_x_seqs(tester, seq_number=3)
    assert len(sb_list) == 5
    counter = 0
    for seqbuddy in sb_list:
        assert type(seqbuddy) == Sb.SeqBuddy
        for record in seqbuddy.records:
            assert record.id == tester.records[counter].id
            assert record.seq == tester.records[counter].seq
            counter += 1


# ######################  '-tb', '--taxonomic_breakdown' ###################### #
def test_taxonomic_breakdown(sb_resources):
    tester = sb_resources.get_one("p g")
    assert Sb.taxonomic_breakdown(tester) == """\
Total: 13

Unknown    11
Eukaryota    2
 |Opisthokonta    2
 | |Metazoa    2
 | | |Eumetazoa    2
 | | | |Ctenophora    2
"""
    assert Sb.taxonomic_breakdown(tester, 7) == """\
Total: 13

Unknown    11
Eukaryota    2
 |Opisthokonta    2
 | |Metazoa    2
 | | |Eumetazoa    2
 | | | |Ctenophora    2
 | | | | |Tentaculata    2
 | | | | | |Lobata    2
"""
    assert Sb.taxonomic_breakdown(tester, -7) == """\
Total: 13

Unknown    11
Eukaryota    2
 |Opisthokonta    2
 | |Metazoa    2
 | | |Eumetazoa    2
 | | | |Ctenophora    2
 | | | | |Tentaculata    2
 | | | | | |Lobata    2
"""
    assert Sb.taxonomic_breakdown(tester, 0) == """\
Total: 13

Unknown    11
Eukaryota    2
 |Opisthokonta    2
 | |Metazoa    2
 | | |Eumetazoa    2
 | | | |Ctenophora    2
 | | | | |Tentaculata    2
 | | | | | |Lobata    2
 | | | | | | |Bolinopsidae    2
 | | | | | | | |Mnemiopsis    2
 | | | | | | | | |leidyi    2
"""
    tester = sb_resources.get_one("p f")
    assert Sb.taxonomic_breakdown(tester) == """\
Total: 13

Unknown    13
"""


# ######################  '-tr6', '--translate6frames' ###################### #
def test_translate6frames(sb_resources, hf):
    tester = Sb.translate6frames(sb_resources.get_one("d f"))
    assert hf.buddy2hash(tester) == '95cf24202007399e6ccd6e6f33ae012e'

    tester = Sb.translate6frames(sb_resources.get_one("d g"))
    assert hf.buddy2hash(tester) == 'a32f067e0c0de5928dbcc211bb961532'


def test_translate6frames_pep_exception(sb_resources):
    with pytest.raises(TypeError):
        Sb.translate6frames(sb_resources.get_one("p f"))


# ######################  '-tr', '--translate' ###################### #
hashes = [('d f', '06893e14839dc0448e6f522c1b8f8957'), ('d g', '78a53e66bb4b8f6c26fa2ae0fb29f0ab'),
          ('d n', 'f3339e0193c10427f017dd8f6bd81d7e'), ('r f', '06893e14839dc0448e6f522c1b8f8957'),
          ('r n', 'f3339e0193c10427f017dd8f6bd81d7e')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_translate(key, next_hash, sb_resources, hf):
    tester = Sb.translate_cds(sb_resources.get_one(key))
    assert hf.buddy2hash(tester) == next_hash


def test_translate_ambig(sb_odd_resources, hf):
    tester = Sb.SeqBuddy(sb_odd_resources['ambiguous_dna'])
    tester = Sb.translate_cds(tester)
    assert hf.buddy2hash(tester) == '648ccc7c3400882be5bf6e8d9781f74e'

    tester = Sb.SeqBuddy(sb_odd_resources['ambiguous_rna'])
    tester = Sb.translate_cds(tester)
    assert hf.buddy2hash(tester) == '648ccc7c3400882be5bf6e8d9781f74e'


def test_translate_edges_and_exceptions(capsys, sb_resources, hf):
    with pytest.raises(TypeError):
        Sb.translate_cds(sb_resources.get_one("p f"))

    tester = Sb.SeqBuddy("ATTCGTTAACGCTAGCGTCG", in_format="raw")
    tester = Sb.translate_cds(tester)
    assert str(tester) == ">raw_input\nIR*R*R\n"
    out, err = capsys.readouterr()
    assert err == "Warning: size mismatch between aa and nucl seqs for raw_input --> 20, 6\n"

    tester = Sb.select_frame(sb_resources.get_one("d g"), 3)
    tester = Sb.translate_cds(tester)
    assert hf.buddy2hash(tester) == "85c3bd973cfb683f222388b1529e787f"
    out, err = capsys.readouterr()
    assert hf.string2hash(err) == "9e2a0b4b03f54c209d3a9111792762df"


# ######################  '-tmd', '--transmembrane_domains' ###################### #
def test_transmembrane_domains_pep(sb_resources, hf, monkeypatch, capsys):
    def mock_runtimeerror(*args, **kwargs):
        raise RuntimeError("%s %s" % (args, kwargs))

    def mock_contenttooshorterror(*args, **kwargs):
        raise urllib.error.ContentTooShortError("%s %s" % (args, kwargs), "content")

    def mock_httperror(*args, **kwargs):
        raise urllib.error.HTTPError(url="http://fake.come", code=503, msg=args, hdrs=kwargs, fp="Bar")

    class MockSudsClient(object):
        def __init__(self):
            self.service = MockServiceSelector()

    class MockServiceSelector(object):
        def __init__(self):
            self.job_id_generator = self.job_id_gen()
            self.current_job_id = next(self.job_id_generator)
            self.result_url = "www.something.com"
            self.numseq_str = "doesn't matter"
            self.errinfo = "Did I fail?"
            self.warninfo = "Also doesn't matter"
            self.job_check = self.status()
            self.job_statuses = ["Queued", "Running", "Finished"]

        def status(self):
            while True:
                for next_status in self.job_statuses:
                    yield next_status

        @staticmethod
        def job_id_gen():
            for job_id in ["rst_MFhyxO", "rst_lE27A5"]:
                yield job_id

        def submitjob(self, *args):
            print(args)
            return [[self.current_job_id, self.result_url, self.numseq_str, self.errinfo, self.warninfo]]

        def checkjob(self, jobid):
            print(jobid)
            return [[next(self.job_check), self.result_url, self.errinfo]]

    def mock_urlretrieve(result_url, filename, reporthook):
        assert result_url
        job_id = os.path.split(filename)[-1].split(".")[0]
        if os.path.isfile(join(work_dir.path, "%s.hashmap" % job_id)):
            os.remove(join(work_dir.path, "%s.hashmap" % job_id))
        shutil.copy(join(hf.resource_path, "topcons", "%s.zip" % job_id), work_dir.path)
        for root, _dirs, _files in br.walklevel(work_dir.path):
            print(root)
            print(_files)
        reporthook(2, 10, 100)
        return

    def mock_hash_ids(seqbuddy):
        hashmap = OrderedDict()
        with open(join(hf.resource_path, "topcons", "%s.hashmap" % suds_client.service.current_job_id),
                  "r", encoding="utf-8") as ifile:
            for line in ifile:
                if line:
                    line = line.strip().split("\t")
                    hashmap[line[0]] = line[1]
        seqbuddy.hash_map = hashmap
        for rec, hash_id in zip(seqbuddy.records, list(hashmap.items())):
            rec.id = hash_id[0]
            rec.name = hash_id[0]
        return seqbuddy

    work_dir = br.TempDir()
    keep_dir = br.TempDir()
    repo_dir = br.TempDir()
    data_dir = repo_dir.subdir("buddy_data")
    suds_client = MockSudsClient()

    monkeypatch.setattr(time, "sleep", lambda _: True)
    monkeypatch.setattr(suds.client, "Client", lambda *_, **__: suds_client)
    monkeypatch.setattr(urllib.request, "urlretrieve", mock_urlretrieve)
    monkeypatch.setattr(br, "TempDir", lambda: work_dir)
    monkeypatch.setattr(Sb, "hash_ids", mock_hash_ids)
    monkeypatch.setattr(br, "config_values", lambda *_: {"data_dir": data_dir})

    tester = sb_resources.get_one("d g")
    Sb.pull_recs(tester, "α[56]")
    Sb.delete_features(tester, "splice|TMD")
    capsys.readouterr()
    tester = Sb.transmembrane_domains(tester)
    assert hf.buddy2hash(tester) == "0e199d23cefa502d8fbcd38fabe145f3"

    suds_client.service.current_job_id = next(suds_client.service.job_id_generator)
    tester = sb_resources.get_one("p g")
    Sb.pull_recs(tester, "α[56]")
    Sb.delete_features(tester, "splice|TMD")
    tester = Sb.transmembrane_domains(tester)
    assert hf.buddy2hash(tester) == "808d237e1e84f8f6857ab8766cfbaefc"

    tester = sb_resources.get_one("p g")
    Sb.pull_recs(tester, "α[56]")
    Sb.delete_features(tester, "splice|TMD")
    capsys.readouterr()
    tester = Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    assert hf.buddy2hash(tester) == "808d237e1e84f8f6857ab8766cfbaefc"

    tester = sb_resources.get_one("p g")
    Sb.pull_recs(tester, "α[56]")
    Sb.delete_features(tester, "splice|TMD")
    tester = Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"], keep_temp=keep_dir.path)
    _root, dirs, files = next(br.walklevel(keep_dir.path))

    assert sorted(dirs) == ['rst_MFhyxO', 'rst_lE27A5']
    assert sorted(files) == sorted(['seqs.tmp', os.path.split(work_dir.path)[-1]])

    with pytest.raises(FileNotFoundError) as err:
        Sb.transmembrane_domains(tester, job_ids=["rst_BLAHHH!!"])
    assert "SeqBuddy does not have the necessary hash-map to process job id 'rst_BLAHHH!!'." in str(err)

    suds_client.service.job_statuses = ["Failed"]
    with pytest.raises(ConnectionError) as err:
        Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    assert "Job failed..." in str(err)

    suds_client.service.job_statuses = ["None"]
    with pytest.raises(ConnectionError) as err:
        Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    assert "The job seems to have been lost by the server." in str(err)

    suds_client.service.job_statuses = ["Finished"]
    monkeypatch.setattr(urllib.request, "urlretrieve", mock_runtimeerror)
    Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    out, err = capsys.readouterr()
    assert "Error: Failed to download TOPCONS job rst_lE27A5 after 5 attempts." in err

    monkeypatch.setattr(urllib.request, "urlretrieve", mock_contenttooshorterror)
    Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    out, err = capsys.readouterr()
    assert "Error: Failed to download TOPCONS job rst_lE27A5 after 5 attempts." in err

    monkeypatch.setattr(urllib.request, "urlretrieve", mock_httperror)
    Sb.transmembrane_domains(tester, job_ids=["rst_lE27A5"])
    out, err = capsys.readouterr()
    assert "Error: Failed to download TOPCONS job rst_lE27A5 after 5 attempts." in err
    os.remove("rst_lE27A5.hashmap")
