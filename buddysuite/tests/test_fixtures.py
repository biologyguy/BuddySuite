# -*- coding: utf-8 -*-

"""
verify that fixtures are working as expected. I.e., test the tests before testing ;)
"""
import os
import pytest
from buddysuite import SeqBuddy, AlignBuddy, PhyloBuddy, DatabaseBuddy


# #################################  -  Helper functions  -  ################################## #
def test_hf_buddy2hash(hf):
    seqbuddy = SeqBuddy.SeqBuddy(">Foo\nATGCGCGCATGCTA")
    assert hf.buddy2hash(seqbuddy) == "0b7482f09f950144574337d95376d644"

    alignbuddy = AlignBuddy.AlignBuddy(">Foo1\nATGCGCGCATGCTA>Foo2\nATGGGCGAATTTTA")
    assert hf.buddy2hash(alignbuddy) == "71b794705f9817006cd053bd20fb1481"

    phylobuddy = PhyloBuddy.PhyloBuddy("(((A, B), C),(D, E));")
    assert hf.buddy2hash(phylobuddy) == "20f81e9f99c673e6cde3afb4b30cc6da"

    dbbuddy = DatabaseBuddy.DbBuddy("Casp9,Inx1")
    with open("temp.del", "w") as ofile:
        ofile.write(str(dbbuddy))
    assert hf.buddy2hash(dbbuddy) == "058055f7f09d0e8bcf8ae75d3ed73a1f"

    with pytest.raises(AttributeError) as err:
        hf.buddy2hash("I'm not a buddy object!")
    assert "Buddy object required" in str(err)


def test_hf_string2hash(hf):
    assert hf.string2hash("Foo bar baz") == "520c28a8ac3459af817a1abfb3bd152e"


def test_hf_res_path(hf):
    assert "unit_test_resources" in hf.resource_path


# #################################  -  SeqBuddy  -  ################################## #
def test_sb_resources_init(hf, sb_resources, capsys):
    assert sb_resources
    for molecule in [('dna', 11), ('rna', 2), ('pep', 8)]:
        assert molecule[0] in sb_resources.resources
        assert len(sb_resources.resources[molecule[0]]) == molecule[1]
        assert molecule[0] in sb_resources.sb_objs

    assert sb_resources.resources['dna']["clustal"] == "%s/Mnemiopsis_cds.clus" % hf.resource_path
    assert type(sb_resources.sb_objs['dna']["clustal"]) == SeqBuddy.SeqBuddy

    for key in ["molecule", "format"]:
        assert key in sb_resources.code_dict
    for key in ["p", "d", "r"]:
        assert key in sb_resources.code_dict["molecule"]
    for key in ["c", "f", "g", "n", "py", "pr", "pss", "psr", "s", "e", "x"]:
        assert key in sb_resources.code_dict["format"]

    for key in ["p", "d", "r", "c", "f", "g", "n", "py", "pr", "pss", "psr", "s", "e", "x"]:
        assert key in sb_resources.single_letter_codes


def test_sb_resources_parse_code(sb_resources):
    codes = sb_resources.parse_code()
    assert sorted(codes["molecule"]) == ['d', 'p', 'r']
    assert sorted(codes["format"]) == ['c', 'e', 'f', 'g', 'n', 'pr', 'psr', 'pss', 'py', 's', 'x']

    codes = sb_resources.parse_code("d")
    assert codes["molecule"] == ['d']
    assert sorted(codes["format"]) == ['c', 'e', 'f', 'g', 'n', 'pr', 'psr', 'pss', 'py', 's', 'x']

    codes = sb_resources.parse_code("c")
    assert sorted(codes["molecule"]) == ['d', 'p', 'r']
    assert codes["format"] == ['c']

    codes = sb_resources.parse_code('d r c f n')
    assert sorted(codes["molecule"]) == ['d', 'r']
    assert sorted(codes["format"]) == ['c', 'f', 'n']


def test_sb_resources_get_key(sb_resources):
    key = sb_resources.get_key("e d")
    assert key["molecule"] == "dna"
    assert key["format"] == "embl"

    with pytest.raises(AttributeError) as err:
        sb_resources.get_key("d e f")
    assert "Only explicit two-component codes are accepted" in str(err)

    with pytest.raises(AttributeError) as err:
        sb_resources.get_key("d z")
    assert "Malformed letter code, 'z' not recognized" in str(err)
    with pytest.raises(AttributeError) as err:
        sb_resources.get_key("d p")
    assert "Malformed letter code, trying to append multiple values to the molecule component of the key" in str(err)


def test_sb_resources_get(sb_resources):
    objs = sb_resources.get("d p pss")
    assert len(objs) == 2
    assert type(objs["d pss"]) == SeqBuddy.SeqBuddy
    assert type(objs["p pss"]) == SeqBuddy.SeqBuddy

    objs = sb_resources.get("d p pss", mode="paths")
    assert len(objs) == 2
    assert "Mnemiopsis_cds.physs" in objs["d pss"]
    assert "Mnemiopsis_pep.physs" in objs["p pss"]

    objs = sb_resources.get("d z pss", mode="objs")
    assert len(objs) == 1
    assert "d pss" in objs

    with pytest.raises(ValueError) as err:
        sb_resources.get('d p pss', mode="foo")
    assert "The 'mode' parameter only accepts 'objs' or 'paths' as input." in str(err)

    assert len(sb_resources.get('z f')) == 3


def test_sb_resources_get_list(sb_resources):
    objs = sb_resources.get_list("d p pss")
    assert len(objs) == 2
    assert type(objs[0]) == SeqBuddy.SeqBuddy
    assert type(objs[1]) == SeqBuddy.SeqBuddy

    objs = sb_resources.get_list("d p pss", mode="paths")
    assert len(objs) == 2
    assert "Mnemiopsis_cds.physs" in objs[0]
    assert "Mnemiopsis_pep.physs" in objs[1]

    objs = sb_resources.get_list("d z pss", mode="objs")
    assert len(objs) == 1
    assert type(objs[0]) == SeqBuddy.SeqBuddy

    with pytest.raises(ValueError) as err:
        sb_resources.get_list('d p pss', mode="foo")
    assert "The 'mode' parameter only accepts 'objs' or 'paths' as input." in str(err)


def test_sb_resources_get_one(sb_resources):
    assert not sb_resources.get_one("d p pss")
    assert not sb_resources.get_one("z pss")

    objs = sb_resources.get_one("d pss")
    assert type(objs) == SeqBuddy.SeqBuddy

    objs = sb_resources.get_one("p pss", mode="paths")
    assert "Mnemiopsis_pep.physs" in objs

    with pytest.raises(ValueError) as err:
        sb_resources.get_one('d pss', mode="foo")
    assert "The 'mode' parameter only accepts 'objs' or 'paths' as input." in str(err)


def test_sb_odd_resources(sb_odd_resources):
    assert len(sb_odd_resources) == 13
    for key in ["blank", "figtree", "unrecognizable", "gibberish", "phylipss_cols", "duplicate",
                "ambiguous_dna", "ambiguous_rna", "blastn", "blastp", "dummy_feats", "cnidaria_pep", "mixed"]:
        assert key in sb_odd_resources


# #################################  -  AlignBuddy  -  ################################ #
def test_alb_resources_init(hf, alb_resources, capsys):
    assert alb_resources
    assert 'dna' in alb_resources.resources
    assert 'pep' in alb_resources.resources
    assert 'rna' in alb_resources.resources

    for files in [('dna', 'single', 9), ('dna', 'multi', 6), ('rna', 'single', 1),
                  ('rna', 'multi', 0), ('pep', 'single', 7), ('pep', 'multi', 6)]:
        assert len(alb_resources.resources[files[0]][files[1]]) == files[2]

    assert alb_resources.resources['dna']['single']['clustal'] == "%s/Mnemiopsis_cds.clus" % hf.resource_path
    assert type(alb_resources.alb_objs['dna']['single']['clustal']) == AlignBuddy.AlignBuddy

    for key in ["molecule", "format", "num_aligns"]:
        assert key in alb_resources.code_dict
    for key in ["p", "d", "r"]:
        assert key in alb_resources.code_dict["molecule"]
    for key in ["o", "m"]:
        assert key in alb_resources.code_dict["num_aligns"]
    for key in ["c", "f", "g", "n", "py", "pr", "pss", "psr", "s"]:
        assert key in alb_resources.code_dict["format"]

    for key in ["p", "d", "r", "o", "m", "c", "f", "g", "n", "py", "pr", "pss", "psr", "s"]:
        assert key in alb_resources.single_letter_codes


def test_alb_resources_parse_code(alb_resources):
    codes = alb_resources.parse_code()
    assert sorted(codes["molecule"]) == ['d', 'p', 'r']
    assert sorted(codes["num_aligns"]) == ['m', 'o']
    assert sorted(codes["format"]) == ['c', 'f', 'g', 'n', 'pr', 'psr', 'pss', 'py', 's']

    codes = alb_resources.parse_code("d")
    assert codes["molecule"] == ['d']
    assert sorted(codes["num_aligns"]) == ['m', 'o']
    assert sorted(codes["format"]) == ['c', 'f', 'g', 'n', 'pr', 'psr', 'pss', 'py', 's']

    codes = alb_resources.parse_code("m")
    assert sorted(codes["molecule"]) == ['d', 'p', 'r']
    assert codes["num_aligns"] == ['m']
    assert sorted(codes["format"]) == ['c', 'f', 'g', 'n', 'pr', 'psr', 'pss', 'py', 's']

    codes = alb_resources.parse_code("c")
    assert sorted(codes["molecule"]) == ['d', 'p', 'r']
    assert sorted(codes["num_aligns"]) == ['m', 'o']
    assert codes["format"] == ['c']

    codes = alb_resources.parse_code('d r o c f n')
    assert sorted(codes["molecule"]) == ['d', 'r']
    assert codes["num_aligns"] == ['o']
    assert sorted(codes["format"]) == ['c', 'f', 'n']


def test_alb_resources_get_key(alb_resources):
    key = alb_resources.get_key("d o f")
    assert key["molecule"] == "dna"
    assert key["num_aligns"] == "single"
    assert key["format"] == "fasta"

    with pytest.raises(AttributeError) as err:
        alb_resources.get_key("d o f p")
    assert "Only explicit three-component codes are accepted" in str(err)

    with pytest.raises(AttributeError) as err:
        alb_resources.get_key("d o z")
    assert "Malformed letter code, 'z' not recognized" in str(err)

    # ToDo: Fix the get_key method (it's not raising AttributeErrors correctly...)
    # with pytest.raises(AttributeError) as err:
    #    alb_resources.get_key("d p o")
    # assert "Malformed letter code, trying to append multiple values to the molecule component of the key" in str(err)


def test_alb_resources_get(alb_resources):
    objs = alb_resources.get("d p pss")
    assert len(objs) == 4
    assert type(objs["o d pss"]) == AlignBuddy.AlignBuddy
    assert type(objs["o p pss"]) == AlignBuddy.AlignBuddy
    assert type(objs["m d pss"]) == AlignBuddy.AlignBuddy
    assert type(objs["m p pss"]) == AlignBuddy.AlignBuddy

    objs = alb_resources.get("d p o pss", mode="paths")
    assert len(objs) == 2
    assert "Mnemiopsis_cds.physs" in objs["o d pss"]
    assert "Mnemiopsis_pep.physs" in objs["o p pss"]

    objs = alb_resources.get("d z pss", mode="objs")
    assert len(objs) == 2
    assert "o d pss" in objs
    assert "m d pss" in objs

    with pytest.raises(ValueError) as err:
        alb_resources.get('d p pss', mode="foo")
    assert "The 'mode' parameter only accepts 'objs' or 'paths' as input." in str(err)

    assert len(alb_resources.get('z py')) == 4


def test_alb_resources_get_list(alb_resources):
    objs = alb_resources.get_list("d p pss")
    assert len(objs) == 4
    assert type(objs[0]) == AlignBuddy.AlignBuddy
    assert type(objs[2]) == AlignBuddy.AlignBuddy

    objs = alb_resources.get_list("d p o pss", mode="paths")
    assert len(objs) == 2
    assert "Mnemiopsis_cds.physs" in objs[0]
    assert "Mnemiopsis_pep.physs" in objs[1]

    objs = alb_resources.get_list("d z pss", mode="objs")
    assert len(objs) == 2
    assert type(objs[0]) == AlignBuddy.AlignBuddy
    assert type(objs[1]) == AlignBuddy.AlignBuddy

    with pytest.raises(ValueError) as err:
        alb_resources.get_list('d p pss', mode="foo")
    assert "The 'mode' parameter only accepts 'objs' or 'paths' as input." in str(err)


def test_alb_resources_get_one(alb_resources):
    assert not alb_resources.get_one("d p pss")
    assert not alb_resources.get_one("d z pss")

    objs = alb_resources.get_one("d m pss")
    assert type(objs) == AlignBuddy.AlignBuddy

    objs = alb_resources.get_one("p m pss", mode="paths")
    assert "Alignments_pep.physs" in objs

    with pytest.raises(ValueError) as err:
        alb_resources.get_one('d m pss', mode="foo")
    assert "The 'mode' parameter only accepts 'objs' or 'paths' as input." in str(err)


def test_alb_odd_resources(alb_odd_resources):
    assert len(alb_odd_resources) == 3
    assert 'blank' in alb_odd_resources

    assert len(alb_odd_resources['dna']) == 1
    assert len(alb_odd_resources['protein']) == 1

    assert len(alb_odd_resources['dna']['single']) == 4
    for key in ["fasta", "phylipss_recs", "phylipss_cols", "ambiguous"]:
        assert key in alb_odd_resources['dna']['single']

    assert len(alb_odd_resources['protein']['single']) == 1
    assert 'phylip' in alb_odd_resources['protein']['single']


# ################################  -  PhyloBuddy  -  ################################# #
