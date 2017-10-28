#!/usr/bin/env python3
# coding=utf-8
""" tests basic functionality of PhyloBuddy class """
import pytest

from PhyloBuddy import PhyloBuddy, _guess_format
import buddy_resources as br


def mock_raiseattribute(*args, **kwargs):
    raise AttributeError("has no attribute 'NodeStyle': %s, %s" % (args, kwargs))


def test_instantiate_phylobuddy_from_file(pb_resources):
    for key, _path in pb_resources.get("o m k n l", "paths").items():
        in_format = pb_resources.parse_code(key, strict=True)
        in_format = pb_resources.single_letter_codes[in_format["format"][0]]
        assert type(PhyloBuddy(_path, _in_format=in_format)) == PhyloBuddy


def test_instantiate_phylobuddy_from_file_guess(pb_resources):
    for key, _path in pb_resources.get("o m k n l", "paths").items():
        assert type(PhyloBuddy(_path)) == PhyloBuddy


def test_instantiate_phylobuddy_from_handle(pb_resources):
    for key, _path in pb_resources.get("o m k n l", "paths").items():
        with open(_path, 'r') as ifile:
            assert type(PhyloBuddy(ifile)) == PhyloBuddy


def test_instantiate_phylobuddy_from_raw(pb_resources):
    for key, _path in pb_resources.get("o m k n l", "paths").items():
        with open(_path, 'r') as ifile:
            assert type(PhyloBuddy(ifile.read())) == PhyloBuddy


def test_instantiate_phylobuddy_from_phylobuddy(pb_resources):
    for key, pb_obj in pb_resources.get("o m k n l").items():
        tester = PhyloBuddy(pb_obj)
        assert type(PhyloBuddy(tester)) == PhyloBuddy


def test_instantiate_phylobuddy_from_list(pb_resources):
    for key, pb_obj in pb_resources.get("o m k n l").items():
        tester = PhyloBuddy(pb_obj)
        assert type(PhyloBuddy(tester.trees)) == PhyloBuddy


def test_empty_file(pb_odd_resources):
    with open(pb_odd_resources['blank'], "r") as ifile:
        with pytest.raises(SystemExit):
            PhyloBuddy(ifile)


def test_guess_error(pb_odd_resources):
    # File path
    with pytest.raises(br.GuessError):
        PhyloBuddy(pb_odd_resources["unrecognizable"])

    with open(pb_odd_resources["unrecognizable"], 'r') as ifile:
        # Raw
        with pytest.raises(br.GuessError):
            PhyloBuddy(ifile.read())

        # Handle
        with pytest.raises(br.GuessError):
            ifile.seek(0)
            PhyloBuddy(ifile)

    # GuessError output
    test_error = br.GuessError("This is a test")
    assert str(test_error) == "This is a test"

    try:
        PhyloBuddy(pb_odd_resources["unrecognizable"])
    except br.GuessError as e:
        assert "Could not automatically determine the format of" in str(e.value) and \
               "\nTry explicitly setting it with the -f flag." in str(e.value)


def test_phylobuddy_edges(pb_odd_resources):
    # If the input list isn't a list of PhyloBuddy objects
    with pytest.raises(TypeError):
        PhyloBuddy(["Foo", "Bar"])

    # Catch figtree metadata
    assert type(PhyloBuddy(pb_odd_resources["figtree"])) == PhyloBuddy
    with open(pb_odd_resources["figtree"], 'r') as ifile:
        tester = PhyloBuddy(ifile)
        assert type(tester) == PhyloBuddy

    # Unsupported output format
    tester.out_format = "foo"
    with pytest.raises(TypeError):
        str(tester)

    # No trees in PhyloBuddy object
    tester.trees = []
    assert str(tester) == "Error: No trees in object.\n"


# ################################################# HELPER FUNCTIONS ################################################# #
hashes = [('m k', '6843a620b725a3a0e0940d4352f2036f'), ('m n', '543d2fc90ca1f391312d6b8fe896c59c'),
          ('m l', '6ce146e635c20ad62e21a1ed6fddbd3a'), ('o k', '4dfed97b2a23b8957ee5141bf4681fe4'),
          ('o n', '77d00fdc512fa09bd1146037d25eafa0'), ('o l', '9b1014be1b38d27f6b7ef73d17003dae')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_str(key, next_hash, pb_resources, hf):
    tester = pb_resources.get_one(key)
    assert hf.string2hash(str(tester)) == next_hash


@pytest.mark.parametrize("key,next_hash", hashes)
def test_write1(key, next_hash, pb_resources, hf):
    temp_file = br.TempFile()
    tester = pb_resources.get_one(key)
    tester.write(temp_file.path)
    out = "{0}\n".format(temp_file.read().rstrip())
    assert hf.string2hash(out) == next_hash


def test_guess_format(pb_resources):
    guessed_format = _guess_format([])
    assert guessed_format == "newick"

    guessed_format = _guess_format(pb_resources.get_one("o n"))
    assert guessed_format == "nexus"

    guessed_format = _guess_format(pb_resources.get_one("o l", mode="paths"))
    assert guessed_format == "nexml"

    with pytest.raises(br.GuessError):
        _guess_format(dict)
