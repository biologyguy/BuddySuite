#!/usr/bin/env python3
# coding=utf-8
""" tests basic functionality of PhyloBuddy class """
import pytest

try:
    from buddysuite import MyFuncs
    from buddysuite.PhyloBuddy import PhyloBuddy
    from buddysuite.buddy_resources import GuessError
except ImportError:
    import MyFuncs
    from PhyloBuddy import PhyloBuddy
    from buddy_resources import GuessError


def test_instantiate_phylobuddy_from_file(pb_resources):
    for key, _path in pb_resources.get("o m k n l", "paths").items():
        in_format = pb_resources.parse_code(key, strict=True)
        in_format = pb_resources.single_letter_codes[in_format["format"][0]]
        assert type(PhyloBuddy(_path, _in_format=in_format)) == PhyloBuddy

"""
@pytest.mark.parametrize("phylo_file", phylo_files)
def test_instantiate_phylobuddy_from_file_guess(phylo_file):
    assert type(Pb.PhyloBuddy(resource(phylo_file))) == Pb.PhyloBuddy


@pytest.mark.parametrize("phylo_file", phylo_files)
def test_instantiate_phylobuddy_from_handle(phylo_file):
    with open(resource(phylo_file), 'r') as ifile:
        assert type(Pb.PhyloBuddy(ifile)) == Pb.PhyloBuddy


@pytest.mark.parametrize("phylo_file", phylo_files)
def test_instantiate_phylobuddy_from_raw(phylo_file):
    with open(resource(phylo_file), 'r') as ifile:
        assert type(Pb.PhyloBuddy(ifile.read())) == Pb.PhyloBuddy


@pytest.mark.parametrize("phylo_file", phylo_files)
def test_instantiate_phylobuddy_from_phylobuddy(phylo_file):
    tester = Pb.PhyloBuddy(resource(phylo_file))
    assert type(Pb.PhyloBuddy(tester)) == Pb.PhyloBuddy


@pytest.mark.parametrize("phylo_file", phylo_files)
def test_instantiate_phylobuddy_from_list(phylo_file):
    tester = Pb.PhyloBuddy(resource(phylo_file))
    assert type(Pb.PhyloBuddy(tester.trees)) == Pb.PhyloBuddy


def test_empty_file():
    with open(resource("blank.fa"), "r") as ifile:
        with pytest.raises(SystemExit):
            Pb.PhyloBuddy(ifile)


def test_guess_error():
    # File path
    with pytest.raises(br.GuessError):
        Pb.PhyloBuddy(resource("unrecognizable.txt"))

    with open(resource("unrecognizable.txt"), 'r') as ifile:
        # Raw
        with pytest.raises(br.GuessError):
            Pb.PhyloBuddy(ifile.read())

        # Handle
        with pytest.raises(br.GuessError):
            ifile.seek(0)
            Pb.PhyloBuddy(ifile)

    # GuessError output
    test_error = br.GuessError("This is a test")
    assert str(test_error) == "This is a test"

    try:
        Pb.PhyloBuddy(resource("unrecognizable.txt"))
    except br.GuessError as e:
        assert "Could not automatically determine the format of" in str(e.value) and \
               "\nTry explicitly setting it with the -f flag." in str(e.value)


def test_stderr(capsys):
    Pb._stderr("Hello std_err", quiet=False)
    out, err = capsys.readouterr()
    assert err == "Hello std_err"

    Pb._stderr("Hello std_err", quiet=True)
    out, err = capsys.readouterr()
    assert err == ""


def test_stdout(capsys):
    Pb._stdout("Hello std_out", quiet=False)
    out, err = capsys.readouterr()
    assert out == "Hello std_out"

    Pb._stdout("Hello std_out", quiet=True)
    out, err = capsys.readouterr()
    assert out == ""


def test_phylobuddy_edges():
    # If the input list isn't a list of PhyloBuddy objects
    with pytest.raises(TypeError):
        Pb.PhyloBuddy(["Foo", "Bar"])

    # Catch figtree metadata
    assert type(Pb.PhyloBuddy(resource("figtree.nexus"))) == Pb.PhyloBuddy
    with open(resource("figtree.nexus"), 'r') as ifile:
        tester = Pb.PhyloBuddy(ifile)
        assert type(tester) == Pb.PhyloBuddy

    # Unsupported output format
    tester.out_format = "foo"
    with pytest.raises(TypeError):
        str(tester)

    # No trees in PhyloBuddy object
    tester.trees = []
    assert str(tester) == "Error: No trees in object.\n"

pb_objects = [Pb.PhyloBuddy(resource(x)) for x in phylo_files]

# ################################################# HELPER FUNCTIONS ################################################# #
hashes = ['6843a620b725a3a0e0940d4352f2036f', '543d2fc90ca1f391312d6b8fe896c59c', '6ce146e635c20ad62e21a1ed6fddbd3a',
          '4dfed97b2a23b8957ee5141bf4681fe4', '77d00fdc512fa09bd1146037d25eafa0', '9b1014be1b38d27f6b7ef73d17003dae']

hashes = [(Pb.make_copy(pb_objects[x]), hashes[x]) for x in range(len(pb_objects))]


@pytest.mark.parametrize("phylobuddy,next_hash", hashes)
def test_str(phylobuddy, next_hash):
    tester = str(phylobuddy)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash


@pytest.mark.parametrize("phylobuddy,next_hash", hashes)
def test_write1(phylobuddy, next_hash):
    temp_file = MyFuncs.TempFile()
    phylobuddy.write(temp_file.path)
    out = "{0}\n".format(temp_file.read().rstrip())
    tester = md5(out.encode()).hexdigest()
    assert tester == next_hash


def test_convert_to_ete():
    tester = Pb.make_copy(pb_objects[0])
    tester.trees[0].seed_node.annotations.add_new("pb_color", '#ff0000')
    ete_tree = Pb._convert_to_ete(tester.trees[0])
    assert ete_tree.pb_color == '#ff0000'


def test_guess_format():
    with pytest.raises(br.GuessError):
        Pb._guess_format(dict)
"""