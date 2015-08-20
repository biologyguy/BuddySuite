import pytest
from hashlib import md5
import os
from io import StringIO
import re
from copy import deepcopy
from MyFuncs import TempFile

try:
    import workshop.PhyloBuddy as Pb
except ImportError:
    import PhyloBuddy as Pb

def phylo_to_hash(_phylobuddy, mode='hash'):
    if mode != "hash":
        return str(_phylobuddy)
    _hash = md5("{0}\n".format(str(_phylobuddy).rstrip()).encode()).hexdigest()
    return _hash

root_dir = os.getcwd()

def resource(file_name):
    return "{0}/unit_test_resources/{1}".format(root_dir, file_name)

phylo_files = ['multi_tree.newick', 'multi_tree.nex', 'multi_tree.xml', 'single_tree.newick', 'single_tree.nex',
            'single_tree.xml']

file_types = ['newick', 'nexus', 'nexml', 'newick', 'nexus', 'nexml']

@pytest.mark.parametrize("phylo_file,file_type", [(phylo_files[x], file_types[x]) for x in range(len(phylo_files))])
def test_instantiate_phylobuddy_from_file(phylo_file, file_type):
    assert type(Pb.PhyloBuddy(resource(phylo_file), _in_format=file_type)) == Pb.PhyloBuddy


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
    with pytest.raises(Pb.GuessError):
        Pb.PhyloBuddy(resource("unrecognizable.txt"))

    with open(resource("unrecognizable.txt"), 'r') as ifile:
        # Raw
        with pytest.raises(Pb.GuessError):
            Pb.PhyloBuddy(ifile.read())

        # Handle
        with pytest.raises(Pb.GuessError):
            ifile.seek(0)
            Pb.PhyloBuddy(ifile)

    # GuessError output
    try:
        Pb.PhyloBuddy(resource("unrecognizable.txt"))
    except Pb.GuessError as e:
        assert "Could not determine format from _input file" in str(e.value) and \
               "\nTry explicitly setting with -f flag." in str(e.value)


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

pb_objects = [Pb.PhyloBuddy(resource(x)) for x in phylo_files]

hashes = ['6843a620b725a3a0e0940d4352f2036f', '543d2fc90ca1f391312d6b8fe896c59c', '6ce146e635c20ad62e21a1ed6fddbd3a',
          '4dfed97b2a23b8957ee5141bf4681fe4', '77d00fdc512fa09bd1146037d25eafa0', '9b1014be1b38d27f6b7ef73d17003dae']

hashes = [(pb_objects[x], hashes[x]) for x in range(len(pb_objects))]

@pytest.mark.parametrize("phylobuddy,next_hash", hashes)
def test_print(phylobuddy, next_hash, capsys):
    phylobuddy.print()
    out, err = capsys.readouterr()
    out = "{0}\n".format(out.rstrip())
    tester = md5(out.encode()).hexdigest()
    assert tester == next_hash


@pytest.mark.parametrize("phylobuddy,next_hash", hashes)
def test_str(phylobuddy, next_hash):
    tester = str(phylobuddy)
    tester = md5(tester.encode()).hexdigest()
    assert tester == next_hash


@pytest.mark.parametrize("phylobuddy,next_hash", hashes)
def test_write1(phylobuddy, next_hash):
    temp_file = TempFile()
    phylobuddy.write(temp_file.path)
    out = "{0}\n".format(temp_file.read().rstrip())
    tester = md5(out.encode()).hexdigest()
    assert tester == next_hash

def test_split_polytomies():
    tester = Pb.PhyloBuddy('(A,(B,C,D));')
    Pb.split_polytomies(tester)
    assert str(tester) in ['(A:1.0,(B:1.0,(C:1.0,D:1.0)):1.0):1.0;\n', '(A:1.0,(B:1.0,(D:1.0,C:1.0)):1.0):1.0;\n',
                           '(A:1.0,(C:1.0,(B:1.0,D:1.0)):1.0):1.0;\n', '(A:1.0,(C:1.0,(D:1.0,B:1.0)):1.0):1.0;\n',
                           '(A:1.0,(D:1.0,(C:1.0,B:1.0)):1.0):1.0;\n', '(A:1.0,(D:1.0,(B:1.0,C:1.0)):1.0):1.0;\n']


# ############################# 'pr', '--prune_taxa' ############################# #
pt_hashes = ['99635c6dbf708f94cf4dfdca87113c44', 'fc03b4f100f038277edf6a9f48913dd0', '001db76033cba463a0f187266855e8d5']
pt_hashes = [(Pb._make_copies(pb_objects[x]), next_hash) for x, next_hash in enumerate(pt_hashes)]
@pytest.mark.parametrize("phylobuddy, next_hash", pt_hashes)
def test_prune_taxa(phylobuddy, next_hash):
    Pb.prune_taxa(phylobuddy, 'fir')
    assert phylo_to_hash(phylobuddy) == next_hash


# ############################# 'li', '--list_ids' ############################# #
li_hashes = ['188aaa14889c438a66bcf6e02ea77eb7', 'f06f615e2bf0c4aa91c982aa8956cf12', 'b99e5b92ec1ab59ad85e94e43131bb78']
li_hashes = [(Pb._make_copies(pb_objects[x]), next_hash) for x, next_hash in enumerate(li_hashes)]
@pytest.mark.parametrize("phylobuddy, next_hash", li_hashes)
def test_list_ids(phylobuddy, next_hash):
    tester = str(Pb.list_ids(phylobuddy))
    assert md5(tester.encode()).hexdigest() == next_hash


# ############################# 'cd', '--calculate_distance' ############################# #
cd_hashes = ['3c49c6a7f06244c0b5d45812f6791519', 'df0d56fe8bdc120a9898cbdf186182f6', '7df609f2e6ee613d3bf3c3d2aae26ad4']
cd_hashes = [(Pb._make_copies(pb_objects[x]), next_hash) for x, next_hash in enumerate(cd_hashes)]
@pytest.mark.parametrize("phylobuddy, next_hash", cd_hashes)
def test_calculate_distance_wrf(phylobuddy, next_hash):
    tester = str(Pb.calculate_distance(phylobuddy, _method='wrf'))
    assert md5(tester.encode()).hexdigest() == next_hash

cd_hashes = ['6d087b86aa9f5bc5013113972173fe0f', 'ec0b22eb92836fc0666ce89e89c17336', '7ef096e3c32dbf898d4b1a035d5c9ad4']
cd_hashes = [(Pb._make_copies(pb_objects[x]), next_hash) for x, next_hash in enumerate(cd_hashes)]
@pytest.mark.parametrize("phylobuddy, next_hash", cd_hashes)
def test_calculate_distance_uwrf(phylobuddy, next_hash):
    tester = str(Pb.calculate_distance(phylobuddy, _method='uwrf'))
    assert md5(tester.encode()).hexdigest() == next_hash

cd_hashes = ['3dba6b10fdd04505b4e4482d926b67d3', 'a36937f19d7398dabe4309046506f1fe', '8d0b3a035015d62916b525f371684bf8']
cd_hashes = [(Pb._make_copies(pb_objects[x]), next_hash) for x, next_hash in enumerate(cd_hashes)]
@pytest.mark.parametrize("phylobuddy, next_hash", cd_hashes)
def test_calculate_distance_ed(phylobuddy, next_hash):
    tester = str(Pb.calculate_distance(phylobuddy, _method='ed'))
    assert md5(tester.encode()).hexdigest() == next_hash