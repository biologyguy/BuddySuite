import pytest
from hashlib import md5
import os
import sys
import argparse
from copy import deepcopy

sys.path.insert(0, "./")
import buddy_resources as br
from MyFuncs import TempFile, TempDir
import PhyloBuddy as Pb
import AlignBuddy as Alb
import MyFuncs

VERSION = Pb.VERSION
WRITE_FILE = MyFuncs.TempFile()


def fmt(prog):
    return br.CustomHelpFormatter(prog)

parser = argparse.ArgumentParser(prog="PhyloBuddy.py", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                 description='''\
\033[1mPhyloBuddy\033[m
Put a little bonsai into your phylogeny.

\033[1mUsage examples\033[m:
PhyloBuddy.py "/path/to/tree_file" -<cmd>
PhyloBuddy.py "/path/to/tree_file" -<cmd> | PhyloBuddy.py -<cmd>
PhyloBuddy.py "(A,(B,C));" -f "raw" -<cmd>
''')

br.flags(parser, ("trees", "Supply file path(s) or raw tree string, If piping trees into PhyloBuddy "
                           "this argument can be left blank."), br.pb_flags, br.pb_modifiers, VERSION)

in_args = parser.parse_args()


def phylo_to_hash(_phylobuddy, mode='hash'):
    if mode != "hash":
        return str(_phylobuddy)
    _hash = md5("{0}\n".format(str(_phylobuddy).rstrip()).encode()).hexdigest()
    return _hash

root_dir = os.getcwd()


def command_line_output_hash(output):
    return md5(output.encode()).hexdigest()


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
    test_error = Pb.GuessError("This is a test")
    assert str(test_error) == "This is a test"

    try:
        Pb.PhyloBuddy(resource("unrecognizable.txt"))
    except Pb.GuessError as e:
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

hashes = [(Pb._make_copy(pb_objects[x]), hashes[x]) for x in range(len(pb_objects))]


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


def test_convert_to_ete():
    tester = Pb._make_copy(pb_objects[0])
    tester.trees[0].seed_node.annotations.add_new("pb_color", '#ff0000')
    ete_tree = Pb._convert_to_ete(tester.trees[0])
    assert ete_tree.pb_color == '#ff0000'


def test_format_to_extension():
    assert Pb._format_to_extension("nexus") == "nex"


def test_guess_format():
    with pytest.raises(Pb.GuessError):
        Pb._guess_format(dict)


# ################################################ MAIN API FUNCTIONS ################################################ #
# ###################### 'ct', '--consensus_tree' ###################### #
hashes = ["f20cbd5aae5971cce8efbda15e4e0b7e", "109c2be4172dc15a7dad477d0041933f",
          "5a813bf69c5f3202842b15c3d7d19846", "59956aecc03e49a4d3cce8505b80b088"]
hashes = [(Pb._make_copy(pb_objects[x]), hashes[x]) for x in range(4)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_consensus_tree(phylobuddy, next_hash):
    tester = Pb.consensus_tree(phylobuddy)
    assert phylo_to_hash(tester) == next_hash


hashes = ["447862ed1ed6e98f2fb535ecce70218b", "253aef020f7b5bc36ed03fe2ad8f0837",
          "77afcc1db9757c9f29a195e7ef380c0c", "59956aecc03e49a4d3cce8505b80b088"]
hashes = [(Pb._make_copy(pb_objects[x]), hashes[x]) for x in range(4)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_consensus_tree_95(phylobuddy, next_hash):
    tester = Pb.consensus_tree(phylobuddy, frequency=0.95)
    assert phylo_to_hash(tester) == next_hash


# ###################### 'dt', '--display_trees' ###################### #
# def test_display_trees():
    # foo = Pb.display_trees(pb_objects[0])


def test_display_trees_error():
    if "DISPLAY" in os.environ:
        display_value = os.environ["DISPLAY"]
        del os.environ["DISPLAY"]
    else:
        display_value = None

    with pytest.raises(SystemError):
        Pb.display_trees(Pb._make_copy(pb_objects[0]))

    if display_value:
        os.environ["DISPLAY"] = display_value


# ###################### 'dis', '--distance' ###################### #
hashes = ['e39a8aadaec77680ad0d9004bab824ea', '3c49c6a7f06244c0b5d45812f6791519', '7df609f2e6ee613d3bf3c3d2aae26ad4']
hashes = [(Pb._make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_distance_wrf(phylobuddy, next_hash):
    tester = str(Pb.distance(phylobuddy, method='wrf'))
    assert md5(tester.encode()).hexdigest() == next_hash

hashes = ['c15d06fc5344da3149e19b134ca31c62', '6d087b86aa9f5bc5013113972173fe0f', '7ef096e3c32dbf898d4b1a035d5c9ad4']
hashes = [(Pb._make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_distance_uwrf(phylobuddy, next_hash):
    tester = str(Pb.distance(phylobuddy, method='uwrf'))
    assert md5(tester.encode()).hexdigest() == next_hash

hashes = ['68942718c8baf4e4bdf5dd2992fbbf9d', '3dba6b10fdd04505b4e4482d926b67d3', '8d0b3a035015d62916b525f371684bf8']
hashes = [(Pb._make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_distance_ed(phylobuddy, next_hash):
    tester = str(Pb.distance(phylobuddy, method='ed'))
    assert md5(tester.encode()).hexdigest() == next_hash


def test_distance_unknown_method():
    with pytest.raises(AttributeError):
        Pb.distance(Pb._make_copy(pb_objects[0]), method='foo')


# ######################  'gt', '--generate_trees' ###################### #
# Hashes are for RAxML version 8.2.3
@pytest.mark.generate_trees
def test_raxml_inputs():
    # Nucleotide
    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))
    tester = Pb.generate_tree(tester, 'raxml')
    assert phylo_to_hash(tester) == '706ba436f8657ef3aee7875217dd07c0'
    # Peptide
    tester = Alb.AlignBuddy(resource("Mnemiopsis_pep.nex"))
    tester = Pb.generate_tree(tester, 'raxml')
    assert phylo_to_hash(tester) == 'fc35569091eeba49ac4dcec7fc6890bf'


@pytest.mark.generate_trees
def test_raxml_multi_param():
    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))
    tester = Pb.generate_tree(tester, 'raxml', '-m GTRCAT -p 112358 -K MK -N 2')
    assert phylo_to_hash(tester) == '2bce58a9c6756fa68fd828a307850d7d'


# PhyML version 20120412
@pytest.mark.generate_trees
def test_phyml_inputs():
    # Nucleotide
    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))
    tester = Pb.generate_tree(tester, 'phyml', '-m GTR --r_seed 12345')
    assert phylo_to_hash(tester) == 'd3a4e7601998885f333ddd714ca764db'
    # Peptide
    tester = Alb.AlignBuddy(resource("Mnemiopsis_pep.nex"))
    tester = Pb.generate_tree(tester, 'phyml', '-m Blosum62 --r_seed 12345')
    assert phylo_to_hash(tester) == '52c7d028341b250bcc867d57a68c794c'


@pytest.mark.generate_trees
def test_phyml_multi_param():
    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))
    tester = Pb.generate_tree(tester, 'phyml', '-m GTR -o tl -b 2 --r_seed 12345')
    assert phylo_to_hash(tester) == '5434f29509eab76dd52dd69d2c0e186f'


@pytest.mark.generate_trees
def test_fasttree_inputs():
    temp_dir = TempDir()
    # Nucleotide
    alignbuddy = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))

    tester = Pb.generate_tree(alignbuddy, 'fasttree', '-seed 12345')
    assert phylo_to_hash(tester) == 'd7f505182dd1a1744b45cc326096f70c'

    alignbuddy = Alb.AlignBuddy(resource("Mnemiopsis_pep.nex"))
    tester = Pb.generate_tree(alignbuddy, 'fasttree', '-seed 12345', keep_temp="%s/new_dir" % temp_dir.path)
    assert phylo_to_hash(tester) == '57eace9bdd2074297cbd2692c1f4cd38'


@pytest.mark.generate_trees
def test_fasttree_multi_param():
    temp_file = TempFile()
    tester = Alb.AlignBuddy(resource("Alignments_cds.phyr"))
    tester = Pb.generate_tree(tester, 'fasttree', '-seed 12345 -wag -fastest -log %s' % temp_file.path)
    assert phylo_to_hash(tester) == '60dc961f90041dd8bc1eab47ffe581c0'


def test_generate_trees_edge_cases():
    temp_file = TempFile()
    tester = Alb.AlignBuddy(resource("Mnemiopsis_cds.nex"))
    with pytest.raises(FileExistsError):
        Pb.generate_tree(tester, "raxml", keep_temp=temp_file.path)

    with pytest.raises(AttributeError):
        Pb.generate_tree(tester, "foo")

    with pytest.raises(FileNotFoundError):
        Pb.generate_tree(tester, "raxml", "-m BINCATLG")

    with pytest.raises(RuntimeError):
        Pb.generate_tree(tester, "fasttree", "-s 12345")

    with pytest.raises(ProcessLookupError):
        backup_path = os.environ["PATH"]
        os.environ["PATH"] = ""
        Pb.generate_tree(tester, "raxml")
        os.environ["PATH"] = backup_path


# ###################### 'li', '--list_ids' ###################### #
li_hashes = ['514675543e958d5177f248708405224d', '229e5d7cd8bb2bfc300fd45ec18e8424', '514675543e958d5177f248708405224d']
li_hashes = [(Pb._make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(li_hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", li_hashes)
def test_list_ids(phylobuddy, next_hash):
    tester = str(Pb.list_ids(phylobuddy))
    assert md5(tester.encode()).hexdigest() == next_hash


# ###################### 'ptr', '--print_trees' ###################### #
hashes = ['1f62f35b8fe64c5eaaf394c272febf4f', 'd5b502ae7e83f0c6ae2211b7d8281b90', '1f62f35b8fe64c5eaaf394c272febf4f',
          '8a5a1f3e409471184da4138ab45fc54e', '370b1c17c1d96ebb866277747f276add', '8a5a1f3e409471184da4138ab45fc54e']
hashes = [(Pb._make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", hashes)
def test_print_trees(phylobuddy, next_hash):
    tester = Pb.trees_to_ascii(phylobuddy)
    assert phylo_to_hash(str(tester)) == next_hash


# ###################### 'pr', '--prune_taxa' ###################### #
pt_hashes = ['99635c6dbf708f94cf4dfdca87113c44', 'fc03b4f100f038277edf6a9f48913dd0', '001db76033cba463a0f187266855e8d5']
pt_hashes = [(Pb._make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(pt_hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", pt_hashes)
def test_prune_taxa(phylobuddy, next_hash):
    Pb.prune_taxa(phylobuddy, 'fir')
    assert phylo_to_hash(phylobuddy) == next_hash


# ######################  'ri', '--rename_ids' ###################### #
ri_hashes = ['6843a620b725a3a0e0940d4352f2036f', '543d2fc90ca1f391312d6b8fe896c59c', '6ce146e635c20ad62e21a1ed6fddbd3a',
             '4dfed97b2a23b8957ee5141bf4681fe4', '77d00fdc512fa09bd1146037d25eafa0', '9b1014be1b38d27f6b7ef73d17003dae']
ri_hashes = [(Pb._make_copy(pb_objects[x]), next_hash) for x, next_hash in enumerate(ri_hashes)]


@pytest.mark.parametrize("phylobuddy, next_hash", ri_hashes)
def test_rename_ids(phylobuddy, next_hash):
    tester = Pb.rename(phylobuddy, 'Mle', 'Phylo')
    assert phylo_to_hash(tester) == next_hash


def test_rename_nodes():
    tester = Pb.PhyloBuddy(resource("tree_with_node_lables.nwk"))
    Pb.rename(tester, "Inner", "Outer")
    assert phylo_to_hash(tester) == "5e02eaa78e25970d7cda0111eef5adba"


# ######################  'su', '--show_unique' ###################### #
def test_show_unique():
    tester = Pb.PhyloBuddy(resource("compare_trees.newick"))
    tester = Pb.show_unique(tester)
    assert phylo_to_hash(tester) == "ea5b0d1fcd7f39cb556c0f5df96281cf"

    with pytest.raises(AssertionError):
        tester = Pb.PhyloBuddy(Pb._make_copy(pb_objects[0]))
        Pb.show_unique(tester)


# ###################### 'sp', '--split_polytomies' ###################### #
def test_split_polytomies():
    tester = Pb.PhyloBuddy('(A,(B,C,D));')
    Pb.split_polytomies(tester)
    assert str(tester) in ['(A:1.0,(B:1.0,(C:1.0,D:1.0)):1.0):1.0;\n', '(A:1.0,(B:1.0,(D:1.0,C:1.0)):1.0):1.0;\n',
                           '(A:1.0,(C:1.0,(B:1.0,D:1.0)):1.0):1.0;\n', '(A:1.0,(C:1.0,(D:1.0,B:1.0)):1.0):1.0;\n',
                           '(A:1.0,(D:1.0,(C:1.0,B:1.0)):1.0):1.0;\n', '(A:1.0,(D:1.0,(B:1.0,C:1.0)):1.0):1.0;\n']


# ################################################# COMMAND LINE UI ################################################## #
# ###################### argparse_init() ###################### #
def test_argparse_init(capsys):
    sys.argv = ['PhyloBuddy.py', resource("compare_trees.newick")]
    temp_in_args, phylobuddy = Pb.argparse_init()
    assert command_line_output_hash(str(phylobuddy)) == "d8e14a2bfc8e9c0ac3c524f5fb478c67"

    sys.argv += ["-f", "foo"]
    with pytest.raises(SystemExit):
        Pb.argparse_init()

    out, err = capsys.readouterr()
    assert "Error: The format 'foo' passed in with the -f flag is not recognized." in err


# ###################### INTERNAL FUNCTIONS ###################### #
def test_print_trees_ui(capsys):
    # Most of this function is covered repeatedly below, so just test the -t flag
    test_in_args = deepcopy(in_args)
    test_in_args.test = True
    test_in_args.list_ids = [True]
    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert err == "*** Test passed ***\n"


def test_in_place_ui(capsys):
    # Some of this function is covered below, so just test an edge
    test_in_args = deepcopy(in_args)
    test_in_args.in_place = True
    test_in_args.trees = [Pb._make_copy(pb_objects[0])]
    test_in_args.screw_formats = "nexus"
    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert "Warning: The -i flag was passed in, but the positional" in err


# ###################### 'ct', '--consensus_tree' ###################### #
def test_consensus_tree_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.consensus_tree = [False]
    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "f20cbd5aae5971cce8efbda15e4e0b7e"

    test_in_args.consensus_tree = [0.9]
    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "447862ed1ed6e98f2fb535ecce70218b"

    test_in_args.consensus_tree = [1.5]
    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "f20cbd5aae5971cce8efbda15e4e0b7e"


# ###################### 'dis', '--distance' ###################### #
def test_distance_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.distance = [False]
    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "7a9f9902ef30ca3d7ac97b63cfdd0b2e"

    test_in_args.distance = ["uwrf"]
    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "6e520f6911aabdf91dfd075d1545bc1e"

    test_in_args.distance = ["ed"]
    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "f0acdd9c751bce77fe14355cc77b69cc"


# ###################### 'gt', '--generate_tree' ###################### #
@pytest.mark.generate_trees
def test_generate_tree_ui1(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.in_format, test_in_args.out_format = "nexus", "newick"
    test_in_args.trees = [resource("Mnemiopsis_cds.nex")]

    test_in_args.generate_tree = [["fasttree", "-seed 12345"]]
    Pb.command_line_ui(test_in_args, [], skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "d7f505182dd1a1744b45cc326096f70c"


@pytest.mark.generate_trees
def test_generate_tree_ui2():
    test_in_args = deepcopy(in_args)
    test_in_args.in_format, test_in_args.out_format = "nexus", "newick"
    test_in_args.trees = [resource("Mnemiopsis_cds.nex")]

    test_in_args.generate_tree = [["fasttree", "-s 12345"]]  # Command doesn't exist in fasttree
    with pytest.raises(SystemExit):
        Pb.command_line_ui(test_in_args, [])


@pytest.mark.generate_trees
def test_generate_tree_ui3():
    test_in_args = deepcopy(in_args)
    test_in_args.in_format, test_in_args.out_format = "nexus", "newick"
    test_in_args.trees = [resource("Mnemiopsis_cds.nex")]

    test_in_args.generate_tree = [["foo"]]
    with pytest.raises(SystemExit):
        Pb.command_line_ui(test_in_args, [])


# ###################### 'li', '--list_ids' ###################### #
def test_list_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.list_ids = [True]

    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "4f0857e0211ff0cd058d0cf7cbaf64d5"

    test_in_args.list_ids = [4]
    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "bf2fbfe1bd52e9b27ae21f5c06e7763a"

    test_in_args.list_ids = [4]
    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy("(, );"), skip_exit=True)
    out, err = capsys.readouterr()
    assert out == "#### tree_1 ####\nNone\n\n"


# ###################### 'ptr', '--print_trees' ###################### #
def test_print_trees_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.print_trees = True

    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "fe340117cb8f573100c00fc897e6c8ce"


# ###################### 'pt', '--prune_taxa' ###################### #
def test_prune_taxa_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.prune_taxa = [["fir"]]

    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "99635c6dbf708f94cf4dfdca87113c44"

    test_in_args.prune_taxa = [["fir", "ovi"]]

    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[1]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "2a385fa95024323fea412fd2b3c3e91f"


# ###################### 'ri', '--rename_ids' ###################### #
def test_rename_ids_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.rename_ids = ['Mle', 'Phylo']

    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "6843a620b725a3a0e0940d4352f2036f"


# ###################### 'sf', '--screw_formats' ###################### #
def test_screw_formats_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = "nexus"

    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "543d2fc90ca1f391312d6b8fe896c59c"


def test_screw_formats_fail(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = "foo"
    Pb.command_line_ui(test_in_args, Pb._make_copy(pb_objects[0]), skip_exit=True)
    foo_out, foo_err = capsys.readouterr()
    assert foo_err == "Error: unknown format 'foo'\n"


def test_screw_formats_inplace_ui(capsys):
    temp_file = TempFile()
    with open(resource("compare_trees.newick"), "r") as ifile:
        temp_file.write(ifile.read())

    test_in_args = deepcopy(in_args)
    test_in_args.screw_formats = "nexus"
    test_in_args.in_place = True
    test_in_args.trees[0] = temp_file.path

    tester = Pb.PhyloBuddy(temp_file.path)
    Pb.command_line_ui(test_in_args, tester, skip_exit=True)
    out, err = capsys.readouterr()
    assert "File over-written at:" in err
    check_file = os.path.isfile("%s.nex" % temp_file.path)
    assert check_file

    test_in_args.trees[0] = "%s.nex" % temp_file.path
    test_in_args.screw_formats = "newick"

    tester = Pb.PhyloBuddy("%s.nex" % temp_file.path)
    Pb.command_line_ui(test_in_args, tester, skip_exit=True)
    nexus_out, nexus_err = capsys.readouterr()
    assert "File over-written at:" in nexus_err
    check_file = os.path.isfile("%s.nwk" % temp_file.path)
    assert check_file


# ###################### 'su', '--show_unique' ###################### #
def test_show_unique_ui(capsys):
    test_in_args = deepcopy(in_args)
    test_in_args.show_unique = True

    Pb.command_line_ui(test_in_args, Pb.PhyloBuddy(resource("compare_trees.newick")), skip_exit=True)
    out, err = capsys.readouterr()
    assert command_line_output_hash(out) == "ea5b0d1fcd7f39cb556c0f5df96281cf"
