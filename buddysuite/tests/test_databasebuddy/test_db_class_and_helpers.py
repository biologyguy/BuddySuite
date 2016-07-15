""" tests basic functionality of DatabaseBuddy class """
import pytest
import datetime

from ... import buddy_resources as br
from ... import DatabaseBuddy as Db


# ##################################################### GLOBALS ###################################################### #
def test_globals():
    assert Db.TRASH_SYNOS == ["t", "tb", "t_bin", "tbin", "trash", "trashbin", "trash-bin", "trash_bin"]
    assert Db.RECORD_SYNOS == ["r", "rec", "recs", "records", "main", "filtered"]
    assert Db.SEARCH_SYNOS == ["st", "search", "search-terms", "search_terms", "terms"]
    assert Db.DATABASES == ["ncbi_nuc", "ncbi_prot", "uniprot", "ensembl"]
    assert Db.RETRIEVAL_TYPES == ["protein", "nucleotide", "gi_num"]
    assert Db.FORMATS == ["ids", "accessions", "summary", "full-summary", "clustal", "embl", "fasta", "fastq",
                          "fastq-sanger", "fastq-solexa", "fastq-illumina", "genbank", "gb", "imgt", "nexus", "phd",
                          "phylip", "seqxml", "sff", "stockholm", "tab", "qual"]
    assert sorted(list(Db.CONFIG)) == ['data_dir', 'diagnostics', 'email', 'user_hash']
    assert type(Db.VERSION) == br.Version
    assert str(Db.VERSION) == """DatabaseBuddy 1.0 (%s)

Public Domain Notice
This is free software; see the source for detailed copying conditions.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
Questions/comments/concerns can be directed to Steve Bond, steve.bond@nih.gov

Contributors:
Stephen Bond    https://github.com/biologyguy
Karl Keat       https://github.com/KarlKeat
Jeremy Labarge  https://github.com/biojerm
""" % datetime.date.today()
    assert Db.GREY == "\033[90m"
    assert Db.RED == "\033[91m"
    assert Db.GREEN == "\033[92m"
    assert Db.YELLOW == "\033[93m"
    assert Db.BLUE == "\033[94m"
    assert Db.MAGENTA == "\033[95m"
    assert Db.CYAN == "\033[96m"
    assert Db.WHITE == "\033[97m"
    assert Db.BOLD == "\033[1m"
    assert Db.UNDERLINE == "\033[4m"
    assert Db.NO_UNDERLINE == "\033[24m"
    assert Db.DEF_FONT == "\033[39m"


# ################################################# HELPER FUNCTIONS ################################################# #
def test_databaseerror():
    with pytest.raises(Db.DatabaseError) as err:
        raise Db.DatabaseError("Raised and err")
    assert "Raised and err" in str(err)


def test_stderr(capsys):
    Db._stderr("Hello from _stderr()")
    out, err = capsys.readouterr()
    assert err == "Hello from _stderr()"

    Db._stderr("Hello from _stderr()", quiet=True)
    out, err = capsys.readouterr()
    assert err == ""


def test_stdout(capsys):
    Db._stdout("Hello from _stdout()")
    out, err = capsys.readouterr()
    assert out == "Hello from _stdout()\033[m"

    Db._stdout("Hello from _stdout()", quiet=True)
    out, err = capsys.readouterr()
    assert out == ""

    Db._stdout("Hello from _stdout()", format_in=Db.RED)
    out, err = capsys.readouterr()
    assert out == "\x1b[91mHello from _stdout()\033[m"

    Db._stdout("Hello from _stdout()", format_in=[Db.UNDERLINE, Db.RED])
    out, err = capsys.readouterr()
    assert out == "\033[4m\x1b[91mHello from _stdout()\033[m"

    Db._stdout("Hello from _stdout()", format_out=Db.RED)
    out, err = capsys.readouterr()
    assert out == "Hello from _stdout()\x1b[91m"

    Db._stdout("Hello from _stdout()", format_out=[Db.UNDERLINE, Db.RED])
    out, err = capsys.readouterr()
    assert out == "Hello from _stdout()\033[4m\x1b[91m"

    with pytest.raises(AttributeError) as err:
        Db._stdout("Hello from _stdout()", format_in=["\033[4m", "\x1b[91"])
    assert "Malformed format_in attribute escape code" in str(err)

    with pytest.raises(AttributeError) as err:
        Db._stdout("Hello from _stdout()", format_out=["\033[4m", "\x1b[91"])
    assert "Malformed format_out attribute escape code" in str(err)


def test_terminal_colors():
    term_colors = Db.terminal_colors()
    for color in [Db.MAGENTA, Db.CYAN, Db.GREEN, Db.RED, Db.YELLOW, Db.GREY, Db.MAGENTA, Db.CYAN]:
        assert next(term_colors) == color


def test_check_database(capsys):
    assert Db.check_database() == Db.DATABASES
    assert Db.check_database('all') == Db.DATABASES
    assert Db.check_database('ALL') == Db.DATABASES
    assert Db.check_database('Foo') == Db.DATABASES
    out, err = capsys.readouterr()
    assert "Warning: 'foo' is not a valid database choice, omitted.\n" in err
    assert "Warning: No valid database choice provided. Setting to default 'all'.\n" in err
    assert Db.check_database("ncbi_prot") == ["ncbi_prot"]
    assert Db.check_database(["foo", "uniprot", "ncbi_prot"]) == ["uniprot", "ncbi_prot"]
    out, err = capsys.readouterr()
    assert "Warning: 'foo' is not a valid database choice, omitted.\n" in err


def test_check_type_protein():
    for _input in ["p", "pr", "prt", "prtn", "prn", "prot", "protn", "protien", "protein"]:
        assert Db.check_type(_input) == "protein"


def test_check_type_nuc():
    for _input in ["n", "ncl", "nuc", "dna", "nt", "gene", "transcript", "nucleotide"]:
        assert Db.check_type(_input) == "nucleotide"


def test_check_type_gi():
    for _input in ["g", "gi", "gn", "gin", "gi_num", "ginum", "gi_number"]:
        assert Db.check_type(_input) == "gi_num"


def test_check_type_default(capsys):
    assert Db.check_type("foo") == "protein"
    out, err = capsys.readouterr()
    assert err == "Warning: 'foo' is not a valid choice for '_type'. Setting to default 'protein'.\n"
