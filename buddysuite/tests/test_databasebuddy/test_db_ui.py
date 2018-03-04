import pytest
import os
import re
import sys
import argparse
from copy import deepcopy
from collections import OrderedDict

import buddy_resources as br
import DatabaseBuddy as Db


def fmt(prog):
    return br.CustomHelpFormatter(prog)


parser = argparse.ArgumentParser(prog="DbBuddy.py", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                 description='''
\033[1mDatabaseBuddy\033[m
Go forth to the servers of sequence, and discover.

\033[1mUsage examples\033[m:
DbBuddy.py -ls (launch empty live session)
DbBuddy.py "<accn1,accn2,accn3,...>" -<cmd>
DbBuddy.py "<search term1, search term2,...>" -<cmd>
DbBuddy.py "<accn1,search term1>" -<cmd>
DbBuddy.py "/path/to/file_of_accns" -<cmd>
''')

br.db_modifiers["database"]["choices"] = Db.DATABASES
br.flags(parser, ("user_input", "Specify accession numbers or search terms, "
                                "either in a file or as a comma separated list"),
         br.db_flags, br.db_modifiers, Db.VERSION)

# This is to allow py.test to work with its own flags
in_args = parser.parse_args([])


def mock_cmdloop(*args):
    print(args)
    return True


class MockUsage(object):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def increment(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        return "".join(self.args)


def mock_filter(_, line, mode):
    print("'%s' filter mocked! %s" % (mode, line))
    return


class OpenPermissionError(object):
    def __init__(self, *args, **kwargs):
        pass

    @staticmethod
    def close():
        raise PermissionError


def mock_fileexistserror(*args, **kwargs):
    raise FileExistsError(args, kwargs)


def mock_keyboardinterrupt(*args, **kwargs):
    raise KeyboardInterrupt(args, kwargs)


def mock_guesserror(*args, **kwargs):
    raise br.GuessError("%s, %s" % (args, kwargs))


def mock_systemexit(*args, **kwargs):
    sys.exit("%s, %s" % (args, kwargs))


# A few real accession numbers to test things out with
ACCNS = ["NP_001287575.1", "ADH10263.1", "XP_005165403.2", "A0A087WX72", "A0A096MTH0", "A0A0A9YFB0",
         "XM_003978475", "ENSAMEG00000011912", "ENSCJAG00000008732", "ENSMEUG00000000523"]


# ###################### argparse_init() ###################### #
def test_argparse_init(capsys, monkeypatch, hf):
    monkeypatch.setattr(sys, "argv", ['DatabaseBuddy.py', "Casp9"])
    temp_in_args, dbbuddy = Db.argparse_init()
    assert hf.string2hash(str(dbbuddy)) == "b61a8e0e0a97f33ec1e85c09391ada64"

    monkeypatch.setattr(sys, "argv", ['DatabaseBuddy.py', "Casp9,Panx3", "Cx43"])
    temp_in_args, dbbuddy = Db.argparse_init()
    assert hf.string2hash(str(dbbuddy)) == "c717f3c1636ab03f0c5f5e86d5e909cb"

    monkeypatch.setattr(sys, "argv", ['DatabaseBuddy.py', "-f"])
    with pytest.raises(SystemExit):
        Db.argparse_init()

    out, err = capsys.readouterr()
    assert "DbBuddy.py: error: unrecognized arguments: -f" in err


def test_liveshell_init(monkeypatch, capsys, hf):
    # Default instantiate
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    assert type(liveshell.tmpdir) == br.TempDir
    assert liveshell.terminal_default == "\033[m\033[40m\033[97m"
    assert liveshell.prompt == '[95m[1mDbBuddy[m[40m[97m[1m>[m[40m[97m '
    assert hf.string2hash(liveshell.doc_leader) == "e71aa4976437bdb0c22eeaacfaea6f9f"
    assert hash(liveshell.dbbuddy) == hash(dbbuddy)
    assert liveshell.crash_file == crash_file
    assert os.path.split(liveshell.history_path)[-1] == "cmd_history"
    assert not liveshell.undo
    assert not liveshell.hash
    assert liveshell.shell_execs == []
    assert type(liveshell.usage) == br.Usage
    out, err = capsys.readouterr()
    assert "Your session is currently unpopulated. Use 'search' to retrieve records." in out

    # Set cmd history path
    tmp_dir = br.TempDir()
    monkeypatch.setitem(Db.CONFIG, "data_dir", tmp_dir.path)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    assert liveshell.history_path == "%s%scmd_history" % (tmp_dir.path, os.sep)

    # Permission error
    os.chmod("%s%scmd_history" % (tmp_dir.path, os.sep), 0o333)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    # Windows does not actually change the permissions with os.chmod...
    if os.name != "nt":
        assert liveshell.history_path == "%s%scmd_history" % (liveshell.tmpdir.path, os.sep)

    # Run initial search
    monkeypatch.setattr(Db, "retrieve_summary", lambda _: True)
    dbbuddy = Db.DbBuddy("Inx15")
    Db.LiveShell(dbbuddy, crash_file)
    assert not dbbuddy.records

    Db.LiveShell(dbbuddy, crash_file)
    dbbuddy = Db.DbBuddy("ENSAMEG00000011912")
    assert len(dbbuddy.records) == 1


def test_liveshell_precmd(monkeypatch):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    assert liveshell.precmd("foo bar line") == "foo bar line"


def test_liveshell_postcmd(monkeypatch):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    assert liveshell.postcmd("STOP!", "foo bar line") == "STOP!"
    assert liveshell.usage.stats['LiveShell'][Db.VERSION.short()]['foo'] == 1


def test_liveshell_dump_session(monkeypatch):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy(_databases="uniprot")
    dbbuddy.server_clients["uniprot"] = Db.UniProtRestClient(dbbuddy)
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    pre_dump = liveshell.crash_file.read()
    liveshell.dump_session()
    assert pre_dump == liveshell.crash_file.read()

    liveshell.dbbuddy.search_terms.append("Blahh")
    liveshell.dump_session()
    assert pre_dump != liveshell.crash_file.read()

    assert liveshell.dbbuddy.server_clients['uniprot'].lock
    assert not liveshell.dbbuddy.server_clients['ensembl']
    assert liveshell.undo


def test_liveshell_default(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    liveshell.default("Dunno...")
    out, err = capsys.readouterr()
    assert '*** Unknown syntax: Dunno...\n\n' in out

    with pytest.raises(SystemExit):
        liveshell.default("exit")
    out, err = capsys.readouterr()
    assert "Goodbye" in out


def test_liveshell_append_slash_if_dir():
    tmp_dir = br.TempDir()
    tmp_dir.subfile("test.txt")
    assert Db.LiveShell._append_slash_if_dir(tmp_dir.path) == "%s%s" % (tmp_dir.path, os.sep)
    assert Db.LiveShell._append_slash_if_dir(tmp_dir.subfiles[0]) == tmp_dir.subfiles[0]


def test_liveshell_get_headings(monkeypatch):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db, "retrieve_summary", lambda _: True)
    dbbuddy = Db.DbBuddy(",".join(ACCNS))
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    dbbuddy.records['XM_003978475'].summary = {'organism': 'velociraptor'}
    assert liveshell.get_headings() == ['ACCN', 'DB', 'Type', 'record', 'organism']


def test_liveshell_filter(monkeypatch, hf, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    liveshell.do_load(hf.dbsave)

    # 'keep' (default)
    capsys.readouterr()
    liveshell.filter("(organism) Mouse")
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    with open("temp.del", "w") as ofile:
        assert hf.string2hash(out) == "8c8bc0d638e981c71c41407337bb134d", ofile.write(out)

    # 'restore'
    liveshell.filter("Phaethon", mode='restore')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "375874594a3748a3e13ad713610882c4"

    # 'remove'
    liveshell.filter("Fragment", mode='remove')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "dbd45b5b47c9052112bdc8bc511081b4"

    # Wrong mode
    with pytest.raises(ValueError) as err:
        liveshell.filter("Fragment", mode='Foo')
    assert "The 'mode' argument in filter() must be 'keep', 'remove', or 'restore', not Foo." in str(err)

    # No search string given at all
    monkeypatch.setattr("builtins.input", lambda _: False)
    liveshell.filter(None)
    out, err = capsys.readouterr()
    assert "Error: you must specify a search string.\n" in out

    # No search string given at first
    monkeypatch.setattr("builtins.input", lambda _: "Casein")
    liveshell.filter(None, mode="remove")
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "21cf3142c80dd611b350283b118a14bf"

    monkeypatch.setattr("builtins.input", lambda _: "Apoptosis")
    liveshell.filter(None, mode="restore")
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "cd47b70b98f663c462fe5af5a2e1f729"

    # Multiple terms
    liveshell.filter('"Baculoviral" "Mitogen"', mode='remove')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "f6ca34607e0c5562edfd54135854b793"

    liveshell.filter("'partial' 'Q[0-9]'", mode='remove')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "c4da39b3c33123ed5347a96f9e75995f"

    # Wonkey quotes given as input
    error_msg = "Error: It appears that you are trying to mix quote types (\" and ') while specifying " \
                "multiple filters. Please pick one or the other.\n\n"
    liveshell.filter("'Foo' \"Bar\"", mode='remove')
    out, err = capsys.readouterr()
    assert error_msg in out

    liveshell.filter('"Foo" \'Bar\'', mode='remove')
    out, err = capsys.readouterr()
    assert error_msg in out


def test_liveshell_do_bash(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    capsys.readouterr()
    tmp_file = br.TempFile()
    liveshell.do_bash("echo 'hello from bash' > %s" % tmp_file.path)
    assert "hello from bash" in tmp_file.read()

    monkeypatch.setattr("builtins.input", lambda _: "echo 'Line from input' > %s" % tmp_file.path)
    liveshell.do_bash(None)
    assert "Line from input" in tmp_file.read()

    liveshell.do_bash("cd /this/path/doesnt/exist")
    out, err = capsys.readouterr()
    assert "-sh: cd: /this/path/doesnt/exist: No such file or directory\n" in out

    tmp_dir = br.TempDir()
    tmp_dir.subfile("foo.txt")
    cwd = os.getcwd()
    liveshell.do_bash("cd %s" % tmp_dir.path)
    out, err = capsys.readouterr()
    assert tmp_dir.path in out
    assert os.path.isfile("foo.txt")
    os.chdir(cwd)


def test_liveshell_do_database(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)

    liveshell.do_database("'NcBi_nuc',\t  \"ENSEMBl,uniprot")
    assert sorted(dbbuddy.databases) == ['ensembl', 'ncbi_nuc', 'uniprot']

    liveshell.do_database("ensembl,all")
    assert sorted(dbbuddy.databases) == ["ensembl", "ncbi_nuc", "ncbi_prot", "uniprot"]

    capsys.readouterr()
    liveshell.do_database("Foo ensembl")
    out, err = capsys.readouterr()
    assert "Invalid database choice(s): foo." in out
    assert dbbuddy.databases == ['ensembl']

    liveshell.do_database("Foo")
    out, err = capsys.readouterr()
    assert "Database search list not changed." in out
    assert dbbuddy.databases == ['ensembl']

    monkeypatch.setattr("builtins.input", lambda _: "'ncbi_nuc', 'ensembl'")
    liveshell.do_database(None)
    assert sorted(dbbuddy.databases) == ['ensembl', 'ncbi_nuc']


def test_liveshell_do_delete(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    capsys.readouterr()

    liveshell.do_delete(None)
    out, err = capsys.readouterr()
    assert "The live session is already empty.\n\n" in out

    dbbuddy.records["Foo"] = "Bar"
    liveshell.do_delete("Foo")
    out, err = capsys.readouterr()
    assert "Sorry, I don't understand what you want to delete." in out

    # Delete failures
    liveshell.do_delete("fail")
    out, err = capsys.readouterr()
    assert "Failures list is already empty.\n\n" in out

    dbbuddy.failures["Foo"] = "Bar"
    monkeypatch.setattr(br, "ask", lambda _, **kwargs: False)
    liveshell.do_delete("fail")
    out, err = capsys.readouterr()
    assert "Aborted...\n" in out
    assert len(dbbuddy.failures) == 1

    monkeypatch.setattr(br, "ask", lambda _, **kwargs: True)
    liveshell.do_delete("fail")
    out, err = capsys.readouterr()
    assert "List of failures removed.\n\n" in out
    assert not dbbuddy.failures

    # Delete searches
    liveshell.do_delete("search")
    out, err = capsys.readouterr()
    assert "Search terms list is already empty.\n\n" in out

    dbbuddy.search_terms.append("Bar")
    monkeypatch.setattr(br, "ask", lambda _, **kwargs: False)
    liveshell.do_delete("terms")
    out, err = capsys.readouterr()
    assert "Aborted...\n" in out
    assert len(dbbuddy.search_terms) == 1

    monkeypatch.setattr(br, "ask", lambda _, **kwargs: True)
    liveshell.do_delete("st")
    out, err = capsys.readouterr()
    assert "Search terms removed.\n\n" in out
    assert not dbbuddy.search_terms

    # Delete trash bin
    liveshell.do_delete("trash")
    out, err = capsys.readouterr()
    assert "Trash bin is already empty.\n" in out

    dbbuddy.trash_bin["Foo"] = "Bar"
    monkeypatch.setattr(br, "ask", lambda _, **kwargs: False)
    liveshell.do_delete("trashbin")
    out, err = capsys.readouterr()
    assert "Aborted...\n" in out
    assert len(dbbuddy.trash_bin) == 1

    monkeypatch.setattr(br, "ask", lambda _, **kwargs: True)
    liveshell.do_delete("tb")
    out, err = capsys.readouterr()
    assert "Trash bin emptied.\n\n" in out
    assert not dbbuddy.trash_bin

    # Delete records
    del dbbuddy.records["Foo"]
    dbbuddy.failures["Foo"] = "Bar"
    liveshell.do_delete("records")
    out, err = capsys.readouterr()
    assert "Records list is already empty.\n" in out

    dbbuddy.records["Foo"] = "Bar"
    monkeypatch.setattr(br, "ask", lambda _, **kwargs: False)
    liveshell.do_delete("main")
    out, err = capsys.readouterr()
    assert "Aborted...\n" in out
    assert len(dbbuddy.records) == 1

    monkeypatch.setattr(br, "ask", lambda _, **kwargs: True)
    liveshell.do_delete("recs")
    out, err = capsys.readouterr()
    assert "All records removed from main list (trash bin is still intact).\n\n" in out
    assert not dbbuddy.records

    # Delete everything
    dbbuddy.search_terms.append("Bar")
    dbbuddy.trash_bin["Foo"] = "Bar"
    dbbuddy.records["Foo"] = "Bar"

    monkeypatch.setattr(br, "ask", lambda _, **kwargs: False)
    liveshell.do_delete("")
    out, err = capsys.readouterr()
    assert "Aborted...\n" in out
    assert len(dbbuddy.failures) == 1
    assert len(dbbuddy.search_terms) == 1
    assert len(dbbuddy.trash_bin) == 1
    assert len(dbbuddy.records) == 1

    monkeypatch.setattr(br, "ask", lambda _, **kwargs: True)
    liveshell.do_delete("all")
    out, err = capsys.readouterr()
    assert "Live session cleared of all data.\n\n" in out
    assert not dbbuddy.failures
    assert not dbbuddy.search_terms
    assert not dbbuddy.trash_bin
    assert not dbbuddy.records


def test_liveshell_do_failures(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)

    liveshell.do_failures("Blahh")
    out, err = capsys.readouterr()
    assert "No failures to report\n\n" in out

    dbbuddy.failures["Foo"] = Db.Failure("Bar", "Fake failure")
    liveshell.do_failures()
    out, err = capsys.readouterr()
    assert "The following failures have occured\n" in out
    assert "Bar\nFake failure" in out


def test_liveshell_do_fetch(monkeypatch, capsys):
    def mock_big_record_no_dl(_dbbuddy):
        _dbbuddy.records["NP_001287575.1"] = Db.Record("NP_001287575.1", _size=5000001)

    def mock_big_record_fetch(_dbbuddy):
        _dbbuddy.records["NP_001287575.1"].record = True

    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    monkeypatch.setattr(Db, "retrieve_summary", lambda _: True)

    dbbuddy = Db.DbBuddy("NP_001287575.1")
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)

    monkeypatch.setattr(Db, "retrieve_summary", mock_big_record_no_dl)
    monkeypatch.setattr(br, "ask", lambda _, **kwargs: False)

    liveshell.do_fetch("Foo")
    assert dbbuddy.records["NP_001287575.1"].size == 5000001
    out, err = capsys.readouterr()
    assert "Aborted...\n\n" in out

    monkeypatch.setattr(Db, "retrieve_sequences", mock_big_record_fetch)
    monkeypatch.setattr(br, "ask", lambda _, **kwargs: True)
    liveshell.do_fetch(None)
    out, err = capsys.readouterr()
    print(out)
    assert "Retrieved 5.0 M residues of sequence data\n\n" in out


def test_liveshell_do_format(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)

    monkeypatch.setattr("builtins.input", lambda _: "Foo")
    liveshell.do_format(None)
    out, err = capsys.readouterr()
    assert "'Foo'" in out
    assert "is not a valid format" in out
    assert dbbuddy.out_format == "summary"

    for frmt in Db.FORMATS:
        liveshell.do_format(frmt)
        out, err = capsys.readouterr()
        assert "Output format changed to" in out
        assert frmt in out
        assert dbbuddy.out_format == frmt


def test_liveshell_do_load(monkeypatch, capsys, hf):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    dbbuddy.server_clients["uniprot"] = Db.UniProtRestClient(dbbuddy)
    dbbuddy.server_clients["uniprot"].http_errors_file.write("Hello!")

    monkeypatch.setattr("builtins.input", lambda _: hf.dbsave)
    liveshell.do_load(None)
    out, err = capsys.readouterr()
    assert "Session loaded from file.\n\n" in out

    assert dbbuddy.server_clients["uniprot"].http_errors_file.read() == ""
    headings = liveshell.get_headings()
    for heading in ['ACCN', 'DB', 'Type', 'record', 'entry_name', 'length', 'TaxId', 'organism',
                    'protein_names', 'comments', 'status', 'name', 'biotype',
                    'object_type', 'strand', 'assembly_name']:
        assert heading in headings, print(heading, headings)

    for heading in headings:
        assert heading in ['ACCN', 'DB', 'Type', 'record', 'entry_name', 'length', 'TaxId', 'organism',
                           'protein_names', 'comments', 'status', 'name', 'biotype',
                           'object_type', 'strand', 'assembly_name']

    monkeypatch.setattr("builtins.input", lambda _: "/no/file/here")
    liveshell.do_load(None)
    out, err = capsys.readouterr()
    assert "Error: Unable to read the provided file. Are you sure it's a saved DbBuddy live session?\n\n" in out


def test_liveshell_do_keep(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    monkeypatch.setattr(Db.LiveShell, "filter", mock_filter)
    liveshell.do_keep(None)
    out, err = capsys.readouterr()
    assert "'keep' filter mocked! None" in out


def test_liveshell_do_quit(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)

    dbbuddy.records["Foo"] = "Bar"
    monkeypatch.setattr(br, "ask", lambda _, **kwargs: False)
    liveshell.do_quit(None)
    out, err = capsys.readouterr()
    assert "Aborted...\n\n" in out

    with pytest.raises(SystemExit):
        monkeypatch.setattr(br, "ask", lambda _, **kwargs: True)
        liveshell.do_quit(None)
    out, err = capsys.readouterr()
    assert "Goodbye" in out


def test_liveshell_do_trash(monkeypatch, capsys):
    def mock_show(_, line, mode):
        print("%s show mocked! %s" % (mode, line))
        return

    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    monkeypatch.setattr(Db.LiveShell, "do_show", mock_show)
    liveshell.do_trash(None)
    out, err = capsys.readouterr()
    assert "trash_bin show mocked! None" in out


def test_liveshell_do_remove(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    monkeypatch.setattr(Db.LiveShell, "filter", mock_filter)
    liveshell.do_remove(None)
    out, err = capsys.readouterr()
    assert "'remove' filter mocked! None" in out


def test_liveshell_do_restore(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    monkeypatch.setattr(Db.LiveShell, "filter", mock_filter)
    liveshell.do_restore(None)
    out, err = capsys.readouterr()
    assert "'restore' filter mocked! None" in out


def test_liveshell_do_save(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    dbbuddy.records["foo"] = "bar"

    # Standard, no problems
    tmp_dir = br.TempDir()
    monkeypatch.setattr("builtins.input", lambda _: "%s/save_dir/save_file1" % tmp_dir.path)
    liveshell.do_save(None)
    out, err = capsys.readouterr()
    assert "Live session saved\n\n" in out
    assert os.path.isfile("%s/save_dir/save_file1.db" % tmp_dir.path)
    with open("%s/save_dir/save_file1.db" % tmp_dir.path, "rb") as ifile:
        assert len(ifile.read()) in [279, 281]  # Different versions of python give different file sizes

    # File exists, abort
    monkeypatch.setattr(br, "ask", lambda _, **kwargs: False)
    liveshell.do_save("%s/save_dir/save_file1" % tmp_dir.path)
    out, err = capsys.readouterr()
    assert "Abort...\n\n" in out

    # PermissionError
    monkeypatch.setattr("builtins.open", OpenPermissionError)
    liveshell.do_save("%s/save_dir/save_file2" % tmp_dir.path)
    out, err = capsys.readouterr()
    assert "Error: You do not have write privileges to create a file in the specified directory.\n\n" in out
    assert not os.path.isfile("%s/save_dir/save_file2.db" % tmp_dir.path)

    def makedirs_permissionerror(*args, **kwargs):
        print("makedirs_permissionerror\nargs: %s\nkwargs: %s" % (args, kwargs))
        raise PermissionError

    monkeypatch.setattr(os, "makedirs", makedirs_permissionerror)
    liveshell.do_save("%s/save_dir/deeper_dir/save_file2" % tmp_dir.path)
    out, err = capsys.readouterr()
    assert "Error: You do not have write privileges to create a directory in the specified path.\n\n" in out
    assert not os.path.isfile("%s/save_dir/deeper_dir/save_file2.db" % tmp_dir.path)


def test_liveshell_do_search(monkeypatch):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    dbbuddy.search_terms += ["Foo", "Bar"]

    monkeypatch.setattr("builtins.input", lambda _: "Panx1, Panx2")
    monkeypatch.setattr(Db, "retrieve_summary", lambda _: True)

    liveshell.do_search(None)
    assert dbbuddy.search_terms == ["Foo", "Bar", "Panx1", "Panx2"]
    assert not dbbuddy.records
    assert not dbbuddy.failures

    mock_buddy = Db.DbBuddy("Cx43, Cx32")
    mock_buddy.records["NewRec1"] = True
    mock_buddy.failures["NewFailure1"] = True

    monkeypatch.setattr(Db, "DbBuddy", lambda _: mock_buddy)
    liveshell.do_search("Blahh")

    assert dbbuddy.search_terms == ["Foo", "Bar", "Panx1", "Panx2", "Cx43", "Cx32"]
    assert "NewRec1" in dbbuddy.records
    assert "NewFailure1" in dbbuddy.failures


def test_liveshell_do_show(monkeypatch, capsys, hf):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)

    liveshell.do_load(hf.dbsave)
    capsys.readouterr()

    # Warning message if nothing to show
    liveshell.do_show(None, "trash_bin")
    out, err = capsys.readouterr()
    assert "Nothing in 'trash bin' to show.\n\n" in out

    # Specify columns and number of records
    liveshell.do_show("ACCN organism 3")
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "615c84b691f49f4aeceb92b2dc211ff4"

    # Large group, say 'no' to display
    monkeypatch.setattr(br, "ask", lambda *_, **kwargs: False)
    liveshell.do_show(None)
    out, err = capsys.readouterr()
    assert "Include an integer value with 'show' to return a specific number of records.\n\n" in out

    # Large group, show it anyway
    monkeypatch.setattr(br, "ask", lambda *_, **kwargs: True)
    liveshell.do_show(None)
    out, err = capsys.readouterr()
    # ENSEMBL order gets messed up, so just sort the characters
    assert hf.string2hash(''.join(sorted(out))) == "276b58c3d12682e9375a71bfeb947f8a"

    # Try sequence format on LiveShell with only summary data
    dbbuddy.out_format = "fasta"
    liveshell.do_show(None)
    out, err = capsys.readouterr()
    assert "Warning: only summary data available; there is nothing to display in fasta format." in out

    # Only some records have full sequence data (patch print to true)
    dbbuddy.records["P00520"].record = True
    dbbuddy.records["Q5R454"].record = True
    monkeypatch.setattr(Db.DbBuddy, "print", lambda *_, **kwargs: True)
    liveshell.do_show(None)
    out, err = capsys.readouterr()
    assert "Warning: 1722 records are only summary data, so will not be displayed in fasta format. " \
           "Use 'fetch' to retrieve all sequence data." in err

    # Raise errors
    def mock_length_error(*args, **kwargs):
        print("mock_length_error\nargs: %s\nkwargs: %s" % (args, kwargs))
        raise ValueError("Sequences must all be the same length")
    dbbuddy.out_format = 'nexus'
    monkeypatch.setattr(Db.DbBuddy, "print", mock_length_error)
    liveshell.do_show(None)
    out, err = capsys.readouterr()
    assert "Error: 'nexus' format does not support sequences of different length." in out

    def mock_qual_score_error(*args, **kwargs):
        print("mock_qual_score_error\nargs: %s\nkwargs: %s" % (args, kwargs))
        raise ValueError("No suitable quality scores found in letter_annotations of SeqRecord")
    dbbuddy.out_format = 'fastq'
    monkeypatch.setattr(Db.DbBuddy, "print", mock_qual_score_error)
    liveshell.do_show(None)
    out, err = capsys.readouterr()
    assert "Error: BioPython requires quality scores to output in 'fastq' format, and this data is not " \
           "currently available to DatabaseBuddy." in out

    def mock_valueerror(*args, **kwargs):
        print("mock_valueerror\nargs: %s\nkwargs: %s" % (args, kwargs))
        raise ValueError("Unknown ValueError")
    monkeypatch.setattr(Db.DbBuddy, "print", mock_valueerror)
    with pytest.raises(ValueError) as err:
        liveshell.do_show(None)
    assert "Unknown ValueError" in str(err)


def test_liveshell_do_sort(monkeypatch, capsys, hf):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)

    liveshell.do_load(hf.dbsave)
    capsys.readouterr()

    # Default sort on accession
    start_accns = [x for x in dbbuddy.records]
    liveshell.do_sort(None)
    accns = [x for x in dbbuddy.records]
    assert start_accns != accns
    assert sorted(start_accns) == accns

    # Default sort reversed
    liveshell.do_sort("rev")
    assert [x for x in dbbuddy.records] == sorted(start_accns, reverse=True)
    liveshell.do_sort(None)
    assert [x for x in dbbuddy.records] != sorted(start_accns, reverse=True)
    liveshell.do_sort("reverse")
    assert [x for x in dbbuddy.records] == sorted(start_accns, reverse=True)

    # Sort on column
    _type = [rec.type for accn, rec in dbbuddy.records.items()]
    liveshell.do_sort("Type")
    after_sort = [rec.type for accn, rec in dbbuddy.records.items()]
    assert after_sort != _type
    assert after_sort == sorted(_type)

    database = [rec.database for accn, rec in dbbuddy.records.items()]
    liveshell.do_sort("DB")
    after_sort = [rec.database for accn, rec in dbbuddy.records.items()]
    assert after_sort != database
    assert after_sort == sorted(database)

    organism = [rec.summary['organism'] for accn, rec in dbbuddy.records.items()]
    liveshell.do_sort("organism")
    after_sort = [rec.summary['organism'] for accn, rec in dbbuddy.records.items()]
    assert after_sort != organism
    assert after_sort == sorted(organism)

    length = [rec.summary['length'] for accn, rec in dbbuddy.records.items()]
    liveshell.do_sort("length")
    after_sort = [rec.summary['length'] for accn, rec in dbbuddy.records.items()]
    assert after_sort != length
    assert after_sort == sorted(length)

    protein_names = [rec.summary['protein_names'] for accn, rec in dbbuddy.records.items()
                     if 'protein_names' in rec.summary]
    liveshell.do_sort("protein_names")
    after_sort = [rec.summary['protein_names'] for accn, rec in dbbuddy.records.items()
                  if 'protein_names' in rec.summary]
    assert after_sort != protein_names
    assert after_sort == sorted(protein_names)

    # Sort on multi-column
    dbbuddy.records["A0A0N8ESW5"].record = True
    dbbuddy.records["XP_011997944.1"].record = True
    liveshell.do_sort("record organism length ACCN")
    capsys.readouterr()
    liveshell.do_show("10")
    out, err = capsys.readouterr()
    with open("temp.del", "w") as ofile:
        assert hf.string2hash(out) == "47747589b285aaccb0ee34891d97bd57", ofile.write(out)


def test_liveshell_do_status(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    capsys.readouterr()

    liveshell.do_status(None)
    out, err = capsys.readouterr()
    assert '''\
############################
### DatabaseBuddy object ###
Databases:    ncbi_nuc, ncbi_prot, uniprot, ensembl
Out format:   summary
Searches:     None
Full Recs:    0
Summary Recs: 0
ACCN only:    0
Trash bin:  0
Failures:     0
############################
''' in out


def test_liveshell_do_write(monkeypatch, capsys, hf):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)

    liveshell.do_load(hf.dbsave)
    capsys.readouterr()
    tmp_dir = br.TempDir()

    # write a summary
    monkeypatch.setattr("builtins.input", lambda _: "%s/save1" % tmp_dir.path)
    monkeypatch.setattr(br, "ask", lambda _: True)
    liveshell.do_write(None)
    assert os.path.isfile("%s/save1" % tmp_dir.path)
    with open("%s/save1" % tmp_dir.path, "r") as ifile:
        assert len(ifile.read()) == 290678
    out, err = capsys.readouterr()
    assert re.search("1724 summary records.*written to.*save1", out)

    # write ids/accns
    dbbuddy.out_format = "ids"
    dbbuddy.records['O14727'].record = Db.Record('O14727', _record=True)
    liveshell.do_write("%s/save2" % tmp_dir.path)
    assert os.path.isfile("%s/save2" % tmp_dir.path)
    with open("%s/save2" % tmp_dir.path, "r") as ifile:
        assert len(ifile.read()) == 22719
    out, err = capsys.readouterr()
    assert re.search("1724 accessions.*written to.*save2", out)

    # Abort summary
    monkeypatch.setattr(br, "ask", lambda _: False)
    liveshell.do_write("%s/save3" % tmp_dir.path)
    assert not os.path.isfile("%s/save3" % tmp_dir.path)
    out, err = capsys.readouterr()
    assert "Abort..." in out

    # Permission error
    dbbuddy.out_format = "fasta"
    monkeypatch.setattr("builtins.open", OpenPermissionError)
    liveshell.do_write("%s/save4" % tmp_dir.path)
    assert not os.path.isfile("%s/save4" % tmp_dir.path)
    out, err = capsys.readouterr()
    assert "Error: You do not have write privileges in the specified directory.\n\n" in out

    # File exists
    monkeypatch.setattr(br, "ask", lambda _: False)
    liveshell.do_write("%s/save2" % tmp_dir.path)
    out, err = capsys.readouterr()
    assert "Abort..." in out
    assert "written" not in out

    # Not a directory
    liveshell.do_write("%s/ghostdir/save5" % tmp_dir.path)
    out, err = capsys.readouterr()
    assert "The specified directory does not exist. Please create it before continuing" in out
    assert "written" not in out


def test_liveshell_do_undo(monkeypatch, capsys, hf):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)

    liveshell.do_undo(None)
    out, err = capsys.readouterr()
    assert "There is currently no undo history (only a single undo is possible).\n\n" in out

    liveshell.do_load(hf.dbsave)

    assert not dbbuddy.trash_bin
    liveshell.do_remove("P00520")
    assert dbbuddy.trash_bin
    liveshell.do_undo(None)
    assert not dbbuddy.trash_bin
    out, err = capsys.readouterr()
    assert "Most recent state reloaded\n\n" in out

    liveshell.do_undo(None)
    out, err = capsys.readouterr()
    assert "There is currently no undo history (only a single undo is possible).\n\n" in out


def test_liveshell_complete_bash(monkeypatch):
    # Note, this doesn't work in Windows
    if os.name == "nt":
        return
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    programs = liveshell.complete_bash("wh")
    for program in ['wheel ', 'whereis ', 'whoami ', 'which ', 'who ']:
        assert program in programs


def test_liveshell_complete_database(monkeypatch):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    assert sorted(liveshell.complete_database("")) == ['ensembl', 'ncbi_nuc', 'ncbi_prot', 'uniprot']


def test_liveshell_complete_delete(monkeypatch):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    assert liveshell.complete_delete("") == ["all", "failures", "search", "trash", "records"]


def test_liveshell_complete_format(monkeypatch):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    assert liveshell.complete_format("f") == ['full-summary', 'fasta', 'fastq', 'fastq-sanger',
                                              'fastq-solexa', 'fastq-illumina']


def test_liveshell_complete_keep_remove_resort_trash_show_sort(monkeypatch, hf):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    liveshell.do_load(hf.dbsave)

    # Keep
    assert liveshell.complete_keep("d") == ['(DB) ']
    assert liveshell.complete_keep("len") == ['(length) ']
    assert liveshell.complete_keep('ac') == ['(ACCN) ']

    # Remove
    assert liveshell.complete_remove("d") == ['(DB) ']
    assert liveshell.complete_remove("len") == ['(length) ']
    assert liveshell.complete_remove('ac') == ['(ACCN) ']

    # Restor
    liveshell.do_remove("Human")
    assert liveshell.complete_restore("d") == ['(DB) ']
    assert liveshell.complete_restore("len") == ['(length) ']
    assert liveshell.complete_restore('ac') == ['(ACCN) ']

    # Trash
    assert liveshell.complete_trash("d") == ['DB ']
    assert liveshell.complete_trash("len") == ['length ']
    assert liveshell.complete_trash('ac') == ['ACCN ']

    # Show
    assert liveshell.complete_show("d") == ['DB ']
    assert liveshell.complete_show("len") == ['length ']
    assert liveshell.complete_show('ac') == ['ACCN ']

    # sort
    assert liveshell.complete_sort("d") == ['DB ']
    assert liveshell.complete_sort("len") == ['length ']
    assert liveshell.complete_sort('ac') == ['ACCN ']
    assert liveshell.complete_sort('re') == ['record ', 'reverse ']


def test_liveshell_complete_load_save_write(monkeypatch):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    tmpdir = br.TempDir()
    os.chdir(tmpdir.path)
    tmpdir.subfile("file.txt")
    tmpdir.subdir("extra_dir")

    # Load
    assert liveshell.complete_load("load fi ", "load fi ", 5, 7) == ['file.txt']
    assert liveshell.complete_load("load ", "load ", 5, 5) == ['extra_dir%s' % os.path.sep, 'file.txt']
    assert not liveshell.complete_load("load ", "load ", 4, 5)

    # Save
    assert liveshell.complete_save("save fi ", "save fi ", 5, 7) == ['file.txt']
    assert liveshell.complete_save("save ", "save ", 5, 5) == ['extra_dir%s' % os.path.sep, 'file.txt']
    assert not liveshell.complete_save("save ", "save ", 4, 5)

    # Save
    assert liveshell.complete_write("write fi ", "write fi ", 6, 8) == ['file.txt']
    assert liveshell.complete_write("write ", "write ", 6, 6) == ['extra_dir%s' % os.path.sep, 'file.txt']
    assert not liveshell.complete_write("write ", "write ", 4, 5)


def test_helps(monkeypatch, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)

    liveshell.help_bash()
    out, err = capsys.readouterr()
    assert "Run bash commands" in out

    liveshell.help_database()
    out, err = capsys.readouterr()
    assert "Reset the database" in out

    liveshell.help_delete()
    out, err = capsys.readouterr()
    assert "Remove records completely" in out

    liveshell.help_failures()
    out, err = capsys.readouterr()
    assert "Print the status of" in out

    liveshell.help_fetch()
    out, err = capsys.readouterr()
    assert "Retrieve full records for" in out

    liveshell.help_format()
    out, err = capsys.readouterr()
    assert "Set the output format" in out

    liveshell.help_keep()
    out, err = capsys.readouterr()
    assert "Further refine your results" in out

    liveshell.help_quit()
    out, err = capsys.readouterr()
    assert "End the live session" in out

    liveshell.help_load()
    out, err = capsys.readouterr()
    assert "Recover the contents of a " in out

    liveshell.help_trash()
    out, err = capsys.readouterr()
    assert "Output the records held in" in out

    liveshell.help_remove()
    out, err = capsys.readouterr()
    assert "Further refine your results" in out

    liveshell.help_restore()
    out, err = capsys.readouterr()
    assert "Return a subset of filtered" in out

    liveshell.help_save()
    out, err = capsys.readouterr()
    assert "Save your live session in DB" in out

    liveshell.help_search()
    out, err = capsys.readouterr()
    assert "Search databases (currently set to" in out

    liveshell.help_show()
    out, err = capsys.readouterr()
    assert "Output the records held in" in out

    liveshell.help_sort()
    out, err = capsys.readouterr()
    assert "Alter the order that records" in out

    liveshell.help_status()
    out, err = capsys.readouterr()
    assert "Display the current state of your Live" in out

    liveshell.help_undo()
    out, err = capsys.readouterr()
    assert "Revert the most recent change to your live session." in out

    liveshell.help_write()
    out, err = capsys.readouterr()
    assert "Send records to a file" in out


# ###################### main() ###################### #
def test_main(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["DatabaseBuddy", "Casp9,Panx3", "-ls"])
    monkeypatch.setattr(Db, "command_line_ui", lambda *_: True)
    assert Db.main()

    monkeypatch.setattr(Db, "command_line_ui", mock_guesserror)
    assert not Db.main()

    monkeypatch.setattr(Db, "command_line_ui", mock_systemexit)
    assert not Db.main()

    monkeypatch.setattr(Db, "command_line_ui", mock_fileexistserror)
    monkeypatch.setattr(br, "send_traceback", lambda *_: True)
    assert not Db.main()


# ######################  loose command line ui helpers ###################### #
@pytest.mark.loose
def test_exit(monkeypatch, capsys):
    class MockExitUsage(object):
        @staticmethod
        def increment(*args):
            print(args)
            return True

        @staticmethod
        def save():
            return True

    monkeypatch.setattr(br, "Usage", MockExitUsage)
    monkeypatch.setattr(Db, "LiveShell", lambda *_: True)
    test_in_args = deepcopy(in_args)

    with pytest.raises(SystemExit):
        Db.command_line_ui(test_in_args, Db.DbBuddy())
    out, err = capsys.readouterr()
    assert "('DatabaseBuddy', '%s', 'LiveShell', 0)" % Db.VERSION.short() in out


@pytest.mark.loose
def test_error(monkeypatch, capsys):
    monkeypatch.setattr(Db, "LiveShell", mock_systemexit)

    test_in_args = deepcopy(in_args)
    test_in_args.live_shell = True

    assert Db.command_line_ui(test_in_args, Db.DbBuddy(), skip_exit=True) is None

    monkeypatch.setattr(Db, "LiveShell", mock_fileexistserror)
    monkeypatch.setattr(br.TempFile, "save", lambda *_: True)
    monkeypatch.setattr(br, "send_traceback", lambda *_: True)
    capsys.readouterr()
    assert Db.command_line_ui(test_in_args, Db.DbBuddy(), skip_exit=True) is None
    out, err = capsys.readouterr()
    assert "can be loaded by launching DatabaseBuddy and using the 'load' command." in err


@pytest.mark.loose
def test_retrieve_accessions(monkeypatch):
    # Don't actually run anything, retrieve_summary() is tested elsewhere
    monkeypatch.setattr(Db, "retrieve_summary", lambda *_: True)
    test_in_args = deepcopy(in_args)
    test_in_args.retrieve_accessions = True
    dbbuddy = Db.DbBuddy()
    with pytest.raises(SystemExit):
        Db.command_line_ui(test_in_args, dbbuddy)
    assert dbbuddy.out_format == "ids"

    test_in_args.out_format = "genbank"
    dbbuddy = Db.DbBuddy(_out_format="genbank")
    with pytest.raises(SystemExit):
        Db.command_line_ui(test_in_args, dbbuddy)
    assert dbbuddy.out_format == "ids"

    for out_format in ["ids", "accessions", "full-summary", "summary"]:
        test_in_args.out_format = out_format
        dbbuddy = Db.DbBuddy(_out_format=out_format)
        with pytest.raises(SystemExit):
            Db.command_line_ui(test_in_args, dbbuddy)
        assert dbbuddy.out_format == out_format


@pytest.mark.loose
def test_retrieve_sequences(monkeypatch, capsys, sb_resources, hf):
    # Don't actually run anything, retrieve_summary() is tested elsewhere
    monkeypatch.setattr(Db, "retrieve_summary", lambda *_: True)
    monkeypatch.setattr(Db, "retrieve_sequences", lambda *_: True)
    test_in_args = deepcopy(in_args)
    test_in_args.retrieve_sequences = True

    dbbuddy = Db.DbBuddy()
    dbbuddy.failures = OrderedDict([("Foo", ValueError("Blahh")), ("Bar", AttributeError("Blahh"))])
    dbbuddy.records = sb_resources.get_one("f d").records
    for rec in dbbuddy.records:
        rec.id = re.sub("α", "", rec.id)
        rec.description = re.sub("α", "", rec.description)
    dbbuddy.records = OrderedDict([(rec.id, rec) for rec in dbbuddy.records])
    for accn, rec in dbbuddy.records.items():
        rec.size = len(rec)
        rec.record = rec

    with pytest.raises(SystemExit):
        Db.command_line_ui(test_in_args, dbbuddy)
    out, err = capsys.readouterr()
    assert dbbuddy.out_format == "fasta"
    assert hf.string2hash(out) == "5438dd48c22fd7b6cf2d3da637d333b6"
    assert err == """\
# ###################### Failures ###################### #
Foo
Blahh

Bar
Blahh
# ############################################## #

"""

    monkeypatch.setattr(br, "ask", lambda _, **kwargs: False)
    for accn, rec in dbbuddy.records.items():
        rec.size *= 100
    with pytest.raises(SystemExit):
        Db.command_line_ui(test_in_args, dbbuddy)
    out, err = capsys.readouterr()
    assert out == '\x1b[91mAborted...\n\n\x1b[m\x1b[40m'


@pytest.mark.loose
def test_guess_db(capsys, hf):
    test_in_args = deepcopy(in_args)
    test_in_args.guess_database = True

    with pytest.raises(SystemExit):
        Db.command_line_ui(test_in_args, Db.DbBuddy())

    out, err = capsys.readouterr()
    assert 'Nothing to return' in out

    with pytest.raises(SystemExit):
        Db.command_line_ui(test_in_args, Db.DbBuddy(",".join(ACCNS) + ",Casp9"))

    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "4b3edb0272b02d8e18ce591304fdea1d"
