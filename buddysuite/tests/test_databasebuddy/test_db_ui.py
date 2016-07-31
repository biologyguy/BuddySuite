import pytest
from unittest import mock
import os
import re

from ... import buddy_resources as br
from ... import DatabaseBuddy as Db


def mock_cmdloop(*args):
    return True


class MockUsage():
    def __init__(self, *args):
        self.args = args

    def increment(self, *args, **kwargs):
        self.args = args
        return "".join(self.args)


def mock_filter(_, line, mode):
    print("'%s' filter mocked! %s" % (mode, line))
    return

# A few real accession numbers to test things out with
ACCNS = ["NP_001287575.1", "ADH10263.1", "XP_005165403.2", "A0A087WX72", "A0A096MTH0", "A0A0A9YFB0",
         "XM_003978475", "ENSAMEG00000011912", "ENSCJAG00000008732", "ENSMEUG00000000523"]


def test_liveshell_init(monkeypatch, capsys, sb_helpers):
    # Default instantiate
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    assert type(liveshell.tmpdir) == br.TempDir
    assert liveshell.terminal_default == "\033[m\033[40m\033[97m"
    assert liveshell.prompt == '[95m[1mDbBuddy[m[40m[97m[1m>[m[40m[97m '
    assert sb_helpers.string2hash(liveshell.doc_leader) == "e71aa4976437bdb0c22eeaacfaea6f9f"
    assert hash(liveshell.dbbuddy) == hash(dbbuddy)
    assert liveshell.crash_file == crash_file
    assert liveshell.history_path.split("/")[-1] == "cmd_history"
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
    assert liveshell.history_path == "%s/cmd_history" % tmp_dir.path

    # Permission error
    os.chmod("%s/cmd_history" % tmp_dir.path, 0o333)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    assert liveshell.history_path == "%s/cmd_history" % liveshell.tmpdir.path

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
    assert liveshell.usage.stats['LiveShell'] == {'1.0': {'foo': 1}}


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
    assert Db.LiveShell._append_slash_if_dir(tmp_dir.path) == "%s/" % tmp_dir.path
    assert Db.LiveShell._append_slash_if_dir(tmp_dir.subfiles[0]) == tmp_dir.subfiles[0]


def test_liveshell_get_headings(monkeypatch):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db, "retrieve_summary", lambda _: True)
    dbbuddy = Db.DbBuddy(",".join(ACCNS))
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    dbbuddy.records['XM_003978475'].summary = {'organism': 'velociraptor'}
    assert liveshell.get_headings() == ['ACCN', 'DB', 'Type', 'record', 'organism']


def test_liveshell_filter(monkeypatch, sb_resources, sb_helpers, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    load_file = "%s/mock_resources/test_databasebuddy_clients/dbbuddy_save.db" % sb_resources.res_path
    liveshell.do_load(load_file)

    # 'keep' (default)
    capsys.readouterr()
    liveshell.filter("(organism) Mouse")
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "9774790626857cd05298b4e9c5e09836"

    # 'restore'
    liveshell.filter("Phaethon", mode='restore')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "836e1b6810b2e349634face7b19d4999"

    # 'remove'
    liveshell.filter("Fragment", mode='remove')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "746d5e86ff1d3b23707977e0e41fd210"

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
    assert sb_helpers.string2hash(out) == "fdcfcc6d32d726cba592e5c9d0bfdf44"

    monkeypatch.setattr("builtins.input", lambda _: "Apoptosis")
    liveshell.filter(None, mode="restore")
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "a3249f5616e3ec863d911638e7f82ed8"

    # Multiple terms
    liveshell.filter('"Baculoviral" "Mitogen"', mode='remove')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "ef0ef9f16687530cadea9a465ff92634"

    liveshell.filter("'partial' 'Q[0-9]'", mode='remove')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "4aa2b9aaf54bbcb874e17621da1a43c5"

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
    assert tmp_file.read() == "hello from bash\n"

    monkeypatch.setattr("builtins.input", lambda _: "echo 'Line from input' > %s" % tmp_file.path)
    liveshell.do_bash(None)
    assert tmp_file.read() == "Line from input\n"

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
    assert dbbuddy.databases == ['ncbi_nuc', 'ensembl', 'uniprot']

    liveshell.do_database("ensembl,all")
    assert dbbuddy.databases == ["ncbi_nuc", "ncbi_prot", "uniprot", "ensembl"]

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
    assert dbbuddy.databases == ['ncbi_nuc', 'ensembl']


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

    for fmt in Db.FORMATS:
        liveshell.do_format(fmt)
        out, err = capsys.readouterr()
        assert "Output format changed to" in out
        assert fmt in out
        assert dbbuddy.out_format == fmt


def test_liveshell_do_load(monkeypatch, capsys, sb_resources):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(Db.LiveShell, "dump_session", lambda _: True)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    dbbuddy.server_clients["uniprot"] = Db.UniProtRestClient(dbbuddy)
    dbbuddy.server_clients["uniprot"].http_errors_file.write("Hello!")

    db_session = "%s/mock_resources/test_databasebuddy_clients/dbbuddy_save.db" % sb_resources.res_path
    monkeypatch.setattr("builtins.input", lambda _: db_session)
    liveshell.do_load(None)
    out, err = capsys.readouterr()
    assert "Session loaded from file.\n\n" in out

    assert dbbuddy.server_clients["uniprot"].http_errors_file.read() == ""
    headings = liveshell.get_headings()
    for heading in ['ACCN', 'DB', 'Type', 'record', 'entry_name', 'length', 'organism-id', 'organism',
                    'protein_names', 'comments', 'gi_num', 'TaxId', 'status', 'name', 'biotype',
                    'object_type', 'strand', 'assembly_name', 'name']:
        assert heading in headings

    for heading in headings:
        assert heading in ['ACCN', 'DB', 'Type', 'record', 'entry_name', 'length', 'organism-id', 'organism',
                           'protein_names', 'comments', 'gi_num', 'TaxId', 'status', 'name', 'biotype',
                           'object_type', 'strand', 'assembly_name', 'name']

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
        assert len(ifile.read()) == 290

    # File exists, abort
    monkeypatch.setattr(br, "ask", lambda _, **kwargs: False)
    liveshell.do_save("%s/save_dir/save_file1" % tmp_dir.path)
    out, err = capsys.readouterr()
    assert "Abort...\n\n" in out

    # Permission errors
    class OpenPermissionError():
        def __init__(self, *args):
            pass

        @staticmethod
        def close():
            raise PermissionError

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


