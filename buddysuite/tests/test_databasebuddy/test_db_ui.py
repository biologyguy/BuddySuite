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
    assert liveshell.prompt == '\033[m\033[40m\033[97m\033[1mDbBuddy>\033[m\033[40m\033[97m '
    assert sb_helpers.string2hash(liveshell.doc_leader) == "225570084f034d82e9205b377c75d8d8"
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
    assert sb_helpers.string2hash(out) == "f4ee93430a0d474c75f3ecbf0a67c2ac"

    # 'restore'
    liveshell.filter("Phaethon", mode='restore')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "90245630ece0198d208cd1f89248ef06"

    # 'remove'
    liveshell.filter("Fragment", mode='remove')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "7de5496b941e3fa48018238a8729b162"

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
    assert sb_helpers.string2hash(out) == "4ab080610de6c2d44bd72d3519197dab"

    monkeypatch.setattr("builtins.input", lambda _: "Apoptosis")
    liveshell.filter(None, mode="restore")
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "41d50e2f8b1b0beccd8bc21712e9c986"

    # Multiple terms
    liveshell.filter('"Baculoviral" "Mitogen"', mode='remove')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "0285e4eb25a37199dda47774e79525fe"

    liveshell.filter("'partial' 'Q[0-9]'", mode='remove')
    liveshell.dbbuddy.print()
    out, err = capsys.readouterr()
    assert sb_helpers.string2hash(out) == "611b36bb7f11689e5bf19e8c1bfdfdbc"

    # Wonkey quotes given as input
    error_msg = "Error: It appears that you are trying to mix quote types (\" and ') while specifying " \
                "multiple filters. Please pick one or the other.\n\n"
    liveshell.filter("'Foo' \"Bar\"", mode='remove')
    out, err = capsys.readouterr()
    assert error_msg in out

    liveshell.filter('"Foo" \'Bar\'', mode='remove')
    out, err = capsys.readouterr()
    assert error_msg in out


def test_liveshell_do_load(monkeypatch, sb_resources, capsys):
    monkeypatch.setattr(Db.LiveShell, "cmdloop", mock_cmdloop)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveShell(dbbuddy, crash_file)
    load_file = "%s/mock_resources/test_databasebuddy_clients/dbbuddy_save.db" % sb_resources.res_path
    capsys.readouterr()
    liveshell.do_load(load_file)
    out, err = capsys.readouterr()
    assert "Session loaded from file." in out
    headings = liveshell.get_headings()
    for heading in ['ACCN', 'DB', 'Type', 'record', 'entry_name', 'length', 'organism-id', 'organism',
                    'protein_names', 'comments', 'gi_num', 'TaxId', 'status', 'name', 'biotype',
                    'object_type', 'strand', 'assembly_name', 'name']:
        assert heading in headings

    for heading in headings:
        assert heading in ['ACCN', 'DB', 'Type', 'record', 'entry_name', 'length', 'organism-id', 'organism',
                           'protein_names', 'comments', 'gi_num', 'TaxId', 'status', 'name', 'biotype',
                           'object_type', 'strand', 'assembly_name', 'name']