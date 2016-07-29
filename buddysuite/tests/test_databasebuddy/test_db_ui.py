import pytest
from unittest import mock
import os

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
