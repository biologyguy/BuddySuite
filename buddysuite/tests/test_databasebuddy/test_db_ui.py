import pytest

from ... import buddy_resources as br
from ... import DatabaseBuddy as Db


class MockUsage():
    def __init__(self, *args):
        self.args = args

    def increment(self, *args, **kwargs):
        self.args = args
        return "".join(self.args)


def mock_cmdloop(*args):
    return True


def test_livesearch_init(monkeypatch, capsys, sb_helpers):
    monkeypatch.setattr(Db.LiveSearch, "cmdloop", mock_cmdloop)
    monkeypatch.setattr(br, "Usage", MockUsage)
    dbbuddy = Db.DbBuddy()
    crash_file = br.TempFile(byte_mode=True)
    liveshell = Db.LiveSearch(dbbuddy, crash_file)
    assert type(liveshell.tmpdir) == br.TempDir
    assert liveshell.terminal_default == "\033[m\033[40m\033[97m"
    assert liveshell.prompt == '\033[m\033[40m\033[97m\033[1mDbBuddy>\033[m\033[40m\033[97m '
    assert sb_helpers.string2hash(liveshell.doc_leader) == "225570084f034d82e9205b377c75d8d8"
    assert hash(liveshell.dbbuddy) == hash(dbbuddy)
    assert liveshell.crash_file == crash_file
    assert liveshell.history_path.split("/")[-1] == "cmd_history"
    assert not liveshell.undo
    assert not liveshell.hash
    assert not liveshell.shell_execs
    #assert liveshell.usage.increment() == ""
