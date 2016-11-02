#!/use/bin/env python3
# coding=utf-8

import sys
import re
import os
import shutil
import BuddySuite as Bs
import buddy_resources as br


def test_version(capsys):
    sys.argv = ['BuddySuite.py', "-v"]
    Bs.main()
    out, err = capsys.readouterr()
    assert len(out.strip().split("\n")) == 4
    for buddy in ["SeqBuddy", "AlignBuddy", "PhyloBuddy", "DatabaseBuddy"]:
        assert buddy in out


def test_tools(capsys):
    sys.argv = ['BuddySuite.py', "-t"]
    Bs.main()
    out, err = capsys.readouterr()
    for buddy in ["SeqBuddy", "AlignBuddy", "PhyloBuddy", "DatabaseBuddy"]:
        assert "### %s" % buddy in out


def test_count(capsys):
    sys.argv = ['BuddySuite.py', "-c"]
    Bs.main()
    out, err = capsys.readouterr()
    count = 0
    for buddy in ["SeqBuddy", "AlignBuddy", "PhyloBuddy", "DatabaseBuddy"]:
        num = re.search("%s: ([0-9]+)" % buddy, out)
        count += int(num.group(1))
    num = re.search("Total: ([0-9]+)", out)
    assert count == int(num.group(1))


def test_setup(capsys, monkeypatch):
    tmp_dir = br.TempDir()
    __init__ = tmp_dir.subfile("__init__")
    monkeypatch.setattr(Bs.buddysuite, "__file__", __init__)
    monkeypatch.setattr(br, "ask", lambda *_: True)
    monkeypatch.setattr("builtins.input", lambda *_: "")
    sys.argv = ['BuddySuite.py', "-s"]
    Bs.main()
    assert os.path.isdir("%s/buddy_data" % tmp_dir.path)
    assert os.path.isfile("%s/buddy_data/config.ini" % tmp_dir.path)
    assert os.path.isfile("%s/buddy_data/buddysuite_usage.json" % tmp_dir.path)
    assert os.path.isfile("%s/buddy_data/cmd_history" % tmp_dir.path)
    out, err = capsys.readouterr()
    assert out == """\
[1mWelcome to BuddySuite![m
Let's configure your installation:

[1mProviding a valid email address is recommended if accessing public databases with BuddySuite.
The maintainers of those resources may attempt to contact you before blocking your IP if you are not adhering to their usage limitations.[m

[1mBuddySuite is able to automatically send anonymized usage statistics and crash reports to the developers as part of the software improvement program.[m

[1mSuccess! You're all set.[m
    Email address:     buddysuite@nih.gov
    Send diagnostics:  True

These choices can be changed at any time by re-running setup.
Enjoy the BuddySuite!

"""
    monkeypatch.setattr(br, "ask", lambda *_: False)
    monkeypatch.setattr("builtins.input", lambda *_: "foo@bar.com")
    Bs.main()
    out, err = capsys.readouterr()
    assert out == """\
[1mWelcome to BuddySuite![m
Let's configure your installation:

[1mProviding a valid email address is recommended if accessing public databases with BuddySuite.
The maintainers of those resources may attempt to contact you before blocking your IP if you are not adhering to their usage limitations.[m

[1mBuddySuite is able to automatically send anonymized usage statistics and crash reports to the developers as part of the software improvement program.[m

[1mSuccess! You're all set.[m
    Email address:     foo@bar.com
    Send diagnostics:  False

These choices can be changed at any time by re-running setup.
Enjoy the BuddySuite!

"""
    monkeypatch.setattr("builtins.input", lambda *_: "not an email address")
    Bs.main()
    out, err = capsys.readouterr()
    assert out == """\
[1mWelcome to BuddySuite![m
Let's configure your installation:

[1mProviding a valid email address is recommended if accessing public databases with BuddySuite.
The maintainers of those resources may attempt to contact you before blocking your IP if you are not adhering to their usage limitations.[m

[1mBuddySuite is able to automatically send anonymized usage statistics and crash reports to the developers as part of the software improvement program.[m

[1mSuccess! You're all set.[m
    Email address:     foo@bar.com
    Send diagnostics:  False

These choices can be changed at any time by re-running setup.
Enjoy the BuddySuite!

"""


def test_uninstall(monkeypatch):
    tmp_dir = br.TempDir()
    tmp_dir.subdir("main_site")
    tmp_dir.subdir("main_site/buddysuite")
    __init__ = tmp_dir.subfile("main_site/buddysuite/__init__")
    monkeypatch.setattr(Bs.buddysuite, "__file__", __init__)

    tmp_dir.subfile("seqbuddy")
    tmp_dir.subfile("alignbuddy")
    tmp_dir.subfile("databasebuddy")
    tmp_dir.subfile("phylobuddy")
    tmp_dir.subfile("sb")
    tmp_dir.subfile("something_else")

    monkeypatch.setattr(br, "ask", lambda *_, **__: True)
    monkeypatch.setattr(br, "config_values", lambda *_: {"shortcuts": ["/not/a/link", "%s/sb" % tmp_dir.path]})
    monkeypatch.setattr(shutil, "which", lambda buddy: "%s/%s" % (tmp_dir.path, buddy))
    sys.argv = ['BuddySuite.py', "-u"]

    root, dirs, files = next(br.walklevel(tmp_dir.path))
    assert dirs == ["main_site"]
    assert sorted(files) == ['alignbuddy', 'databasebuddy', 'phylobuddy', 'sb', 'seqbuddy', 'something_else']

    Bs.main()
    root, dirs, files = next(br.walklevel(tmp_dir.path))
    assert dirs == []
    assert files == ["something_else"]

    # Run one more time to ensure that the rmtree FileNotFoundError raises
    Bs.main()
    root, dirs, files = next(br.walklevel(tmp_dir.path))
    assert dirs == []
    assert files == ["something_else"]
