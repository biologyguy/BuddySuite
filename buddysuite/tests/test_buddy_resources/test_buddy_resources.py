#!/use/bin/env python3
# coding=utf-8

import pytest
import os
import sys
import io
import builtins
import re
import ftplib
import urllib.request
import argparse
import json
from hashlib import md5
from time import sleep
from unittest import mock
from ... import AlignBuddy as Alb
from ... import buddy_resources as br
from pkg_resources import DistributionNotFound
from configparser import ConfigParser
if os.name == "nt":
    import msvcrt

# Globals
TEMP_DIR = br.TempDir()
RESOURCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '..{0}unit_test_resources{0}'.format(os.path.sep))


# Mock resources
class MockLocation(object):
    def __init__(self):
        self.start = 0
        self.end = 1


def mock_valueerror(*args, **kwargs):
    raise ValueError(args, kwargs)


def mock_permissionerror(*args, **kwargs):
    raise PermissionError(args, kwargs)


def mock_unicodeencodeerror(*args, **kwargs):
    raise UnicodeEncodeError(str(args), str(kwargs), 1, 2, "baz")


# Tests
def string2hash(_input):
    return md5(_input.encode("utf-8")).hexdigest()


def test_parse_phylip_format():
    for _format in ["phylip", "phylipis", "phylip-strict", "phylip-interleaved-strict"]:
        assert br.parse_format(_format) == "phylip"

    for _format in ["phylipi", "phylip-relaxed", "phylip-interleaved", "phylipr"]:
        assert br.parse_format(_format) == "phylip-relaxed"

    for _format in ["phylips", "phylipsr", "phylip-sequential", "phylip-sequential-relaxed"]:
        assert br.parse_format(_format) == "phylipsr"

    for _format in ["phylipss", "phylip-sequential-strict"]:
        assert br.parse_format(_format) == "phylipss"

    with pytest.raises(TypeError) as e:
        br.parse_format("foo")
    assert "Format type 'foo' is not recognized/supported" in str(e)


def test_timer():
    timer = br.Timer()
    timer.start()
    sleep(1)
    assert timer.end() == '1 sec'


def test_runtime():
    timer = br.RunTime('x ', ' y', "stdout")
    assert timer.out_type == "stdout"
    assert timer.prefix == "x "
    assert timer.postfix == " y"
    assert not timer.running_process
    timer.start()
    assert timer.running_process
    sleep(1.1)
    timer.end()
    assert not timer.running_process


def test_dynamicprint_init():
    printer = br.DynamicPrint()
    assert printer._last_print == ""
    assert printer._next_print == ""
    assert printer.out_type == sys.stdout
    assert not printer.quiet

    printer = br.DynamicPrint(out_type="stderr")
    assert printer.out_type == sys.stderr


def test_dynamicprint_write(capsys):
    printer = br.DynamicPrint()
    printer.write("Hello")
    printer.new_line(2)
    printer.write("foo")
    printer.clear()
    printer.write("bar")
    out, err = capsys.readouterr()
    assert out == "\r\rHello\n\n\r\rfoo\r   \r\r\rbar"
    assert err == ""

    printer = br.DynamicPrint(out_type="stderr")
    printer.write("Hello")
    out, err = capsys.readouterr()
    assert out == ""
    assert err == "\r\rHello"

    printer = br.DynamicPrint(quiet=True)
    printer.write("Hello")
    out, err = capsys.readouterr()
    assert out == ""
    assert err == ""


def test_pretty_time():
    assert br.pretty_time(1) == '1 sec'
    assert br.pretty_time(60) == '1 min, 0 sec'
    assert br.pretty_time(-1) == '-1 sec'
    assert br.pretty_time(3600) == '1 hrs, 0 min, 0 sec'
    assert br.pretty_time(3740) == '1 hrs, 2 min, 20 sec'
    assert br.pretty_time(100000) == '1 days, 3 hrs, 46 min, 40 sec'


def test_pretty_number():
    # short
    num = 100
    assert br.pretty_number(num, mode='short') == '100'
    num *= 100
    assert br.pretty_number(num, mode='short') == '10.0 K'
    num += 1
    assert br.pretty_number(num, mode='short', precision=3) == '10.001 K'
    assert br.pretty_number(num, mode='short', precision=2) == '10.0 K'
    num -= 1
    num *= 100
    assert br.pretty_number(num, mode='short') == '1.0 M'
    suffixes = ['G', 'T', 'P', 'E', 'Z', 'Y']
    for suffix in suffixes:
        num *= 1000
        assert br.pretty_number(num, mode='short') == '1.0 %s' % suffix

    # medium
    suffixes = ['Kilo', 'Mega', 'Giga', 'Tera', 'Peta', 'Exa', 'Zetta', 'Yotta']
    num = 1
    assert br.pretty_number(num, mode='medium') == '1'
    num *= 1000
    for suffix in suffixes:
        assert br.pretty_number(num, mode='medium') == '1.0 %s' % suffix
        num *= 1000

    # long
    suffixes = ['Thousand', 'Million', 'Billion', 'Trillion', 'Quadrillion', 'Quintillion', 'Sextillion', 'Septillion']
    num = 1
    assert br.pretty_number(num, mode='long') == '1'
    num *= 1000
    for suffix in suffixes:
        assert br.pretty_number(num, mode='long') == '1.0 %s' % suffix
        num *= 1000

    # invalid
    with pytest.raises(ValueError):
        br.pretty_number(100, mode='xxl')


def test_usable_cpu_count(monkeypatch):
    cpu_func = mock.Mock(return_value=10)
    monkeypatch.setattr(br, 'cpu_count', cpu_func)
    assert br.usable_cpu_count() == 7

    cpu_func = mock.Mock(return_value=7)
    monkeypatch.setattr(br, 'cpu_count', cpu_func)
    assert br.usable_cpu_count() == 5

    cpu_func = mock.Mock(return_value=3)
    monkeypatch.setattr(br, 'cpu_count', cpu_func)
    assert br.usable_cpu_count() == 2

    cpu_func = mock.Mock(return_value=1)
    monkeypatch.setattr(br, 'cpu_count', cpu_func)
    assert br.usable_cpu_count() == 1

# skipping run_multicore function for now


def test_run_multicore_function(monkeypatch, hf):
    temp_file = br.TempFile()
    temp_path = temp_file.path
    monkeypatch.setattr(br, "time", lambda: 1)

    monkeypatch.setattr(br, "cpu_count", mock.Mock(return_value=4))
    monkeypatch.setattr(br, "usable_cpu_count", mock.Mock(return_value=4))

    testfunc = mock.Mock()
    testfunc.configure_mock(**{"__name__": "testfunc"})

    nums = range(1, 5)

    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, lambda *_: True, func_args=False,
                                  max_processes=0, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert hf.string2hash(output) in ["107696d60ee9b932ecaffad7c97a609f", "c2adacdaf0de6526c707564068a3460a"]

    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, lambda *_: True, func_args=["Foo"],
                                  max_processes=5, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert hf.string2hash(output) in ["107696d60ee9b932ecaffad7c97a609f", "c2adacdaf0de6526c707564068a3460a"]

    with open(temp_path, "w") as output:
        br.run_multicore_function({"a": 1, "b": 2, "c": 3, "d": 4}, lambda *_: True, func_args=False,
                                  max_processes=-4, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert hf.string2hash(output) in ["cf5afec941a4b854ed78f01d2753009d", "b9a2268fefae3786a611f5e699fd6200"]

    with pytest.raises(AttributeError) as err:
        br.run_multicore_function(nums, lambda *_: True, func_args="Foo", max_processes=4, quiet=False,
                                  out_type=sys.stdout)
    assert "The arguments passed into the multi-thread function must be provided" in str(err)

    class MockTime(object):
        def __init__(self):
            self.timer = 0

        def time(self):
            self.timer += 1 if self.timer < 3 else 0
            return self.timer

    timer = MockTime()
    monkeypatch.setattr(br, "time", timer.time)
    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, lambda *_: True, func_args=False,
                                  max_processes=1, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert hf.string2hash(output) in ["08f48aaa5d64ac48ea15a4aa63d75141", 'b9a2268fefae3786a611f5e699fd6200']

    timer = MockTime()
    monkeypatch.setattr(br, "time", timer.time)
    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, lambda *_: True, func_args=False,
                                  max_processes=0, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert hf.string2hash(output) in ["41cc02db5a989591601308d0657544a8", "c2adacdaf0de6526c707564068a3460a"]


# ######################################  TempDir  ###################################### #
def test_tempdir_init():
    test_dir = br.TempDir()
    assert os.path.exists(test_dir.path)
    assert not len(test_dir.subdirs)
    assert not len(test_dir.subfiles)


def test_tempdir_subdirs():
    test_dir = br.TempDir()
    subdir = test_dir.subdir("test")
    assert os.path.exists(subdir)
    assert os.path.exists("{0}/test".format(test_dir.path))

    with pytest.raises(ValueError):
        test_dir.del_subdir("dasdfs")

    test_dir.del_subdir("test")
    assert not os.path.exists(subdir)
    assert not os.path.exists("{0}/test".format(test_dir.path))

    subdir = test_dir.subdir()
    assert os.path.exists(subdir)
    test_dir.del_subdir(subdir)


def test_tempdir_subfiles():
    test_dir = br.TempDir()
    subfile = test_dir.subfile("testfile")
    assert os.path.exists(subfile)
    assert os.path.exists("{0}/testfile".format(test_dir.path))
    with open(subfile, 'w') as file_to_write:
        file_to_write.write("hello world")
    with open(subfile, 'r') as file_to_write:
        assert file_to_write.read() == "hello world"
    with open("{0}/testfile".format(test_dir.path), 'r') as file_to_write:
        assert file_to_write.read() == "hello world"

    with pytest.raises(ValueError):
        test_dir.del_subfile("dasdfs")

    test_dir.del_subfile("testfile")
    assert not os.path.exists(subfile)
    assert not os.path.exists("{0}/testfile".format(test_dir.path))

    subfile = test_dir.subfile()
    assert os.path.exists(subfile)
    test_dir.del_subfile(subfile)


def test_tempdir_save():
    save_dir = br.TempDir()
    save_dir.subfile("testfile")
    assert save_dir.save("%s/fakedir" % save_dir.path)
    assert os.path.exists("%s/fakedir" % save_dir.path)
    assert os.path.exists("%s/fakedir/testfile" % save_dir.path)
    assert not save_dir.save("%s/fakedir" % save_dir.path)


def test_tempdir_make_dir_generator():
    test_dir = br.TempDir()
    dir_path = test_dir.path
    for _ in test_dir._make_dir():
        pass
    assert not os.path.isdir(dir_path)
    os.makedirs(dir_path)


# ######################################  TempFile  ###################################### #
def test_tempfile():
    test_file = br.TempFile()
    assert os.path.exists(test_file.path)

    test_file.open(mode='w')
    test_file.open(mode='w')

    assert isinstance(test_file.get_handle(), io.TextIOWrapper)

    assert test_file.write("hello world")
    assert test_file.read() == "hello world"
    test_file.close()

    test_file.open("w")
    assert not test_file.write("fail", mode='r')
    test_file.close()

    test_file.clear()
    assert test_file.read() == ""
    test_file.write("hello world")
    test_file.save("{0}/temp".format(TEMP_DIR.path))
    assert os.path.exists("{0}/temp".format(TEMP_DIR.path))
    assert open("{0}/temp".format(TEMP_DIR.path), 'r').read() == "hello world"


def test_safetyvalve():
    valve = br.SafetyValve()
    with pytest.raises(RuntimeError):
        while True:
            valve.step()

    state = 0
    for x in range(20):
        valve.test(state)
        state += 1

    with pytest.raises(RuntimeError):
        while True:
            valve.test(state)


def test_walklevel():
    tmp_dir = br.TempDir()
    tmp_dir.subdir("mydir")
    tmp_dir.subdir("mydir{0}subdir".format(os.path.sep))
    tmp_dir.subdir("mydir{0}subdir{0}subsubdir".format(os.path.sep))

    tmp_dir.subfile("myfile.txt")
    tmp_dir.subfile("mydir{0}subfile.txt".format(os.path.sep))
    tmp_dir.subfile("mydir{0}subdir{0}subsubfile.txt".format(os.path.sep))

    walker = br.walklevel(tmp_dir.path)
    root, dirs, files = next(walker)
    assert root == tmp_dir.path
    assert dirs == ["mydir"]
    assert files == ["myfile.txt"]

    root, dirs, files = next(walker)
    assert root == "%s%smydir" % (tmp_dir.path, os.path.sep)
    assert dirs == ["subdir"]
    assert files == ["subfile.txt"]

    with pytest.raises(StopIteration):
        next(walker)

    walker = br.walklevel(tmp_dir.path)
    counter = 0
    for _ in walker:
        counter += 1
    assert counter == 2

    walker = br.walklevel(tmp_dir.path, level=2)
    counter = 0
    for _ in walker:
        counter += 1
    assert counter == 3


def test_copydir():
    tmp_path = TEMP_DIR.path
    os.makedirs('{0}{1}fakedir'.format(tmp_path, os.path.sep))
    os.makedirs('{0}{1}fakedir{1}fakesub'.format(tmp_path, os.path.sep))
    os.makedirs('{0}{1}fakedir{1}fakesub{1}empty'.format(tmp_path, os.path.sep))
    os.makedirs('{0}{1}fakedir{1}fakesub{1}subsub'.format(tmp_path, os.path.sep))
    open("{0}{1}fakedir{1}fakefile".format(tmp_path, os.path.sep), 'w+').close()
    open("{0}{1}fakedir{1}fakesub{1}fakesubfile".format(tmp_path, os.path.sep), 'w+').close()
    open("{0}{1}fakedir{1}fakesub{1}subsub{1}subsubfile".format(tmp_path, os.path.sep), 'w+').close()

    br.copydir('{0}{1}fakedir'.format(tmp_path, os.path.sep), '{0}{1}fakecopy'.format(tmp_path, os.path.sep))

    for x in os.listdir('{0}{1}fakecopy{1}'.format(tmp_path, os.path.sep)):
        assert x in ["fakefile", "fakesubfile", "subsubfile"]


def test_windows_ask(monkeypatch):
    if os.name == "nt":
        import msvcrt

        class MockMsvcrt(object):
            def __init__(self, outcome=True):
                self.run = outcome
                self.output_chars = ["y", "e", "s", "\r"] if outcome else ["n", "o", "\r"]
                self.output = self._output()

            @staticmethod
            def kbhit():
                return True

            def _output(self):
                for _chr in self.output_chars:
                    yield _chr.encode()

            def getche(self):
                return next(self.output)

        mock_msvcrt = MockMsvcrt(outcome=True)
        monkeypatch.setattr(msvcrt, "getche", mock_msvcrt.getche)
        monkeypatch.setattr(msvcrt, "kbhit", mock_msvcrt.kbhit)
        assert br.ask("test")

        mock_msvcrt = MockMsvcrt(outcome=False)
        monkeypatch.setattr(msvcrt, "getche", mock_msvcrt.getche)
        monkeypatch.setattr(msvcrt, "kbhit", mock_msvcrt.kbhit)
        assert not br.ask("test")

        # Timeout possitive
        monkeypatch.setattr(msvcrt, "kbhit", lambda *_: False)
        assert br.ask("test", timeout=1)
        assert not br.ask("test", timeout=1, default="no")


def test_ask_unix(monkeypatch):
    if os.name != "nt":
        def wait(*args, **kwargs):
            print(args, kwargs)
            sleep(2)
            return 'yes'

        monkeypatch.setattr(builtins, "input", wait)
        assert not br.ask("test", timeout=1)

        fake_input = mock.Mock(return_value="yes")
        monkeypatch.setattr(builtins, "input", fake_input)
        assert br.ask("test")

        fake_input = mock.Mock(return_value="no")
        monkeypatch.setattr(builtins, "input", fake_input)
        assert not br.ask("test")

        fake_input = mock.Mock(return_value="abort")
        monkeypatch.setattr(builtins, "input", fake_input)
        assert not br.ask("test")

        fake_input = mock.Mock(side_effect=["dkjsfaksd", "fsdjgdfgdf", "no"])
        monkeypatch.setattr(builtins, "input", fake_input)
        assert not br.ask("test", timeout=1)

        fake_input = mock.Mock(return_value="")
        monkeypatch.setattr(builtins, "input", fake_input)
        assert br.ask("test", default="yes")

        fake_input = mock.Mock(return_value="")
        monkeypatch.setattr(builtins, "input", fake_input)
        assert not br.ask("test", default="no")


def test_guesserror():
    with pytest.raises(br.GuessError):
        error = br.GuessError("test")
        assert str(error) == "test"
        raise error


def test_phyliperror():
    with pytest.raises(br.PhylipError):
        error = br.PhylipError("test")
        assert str(error) == "test"
        raise error


def test_contributor():
    contributor = br.Contributor("Bud", "Suite", "D", commits=10, github="buddysuite")
    assert contributor.name() == "Bud D Suite"
    assert str(contributor) == "Bud Suite, buddysuite"

    contributor = br.Contributor("Bud", "Suite")
    assert contributor.name() == "Bud Suite"
    assert str(contributor) == "Bud Suite"


def test_customhelpformatter(capsys, hf):
    oldparser = argparse.ArgumentParser(prog="SeqBuddy.py")
    oldparser.add_argument('-v', '--version', help='Show module version #s', action='store')
    oldparser.print_help()
    oldout, olderr = capsys.readouterr()

    parser = argparse.ArgumentParser(prog="SeqBuddy.py", formatter_class=br.CustomHelpFormatter, add_help=False)
    parser.add_argument('-h', '--help', help='show this help message and exit', action='help')
    parser.add_argument('-v', '--version', help='Show module version #s', action='store')
    parser.print_help()
    out, err = capsys.readouterr()
    assert out != oldout
    assert hf.string2hash(oldout) == "5964e3c404d744898171d953c4c6c5f0"
    assert hf.string2hash(out) == "790871bff8c3811c3f7bce6f9e19bb05"
    assert err == "" and olderr == ""


def test_usage(monkeypatch):
    class FakeFTP:
        def __init__(self, *args, **kwargs):
            return

        @staticmethod
        def storlines(*args, **kwargs):
            raise RuntimeError

    config = mock.Mock(return_value={"email": "buddysuite@nih.gov", "diagnostics": True, "user_hash": "ABCDEF",
                                     "data_dir": TEMP_DIR.path})
    monkeypatch.setattr(br, "config_values", config)
    monkeypatch.setattr(br, "FTP", FakeFTP)
    usage = br.Usage()
    usage.stats["last_upload"] = "2015-01-01"
    with pytest.raises(RuntimeError):
        usage.save(send_report=True)

    usage.increment("seqbuddy", "1.3", "usage_test", "10MB")
    usage.increment("seqbuddy", "1.3", "other", "15MB")
    usage.increment("seqbuddy", "1.3", "usage_test", "3MB")
    usage.increment("seqbuddy", "1.4", "usage_test", "5MB")

    usage.save(send_report=False)

    with open(usage.usage_file_path, "r") as usage_file:
        contents = usage_file.read()
        assert "\"seqbuddy\": " in contents
        assert "\"1.4\": " in contents
        assert "\"sizes\": [\"5MB\"]" in contents
        assert "\"usage_test\": 1" in contents
        assert "\"1.3\": " in contents
        assert "\"other\": 1" in contents
        assert "\"sizes\": [\"10MB\", \"15MB\", \"3MB\"]" in contents
        assert "\"usage_test\": 2" in contents
        assert "\"user_hash\": \"ABCDEF\"" in contents

    def raise_ftp_errors(*args, **kwargs):
        print(args, kwargs)
        raise ftplib.error_perm

    # Gracefully handling FTP errors
    FakeFTP.storlines = raise_ftp_errors
    monkeypatch.setattr(br, "FTP", FakeFTP)

    usage = br.Usage()
    usage.stats["last_upload"] = "2015-01-01"
    usage.save(send_report=True)
    with open(usage.usage_file_path, "r") as ifile:
        current_content = ifile.read()
    assert current_content != ""

    monkeypatch.setattr(json, "dump", mock_permissionerror)
    usage.config["diagnostics"] = False
    usage.save()
    with open(usage.usage_file_path, "r") as ifile:
        current_content = ifile.read()
    assert current_content == ""

    monkeypatch.setattr(os.path, "isfile", mock_valueerror)
    usage.config["diagnostics"] = True
    usage = br.Usage()
    assert usage.usage_file_path == usage.tmpfile.path


def test_version():
    contributors = list()
    contributors.append(br.Contributor("Bud", "Suite", "D", commits=10, github="buddysuite"))
    contributors.append(br.Contributor("Sweet", "Water", commits=5, github="sweetwater"))
    version = br.Version("BudddySuite", "3", "5", contributors, release_date={"day": 13, "month": 7, "year": 2016})
    assert version.short() == "3.5"
    assert version.contributors_string() == "Bud D Suite  buddysuite\nSweet Water  sweetwater"
    version_string = re.sub("[\n| ]", "", str(version))
    assert version_string == "BudddySuite3.5(2016-07-13)PublicDomainNoticeThisisfreesoftware;seethesourcefordetailed" \
                             "copyingconditions.ThereisNOwarranty;notevenforMERCHANTABILITYorFITNESSFORAPARTICULAR" \
                             "PURPOSE.Questions/comments/concernscanbedirectedtoSteveBond,steve.bond@nih.gov" \
                             "Contributors:BudDSuitebuddysuiteSweetWatersweetwater"


def test_config_values(monkeypatch):
    fake_config = br.TempFile()
    fake_config.write("[DEFAULT]\nuser_hash = ABCDEFG\ndiagnostics = True\nemail = buddysuite@mockmail.com"
                      "\nshortcuts = /usr/local/sb,/usr/local/alb")
    fake_config.close()
    config_path = fake_config.path

    monkeypatch.setattr(br, "resource_filename", lambda *_: config_path)
    options = br.config_values()
    assert options["user_hash"] == "ABCDEFG"
    assert options["diagnostics"]
    assert options["email"] == "buddysuite@mockmail.com"
    assert options["shortcuts"] == ['/usr/local/sb', '/usr/local/alb']

    def mock_keyerror(*args, **kwargs):
        raise KeyError(args, kwargs)

    monkeypatch.setattr(ConfigParser, "getboolean", mock_keyerror)
    monkeypatch.setattr(ConfigParser, "get", mock_keyerror)
    options = br.config_values()
    assert options["user_hash"] == "hashless"
    assert not options["diagnostics"]
    assert options["email"] == "buddysuite@nih.gov"

    def mock_distributionerror(*args, **kwargs):
        raise DistributionNotFound(args, kwargs)

    monkeypatch.setattr(br, "resource_filename", mock_distributionerror)
    options = br.config_values()
    assert not options["data_dir"]


def test_error_report(monkeypatch):
    class FakeFTP:
        def __init__(self, *args, **kwargs):
            return

        @staticmethod
        def storlines(*args, **kwargs):
            raise RuntimeError  # If a runtime error is raised, the file was "sent"

    fake_error = "ABCD"
    error_hash = b'cb08ca4a7bb5f9683c19133a84872ca7'
    fake_raw_output = io.BytesIO(b'{\"%b\": [1.1, 1.2]}' % error_hash)
    mock_json = mock.Mock(return_value=fake_raw_output)
    monkeypatch.setattr(urllib.request, "urlopen", mock_json)
    config = mock.Mock(return_value={"email": "buddysuite@nih.gov", "diagnostics": True, "user_hash": "hashless",
                                     "data_dir": False})

    monkeypatch.setattr(br, "FTP", FakeFTP)
    monkeypatch.setattr(br, "config_values", config)

    br.error_report(fake_error, "test", "test", br.Version("BuddySuite", 3, 5, _contributors=[]))  # Known bug

    fake_error = "WXYZ"

    fake_raw_output = io.BytesIO(b'{\"%b\": [1.1, 1.2]}' % error_hash)  # Needs to be reset every time
    mock_json = mock.Mock(return_value=fake_raw_output)
    monkeypatch.setattr(urllib.request, "urlopen", mock_json)

    with pytest.raises(RuntimeError):  # Unknown error, diagnostics true
        br.error_report(fake_error, "test", "test", br.Version("BuddySuite", 3, 5, _contributors=[]))

    fake_raw_output = io.BytesIO(b'{\"%b\": [1.1, 1.2]}' % error_hash)
    mock_json = mock.Mock(return_value=fake_raw_output)
    monkeypatch.setattr(urllib.request, "urlopen", mock_json)

    def raise_ftp_errors(*args, **kwargs):
        print(args, kwargs)
        raise ftplib.error_perm

    # Gracefully handling FTP errors
    FakeFTP.storlines = raise_ftp_errors
    monkeypatch.setattr(br, "FTP", FakeFTP)
    br.error_report(fake_error, "test", "test", br.Version("BuddySuite", 3, 5, _contributors=[]))


def test_flags(capsys, hf):
    contributors = list()
    contributors.append(br.Contributor("Bud", "Suite", "D", commits=10, github="buddysuite"))
    contributors.append(br.Contributor("Sweet", "Water", commits=5, github="sweetwater"))
    version = br.Version("BudddySuite", "3", "5", contributors, release_date={"day": 13, "month": 7, "year": 2016})

    parser = argparse.ArgumentParser(prog="SeqBuddy.py", add_help=False)
    pos_dict = ("sequence", "Supply file path(s) or raw sequence. If piping sequences "
                            "into SeqBuddy this argument can be left blank.")
    flag_dict = {"annotate": {"flag": "ano",
                              "nargs": "*",
                              "metavar": "args",
                              "help": "Add a feature (annotation) to selected sequences. "
                                      "Args: <name>, <location (start1-end1,start2-end2...)>, [strand (+|-)], "
                                      "[qualifier (foo=bar) [qualifier]], [regex_pattern [regex_[pattern]]}"},
                 "ave_seq_length": {"flag": "asl",
                                    "action": "append",
                                    "nargs": "?",
                                    "metavar": "'clean'",
                                    "help": "Calculate average sequence length. Specify 'clean' to remove gaps etc "
                                            "first"}}
    mod_dict = {"alpha": {"flag": "a",
                          "action": "store",
                          "help": "If you want the file read with a specific alphabet"},
                "in_format": {"flag": "f",
                              "action": "store",
                              "help": "If SeqBuddy can't guess the file format, try specifying it directly"}}

    br.flags(parser, _positional=pos_dict, _flags=flag_dict, _modifiers=mod_dict, version=version)
    parser.print_help()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "9ea2d4ac842b1712687034ea0abf497b"


def test_parse_format():
    assert br.parse_format("CLUSTAL") == "clustal"
    assert br.parse_format("clustal") == "clustal"
    assert br.parse_format("phylip") == "phylip"
    assert br.parse_format("phylip-interleaved-strict") == "phylip"
    assert br.parse_format("phylipr") == "phylip-relaxed"
    assert br.parse_format("phylips") == "phylipsr"
    assert br.parse_format("phylipss") == "phylipss"

    with pytest.raises(TypeError):
        br.parse_format("buddy")


def test_preparse_flags():
    sys.argv = ['buddy_resources.py', "-v", "-foo", "blahh", "-c", "-ns", "57684", "--blast", "--bar"]
    br.preparse_flags()
    print(sys.argv)
    assert sys.argv == ['buddy_resources.py', '-v', ' -foo', 'blahh', '-c', '-ns', "57684", '--blast', " --bar"]


def test_phylip_sequential_out(alb_resources, sb_resources):
    buddy = alb_resources.get_one("o d n")
    output = br.phylip_sequential_out(buddy)
    assert string2hash(output) == '0379295eb39370bdba17c848ec9a8b73'

    cloned_rec = buddy.alignments[0][3]
    buddy.alignments[0].append(cloned_rec)
    with pytest.raises(br.PhylipError):
        br.phylip_sequential_out(buddy)

    buddy = alb_resources.get_one("o d n")
    buddy = Alb.rename(buddy, "Mle", "M")
    output = br.phylip_sequential_out(buddy, relaxed=False)
    assert string2hash(output) == '830f75901a9e69a91679629613dc0a57'

    buddy = Alb.rename(buddy, "M", "Mleeeeeeeeeeeeeeeee")
    print(buddy.alignments[0])
    with pytest.raises(br.PhylipError):
        br.phylip_sequential_out(buddy, relaxed=False)

    buddy = sb_resources.get_one("d f")
    with pytest.raises(br.PhylipError):
        br.phylip_sequential_out(buddy, _type="seq")


def test_phylip_sequential_read(alb_odd_resources, hf, capsys):
    records = br.phylip_sequential_read(open("{0}Mnemiopsis_cds.physr".format(RESOURCE_PATH),
                                             "r", encoding="utf-8").read())
    buddy = Alb.AlignBuddy(records, out_format="phylipsr")
    assert hf.buddy2hash(buddy) == "c5fb6a5ce437afa1a4004e4f8780ad68", buddy.write("temp.del")

    records = br.phylip_sequential_read(open("{0}Mnemiopsis_cds.physs".format(RESOURCE_PATH),
                                             "r", encoding="utf-8").read(), relaxed=False)
    buddy = Alb.AlignBuddy(records, out_format="phylipss")
    assert hf.buddy2hash(buddy) == "4c0c1c0c63298786e6fb3db1385af4d5"

    with open(alb_odd_resources['dna']['single']['phylipss_cols'], "r", encoding="utf-8") as ifile:
            records = ifile.read()
    with pytest.raises(br.PhylipError) as err:
        br.phylip_sequential_read(records)
    assert "Malformed Phylip --> Less sequence found than expected" in str(err)

    with open(alb_odd_resources['dna']['single']['phylipss_recs'], "r", encoding="utf-8") as ifile:
            records = ifile.read()
    with pytest.raises(br.PhylipError) as err:
        br.phylip_sequential_read(records)
    assert "Malformed Phylip --> 9 sequences expected, 4 found." in str(err)

    capsys.readouterr()

    records = """  3 15
Mle-Panxα4  M--VIE---------A
Mle-Panxα8  M--VLE---------A
Mle-Panxα6  M--LLE----------A
"""
    with pytest.raises(br.PhylipError) as err:
        br.phylip_sequential_read(records)
    assert "Malformed Phylip --> Sequence Mle-Panxα4 has 16 columns, 15 expected." in str(err)

    records = """  3 15
Mle-Panxα4  M--VIE--------A
Mle-Panxα8  M--VLE--------A
Mle-Panxα8  M--LLE--------A
"""
    with pytest.raises(br.PhylipError) as err:
        br.phylip_sequential_read(records)
    assert "Malformed Phylip --> Repeat ID Mle-Panxα8." in str(err)

    records = """  3 15
Mle-Panxα4M--VIE--------A
Mle-Panxα8M--VLE--------A
Mle-Panxα8M--LLE--------A
"""
    with pytest.raises(br.PhylipError) as err:
        br.phylip_sequential_read(records, relaxed=False)
    assert "Malformed Phylip --> Repeat id 'Mle-Panxα8' after strict truncation. " in str(err)


def test_replacements():
    input_str = "This test is A string with numbers (12345) and This [CHARS] is a test"
    assert br.replacements(input_str, "numbers", "integers") == "This test is A string with integers (12345) and " \
                                                                "This [CHARS] is a test"
    assert br.replacements(input_str, "This", "some") == "some test is A string with numbers (12345) and some " \
                                                         "[CHARS] is a test"
    assert br.replacements(input_str, "This", "some", -1) == "This test is A string with numbers (12345) and " \
                                                             "some [CHARS] is a test"
    assert br.replacements(input_str, "This", "some", -2) == "some test is A string with numbers (12345) and " \
                                                             "some [CHARS] is a test"
    assert br.replacements(input_str, "This", "some", -3) == "some test is A string with numbers (12345) and " \
                                                             "some [CHARS] is a test"
    assert br.replacements(input_str, "(?:(?:Thi)|(?:tes))", "foo", -3) == "This foot is A string with numbers " \
                                                                           "(12345) and foos [CHARS] is a foot"
    assert br.replacements(input_str, r'\(([0-9]+)\)', r'\1') == "This test is A string with numbers 12345 and This " \
                                                                 "[CHARS] is a test"
    query = '(.+)\(([0-9]+)\)(.+)\[(CHARS)\](.+)'
    assert br.replacements(input_str, query, r'\1\2\3\4\5') == "This test is A string with numbers 12345 and This " \
                                                               "CHARS is a test"
    assert br.replacements(input_str, '([Tt].{3}).*?(.{3}[Tt])', r'\1\2', -2) == "This test is A strin with numbers " \
                                                                                 "(12345) and This a test"
    with pytest.raises(AttributeError) as err:
        br.replacements(input_str, '(.)(.).{3}', r'\1\2\3')
    assert "There are more replacement match values specified than query parenthesized groups" in str(err)


def test_send_traceback(capsys, monkeypatch):
    donothing = mock.Mock(return_value=0)
    monkeypatch.setattr(br, "error_report", donothing)
    br.send_traceback("test", "test", "RuntimeError\nTraceback (most recent call last)\n\t1 raise "
                                      "RuntimeError(\"Something broke!\"\nRuntimeError: Something broke!", 1.2)
    out, err = capsys.readouterr()
    assert str(out) == "\033[mtest::test has crashed with the following traceback:\033[91m\n\nstr: RuntimeError\n" \
                       "Traceback (most recent call last)\n\t1 raise RuntimeError(\"Something broke!\"\nRuntimeError:" \
                       " Something broke!\n\n\n\n\033[m\n"
    try:
        raise TypeError
    except TypeError as err:
        br.send_traceback("test2", "test2", err, 1.3)
        out, err = capsys.readouterr()
        assert "test2::test2 has crashed with the following traceback:" in out
        assert re.search('File ".*test_buddy_resources.py", line', out)
        assert "raise TypeError" in out


def test_shift_features(sb_resources, hf):
    buddy = sb_resources.get_one("d g")
    buddy.records = [buddy.records[0]]
    features = buddy.records[0].features
    shifted_features = br.shift_features(features, 10, len(buddy.records[0]))
    buddy.records[0].features = shifted_features
    assert hf.buddy2hash(buddy) == "e494bf6dc9c6c73c509e66ffc3db57a9"

    buddy = sb_resources.get_one("d g")
    buddy.records = [buddy.records[0]]
    features = buddy.records[0].features
    shifted_features = br.shift_features(features, -10, len(buddy.records[0]))
    buddy.records[0].features = shifted_features
    assert hf.buddy2hash(buddy) == "c2b852ded3b2829f4aaa7f98f0a9f7f3"

    buddy = sb_resources.get_one("d g")
    buddy.records = [buddy.records[0]]
    features = buddy.records[0].features
    shifted_features = br.shift_features(features, 1160, len(buddy.records[0]))
    buddy.records[0].features = shifted_features
    assert hf.buddy2hash(buddy) == "35e0b3b0079ee41591f5f28ac27f039d"

    buddy = sb_resources.get_one("d g")
    buddy.records = [buddy.records[0]]
    features = buddy.records[0].features
    shifted_features = br.shift_features(features, 1400, len(buddy.records[0]))
    buddy.records[0].features = shifted_features
    assert hf.buddy2hash(buddy) == "750a1dd925dd9cef7e3ec0760fd9de9b"

    buddy = sb_resources.get_one("d g")
    buddy.records = [buddy.records[0]]
    features = buddy.records[0].features
    features[0].location = str("Hello World!")
    with pytest.raises(TypeError):
        br.shift_features(features, 10, len(buddy.records[0]))


def test_ungap_feature_ends_simple(alb_resources):
    rec = alb_resources.get_one("o d g").records()[0]
    feature = rec.features[1]
    feature.location._start = -4
    feature = br.ungap_feature_ends(feature, rec)
    assert feature.location.start == 9

    feature.location._start = 77
    feature.location._end = -4
    feature = br.ungap_feature_ends(feature, rec)
    assert feature.location.start == 9
    assert feature.location.end == 77

    feature.location._end = 238
    feature = br.ungap_feature_ends(feature, rec)
    assert feature.location.end == 225

    feature.location.strand = -1
    feature.location._start = -4
    feature = br.ungap_feature_ends(feature, rec)
    assert feature.location.start == 9

    feature.location._end = 238
    feature = br.ungap_feature_ends(feature, rec)
    assert feature.location.end == 225


def test_ungap_feature_ends_compound(alb_resources):
    rec = alb_resources.get_one("o d g").records()[0]
    feature = rec.features[0]
    feature.location.parts[0]._start = 0
    feature.location.parts[1]._end = 238
    print(feature)
    feature = br.ungap_feature_ends(feature, rec)
    assert feature.location.parts[0].start == 9
    assert feature.location.parts[1].end == 225
    print(feature)


def test_ungap_feature_ends_error(alb_resources):
    rec = alb_resources.get_one("o d g").records()[0]
    feature = rec.features[1]
    feature.location = MockLocation()
    with pytest.raises(TypeError) as err:
        br.ungap_feature_ends(feature, rec)
    assert "FeatureLocation or CompoundLocation object required." in str(err)


def test_old2new_simple(alb_resources, sb_resources):
    # Pulling out Mle-Panxα9 from both buddy objects
    align_rec = alb_resources.get_one("o d g").records()[0]
    seq_rec = sb_resources.get_one("d g").records[0]
    feature = br._old2new(seq_rec.features[1], seq_rec, align_rec)
    assert str(feature.location) == str(align_rec.features[1].location)

    feature.location._start = 45
    feature.location._end = 35
    feature = br._old2new(feature, seq_rec, align_rec)
    assert str(feature.location) == str(align_rec.features[1].location)

    feature.location._start = 0
    feature.location._end = 0
    feature_str = str(feature)
    feature = br._old2new(feature, seq_rec, align_rec)
    assert feature_str == str(feature)

    # Pulling out Mle-Panxα1 to test case where a feature doesn't change between old and new
    align_rec = alb_resources.get_one("o d g").records()[2]
    seq_rec = sb_resources.get_one("d g").records[2]
    feature = seq_rec.features[1]
    feature.location._start = 0
    feature.location._end = 1
    feature_str = str(feature)
    feature = br._old2new(feature, seq_rec, align_rec)
    assert feature_str == str(feature)


def test_old2new_compound(alb_resources, sb_resources):
    align_rec = alb_resources.get_one("o d g").records()[0]
    seq_rec = sb_resources.get_one("d g").records[0]
    feature = br._old2new(seq_rec.features[0], seq_rec, align_rec)
    assert str(feature.location) == str(align_rec.features[0].location)

    part = feature.location.parts[0]
    part._start, part._end = [0, 45]
    feature.location.parts = [part]
    feature = br._old2new(feature, seq_rec, align_rec)
    assert str(feature.location) == str(align_rec.features[0].location.parts[0])

    seq_rec = sb_resources.get_one("d g").records[0]
    feature = seq_rec.features[0]
    feature.location.parts = []
    assert not br._old2new(feature, seq_rec, align_rec)


def test_old2new_error(alb_resources):
    rec = alb_resources.get_one("o d g").records()[0]
    feature = rec.features[1]
    feature.location = MockLocation()
    with pytest.raises(TypeError) as err:
        br._old2new(feature, "seq_rec", "align_rec")
    assert "FeatureLocation or CompoundLocation object required." in str(err)


def test_remap_gapped_features(alb_resources, sb_resources):
    align_recs = alb_resources.get_one("o d g").records()
    seq_recs = sb_resources.get_one("d g").records

    new_recs = br.remap_gapped_features(seq_recs, align_recs)
    align_str, new_str = "", ""
    for align_rec, new_rec in zip(align_recs, new_recs):
        for align_feat, new_feat in zip(align_rec.features, new_rec.features):
            align_str += "%s\n" % align_feat
            new_str += "%s\n" % new_feat
    assert align_str == new_str


def test_stderr(capsys):
    br._stderr("Hello std_err", quiet=False)
    out, err = capsys.readouterr()
    assert err == "Hello std_err"

    br._stderr("Hello std_err", quiet=True)
    out, err = capsys.readouterr()
    assert err == ""


def test_stdout(capsys):
    br._stdout("Hello std_out", quiet=False)
    out, err = capsys.readouterr()
    assert out == "Hello std_out"

    br._stdout("Hello std_out", quiet=True)
    out, err = capsys.readouterr()
    assert out == ""


def test_std_errors(capfd, monkeypatch):
    monkeypatch.setattr(sys.stderr, "write", mock_unicodeencodeerror)
    br._stderr("Hello std_err α", quiet=False)
    out, err = capfd.readouterr()
    assert err == "Hello std_err α"

    monkeypatch.setattr(sys.stdout, "write", mock_unicodeencodeerror)
    br._stdout("Hello std_out α", quiet=False)
    out, err = capfd.readouterr()
    assert out == "Hello std_out α"
