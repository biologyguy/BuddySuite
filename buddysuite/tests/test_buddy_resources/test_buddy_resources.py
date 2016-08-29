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
from hashlib import md5
from time import sleep
from unittest import mock
from ... import AlignBuddy as Alb
from ... import SeqBuddy as Sb
from ... import buddy_resources as br

# Globals
temp_dir = br.TempDir()
RESOURCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../unit_test_resources')


# Mock resources
class MockLocation(object):
    def __init__(self):
        self.start = 0
        self.end = 1


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
    temp_file_path = temp_dir.subfile("runtime")
    with open(temp_file_path, "w") as temp_file:
        timer = br.RunTime('x ', ' y', temp_file)
        timer.start()
        sleep(3)
        timer.end()
    with open(temp_file_path, "r") as temp_file:
        out = temp_file.read()
        assert re.search('\n\nx 0 sec y\n         \nx 1 sec y\n         \nx 2 sec y\n',
                         out, re.MULTILINE)


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


def test_run_multicore_function(monkeypatch, sb_helpers):
    temp_file = br.TempFile()
    temp_path = temp_file.path
    faketime = mock.Mock(return_value=1)
    monkeypatch.setattr(br, "time", faketime)

    monkeypatch.setattr(br, "cpu_count", mock.Mock(return_value=4))
    monkeypatch.setattr(br, "usable_cpu_count", mock.Mock(return_value=4))

    testfunc = mock.Mock()
    testfunc.configure_mock(**{"__name__": "testfunc"})
    nums = range(1, 21)

    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, testfunc, func_args=False, max_processes=0, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert sb_helpers.string2hash(output) == "5caff0f554558b6a5b972f25be4e3568"

    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, testfunc, func_args=False, max_processes=4, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert sb_helpers.string2hash(output) == "5caff0f554558b6a5b972f25be4e3568"

    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, testfunc, func_args=False, max_processes=4, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert sb_helpers.string2hash(output) == "5caff0f554558b6a5b972f25be4e3568"

    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, testfunc, func_args=False, max_processes=4, quiet=True, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert output == ""

    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, testfunc, func_args=False, max_processes=400, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert sb_helpers.string2hash(output) == "5caff0f554558b6a5b972f25be4e3568"

    with pytest.raises(AttributeError):
        br.run_multicore_function(nums, testfunc, func_args="Test", max_processes=4, quiet=False, out_type=sys.stdout)

    input_args = ["test_string", 4]
    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, testfunc, input_args, max_processes=4, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert sb_helpers.string2hash(output) == "5caff0f554558b6a5b972f25be4e3568"

    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, testfunc, func_args=False, max_processes=2, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert sb_helpers.string2hash(output) == "3f2495fdec6684e0a602579806fc421d"

    nums = [0, 1, 2]
    with open(temp_path, "w") as output:
        br.run_multicore_function(nums, testfunc, func_args=False, max_processes=4, quiet=False, out_type=output)
    with open(temp_path, "r") as out:
        output = out.read()
        assert sb_helpers.string2hash(output) == "b9cc1f9ac03116b3f1727ad2c5ea5a05"


def test_tempdir():
    test_dir = br.TempDir()
    assert os.path.exists(test_dir.path)

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

    save_dir = br.TempDir()
    test_dir.subfile("testfile")
    assert test_dir.save("%s/fakedir" % save_dir.path)
    assert os.path.exists("%s/fakedir" % save_dir.path)
    assert os.path.exists("%s/fakedir/testfile" % save_dir.path)
    assert not test_dir.save("%s/fakedir" % save_dir.path)


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
    test_file.save("{0}/temp".format(temp_dir.path))
    assert os.path.exists("{0}/temp".format(temp_dir.path))
    assert open("{0}/temp".format(temp_dir.path), 'r').read() == "hello world"


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

# Skipped walklevel because idk what it's for


def test_copydir():
    tmp_path = temp_dir.path
    os.makedirs('{0}/fakedir'.format(tmp_path))
    os.makedirs('{0}/fakedir/fakesub'.format(tmp_path))
    os.makedirs('{0}/fakedir/fakesub/empty'.format(tmp_path))
    os.makedirs('{0}/fakedir/fakesub/subsub'.format(tmp_path))
    open("{0}/fakedir/fakefile".format(tmp_path), 'w+').close()
    open("{0}/fakedir/fakesub/fakesubfile".format(tmp_path), 'w+').close()
    open("{0}/fakedir/fakesub/subsub/subsubfile".format(tmp_path), 'w+').close()

    br.copydir('{0}/fakedir'.format(tmp_path), '{0}/fakecopy'.format(tmp_path))

    for x in os.listdir('{0}/fakecopy/'.format(tmp_path)):
        assert x in ["fakefile", "fakesubfile", "subsubfile"]


def test_ask(monkeypatch):

    def wait(*args, **kwargs):
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


def test_customhelpformatter(capsys, sb_helpers):
    oldparser = argparse.ArgumentParser()
    oldparser.add_argument('-v', '--version', help='Show module version #s', action='store')
    oldparser.print_help()
    oldout, olderr = capsys.readouterr()

    parser = argparse.ArgumentParser(formatter_class=br.CustomHelpFormatter, add_help=False)
    parser.add_argument('-h', '--help', help='show this help message and exit', action='help')
    parser.add_argument('-v', '--version', help='Show module version #s', action='store')
    parser.print_help()
    out, err = capsys.readouterr()
    assert out != oldout
    assert sb_helpers.string2hash(oldout) == "1429ec82e7bfffa2189a47404e12e6bd"
    assert sb_helpers.string2hash(out) == "c99e5e1197b0a778a73025d1ab623c1c"
    assert err == "" and olderr == ""

def test_usage(monkeypatch):
    class FakeFTP:
        def __init__(self, *args, **kwargs):
            return

        @staticmethod
        def storlines(*args, **kwargs):
            raise RuntimeError

    config = mock.Mock(return_value={"email": "buddysuite@nih.gov", "diagnostics": True, "user_hash": "ABCDEF",
                                     "data_dir": temp_dir.path})
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
        raise ftplib.error_perm

    # Gracefully handling FTP errors
    FakeFTP.storlines = raise_ftp_errors
    monkeypatch.setattr(br, "FTP", FakeFTP)

    usage = br.Usage()
    usage.stats["last_upload"] = "2015-01-01"
    usage.save(send_report=True)


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
    fake_config.write("[DEFAULT]\nuser_hash = ABCDEFG\ndiagnostics = True\nemail = buddysuite@mockmail.com")
    fake_config.close()
    config_path = fake_config.path

    monkeypatch.setattr(br, "resource_filename", mock.Mock(return_value=config_path))
    options = br.config_values()
    assert options["user_hash"] == "ABCDEFG"
    assert options["diagnostics"]
    assert options["email"] == "buddysuite@mockmail.com"


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
        raise ftplib.error_perm

    # Gracefully handling FTP errors
    FakeFTP.storlines = raise_ftp_errors
    monkeypatch.setattr(br, "FTP", FakeFTP)
    br.error_report(fake_error, "test", "test", br.Version("BuddySuite", 3, 5, _contributors=[]))


def test_flags(capsys, sb_helpers):
    contributors = list()
    contributors.append(br.Contributor("Bud", "Suite", "D", commits=10, github="buddysuite"))
    contributors.append(br.Contributor("Sweet", "Water", commits=5, github="sweetwater"))
    version = br.Version("BudddySuite", "3", "5", contributors, release_date={"day": 13, "month": 7, "year": 2016})

    parser = argparse.ArgumentParser(add_help=False)
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
    assert sb_helpers.string2hash(out) == "08c59420d1a07f528c1d80ed860c511a"


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


def test_phylip_sequential_read(alb_helpers):
    records = br.phylip_sequential_read(open("{0}/Mnemiopsis_cds.physr".format(RESOURCE_PATH), "r").read())
    buddy = Alb.AlignBuddy(records, out_format="phylipsr")
    assert alb_helpers.align2hash(buddy) == "c5fb6a5ce437afa1a4004e4f8780ad68"

    records = br.phylip_sequential_read(open("{0}/Mnemiopsis_cds.physs".format(RESOURCE_PATH), "r").read(),
                                        relaxed=False)
    buddy = Alb.AlignBuddy(records, out_format="phylipss")
    assert alb_helpers.align2hash(buddy) == "4c0c1c0c63298786e6fb3db1385af4d5"


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


def test_shift_features(sb_resources, sb_helpers):
    buddy = sb_resources.get_one("d g")
    buddy.records = [buddy.records[0]]
    features = buddy.records[0].features
    shifted_features = br.shift_features(features, 10, len(buddy.records[0]))
    buddy.records[0].features = shifted_features
    assert sb_helpers.features2hash(buddy) == "2929e27c194fbb4a530023faa602d611"

    buddy = sb_resources.get_one("d g")
    buddy.records = [buddy.records[0]]
    features = buddy.records[0].features
    shifted_features = br.shift_features(features, -10, len(buddy.records[0]))
    buddy.records[0].features = shifted_features
    assert sb_helpers.features2hash(buddy) == "5918b48a9ec783b4010916ec517b66a6"

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


def test_old2new_error(alb_resources, sb_resources):
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
