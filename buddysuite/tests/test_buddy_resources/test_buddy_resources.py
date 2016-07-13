#!/use/bin/env python3
# coding=utf-8

import pytest
import os
import io
import builtins
import re
import ftplib
import urllib.request
from hashlib import md5
from time import sleep
from unittest import mock
from ... import AlignBuddy as Alb
from ... import SeqBuddy as Sb
from ... import buddy_resources as br

# Globals
temp_dir = br.TempDir()

def string2hash(_input):
    return md5(_input.encode("utf-8")).hexdigest()


def test_parse_format():
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
        sleep(2)
        timer.end()
    with open(temp_file_path, "r") as temp_file:
        out = temp_file.read()
        print(repr(out))
        assert out == '\n\nx 0 sec y\n         \nx 1 sec y\n         \nx 2 sec y\n'

""" Test works but throws a ValueError (still passes) "I/O operation on closed file."
def test_dynamicprint():
    temp_file_path = temp_dir.subfile("dynamicprint")
    with open(temp_file_path, "w") as temp_file:
        printer = br.DynamicPrint(out_type=temp_file)
        printer.write("Hello")
        printer.write(" World")
        printer.new_line(number=2)
        printer.write("I am")
        printer.clear()
        printer.write("buddysuite")

    with open(temp_file_path, "r") as temp_file:
        out = temp_file.read()
        print(repr(out))
        assert out == '\n\nHello\n     \n World\n\n\n\nI am\n    \n\n\nbuddysuite'
"""


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

    test_dir.subfile("testfile")
    assert test_dir.save("fakedir")
    assert os.path.exists("fakedir")
    assert os.path.exists("fakedir/testfile")
    assert not test_dir.save("fakedir")
    os.remove("fakedir/testfile")
    os.removedirs("fakedir")
    assert not os.path.exists("fakedir/testfile")
    assert not os.path.exists("fakedir")


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

    assert ["fakefile", "fakesubfile", "subsubfile"] == os.listdir('{0}/fakecopy/'.format(tmp_path))


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

# Skipped CustomHelpFormatter
'''
def test_usage(monkeypatch):
    config = mock.Mock(return_value={"email": "buddysuite@nih.gov", "diagnostics": False, "user_hash": "hashless",
                                     "data_dir": False})
    monkeypatch.setattr(br, "config_values", config)
    usage = br.Usage()
'''


def test_version():
    contributors = list()
    contributors.append(br.Contributor("Bud", "Suite", "D", commits=10, github="buddysuite"))
    contributors.append(br.Contributor("Sweet", "Water", commits=5, github="sweetwater"))
    version = br.Version("BudddySuite", "3", "5", contributors)
    print(version.contributors_string())
    print(version)
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

        def storlines(self, *args, **kwargs):
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

# Skip flags


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


def test_phylip_sequential_out():
    buddy = Alb.AlignBuddy("unit_test_resources/Mnemiopsis_cds.nex", in_format="nexus")
    output = br.phylip_sequential_out(buddy)
    assert string2hash(output) == '0379295eb39370bdba17c848ec9a8b73'

    cloned_rec = buddy.alignments[0][3]
    buddy.alignments[0].append(cloned_rec)
    with pytest.raises(br.PhylipError):
        br.phylip_sequential_out(buddy)

    buddy = Alb.AlignBuddy("unit_test_resources/Mnemiopsis_cds.nex", in_format="nexus")
    buddy = Alb.rename(buddy, "Mle", "M")
    output = br.phylip_sequential_out(buddy, relaxed=False)
    assert string2hash(output) == '830f75901a9e69a91679629613dc0a57'

    buddy = Alb.rename(buddy, "M", "Mleeeeeeeeeeeeeeeee")
    print(buddy.alignments[0])
    with pytest.raises(br.PhylipError):
        br.phylip_sequential_out(buddy, relaxed=False)

    buddy = Sb.SeqBuddy("unit_test_resources/Mnemiopsis_cds.fa", in_format="fasta")
    with pytest.raises(br.PhylipError):
        br.phylip_sequential_out(buddy, _type="seq")