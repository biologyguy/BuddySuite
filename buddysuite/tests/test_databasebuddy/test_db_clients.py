import pytest
from unittest import mock
from urllib.error import HTTPError, URLError

from ... import buddy_resources as br
from ... import DatabaseBuddy as Db

# A few real accession numbers to test things out with
ACCNS = ["NP_001287575.1", "ADH10263.1", "XP_005165403.2", "A0A087WX72", "A0A096MTH0", "A0A0A9YFB0",
         "XM_003978475", "ENSAMEG00000011912", "ENSCJAG00000008732", "ENSMEUG00000000523"]


# ################################################# Database Clients ################################################# #
def test_uniprotrestclient_init():
    dbbuddy = Db.DbBuddy(", ".join(ACCNS[3:6]))
    client = Db.UniProtRestClient(dbbuddy)
    assert hash(dbbuddy) == hash(client.dbbuddy)
    assert client.server == 'http://www.uniprot.org/uniprot'
    assert type(client.temp_dir) == br.TempDir
    assert client.http_errors_file == "%s/errors.txt" % client.temp_dir.path
    assert client.results_file == "%s/results.txt" % client.temp_dir.path
    assert client.max_url == 1000


def test_uniprotrestclient_query_uniprot(monkeypatch, capsys):
    tmp_file = br.TempFile(byte_mode=True)
    tmp_file.write('''A8XEF9
O61786
A0A0H5SBJ0
'''.encode("utf-8"))

    def mock_urlopen_return_handle(*args):
        return tmp_file.get_handle("r")

    def mock_urlopen_raise_httperror(*args):
        raise HTTPError("101", "Fake HTTPError from Mock", "Foo", "Bar", "Baz")

    def mock_urlopen_raise_urlerror(*args):
        raise URLError("Fake URLError from Mock")

    def mock_urlopen_raise_keyboardinterrupt(*args):
        raise KeyboardInterrupt()

    dbbuddy = Db.DbBuddy(", ".join(ACCNS[3:6]))
    client = Db.UniProtRestClient(dbbuddy)
    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_return_handle):
        client.query_uniprot("inx15", {"format": "list"})

    with open(client.results_file, "r") as ifile:
        assert ifile.read() == '''# Search: inx15
A8XEF9
O61786
A0A0H5SBJ0
//
'''
    # Also make sure request_params can come in as a list
    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_return_handle):
        client.query_uniprot("inx15", [{"format": "list"}])

    # Errors
    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_raise_httperror):
        client.query_uniprot("inx15", [{"format": "list"}])
    with open(client.http_errors_file, "r") as ifile:
        assert ifile.read() == "inx15\nHTTP Error Fake HTTPError from Mock: Foo\n//\n"

    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_raise_urlerror):
        client.query_uniprot("inx15", [{"format": "list"}])
    with open(client.http_errors_file, "r") as ifile:
        assert "<urlopen error Fake URLError from Mock>" in ifile.read()

    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_raise_keyboardinterrupt):
        client.query_uniprot("inx15", [{"format": "list"}])
    out, err = capsys.readouterr()
    assert "\n\tUniProt query interrupted by user\n" in err
