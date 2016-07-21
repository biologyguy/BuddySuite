import pytest
from unittest import mock
from urllib.error import HTTPError, URLError
from collections import OrderedDict

from ... import buddy_resources as br
from ... import DatabaseBuddy as Db

# A few real accession numbers to test things out with
ACCNS = ["NP_001287575.1", "ADH10263.1", "XP_005165403.2", "A0A087WX72", "A0A096MTH0", "A0A0A9YFB0",
         "XM_003978475", "ENSAMEG00000011912", "ENSCJAG00000008732", "ENSMEUG00000000523"]


# Mock functions and classes
def mock_urlopen_handle_uniprot_ids(*args):
    tmp_file = br.TempFile(byte_mode=True)
    tmp_file.write('''A8XEF9
O61786
A0A0H5SBJ0
'''.encode("utf-8"))
    return tmp_file.get_handle("r")


def mock_urlopen_uniprot_count_hits(*args):
    tmp_file = br.TempFile(byte_mode=True)
    tmp_file.write('''# Search: (inx15)+OR+(inx16)
O61787
A0A0V1AZ11
A8XEF9
A8XEF8
A0A0B2VB60
A0A0V0W5E2
O61786
A0A0H5SBJ0
E3MGD6
E3MGD5
//'''.encode("utf-8"))
    return tmp_file.get_handle("r")


def mock_urlopen_uniprot_summary(*args):
    tmp_file = br.TempFile(byte_mode=True)
    tmp_file.write('''# Search: inx15
A8XEF9	A8XEF9_CAEBR	381	6238	Caenorhabditis briggsae	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
O61786	O61786_CAEEL	382	6239	Caenorhabditis elegans	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
A0A0H5SBJ0	A0A0H5SBJ0_BRUMA	129	6279	Brugia malayi (Filarial nematode worm)	Innexin	Function (1); Sequence similarities (1); Subcellular location (1)
E3MGD6	E3MGD6_CAERE	384	31234	Caenorhabditis remanei (Caenorhabditis vulgaris)	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
//
# Search: inx16
O61787	INX16_CAEEL	372	6239	Caenorhabditis elegans	Innexin-16 (Protein opu-16)	Function (1); Sequence similarities (1); Subcellular location (1)
A0A0V1AZ11	A0A0V1AZ11_TRISP	406	6334	Trichinella spiralis (Trichina worm)	Innexin	Caution (1); Function (1); Sequence similarities (1); Subcellular location (2)
A8XEF8	A8XEF8_CAEBR	374	6238	Caenorhabditis briggsae	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
A0A0B2VB60	A0A0B2VB60_TOXCA	366	6265	Toxocara canis (Canine roundworm)	Innexin	Caution (2); Function (1); Sequence similarities (1); Subcellular location (1)
A0A0V0W5E2	A0A0V0W5E2_9BILA	410	92179	Trichinella sp. T6	Innexin	Caution (2); Function (1); Sequence similarities (1); Subcellular location (1)
//'''.encode("utf-8"))
    return tmp_file.get_handle("r")


def mock_urlopen_raise_httperror(*args):
    raise HTTPError("101", "Fake HTTPError from Mock", "Foo", "Bar", "Baz")


def mock_urlopen_raise_urlerror(*args):
    raise URLError("Fake URLError from Mock")


def mock_urlopen_raise_keyboardinterrupt(*args):
    raise KeyboardInterrupt()


# ################################################# Database Clients ################################################# #
# Generic
def test_client_init():
    dbbuddy = Db.DbBuddy(", ".join(ACCNS[3:6]))
    client = Db.GenericClient(dbbuddy)
    assert hash(dbbuddy) == hash(client.dbbuddy)
    assert type(client.http_errors_file) == br.TempFile
    assert type(client.results_file) == br.TempFile
    assert client.max_url == 1000
    with client.lock:
        assert True


def test_client_parse_error_file():
    dbbuddy = Db.DbBuddy()
    client = Db.GenericClient(dbbuddy)

    assert not client._parse_error_file()
    assert not dbbuddy.failures
    client.http_errors_file.write("Casp9\n%s\n//\n" % HTTPError("101", "Fake HTTPError from Mock", "Foo", "Bar", "Baz"))
    client.http_errors_file.write("Inx1\n%s\n//\n" % URLError("Fake URLError from Mock"))

    assert client._parse_error_file() == '''Casp9
HTTP Error Fake HTTPError from Mock: Foo

Inx1
<urlopen error Fake URLError from Mock>

'''
    assert len(dbbuddy.failures) == 2

    # Repeat to make sure that the same error is not added again
    client.http_errors_file.write("Inx1\n%s\n//\n" % URLError("Fake URLError from Mock"))
    assert not client._parse_error_file()
    assert len(dbbuddy.failures) == 2


# UniProt
def test_uniprotrestclient_init():
    dbbuddy = Db.DbBuddy(", ".join(ACCNS[3:6]))
    client = Db.UniProtRestClient(dbbuddy)
    assert hash(dbbuddy) == hash(client.dbbuddy)
    assert client.server == 'http://www.uniprot.org/uniprot'
    assert type(client.http_errors_file) == br.TempFile
    assert type(client.results_file) == br.TempFile
    assert client.max_url == 1000


def test_uniprotrestclient_query_uniprot(capsys):
    dbbuddy = Db.DbBuddy()
    client = Db.UniProtRestClient(dbbuddy)
    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_handle_uniprot_ids):
        client.query_uniprot("inx15", {"format": "list"})

    assert client.results_file.read() == '''# Search: inx15
A8XEF9
O61786
A0A0H5SBJ0
//
'''
    # Also make sure request_params can come in as a list
    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_handle_uniprot_ids):
        client.query_uniprot("inx15", [{"format": "list"}])

    # Errors
    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_raise_httperror):
        client.query_uniprot("inx15", [{"format": "list"}])
    assert client.http_errors_file.read() == "inx15\nHTTP Error Fake HTTPError from Mock: Foo\n//\n"

    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_raise_urlerror):
        client.query_uniprot("inx15", [{"format": "list"}])
    assert "<urlopen error Fake URLError from Mock>" in client.http_errors_file.read()

    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_raise_keyboardinterrupt):
        client.query_uniprot("inx15", [{"format": "list"}])
    out, err = capsys.readouterr()
    assert "\n\tUniProt query interrupted by user\n" in err

    params = {"format": "tab", "columns": "id,entry name,length,organism-id,organism,protein names,comments"}
    client.query_uniprot("ABXEF9", params)


def test_uniprotrestclient_count_hits(capsys):
    dbbuddy = Db.DbBuddy("inx15,inx16")
    client = Db.UniProtRestClient(dbbuddy)
    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_uniprot_count_hits):
        assert client.count_hits() == 10

        for indx in range(10):
            client.dbbuddy.search_terms.append("a" * 110)
        assert client.count_hits() == 20

    with mock.patch('buddysuite.DatabaseBuddy.urlopen', mock_urlopen_raise_httperror):
        assert client.count_hits() == 0
        out, err = capsys.readouterr()
        assert "The following errors were encountered while querying UniProt with count_hits():" in err

    with pytest.raises(ValueError) as err:
        client.dbbuddy.search_terms[0] = "a" * 1001
        client.count_hits()
    assert "Search term exceeds size limit of 1000 characters." in str(err)


def test_uniprotrestclient_search_proteins(monkeypatch, capsys):
    def patch_query_uniprot_multi(*args, **kwargs):
        client1.results_file.write('''# Search: inx15
A8XEF9	A8XEF9_CAEBR	381	6238	Caenorhabditis briggsae	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
O61786	O61786_CAEEL	382	6239	Caenorhabditis elegans	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
A0A0H5SBJ0	A0A0H5SBJ0_BRUMA	129	6279	Brugia malayi (Filarial nematode worm)	Innexin	Function (1); Sequence similarities (1); Subcellular location (1)
E3MGD6	E3MGD6_CAERE	384	31234	Caenorhabditis remanei (Caenorhabditis vulgaris)	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
//
# Search: inx16
O61787	INX16_CAEEL	372	6239	Caenorhabditis elegans	Innexin-16 (Protein opu-16)	Function (1); Sequence similarities (1); Subcellular location (1)
A0A0V1AZ11	A0A0V1AZ11_TRISP	406	6334	Trichinella spiralis (Trichina worm)	Innexin	Caution (1); Function (1); Sequence similarities (1); Subcellular location (2)
A8XEF8	A8XEF8_CAEBR	374	6238	Caenorhabditis briggsae	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
A0A0B2VB60	A0A0B2VB60_TOXCA	366	6265	Toxocara canis (Canine roundworm)	Innexin	Caution (2); Function (1); Sequence similarities (1); Subcellular location (1)
A0A0V0W5E2	A0A0V0W5E2_9BILA	410	92179	Trichinella sp. T6	Innexin	Caution (2); Function (1); Sequence similarities (1); Subcellular location (1)
//''', "w")
        return

    def patch_query_uniprot_single(*args, **kwargs):
        client2.results_file.write('''# Search: inx15
A8XEF9	A8XEF9_CAEBR	381	6238	Caenorhabditis briggsae	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
O61786	O61786_CAEEL	382	6239	Caenorhabditis elegans	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
A0A0H5SBJ0	A0A0H5SBJ0_BRUMA	129	6279	Brugia malayi (Filarial nematode worm)	Innexin
E3MGD6	E3MGD6_CAERE	384	31234	Caenorhabditis remanei (Caenorhabditis vulgaris)	Innexin
//''', "w")
        return

    monkeypatch.setattr(Db.UniProtRestClient, "count_hits", lambda _: 0)
    dbbuddy = Db.DbBuddy("inx15,inx16")
    client1 = Db.UniProtRestClient(dbbuddy)
    client1.search_proteins()
    out, err = capsys.readouterr()
    assert "Uniprot returned no results\n\n" in err

    monkeypatch.setattr(Db.UniProtRestClient, "count_hits", lambda _: 9)
    monkeypatch.setattr(br, "run_multicore_function", patch_query_uniprot_multi)
    client1.search_proteins()
    out, err = capsys.readouterr()
    assert "Retrieving summary data for 9 records from UniProt\n" in err
    assert "Querying UniProt with 2 search terms (Ctrl+c to abort)\n" in err
    assert len(dbbuddy.records) == 9

    monkeypatch.setattr(Db.UniProtRestClient, "query_uniprot", patch_query_uniprot_single)
    dbbuddy = Db.DbBuddy("inx15")
    client2 = Db.UniProtRestClient(dbbuddy)
    client2.http_errors_file.write("inx15\n%s\n//\n" % URLError("Fake URLError from Mock"))
    client2.search_proteins()
    out, err = capsys.readouterr()
    assert "Querying UniProt with the search term 'inx15'...\n" in err
    assert "The following errors were encountered while querying UniProt with search_proteins():" in err
    assert len(dbbuddy.records) == 4


def test_uniprotrestclient_fetch_proteins(monkeypatch, capsys, sb_resources, sb_helpers):
    def patch_query_uniprot_search(*args, **kwargs):
        client.results_file.write('''# Search: inx15
A8XEF9	A8XEF9_CAEBR	381	6238	Caenorhabditis briggsae	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
O61786	O61786_CAEEL	382	6239	Caenorhabditis elegans	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
A0A0H5SBJ0	A0A0H5SBJ0_BRUMA	129	6279	Brugia malayi (Filarial nematode worm)	Innexin	Function (1); Sequence similarities (1); Subcellular location (1)
E3MGD6	E3MGD6_CAERE	384	31234	Caenorhabditis remanei (Caenorhabditis vulgaris)	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
//
# Search: inx16
O61787	INX16_CAEEL	372	6239	Caenorhabditis elegans	Innexin-16 (Protein opu-16)	Function (1); Sequence similarities (1); Subcellular location (1)
A0A0V1AZ11	A0A0V1AZ11_TRISP	406	6334	Trichinella spiralis (Trichina worm)	Innexin	Caution (1); Function (1); Sequence similarities (1); Subcellular location (2)
A8XEF8	A8XEF8_CAEBR	374	6238	Caenorhabditis briggsae	Innexin	Function (1); Sequence similarities (1); Subcellular location (2)
A0A0B2VB60	A0A0B2VB60_TOXCA	366	6265	Toxocara canis (Canine roundworm)	Innexin	Caution (2); Function (1); Sequence similarities (1); Subcellular location (1)
A0A0V0W5E2	A0A0V0W5E2_9BILA	410	92179	Trichinella sp. T6	Innexin	Caution (2); Function (1); Sequence similarities (1); Subcellular location (1)
//''', "w")
        return

    def patch_query_uniprot_fetch(*args, **kwargs):
        with open("%s/mock_resources/test_databasebuddy_clients/uniprot_fetch.txt" % sb_resources.res_path, "r") \
                as ifile:
            client.results_file.write(ifile.read(), "w")
        return

    def patch_query_uniprot_fetch_nothing(*args, **kwargs):
        client.results_file.write("# Search: A8XEF9,O61786,A0A0H5SBJ0,E3MGD6,O61787,A0A0V1AZ11,A8XEF8,A0A0B2VB60,"
                                  "A0A0V0W5E2\n//\n//", "w")
        return

    dbbuddy = Db.DbBuddy("inx15,inx16")
    client = Db.UniProtRestClient(dbbuddy)
    client.fetch_proteins()

    out, err = capsys.readouterr()
    assert err == out == ""

    # Test a single call to query_uniprot
    monkeypatch.setattr(Db.UniProtRestClient, "query_uniprot", patch_query_uniprot_search)
    client.search_proteins()
    monkeypatch.setattr(Db.UniProtRestClient, "query_uniprot", patch_query_uniprot_fetch)
    client.fetch_proteins()
    out, err = capsys.readouterr()
    assert "Requesting 9 full records from UniProt..." in err

    # Test multicore call to query_uniprot
    monkeypatch.setattr(br, "run_multicore_function", patch_query_uniprot_fetch)
    for accn, rec in client.dbbuddy.records.items():
        rec.record = None
    client.dbbuddy.records["a" * 999] = Db.Record("a" * 999, _database="uniprot")
    client.fetch_proteins()
    out, err = capsys.readouterr()
    assert "Requesting 10 full records from UniProt..." in err
    assert sb_helpers.string2hash(str(client.dbbuddy.records["A8XEF9"].record.seq)) == "04f13629336cf6cdd5859c8913b742a5"

    # Some edge cases
    monkeypatch.setattr(Db.UniProtRestClient, "query_uniprot", patch_query_uniprot_fetch_nothing)
    client.http_errors_file.write("inx15\n%s\n//\n" % URLError("Fake URLError from Mock"))

    client.dbbuddy.records = OrderedDict([("a" * 999, Db.Record("a" * 999, _database="uniprot"))])
    client.fetch_proteins()
    out, err = capsys.readouterr()
    assert "Requesting 1 full records from UniProt..." in err
    assert "No sequences returned\n\n" in err
    assert "The following errors were encountered while querying UniProt with fetch_proteins():" in err
    assert sb_helpers.string2hash(str(client.dbbuddy.records["a" * 999])) == "670bf9c6ae5832b42841798d882a7276"

    with pytest.raises(ValueError) as err:
        client.dbbuddy.records["a" * 1001] = Db.Record("a" * 1001, _database="uniprot")
        client.fetch_proteins()
    assert "The provided accession or search term is too long (>1000)." in str(err)


# NCBI
def test_ncbiclient_init():
    dbbuddy = Db.DbBuddy(", ".join(ACCNS[:3]))
    client = Db.NCBIClient(dbbuddy)
    assert client.Entrez.email == br.config_values()['email']
    assert client.Entrez.tool == "buddysuite"
    assert hash(dbbuddy) == hash(client.dbbuddy)
    assert type(client.http_errors_file) == br.TempFile
    assert type(client.results_file) == br.TempFile
    assert client.max_url == 1000
    assert client.max_attempts == 5


def test_ncbiclient_split_for_url():
    pass
