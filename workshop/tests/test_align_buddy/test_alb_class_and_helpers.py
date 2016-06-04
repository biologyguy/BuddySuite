#!/usr/bin/env python3
# coding=utf-8
""" tests basic functionality of AlignBuddy class """
import pytest
import MyFuncs
import io
from AlignBuddy import AlignBuddy, guess_alphabet, guess_format
from buddy_resources import GuessError, parse_format, PhylipError
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def test_instantiate_alignbuddy_from_file(alb_resources):
    for key, _path in alb_resources.get(mode="paths").items():
        in_format = alb_resources.parse_code(key, strict=True)
        in_format = alb_resources.single_letter_codes[in_format["format"][0]]
        assert type(AlignBuddy(_path, in_format=in_format)) == AlignBuddy


def test_instantiate_alignbuddy_from_file_guess(alb_resources):
    for _path in alb_resources.get_list(mode="paths"):
        assert type(AlignBuddy(_path)) == AlignBuddy


def test_instantiate_alignbuddy_from_handle(alb_resources):
    for _path in alb_resources.get_list(mode="paths"):
        with open(_path, 'r') as ifile:
            assert type(AlignBuddy(ifile)) == AlignBuddy


def test_instantiate_alignbuddy_from_raw(alb_resources):
    for _path in alb_resources.get_list(mode="paths"):
        with open(_path, 'r') as ifile:
            assert type(AlignBuddy(ifile.read())) == AlignBuddy


def test_instantiate_alignbuddy_from_alignbuddy(alb_resources):
    for alignbuddy in alb_resources.get_list():
        assert type(AlignBuddy(alignbuddy)) == AlignBuddy


def test_instantiate_alignbuddy_from_list(alb_resources):
    for alignbuddy in alb_resources.get_list():
        assert type(AlignBuddy(alignbuddy.alignments)) == AlignBuddy

    with pytest.raises(TypeError):  # When non-MultipleSeqAlignment objects are in the .alignments list
        alignbuddy = alb_resources.get_one("p s m")
        alignbuddy.alignments.append("Dummy string object")
        AlignBuddy(alignbuddy.alignments)


def test_instantiation_alignbuddy_errors(alignment_bad_resources):
    with pytest.raises(GuessError) as e:
        AlignBuddy(alignment_bad_resources["dna"]["single"]["fasta"])
    assert "Could not determine format from _input file" in str(e)

    tester = open(alignment_bad_resources["dna"]["single"]["fasta"], "r")
    with pytest.raises(GuessError) as e:
        AlignBuddy(tester.read())
    assert "Could not determine format from raw" in str(e)

    tester.seek(0)
    with pytest.raises(GuessError) as e:
        AlignBuddy(tester)
    assert "Could not determine format from input file-like object" in str(e)


def test_empty_file(alignment_bad_resources):
    with open(alignment_bad_resources["blank"], "r") as ifile:
        with pytest.raises(GuessError) as e:
            AlignBuddy(ifile)
        assert "Empty file" in str(e)


def test_throws_errors_on_invalid_files(alignment_bad_resources):
    """ expect AlignBuddy to raise errors on invalid filesr """
    with pytest.raises(GuessError):
        AlignBuddy(alignment_bad_resources['dna']['single']['fasta'])


# ##################### AlignBuddy methods ###################### ##
def test_set_format(alb_resources):
    tester = alb_resources.get_list("o d g")[0]
    tester.set_format("fasta")
    assert tester.out_format == "fasta"


def test_records(alb_resources):
    tester = alb_resources.get_list("m p py")[0]
    assert len(tester.records()) == 29


def test_records_iter(alb_resources):
    tester = alb_resources.get_list("m p py")[0]
    counter = 0
    for rec in tester.records_iter():
        assert type(rec) == SeqRecord
        counter += 1
    assert counter == 29


def test_lengths_single(alb_resources):
    for alignbuddy in alb_resources.get_list("o p g py pr pss psr"):
        assert alignbuddy.lengths()[0] == 681


def test_lengths_multi(alb_resources):
    for alignbuddy in alb_resources.get_list("m p py pr pss psr"):
        assert alignbuddy.lengths()[1] == 480

hashes = [('o p g', 'bf8485cbd30ff8986c2f50b677da4332'), ('o p n', '17ff1b919cac899c5f918ce8d71904f6'),
          ('o p py', '968ed9fa772e65750f201000d7da670f'), ('o p pr', 'ce423d5b99d5917fbef6f3b47df40513'),
          ('o p pss', '4bd927145de635c429b2917e0a1db176'), ('o p psr', '8ff80c7f0b8fc7f237060f94603c17be'),
          ('o p s', 'c0dce60745515b31a27de1f919083fe9'),

          ('o d c', '3c937c9fec251a42f0994caabb64420c'), ('o d f', '98a3a08389284461ea9379c217e99770'),
          ('o d g', '2a42c56df314609d042bdbfa742871a3'), ('o d n', 'cb1169c2dd357771a97a02ae2160935d'),
          ('o d py', '503e23720beea201f8fadf5dabda75e4'), ('o d pr', '52c23bd793c9761b7c0f897d3d757c12'),
          ('o d pss', '4c0c1c0c63298786e6fb3db1385af4d5'), ('o d psr', 'c5fb6a5ce437afa1a4004e4f8780ad68'),
          ('o d s', '228e36a30e8433e4ee2cd78c3290fa6b'),

          ('o r n', 'f3bd73151645359af5db50d2bdb6a33d'),

          ('m p c', 'f0e20a55f679ee492bb0b3be444b46f9'), ('m p s', '3fd5805f61777f7f329767c5f0fb7467'),
          ('m p py', '2a77f5761d4f51b88cb86b079e564e3b'), ('m p pr', '3fef9a05058a5259ebd517d1500388d4'),
          ('m p pss', 'eb82cda31fcb2cf00e11d7e910fde695'), ('m p psr', 'a16c6e617e5a88fef080eea54e54e8a8'),

          ('m d c', '058ef1525cfc1364f26dd5a5bd6b97fb'), ('m d s', 'ae352b908be94738d6d9cd54770e5b5d'),
          ('m d py', '42679a32ebd93b628303865f68b0293d'), ('m d pr', '22c0f0c8f014a34be8edd394bf477a2d'),
          ('m d pss', 'c789860da8f0b59e0adc7bde6342b4b0'), ('m d psr', '28b2861275e0a488042cff35393ac36d')]


@pytest.mark.parametrize('key,next_hash', hashes)
def test_str(alb_resources, helpers, key, next_hash):
    tester = str(alb_resources.get_one(key))
    assert helpers.string2hash(tester) == next_hash, open("error_files/%s" % next_hash, "w").write(tester)


@pytest.mark.parametrize('key,next_hash', hashes)
def test_write1(alb_resources, helpers, key, next_hash):
    temp_file = MyFuncs.TempFile()
    alignbuddy = alb_resources.get_one(key)
    alignbuddy.write(temp_file.path)
    tester_hash = helpers.string2hash(temp_file.read())
    assert tester_hash == next_hash, alignbuddy.write("error_files/%s" % next_hash)

hashes = [('m p c', '9c6773e7d24000f8b72dd9d25620cff1'), ('m p s', '9c6773e7d24000f8b72dd9d25620cff1'),
          ('m p py', '1f172a3beef76e8e3d42698bb2c3c87d'), ('m p pr', '3fef9a05058a5259ebd517d1500388d4'),
          ('m p pss', '1f172a3beef76e8e3d42698bb2c3c87d'), ('m p psr', '3fef9a05058a5259ebd517d1500388d4')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_write2(alb_resources, helpers, key, next_hash):
    alignbuddy = alb_resources.get_one(key)
    temp_file = MyFuncs.TempFile()
    alignbuddy.write(temp_file.path, out_format="phylipr")
    out = temp_file.read()
    assert helpers.string2hash(out) == next_hash, alignbuddy.write("error_files/%s" % next_hash)


def test_write3(alb_resources, helpers):  # Unloopable components
    tester = alb_resources.get_one("m p py")
    tester.set_format("fasta")
    with pytest.raises(ValueError):
        str(tester)

    tester.alignments = []
    assert str(tester) == "AlignBuddy object contains no alignments.\n"

    tester = alb_resources.get_one("o d pr")
    tester.set_format("phylipi")
    assert helpers.align_to_hash(tester) == "52c23bd793c9761b7c0f897d3d757c12"

    tester = AlignBuddy("%s/Mnemiopsis_cds_hashed_ids.nex" % helpers.resource_path)
    tester.set_format("phylip-strict")
    assert helpers.align_to_hash(tester) == "16b3397d6315786e8ad8b66e0d9c798f"


# ################################################# HELPER FUNCTIONS ################################################# #
def test_guess_error(alignment_bad_resources):
    # File path
    with pytest.raises(GuessError):
        unrecognizable = alignment_bad_resources['protein']['single']['phylip']
        AlignBuddy(unrecognizable)

    with open(unrecognizable, 'r') as ifile:
        # Raw
        with pytest.raises(GuessError) as e:
            AlignBuddy(ifile.read())
        assert "Could not determine format from raw input" in str(e)

        # Handle
        with pytest.raises(GuessError) as e:
            ifile.seek(0)
            AlignBuddy(ifile)
        assert "Could not determine format from input file-like object" in str(e)

    # GuessError output
    try:
        AlignBuddy(unrecognizable)
    except GuessError as e:
        assert "Could not determine format from _input file" in str(e) and \
               "\nTry explicitly setting with -f flag." in str(e)


def test_guess_alphabet(alb_resources):
    for alb in alb_resources.get_list("d"):
        assert guess_alphabet(alb) == IUPAC.ambiguous_dna
    for alb in alb_resources.get_list("p"):
        assert guess_alphabet(alb) == IUPAC.protein
    for alb in alb_resources.get_list("r"):
        assert guess_alphabet(alb) == IUPAC.ambiguous_rna

    assert not guess_alphabet(AlignBuddy("", in_format="fasta"))


def test_guess_format(alb_resources, alignment_bad_resources):
    assert guess_format(["dummy", "list"]) == "stockholm"

    for key, obj in alb_resources.get().items():
        assert guess_format(obj) == parse_format(alb_resources.get_key(key)["format"])

    for key, path in alb_resources.get(mode="paths").items():
        assert guess_format(path) == parse_format(alb_resources.get_key(key)["format"])
        with open(path, "r") as ifile:
            assert guess_format(ifile) == parse_format(alb_resources.get_key(key)["format"])
            ifile.seek(0)
            string_io = io.StringIO(ifile.read())
        assert guess_format(string_io) == parse_format(alb_resources.get_key(key)["format"])

    guess_format(alignment_bad_resources['blank']) == "empty file"
    assert not guess_format(alignment_bad_resources['dna']['single']['phylipss_recs'])
    assert not guess_format(alignment_bad_resources['dna']['single']['phylipss_cols'])

    with pytest.raises(GuessError) as e:
        guess_format({"Dummy dict": "Type not recognized by guess_format()"})
    assert "Unsupported _input argument in guess_format()" in str(e)

"""
def test_make_copy(alb_resources):
    for alb in alb_resources.get_list():
        tester = make_copy(alb)
        align_to_hash(tester) == align_to_hash(alb)


def test_stderr(capsys):
    _stderr("Hello std_err", quiet=False)
    out, err = capsys.readouterr()
    assert err == "Hello std_err"

    _stderr("Hello std_err", quiet=True)
    out, err = capsys.readouterr()
    assert err == ""


def test_stdout(capsys):
    _stdout("Hello std_out", quiet=False)
    out, err = capsys.readouterr()
    assert out == "Hello std_out"

    _stdout("Hello std_out", quiet=True)
    out, err = capsys.readouterr()
    assert out == ""


# ToDo: def test_feature_remapper()
"""
