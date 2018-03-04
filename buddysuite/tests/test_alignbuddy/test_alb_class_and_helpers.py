#!/usr/bin/env python3
# coding=utf-8
""" tests basic functionality of AlignBuddy class """
import pytest
import io
import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import AlignIO

import buddy_resources as br
from AlignBuddy import AlignBuddy, guess_alphabet, make_copy
from buddy_resources import GuessError, parse_format


def mock_valueerror(*args, **kwargs):
    raise ValueError(args, kwargs)


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
        with open(_path, 'r', encoding="utf-8") as ifile:
            assert type(AlignBuddy(ifile)) == AlignBuddy


def test_instantiate_alignbuddy_from_raw(alb_resources):
    for _path in alb_resources.get_list(mode="paths"):
        with open(_path, 'r', encoding="utf-8") as ifile:
            output = ifile.read()
            output = br.utf_encode(output)
            print(output)
            assert type(AlignBuddy(output)) == AlignBuddy


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


def test_instantiation_alignbuddy_errors(alb_odd_resources):
    with pytest.raises(GuessError) as e:
        AlignBuddy(alb_odd_resources["dna"]["single"]["fasta"])
    assert "Could not determine format from _input file" in str(e)

    tester = open(alb_odd_resources["dna"]["single"]["fasta"], "r", encoding="utf-8")
    with pytest.raises(GuessError) as e:
        AlignBuddy(tester.read())
    assert "Could not determine format from raw" in str(e)

    tester.seek(0)
    with pytest.raises(GuessError) as e:
        AlignBuddy(tester)
    assert "Could not determine format from input file-like object" in str(e)


def test_empty_file(alb_odd_resources):
    with open(alb_odd_resources["blank"], "r", encoding="utf-8") as ifile:
        with pytest.raises(GuessError) as e:
            AlignBuddy(ifile)
        assert "Empty file" in str(e)


def test_throws_errors_on_invalid_files(alb_odd_resources):
    """ expect AlignBuddy to raise errors on invalid filesr """
    with pytest.raises(GuessError):
        AlignBuddy(alb_odd_resources['dna']['single']['fasta'])


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


def test_records_dict(alb_resources, hf):
    alignbuddy = alb_resources.get_one("o p g")
    alb_dict = alignbuddy.records_dict()
    assert hf.string2hash(str(alb_dict)) == "a11e822d85aa7dc43afad3eda4f1708d"


def test_lengths_single(alb_resources):
    for alignbuddy in alb_resources.get_list("o p g py pr pss psr"):
        assert alignbuddy.lengths()[0] == 681


def test_lengths_multi(alb_resources):
    for alignbuddy in alb_resources.get_list("m p py pr pss psr"):
        assert alignbuddy.lengths()[1] == 480

hashes = [('o p g', '46388b175b31b81f47199ae6327768af'), ('o p n', '17ff1b919cac899c5f918ce8d71904f6'),
          ('o p py', '968ed9fa772e65750f201000d7da670f'), ('o p pr', 'ce423d5b99d5917fbef6f3b47df40513'),
          ('o p pss', '4bd927145de635c429b2917e0a1db176'), ('o p psr', '8ff80c7f0b8fc7f237060f94603c17be'),
          ('o p s', 'c0dce60745515b31a27de1f919083fe9'),

          ('o d c', '3c937c9fec251a42f0994caabb64420c'), ('o d f', '98a3a08389284461ea9379c217e99770'),
          ('o d g', '842d9c911a33c0fd0484383eabefb0fe'), ('o d n', 'cb1169c2dd357771a97a02ae2160935d'),
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
def test_str(alb_resources, hf, key, next_hash):
    tester = str(alb_resources.get_one(key))
    assert hf.string2hash(tester) == next_hash, open("error_files/%s" % next_hash, "w", encoding="utf-8").write(tester)


def test_str2(alb_resources, hf, capsys, monkeypatch):
    alignbuddy = alb_resources.get_one("m p c")
    alignbuddy.alignments[0] = []
    assert hf.string2hash(str(alignbuddy)) == "cd017eee573d2cf27eb7c161babe0cad"

    alignbuddy = alb_resources.get_one("o d f")
    alignbuddy.alignments[0] = []
    assert str(alignbuddy) == "AlignBuddy object contains no alignments.\n"

    alignbuddy = alb_resources.get_one("o p g")
    for rec in alignbuddy.records():
        rec.annotations['organism'] = ". . "
    assert hf.string2hash(str(alignbuddy)) == "46388b175b31b81f47199ae6327768af"

    alignbuddy = alb_resources.get_one("o p g")
    for rec in alignbuddy.records():
        del rec.annotations['organism']
    assert hf.string2hash(str(alignbuddy)) == "46388b175b31b81f47199ae6327768af"

    alignbuddy = alb_resources.get_one("m p c")
    alignbuddy.set_format("genbank")
    with pytest.raises(ValueError) as err:
        str(alignbuddy)
    assert "genbank format does not support multiple alignments in one file." in str(err)

    alignbuddy = alb_resources.get_one("o p g")
    alignbuddy.set_format("phylipsr")
    assert hf.string2hash(str(alignbuddy)) == "8ff80c7f0b8fc7f237060f94603c17be"

    alignbuddy = alb_resources.get_one("o p py")
    alignbuddy.set_format("phylipss")
    assert hf.string2hash(str(alignbuddy)) == "4bd927145de635c429b2917e0a1db176"

    alignbuddy = alb_resources.get_one("o p g")
    alignbuddy.set_format("phylip")
    assert hf.string2hash(str(alignbuddy)) == "ce423d5b99d5917fbef6f3b47df40513"
    out, err = capsys.readouterr()
    assert "Warning: Phylip format returned a 'repeat name' error, probably due to truncation." in err

    alignbuddy = alb_resources.get_one("o p py")
    rec = alignbuddy.records()[0]
    rec.seq = Seq(str(rec.seq)[:-2], alphabet=rec.seq.alphabet)
    assert hf.string2hash(str(alignbuddy)) == "9337bba9fb455f1e6257cc236a663001"
    out, err = capsys.readouterr()
    assert "Warning: Alignment format detected but sequences are different lengths." in err

    monkeypatch.setattr(AlignIO, "write", mock_valueerror)
    alignbuddy = alb_resources.get_one("o p py")
    with pytest.raises(ValueError):
        str(alignbuddy)


@pytest.mark.parametrize('key,next_hash', hashes)
def test_write1(alb_resources, hf, key, next_hash):
    temp_file = br.TempFile()
    alignbuddy = alb_resources.get_one(key)
    alignbuddy.write(temp_file.path)
    with open(temp_file.path, "r", encoding="utf-8") as ifile:
        output = ifile.read()
        tester_hash = hf.string2hash(output)
    assert tester_hash == next_hash

hashes = [('m p c', '9c6773e7d24000f8b72dd9d25620cff1'), ('m p s', '9c6773e7d24000f8b72dd9d25620cff1'),
          ('m p py', '1f172a3beef76e8e3d42698bb2c3c87d'), ('m p pr', '3fef9a05058a5259ebd517d1500388d4'),
          ('m p pss', '1f172a3beef76e8e3d42698bb2c3c87d'), ('m p psr', '3fef9a05058a5259ebd517d1500388d4')]


@pytest.mark.parametrize("key,next_hash", hashes)
def test_write2(alb_resources, hf, key, next_hash):
    temp_file = br.TempFile()
    alignbuddy = alb_resources.get_one(key)
    alignbuddy.write(temp_file.path, out_format="phylipr")
    with open(temp_file.path, "r", encoding="utf-8") as ifile:
        tester_hash = hf.string2hash(ifile.read())
    assert tester_hash == next_hash


def test_write3(alb_resources, hf):  # Unloopable components
    tester = alb_resources.get_one("m p py")
    tester.set_format("fasta")
    with pytest.raises(ValueError):
        str(tester)

    tester.alignments = []
    assert str(tester) == "AlignBuddy object contains no alignments.\n"

    tester = alb_resources.get_one("o d pr")
    tester.set_format("phylipi")
    assert hf.buddy2hash(tester) == "52c23bd793c9761b7c0f897d3d757c12"

    tester = AlignBuddy("%s/Mnemiopsis_cds_hashed_ids.nex" % hf.resource_path)
    tester.set_format("phylip-strict")
    assert hf.buddy2hash(tester) == "16b3397d6315786e8ad8b66e0d9c798f"


# ################################################# HELPER FUNCTIONS ################################################# #
def test_guess_error(alb_odd_resources):
    # File path
    with pytest.raises(GuessError):
        unrecognizable = alb_odd_resources['protein']['single']['phylip']
        AlignBuddy(unrecognizable)

    with open(unrecognizable, 'r', encoding="utf-8") as ifile:
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


def test_make_copy(alb_resources, hf):
    for alb in alb_resources.get_list():
        tester = make_copy(alb)
        assert hf.buddy2hash(tester) == hf.buddy2hash(alb)

# ToDo: def test_feature_remapper()
