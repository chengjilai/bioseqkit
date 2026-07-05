"""Tests for reverse complement and six-frame translation."""

from bioseqkit.transform import (
    back_transcribe,
    reverse_complement,
    six_frame_translation,
    transcribe,
    translate,
)


def test_reverse_complement_known():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("ACGTN") == "NACGT"


def test_transcribe_and_back():
    assert transcribe("ATGC") == "AUGC"
    assert transcribe("atgc") == "augc"
    assert back_transcribe("AUGC") == "ATGC"
    # round trip
    assert back_transcribe(transcribe("ACGTACGT")) == "ACGTACGT"


def test_reverse_complement_lowercase():
    assert reverse_complement("atgc") == "gcat"


def test_translate_start_codon():
    # ATG GCC TAA -> M A *
    assert translate("ATGGCCTAA") == "MA*"


def test_translate_incomplete_codon_dropped():
    assert translate("ATGGC") == "M"


def test_six_frame_count_and_names():
    frames = six_frame_translation("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    assert len(frames) == 6
    assert [f.name for f in frames] == ["+1", "+2", "+3", "-1", "-2", "-3"]


def test_six_frame_reading_frame_coordinates():
    seq = "AAATGGCCTAA"  # +3 frame (offset 2) starts ATG GCC TAA
    frames = {f.name: f.protein for f in six_frame_translation(seq)}
    assert frames["+3"].startswith("MA*")
    # forward frame +1 uses offset 0
    assert frames["+1"] == translate(seq)
