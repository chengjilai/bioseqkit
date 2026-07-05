"""Tests for FASTA/FASTQ parsing (io module)."""

import io as _io

import pytest

from bioseqkit.io import parse_fasta, parse_fastq, write_fasta


def test_single_fasta(single_fasta):
    records = list(parse_fasta(str(single_fasta)))
    assert len(records) == 1
    rec = records[0]
    assert rec.id == "seq1"
    assert rec.description == "example single"
    assert rec.sequence == "ACGTACGTAC"
    assert len(rec) == 10


def test_multiline_fasta(multiline_fasta):
    records = list(parse_fasta(str(multiline_fasta)))
    assert [r.id for r in records] == ["seq2", "seq3"]
    assert records[0].sequence == "ACGTACGTAC"  # multi-line joined
    assert records[1].sequence == "GGGGCCCC"


def test_gzip_fasta(gz_fasta):
    records = list(parse_fasta(str(gz_fasta)))
    assert len(records) == 1
    assert records[0].sequence == "ACGTACGTAC"


def test_blank_lines_and_stream():
    text = ">a\n\nACG\n\nTAC\n\n>b\nGGG\n"
    records = list(parse_fasta(_io.StringIO(text)))
    assert records[0].sequence == "ACGTAC"
    assert records[1].sequence == "GGG"


def test_sequence_before_header_raises():
    with pytest.raises(ValueError):
        list(parse_fasta(_io.StringIO("ACGT\n>a\nACGT\n")))


def test_illegal_characters_are_preserved_not_crashing():
    # Records with unexpected/illegal characters (digits, gaps, ambiguity codes,
    # lowercase) must be parsed without error; validation is left to the caller.
    text = ">weird desc\nACGT-N*12\nacgtRYKM\n"
    records = list(parse_fasta(_io.StringIO(text)))
    assert len(records) == 1
    assert records[0].sequence == "ACGT-N*12acgtRYKM"
    # downstream stats should still work on such input
    from bioseqkit.stats import gc_content

    assert 0.0 <= gc_content(records[0].sequence) <= 1.0


def test_fastq_parsing(fastq_file):
    records = list(parse_fastq(str(fastq_file)))
    assert len(records) == 2
    assert records[0].id == "read1"
    assert records[0].sequence == "ACGTACGT"
    assert records[0].quality == "IIIIIIII"
    assert records[0].phred_scores()[0] == 40
    assert records[1].phred_scores() == [0, 0, 0, 0, 2, 2, 2, 2]


def test_fastq_length_mismatch_raises():
    bad = "@r\nACGT\n+\nII\n"
    with pytest.raises(ValueError):
        list(parse_fastq(_io.StringIO(bad)))


def test_write_fasta_roundtrip(multiline_fasta):
    records = list(parse_fasta(str(multiline_fasta)))
    buf = _io.StringIO()
    n = write_fasta(records, buf, width=4)
    assert n == 2
    reparsed = list(parse_fasta(_io.StringIO(buf.getvalue())))
    assert [r.sequence for r in reparsed] == [r.sequence for r in records]
