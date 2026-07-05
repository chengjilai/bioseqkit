"""Tests for FAI-like indexing and random access."""

from pathlib import Path

import pytest

from bioseqkit.index import FaidxIndex, build_faidx, fetch, parse_region


@pytest.fixture
def indexed_fasta(tmp_path: Path) -> Path:
    p = tmp_path / "genome.fa"
    p.write_text(">chr1 test\n" + "ACGT" * 10 + "\n" + "ACGT" * 10 + "\n>chr2\nGGGGCCCCTTTTAAAA\n")
    return p


def test_build_and_fetch_whole(indexed_fasta):
    index = build_faidx(str(indexed_fasta))
    assert index.names() == ["chr1", "chr2"]
    assert index.records["chr1"].length == 80
    assert index.fetch("chr2") == "GGGGCCCCTTTTAAAA"


def test_fetch_region(indexed_fasta):
    index = build_faidx(str(indexed_fasta))
    # 1-based inclusive: bases 1-4 of chr1 == "ACGT"
    assert index.fetch("chr1", 1, 4) == "ACGT"
    # cross line boundary (line width 40); bases 39-42 span the two lines
    assert index.fetch("chr1", 39, 42) == "GTAC"


def test_index_write_load_roundtrip(indexed_fasta):
    index = build_faidx(str(indexed_fasta))
    path = index.write()
    assert Path(path).exists()
    loaded = FaidxIndex.load(str(indexed_fasta))
    assert loaded.fetch("chr1", 1, 4) == "ACGT"


def test_parse_region():
    assert parse_region("chr1:1000-2000") == ("chr1", 1000, 2000)
    assert parse_region("chr1") == ("chr1", None, None)
    assert parse_region("chr1:5") == ("chr1", 5, 5)


def test_fetch_helper(indexed_fasta):
    assert fetch(str(indexed_fasta), "chr2:1-4") == "GGGG"


def test_out_of_bounds(indexed_fasta):
    index = build_faidx(str(indexed_fasta))
    with pytest.raises(ValueError):
        index.fetch("chr1", 1, 999)
