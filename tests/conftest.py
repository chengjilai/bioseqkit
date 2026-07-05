"""Shared pytest fixtures and sample data."""

import gzip
from pathlib import Path

import pytest

SINGLE_FASTA = ">seq1 example single\nACGTACGTAC\n"

MULTILINE_FASTA = ">seq2 multi line\nACGT\nACGT\nAC\n>seq3\nGGGGCCCC\n"

FASTQ = (
    "@read1 first\n"
    "ACGTACGT\n"
    "+\n"
    "IIIIIIII\n"
    "@read2\n"
    "GGGGCCCC\n"
    "+\n"
    "!!!!####\n"
)


@pytest.fixture
def single_fasta(tmp_path: Path) -> Path:
    p = tmp_path / "single.fa"
    p.write_text(SINGLE_FASTA)
    return p


@pytest.fixture
def multiline_fasta(tmp_path: Path) -> Path:
    p = tmp_path / "multi.fa"
    p.write_text(MULTILINE_FASTA)
    return p


@pytest.fixture
def gz_fasta(tmp_path: Path) -> Path:
    p = tmp_path / "single.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(SINGLE_FASTA)
    return p


@pytest.fixture
def fastq_file(tmp_path: Path) -> Path:
    p = tmp_path / "reads.fq"
    p.write_text(FASTQ)
    return p
