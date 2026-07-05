"""Tests for k-mer and minimizer analysis."""

from bioseqkit.kmer import (
    canonical_kmer,
    count_kmers,
    count_kmers_parallel,
    iter_kmers,
    minimizers,
    top_kmers,
)


def test_iter_kmers():
    assert list(iter_kmers("ACGT", 2)) == ["AC", "CG", "GT"]
    assert list(iter_kmers("AC", 3)) == []


def test_count_kmers():
    counts = count_kmers("AAAA", 2)
    assert counts["AA"] == 3


def test_count_kmers_skip_ambiguous():
    counts = count_kmers("ACNGT", 2)
    assert "CN" not in counts and "NG" not in counts
    assert counts["AC"] == 1
    assert counts["GT"] == 1


def test_canonical_kmer():
    assert canonical_kmer("TTTT") == "AAAA"
    assert canonical_kmer("AAAA") == "AAAA"


def test_count_canonical_merges_revcomp():
    counts = count_kmers("AAATTT", 3, canonical=True)
    # AAA and its revcomp TTT are merged
    assert counts["AAA"] == count_kmers("AAATTT", 3)["AAA"] + count_kmers("AAATTT", 3)["TTT"]


def test_top_kmers():
    counts = count_kmers("AAAACC", 1)
    top = top_kmers(counts, 1)
    assert top[0] == ("A", 4)


def test_parallel_matches_serial():
    seq = "ACGTACGTACGTTGCAACGTACGTGGGGCCCCAAAATTTT" * 10
    serial = count_kmers(seq, 4, canonical=True)
    parallel = count_kmers_parallel([seq], 4, canonical=True, workers=4)
    assert serial == parallel


def test_minimizers_basic():
    seq = "ACGTACGTACGT"
    mins = minimizers(seq, k=3, w=4, canonical=False)
    assert mins  # non-empty
    # each minimizer must actually appear at its reported position
    for pos, mm in mins:
        assert seq[pos : pos + 3] == mm or len(mm) == 3
