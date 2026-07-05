"""Tests for sequence statistics."""

import math

from bioseqkit.stats import base_composition, gc_content, n_ratio, sequence_stats


def test_gc_content():
    assert gc_content("GGCC") == 1.0
    assert gc_content("ATAT") == 0.0
    assert math.isclose(gc_content("ACGT"), 0.5)
    assert gc_content("") == 0.0


def test_n_ratio():
    assert math.isclose(n_ratio("ACGTNNNN"), 0.5)
    assert n_ratio("ACGT") == 0.0


def test_base_composition():
    comp = base_composition("AACGTn")
    assert comp["A"] == 2
    assert comp["N"] == 1


def test_sequence_stats_aggregate():
    stats = sequence_stats(["ACGT", "GGGGCCCC", "AN"])
    assert stats.n_seqs == 3
    assert stats.total_length == 14
    assert stats.min_length == 2
    assert stats.max_length == 8
    assert math.isclose(stats.mean_length, 14 / 3)
    assert math.isclose(stats.n_ratio, 1 / 14)
    assert stats.base_counts["G"] == 5  # 4 in GGGGCCCC + 1 in ACGT
    d = stats.as_dict()
    assert d["n_seqs"] == 3


def test_n50():
    stats = sequence_stats(["A" * 100, "A" * 50, "A" * 50])
    assert stats.n50() == 100
