"""Sequence statistics.

Length distribution, GC content, N-base ratio and base-composition matrix,
computed with the Python standard library only.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable

__all__ = [
    "SeqStats",
    "gc_content",
    "n_ratio",
    "base_composition",
    "sequence_stats",
]


def gc_content(sequence: str) -> float:
    """Fraction of G/C bases (case-insensitive). Returns 0.0 for empty input."""
    if not sequence:
        return 0.0
    gc = 0
    for base in sequence:
        if base in "GCgc":
            gc += 1
    return gc / len(sequence)


def n_ratio(sequence: str) -> float:
    """Fraction of ambiguous ``N`` bases (case-insensitive)."""
    if not sequence:
        return 0.0
    n = sequence.count("N") + sequence.count("n")
    return n / len(sequence)


def base_composition(sequence: str) -> dict[str, int]:
    """Return a per-base count dictionary (upper-cased keys)."""
    return dict(Counter(sequence.upper()))


@dataclass
class SeqStats:
    """Aggregate statistics over a collection of sequences."""

    n_seqs: int = 0
    total_length: int = 0
    min_length: int = 0
    max_length: int = 0
    lengths: list[int] = field(default_factory=list)
    gc_content: float = 0.0
    n_ratio: float = 0.0
    base_counts: dict[str, int] = field(default_factory=dict)

    @property
    def mean_length(self) -> float:
        return self.total_length / self.n_seqs if self.n_seqs else 0.0

    def n50(self) -> int:
        """Return the N50 of the length distribution."""
        if not self.lengths:
            return 0
        half = self.total_length / 2
        acc = 0
        for length in sorted(self.lengths, reverse=True):
            acc += length
            if acc >= half:
                return length
        return 0

    def as_dict(self) -> dict[str, object]:
        return {
            "n_seqs": self.n_seqs,
            "total_length": self.total_length,
            "min_length": self.min_length,
            "max_length": self.max_length,
            "mean_length": round(self.mean_length, 3),
            "n50": self.n50(),
            "gc_content": round(self.gc_content, 6),
            "n_ratio": round(self.n_ratio, 6),
            "base_counts": self.base_counts,
        }


def sequence_stats(sequences: Iterable[str]) -> SeqStats:
    """Compute aggregate :class:`SeqStats` over an iterable of sequences."""
    stats = SeqStats()
    gc = 0
    n_bases = 0
    counts: Counter[str] = Counter()
    for seq in sequences:
        length = len(seq)
        stats.n_seqs += 1
        stats.total_length += length
        stats.lengths.append(length)
        upper = seq.upper()
        counts.update(upper)
        gc += upper.count("G") + upper.count("C")
        n_bases += upper.count("N")
    if stats.lengths:
        stats.min_length = min(stats.lengths)
        stats.max_length = max(stats.lengths)
    if stats.total_length:
        stats.gc_content = gc / stats.total_length
        stats.n_ratio = n_bases / stats.total_length
    stats.base_counts = dict(counts)
    return stats
