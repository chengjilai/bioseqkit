"""k-mer analysis: counting, top-k, canonical k-mers, parallel counting and minimizers."""

from __future__ import annotations

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from typing import Iterable, Iterator

from bioseqkit.transform import reverse_complement

__all__ = [
    "iter_kmers",
    "count_kmers",
    "count_kmers_parallel",
    "top_kmers",
    "canonical_kmer",
    "minimizers",
]


def iter_kmers(sequence: str, k: int) -> Iterator[str]:
    """Yield successive k-mers of length ``k`` from ``sequence``."""
    if k <= 0:
        raise ValueError("k must be a positive integer")
    for i in range(len(sequence) - k + 1):
        yield sequence[i : i + k]


def canonical_kmer(kmer: str) -> str:
    """Return the lexicographically smaller of a k-mer and its reverse complement."""
    rc = reverse_complement(kmer)
    return kmer if kmer <= rc else rc


def count_kmers(
    sequence: str,
    k: int,
    canonical: bool = False,
    skip_ambiguous: bool = True,
) -> Counter[str]:
    """Count k-mers in a single sequence.

    If ``canonical`` is True, a k-mer and its reverse complement are merged.
    If ``skip_ambiguous`` is True, k-mers containing bases outside ``ACGT`` are
    skipped.
    """
    counts: Counter[str] = Counter()
    seq = sequence.upper()
    valid = set("ACGT")
    for kmer in iter_kmers(seq, k):
        if skip_ambiguous and not set(kmer) <= valid:
            continue
        counts[canonical_kmer(kmer) if canonical else kmer] += 1
    return counts


def _count_chunk(args: tuple[str, int, bool, bool]) -> Counter[str]:
    sequence, k, canonical, skip_ambiguous = args
    return count_kmers(sequence, k, canonical, skip_ambiguous)


def _chunk_sequence(sequence: str, n_chunks: int, k: int) -> list[str]:
    """Split a sequence into ``n_chunks`` overlapping chunks (overlap = k - 1).

    The overlap guarantees k-mers spanning chunk boundaries are still counted.
    """
    n = len(sequence)
    if n_chunks <= 1 or n < 2 * k:
        return [sequence]
    size = max(k, n // n_chunks)
    chunks: list[str] = []
    start = 0
    while start < n:
        end = min(n, start + size)
        chunks.append(sequence[start : min(n, end + k - 1)])
        start = end
    return chunks


def count_kmers_parallel(
    sequences: Iterable[str],
    k: int,
    canonical: bool = False,
    skip_ambiguous: bool = True,
    workers: int = 4,
) -> Counter[str]:
    """Count k-mers across sequences using a process pool, then merge counts.

    Each input sequence is split into overlapping chunks that are distributed
    to worker processes; the per-chunk counters are summed into a single result.
    """
    chunks: list[str] = []
    for seq in sequences:
        chunks.extend(_chunk_sequence(seq.upper(), workers, k))
    if not chunks:
        return Counter()
    if workers <= 1 or len(chunks) == 1:
        total: Counter[str] = Counter()
        for chunk in chunks:
            total += count_kmers(chunk, k, canonical, skip_ambiguous)
        return total

    tasks = [(chunk, k, canonical, skip_ambiguous) for chunk in chunks]
    total = Counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for partial in pool.map(_count_chunk, tasks):
            total += partial
    return total


def top_kmers(counts: Counter[str], n: int = 10) -> list[tuple[str, int]]:
    """Return the ``n`` most common (k-mer, count) pairs."""
    return counts.most_common(n)


def minimizers(sequence: str, k: int, w: int, canonical: bool = True) -> list[tuple[int, str]]:
    """Compute (position, minimizer) pairs over a sliding window.

    For each window of ``w`` consecutive k-mers, the lexicographically smallest
    (canonical, by default) k-mer is selected. Consecutive duplicate minimizers
    are collapsed, mirroring minimap2/Mash behaviour.
    """
    if k <= 0 or w <= 0:
        raise ValueError("k and w must be positive integers")
    seq = sequence.upper()
    kmers = list(iter_kmers(seq, k))
    if not kmers:
        return []
    keyed = [canonical_kmer(km) if canonical else km for km in kmers]
    result: list[tuple[int, str]] = []
    last: tuple[int, str] | None = None
    n_windows = max(1, len(kmers) - w + 1)
    for start in range(n_windows):
        window = keyed[start : start + w]
        best_val = min(window)
        best_pos = start + window.index(best_val)
        picked = (best_pos, best_val)
        if picked != last:
            result.append(picked)
            last = picked
    return result
