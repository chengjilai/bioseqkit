"""Micro-benchmarks executed via pytest-benchmark.

These are skipped automatically when the ``pytest-benchmark`` plugin is not
installed, so the ordinary test suite remains dependency-free. Run them with::

    pytest tests/test_benchmark.py --benchmark-only
"""

import random

import pytest

import bioseqkit as bsk

pytest.importorskip("pytest_benchmark")


@pytest.fixture(scope="module")
def sequence() -> str:
    rng = random.Random(0)
    return "".join(rng.choice("ACGT") for _ in range(200_000))


def test_count_kmers_benchmark(benchmark, sequence):
    result = benchmark(bsk.count_kmers, sequence, 8, canonical=True)
    assert sum(result.values()) > 0


def test_faidx_build_benchmark(benchmark, tmp_path, sequence):
    fa = tmp_path / "g.fa"
    with open(fa, "w") as fh:
        fh.write(">chr1\n")
        for i in range(0, len(sequence), 70):
            fh.write(sequence[i : i + 70] + "\n")
    index = benchmark(bsk.build_faidx, str(fa))
    assert index.records["chr1"].length == len(sequence)
