"""Benchmark and complexity-analysis driver for bioseqkit.

Measures wall-clock time and empirical scaling of the performance-critical
operations across a range of input sizes:

* k-mer counting        -- expected O(n k)
* parallel k-mer counting (2 / 4 workers)
* FAI-like index build  -- expected O(n)
* random sub-sequence fetch -- expected O(1) per query

Results are printed as a table, written to ``benchmarks/results.csv`` and, if
matplotlib is available, plotted to ``benchmarks/scalability.png``.

Run with:  ``python benchmarks/bench.py``  or  ``python benchmarks/bench.py --quick``
"""

from __future__ import annotations

import argparse
import csv
import random
import tempfile
import time
from pathlib import Path

import bioseqkit as bsk

HERE = Path(__file__).resolve().parent


def random_sequence(n: int, seed: int = 42) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_fasta(seq: str, path: Path, width: int = 70) -> None:
    with open(path, "w") as fh:
        fh.write(">chr1\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i : i + width] + "\n")


def time_call(fn, repeats: int = 1) -> float:
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def run(sizes: list[int], k: int = 8) -> list[dict]:
    rows: list[dict] = []
    for n in sizes:
        seq = random_sequence(n)
        serial = time_call(lambda: bsk.count_kmers(seq, k, canonical=True))
        par2 = time_call(lambda: bsk.count_kmers_parallel([seq], k, canonical=True, workers=2))
        par4 = time_call(lambda: bsk.count_kmers_parallel([seq], k, canonical=True, workers=4))

        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "g.fa"
            write_fasta(seq, fa)
            build = time_call(lambda: bsk.build_faidx(str(fa)))
            idx = bsk.build_faidx(str(fa))
            mid = max(1, n // 2)
            end = min(n, mid + 100)
            n_q = 1000
            fetch_total = time_call(lambda: [idx.fetch("chr1", mid, end) for _ in range(n_q)])
            fetch_per_query_ms = fetch_total / n_q * 1e3

        rows.append(
            {
                "n": n,
                "kmer_serial_s": round(serial, 4),
                "kmer_2w_s": round(par2, 4),
                "kmer_4w_s": round(par4, 4),
                "speedup_4w": round(serial / par4, 2) if par4 else 0.0,
                "faidx_build_s": round(build, 5),
                "fetch_ms_per_query": round(fetch_per_query_ms, 4),
            }
        )
        print(
            f"n={n:>9,d}  serial={serial:7.3f}s  2w={par2:7.3f}s  4w={par4:7.3f}s  "
            f"(x{rows[-1]['speedup_4w']})  build={build:7.4f}s  fetch={fetch_per_query_ms:.4f}ms"
        )
    return rows


def save_csv(rows: list[dict], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {path}")


def save_plot(rows: list[dict], path: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping plot")
        return

    ns = [r["n"] for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    axes[0].plot(ns, [r["kmer_serial_s"] for r in rows], "o-", label="serial")
    axes[0].plot(ns, [r["kmer_2w_s"] for r in rows], "s-", label="2 workers")
    axes[0].plot(ns, [r["kmer_4w_s"] for r in rows], "^-", label="4 workers")
    axes[0].set_xlabel("sequence length (bp)")
    axes[0].set_ylabel("time (s)")
    axes[0].set_title("k-mer counting scalability (k=8)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(ns, [r["faidx_build_s"] for r in rows], "o-", color="seagreen", label="faidx build")
    axes[1].set_xlabel("sequence length (bp)")
    axes[1].set_ylabel("build time (s)", color="seagreen")
    axes[1].set_title("Index build O(n) & O(1) fetch")
    ax2 = axes[1].twinx()
    ax2.plot(ns, [r["fetch_ms_per_query"] for r in rows], "s--", color="indianred", label="fetch/query")
    ax2.set_ylabel("fetch (ms/query)", color="indianred")
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"wrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="smaller sizes for a fast run")
    args = parser.parse_args()

    sizes = [50_000, 200_000, 500_000] if args.quick else [
        100_000,
        500_000,
        1_000_000,
        2_000_000,
        4_000_000,
    ]
    rows = run(sizes)
    save_csv(rows, HERE / "results.csv")
    save_plot(rows, HERE / "scalability.png")


if __name__ == "__main__":
    main()
