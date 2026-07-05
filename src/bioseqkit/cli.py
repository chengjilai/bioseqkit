"""Command-line interface for bioseqkit.

Implemented with :mod:`argparse` (standard library) to keep the core package
dependency-free. Sub-commands: ``stats``, ``revcomp``, ``translate``,
``kmer``, ``minimizer``, ``index`` and ``fetch``.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Sequence

from bioseqkit import __version__
from bioseqkit.index import build_faidx, fetch as fetch_region
from bioseqkit.io import FastaRecord, parse_fasta, write_fasta
from bioseqkit.kmer import count_kmers, count_kmers_parallel, minimizers, top_kmers
from bioseqkit.stats import sequence_stats
from bioseqkit.transform import reverse_complement, six_frame_translation


def _read_sequences(path: str) -> list[FastaRecord]:
    return list(parse_fasta(path))


def cmd_stats(args: argparse.Namespace) -> int:
    records = _read_sequences(args.input)
    stats = sequence_stats(rec.sequence for rec in records)
    print(json.dumps(stats.as_dict(), indent=2))
    return 0


def cmd_revcomp(args: argparse.Namespace) -> int:
    out = [
        FastaRecord(rec.id, (rec.description + " revcomp").strip(), reverse_complement(rec.sequence))
        for rec in _read_sequences(args.input)
    ]
    write_fasta(out, sys.stdout)
    return 0


def cmd_translate(args: argparse.Namespace) -> int:
    out: list[FastaRecord] = []
    for rec in _read_sequences(args.input):
        for frame in six_frame_translation(rec.sequence):
            out.append(FastaRecord(f"{rec.id}_frame{frame.name}", "", frame.protein))
    write_fasta(out, sys.stdout)
    return 0


def cmd_kmer(args: argparse.Namespace) -> int:
    seqs = [rec.sequence for rec in _read_sequences(args.input)]
    if args.threads and args.threads > 1:
        counts = count_kmers_parallel(seqs, args.k, canonical=args.canonical, workers=args.threads)
    else:
        counts = sum(
            (count_kmers(s, args.k, canonical=args.canonical) for s in seqs),
            start=type(count_kmers("", args.k))(),
        )
    for kmer, count in top_kmers(counts, args.top):
        print(f"{kmer}\t{count}")
    return 0


def cmd_minimizer(args: argparse.Namespace) -> int:
    for rec in _read_sequences(args.input):
        for pos, mm in minimizers(rec.sequence, args.k, args.w, canonical=not args.no_canonical):
            print(f"{rec.id}\t{pos}\t{mm}")
    return 0


def cmd_index(args: argparse.Namespace) -> int:
    index = build_faidx(args.input)
    path = index.write()
    print(f"Wrote index: {path} ({len(index.records)} sequences)")
    return 0


def cmd_fetch(args: argparse.Namespace) -> int:
    seq = fetch_region(args.input, args.region)
    name = args.region
    print(f">{name}")
    for i in range(0, len(seq), 70):
        print(seq[i : i + 70])
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bioseqkit", description=__doc__.split("\n")[0])
    parser.add_argument("--version", action="version", version=f"bioseqkit {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("stats", help="sequence statistics (JSON)")
    p.add_argument("input", help="FASTA/FASTA.gz file")
    p.set_defaults(func=cmd_stats)

    p = sub.add_parser("revcomp", help="reverse complement")
    p.add_argument("input")
    p.set_defaults(func=cmd_revcomp)

    p = sub.add_parser("translate", help="six-frame translation")
    p.add_argument("input")
    p.set_defaults(func=cmd_translate)

    p = sub.add_parser("kmer", help="k-mer frequency analysis")
    p.add_argument("input")
    p.add_argument("-k", type=int, default=5, help="k-mer size (default 5)")
    p.add_argument("--top", type=int, default=10, help="report top-N k-mers")
    p.add_argument("--canonical", action="store_true", help="merge reverse-complement k-mers")
    p.add_argument("-t", "--threads", type=int, default=1, help="parallel worker processes")
    p.set_defaults(func=cmd_kmer)

    p = sub.add_parser("minimizer", help="minimizer sampling")
    p.add_argument("input")
    p.add_argument("-k", type=int, default=15, help="k-mer size (default 15)")
    p.add_argument("-w", type=int, default=10, help="window size (default 10)")
    p.add_argument("--no-canonical", action="store_true", help="disable canonical k-mers")
    p.set_defaults(func=cmd_minimizer)

    p = sub.add_parser("index", help="build a FAI-like index")
    p.add_argument("input")
    p.set_defaults(func=cmd_index)

    p = sub.add_parser("fetch", help="fetch a sub-sequence (chr:start-end)")
    p.add_argument("input")
    p.add_argument("region", help="region string, e.g. chr1:1000-2000")
    p.set_defaults(func=cmd_fetch)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
