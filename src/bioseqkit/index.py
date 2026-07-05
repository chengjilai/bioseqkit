"""FAI-like FASTA indexing for random access.

Mirrors the ``samtools faidx`` (``*.fai``) format so that arbitrary
sub-sequences can be fetched without reading the whole file. Each index line
holds: name, sequence length, byte offset of the first base, number of bases
per line and number of bytes per line (including the newline).
"""

from __future__ import annotations

import os
from dataclasses import dataclass

__all__ = ["FaidxRecord", "FaidxIndex", "build_faidx", "fetch", "parse_region"]


@dataclass(frozen=True)
class FaidxRecord:
    name: str
    length: int
    offset: int
    linebases: int
    linewidth: int

    def to_line(self) -> str:
        return f"{self.name}\t{self.length}\t{self.offset}\t{self.linebases}\t{self.linewidth}"

    @classmethod
    def from_line(cls, line: str) -> "FaidxRecord":
        name, length, offset, linebases, linewidth = line.rstrip("\n").split("\t")
        return cls(name, int(length), int(offset), int(linebases), int(linewidth))


class FaidxIndex:
    """An in-memory FASTA index bound to a plain-text FASTA file."""

    def __init__(self, fasta_path: str, records: dict[str, FaidxRecord]):
        self.fasta_path = fasta_path
        self.records = records

    def names(self) -> list[str]:
        return list(self.records)

    def write(self, fai_path: str | None = None) -> str:
        path = fai_path or self.fasta_path + ".fai"
        with open(path, "w") as fh:
            for rec in self.records.values():
                fh.write(rec.to_line() + "\n")
        return path

    @classmethod
    def load(cls, fasta_path: str, fai_path: str | None = None) -> "FaidxIndex":
        path = fai_path or fasta_path + ".fai"
        records: dict[str, FaidxRecord] = {}
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    rec = FaidxRecord.from_line(line)
                    records[rec.name] = rec
        return cls(fasta_path, records)

    def fetch(self, name: str, start: int | None = None, end: int | None = None) -> str:
        """Fetch a sub-sequence using 1-based inclusive ``start``/``end``.

        With no coordinates the whole sequence is returned.
        """
        if name not in self.records:
            raise KeyError(f"Sequence {name!r} not found in index")
        rec = self.records[name]
        s = 0 if start is None else start - 1
        e = rec.length if end is None else end
        if s < 0 or e > rec.length or s > e:
            raise ValueError(f"Region out of bounds for {name} (length {rec.length})")
        newline_bytes = rec.linewidth - rec.linebases
        with open(self.fasta_path, "rb") as fh:
            start_line, start_col = divmod(s, rec.linebases)
            byte_start = rec.offset + start_line * rec.linewidth + start_col
            fh.seek(byte_start)
            n_bases = e - s
            # Read enough bytes to cover the requested bases plus newlines.
            n_lines = (start_col + n_bases) // rec.linebases
            n_read = n_bases + n_lines * newline_bytes + rec.linewidth
            raw = fh.read(n_read)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii")
        return seq[:n_bases]


def build_faidx(fasta_path: str) -> FaidxIndex:
    """Scan a plain-text FASTA file and build a :class:`FaidxIndex`.

    Raises :class:`ValueError` for gzip files (which are not seekable by base)
    or for records with inconsistent line lengths.
    """
    if fasta_path.endswith(".gz"):
        raise ValueError("faidx requires an uncompressed FASTA file")
    records: dict[str, FaidxRecord] = {}
    with open(fasta_path, "rb") as fh:
        name: str | None = None
        length = 0
        offset = 0
        linebases = 0
        linewidth = 0
        line_lengths: list[int] = []
        pos = 0

        def flush() -> None:
            nonlocal name
            if name is not None:
                _validate_lines(name, line_lengths)
                records[name] = FaidxRecord(name, length, offset, linebases, linewidth)

        for raw in fh:
            if raw.startswith(b">"):
                flush()
                header = raw[1:].decode("ascii", "replace").strip()
                name = header.split()[0] if header else ""
                length = 0
                linebases = 0
                linewidth = 0
                line_lengths = []
                offset = pos + len(raw)
            else:
                stripped = raw.rstrip(b"\r\n")
                if name is not None:
                    if linebases == 0:
                        linebases = len(stripped)
                        linewidth = len(raw)
                    line_lengths.append(len(stripped))
                    length += len(stripped)
            pos += len(raw)
        flush()
    return FaidxIndex(fasta_path, records)


def _validate_lines(name: str, line_lengths: list[int]) -> None:
    if len(line_lengths) > 1:
        body = line_lengths[:-1]
        if len(set(body)) > 1:
            raise ValueError(f"Inconsistent line lengths in record {name!r}; cannot index")
        if line_lengths[-1] > body[0]:
            raise ValueError(f"Last line longer than others in record {name!r}")


def parse_region(region: str) -> tuple[str, int | None, int | None]:
    """Parse a ``chr:start-end`` region string (1-based, inclusive).

    ``chr`` alone returns the whole sequence.
    """
    if ":" not in region:
        return region, None, None
    name, span = region.rsplit(":", 1)
    span = span.replace(",", "")
    if "-" in span:
        start_s, end_s = span.split("-", 1)
        return name, int(start_s), int(end_s)
    return name, int(span), int(span)


def fetch(fasta_path: str, region: str) -> str:
    """Convenience: build/load index and fetch a region in one call."""
    fai_path = fasta_path + ".fai"
    index = FaidxIndex.load(fasta_path) if os.path.exists(fai_path) else build_faidx(fasta_path)
    name, start, end = parse_region(region)
    return index.fetch(name, start, end)
