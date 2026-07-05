"""Sequence file I/O.

Pure-Python, streaming FASTA/FASTQ parsers implemented with the
iterator/generator pattern so that arbitrarily large files can be processed
with constant memory. Both plain-text and gzip-compressed files are supported
transparently (detected by the ``.gz`` extension or the gzip magic bytes).
"""

from __future__ import annotations

import gzip
import io as _io
from dataclasses import dataclass
from typing import IO, Iterable, Iterator

__all__ = [
    "FastaRecord",
    "FastqRecord",
    "open_text",
    "parse_fasta",
    "parse_fastq",
    "write_fasta",
]


@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA record."""

    id: str
    description: str
    sequence: str

    def __len__(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class FastqRecord:
    """A single FASTQ record, carrying Phred quality string."""

    id: str
    description: str
    sequence: str
    quality: str

    def __len__(self) -> int:
        return len(self.sequence)

    def phred_scores(self, offset: int = 33) -> list[int]:
        """Decode the ASCII quality string into integer Phred scores."""
        return [ord(c) - offset for c in self.quality]


def _is_gzip(path: str) -> bool:
    if path.endswith(".gz"):
        return True
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def open_text(path: str) -> IO[str]:
    """Open ``path`` as a text stream, transparently handling gzip files."""
    if _is_gzip(path):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _split_header(header: str) -> tuple[str, str]:
    header = header.strip()
    if not header:
        return "", ""
    parts = header.split(None, 1)
    seq_id = parts[0]
    description = parts[1] if len(parts) > 1 else ""
    return seq_id, description


def parse_fasta(source: str | IO[str]) -> Iterator[FastaRecord]:
    """Stream :class:`FastaRecord` objects from a FASTA file or text stream.

    Blank lines are ignored. Sequence lines are concatenated so multi-line
    records are handled. Raises :class:`ValueError` on malformed input.
    """
    handle, own = _as_handle(source)
    try:
        seq_id: str | None = None
        description = ""
        chunks: list[str] = []
        started = False
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                started = True
                if seq_id is not None:
                    yield FastaRecord(seq_id, description, "".join(chunks))
                seq_id, description = _split_header(line[1:])
                chunks = []
            else:
                if not started:
                    raise ValueError("FASTA sequence data before any '>' header")
                chunks.append(line.strip())
        if seq_id is not None:
            yield FastaRecord(seq_id, description, "".join(chunks))
    finally:
        if own:
            handle.close()


def parse_fastq(source: str | IO[str]) -> Iterator[FastqRecord]:
    """Stream :class:`FastqRecord` objects from a FASTQ file or text stream."""
    handle, own = _as_handle(source)
    try:
        it = iter(handle)
        while True:
            header = _next_nonblank(it)
            if header is None:
                break
            if not header.startswith("@"):
                raise ValueError(f"FASTQ header must start with '@': {header!r}")
            seq_line = _require(it, "sequence")
            plus = _require(it, "'+' separator")
            if not plus.startswith("+"):
                raise ValueError(f"FASTQ separator must start with '+': {plus!r}")
            qual_line = _require(it, "quality")
            if len(seq_line) != len(qual_line):
                raise ValueError(
                    "FASTQ sequence and quality length mismatch "
                    f"({len(seq_line)} vs {len(qual_line)})"
                )
            seq_id, description = _split_header(header[1:])
            yield FastqRecord(seq_id, description, seq_line, qual_line)
    finally:
        if own:
            handle.close()


def write_fasta(records: Iterable[FastaRecord], handle: IO[str], width: int = 70) -> int:
    """Write records to a text handle, wrapping sequences at ``width`` columns.

    Returns the number of records written. ``width <= 0`` disables wrapping.
    """
    n = 0
    for rec in records:
        header = rec.id if not rec.description else f"{rec.id} {rec.description}"
        handle.write(f">{header}\n")
        seq = rec.sequence
        if width and width > 0:
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")
        else:
            handle.write(seq + "\n")
        n += 1
    return n


def _as_handle(source: str | IO[str]) -> tuple[IO[str], bool]:
    if isinstance(source, str):
        return open_text(source), True
    if isinstance(source, (_io.TextIOBase,)) or hasattr(source, "read"):
        return source, False
    raise TypeError(f"Unsupported source type: {type(source)!r}")


def _next_nonblank(it: Iterator[str]) -> str | None:
    for raw in it:
        line = raw.rstrip("\r\n")
        if line.strip():
            return line
    return None


def _require(it: Iterator[str], what: str) -> str:
    try:
        return next(it).rstrip("\r\n")
    except StopIteration as exc:  # noqa: F841
        raise ValueError(f"Truncated FASTQ record: missing {what} line") from None
