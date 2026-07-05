"""Sequence transformations: reverse complement and translation.

Implements DNA reverse complement and the standard genetic code, including
six-frame translation (three forward frames + three reverse-complement frames).
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "reverse_complement",
    "complement",
    "transcribe",
    "back_transcribe",
    "translate",
    "six_frame_translation",
    "Frame",
    "CODON_TABLE",
]

_COMPLEMENT = str.maketrans(
    "ACGTUNacgtunRYSWKMBDHVryswkmbdhv",
    "TGCAANtgcaanYRSWMKVHDByrswmkvhdb",
)

# Standard genetic code (NCBI transl_table=1).
CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def complement(sequence: str) -> str:
    """Return the base-wise complement (IUPAC aware, preserves length)."""
    return sequence.translate(_COMPLEMENT)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return sequence.translate(_COMPLEMENT)[::-1]


_DNA_TO_RNA = str.maketrans("Tt", "Uu")
_RNA_TO_DNA = str.maketrans("Uu", "Tt")


def transcribe(sequence: str) -> str:
    """Transcribe a DNA sequence into RNA by replacing ``T`` with ``U``.

    Case is preserved. This mirrors the coding (sense) strand convention, where
    the RNA has the same sequence as the template's complement.
    """
    return sequence.translate(_DNA_TO_RNA)


def back_transcribe(sequence: str) -> str:
    """Reverse-transcribe an RNA sequence into DNA by replacing ``U`` with ``T``."""
    return sequence.translate(_RNA_TO_DNA)


def translate(sequence: str, unknown: str = "X") -> str:
    """Translate a nucleotide sequence in reading frame 0.

    ``U`` is treated as ``T``. Incomplete trailing codons are dropped.
    Unknown codons map to ``unknown``.
    """
    seq = sequence.upper().replace("U", "T")
    protein: list[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        protein.append(CODON_TABLE.get(codon, unknown))
    return "".join(protein)


@dataclass(frozen=True)
class Frame:
    """One of the six reading frames of a sequence."""

    strand: str  # '+' or '-'
    offset: int  # 0, 1 or 2
    protein: str

    @property
    def name(self) -> str:
        return f"{self.strand}{self.offset + 1}"


def six_frame_translation(sequence: str, unknown: str = "X") -> list[Frame]:
    """Return all six reading-frame translations of ``sequence``.

    Frames ``+1/+2/+3`` translate the forward strand at offsets 0/1/2;
    frames ``-1/-2/-3`` translate the reverse complement at offsets 0/1/2.
    """
    frames: list[Frame] = []
    forward = sequence.upper().replace("U", "T")
    reverse = reverse_complement(forward)
    for strand, seq in (("+", forward), ("-", reverse)):
        for offset in range(3):
            frames.append(Frame(strand, offset, translate(seq[offset:], unknown)))
    return frames
