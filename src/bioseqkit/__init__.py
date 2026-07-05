"""bioseqkit: a lightweight, dependency-free biological sequence toolkit.

The package provides pure-Python FASTA/FASTQ parsing, sequence statistics,
transformations (reverse complement, six-frame translation), k-mer analysis,
minimizer sampling and FAI-like random-access indexing.
"""

from bioseqkit.io import (
    FastaRecord,
    FastqRecord,
    parse_fasta,
    parse_fastq,
    write_fasta,
)
from bioseqkit.stats import (
    SeqStats,
    base_composition,
    gc_content,
    n_ratio,
    sequence_stats,
)
from bioseqkit.transform import (
    back_transcribe,
    reverse_complement,
    six_frame_translation,
    transcribe,
    translate,
)
from bioseqkit.kmer import (
    canonical_kmer,
    count_kmers,
    count_kmers_parallel,
    minimizers,
    top_kmers,
)
from bioseqkit.index import FaidxIndex, build_faidx, fetch

__version__ = "0.1.0"

__all__ = [
    "FastaRecord",
    "FastqRecord",
    "parse_fasta",
    "parse_fastq",
    "write_fasta",
    "SeqStats",
    "sequence_stats",
    "gc_content",
    "n_ratio",
    "base_composition",
    "reverse_complement",
    "translate",
    "transcribe",
    "back_transcribe",
    "six_frame_translation",
    "count_kmers",
    "count_kmers_parallel",
    "top_kmers",
    "canonical_kmer",
    "minimizers",
    "FaidxIndex",
    "build_faidx",
    "fetch",
    "__version__",
]
