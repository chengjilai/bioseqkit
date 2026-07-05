# bioseqkit

[![CI](https://github.com/chengjilai/bioseqkit/actions/workflows/ci.yml/badge.svg)](https://github.com/chengjilai/bioseqkit/actions/workflows/ci.yml)
[![Docs](https://github.com/chengjilai/bioseqkit/actions/workflows/docs.yml/badge.svg)](https://chengjilai.github.io/bioseqkit/)
[![PyPI](https://img.shields.io/pypi/v/bioseqkit.svg)](https://pypi.org/project/bioseqkit/)
[![Python](https://img.shields.io/pypi/pyversions/bioseqkit.svg)](https://pypi.org/project/bioseqkit/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A lightweight, **dependency-free** biological sequence processing toolkit built
from scratch in pure Python. `bioseqkit` implements FASTA/FASTQ parsing,
sequence statistics, transformations, k-mer / minimizer analysis and FAI-like
random-access indexing, exposed both as a Python API and a command-line tool.

The project is a teaching implementation for **BIO2502 (Programming Languages
for Biological Computing)**: it deliberately re-implements the low-level I/O,
streaming and indexing logic instead of relying on Biopython, so the core
design patterns of bioinformatics data handling are made explicit.

## Features

- **Streaming FASTA/FASTQ parsers** (`io`) — generator based, constant memory,
  transparent gzip support, Phred quality decoding.
- **Statistics** (`stats`) — length distribution, N50, GC content, N-base
  ratio, base-composition matrix.
- **Transformations** (`transform`) — reverse complement (IUPAC aware) and
  six-frame translation with the standard genetic code.
- **k-mer analysis** (`kmer`) — counting, top-k, canonical k-mers,
  **multi-process** parallel counting, and **minimizer** sampling.
- **FAI-like indexing** (`index`) — `samtools faidx`-compatible index for
  `chr:start-end` random access without scanning the whole file.
- **CLI** (`cli`) — `stats`, `revcomp`, `translate`, `kmer`, `minimizer`,
  `index`, `fetch`.
- **NCBI download** (`entrez`) — fetch reference sequences via E-utilities
  (standard-library HTTP only).

## Project layout

```
bioseqkit/
├── pyproject.toml          # src-layout, PEP 621 metadata, console script
├── README.md
├── LICENSE
├── src/bioseqkit/
│   ├── __init__.py         # public API
│   ├── io.py               # FASTA/FASTQ parsers
│   ├── stats.py            # sequence statistics
│   ├── transform.py        # revcomp + six-frame translation
│   ├── kmer.py             # k-mer / minimizer (serial + parallel)
│   ├── index.py            # FAI-like random-access index
│   ├── entrez.py           # NCBI download helper
│   └── cli.py              # argparse CLI
├── tests/                  # pytest suite (io/stats/transform/kmer/index/cli)
├── benchmarks/             # benchmark driver + complexity analysis
├── workflow/               # Snakemake pipeline (Snakefile)
├── config/config.yaml      # workflow configuration
├── examples/
│   ├── demo.ipynb          # Jupyter demo (stats, GC, k-mer spectrum, ...)
│   └── example_data/sample.fa
├── docs/                   # Sphinx documentation
├── flake.nix               # Nix flake: reproducible dev shell + Docker/Apptainer images
├── flake.lock              # pinned nixpkgs revision (bit-reproducible builds)
└── .github/workflows/ci.yml
```

## Installation

Requires Python >= 3.10. The core package has **no runtime dependencies**.

```bash
pip install bioseqkit

# with optional extras (plots for the notebook / NCBI download / docs)
pip install "bioseqkit[viz,net,docs]"
```

### Reproducible builds with Nix

A [Nix flake](flake.nix) is provided as a single, bit-for-bit reproducible
source of truth. Because `flake.lock` pins the exact `nixpkgs` revision, every
build — the package, the dev shell and the container images — is fully
reproducible. From the one flake you can:

```bash
nix run   .                 # run the bioseqkit CLI
nix build .#default         # build the Python package (runs the test suite)
nix develop                 # enter a reproducible dev shell (python + uv + typst + pytest ...)

nix build .#docker          # build an OCI/Docker image  -> docker load < result
nix build .#apptainer       # build an Apptainer/Singularity (.sif) image
```

This supersedes a hand-written Dockerfile/Apptainer recipe: both container
images are *derived* from the same pinned dependency graph, so they can never
drift from the tested build.

## Command-line usage

```bash
bioseqkit stats    examples/example_data/sample.fa      # JSON statistics
bioseqkit revcomp  examples/example_data/sample.fa      # reverse complement
bioseqkit translate examples/example_data/sample.fa     # six-frame translation
bioseqkit kmer     examples/example_data/sample.fa -k 5 --top 10 --canonical
bioseqkit kmer     examples/example_data/sample.fa -k 5 -t 4   # parallel
bioseqkit minimizer examples/example_data/sample.fa -k 15 -w 10
bioseqkit index    examples/example_data/sample.fa      # write *.fai
bioseqkit fetch    examples/example_data/sample.fa seq2:1-16
```

## Understanding the output

**`stats`** prints a JSON object; each field means:

| Field | Meaning |
|---|---|
| `n_seqs` | number of sequences in the file |
| `total_length` | sum of all sequence lengths (bp) |
| `min_length` / `max_length` | shortest / longest sequence |
| `mean_length` | average sequence length |
| `n50` | length such that sequences ≥ this length cover ≥ 50% of `total_length` (assembly contiguity metric) |
| `gc_content` | fraction of G/C bases in `[0, 1]` — useful for species/GC-bias assessment |
| `n_ratio` | fraction of ambiguous `N` bases — a data-quality indicator |
| `base_counts` | per-base counts (`A/C/G/T/N/...`), the base-composition matrix |

**`kmer`** prints a tab-separated `kmer<TAB>count` table, one line per k-mer,
sorted by descending frequency (top-`--top`). With `--canonical`, a k-mer and
its reverse complement are counted together, so results are strand-independent.

**`revcomp` / `translate`** emit FASTA to stdout. `translate` produces six
records per input sequence, suffixed `_frame+1..+3` (forward strand, offsets
0/1/2) and `_frame-1..-3` (reverse-complement strand); `*` denotes a stop codon
and `X` an untranslatable codon.

**`minimizer`** prints `seq_id<TAB>position<TAB>minimizer` — the sampled k-mers
and their 0-based positions along the sequence.

**`index`** writes a `<file>.fai` (name, length, offset, bases-per-line,
bytes-per-line); **`fetch`** prints the requested sub-sequence as FASTA using
1-based inclusive coordinates.

## Reproducible pipeline (Snakemake)

A [Snakemake](workflow/Snakefile) workflow chains the CLI into a
"data acquisition → processing → analysis" pipeline
(`fetch_input → stats + kmer + index`):

```bash
snakemake -c1                                   # run offline on the bundled example
snakemake -n                                    # dry run: show the job DAG
snakemake -c4 --config accession=NC_012920.1 k=6  # download human chrM from NCBI, then analyse
```

Outputs are written to `results/` (`stats.json`, `kmers.tsv`, `input.fa.fai`).
Configure input, accession and parameters in [`config/config.yaml`](config/config.yaml).

## Benchmarks & complexity

```bash
python benchmarks/bench.py            # sweep 100 kbp .. 4 Mbp, write results.csv + scalability.png
```

See [`benchmarks/README.md`](benchmarks/README.md) for the full time/space
complexity table and measured scaling. In short: k-mer counting is O(n·k) and
scales to ~3× on four workers for inputs ≥ 1 Mbp; the FAI index builds in O(n)
and answers random `fetch` queries in tens of microseconds, independent of file
size. Micro-benchmarks are also runnable via `pytest tests/test_benchmark.py
--benchmark-only`.

## Python API

```python
import bioseqkit as bsk

for rec in bsk.parse_fasta("examples/example_data/sample.fa"):
    print(rec.id, len(rec), bsk.gc_content(rec.sequence))

print(bsk.reverse_complement("ATGC"))            # -> GCAT
print(bsk.translate("ATGGCCTAA"))                # -> MA*

counts = bsk.count_kmers("ACGTACGTACGT", k=3, canonical=True)
print(bsk.top_kmers(counts, 3))

idx = bsk.build_faidx("examples/example_data/sample.fa")
print(idx.fetch("seq2", 1, 16))
```

## Testing

```bash
uv run --with pytest pytest -q      # 39 tests
```

Continuous integration (GitHub Actions) runs `ruff` linting, the `pytest`
suite with coverage on Python 3.10–3.12 for every push, and builds & deploys
the Sphinx documentation to GitHub Pages.

## Data sources

- NCBI Nucleotide: <https://www.ncbi.nlm.nih.gov/nucleotide/>
- UCSC Genome Browser: <https://genome.ucsc.edu/>

The bundled `examples/example_data/sample.fa` is a small synthetic sequence for
offline testing; `demo.ipynb` will download real data from NCBI when a network
connection is available and fall back to the bundled file otherwise.

## License

MIT — see [LICENSE](LICENSE).
