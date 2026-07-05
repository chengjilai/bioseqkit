# bioseqkit

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

Continuous integration (GitHub Actions) runs `ruff` linting and the `pytest`
suite on Python 3.10–3.12 for every push.

## Data sources

- NCBI Nucleotide: <https://www.ncbi.nlm.nih.gov/nucleotide/>
- UCSC Genome Browser: <https://genome.ucsc.edu/>

The bundled `examples/example_data/sample.fa` is a small synthetic sequence for
offline testing; `demo.ipynb` will download real data from NCBI when a network
connection is available and fall back to the bundled file otherwise.

## License

MIT — see [LICENSE](LICENSE).
