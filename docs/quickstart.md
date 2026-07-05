# Quickstart

Install the package and run the CLI:

```bash
pip install bioseqkit
bioseqkit stats examples/example_data/sample.fa
```

Use the Python API:

```python
import bioseqkit as bsk

records = list(bsk.parse_fasta("examples/example_data/sample.fa"))
stats = bsk.sequence_stats(r.sequence for r in records)
print(stats.as_dict())

print(bsk.reverse_complement("ATGC"))          # GCAT
print(bsk.count_kmers("ACGTACGT", 3).most_common(3))
```
