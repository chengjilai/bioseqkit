# Benchmarks & complexity analysis

Run the driver:

```bash
python benchmarks/bench.py           # full sweep (100 kbp .. 4 Mbp)
python benchmarks/bench.py --quick   # fast sweep
```

It writes [`results.csv`](results.csv) and, with matplotlib installed, the
scalability figure [`scalability.png`](scalability.png).

## Asymptotic complexity

| Operation | Time | Space | Notes |
|---|---|---|---|
| FASTA/FASTQ parse | O(n) | O(1)\* | streaming generator; \*O(longest record) |
| GC / N / composition | O(n) | O(1) | single pass, fixed alphabet |
| reverse complement | O(n) | O(n) | one `str.translate` + reverse |
| six-frame translation | O(n) | O(n) | 6 frames, one codon table |
| k-mer counting | O(n·k) | O(min(4^k, n)) | sliding window; hash table of distinct k-mers |
| parallel k-mer counting | O(n·k / p + m) | O(n) | p workers, m = merge of partial counters |
| minimizer sampling | O(n·w) | O(n) | naive window minimum (O(n) with a monotonic deque) |
| faidx build | O(n) | O(s) | one byte-wise pass; s = number of sequences |
| faidx fetch | O(L) | O(L) | one `seek` + read of an L-base region; independent of file size |

## Measured results

Representative run (random sequence, k = 8 canonical, 4 physical cores):

| n (bp) | serial | 2 workers | 4 workers | speed-up | faidx build | fetch/query |
|---:|---:|---:|---:|---:|---:|---:|
| 100 k | 0.10 s | 0.23 s | 0.10 s | 1.1× | 0.5 ms | 0.016 ms |
| 500 k | 0.52 s | 0.34 s | 0.23 s | 2.3× | 2.0 ms | 0.034 ms |
| 1 M | 0.98 s | 0.66 s | 0.40 s | 2.5× | 3.6 ms | 0.056 ms |
| 2 M | 2.12 s | 1.15 s | 0.66 s | 3.2× | 7.1 ms | 0.100 ms |
| 4 M | 3.85 s | 2.17 s | 1.23 s | 3.1× | 14.3 ms | 0.196 ms |

## Interpretation

- **k-mer counting is linear in n** (time roughly doubles as n doubles),
  matching the expected O(n·k).
- **Parallel scaling improves with input size**: at 100 kbp the process
  start-up and counter-merge overhead dominates (no gain), but from ~1 Mbp the
  4-worker run approaches a 3× speed-up. This is the classic Amdahl/overhead
  regime — parallelism pays off only once the per-chunk work outweighs the
  fixed cost of spawning workers and merging `Counter`s.
- **Index build is linear** (≈3.6 µs/kbp) and the per-query fetch stays in the
  tens-to-hundreds of microseconds, confirming that a query cost is a single
  seek-and-read rather than a scan of the file.
