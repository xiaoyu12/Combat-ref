# Optional benchmark

Run from the repository root after installing `combatrefql` and optional
Bioconductor package `sva`:

```r
source("benchmark/benchmark.R")
```

On 2026-08-04 with seed 20260804, R 4.6.1, edgeR 4.10.1, and sva 3.60.0,
identical 2,000-gene by 24-sample input with three batches produced:

| Method | Elapsed seconds | Approximate result size (MiB) |
|---|---:|---:|
| combatrefql | 0.379 | 1.414 |
| sva::ComBat_seq | 0.422 | 0.506 |

Elapsed time is one local run, not a stable estimate. Result object size is not
peak memory. The script reports both transparently and makes no superiority
claim.
