# ComBat-refQL 0.4.0 validation report

## Simulation

The fast subset produced 4 scenario-method rows with no errors. The full study
produced 576 rows: 144 scenarios evaluated with ComBat-refQL, raw-count edgeR,
pinned ComBat-ref, and pinned ComBat-seq. Every scenario completed.

| Method | Mean delta RMSE | Mean log-dispersion RMSE | Biology LFC RMSE | Mean absolute change |
|---|---:|---:|---:|---:|
| ComBat-refQL | 0.200 | 0.200 | 0.566 | 8.52 |
| ComBat-ref | — | — | 0.517 | 12.25 |
| ComBat-seq | — | — | 0.580 | 12.85 |
| Raw edgeR design input | — | — | 0.575 | 0.00 |

For ComBat-refQL, mean absolute log-dispersion error was 0.332, 0.161, and
0.107 for true source/reference ratios 0.5, 1, and 2. The simulations support
the implemented uncertainty moderation and dispersion-ratio recovery in these
scenarios; they do not establish general superiority.

Raw-count edgeR design modelling was the only workflow evaluated for
inferential FDR and power. Those scenario-level values remain in the tidy full
results.

## GFRN

The analysis retained 14,724 genes and 47 samples. The ComBat-refQL reference
was batch 3 and its observed counts remained exactly unchanged.

| Method | Conditional batch R2 | Mean-alignment RMSE | Dispersion-alignment score |
|---|---:|---:|---:|
| Unadjusted | 0.549 | 1.191 | 0.270 |
| ComBat-refQL | 0.102 | 0.127 | 0.178 |
| ComBat-ref | 0.001 | 0.073 | 0.279 |
| ComBat-seq | 0.194 | 0.221 | 0.470 |

## NASA GeneLab

The upstream GLDS-168 exclusion left 18,112 genes and 35 samples. ComBat-refQL
selected GLDS-137/RR3 as reference and preserved its counts exactly.

| Method | Conditional batch R2 | Mean-alignment RMSE | Dispersion-alignment score |
|---|---:|---:|---:|
| Unadjusted | 0.839 | 1.382 | 0.530 |
| ComBat-refQL | 0.224 | 0.288 | 0.114 |
| ComBat-ref | 0.067 | 0.208 | 0.140 |
| ComBat-seq | 0.086 | 0.313 | 0.154 |

Lower values indicate less residual batch association or closer cross-batch
alignment for the named metric. The datasets lack complete ground truth, so
these results describe empirical structure rather than correctness or method
ranking. Stronger batch removal alone does not demonstrate better preservation
of biology.

## Reproducibility

```sh
Rscript validation/simulation.R fast
Rscript validation/simulation.R full
Rscript validation/real_data/gfrn/run.R
Rscript validation/real_data/nasa/run.R
```

The implementation and validation environment is recorded in
`validation/session-info.txt`.
