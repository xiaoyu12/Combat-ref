# Validation

Validation is separate from package code and ordinary `R CMD check`. Unit tests
remain offline and fast. Public datasets and comparator source files are cached
under ignored `validation/.cache/`.

## Simulation

```sh
Rscript validation/simulation.R fast
Rscript validation/simulation.R full
```

The fast run evaluates two deterministic scenarios. The full run evaluates 144
scenarios and adds ComBat-ref and ComBat-seq from pinned public commits. Raw
counts analyzed with an edgeR batch-and-biology design provide the inference
baseline. Outputs are `results-fast.tsv` and `results-full.tsv`.

## Real data

```sh
Rscript validation/real_data/gfrn/run.R
Rscript validation/real_data/nasa/run.R
```

The GFRN workflow downloads the `signature_data.rds` input used by ComBat-ref,
selects GFP controls and HER2, EGFR, and KRAS perturbations, and applies the
upstream per-batch expression filter. The NASA workflow downloads the upstream
`nasa.RData`, excludes GLDS-168 as in the published script, and analyzes the
remaining GeneLab liver samples. Each workflow writes:

- `metrics.tsv`: quantitative PCA, mean, dispersion, and invariance metrics;
- `pca-scores.tsv`: tidy coordinates for reproducible plotting;
- `pca.svg`: a compact vector before/after comparison.

Inputs are pinned to ComBat-ref commit
`f0c6d3f313e9f1f29de2c8c8f12b98113443be55`. ComBat-seq is pinned to
`2673dd8f23052ea4354ffa5e7600fabcd3f08269`. The scripts do not commit the
downloaded data.

See [VALIDATION_PLAN_0.4.0.md](VALIDATION_PLAN_0.4.0.md) for the preregistered
design and [VALIDATION_REPORT.md](VALIDATION_REPORT.md) for results and
interpretation.
