# ComBat-refQL

Reference-batch correction for raw RNA-seq counts. Version 0.1.0 uses TMMwsp
normalization, negative-binomial counterfactual prediction, and deterministic
mid-P mapping. The optional QL pathway is diagnostic and experimental.

## R API

```r
fit <- combat_ref(
  counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference_batch = "batch_A"
)
fit$adjusted_counts
fit$diagnostics$input_actions
```

`verbosity = "auto"` is quiet in scripts and normal interactively. Choose
`quiet`, `normal`, or `verbose` explicitly when needed. The older `verbose`
logical remains accepted.

### Input preparation

Lossless repairs are automatic, announced in normal mode, and recorded in
`diagnostics$input_actions`: data frames to numeric matrices; numeric character
columns; a unique first gene-ID column; integer storage to double; batch/group
factors; one-level group to `NULL`; named sample reordering; unambiguous
transposition; and case-insensitive method spellings.

Lossy changes require explicit opt-in:

```r
combat_ref(..., fractional_counts = "round")
combat_ref(..., duplicate_genes = "sum")
```

Negative, missing, non-finite, unparseable, ambiguous, normalized, or undeclared
fractional inputs remain errors. Declare `input_scale = "TPM"` (or CPM, logCPM,
RPKM, FPKM) to obtain an explicit rejection; only `"counts"` is accepted.

## CLI

```bash
combat-ref \
  --counts counts.tsv \
  --metadata metadata.csv \
  --batch-column batch \
  --group-column condition \
  --output adjusted.tsv \
  --diagnostics-json diagnostics.json \
  --gene-status gene_status.tsv
```

Normal mode prints input, configuration, design, normalization, reference,
correction, warnings, and output summaries. Add `--verbose` for timings,
`--quiet` for pipelines, `--no-color` for plain text, `--log-file` for a run
log, or `--validate-only` to check inputs/design without fitting.

Stable exit codes are 0 success, 2 input/format, 3 design, 4 fit/mapping, and 5
output writing. TSV, CSV, clearly delimited TXT, and gzipped TSV/CSV are
supported.

## Statistical details

Every within-batch dispersion score uses an intercept-bearing local design.
Missing levels and non-estimable columns are reduced explicitly and recorded in
`reference_local_designs`. Automatic reference ties are resolved by batch name.

The ordinary fit uses tagwise NB dispersion. With edgeR >= 4.2.0, QL fitting
uses `glmQLFit(..., dispersion = NULL, legacy = FALSE)`, which estimates a
global common NB dispersion. The global QL input, tagwise baseline dispersion,
and moderated QL dispersion are stored separately. QL dispersion is never used
in NB quantile functions.

For differential expression, analyze raw counts directly with an edgeR QL
model containing batch and biological group; corrected counts do not propagate
correction uncertainty.

## Citation and license

The original ComBat-ref method is described by Zhang X. “Highly effective batch
effect correction method for RNA-seq count data.” *Computational and Structural
Biotechnology Journal* (2025). DOI: 10.1016/j.csbj.2024.12.010; PMID: 39802213.

The newly written package implementation is MIT licensed. This does not invent
or alter licensing terms for upstream code or data.
