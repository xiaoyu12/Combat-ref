# Validation report — 0.1.0

## Status

This release has deterministic unit and integration coverage, but it has not
passed a broad simulation or external real-data validation gate. It must not be
described as validated, superior, or production ready.

## Reproducible smoke simulation

`validation/simulation.R` generates 200 NB genes across two balanced batches
and two biological groups with seed 104729. Batch B has a 1.8-fold simulated
mean multiplier; half the genes have a two-fold group effect. The reference is
fixed to batch A.

| Metric | Result |
|---|---:|
| Raw median B/A mean ratio | 1.7112 |
| Adjusted median B/A mean ratio | 1.6447 |
| Unsupported genes | 0 |

The correction moved the simple batch ratio toward one only modestly. This
single smoke scenario is insufficient for type-I error, FDR, power, signal
preservation, reference stability, runtime, or memory claims.

## Gate for 0.2.0

Before 0.2.0, run the predefined multi-scenario study covering 2–5 batches,
balanced/unbalanced groups, depth/composition shifts, zeros, outliers,
confounding failures, null effects, and sparse/widespread DE. Compare direct
edgeR QL analysis of raw counts, legacy, repaired baseline, and QL-assisted
diagnostics. Predefine acceptable false-positive, FDR, bias, and biological
signal-preservation margins before reading results.
