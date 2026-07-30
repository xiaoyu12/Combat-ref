# Changelog

## 0.1.0 pre-release corrections

| Area | Previous behavior | New behavior | Consequence | Files/tests |
|---|---|---|---|---|
| Local dispersion design | Removed batch columns and could omit an intercept | Dedicated intercept-bearing, reduced estimable design | Valid baselines and auditable fallbacks | `build_design.R`, statistical follow-up tests |
| QL NB dispersion | Supplied tagwise dispersion | edgeR v4 estimates one global common dispersion | Correct QL variance decomposition | `combat_ref_baseline.R`, QL regression test |
| Reference diagnostics | One-hot vector called frequency | Frequency `NULL`; selected batch/method named | No fabricated uncertainty | `combat_ref.R` |
| Convergence | Supported implied converged | edgeR/numerical checks or `NA` | Honest gene status | `combat_ref.R` |
| Count storage | Forced 32-bit integer | Integer-valued numeric double | No overflow above 2^31-1 | `quantile_mapping.R`, overflow test |
| Input handling | Mostly rejected coercible inputs | Safe repairs with structured audit | Reproducible conversion history | `prepare_inputs.R`, preparation tests |
| Lossy repair | Unavailable | Explicit round/sum policies | Never silent | preparation tests |
| CLI | Minimal success messages | Structured default report and machine outputs | Useful interactive and pipeline modes | CLI script/tests |
