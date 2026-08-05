# Validation resources

`simulation.R` runs the fast subset or the full 0.4.0 factorial study. The plan
fixes scenarios and metrics; the report summarizes the resulting tidy TSV
files. Full runs fetch the two upstream comparator implementations from pinned
public commits into a temporary directory.

```sh
Rscript validation/simulation.R fast
Rscript validation/simulation.R full
```
