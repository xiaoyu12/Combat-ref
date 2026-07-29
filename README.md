# ComBat-refQL

**Development status:** ComBat-refQL is an in-progress refactoring and extension of the original ComBat-ref research implementation. The quasi-likelihood extension is not yet complete, validated, or packaged for routine use.

## Overview

ComBat-refQL is being developed as an R package for reference-batch RNA-seq batch correction. The project starts from a script-oriented ComBat-ref research codebase and aims to provide a documented, tested, installable library that can be called directly from an R session.

The long-term goal is to preserve the intended statistical behavior of ComBat-ref while investigating whether parts of the algorithm can benefit from the edgeR quasi-likelihood framework. Any corrected counts produced by this project must remain suitable for the downstream RNA-seq analyses for which ComBat-ref-style adjustment is intended.

## Motivation

RNA-seq batch correction must separate technical batch effects from biological signal without producing corrected counts that are inconsistent with downstream count-based analysis. The current repository is useful research code, but it is organized around terminal scripts, global file paths, simulation workflows, and analysis-specific outputs.

ComBat-refQL is intended to make the method easier to use and validate by:

- exposing a small R-facing API;
- preserving sample and gene identifiers through correction;
- separating reusable statistical code from simulations and real-data analyses;
- establishing regression tests before changing the original algorithm;
- documenting where negative-binomial modeling, dispersion estimation, quasi-likelihood fitting, and count adjustment enter the workflow.

## Relationship to ComBat-ref

This repository is derived from the original ComBat-ref implementation by Xiaoyu Zhang and collaborators. The original method, publication, authorship, citation, and acknowledgements are preserved below in [References and Acknowledgements](#references-and-acknowledgements).

ComBat-refQL is not presented as the official ComBat-ref implementation. It is a derivative development repository that is being refactored toward an R package and used to investigate a quasi-likelihood extension.

## Proposed Quasi-Likelihood Extension

edgeR quasi-likelihood modeling is primarily an inference framework. It improves uncertainty modeling for differential expression by estimating quasi-likelihood variability around negative-binomial mean-variance relationships. It does not, by itself, define a batch-adjusted count matrix.

For this project, the following pieces must remain distinct:

1. **Negative-binomial model fitting:** modeling RNA-seq counts with mean and dispersion parameters.
2. **Dispersion estimation:** estimating biological or technical variability used by the negative-binomial model.
3. **Quasi-likelihood fitting and uncertainty estimation:** estimating additional uncertainty for inference, for example through edgeR-style QL methods.
4. **Count transformation or adjustment:** the explicit algorithm that maps observed counts to corrected counts while preserving the intended reference-batch behavior.

Replacing `glmFit()` with `glmQLFit()` alone does not define a valid count-adjustment method. The exact role of quasi-likelihood quantities in the ComBat-refQL correction algorithm must be formally specified, implemented as a separate pathway, and validated against baseline ComBat-ref behavior and simulation targets.

## Planned R-Package Interface

The intended package name is `combatrefql`, and the main planned user-facing function is `combat_ref_ql()`.

Planned interface, not yet implemented:

```r
library(combatrefql)

adjusted <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference_batch = "batch1"
)
```

The package interface should eventually accept count matrices and metadata vectors directly, rather than requiring command-line arguments and fixed output files.

## Installation

ComBat-refQL is not yet structured as an installable R package. For now, this repository should be treated as a development workspace.

Planned development installation, once package metadata exists:

```r
# install.packages("pak")
pak::pak("Thokas99/ComBat-refQL")
```

This command is provisional and should not be considered supported until `DESCRIPTION`, `NAMESPACE`, documentation, and tests are in place.

## Minimal Usage Example

Current script-oriented workflow:

```sh
Rscript batch_correct.R gfrn.mat gfrn_samples.csv
```

This reads a count matrix and sample metadata file, runs the current ComBat-ref correction code, and writes a corrected matrix. This is legacy behavior from the research implementation and will be retained as baseline behavior until regression tests exist.

Planned R workflow:

```r
adjusted <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference_batch = "batch1"
)
```

## Development Roadmap

1. Establish baseline behavior and regression tests for the current implementation.
2. Create a minimal R-package skeleton.
3. Extract the original ComBat-ref algorithm into documented functions without changing results.
4. Replace script-level global state with explicit arguments and structured return values.
5. Add input validation for counts, metadata, reference batches, and design matrices.
6. Implement the QL method as a separate internal pathway, not a mechanical replacement.
7. Validate original and QL implementations side by side on small simulations and known fixtures.
8. Optimize performance only after correctness and reproducibility are established.
9. Add vignettes, package checks, continuous integration, and release documentation.

## Reproducibility and Validation

Before any statistical changes are made, the current script output should be captured with small reproducible fixtures. Validation should include:

- preservation of count matrix dimensions and dimnames;
- validation of batch, group, and reference-batch inputs;
- deterministic behavior on fixed simulated datasets;
- equivalence with the original implementation where expected;
- behavior with zero-count genes, low-count genes, missing values, and confounded designs;
- downstream checks for type-I error, power, bias, false-positive rate, and batch-removal metrics.

Claims about improved robustness, accuracy, power, or false-positive control require explicit validation and should not be inferred from using edgeR QL machinery alone.

## License

No repository license file is currently present in this checkout. Licensing and redistribution terms for derivative package work need manual review before release or publication.

## References and Acknowledgements

This repository is derived from the original ComBat-ref research implementation. Upstream attribution is retained here prominently because the original method and results are not the work of this derivative package project.

Original ComBat-ref publication:

- Zhang X. and collaborators. ComBat-ref for RNA-seq batch correction. PubMed: <https://pubmed.ncbi.nlm.nih.gov/39802213/>

Original author acknowledgement retained from the upstream README:

- Xiaoyu Zhang

Additional citations for edgeR and related statistical methods will be added as the package documentation is formalized.
