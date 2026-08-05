# ComBat-refQL 0.4.0 statistical specification

## Status

This document defines the intended statistical model for ComBat-refQL 0.4.0.
It separates the mathematical decisions from their R implementation. The
implemented package documentation must record any justified deviation from
this specification.

## 1. Objective

ComBat-refQL aligns supported non-reference RNA-seq count distributions to an
observed reference batch while retaining the reference counts. Version 0.4.0
combines:

- QL-aware empirical-Bayes moderation of gene-specific batch mean contrasts;
- hierarchical batch-specific NB dispersion estimation;
- deterministic source-to-reference NB mid-P transport.

The corrected matrix supports batch-aligned exploration, visualization,
integration, and count-scale downstream operations. Differential-expression
inference remains anchored to the original counts with batch and biological
variables represented in the statistical design.

## 2. Estimand

The method estimates gene-specific relative-abundance batch effects after
TMMwsp library/composition normalization.

A multiplier shared by nearly all genes is structurally difficult to separate
from library size after compositional normalization. Absolute global RNA-output
shifts require an external anchor such as spike-ins. ComBat-refQL therefore
targets residual gene-specific mean effects and batch-associated dispersion
differences on the normalized count scale.

## 3. Mean model

For gene \(g=1,\ldots,G\) and sample \(j=1,\ldots,n\), let

\[
Y_{gj}
\]

denote the observed raw count. The model mean satisfies

\[
\log\mu_{gj}
=
\log L_j+X_j\beta_g+\gamma_{g,b(j)},
\]

where:

- \(L_j\) is the TMMwsp effective library size;
- \(X_j\) contains biological groups and supplied covariates;
- \(\beta_g\) contains biological/covariate coefficients;
- \(b(j)\) is the observed batch;
- \(\gamma_{g,b(j)}\) is the gene-specific batch coefficient.

The design must be full rank and retain positive global residual information.
Correction is defined within the biological and covariate support represented
across the source and reference batches.

## 4. Reference contrasts

Let \(r\) denote the reference batch. For each non-reference batch \(b\), define

\[
\delta_{gb}=\gamma_{gb}-\gamma_{gr}.
\]

The unmoderated estimator and its QL-scaled variance are

\[
\widehat\delta_{gb}=c_b^T\widehat\gamma_g
\]

and

\[
s^2_{gb,\mathrm{QL}}
=
\widetilde\psi_g
c_b^T\operatorname{Cov}_0(\widehat\gamma_g)c_b,
\]

where \(c_b\) is the source-versus-reference contrast and
\(\widetilde\psi_g\) is the moderated QL scale used by the final edgeR fit.

## 5. QL-aware empirical-Bayes moderation

For each non-reference batch, use the heteroscedastic observation model

\[
\widehat\delta_{gb}\mid\delta_{gb}
\approx
N(\delta_{gb},s^2_{gb,\mathrm{QL}})
\]

and the single normal prior

\[
\delta_{gb}\sim N(m_b,\tau_b^2).
\]

Estimate \(m_b\) and \(\tau_b^2\ge0\) by marginal maximum likelihood or REML:

\[
\ell_b(m_b,\tau_b^2)
=
-\frac12\sum_{g\in\mathcal G_b}
\left[
\log(s^2_{gb,\mathrm{QL}}+\tau_b^2)
+
\frac{(\widehat\delta_{gb}-m_b)^2}
{s^2_{gb,\mathrm{QL}}+\tau_b^2}
\right],
\]

where \(\mathcal G_b\) contains genes with valid contrasts and finite positive
uncertainty.

The posterior mean is

\[
\widetilde\delta_{gb}
=
w_{gb}\widehat\delta_{gb}
+
(1-w_{gb})m_b,
\]

with

\[
w_{gb}
=
\frac{\tau_b^2}{\tau_b^2+s^2_{gb,\mathrm{QL}}}.
\]

Thus precise contrasts receive stronger data weight and uncertain contrasts
receive stronger prior weight.

The counterfactual reference mean is

\[
\mu^{(r)}_{gj}
=
\mu_{gj}\exp(-\widetilde\delta_{gb})
\]

for samples in source batch \(b\).

## 6. NB dispersion and QL dispersion

The NB dispersion defines a count distribution through

\[
\operatorname{Var}_{NB}(Y_{gj})
=
\mu_{gj}+\phi_{gb}\mu_{gj}^2.
\]

QL adds a separate residual scale through a mean-variance relationship such as

\[
\operatorname{Var}_{QL}(Y_{gj})
=
\psi_g
\left(\mu_{gj}+\phi_{gb}\mu_{gj}^2\right).
\]

QL specifies mean and variance behavior rather than a unique count CDF.
Therefore:

- \(\phi_{gb}\) is a valid NB distribution parameter;
- QL quantities govern uncertainty and moderation;
- QL quasi-dispersions do not enter `pnbinom()`, `dnbinom()`, or `qnbinom()`.

## 7. Hierarchical batch-specific dispersion

Version 0.4.0 models

\[
\log\phi_{gb}=f(A_g)+u_g+d_b,
\]

where:

- \(A_g\) is average abundance;
- \(f(A_g)\) is a shared smooth abundance-dispersion trend;
- \(u_g\) is a moderated gene-level deviation shared across batches;
- \(d_b\) is a batch-wide log-dispersion shift.

The reference constraint is

\[
d_r=0.
\]

Define

\[
\phi_{gr}=\exp\{f(A_g)+u_g\}
\]

and

\[
\kappa_b=\exp(d_b).
\]

Then

\[
\phi_{gb}=\phi_{gr}\kappa_b
\]

and

\[
\frac{\phi_{gb}}{\phi_{gr}}=\kappa_b.
\]

The model retains gene-specific absolute dispersion while representing the
source-versus-reference technical dispersion difference through a batch-wide
multiplier. Version 0.4.0 sets the gene-by-batch interaction to zero:

\[
h_{gb}=0.
\]

## 8. Effective residual information

Pooling genes supports estimation of shared trends and batch shifts while the
sample design determines residual information.

For an initial design-only diagnostic,

\[
H^{(0)}=X(X^TX)^{-1}X^T
\]

and

\[
\nu_b^{(0)}
=
\sum_{j:b(j)=b}(1-h^{(0)}_{jj}).
\]

For the fitted NB GLM, gene-specific working weights define

\[
H_g
=
W_g^{1/2}X
(X^TW_gX)^{-1}
X^TW_g^{1/2}
\]

and

\[
\nu_{gb}
=
\sum_{j:b(j)=b}(1-h_{gjj}).
\]

The dispersion precision contribution must satisfy:

\[
p_{gb}=0\quad\text{when}\quad\nu_{gb}=0
\]

and increase monotonically with \(\nu_{gb}\). A preferred form is

\[
p_{gb}
=
\left[
\psi_1(\nu_{gb}/2)
\right]^{-1},
\]

with \(\psi_1\) denoting the trigamma function. For moderate information,

\[
p_{gb}\approx\nu_{gb}/2.
\]

A saturated batch contributes zero information to its dispersion shift. A
weakly informed batch contributes proportionally less information. Zero
dispersion information produces a neutral source shift with a clear
estimability diagnostic while leaving an identifiable mean correction
available.

## 9. Distributional transport

For a supported count \(Y_{gj}=y\) from source batch \(b\), the source model is

\[
Y_{gj}\sim NB(\mu_{gj},\phi_{gb}).
\]

The reference-counterfactual target is

\[
Y^{(r)}_{gj}\sim NB(\mu^{(r)}_{gj},\phi_{gr}).
\]

Define the deterministic mid-P coordinate

\[
q_{gj}
=
F_b(y-1;\mu_{gj},\phi_{gb})
+
\frac12
P_b(Y=y;\mu_{gj},\phi_{gb}).
\]

The transported integer count is

\[
Y^*_{gj}
=
F_r^{-1}
(q_{gj};\mu^{(r)}_{gj},\phi_{gr}).
\]

The transformation aligns source counts to the reference mean-dispersion
relationship while preserving the reference observations.

## 10. Zero policy

Under the default model-consistent policy, an observed zero receives

\[
q_{gj}
=
\frac12P_b(Y=0;\mu_{gj},\phi_{gb}).
\]

The target quantile determines whether the transported value remains zero.

An explicit compatibility policy may retain observed zeros as zeros. That
policy adds a structural constraint to the transport and must be labelled in
the result diagnostics.

## 11. Reference invariance

For every reference-batch sample,

\[
Y^*_{gj}=Y_{gj}.
\]

This is an exact implementation invariant rather than an approximate fitted
property. The reference batch defines the target distribution and its observed
counts remain untouched.

## 12. Supported genes

Genes require sufficient count support for fitting and transport under the
package's documented support rule. Unsupported genes retain their original
counts and receive a structured status explaining the retained result.

Numerical mapping failure for a supported gene triggers gene-level rollback to
the original counts and a structured diagnostic.

## 13. Reference selection

An externally justified reference is preferred. When the reference is selected
automatically, the package reports:

- every candidate score;
- the selected candidate;
- the second-best candidate;
- the score margin;
- the biological terms retained in each local reference model;
- each candidate's effective residual information.

Automatic selection remains a heuristic and does not alter the exact
reference-invariance requirement.

## 14. Identifiability and interpretation

The mean design must separate biological variables from batch. Full
confounding leaves the corresponding effects unidentified.

Batch-dispersion estimation draws information from replicated residual
variation after accounting for locally estimable biological terms. Gene
pooling strengthens shared dispersion quantities while preserving the residual
information limits imposed by the sample design.

The batch-wide multiplier model assumes

\[
\phi_{gb}/\phi_{gr}=\kappa_b
\]

across genes. It represents a shared technical variance shift while the
abundance trend and moderated gene deviation represent gene-specific absolute
dispersion.

## 15. Default and ablation configurations

The 0.4.0 default is:

```r
mean_shrinkage = "ql_eb"
dispersion_model = "batch_scaled"
zero_policy = "model"
```

The following configurations support validation:

| Configuration | Mean correction | Dispersion correction |
|---|---|---|
| historical comparator | unmoderated | shared |
| A only | QL-EB | shared |
| B only | unmoderated | batch-specific hierarchical |
| full 0.4.0 | QL-EB | batch-specific hierarchical |

These ablations establish the contribution of each component without dividing
the method into separate release milestones.

## 16. Validation claims

Validation must separately establish:

1. recovery of simulated mean batch contrasts;
2. recovery of simulated dispersion multipliers;
3. preservation of biological contrasts;
4. minimal transformation under the no-batch-effect condition;
5. exact reference invariance;
6. stable behavior as residual information changes;
7. calibrated diagnostics for unsupported and weakly informed designs.

Claims in the README, vignette, or release notes must follow the completed
validation evidence.

## 17. Scope beyond 0.4.0

Mixture priors, gene-by-batch dispersion interactions, fully QL-native count
families, randomized PIT transport, semiparametric residual transport, and
distributional-medoid reference selection belong to later methodological work.
They are not part of the 0.4.0 model.
