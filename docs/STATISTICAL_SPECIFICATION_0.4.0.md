# ComBat-refQL 0.4.0 statistical specification

## Estimand and identifiability

ComBat-refQL estimates how each supported gene differs between every source
batch and a selected reference after accounting for biological group, supplied
covariates, and effective library size. It transports source observations to
the fitted reference-batch distribution. Batch effects that are perfectly
confounded with biology are not identifiable; the package requires a full-rank
design with positive residual degrees of freedom.

## Mean model and contrasts

For gene $g$ and sample $j$,

$$
\log\mu_{gj}=\log L_j+X_j\beta_g+\gamma_{g,b(j)},
$$

where $L_j$ is the TMMwsp effective library size. The design contains one
coefficient per batch plus estimable biological and covariate terms. For
reference $r$, the exact source contrast is

$$\delta_{gb}=\gamma_{gb}-\gamma_{gr}.$$

The same contrast vector is used for coefficient extraction and covariance.
Changing factor order or an equivalent parameterization does not change the
estimand.

## QL-aware empirical Bayes

The final edgeR v4 fit supplies a raw contrast and its moderated QL variance:

$$\widehat\delta_{gb}\mid\delta_{gb}\approx
N(\delta_{gb},s^2_{gb,\mathrm{QL}}).$$

For each non-reference batch, 0.4.0 uses one heteroscedastic normal prior,

$$\delta_{gb}\sim N(m_b,\tau_b^2).$$

The implementation profiles the marginal likelihood over $\tau_b^2\ge0$ and
computes the precision-weighted $m_b$. The posterior mean is

$$
\widetilde\delta_{gb}=w_{gb}\widehat\delta_{gb}+(1-w_{gb})m_b,
\qquad
w_{gb}=\frac{\tau_b^2}{\tau_b^2+s^2_{gb,\mathrm{QL}}}.
$$

The boundary $\tau_b^2=0$ gives complete pooling. Invalid standard errors and
fewer than two supported genes receive explicit statuses. With
`mean_shrinkage = "none"`, the raw contrast is used.

## NB dispersion and QL dispersion

NB dispersion determines the count mean-variance relationship and is the only
dispersion used in `dnbinom()`, `pnbinom()`, and `qnbinom()`. edgeR posterior QL
dispersion scales coefficient uncertainty. It never enters an NB CDF.

Version 0.4.0 models

$$
\log\phi_{gb}=f(A_g)+u_g+d_b,\qquad d_r=0,
$$

or equivalently $\phi_{gb}=\phi_{gr}\exp(d_b)$. There is no gene-by-batch
dispersion interaction. edgeR `estimateDisp()` provides the moderated
gene-level baseline. With preliminary fitted means held fixed, each batch-wide
log multiplier is estimated by bounded, effective-information-weighted NB
maximum likelihood. The raw multipliers are then centred on the reference.

The final mean model is refitted by `glmQLFit(..., legacy = FALSE)` with the
gene-by-sample NB dispersion matrix. `dispersion_model = "shared"` instead uses
the preliminary shared dispersion and fixes every $d_b$ to zero.

## Effective residual information

The design-only diagnostic is

$$H^{(0)}=X(X^TX)^{-1}X^T,\qquad
\nu_b^{(0)}=\sum_{j:b(j)=b}(1-h^{(0)}_{jj}).$$

After preliminary fitting, gene-level working weights define

$$
H_g=W_g^{1/2}X(X^TW_gX)^{-1}X^TW_g^{1/2},\qquad
\nu_{gb}=\sum_{j:b(j)=b}(1-h_{gjj}).
$$

Dispersion likelihood contributions use

$$p_{gb}=\{\psi_1(\nu_{gb}/2)\}^{-1}.$$

This precision increases continuously with information and is zero at
$\nu_{gb}=0$. A source batch with zero total precision receives $d_b=0$ and
`not_estimable`; its identifiable mean correction remains active. A reference
with zero total precision raises `combatrefql_statistical_error`, with
`dispersion_model = "shared"` named as the explicit fallback.

## Counterfactual and transport

For a source sample,

$$\mu^{(r)}_{gj}=\mu_{gj}\exp(-\widetilde\delta_{gb}).$$

Given observed $y$, source $(\mu_{gj},\phi_{gb})$, and target
$(\mu^{(r)}_{gj},\phi_{gr})$, deterministic mid-P transport calculates

$$
q=F_b(y-1)+\tfrac12P_b(Y=y),\qquad y^*=F_r^{-1}(q).
$$

Under the default `zero_policy = "model"`, $q=\tfrac12P_b(Y=0)$ for zero.
`zero_policy = "preserve"` sets adjusted zeros to zero for compatibility.
Probabilities are clipped at machine-safe tails. Invalid gene mappings roll
back the complete gene row. Chunk size changes allocation only, not results.

Reference samples bypass transport and are checked for exact invariance.

## Supported conditions and fallbacks

Genes must have non-zero support in every batch. Unsupported genes return
unchanged with a recorded reason. The design must be full rank with positive
global residual degrees of freedom. Explicit references are preferred when
study design or quality control identifies a target. Automatic selection fits
audited within-batch biological designs, selects the smallest valid NB score,
and records the best score, second-best candidate, and margin.

The formal configurations are `none + shared` (0.3.1-style), `ql_eb + shared`
(A-only), `none + batch_scaled` (B-only), and `ql_eb + batch_scaled` (0.4.0
default). Randomized transport, mixture priors, gene-by-batch dispersion
interactions, and automatic tuning are outside the 0.4.0 model.
