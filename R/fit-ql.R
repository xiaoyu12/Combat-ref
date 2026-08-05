fit_ql_model <- function(dge, design, dispersion = NULL) {
  fit <- tryCatch(
    edgeR::glmQLFit(dge, design, dispersion = dispersion, robust = TRUE,
      legacy = FALSE, prior.count = 0, keep.unit.mat = TRUE),
    error = function(error) abort_combatrefql("combatrefql_ql_error",
      sprintf("edgeR v4 quasi-likelihood fitting failed: %s", conditionMessage(error)),
      "*" = "Check replication, design rank, and count scale.", stage = "ql_fit")
  )
  if (any(!is.finite(fit$dispersion)) || any(fit$dispersion < 0)) {
    abort_combatrefql("combatrefql_internal_error",
      "edgeR did not return valid negative-binomial dispersions.", stage = "ql_fit")
  }
  fit
}

contrast_covariance <- function(fit, dispersion, source, reference) {
  contrast <- batch_contrast(fit$design, source, reference)
  if (is.null(dim(dispersion))) dispersion <- matrix(dispersion,
    nrow(fit$fitted.values), ncol(fit$fitted.values))
  prior_weights <- fit$weights
  if (is.null(prior_weights)) prior_weights <- matrix(1,
    nrow(fit$fitted.values), ncol(fit$fitted.values))
  average_ql <- as.numeric(fit$average.ql.dispersion)
  working_dispersion <- dispersion / average_ql
  coefficients <- fit$unshrunk.coefficients %||% fit$coefficients
  raw <- variance <- rep(NA_real_, nrow(fit$fitted.values))
  for (g in seq_len(nrow(fit$fitted.values))) {
    mu <- fit$fitted.values[g, ]
    w <- prior_weights[g, ] * mu / (1 + working_dispersion[g, ] * mu)
    covariance <- try(solve(crossprod(fit$design * sqrt(w))), silent = TRUE)
    if (inherits(covariance, "try-error")) next
    raw[g] <- sum(contrast * coefficients[g, ])
    variance[g] <- fit$s2.post[g] *
      as.numeric(crossprod(contrast, covariance %*% contrast))
  }
  data.frame(gene = rownames(fit$fitted.values), source_batch = source,
    reference_batch = reference, raw_delta = raw,
    ql_variance = variance, ql_se = sqrt(variance),
    stringsAsFactors = FALSE)
}
