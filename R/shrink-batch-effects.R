fit_normal_eb <- function(delta, variance) {
  valid <- is.finite(delta) & is.finite(variance) & variance > 0
  if (sum(valid) < 2L) return(list(mean = NA_real_, variance = NA_real_,
    weight = rep(1, length(delta)), posterior = delta,
    valid = valid, status = "insufficient_genes", convergence = FALSE))
  y <- delta[valid]
  s2 <- variance[valid]
  profile <- function(tau2) {
    precision <- 1 / (s2 + tau2)
    mean <- sum(precision * y) / sum(precision)
    0.5 * sum(log(s2 + tau2) + (y - mean)^2 / (s2 + tau2))
  }
  upper <- max(stats::var(y) * 100, max(s2), 1e-8)
  opt <- stats::optimize(profile, c(0, upper))
  tau2 <- if (profile(0) <= opt$objective + sqrt(.Machine$double.eps)) 0 else
    opt$minimum
  precision <- 1 / (s2 + tau2)
  prior_mean <- sum(precision * y) / sum(precision)
  weight <- rep(NA_real_, length(delta))
  weight[valid] <- tau2 / (tau2 + s2)
  posterior <- delta
  posterior[valid] <- weight[valid] * y + (1 - weight[valid]) * prior_mean
  list(mean = prior_mean, variance = tau2, weight = weight,
    posterior = posterior, valid = valid,
    status = if (tau2 == 0) "complete_pooling" else "estimated",
    convergence = TRUE)
}

shrink_batch_effects <- function(fit, dispersion, batch, reference) {
  rows <- lapply(setdiff(levels(batch), reference), function(source) {
    contrast <- contrast_covariance(fit, dispersion, source, reference)
    eb <- fit_normal_eb(contrast$raw_delta, contrast$ql_variance)
    transform(contrast, prior_mean = eb$mean, prior_variance = eb$variance,
      shrinkage_weight = eb$weight, posterior_delta = eb$posterior,
      valid = eb$valid, status = eb$status, converged = eb$convergence)
  })
  diagnostics <- do.call(rbind, rows)
  sources <- unique(diagnostics$source_batch)
  posterior <- lapply(sources, function(source) stats::setNames(
    diagnostics$posterior_delta[diagnostics$source_batch == source],
    diagnostics$gene[diagnostics$source_batch == source]))
  names(posterior) <- sources
  list(posterior = posterior, diagnostics = diagnostics)
}
