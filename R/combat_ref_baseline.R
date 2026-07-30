fit_combat_ref_baseline <- function(counts, batch, group, covariates, reference_batch,
                                    lib_size, normalization, fit_method,
                                    quantile_method, keep_model) {
  keep <- Reduce(`&`, lapply(levels(batch), function(b) rowSums(counts[, batch == b, drop = FALSE]) > 0))
  dge_all <- edgeR::DGEList(counts, lib.size = lib_size)
  raw_lib_size <- dge_all$samples$lib.size
  dge_all <- edgeR::normLibSizes(dge_all, method = normalization)
  effective_lib_size <- dge_all$samples$lib.size * dge_all$samples$norm.factors
  dge <- dge_all[keep, , keep.lib.sizes = TRUE]
  design <- build_combat_design(batch, group, covariates)

  batch_scores <- vapply(levels(batch), function(b) {
    idx <- batch == b
    local_design <- design[idx, !startsWith(colnames(design), "batch"), drop = FALSE]
    if (!ncol(local_design)) local_design <- matrix(1, sum(idx), 1L)
    if (qr(local_design)$rank < ncol(local_design) || nrow(local_design) <= ncol(local_design)) {
      local_design <- matrix(1, sum(idx), 1L)
    }
    local <- dge[, idx, keep.lib.sizes = TRUE]
    local <- edgeR::estimateDisp(local, local_design, robust = TRUE)
    local$common.dispersion
  }, numeric(1))
  reference_batch <- reference_batch %||% names(which.min(batch_scores))

  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  fit <- edgeR::glmFit(dge, design, dispersion = dge$tagwise.dispersion, prior.count = 1e-4)
  observed_mu <- pmax(fit$fitted.values, 1e-8)
  reference_design <- counterfactual_reference_design(design, batch, reference_batch)
  offset <- fit$offset
  if (is.null(dim(offset))) offset <- matrix(offset, nrow(dge), ncol(dge), byrow = TRUE)
  reference_mu <- exp(offset + fit$coefficients %*% t(reference_design))

  adjusted <- counts[keep, , drop = FALSE]
  ref_dispersion <- batch_scores[[reference_batch]]
  tail_events <- 0L
  for (b in setdiff(levels(batch), reference_batch)) {
    idx <- batch == b
    mapped <- match_quantiles(
      adjusted[, idx, drop = FALSE], observed_mu[, idx, drop = FALSE],
      batch_scores[[b]], reference_mu[, idx, drop = FALSE], ref_dispersion,
      method = quantile_method
    )
    tail_events <- tail_events + attr(mapped, "tail_events")
    adjusted[, idx] <- mapped
  }
  whole <- counts
  whole[keep, ] <- adjusted
  ql <- NULL
  if (fit_method == "ql_assisted") {
    ql <- edgeR::glmQLFit(dge, design, dispersion = dge$tagwise.dispersion, robust = TRUE)
  }
  model <- if (keep_model) list(
    design = design, reference_design = reference_design,
    coefficients = fit$coefficients, fitted_means = observed_mu,
    reference_means = reference_mu, nb_dispersion = dge$tagwise.dispersion,
    ql_dispersion = if (is.null(ql)) NULL else ql$s2.post
  ) else NULL
  list(
    counts = whole, keep = keep, reference_batch = reference_batch,
    batch_scores = batch_scores, raw_lib_size = raw_lib_size,
    normalization_factors = dge_all$samples$norm.factors,
    effective_lib_size = effective_lib_size, design = design,
    residual_df = nrow(design) - ncol(design), tail_events = tail_events,
    ql_summary = if (is.null(ql)) NULL else summary(ql$s2.post), model = model
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
