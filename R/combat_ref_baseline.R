fit_combat_ref_baseline <- function(counts, batch, group, covariates, reference_batch,
                                    lib_size, normalization, fit_method,
                                    quantile_method, keep_model) {
  timing <- list(); tick <- proc.time()[["elapsed"]]
  keep <- Reduce(`&`, lapply(levels(batch), function(b)
    rowSums(counts[, batch == b, drop = FALSE]) > 0))
  dge_all <- edgeR::DGEList(counts, lib.size = lib_size)
  raw_lib_size <- dge_all$samples$lib.size
  dge_all <- edgeR::normLibSizes(dge_all, method = normalization)
  effective_lib_size <- dge_all$samples$lib.size * dge_all$samples$norm.factors
  dge <- dge_all[keep, , keep.lib.sizes = TRUE]
  design <- build_combat_design(batch, group, covariates)
  timing$normalization <- proc.time()[["elapsed"]] - tick; tick <- proc.time()[["elapsed"]]

  local_designs <- lapply(levels(batch), function(b) {
    idx <- batch == b
    info <- build_within_batch_design(group, covariates, idx)
    score <- tryCatch({
      if (info$residual_df <= 0L) stop("no residual degrees of freedom")
      local <- dge[, idx, keep.lib.sizes = TRUE]
      estimate <- edgeR::estimateDisp(local, info$design, robust = TRUE)$common.dispersion
      if (length(estimate) != 1L || !is.finite(estimate) || estimate <= 0) {
        stop("non-finite or non-positive common dispersion")
      }
      estimate
    }, error = function(e) {
      info$error <<- conditionMessage(e)
      NA_real_
    })
    info$score <- score
    info$batch <- b
    info$n <- sum(idx)
    info
  })
  names(local_designs) <- levels(batch)
  batch_scores <- vapply(local_designs, `[[`, numeric(1), "score")
  if (all(!is.finite(batch_scores))) {
    cli::cli_abort("Automatic reference selection failed: all batch dispersion scores are invalid. Supply {.arg reference_batch} or add within-batch replication.")
  }
  if (is.null(reference_batch)) {
    valid_scores <- batch_scores[is.finite(batch_scores)]
    reference_batch <- sort(names(valid_scores)[valid_scores == min(valid_scores)])[1L]
  }
  warnings <- unique(unlist(lapply(local_designs, function(x) {
    c(if (isTRUE(x$fallback)) paste0("Batch ", x$batch, " local design fallback: ", x$reason),
      if (!is.null(x$error)) paste0("Batch ", x$batch, " reference score failed: ", x$error))
  })))
  timing$reference <- proc.time()[["elapsed"]] - tick; tick <- proc.time()[["elapsed"]]

  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  fit <- edgeR::glmFit(dge, design, dispersion = dge$tagwise.dispersion,
                       prior.count = 1e-4)
  observed_mu <- pmax(fit$fitted.values, 1e-8)
  reference_design <- counterfactual_reference_design(design, batch, reference_batch)
  offset <- fit$offset
  if (is.null(dim(offset))) offset <- matrix(offset, nrow(dge), ncol(dge), byrow = TRUE)
  reference_mu <- exp(offset + fit$coefficients %*% t(reference_design))
  finite_fit <- apply(is.finite(fit$coefficients), 1L, all) &
    apply(is.finite(observed_mu) & observed_mu > 0, 1L, all) &
    is.finite(dge$tagwise.dispersion) & dge$tagwise.dispersion > 0
  if (!is.null(fit$failed)) finite_fit <- finite_fit & !fit$failed
  timing$model <- proc.time()[["elapsed"]] - tick; tick <- proc.time()[["elapsed"]]

  adjusted <- counts[keep, , drop = FALSE]
  mapping_valid <- rep(TRUE, sum(keep))
  ref_dispersion <- batch_scores[[reference_batch]]
  if (!is.finite(ref_dispersion)) {
    cli::cli_abort("The requested reference batch has an invalid dispersion score.")
  }
  tail_events <- 0L
  for (b in setdiff(levels(batch), reference_batch)) {
    if (!is.finite(batch_scores[[b]])) {
      cli::cli_abort("Batch {.val {b}} has no valid dispersion score for count mapping.")
    }
    idx <- batch == b
    mapped <- match_quantiles(
      adjusted[, idx, drop = FALSE], observed_mu[, idx, drop = FALSE],
      batch_scores[[b]], reference_mu[, idx, drop = FALSE], ref_dispersion,
      method = quantile_method
    )
    tail_events <- tail_events + attr(mapped, "tail_events")
    mapping_valid <- mapping_valid & attr(mapped, "mapping_valid")
    adjusted[, idx] <- mapped
  }
  whole <- counts; whole[keep, ] <- adjusted
  timing$mapping <- proc.time()[["elapsed"]] - tick

  ql <- NULL; ql_nb_dispersion <- NULL
  if (fit_method == "ql_assisted") {
    if (utils::packageVersion("edgeR") < "4.2.0") {
      cli::cli_abort("QL-assisted fitting requires edgeR >= 4.2.0.")
    }
    ql <- edgeR::glmQLFit(dge, design, dispersion = NULL, robust = TRUE,
                          legacy = FALSE)
    ql_nb_dispersion <- ql$dispersion
  }
  model <- if (keep_model) list(
    design = design, reference_design = reference_design,
    coefficients = fit$coefficients, fitted_means = observed_mu,
    reference_means = reference_mu,
    baseline_nb_tagwise_dispersion = dge$tagwise.dispersion,
    ql_nb_global_dispersion = ql_nb_dispersion,
    ql_dispersion = if (is.null(ql)) NULL else ql$s2.post
  ) else NULL
  list(
    counts = whole, keep = keep, reference_batch = reference_batch,
    batch_scores = batch_scores, local_designs = local_designs,
    raw_lib_size = raw_lib_size,
    normalization_factors = dge_all$samples$norm.factors,
    effective_lib_size = effective_lib_size, design = design,
    residual_df = nrow(design) - ncol(design), tail_events = tail_events,
    fit_valid = finite_fit, mapping_valid = mapping_valid,
    ql_nb_global_dispersion = ql_nb_dispersion,
    ql_summary = if (is.null(ql)) NULL else summary(ql$s2.post),
    warnings = warnings, timing = timing, model = model
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
