#' Correct RNA-seq batch effects toward a reference batch
#'
#' Fits a negative-binomial model with one named coefficient per batch, then
#' predicts each sample under its observed and counterfactual reference-batch
#' designs. Count mapping always uses NB dispersion. The `"ql_assisted"` method
#' adds edgeR QL diagnostics but does not treat QL dispersion as a distribution.
#'
#' @param counts Raw non-negative integer matrix, genes by samples.
#' @param batch Batch vector aligned with columns of `counts`.
#' @param group Optional biological-group vector.
#' @param covariates Optional sample-aligned data frame or matrix.
#' @param reference_batch Optional explicit reference batch.
#' @param normalization edgeR normalization method; `"TMMwsp"` is the default.
#' @param fit_method Repaired NB baseline or experimental QL-assisted diagnostics.
#' @param quantile_method Deterministic mid-P, upper-CDF, or randomized mapping.
#' @param lib_size Optional positive library sizes, applied before normalization.
#' @param seed Seed required for randomized mapping.
#' @param verbose Print concise progress messages.
#' @param keep_model Retain fitted coefficients and mean matrices.
#' @return A `combatref_result` object.
#' @export
combat_ref <- function(counts, batch, group = NULL, covariates = NULL,
                       reference_batch = NULL,
                       normalization = c("TMMwsp", "TMM", "none"),
                       fit_method = c("baseline", "ql_assisted"),
                       quantile_method = c("midp", "upper", "randomized"),
                       lib_size = NULL, seed = NULL, verbose = interactive(),
                       keep_model = FALSE) {
  call <- match.call()
  normalization <- match.arg(normalization)
  fit_method <- match.arg(fit_method)
  quantile_method <- match.arg(quantile_method)
  if (quantile_method == "randomized" && is.null(seed)) {
    cli::cli_abort("{.arg seed} is required for randomized quantile mapping.")
  }
  if (!is.null(seed)) set.seed(seed)
  validated <- validate_combat_inputs(counts, batch, group, covariates,
                                      reference_batch, lib_size)
  if (isTRUE(verbose)) cli::cli_alert_info("Fitting {fit_method} ComBat-ref model")
  fit <- suppressMessages(fit_combat_ref_baseline(
    validated$counts, validated$batch, validated$group, validated$covariates,
    validated$reference_batch, validated$lib_size, normalization, fit_method,
    quantile_method, keep_model
  ))
  unsupported <- !fit$keep
  gene_status <- data.frame(
    gene = rownames(validated$counts),
    status = ifelse(unsupported, "unsupported_zero_batch", "adjusted"),
    reason = ifelse(unsupported, "all zero in at least one batch", NA_character_),
    adjusted = !unsupported,
    zero_batches = vapply(seq_len(nrow(validated$counts)), function(i)
      paste(levels(validated$batch)[vapply(levels(validated$batch), function(b)
        all(validated$counts[i, validated$batch == b] == 0), logical(1))], collapse = ","), character(1)),
    low_information = rowSums(validated$counts) == 0,
    fit_converged = !unsupported,
    mapping_converged = !unsupported,
    stringsAsFactors = FALSE
  )
  diagnostics <- list(
    n_genes_input = nrow(validated$counts), n_genes_adjusted = sum(fit$keep),
    n_genes_filtered = sum(unsupported), n_samples = ncol(validated$counts),
    batch_sizes = table(validated$batch),
    group_sizes = if (is.null(validated$group)) NULL else table(validated$group),
    design_rank = qr(fit$design)$rank, design_columns = ncol(fit$design),
    condition_number = kappa(fit$design), residual_df = fit$residual_df,
    raw_library_sizes = fit$raw_lib_size,
    normalization_factors = fit$normalization_factors,
    effective_library_sizes = fit$effective_lib_size,
    reference_scores = fit$batch_scores,
    reference_selection_frequency = stats::setNames(as.numeric(names(fit$batch_scores) == fit$reference_batch), names(fit$batch_scores)),
    reference_stable = NA, numerical_tail_events = fit$tail_events,
    ql_dispersion_summary = fit$ql_summary,
    package_version = as.character(utils::packageVersion("combatrefql")),
    edgeR_version = as.character(utils::packageVersion("edgeR"))
  )
  new_combatref_result(
    fit$counts, fit$reference_batch, normalization, fit_method, quantile_method,
    validated$batch, validated$group, validated$covariates, gene_status,
    rownames(validated$counts)[unsupported], diagnostics, fit$model, call
  )
}

#' Experimental QL-assisted ComBat-ref diagnostics
#' @param ... Arguments passed to [combat_ref()].
#' @export
combat_ref_ql <- function(...) combat_ref(..., fit_method = "ql_assisted")

format_named_counts <- function(x) paste0(names(x), " (n = ", as.integer(x), ")", collapse = ", ")
