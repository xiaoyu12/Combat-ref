#' Correct RNA-seq batch effects toward a reference batch
#'
#' @param counts Raw count matrix or safely coercible data frame, genes by samples.
#' @param batch Batch vector aligned with samples.
#' @param group Optional biological-group vector.
#' @param covariates Optional sample-aligned data frame or matrix.
#' @param reference_batch Optional explicit reference batch.
#' @param normalization edgeR normalization: `"TMMwsp"`, `"TMM"`, or `"none"`.
#' @param fit_method `"baseline"` or experimental `"ql_assisted"`.
#' @param quantile_method `"midp"`, `"upper"`, or `"randomized"`.
#' @param fractional_counts `"error"` or explicit lossy `"round"`.
#' @param duplicate_genes `"error"` or explicit aggregation by `"sum"`.
#' @param input_scale Must be `"counts"`; declared normalized scales are rejected.
#' @param lib_size Optional positive library sizes, applied before normalization.
#' @param seed Seed required for randomized mapping.
#' @param verbosity Messaging level: auto, quiet, normal, or verbose.
#' @param verbose Deprecated compatibility switch; use `verbosity`.
#' @param keep_model Retain fitted coefficients and mean matrices.
#' @return A `combatref_result`; input repairs are in `diagnostics$input_actions`.
#' @export
combat_ref <- function(counts, batch, group = NULL, covariates = NULL,
                       reference_batch = NULL, normalization = "TMMwsp",
                       fit_method = "baseline", quantile_method = "midp",
                       fractional_counts = c("error", "round"),
                       duplicate_genes = c("error", "sum"),
                       input_scale = "counts",
                       lib_size = NULL, seed = NULL,
                       verbosity = c("auto", "quiet", "normal", "verbose"),
                       verbose = NULL, keep_model = FALSE) {
  call <- match.call(); started <- proc.time()[["elapsed"]]
  verbosity <- match.arg(verbosity)
  if (!is.null(verbose)) verbosity <- if (isTRUE(verbose)) "verbose" else "quiet"
  if (verbosity == "auto") verbosity <- if (interactive()) "normal" else "quiet"
  fractional_counts <- match.arg(fractional_counts)
  duplicate_genes <- match.arg(duplicate_genes)
  normal <- verbosity %in% c("normal", "verbose")

  method_actions <- empty_input_actions()
  nm <- canonical_method(normalization, c("TMMwsp", "TMM", "none"), "normalization")
  fm <- canonical_method(fit_method, c("baseline", "ql_assisted"), "fit_method")
  qm <- canonical_method(quantile_method, c("midp", "upper", "randomized"), "quantile_method")
  normalization <- nm$value; fit_method <- fm$value; quantile_method <- qm$value
  method_actions <- rbind(nm$actions, fm$actions, qm$actions)
  if (quantile_method == "randomized" && is.null(seed)) cli::cli_abort("{.arg seed} is required for randomized quantile mapping.")
  if (!is.null(seed)) set.seed(seed)

  prepared <- prepare_combat_inputs(counts, batch, group, covariates,
    reference_batch, fractional_counts, duplicate_genes, input_scale)
  prepared$actions <- rbind(prepared$actions, method_actions)
  if (normal && nrow(prepared$actions)) {
    apply(prepared$actions, 1L, function(x) cli::cli_alert_info(x[["reason"]]))
  }
  validated <- validate_model_input(prepared$counts, prepared$batch,
    prepared$group, prepared$covariates, prepared$reference_batch, lib_size)
  if (normal) cli::cli_alert_info("Fitting {fit_method} model for {nrow(validated$counts)} genes x {ncol(validated$counts)} samples")
  captured_warnings <- character()
  fit <- withCallingHandlers(
    fit_combat_ref_baseline(validated$counts, validated$batch, validated$group,
      validated$covariates, validated$reference_batch, validated$lib_size,
      normalization, fit_method, quantile_method, keep_model),
    warning = function(w) { captured_warnings <<- c(captured_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  unsupported <- !fit$keep
  fit_state <- mapping_state <- rep(NA, nrow(validated$counts))
  fit_state[fit$keep] <- fit$fit_valid
  mapping_state[fit$keep] <- fit$mapping_valid
  supported_index <- which(fit$keep)
  reason <- rep(NA_character_, nrow(validated$counts))
  reason[unsupported] <- "all zero in at least one batch"
  reason[supported_index[!fit$fit_valid]] <- "non-finite coefficient, mean, or NB dispersion"
  reason[supported_index[!fit$mapping_valid]] <- "non-finite, negative, or non-integer mapped count"
  gene_status <- data.frame(
    gene = rownames(validated$counts),
    status = ifelse(unsupported, "unsupported_zero_batch",
      ifelse(fit_state == FALSE | mapping_state == FALSE, "failed", "adjusted")),
    reason = reason, adjusted = !unsupported & fit_state != FALSE & mapping_state != FALSE,
    zero_batches = vapply(seq_len(nrow(validated$counts)), function(i)
      paste(levels(validated$batch)[vapply(levels(validated$batch), function(b)
        all(validated$counts[i, validated$batch == b] == 0), logical(1))], collapse = ","), character(1)),
    low_information = rowSums(validated$counts) == 0,
    fit_converged = fit_state, mapping_converged = mapping_state,
    stringsAsFactors = FALSE
  )
  local_table <- do.call(rbind, lapply(fit$local_designs, function(x) data.frame(
    batch = x$batch, n = x$n, score = x$score, formula = x$formula,
    columns = paste(colnames(x$design), collapse = ","),
    columns_dropped = paste(x$columns_dropped, collapse = ","),
    fallback = x$fallback, reason = x$reason, rank = x$rank,
    residual_df = x$residual_df, stringsAsFactors = FALSE)))
  warnings <- unique(c(prepared$warnings, fit$warnings, captured_warnings))
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
    reference_scores = fit$batch_scores, reference_local_designs = local_table,
    reference_selected = fit$reference_batch,
    reference_selection_method = if (is.null(reference_batch)) "minimum valid robust within-batch common dispersion" else "explicit",
    reference_selection_frequency = NULL, reference_stable = NA,
    numerical_tail_events = fit$tail_events,
    baseline_nb_tagwise_dispersion_summary = summary(if (is.null(fit$model)) numeric() else fit$model$baseline_nb_tagwise_dispersion),
    ql_nb_global_dispersion = fit$ql_nb_global_dispersion,
    ql_dispersion_summary = fit$ql_summary,
    input_actions = prepared$actions, warnings = warnings,
    timing_seconds = c(fit$timing, total = proc.time()[["elapsed"]] - started),
    package_version = as.character(utils::packageVersion("combatrefql")),
    edgeR_version = as.character(utils::packageVersion("edgeR"))
  )
  if (normal) {
    cli::cli_alert_success("Selected reference batch: {fit$reference_batch}")
    cli::cli_alert_success("Adjusted {sum(fit$keep)} supported genes; {sum(unsupported)} unchanged")
    if (length(warnings)) for (w in warnings) cli::cli_alert_warning(w)
  }
  new_combatref_result(fit$counts, fit$reference_batch, normalization, fit_method,
    quantile_method, validated$batch, validated$group, validated$covariates,
    gene_status, rownames(validated$counts)[unsupported], diagnostics,
    fit$model, call)
}

#' Experimental QL-assisted ComBat-ref diagnostics
#' @param ... Arguments passed to [combat_ref()].
#' @export
combat_ref_ql <- function(...) combat_ref(..., fit_method = "ql_assisted")

format_named_counts <- function(x) paste0(names(x), " (n = ", as.integer(x), ")", collapse = ", ")
