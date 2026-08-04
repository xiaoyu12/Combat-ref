#' Fit ComBat-refQL reference-batch correction
#'
#' Fits robust edgeR v4 quasi-likelihood counterfactual means and transports
#' non-reference counts with deterministic negative-binomial mid-P quantiles.
#' Quasi-dispersion is never used as an NB distribution parameter.
#'
#' @param counts Numeric raw-count matrix with genes in rows and samples in columns.
#' @param batch Sample batch vector, optionally named for automatic reordering.
#' @param group Optional biological group vector.
#' @param covariates Optional sample-by-covariate data frame or matrix.
#' @param reference Optional reference batch. The minimum-dispersion batch is selected when `NULL`.
#' @param normalization edgeR library normalization method. Only `"TMMwsp"` is supported.
#' @param fractional_counts Round with an audit record, or reject fractional values.
#' @param chunk_size Positive number of genes transported per chunk.
#' @param verbosity User output level.
#' @return A validated [CombatRefQLFit] object.
#' @export
combat_ref_ql <- function(counts, batch, group = NULL, covariates = NULL,
                          reference = NULL, normalization = "TMMwsp",
                          fractional_counts = c("round", "error"), chunk_size = 5000L,
                          verbosity = c("normal", "quiet", "verbose")) {
  call <- match.call(); fractional_counts <- match.arg(fractional_counts); verbosity <- match.arg(verbosity)
  if (!identical(normalization, "TMMwsp")) input_error("{.arg normalization} must be {.val TMMwsp} in version 0.2.1.")
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1 || chunk_size != as.integer(chunk_size)) input_error("{.arg chunk_size} must be one positive integer.")
  chunk_size <- as.integer(chunk_size)
  emit <- function(type, text) { if (verbosity == "quiet") return(invisible()); if (type == "success") cli::cli_alert_success(text) else cli::cli_text(text) }
  total_start <- proc.time()[["elapsed"]]; emit("text", "ComBat-refQL 0.2.1")
  start <- proc.time()[["elapsed"]]
  prepared <- prepare_inputs(counts, batch, group, covariates, reference, fractional_counts)
  timings <- timing_row("preparation", start)
  if (prepared$input_actions$affected_values[[1L]]) emit("text", prepared$input_actions$message[[1L]])
  emit("success", sprintf("Prepared %s genes x %s samples", nrow(prepared$counts), ncol(prepared$counts)))
  start <- proc.time()[["elapsed"]]; dge <- edgeR::normLibSizes(edgeR::DGEList(prepared$counts), method = normalization)
  timings <- rbind(timings, timing_row("normalization", start)); emit("success", "Normalized with TMMwsp")
  design <- build_design(prepared$batch, prepared$group, prepared$covariates)
  keep <- Reduce(`&`, lapply(levels(prepared$batch), function(x) rowSums(prepared$counts[, prepared$batch == x, drop = FALSE]) > 0))
  if (sum(keep) < 2L) abort_combatrefql("combatrefql_ql_error", "Fewer than two genes have count support in every batch.", stage = "ql_fit", affected_genes = sum(keep))
  reference_progress <- function(candidate) if (verbosity == "verbose" && !interactive()) emit("text", sprintf("Scored reference candidate: %s", candidate))
  start <- proc.time()[["elapsed"]]; selected <- select_reference(dge, prepared$batch, prepared$reference, reference_progress)
  timings <- rbind(timings, timing_row("reference", start)); emit("success", sprintf("Selected reference: %s", selected$reference))
  start <- proc.time()[["elapsed"]]; fit <- fit_ql_model(dge[keep, , keep.lib.sizes = TRUE], design$matrix)
  timings <- rbind(timings, timing_row("ql_fit", start)); emit("success", "Fitted robust edgeR v4 QL model")
  start <- proc.time()[["elapsed"]]
  mapping_steps <- length(setdiff(levels(prepared$batch), selected$reference)) * ceiling(sum(keep) / chunk_size)
  progress_id <- if (verbosity != "quiet" && interactive()) cli::cli_progress_bar("Mapping batches", total = mapping_steps) else NULL
  completed_steps <- 0L
  mapping_progress <- function(source, genes) {
    completed_steps <<- completed_steps + 1L
    if (!is.null(progress_id)) cli::cli_progress_update(id = progress_id, set = completed_steps)
    else if (verbosity == "verbose" && (completed_steps == 1L || completed_steps == mapping_steps || completed_steps %% max(1L, ceiling(mapping_steps / 5)) == 0L)) emit("text", sprintf("Mapped batch %s: %s/%s gene chunks", source, completed_steps, mapping_steps))
  }
  transported <- transport_counts(prepared$counts, keep, prepared$batch, selected$reference, design$matrix, fit, chunk_size, mapping_progress)
  if (!is.null(progress_id)) cli::cli_progress_done(progress_id)
  timings <- rbind(timings, timing_row("mapping", start))
  gene_status <- make_gene_status(prepared$counts, prepared$batch, keep,
    transported$mapping_valid)
  outcomes <- make_outcomes(gene_status)
  mapping_warnings <- if (nrow(transported$failures)) data.frame(
    stage = "mapping", class = "combatrefql_mapping_error",
    message = sprintf("NB transport failed for %s genes; complete original rows were retained.",
      sum(gene_status$status == "failed")), stringsAsFactors = FALSE) else empty_warnings()
  if (nrow(mapping_warnings)) emit("text", mapping_warnings$message[[1L]])
  sample_data <- data.frame(sample = colnames(prepared$counts), batch = as.character(prepared$batch), group = if (is.null(prepared$group)) NA_character_ else as.character(prepared$group), raw_library_size = colSums(prepared$counts), normalization_factor = dge$samples$norm.factors, effective_library_size = dge$samples$lib.size * dge$samples$norm.factors, stringsAsFactors = FALSE)
  specification <- data.frame(method = "QL-fitted counterfactual means with negative-binomial mid-P quantile transport", normalization = normalization, transport = "negative-binomial mid-P", preserve_zeros = TRUE, reference_batch = selected$reference, reference_selection = if (is.null(prepared$reference)) "automatic" else "explicit", chunk_size = chunk_size, stringsAsFactors = FALSE)
  ql <- data.frame(edgeR_pipeline = "v4", robust = TRUE, legacy = FALSE, nb_dispersion = as.numeric(fit$dispersion), average_ql_dispersion = as.numeric(fit$average.ql.dispersion), median_posterior_ql_dispersion = stats::median(fit$s2.post), stringsAsFactors = FALSE)
  model <- list(coefficients = fit$coefficients, nb_dispersion = as.numeric(fit$dispersion), ql_dispersion = fit$s2.prior, posterior_ql_dispersion = fit$s2.post, design = data.frame(samples = nrow(design$matrix), coefficients = ncol(design$matrix), rank = design$rank, residual_df = design$residual_df, condition_number = design$condition_number), design_matrix = design$matrix)
  timings <- rbind(timings, data.frame(stage = "total", seconds = proc.time()[["elapsed"]] - total_start))
  diagnostics <- list(batches = selected$batches, mapping = transported$mapping,
    mapping_failures = transported$failures, timing = timings,
    input_actions = prepared$input_actions, warnings = mapping_warnings,
    outcomes = outcomes, ql = ql)
  result <- CombatRefQLFit(counts = transported$counts, specification = specification, samples = sample_data, gene_status = gene_status, diagnostics = diagnostics, call = call, model = model)
  emit("success", "Constructed and validated result"); emit("text", sprintf("Completed in %.1f seconds", timings$seconds[timings$stage == "total"])); result
}
