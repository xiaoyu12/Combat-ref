#' Fit ComBat-refQL reference-batch adjustment
#'
#' Fits model-based reference-counterfactual means with edgeR v4 and adjusts
#' non-reference counts by deterministic negative-binomial mid-P transport.
#'
#' @param counts Numeric raw-count matrix with genes in rows and samples in columns.
#' @param batch Sample batch vector, optionally named for automatic reordering.
#' @param group Optional biological group vector.
#' @param covariates Optional sample-by-covariate data frame or matrix.
#' @param reference Optional reference batch. When `NULL`, the minimum-dispersion
#'   batch is selected from biologically adjusted within-batch QL fits.
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
  call <- match.call()
  fractional_counts <- match.arg(fractional_counts)
  verbosity <- match.arg(verbosity)
  if (!identical(normalization, "TMMwsp"))
    input_error("{.arg normalization} must be {.val TMMwsp}.")
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1 ||
      chunk_size != as.integer(chunk_size))
    input_error("{.arg chunk_size} must be one positive integer.")
  chunk_size <- as.integer(chunk_size)
  reporter <- new_reporter(verbosity)
  total_start <- proc.time()[["elapsed"]]
  reporter$start()

  start <- proc.time()[["elapsed"]]
  prepared <- prepare_inputs(counts, batch, group, covariates, reference, fractional_counts)
  timings <- timing_row("preparation", start)
  if (prepared$input_actions$affected_values[[1L]])
    reporter$info(prepared$input_actions$message[[1L]])
  reporter$success(sprintf("Prepared %s genes x %s samples",
    nrow(prepared$counts), ncol(prepared$counts)))

  start <- proc.time()[["elapsed"]]
  dge <- normalize_counts(prepared, normalization)
  timings <- rbind(timings, timing_row("normalization", start))
  reporter$success("Normalized with TMMwsp")

  design <- build_design(prepared$batch, prepared$group, prepared$covariates)
  keep <- assess_support(prepared$counts, prepared$batch)
  start <- proc.time()[["elapsed"]]
  selected <- select_reference(dge, prepared$batch, prepared$group,
    prepared$covariates, prepared$reference, reporter$reference)
  timings <- rbind(timings, timing_row("reference", start))
  reporter$success(sprintf("Selected reference: %s", selected$reference))

  start <- proc.time()[["elapsed"]]
  fit <- fit_ql_model(dge[keep, , keep.lib.sizes = TRUE], design$matrix)
  timings <- rbind(timings, timing_row("ql_fit", start))
  reporter$success("Fitted robust edgeR v4 QL model")

  start <- proc.time()[["elapsed"]]
  mapping_steps <- length(setdiff(levels(prepared$batch), selected$reference)) * ceiling(sum(keep) / chunk_size)
  reporter$start_mapping(mapping_steps)
  transported <- transport_counts(prepared$counts, keep, prepared$batch,
    selected$reference, design$matrix, fit, chunk_size, reporter$mapping)
  timings <- rbind(timings, timing_row("mapping", start))
  timings <- rbind(timings, data.frame(stage = "total", seconds = proc.time()[["elapsed"]] - total_start))
  built <- build_combatrefql_result(prepared, dge, design, selected, fit,
    transported, keep, timings, call, normalization, chunk_size)
  if (nrow(built$warnings)) reporter$warning(built$warnings$message[[1L]])
  reporter$success("Constructed and validated result")
  reporter$finish(timings$seconds[timings$stage == "total"])
  built$result
}
