#' Fit the baseline ComBat-ref correction
#'
#' `combat_ref()` exposes the original ComBat-ref reference-batch correction
#' algorithm through an R-facing package API. The input matrix must contain genes
#' in rows and samples in columns. The current implementation preserves the
#' legacy negative-binomial fitting, dispersion estimation, automatic
#' minimum-dispersion reference-batch selection, filtering, zero handling, and
#' negative-binomial quantile mapping.
#'
#' This function does not implement quasi-likelihood correction.
#'
#' @param counts Numeric count matrix-like object with genes in rows and samples
#'   in columns. Row and column names are required and preserved.
#' @param batch Batch vector aligned exactly with columns of `counts`.
#' @param group Optional biological-group vector aligned exactly with columns of
#'   `counts`.
#' @param reference_batch Optional reference batch. `NULL` reproduces the legacy
#'   automatic minimum-dispersion selection.
#' @param lib_size Optional positive finite library sizes aligned with samples.
#'   `NULL` reproduces the legacy edgeR library-size behavior.
#' @param verbose Logical. If `TRUE`, print concise progress messages with
#'   `cli`. Defaults to `interactive()`.
#' @param keep_model Logical. If `TRUE`, include selected intermediate model
#'   quantities under the `model` field. Large edgeR fit objects are not
#'   returned.
#'
#' @return A `combatref_result` object. Use `$adjusted_counts` to access the
#'   corrected count matrix.
#' @export
#' @examples
#' counts <- matrix(
#'   c(10, 12, 30, 33, 8, 9, 28, 31,
#'     4, 5, 14, 16, 3, 4, 15, 17,
#'     50, 55, 120, 130, 48, 52, 125, 135,
#'     1, 2, 4, 5, 1, 2, 5, 6),
#'   nrow = 4,
#'   byrow = TRUE,
#'   dimnames = list(paste0("gene", 1:4), paste0("sample", 1:8))
#' )
#' batch <- rep(c("batch_A", "batch_B"), each = 2, times = 2)
#' group <- rep(c("ctrl", "treated"), each = 4)
#'
#' fit <- combat_ref(counts, batch = batch, group = group)
#' adjusted_counts <- fit$adjusted_counts
combat_ref <- function(counts, batch, group = NULL, reference_batch = NULL,
                       lib_size = NULL, verbose = interactive(),
                       keep_model = FALSE) {
  call <- match.call()
  if (isTRUE(verbose)) {
    cli::cli_alert_info("Validating input data")
  }
  validated <- validate_combat_inputs(
    counts = counts,
    batch = batch,
    group = group,
    reference_batch = reference_batch,
    lib_size = lib_size
  )
  batch_sizes <- table(validated$batch)
  group_sizes <- if (is.null(validated$group)) NULL else table(validated$group)

  if (isTRUE(verbose)) {
    cli::cli_alert_success(
      "Accepted {nrow(validated$counts)} genes across {ncol(validated$counts)} samples"
    )
    cli::cli_alert_info("Detected {length(batch_sizes)} batches: {format_named_counts(batch_sizes)}")
    if (!is.null(group_sizes)) {
      cli::cli_alert_info("Biological groups: {format_named_counts(group_sizes)}")
    }
    if (is.null(validated$reference_batch)) {
      cli::cli_alert_info("Selecting the reference batch automatically")
    } else {
      cli::cli_alert_info("Using requested reference batch: {validated$reference_batch}")
    }
    cli::cli_alert_info("Fitting the baseline ComBat-ref negative-binomial model")
  }

  utils::capture.output({
    adjusted <- ComBat_ref(
      counts = validated$counts,
      batch = validated$batch,
      group = validated$group,
      reference_batch = validated$reference_batch,
      lib_size = validated$lib_size,
      genewise.disp = FALSE
    )
  })
  selected_reference <- attr(adjusted, "reference_batch")
  model <- attr(adjusted, "model")
  attr(adjusted, "reference_batch") <- NULL
  attr(adjusted, "model") <- NULL

  removed_genes <- rownames(validated$counts)[
    vapply(seq_len(nrow(validated$counts)), function(i) {
      any(tapply(validated$counts[i, ], validated$batch, function(x) all(x == 0)))
    }, logical(1))
  ]

  diagnostics <- list(
    n_genes_input = nrow(validated$counts),
    n_genes_adjusted = nrow(adjusted),
    n_genes_filtered = length(removed_genes),
    n_samples = ncol(validated$counts),
    batch_sizes = batch_sizes,
    group_sizes = group_sizes,
    automatic_reference = is.null(validated$reference_batch),
    reference_batch = selected_reference
  )

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Selected reference batch: {selected_reference}")
    if (length(removed_genes) > 0) {
      cli::cli_alert_info(
        "Filtered {length(removed_genes)} genes that were all zero within at least one batch"
      )
    }
    cli::cli_alert_info("Mapping non-reference counts to the reference distribution")
    cli::cli_alert_success("Batch correction completed")
  }

  new_combatref_result(
    adjusted_counts = adjusted,
    reference_batch = selected_reference,
    batch = validated$batch,
    group = validated$group,
    removed_genes = removed_genes,
    diagnostics = diagnostics,
    model = if (isTRUE(keep_model)) model else NULL,
    call = call
  )
}

format_named_counts <- function(x) {
  paste0(names(x), " (n = ", as.integer(x), ")", collapse = ", ")
}
