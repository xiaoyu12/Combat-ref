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
#' @param return_model Logical. If `TRUE`, include selected internal matrices in
#'   `diagnostics`. Large edgeR fit objects are not returned.
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
                       lib_size = NULL, return_model = FALSE) {
  call <- match.call()
  validated <- validate_combat_inputs(
    counts = counts,
    batch = batch,
    group = group,
    reference_batch = reference_batch,
    lib_size = lib_size
  )

  adjusted <- ComBat_ref(
    counts = validated$counts,
    batch = validated$batch,
    group = validated$group,
    reference_batch = validated$reference_batch,
    lib_size = validated$lib_size,
    genewise.disp = FALSE
  )
  selected_reference <- attr(adjusted, "reference_batch")
  attr(adjusted, "reference_batch") <- NULL

  removed_genes <- rownames(validated$counts)[
    vapply(seq_len(nrow(validated$counts)), function(i) {
      any(tapply(validated$counts[i, ], validated$batch, function(x) all(x == 0)))
    }, logical(1))
  ]

  diagnostics <- list(
    n_genes = nrow(validated$counts),
    n_samples = ncol(validated$counts),
    n_removed_genes = length(removed_genes)
  )
  if (return_model) {
    diagnostics$design <- validated$design
  }

  new_combatref_result(
    adjusted_counts = adjusted,
    reference_batch = selected_reference,
    batch = validated$batch,
    group = validated$group,
    design = validated$design,
    dispersions = NULL,
    removed_genes = removed_genes,
    diagnostics = diagnostics,
    call = call
  )
}
