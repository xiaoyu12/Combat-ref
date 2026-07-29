new_combatref_result <- function(adjusted_counts, reference_batch, batch, group,
                                 removed_genes, diagnostics, model = NULL,
                                 call) {
  out <- list(
    adjusted_counts = adjusted_counts,
    reference_batch = reference_batch,
    batch = batch,
    group = group,
    removed_genes = removed_genes,
    diagnostics = diagnostics,
    call = call
  )
  if (!is.null(model)) {
    out$model <- model
  }
  structure(out, class = "combatref_result")
}

#' Print a ComBat-ref result
#'
#' @param x A `combatref_result` object.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
print.combatref_result <- function(x, ...) {
  cli::cli_text("<ComBat-ref result>")
  cli::cli_text("Genes:          {nrow(x$adjusted_counts)}")
  cli::cli_text("Samples:        {ncol(x$adjusted_counts)}")
  cli::cli_text("Batches:        {length(x$diagnostics$batch_sizes)}")
  cli::cli_text("Reference:      {x$reference_batch}")
  cli::cli_text("Filtered genes: {length(x$removed_genes)}")
  invisible(x)
}

#' Summarize a ComBat-ref result
#'
#' @param object A `combatref_result` object.
#' @param ... Unused.
#' @return A list with compact summary fields.
#' @export
summary.combatref_result <- function(object, ...) {
  out <- list(
    n_genes = object$diagnostics$n_genes_adjusted,
    n_samples = object$diagnostics$n_samples,
    reference_batch = object$reference_batch,
    removed_genes = object$removed_genes,
    diagnostics = object$diagnostics
  )
  class(out) <- "summary.combatref_result"
  out
}

#' @export
print.summary.combatref_result <- function(x, ...) {
  cli::cli_text("<ComBat-ref summary>")
  cli::cli_text("Genes:          {x$n_genes}")
  cli::cli_text("Samples:        {x$n_samples}")
  cli::cli_text("Reference:      {x$reference_batch}")
  cli::cli_text("Filtered genes: {length(x$removed_genes)}")
  invisible(x)
}
