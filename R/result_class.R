new_combatref_result <- function(adjusted_counts, reference_batch, batch, group,
                                 design, dispersions, removed_genes,
                                 diagnostics, call) {
  structure(
    list(
      adjusted_counts = adjusted_counts,
      reference_batch = reference_batch,
      batch = batch,
      group = group,
      design = design,
      dispersions = dispersions,
      removed_genes = removed_genes,
      diagnostics = diagnostics,
      call = call
    ),
    class = "combatref_result"
  )
}

#' Print a ComBat-ref result
#'
#' @param x A `combatref_result` object.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
print.combatref_result <- function(x, ...) {
  cat("ComBat-ref result\n")
  cat("  genes:", nrow(x$adjusted_counts), "\n")
  cat("  samples:", ncol(x$adjusted_counts), "\n")
  cat("  reference batch:", x$reference_batch, "\n")
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
    n_genes = nrow(object$adjusted_counts),
    n_samples = ncol(object$adjusted_counts),
    reference_batch = object$reference_batch,
    removed_genes = object$removed_genes,
    diagnostics = object$diagnostics
  )
  class(out) <- "summary.combatref_result"
  out
}

#' @export
print.summary.combatref_result <- function(x, ...) {
  cat("ComBat-ref summary\n")
  cat("  genes:", x$n_genes, "\n")
  cat("  samples:", x$n_samples, "\n")
  cat("  reference batch:", x$reference_batch, "\n")
  cat("  removed genes:", length(x$removed_genes), "\n")
  invisible(x)
}
