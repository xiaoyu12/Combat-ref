select_reference <- function(dge, batch, explicit = NULL, progress = function(...) NULL) {
  scores <- vapply(levels(batch), function(candidate) {
    progress(candidate)
    idx <- batch == candidate
    local <- dge[, idx, keep.lib.sizes = TRUE]
    keep <- rowSums(local$counts) > 0
    if (sum(keep) < 2L) return(Inf)
    fit <- try(edgeR::glmQLFit(local[keep, , keep.lib.sizes = TRUE],
      matrix(1, sum(idx), 1L), dispersion = NULL, robust = TRUE, legacy = FALSE), silent = TRUE)
    if (inherits(fit, "try-error")) Inf else as.numeric(fit$dispersion)
  }, numeric(1))
  if (all(!is.finite(scores))) {
    abort_combatrefql("combatrefql_reference_error",
      "No batch could be scored as a reference.",
      "*" = "Ensure each batch has replicated, non-zero raw counts.", stage = "reference")
  }
  selected <- explicit %||% names(which.min(scores))
  data <- data.frame(batch = names(scores), samples = as.integer(table(batch)[names(scores)]),
    reference_score = unname(scores), selected = names(scores) == selected,
    selection_method = if (is.null(explicit)) "automatic_minimum_nb_dispersion" else "explicit",
    stringsAsFactors = FALSE)
  list(reference = selected, batches = data)
}
