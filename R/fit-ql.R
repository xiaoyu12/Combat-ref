fit_ql_model <- function(dge, design) {
  fit <- tryCatch(
    edgeR::glmQLFit(dge, design, dispersion = NULL, robust = TRUE, legacy = FALSE),
    error = function(error) abort_combatrefql("combatrefql_ql_error",
      sprintf("edgeR v4 quasi-likelihood fitting failed: %s", conditionMessage(error)),
      "*" = "Check replication, design rank, and count scale.", stage = "ql_fit")
  )
  if (length(fit$dispersion) != 1L || !is.finite(fit$dispersion) || fit$dispersion < 0) {
    abort_combatrefql("combatrefql_internal_error",
      "edgeR did not return one valid global negative-binomial dispersion.", stage = "ql_fit")
  }
  fit
}
