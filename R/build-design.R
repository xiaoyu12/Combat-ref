build_design <- function(batch, group = NULL, covariates = NULL) {
  data <- data.frame(batch = batch)
  rhs <- "0 + batch"
  if (!is.null(group)) {
    data$group <- group
    rhs <- paste(rhs, "+ group")
  }
  if (!is.null(covariates)) {
    names(covariates) <- make.names(names(covariates), unique = TRUE)
    data <- cbind(data, covariates)
    rhs <- paste(rhs, "+", paste(names(covariates), collapse = " + "))
  }
  design <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data)
  colnames(design)[seq_len(nlevels(batch))] <- paste0("batch", levels(batch))
  rank <- qr(design)$rank
  residual_df <- nrow(design) - rank
  if (rank < ncol(design) || residual_df <= 0L) {
    abort_combatrefql("combatrefql_design_error",
      "The batch, group, and covariate design is rank deficient or has no residual degrees of freedom.",
      "*" = "Remove confounded or redundant terms, or add replicates.",
      stage = "design", rank = rank, residual_df = residual_df)
  }
  list(matrix = design, rank = rank, residual_df = residual_df,
       condition_number = kappa(design))
}

reference_design <- function(design, batch, reference, columns) {
  out <- design
  out[, columns] <- 0
  out[, paste0("batch", reference)] <- 1
  out
}

