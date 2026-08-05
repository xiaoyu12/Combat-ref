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

batch_contrast <- function(design, source, reference) {
  contrast <- stats::setNames(numeric(ncol(design)), colnames(design))
  source_column <- match(paste0("batch", source), names(contrast))
  reference_column <- match(paste0("batch", reference), names(contrast))
  if (!is.na(source_column)) contrast[source_column] <- 1
  if (!is.na(reference_column)) contrast[reference_column] <- -1
  if (is.na(source_column) && is.na(reference_column))
    abort_combatrefql("combatrefql_internal_error",
      "The fitted design does not identify the requested batch contrast.",
      stage = "contrast")
  contrast
}

effective_residual_information <- function(design, batch, fitted = NULL,
                                           dispersion = NULL, weights = NULL) {
  batches <- levels(batch)
  design_hat <- design %*% solve(crossprod(design), t(design))
  design_only <- vapply(batches, function(b)
    sum(1 - diag(design_hat)[batch == b]), numeric(1))
  if (is.null(fitted)) return(data.frame(batch = batches,
    design_residual_information = design_only))

  if (is.null(dim(dispersion))) dispersion <- matrix(dispersion,
    nrow(fitted), ncol(fitted))
  if (is.null(weights)) weights <- matrix(1, nrow(fitted), ncol(fitted))
  info <- matrix(0, nrow(fitted), length(batches),
    dimnames = list(rownames(fitted), batches))
  for (g in seq_len(nrow(fitted))) {
    working <- weights[g, ] * fitted[g, ] /
      (1 + dispersion[g, ] * fitted[g, ])
    root <- sqrt(pmax(working, 0))
    xw <- design * root
    inverse <- tryCatch(solve(crossprod(xw)), error = function(error) NULL)
    if (is.null(inverse)) next
    leverage <- rowSums((xw %*% inverse) * xw)
    info[g, ] <- vapply(batches, function(b)
      sum(1 - leverage[batch == b]), numeric(1))
  }
  info[info < sqrt(.Machine$double.eps)] <- 0
  precision <- ifelse(info > 0, 1 / trigamma(info / 2), 0)
  list(
    design_only = data.frame(batch = batches,
      design_residual_information = design_only),
    weighted = info,
    precision = precision,
    formula = "inverse_trigamma(weighted_residual_information / 2)"
  )
}
