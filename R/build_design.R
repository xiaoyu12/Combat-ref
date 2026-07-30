build_combat_design <- function(batch, group = NULL, covariates = NULL) {
  data <- data.frame(batch = factor(batch))
  terms <- "0 + batch"
  if (!is.null(group) && nlevels(factor(group)) > 1L) {
    data$group <- factor(group)
    terms <- paste(terms, "+ group")
  }
  if (!is.null(covariates)) {
    covariates <- as.data.frame(covariates, check.names = FALSE)
    names(covariates) <- make.names(names(covariates), unique = TRUE)
    data <- cbind(data, covariates)
    terms <- paste(terms, "+", paste(names(covariates), collapse = " + "))
  }
  design <- stats::model.matrix(stats::as.formula(paste("~", terms)), data)
  colnames(design)[seq_len(nlevels(data$batch))] <- paste0("batch", levels(data$batch))
  design
}

counterfactual_reference_design <- function(design, batch, reference_batch) {
  out <- design
  batch_columns <- paste0("batch", levels(batch))
  out[, batch_columns] <- 0
  out[, paste0("batch", reference_batch)] <- 1
  out
}

build_within_batch_design <- function(group = NULL, covariates = NULL,
                                      sample_index) {
  n <- sum(sample_index)
  data <- data.frame(row.names = seq_len(n))
  missing_levels <- character()
  terms <- character()
  if (!is.null(group)) {
    global <- levels(factor(group))
    local <- droplevels(factor(group[sample_index], levels = global))
    missing_levels <- setdiff(global, levels(local))
    if (nlevels(local) > 1L) {
      data$group <- local
      terms <- c(terms, "group")
    }
  }
  if (!is.null(covariates)) {
    local_covariates <- as.data.frame(covariates[sample_index, , drop = FALSE],
                                      check.names = FALSE)
    local_covariates[] <- lapply(local_covariates, function(x) {
      if (is.factor(x) || is.character(x)) droplevels(factor(x)) else x
    })
    keep <- vapply(local_covariates, function(x) {
      if (is.factor(x)) nlevels(x) > 1L else length(unique(x)) > 1L
    }, logical(1))
    data <- cbind(data, local_covariates[, keep, drop = FALSE])
    terms <- c(terms, names(local_covariates)[keep])
  }
  formula <- if (length(terms)) paste("~", paste(terms, collapse = " + ")) else "~ 1"
  full <- stats::model.matrix(stats::as.formula(formula), data)
  keep_columns <- seq_len(ncol(full))
  dropped <- character()
  while (qr(full[, keep_columns, drop = FALSE])$rank < length(keep_columns) ||
         n <= length(keep_columns)) {
    if (length(keep_columns) == 1L) break
    dropped <- c(colnames(full)[utils::tail(keep_columns, 1L)], dropped)
    keep_columns <- utils::head(keep_columns, -1L)
  }
  design <- full[, keep_columns, drop = FALSE]
  rank <- qr(design)$rank
  residual_df <- nrow(design) - ncol(design)
  fallback <- length(dropped) > 0L || length(missing_levels) > 0L
  reasons <- c(
    if (length(missing_levels)) paste("absent local levels:", paste(missing_levels, collapse = ", ")),
    if (length(dropped)) paste("dropped non-estimable columns:", paste(dropped, collapse = ", ")),
    if (residual_df <= 0L) "no residual degrees of freedom"
  )
  list(
    design = design,
    formula = formula,
    columns_dropped = dropped,
    missing_levels = missing_levels,
    fallback = fallback,
    reason = if (length(reasons)) paste(reasons, collapse = "; ") else NA_character_,
    rank = rank,
    residual_df = residual_df
  )
}
