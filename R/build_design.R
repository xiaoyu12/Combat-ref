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
