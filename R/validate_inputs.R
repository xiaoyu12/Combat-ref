validate_combat_inputs <- function(counts, batch, group = NULL, covariates = NULL,
                                   reference_batch = NULL, lib_size = NULL) {
  if (!is.matrix(counts)) {
    counts <- as.matrix(counts)
  }
  if (!is.numeric(counts)) {
    cli::cli_abort(c(
      "x" = "{.arg counts} must contain raw, non-negative integer counts.",
      "i" = "TPM, CPM, logCPM and other normalized expression values are not supported."
    ))
  }
  if (anyNA(counts) || any(!is.finite(counts))) {
    cli::cli_abort("{.arg counts} must not contain missing or non-finite values.")
  }
  if (any(counts < 0)) {
    cli::cli_abort("{.arg counts} must not contain negative values.")
  }
  if (any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
    cli::cli_abort(c(
      "x" = "{.arg counts} must contain raw, non-negative integer counts.",
      "i" = "TPM, CPM, logCPM and other normalized expression values are not supported."
    ))
  }
  if (is.null(rownames(counts)) || anyNA(rownames(counts)) || any(rownames(counts) == "")) {
    cli::cli_abort("{.arg counts} must have non-missing gene names in {.code rownames(counts)}.")
  }
  if (anyDuplicated(rownames(counts))) {
    cli::cli_abort("{.arg counts} must not have duplicated gene names.")
  }
  if (is.null(colnames(counts)) || anyNA(colnames(counts)) || any(colnames(counts) == "")) {
    cli::cli_abort("{.arg counts} must have non-missing sample names in {.code colnames(counts)}.")
  }
  if (anyDuplicated(colnames(counts))) {
    cli::cli_abort("{.arg counts} must not have duplicated sample names.")
  }
  if (any(colSums(counts) == 0)) {
    cli::cli_abort("{.arg counts} must not contain all-zero samples.")
  }

  batch <- validate_sample_vector(batch, counts, "batch")
  if (anyNA(batch)) {
    cli::cli_abort("{.arg batch} must not contain missing values.")
  }
  batch <- factor(batch)
  if (nlevels(batch) < 2) {
    cli::cli_abort("{.arg batch} must contain at least two batches.")
  }
  if (any(tabulate(batch) < 2)) {
    cli::cli_abort("Each batch must contain at least two samples.")
  }

  if (!is.null(group)) {
    group <- validate_sample_vector(group, counts, "group")
    if (anyNA(group)) {
      cli::cli_abort("{.arg group} must not contain missing values.")
    }
    group <- factor(group)
    if (nlevels(group) == 1L) group <- NULL
  }

  if (!is.null(covariates)) {
    if (!is.data.frame(covariates) && !is.matrix(covariates)) {
      cli::cli_abort("{.arg covariates} must be a data frame or matrix aligned with samples.")
    }
    covariates <- as.data.frame(covariates, check.names = FALSE)
    if (nrow(covariates) != ncol(counts)) cli::cli_abort("{.arg covariates} must have one row per sample.")
    if (!is.null(rownames(covariates)) && !identical(rownames(covariates), colnames(counts))) {
      cli::cli_abort("Row names on {.arg covariates} must exactly match {.code colnames(counts)}.")
    }
    if (anyNA(covariates)) cli::cli_abort("{.arg covariates} must not contain missing values.")
  }

  if (!is.null(reference_batch) && !(reference_batch %in% levels(batch))) {
    cli::cli_abort("{.arg reference_batch} must be one of the observed batch levels.")
  }

  if (!is.null(lib_size)) {
    lib_size <- validate_sample_vector(lib_size, counts, "lib_size")
    if (!is.numeric(lib_size) || anyNA(lib_size) || any(!is.finite(lib_size)) || any(lib_size <= 0)) {
      cli::cli_abort("{.arg lib_size} must be finite, positive, and aligned with samples.")
    }
    lib_size <- as.numeric(lib_size)
  }

  design <- build_combat_design(batch = batch, group = group, covariates = covariates)
  if (qr(design)$rank < ncol(design)) {
    cli::cli_abort("The batch and group design is rank deficient or perfectly confounded.")
  }

  list(
    counts = counts,
    batch = batch,
    group = group,
    covariates = covariates,
    reference_batch = reference_batch,
    lib_size = lib_size,
    design = design
  )
}

validate_sample_vector <- function(x, counts, name) {
  if (is.null(x)) {
    cli::cli_abort("{.arg {name}} must not be {.code NULL}.")
  }
  if (is.matrix(x) || is.data.frame(x)) {
    cli::cli_abort("{.arg {name}} must be a vector aligned with columns of {.arg counts}.")
  }
  if (length(x) != ncol(counts)) {
    cli::cli_abort("{.arg {name}} must have length equal to {.code ncol(counts)}.")
  }

  nm <- names(x)
  if (!is.null(nm)) {
    if (anyNA(nm) || any(nm == "") || anyDuplicated(nm)) {
      cli::cli_abort("Names on {.arg {name}} must be non-missing and unique.")
    }
    if (!identical(nm, colnames(counts))) {
      cli::cli_abort("Names on {.arg {name}} must exactly match {.code colnames(counts)} in the same order.")
    }
  }
  x
}
