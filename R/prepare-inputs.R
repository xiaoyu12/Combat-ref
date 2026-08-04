prepare_inputs <- function(counts, batch, group, covariates, reference,
                           fractional_counts) {
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (!is.matrix(counts) || !is.numeric(counts)) {
    input_error("{.arg counts} must be a numeric matrix of raw counts.",
                "*" = "Supply untransformed gene-by-sample counts.")
  }
  if (!length(counts) || any(dim(counts) == 0L)) input_error("{.arg counts} must not be empty.")
  if (anyNA(counts) || any(!is.finite(counts))) input_error("{.arg counts} contains missing or infinite values.")
  if (any(counts < 0)) input_error("{.arg counts} contains negative values.")
  validate_ids(rownames(counts), "gene", "rownames(counts)")
  validate_ids(colnames(counts), "sample", "colnames(counts)")
  if (any(colSums(counts) == 0)) input_error("{.arg counts} contains an all-zero sample.")

  delta <- round(counts) - counts
  fractional <- abs(delta) > sqrt(.Machine$double.eps)
  n_fractional <- sum(fractional)
  if (n_fractional && fractional_counts == "error") {
    input_error(sprintf("{.arg counts} contains %s fractional values.", n_fractional),
                "*" = "Use {.code fractional_counts = \"round\"} to audit and round them.",
                affected_values = n_fractional)
  }
  library_sizes <- colSums(counts)
  equal_library_sizes <- length(library_sizes) > 1L && mean(library_sizes) > 0 &&
    (max(library_sizes) - min(library_sizes)) / mean(library_sizes) < 1e-6
  dense_small_expression <- max(counts) < 100 && mean(counts == 0) < 0.01
  if (n_fractional && mean(fractional) > 0.9 &&
      (equal_library_sizes || dense_small_expression)) {
    input_error("{.arg counts} appears to be normalized expression rather than raw counts.",
                "*" = "Supply untransformed read counts; TPM, CPM, and log-expression are unsupported.")
  }
  action <- data.frame(
    stage = "preparation", action = "fractional_rounding",
    severity = if (n_fractional) "info" else "none", lossy = n_fractional > 0,
    affected_values = n_fractional,
    affected_genes = if (n_fractional) sum(rowSums(fractional) > 0) else 0L,
    maximum_absolute_change = if (n_fractional) max(abs(delta[fractional])) else 0,
    total_count_change = if (n_fractional) sum(delta[fractional]) else 0,
    message = if (n_fractional) sprintf("Rounded %s fractional values across %s genes", n_fractional, sum(rowSums(fractional) > 0)) else "No fractional values detected",
    stringsAsFactors = FALSE
  )
  counts <- round(counts)
  storage.mode(counts) <- "double"

  batch <- align_vector(batch, counts, "batch", required = TRUE)
  if (anyNA(batch) || any(batch == "")) input_error("{.arg batch} contains missing or empty values.")
  batch <- factor(batch)
  if (nlevels(batch) < 2L) input_error("{.arg batch} must contain at least two batches.")
  if (any(table(batch) < 2L)) input_error("Every batch must contain at least two samples.")
  group <- align_vector(group, counts, "group", required = FALSE)
  if (!is.null(group)) {
    if (anyNA(group) || any(group == "")) input_error("{.arg group} contains missing or empty values.")
    group <- factor(group)
    if (nlevels(group) == 1L) group <- NULL
  }
  covariates <- align_covariates(covariates, counts)
  if (!is.null(reference)) {
    if (length(reference) != 1L || is.na(reference) || !as.character(reference) %in% levels(batch)) {
      abort_combatrefql("combatrefql_reference_error",
        "{.arg reference} must name exactly one observed batch.",
        "*" = sprintf("Choose one of: %s.", paste(levels(batch), collapse = ", ")),
        stage = "reference", batch = reference)
    }
    reference <- as.character(reference)
  }
  list(counts = counts, batch = batch, group = group, covariates = covariates,
       reference = reference, input_actions = action)
}

validate_ids <- function(x, type, location) {
  if (is.null(x) || anyNA(x) || any(x == "") || anyDuplicated(x)) {
    input_error(sprintf("%s must contain unique, non-missing %s identifiers.", location, type))
  }
}

align_vector <- function(x, counts, name, required) {
  if (is.null(x)) {
    if (required) input_error(sprintf("%s must not be NULL.", name))
    return(NULL)
  }
  if (is.matrix(x) || is.data.frame(x) || length(x) != ncol(counts)) {
    abort_combatrefql("combatrefql_alignment_error",
      sprintf("%s must be a vector with one value per sample.", name),
      "*" = "Align it to {.code colnames(counts)}.", stage = "alignment")
  }
  nm <- names(x)
  if (!is.null(nm)) {
    if (anyNA(nm) || any(nm == "") || anyDuplicated(nm) || !setequal(nm, colnames(counts))) {
      abort_combatrefql("combatrefql_alignment_error",
        sprintf("Names on %s do not uniquely match colnames(counts).", name),
        "*" = "Use the same sample identifiers; order may differ.", stage = "alignment")
    }
    x <- x[match(colnames(counts), nm)]
  }
  unname(x)
}

align_covariates <- function(x, counts) {
  if (is.null(x)) return(NULL)
  if (!is.data.frame(x) && !is.matrix(x)) input_error("{.arg covariates} must be a data frame or matrix.")
  x <- as.data.frame(x, check.names = FALSE)
  if (nrow(x) != ncol(counts) || anyNA(x)) input_error("{.arg covariates} must have one complete row per sample.")
  if (!is.null(rownames(x))) {
    if (anyDuplicated(rownames(x)) || !setequal(rownames(x), colnames(counts))) {
      abort_combatrefql("combatrefql_alignment_error", "Covariate row names do not match sample identifiers.", stage = "alignment")
    }
    x <- x[match(colnames(counts), rownames(x)), , drop = FALSE]
  }
  x
}

normalize_counts <- function(inputs, method) {
  edgeR::normLibSizes(edgeR::DGEList(inputs$counts), method = method)
}

assess_support <- function(counts, batch) {
  supported <- Reduce(`&`, lapply(levels(batch), function(value) {
    rowSums(counts[, batch == value, drop = FALSE]) > 0
  }))
  if (sum(supported) < 2L) {
    abort_combatrefql("combatrefql_ql_error",
      "Fewer than two genes have count support in every batch.",
      stage = "ql_fit", affected_genes = sum(supported))
  }
  supported
}
