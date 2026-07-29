validate_combat_inputs <- function(counts, batch, group = NULL, reference_batch = NULL,
                                   lib_size = NULL) {
  if (!is.matrix(counts)) {
    counts <- as.matrix(counts)
  }
  if (!is.numeric(counts)) {
    stop("`counts` must be numeric.", call. = FALSE)
  }
  if (anyNA(counts) || any(!is.finite(counts))) {
    stop("`counts` must not contain missing or non-finite values.", call. = FALSE)
  }
  if (any(counts < 0)) {
    stop("`counts` must not contain negative values.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || anyNA(rownames(counts)) || any(rownames(counts) == "")) {
    stop("`counts` must have non-missing gene names in `rownames(counts)`.", call. = FALSE)
  }
  if (anyDuplicated(rownames(counts))) {
    stop("`counts` must not have duplicated gene names.", call. = FALSE)
  }
  if (is.null(colnames(counts)) || anyNA(colnames(counts)) || any(colnames(counts) == "")) {
    stop("`counts` must have non-missing sample names in `colnames(counts)`.", call. = FALSE)
  }
  if (anyDuplicated(colnames(counts))) {
    stop("`counts` must not have duplicated sample names.", call. = FALSE)
  }
  if (any(colSums(counts) == 0)) {
    stop("`counts` must not contain all-zero samples.", call. = FALSE)
  }

  batch <- validate_sample_vector(batch, counts, "batch")
  if (anyNA(batch)) {
    stop("`batch` must not contain missing values.", call. = FALSE)
  }
  batch <- factor(batch)
  if (nlevels(batch) < 2) {
    stop("`batch` must contain at least two batches.", call. = FALSE)
  }
  if (any(tabulate(batch) < 2)) {
    stop("Each batch must contain at least two samples.", call. = FALSE)
  }

  if (!is.null(group)) {
    group <- validate_sample_vector(group, counts, "group")
    if (anyNA(group)) {
      stop("`group` must not contain missing values.", call. = FALSE)
    }
    group <- factor(group)
  }

  if (!is.null(reference_batch) && !(reference_batch %in% levels(batch))) {
    stop("`reference_batch` must be one of the observed batch levels.", call. = FALSE)
  }

  if (!is.null(lib_size)) {
    lib_size <- validate_sample_vector(lib_size, counts, "lib_size")
    if (!is.numeric(lib_size) || anyNA(lib_size) || any(!is.finite(lib_size)) || any(lib_size <= 0)) {
      stop("`lib_size` must be finite, positive, and aligned with samples.", call. = FALSE)
    }
    lib_size <- as.numeric(lib_size)
  }

  design <- build_combat_design(batch = batch, group = group)
  if (qr(design)$rank < ncol(design)) {
    stop("The batch and group design is rank deficient or perfectly confounded.", call. = FALSE)
  }

  list(
    counts = counts,
    batch = batch,
    group = group,
    reference_batch = reference_batch,
    lib_size = lib_size,
    design = design
  )
}

validate_sample_vector <- function(x, counts, name) {
  if (is.null(x)) {
    stop(sprintf("`%s` must not be NULL.", name), call. = FALSE)
  }
  if (is.matrix(x) || is.data.frame(x)) {
    stop(sprintf("`%s` must be a vector aligned with columns of `counts`.", name), call. = FALSE)
  }
  if (length(x) != ncol(counts)) {
    stop(sprintf("`%s` must have length equal to `ncol(counts)`.", name), call. = FALSE)
  }

  nm <- names(x)
  if (!is.null(nm)) {
    if (anyNA(nm) || any(nm == "") || anyDuplicated(nm)) {
      stop(sprintf("Names on `%s` must be non-missing and unique.", name), call. = FALSE)
    }
    if (!identical(nm, colnames(counts))) {
      stop(sprintf("Names on `%s` must exactly match `colnames(counts)` in the same order.", name), call. = FALSE)
    }
  }
  x
}
