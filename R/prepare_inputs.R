empty_input_actions <- function() data.frame(
  field = character(), original_class = character(), resulting_class = character(),
  action = character(), reason = character(), affected = character(),
  lossy = logical(), stringsAsFactors = FALSE
)

input_action <- function(field, original, resulting, action, reason,
                         affected, lossy = FALSE) data.frame(
  field = field, original_class = original, resulting_class = resulting,
  action = action, reason = reason, affected = as.character(affected),
  lossy = lossy, stringsAsFactors = FALSE
)

is_default_rownames <- function(x) {
  rn <- rownames(x)
  is.null(rn) || identical(rn, as.character(seq_len(nrow(x))))
}

lossless_numeric <- function(x, field) {
  if (is.numeric(x)) return(as.numeric(x))
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) cli::cli_abort("{.arg {field}} contains a non-numeric column that cannot be repaired safely.")
  out <- suppressWarnings(as.numeric(x))
  if (anyNA(out) || any(!is.finite(out))) {
    cli::cli_abort("{.arg {field}} contains values that cannot be parsed losslessly as finite numeric counts.")
  }
  out
}

prepare_counts_input <- function(counts, batch = NULL,
                                 fractional_counts = c("error", "round"),
                                 duplicate_genes = c("error", "sum"),
                                 input_scale = "counts") {
  scale <- tolower(input_scale[[1L]])
  declared <- attr(counts, "expression_scale")
  if (!is.null(declared)) scale <- tolower(as.character(declared[[1L]]))
  if (scale %in% c("tpm", "cpm", "logcpm", "rpkm", "fpkm", "log", "normalized")) {
    cli::cli_abort("Declared {.arg input_scale} {.val {scale}} is not raw counts; TPM, CPM, logCPM, RPKM, and FPKM are unsafe inputs.")
  }
  if (scale != "counts") cli::cli_abort("Unknown {.arg input_scale} {.val {input_scale}}; use raw {.val counts} only.")
  fractional_counts <- match.arg(fractional_counts)
  duplicate_genes <- match.arg(duplicate_genes)
  actions <- empty_input_actions(); warnings <- character()
  original_class <- paste(class(counts), collapse = "/")

  if (is.data.frame(counts)) {
    source_rownames <- rownames(counts)
    if (ncol(counts) < 1L) cli::cli_abort("{.arg counts} has no columns.")
    if (is_default_rownames(counts) && ncol(counts) > 1L) {
      ids <- counts[[1L]]
      remaining <- counts[-1L]
      usable_ids <- !anyNA(ids) && all(nzchar(as.character(ids))) && !anyDuplicated(ids)
      numeric_remaining <- tryCatch({
        remaining[] <- lapply(seq_along(remaining), function(i)
          lossless_numeric(remaining[[i]], paste0("counts column ", names(remaining)[i])))
        TRUE
      }, error = function(e) FALSE)
      if (usable_ids && numeric_remaining) {
        rownames(remaining) <- as.character(ids)
        source_rownames <- rownames(remaining)
        actions <- rbind(actions, input_action(
          "counts", original_class, "data.frame", "first_column_to_gene_ids",
          paste0("used `", names(counts)[1L], "` as row names"), length(ids)))
        counts <- remaining
      }
    }
    converted <- lapply(seq_along(counts), function(i)
      lossless_numeric(counts[[i]], paste0("counts column ", names(counts)[i])))
    names(converted) <- names(counts)
    counts <- do.call(cbind, converted)
    colnames(counts) <- names(converted)
    if (!is.null(source_rownames)) rownames(counts) <- source_rownames
    actions <- rbind(actions, input_action(
      "counts", original_class, "matrix/array", "data_frame_to_numeric_matrix",
      "all expression columns converted losslessly", paste(dim(counts), collapse = " x ")))
  } else if (!is.matrix(counts)) {
    cli::cli_abort("{.arg counts} must be a matrix or data frame.")
  }

  if (!is.numeric(counts)) {
    converted <- suppressWarnings(matrix(as.numeric(counts), nrow(counts), ncol(counts),
      dimnames = dimnames(counts)))
    if (anyNA(converted) || any(!is.finite(converted))) {
      cli::cli_abort("{.arg counts} contains values that cannot be parsed losslessly as numeric.")
    }
    counts <- converted
    actions <- rbind(actions, input_action("counts", original_class, "matrix/array",
      "character_to_numeric", "every value parsed without new missing values", length(counts)))
  }
  if (is.integer(counts)) {
    storage.mode(counts) <- "double"
    actions <- rbind(actions, input_action("counts", "integer matrix", "double matrix",
      "integer_to_double", "edgeR uses numeric doubles", length(counts)))
  }
  if (!is.null(batch) && length(batch) == nrow(counts) && length(batch) != ncol(counts)) {
    counts <- t(counts)
    actions <- rbind(actions, input_action("counts", original_class, "matrix/array",
      "transpose", "batch length matched rows but not columns", paste(dim(counts), collapse = " x ")))
  }
  if (anyNA(counts) || any(!is.finite(counts))) cli::cli_abort("{.arg counts} contains NA, NaN, or infinite values; automatic repair is unsafe.")
  if (any(counts < 0)) cli::cli_abort("{.arg counts} contains negative values; raw counts must be non-negative.")
  fractional <- abs(counts - round(counts)) > sqrt(.Machine$double.eps)
  if (any(fractional)) {
    if (fractional_counts == "error") {
      cli::cli_abort("{.arg counts} contains {sum(fractional)} fractional values. Use {.code fractional_counts = 'round'} only when rounding is scientifically justified.")
    }
    difference <- abs(counts[fractional] - round(counts[fractional]))
    counts <- round(counts)
    reason <- sprintf("rounded %d values (%.3f%%); median |difference| %.6g; max %.6g",
      length(difference), 100 * length(difference) / length(counts),
      stats::median(difference), max(difference))
    actions <- rbind(actions, input_action("counts", "double matrix", "double matrix",
      "round_fractional_counts", reason, length(difference), TRUE))
    warnings <- c(warnings, paste("Fractional counts were rounded:", reason,
      "Do not use this for TPM, CPM, logCPM, RPKM, or FPKM."))
  }
  if (is.null(rownames(counts)) || anyNA(rownames(counts)) || any(rownames(counts) == "")) {
    cli::cli_abort("{.arg counts} needs unambiguous gene identifiers in row names or a unique first data-frame column.")
  }
  if (anyDuplicated(rownames(counts))) {
    if (duplicate_genes == "error") {
      cli::cli_abort("{.arg counts} has duplicated gene identifiers. Use {.code duplicate_genes = 'sum'} to aggregate exact identifiers explicitly.")
    }
    ids <- rownames(counts); unique_ids <- unique(ids)
    duplicate_ids <- sum(duplicated(ids) | duplicated(ids, fromLast = TRUE))
    counts <- rowsum(counts, factor(ids, levels = unique_ids), reorder = FALSE)
    actions <- rbind(actions, input_action("counts", "matrix/array", "matrix/array",
      "sum_duplicate_genes", sprintf("aggregated %d duplicated rows into exact identifiers", duplicate_ids),
      duplicate_ids, TRUE))
    warnings <- c(warnings, sprintf("Aggregated %d duplicated gene rows by exact identifier.", duplicate_ids))
  }
  if (is.null(colnames(counts)) || anyNA(colnames(counts)) || any(colnames(counts) == "")) cli::cli_abort("{.arg counts} must have sample names in column names.")
  if (anyDuplicated(colnames(counts))) cli::cli_abort("{.arg counts} has duplicated sample identifiers; automatic repair is unsafe.")
  storage.mode(counts) <- "double"
  list(value = counts, actions = actions, warnings = warnings)
}

prepare_sample_vector <- function(x, counts, field, allow_null = FALSE) {
  actions <- empty_input_actions()
  if (is.null(x)) {
    if (allow_null) return(list(value = NULL, actions = actions, warnings = character()))
    cli::cli_abort("{.arg {field}} must not be NULL.")
  }
  if (is.matrix(x) || is.data.frame(x)) cli::cli_abort("{.arg {field}} must be a sample vector.")
  if (length(x) != ncol(counts)) cli::cli_abort("{.arg {field}} length must equal the number of samples after orientation checks.")
  if (!is.null(names(x)) && !identical(names(x), colnames(counts))) {
    if (setequal(names(x), colnames(counts)) && !anyDuplicated(names(x))) {
      x <- x[colnames(counts)]
      actions <- rbind(actions, input_action(field, class(x)[1], class(x)[1],
        "reorder_by_sample_id", "names matched count samples in a different order", length(x)))
    } else cli::cli_abort("Names on {.arg {field}} do not match count sample identifiers.")
  }
  if (anyNA(x)) cli::cli_abort("{.arg {field}} contains missing labels.")
  original <- class(x)[1]; x <- droplevels(factor(x))
  actions <- rbind(actions, input_action(field, original, "factor", "convert_to_factor",
    paste0("converted `", field, "` to a factor"), length(x)))
  if (field == "group" && nlevels(x) == 1L) {
    actions <- rbind(actions, input_action(field, "factor", "NULL", "drop_one_level_group",
      "a one-level group adds no estimable information", length(x)))
    x <- NULL
  }
  list(value = x, actions = actions, warnings = character())
}

prepare_covariates <- function(covariates, counts) {
  actions <- empty_input_actions()
  if (is.null(covariates)) return(list(value = NULL, actions = actions, warnings = character()))
  if (!is.data.frame(covariates) && !is.matrix(covariates)) cli::cli_abort("{.arg covariates} must be a matrix or data frame.")
  covariates <- as.data.frame(covariates, check.names = FALSE)
  if (nrow(covariates) != ncol(counts)) cli::cli_abort("{.arg covariates} must have one row per sample.")
  if (!is.null(rownames(covariates)) && !identical(rownames(covariates), colnames(counts))) {
    if (setequal(rownames(covariates), colnames(counts)) && !anyDuplicated(rownames(covariates))) {
      covariates <- covariates[colnames(counts), , drop = FALSE]
      actions <- rbind(actions, input_action("covariates", "data.frame", "data.frame",
        "reorder_by_sample_id", "row names matched in a different order", nrow(covariates)))
    } else cli::cli_abort("Covariate row names do not match count sample identifiers.")
  }
  if (anyNA(covariates)) cli::cli_abort("{.arg covariates} contains missing values.")
  list(value = covariates, actions = actions, warnings = character())
}

canonical_method <- function(x, choices, field) {
  key <- gsub("[-_]", "", tolower(x[[1L]]))
  choice_keys <- gsub("[-_]", "", tolower(choices))
  match <- match(key, choice_keys)
  if (is.na(match)) cli::cli_abort("Unknown {.arg {field}} {.val {x}}; choose one of {paste(choices, collapse = ', ')}.")
  value <- choices[[match]]
  actions <- empty_input_actions()
  if (!identical(x[[1L]], value)) actions <- input_action(field, "character", "character",
    "canonicalize_method", paste0("converted `", x[[1L]], "` to `", value, "`"), 1L)
  list(value = value, actions = actions, warnings = character())
}

prepare_combat_inputs <- function(counts, batch, group = NULL, covariates = NULL,
                                  reference_batch = NULL,
                                  fractional_counts = "error",
                                  duplicate_genes = "error",
                                  input_scale = "counts") {
  pc <- prepare_counts_input(counts, batch, fractional_counts, duplicate_genes, input_scale)
  pb <- prepare_sample_vector(batch, pc$value, "batch")
  pg <- prepare_sample_vector(group, pc$value, "group", allow_null = TRUE)
  pv <- prepare_covariates(covariates, pc$value)
  reference <- if (is.null(reference_batch)) NULL else as.character(reference_batch[[1L]])
  list(counts = pc$value, batch = pb$value, group = pg$value,
       covariates = pv$value, reference_batch = reference,
       actions = do.call(rbind, list(pc$actions, pb$actions, pg$actions, pv$actions)),
       warnings = unique(c(pc$warnings, pb$warnings, pg$warnings, pv$warnings)))
}
