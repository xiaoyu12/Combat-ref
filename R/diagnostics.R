make_gene_status <- function(counts, batch, keep) {
  zero_support <- vapply(levels(batch), function(x) rowSums(counts[, batch == x, drop = FALSE]) == 0, logical(nrow(counts)))
  if (is.null(dim(zero_support))) zero_support <- matrix(zero_support, ncol = 1L)
  zero_labels <- apply(zero_support, 1L, function(x) paste(levels(batch)[x], collapse = ","))
  data.frame(
    gene = rownames(counts), status = ifelse(keep, "adjusted", "unsupported"),
    reason = ifelse(keep, NA_character_, "all zero in at least one batch"),
    adjusted = keep, zero_batches = zero_labels, fit_valid = keep,
    mapping_valid = keep, stringsAsFactors = FALSE
  )
}

empty_warnings <- function() data.frame(stage = character(), class = character(), message = character())

timing_row <- function(stage, start) data.frame(stage = stage, seconds = proc.time()[["elapsed"]] - start)

`%||%` <- function(x, y) if (is.null(x)) y else x

