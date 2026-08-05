make_gene_status <- function(counts, batch, keep, mapping_valid) {
  zero_support <- vapply(levels(batch), function(x) rowSums(counts[, batch == x, drop = FALSE]) == 0, logical(nrow(counts)))
  if (is.null(dim(zero_support))) zero_support <- matrix(zero_support, ncol = 1L)
  zero_labels <- apply(zero_support, 1L, function(x) paste(levels(batch)[x], collapse = ","))
  status <- rep("unsupported", nrow(counts))
  status[keep] <- ifelse(mapping_valid, "adjusted", "failed")
  reason <- rep("all zero in at least one batch", nrow(counts))
  reason[status == "adjusted"] <- NA_character_
  reason[status == "failed"] <- "negative-binomial transport failed; original row retained"
  data.frame(gene = rownames(counts), status = status, reason = reason,
    adjusted = status == "adjusted", zero_batches = zero_labels,
    fit_valid = keep, mapping_valid = keep & status == "adjusted",
    stringsAsFactors = FALSE)
}

make_outcomes <- function(gene_status) {
  levels <- c("adjusted", "unsupported", "failed")
  genes <- tabulate(match(gene_status$status, levels), nbins = length(levels))
  data.frame(outcome = levels, genes = genes,
    proportion = genes / nrow(gene_status), stringsAsFactors = FALSE)
}

build_diagnostics <- function(selected, transported, timings, input_actions,
                              warnings, gene_status, ql, shrinkage, dispersion) {
  list(
    batches = selected$batches,
    mapping = transported$mapping,
    mapping_failures = transported$failures,
    timing = timings,
    input_actions = input_actions,
    warnings = warnings,
    outcomes = make_outcomes(gene_status),
    ql = ql,
    batch_contrasts = shrinkage$diagnostics,
    dispersion = dispersion$diagnostics
  )
}

empty_warnings <- function() data.frame(stage = character(), class = character(), message = character())

timing_row <- function(stage, start) data.frame(stage = stage, seconds = proc.time()[["elapsed"]] - start)

`%||%` <- function(x, y) if (is.null(x)) y else x
