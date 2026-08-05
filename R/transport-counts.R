midp_transport <- function(counts, observed_mu, reference_mu,
                           source_dispersion, target_dispersion) {
  valid_means <- is.finite(observed_mu) & observed_mu > 0 &
    is.finite(reference_mu) & reference_mu > 0
  if (!all(valid_means)) {
    abort_combatrefql("combatrefql_mapping_error",
      "Negative-binomial transport received invalid fitted means.",
      "*" = "Inspect the affected genes for numerical overflow or a non-estimable fit.",
      stage = "mapping", affected_genes = sum(rowSums(!valid_means) > 0))
  }
  if (any(!is.finite(source_dispersion)) || any(source_dispersion < 0) ||
      any(!is.finite(target_dispersion)) || any(target_dispersion < 0))
    abort_combatrefql("combatrefql_mapping_error",
      "Negative-binomial transport received invalid dispersions.",
      stage = "mapping")
  source_size <- ifelse(source_dispersion == 0, Inf, 1 / source_dispersion)
  target_size <- ifelse(target_dispersion == 0, Inf, 1 / target_dispersion)
  probability <- stats::pnbinom(counts - 1, mu = observed_mu,
    size = source_size) + 0.5 * stats::dnbinom(counts, mu = observed_mu,
      size = source_size)
  clipped <- probability <= .Machine$double.xmin |
    probability >= 1 - .Machine$double.eps
  probability <- pmin(pmax(probability, .Machine$double.xmin), 1 - .Machine$double.eps)
  mapped <- stats::qnbinom(probability, mu = reference_mu, size = target_size)
  storage.mode(mapped) <- "double"
  valid <- rowSums(is.na(mapped) | !is.finite(mapped) | mapped < 0 |
    abs(mapped - round(mapped)) > sqrt(.Machine$double.eps)) == 0
  list(counts = mapped, probability = probability, valid = valid,
    tail_events = sum(clipped, na.rm = TRUE),
    zero_to_positive = sum(counts == 0 & mapped > 0, na.rm = TRUE),
    positive_to_zero = sum(counts > 0 & mapped == 0, na.rm = TRUE))
}

transport_counts <- function(counts, keep, batch, reference, fit,
                             shrinkage, dispersion, chunk_size,
                             progress = function(...) NULL) {
  adjusted <- counts
  mapping <- list()
  supported_genes <- which(keep)
  mapping_valid <- rep(TRUE, length(supported_genes))
  failures <- list()
  for (source in setdiff(levels(batch), reference)) {
    started <- proc.time()[["elapsed"]]
    source_valid <- rep(TRUE, length(supported_genes))
    sample_idx <- which(batch == source)
    genes <- seq_along(supported_genes)
    chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))
    tail_events <- 0L
    zero_to_positive <- positive_to_zero <- 0L
    for (chunk in chunks) {
      progress(source, length(chunk))
      observed_mu <- fit$fitted.values[chunk, sample_idx, drop = FALSE]
      delta <- shrinkage$posterior[[source]][chunk]
      reference_mu <- observed_mu * exp(-delta)
      source_phi <- dispersion$dispersion[chunk, sample_idx, drop = FALSE]
      target_phi <- matrix(dispersion$reference[chunk], length(chunk),
        length(sample_idx))
      target_genes <- supported_genes[chunk]
      valid_rows <- rowSums(!is.finite(observed_mu) | observed_mu <= 0 |
        !is.finite(reference_mu) | reference_mu <= 0) == 0
      if (any(!valid_rows)) {
        mapping_valid[chunk[!valid_rows]] <- FALSE
        source_valid[chunk[!valid_rows]] <- FALSE
        failures[[length(failures) + 1L]] <- data.frame(
          gene = rownames(counts)[target_genes[!valid_rows]], source_batch = source,
          class = "combatrefql_mapping_error",
          message = "Counterfactual prediction produced invalid fitted means.",
          stringsAsFactors = FALSE)
      }
      if (!any(valid_rows)) next
      mapped <- tryCatch(
        midp_transport(counts[target_genes[valid_rows], sample_idx, drop = FALSE],
          observed_mu[valid_rows, , drop = FALSE],
          reference_mu[valid_rows, , drop = FALSE],
          source_phi[valid_rows, , drop = FALSE],
          target_phi[valid_rows, , drop = FALSE]),
        combatrefql_mapping_error = function(error) error
      )
      if (inherits(mapped, "combatrefql_mapping_error")) {
        mapping_valid[chunk[valid_rows]] <- FALSE
        source_valid[chunk[valid_rows]] <- FALSE
        failures[[length(failures) + 1L]] <- data.frame(
          gene = rownames(counts)[target_genes[valid_rows]], source_batch = source,
          class = "combatrefql_mapping_error", message = conditionMessage(mapped),
          stringsAsFactors = FALSE)
      } else {
        valid_targets <- target_genes[valid_rows]
        if (any(!mapped$valid)) {
          failed_positions <- chunk[which(valid_rows)[!mapped$valid]]
          mapping_valid[failed_positions] <- FALSE
          source_valid[failed_positions] <- FALSE
          failures[[length(failures) + 1L]] <- data.frame(
            gene = rownames(counts)[valid_targets[!mapped$valid]],
            source_batch = source, class = "combatrefql_mapping_error",
            message = "NB transport produced invalid adjusted counts.",
            stringsAsFactors = FALSE)
        }
        adjusted[valid_targets[mapped$valid], sample_idx] <-
          mapped$counts[mapped$valid, , drop = FALSE]
        tail_events <- tail_events + mapped$tail_events
        zero_to_positive <- zero_to_positive + mapped$zero_to_positive
        positive_to_zero <- positive_to_zero + mapped$positive_to_zero
      }
    }
    mapping[[source]] <- data.frame(
      source_batch = source, reference_batch = reference,
      genes_attempted = sum(keep),
      genes_adjusted = sum(source_valid), genes_failed = sum(!source_valid),
      source_dispersion_min = min(dispersion$dispersion[, sample_idx]),
      source_dispersion_max = max(dispersion$dispersion[, sample_idx]),
      target_dispersion_min = min(dispersion$reference),
      target_dispersion_max = max(dispersion$reference),
      zero_to_positive = zero_to_positive,
      positive_to_zero = positive_to_zero, tail_events = tail_events,
      reference_invariant = identical(adjusted[, batch == reference, drop = FALSE],
        counts[, batch == reference, drop = FALSE]),
      seconds = proc.time()[["elapsed"]] - started, stringsAsFactors = FALSE)
  }
  failed_genes <- supported_genes[!mapping_valid]
  if (length(failed_genes)) adjusted[failed_genes, ] <- counts[failed_genes, ]
  if (!identical(adjusted[, batch == reference, drop = FALSE],
      counts[, batch == reference, drop = FALSE]))
    abort_combatrefql("combatrefql_internal_error",
      "Reference-batch counts changed during transport.", stage = "mapping")
  list(counts = adjusted, mapping = do.call(rbind, mapping),
       mapping_valid = mapping_valid,
       failures = if (length(failures)) do.call(rbind, failures) else
         data.frame(gene = character(), source_batch = character(),
           class = character(), message = character()))
}
