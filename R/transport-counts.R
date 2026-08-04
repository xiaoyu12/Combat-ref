midp_transport <- function(counts, observed_mu, reference_mu, nb_dispersion,
                           preserve_zero = TRUE) {
  valid_means <- is.finite(observed_mu) & observed_mu > 0 &
    is.finite(reference_mu) & reference_mu > 0
  if (!all(valid_means)) {
    abort_combatrefql("combatrefql_mapping_error",
      "Negative-binomial transport received invalid fitted means.",
      "*" = "Inspect the affected genes for numerical overflow or a non-estimable fit.",
      stage = "mapping", affected_genes = sum(rowSums(!valid_means) > 0))
  }
  size <- if (nb_dispersion == 0) Inf else 1 / nb_dispersion
  probability <- stats::pnbinom(counts - 1, mu = observed_mu, size = size) +
    0.5 * stats::dnbinom(counts, mu = observed_mu, size = size)
  if (preserve_zero) probability[counts == 0] <- 0
  clipped <- probability <= .Machine$double.xmin |
    probability >= 1 - .Machine$double.eps
  probability <- pmin(pmax(probability, .Machine$double.xmin), 1 - .Machine$double.eps)
  mapped <- stats::qnbinom(probability, mu = reference_mu, size = size)
  if (preserve_zero) mapped[counts == 0] <- 0
  storage.mode(mapped) <- "double"
  valid <- rowSums(is.na(mapped) | !is.finite(mapped) | mapped < 0 |
    abs(mapped - round(mapped)) > sqrt(.Machine$double.eps)) == 0
  list(counts = mapped, valid = valid, tail_events = sum(clipped, na.rm = TRUE))
}

transport_counts <- function(counts, keep, batch, reference, design, fit,
                             chunk_size, progress = function(...) NULL) {
  adjusted <- counts
  nb_dispersion <- as.numeric(fit$dispersion)
  batch_columns <- paste0("batch", levels(batch))
  mapping <- list()
  supported_genes <- which(keep)
  mapping_valid <- rep(TRUE, length(supported_genes))
  failures <- list()
  for (source in setdiff(levels(batch), reference)) {
    started <- proc.time()[["elapsed"]]
    source_valid <- rep(TRUE, length(supported_genes))
    sample_idx <- which(batch == source)
    cf <- reference_design(design[sample_idx, , drop = FALSE], batch, reference, batch_columns)
    observed_design <- design[sample_idx, , drop = FALSE]
    genes <- seq_along(supported_genes)
    chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))
    tail_events <- 0L
    for (chunk in chunks) {
      progress(source, length(chunk))
      chunk_coefficients <- fit$coefficients[chunk, , drop = FALSE]
      chunk_offset <- fit$offset[chunk, sample_idx, drop = FALSE]
      observed_mu <- exp(chunk_offset + chunk_coefficients %*% t(observed_design))
      reference_mu <- exp(chunk_offset + chunk_coefficients %*% t(cf))
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
          reference_mu[valid_rows, , drop = FALSE], nb_dispersion),
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
      }
    }
    mapping[[source]] <- data.frame(
      source_batch = source, reference_batch = reference,
      genes_attempted = sum(keep),
      genes_adjusted = sum(source_valid), genes_failed = sum(!source_valid),
      tail_events = tail_events,
      seconds = proc.time()[["elapsed"]] - started, stringsAsFactors = FALSE)
  }
  failed_genes <- supported_genes[!mapping_valid]
  if (length(failed_genes)) adjusted[failed_genes, ] <- counts[failed_genes, ]
  list(counts = adjusted, mapping = do.call(rbind, mapping),
       mapping_valid = mapping_valid,
       failures = if (length(failures)) do.call(rbind, failures) else
         data.frame(gene = character(), source_batch = character(),
           class = character(), message = character()))
}
