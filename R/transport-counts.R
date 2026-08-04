midp_transport <- function(counts, observed_mu, reference_mu, nb_dispersion,
                           preserve_zero = TRUE) {
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
  list(counts = mapped, tail_events = sum(clipped))
}

transport_counts <- function(counts, keep, batch, reference, design, fit,
                             chunk_size, progress = function(...) NULL) {
  adjusted <- counts
  observed_mu <- fit$fitted.values
  nb_dispersion <- as.numeric(fit$dispersion)
  batch_columns <- paste0("batch", levels(batch))
  mapping <- list()
  for (source in setdiff(levels(batch), reference)) {
    started <- proc.time()[["elapsed"]]
    sample_idx <- which(batch == source)
    cf <- reference_design(design[sample_idx, , drop = FALSE], batch, reference, batch_columns)
    reference_mu <- exp(fit$offset[, sample_idx, drop = FALSE] + fit$coefficients %*% t(cf))
    genes <- seq_len(sum(keep))
    chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))
    tail_events <- 0L
    for (chunk in chunks) {
      progress(source, length(chunk))
      mapped <- midp_transport(
        adjusted[keep, sample_idx, drop = FALSE][chunk, , drop = FALSE],
        observed_mu[chunk, sample_idx, drop = FALSE],
        reference_mu[chunk, , drop = FALSE], nb_dispersion)
      target_genes <- which(keep)[chunk]
      adjusted[target_genes, sample_idx] <- mapped$counts
      tail_events <- tail_events + mapped$tail_events
    }
    mapping[[source]] <- data.frame(
      source_batch = source, reference_batch = reference,
      genes_attempted = sum(keep), genes_adjusted = sum(keep),
      tail_events = tail_events,
      seconds = proc.time()[["elapsed"]] - started, stringsAsFactors = FALSE)
  }
  list(counts = adjusted, mapping = do.call(rbind, mapping))
}
