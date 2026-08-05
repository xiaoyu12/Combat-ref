estimate_batch_dispersion <- function(dge, design, batch, reference,
                                      preliminary, information) {
  genes <- rownames(dge$counts)
  batches <- levels(batch)
  estimated <- tryCatch(edgeR::estimateDisp(dge, design = design,
      robust = TRUE, trend.method = "locfit"), error = function(error)
        abort_combatrefql("combatrefql_statistical_error",
          sprintf("Batch-dispersion baseline estimation failed: %s",
            conditionMessage(error)), stage = "dispersion"))
    baseline <- estimated$tagwise.dispersion
    fallback <- estimated$trended.dispersion %||% estimated$common.dispersion
    if (length(fallback) == 1L) fallback <- rep(fallback, length(baseline))
    baseline[!is.finite(baseline) | baseline <= 0] <-
      fallback[!is.finite(baseline) | baseline <= 0]
    baseline <- pmin(pmax(as.numeric(baseline), 1e-8), 4)
    names(baseline) <- genes
    total_precision <- colSums(information$precision)
    if (!is.finite(total_precision[reference]) || total_precision[reference] <= 0)
      abort_combatrefql("combatrefql_statistical_error",
        "The reference batch has zero effective dispersion information.",
        "*" = "Choose a reference batch with replicated residual information.",
        stage = "dispersion", reference_batch = reference)

    raw <- stats::setNames(numeric(length(batches)), batches)
    status <- stats::setNames(rep("estimated", length(batches)), batches)
    convergence <- stats::setNames(rep(TRUE, length(batches)), batches)
    iterations <- stats::setNames(integer(length(batches)), batches)
    for (b in batches) {
      if (!is.finite(total_precision[b]) || total_precision[b] <= 0) {
        raw[b] <- NA_real_
        status[b] <- "not_estimable"
        next
      }
      idx <- batch == b
      objective <- function(log_multiplier) {
        phi <- baseline * exp(log_multiplier)
        ll <- rowMeans(stats::dnbinom(dge$counts[, idx, drop = FALSE],
          mu = preliminary$fitted.values[, idx, drop = FALSE],
          size = 1 / phi, log = TRUE))
        -sum(information$precision[, b] * ll, na.rm = TRUE)
      }
      opt <- stats::optim(0, objective, method = "Brent",
        lower = -log(10), upper = log(10))
      raw[b] <- opt$par
      convergence[b] <- opt$convergence == 0L
      iterations[b] <- unname(opt$counts[["function"]])
      if (!convergence[b]) status[b] <- "not_converged"
    }
    raw[status == "not_estimable"] <- raw[reference]
    shifts <- raw - raw[reference]
    shifts[reference] <- 0
    shifts[status == "not_estimable"] <- 0
    max_multiplier <- max(exp(shifts))
  baseline <- pmax(pmin(baseline * exp(raw[reference]),
    4 / max_multiplier), 1e-8 / min(exp(shifts)))

  multipliers <- exp(shifts)
  phi_reference <- baseline
  dispersion <- outer(phi_reference, multipliers[as.character(batch)], `*`)
  dimnames(dispersion) <- dimnames(dge$counts)
  rows <- lapply(batches, function(b) data.frame(
    gene = genes, batch = b, reference_batch = reference,
    baseline_dispersion = baseline,
    reference_dispersion = phi_reference,
    source_dispersion = baseline * multipliers[b],
    d_batch = unname(shifts[b]), multiplier = unname(multipliers[b]),
    design_residual_information = information$design_only$
      design_residual_information[match(b, information$design_only$batch)],
    weighted_residual_information = information$weighted[, b],
    precision_weight = information$precision[, b],
    precision_formula = information$formula,
    estimable = status[b] != "not_estimable", status = unname(status[b]),
    converged = unname(convergence[b]), iterations = unname(iterations[b]),
    stringsAsFactors = FALSE))
  list(dispersion = dispersion, reference = phi_reference,
    shifts = shifts, multipliers = multipliers,
    diagnostics = do.call(rbind, rows))
}
