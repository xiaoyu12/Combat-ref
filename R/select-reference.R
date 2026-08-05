build_local_reference_design <- function(n, group = NULL, covariates = NULL) {
  data <- data.frame(row.names = seq_len(n))
  candidates <- character()
  dropped <- character()
  missing <- character()

  add_candidate <- function(name, value) {
    if (is.character(value)) value <- factor(value)
    if (is.factor(value)) {
      absent <- setdiff(levels(value), unique(as.character(value)))
      if (length(absent)) missing <<- c(missing, paste0(name, ":", absent))
      value <- droplevels(value)
    }
    if (length(unique(value)) < 2L) {
      dropped <<- c(dropped, name)
      return()
    }
    data[[name]] <<- value
    candidates <<- c(candidates, name)
  }

  if (!is.null(group)) add_candidate("group", group)
  if (!is.null(covariates)) {
    names(covariates) <- make.names(names(covariates), unique = TRUE)
    for (name in names(covariates)) add_candidate(name, covariates[[name]])
  }

  accepted <- character()
  design <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
  for (term in candidates) {
    proposed <- stats::model.matrix(stats::reformulate(c(accepted, term)), data)
    rank <- qr(proposed)$rank
    if (rank < ncol(proposed) || n - rank <= 0L) {
      dropped <- c(dropped, term)
    } else {
      accepted <- c(accepted, term)
      design <- proposed
    }
  }

  fallback <- !length(accepted)
  list(
    matrix = design,
    formula = if (length(accepted))
      paste("~ 1 +", paste(accepted, collapse = " + ")) else "~ 1",
    dropped_columns = paste(dropped, collapse = ","),
    missing_levels = paste(missing, collapse = ","),
    fallback = fallback,
    fallback_reason = if (fallback) {
      if (length(dropped)) "no biological term was locally estimable" else "no biological terms supplied"
    } else NA_character_,
    residual_df = n - qr(design)$rank
  )
}

select_reference <- function(dge, batch, group = NULL, covariates = NULL,
                             explicit = NULL, progress = function(...) NULL) {
  if (!is.null(explicit)) {
    candidates <- levels(batch)
    return(list(
      reference = explicit,
      batches = data.frame(
        batch = candidates,
        samples = as.integer(table(batch)[candidates]),
        reference_score = NA_real_,
        selected = candidates == explicit,
        selection_method = "explicit",
        local_formula = NA_character_,
        local_columns = NA_character_,
        dropped_columns = NA_character_,
        missing_levels = NA_character_,
        fallback = NA,
        fallback_reason = NA_character_,
        residual_df = NA_integer_,
        effective_residual_information = NA_real_,
        second_best = FALSE,
        score_margin = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }
  rows <- lapply(levels(batch), function(candidate) {
    progress(candidate)
    idx <- batch == candidate
    local_design <- build_local_reference_design(sum(idx),
      if (is.null(group)) NULL else group[idx],
      if (is.null(covariates)) NULL else covariates[idx, , drop = FALSE]
    )
    local <- dge[, idx, keep.lib.sizes = TRUE]
    keep <- rowSums(local$counts) > 0
    fit <- if (sum(keep) < 2L) NULL else try(
      edgeR::glmQLFit(local[keep, , keep.lib.sizes = TRUE], local_design$matrix,
        dispersion = NULL, robust = TRUE, legacy = FALSE), silent = TRUE)
    score <- if (is.null(fit) || inherits(fit, "try-error") ||
      length(fit$dispersion) != 1L || !is.finite(fit$dispersion)) Inf else
      as.numeric(fit$dispersion)
    data.frame(
      batch = candidate, samples = sum(idx), reference_score = score,
      local_formula = local_design$formula,
      local_columns = paste(colnames(local_design$matrix), collapse = ","),
      dropped_columns = local_design$dropped_columns,
      missing_levels = local_design$missing_levels,
      fallback = local_design$fallback,
      fallback_reason = local_design$fallback_reason,
      residual_df = local_design$residual_df,
      effective_residual_information = local_design$residual_df,
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, rows)
  if (all(!is.finite(data$reference_score))) {
    abort_combatrefql("combatrefql_reference_error",
      "No batch could be scored as a reference.",
      "*" = "Ensure each batch has replicated, non-zero raw counts.", stage = "reference")
  }
  selected <- data$batch[which.min(data$reference_score)]
  ordered <- order(data$reference_score)
  second <- if (length(ordered) > 1L) data$batch[ordered[2L]] else NA_character_
  margin <- if (length(ordered) > 1L)
    data$reference_score[ordered[2L]] - data$reference_score[ordered[1L]] else
    NA_real_
  data$selected <- data$batch == selected
  data$second_best <- data$batch == second
  data$score_margin <- margin
  data$selection_method <-
    "automatic_biologically_adjusted_minimum_nb_dispersion"
  data <- data[c("batch", "samples", "reference_score", "selected",
    "selection_method", "local_formula", "local_columns", "dropped_columns",
    "missing_levels", "fallback", "fallback_reason", "residual_df",
    "effective_residual_information", "second_best", "score_margin")]
  list(reference = selected, batches = data)
}
