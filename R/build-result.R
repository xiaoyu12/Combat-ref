build_combatrefql_result <- function(prepared, dge, design, selected, fit,
                                     shrinkage, dispersion, transported, keep,
                                     timings, call, normalization,
                                     chunk_size) {
  gene_status <- make_gene_status(prepared$counts, prepared$batch, keep,
    transported$mapping_valid)
  mapping_warnings <- if (nrow(transported$failures)) data.frame(
    stage = "mapping", class = "combatrefql_mapping_error",
    message = sprintf("NB transport failed for %s genes; complete original rows were retained.",
      sum(gene_status$status == "failed")), stringsAsFactors = FALSE) else
    empty_warnings()
  specification <- data.frame(
    method = "ComBat-refQL integrated reference-batch adjustment",
    normalization = normalization, transport = "negative-binomial mid-P",
    zero_handling = "model-consistent mid-P", reference_batch = selected$reference,
    reference_selection = if (is.null(prepared$reference)) "automatic" else "explicit",
    chunk_size = chunk_size, stringsAsFactors = FALSE)
  samples <- data.frame(
    sample = colnames(prepared$counts), batch = as.character(prepared$batch),
    group = if (is.null(prepared$group)) NA_character_ else as.character(prepared$group),
    raw_library_size = colSums(prepared$counts),
    normalization_factor = dge$samples$norm.factors,
    effective_library_size = dge$samples$lib.size * dge$samples$norm.factors,
    stringsAsFactors = FALSE)
  ql <- data.frame(
    edgeR_pipeline = "v4", robust = TRUE, legacy = FALSE,
    nb_dispersion_min = min(dispersion$dispersion),
    nb_dispersion_max = max(dispersion$dispersion),
    average_ql_dispersion = as.numeric(fit$average.ql.dispersion),
    median_posterior_ql_dispersion = stats::median(fit$s2.post),
    stringsAsFactors = FALSE)
  model <- list(
    coefficients = fit$unshrunk.coefficients %||% fit$coefficients,
    nb_dispersion = dispersion$reference,
    reference_dispersion = dispersion$reference,
    batch_dispersion_shifts = dispersion$shifts,
    batch_dispersion_multipliers = dispersion$multipliers,
    ql_dispersion = fit$s2.prior, posterior_ql_dispersion = fit$s2.post,
    design = data.frame(samples = nrow(design$matrix),
      coefficients = ncol(design$matrix), rank = design$rank,
      residual_df = design$residual_df,
      condition_number = design$condition_number),
    design_matrix = design$matrix)
  diagnostics <- build_diagnostics(selected, transported, timings,
    prepared$input_actions, mapping_warnings, gene_status, ql, shrinkage,
    dispersion)
  list(
    result = CombatRefQLFit(counts = transported$counts,
      specification = specification, samples = samples,
      gene_status = gene_status, diagnostics = diagnostics, call = call,
      model = model),
    warnings = mapping_warnings
  )
}
