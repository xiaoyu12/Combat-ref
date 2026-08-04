build_combatrefql_result <- function(prepared, dge, design, selected, fit,
                                     transported, keep, timings, call,
                                     normalization, chunk_size) {
  gene_status <- make_gene_status(prepared$counts, prepared$batch, keep,
    transported$mapping_valid)
  mapping_warnings <- if (nrow(transported$failures)) data.frame(
    stage = "mapping", class = "combatrefql_mapping_error",
    message = sprintf("NB transport failed for %s genes; complete original rows were retained.",
      sum(gene_status$status == "failed")), stringsAsFactors = FALSE) else
    empty_warnings()
  specification <- data.frame(
    method = "QL-fitted counterfactual means with negative-binomial mid-P quantile transport",
    normalization = normalization, transport = "negative-binomial mid-P",
    preserve_zeros = TRUE, reference_batch = selected$reference,
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
    nb_dispersion = as.numeric(fit$dispersion),
    average_ql_dispersion = as.numeric(fit$average.ql.dispersion),
    median_posterior_ql_dispersion = stats::median(fit$s2.post),
    stringsAsFactors = FALSE)
  model <- list(
    coefficients = fit$coefficients, nb_dispersion = as.numeric(fit$dispersion),
    ql_dispersion = fit$s2.prior, posterior_ql_dispersion = fit$s2.post,
    design = data.frame(samples = nrow(design$matrix),
      coefficients = ncol(design$matrix), rank = design$rank,
      residual_df = design$residual_df,
      condition_number = design$condition_number),
    design_matrix = design$matrix)
  diagnostics <- build_diagnostics(selected, transported, timings,
    prepared$input_actions, mapping_warnings, gene_status, ql)
  list(
    result = CombatRefQLFit(counts = transported$counts,
      specification = specification, samples = samples,
      gene_status = gene_status, diagnostics = diagnostics, call = call,
      model = model),
    warnings = mapping_warnings
  )
}
