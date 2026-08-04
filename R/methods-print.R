#' @export
S7::method(print, CombatRefQLFit) <- function(x, ...) {
  outcomes <- stats::setNames(x@diagnostics$outcomes$genes, x@diagnostics$outcomes$outcome)
  total <- x@diagnostics$timing$seconds[x@diagnostics$timing$stage == "total"]
  cli::cli_h1("ComBat-refQL fit"); cli::cli_h2("Data")
  cli::cli_text("  {format(nrow(x@counts), big.mark = ',')} genes x {ncol(x@counts)} samples")
  cli::cli_text("  {length(unique(x@samples$batch))} batches - {length(unique(stats::na.omit(x@samples$group)))} biological groups")
  cli::cli_h2("Model"); cli::cli_text("  QL GLM:         edgeR v4")
  cli::cli_text("  Normalization:  {x@specification$normalization}")
  cli::cli_text("  Transport:      negative-binomial mid-P")
  selection <- if (x@specification$reference_selection == "automatic") "automatically selected" else "explicit"
  cli::cli_text("  Reference:      {x@specification$reference_batch} - {selection}")
  cli::cli_h2("Result"); cli::cli_alert_success("{outcomes[['adjusted']] %||% 0} adjusted")
  cli::cli_alert_info("{outcomes[['unsupported']] %||% 0} unchanged")
  cli::cli_alert_success("{outcomes[['failed']] %||% 0} failed")
  if (nrow(x@diagnostics$warnings)) cli::cli_text("  {nrow(x@diagnostics$warnings)} stored warnings")
  cli::cli_text("Runtime: {format(round(total, 1), nsmall = 1)} seconds")
  cli::cli_text("Use summary(x) for diagnostics."); invisible(x)
}

#' @export
S7::method(print, CombatRefQLSummary) <- function(x, ...) {
  cli::cli_h1("ComBat-refQL summary"); cli::cli_h2("Design"); design <- x@design[1L, ]
  cli::cli_text("  Samples:          {design$samples}"); cli::cli_text("  Coefficients:     {design$coefficients}")
  cli::cli_text("  Rank:             {design$rank}"); cli::cli_text("  Residual df:      {design$residual_df}")
  cli::cli_text("  Condition number: {format(round(design$condition_number, 2), nsmall = 2)}")
  cli::cli_h2("Reference"); cli::cli_text("  batch  samples  score  selected")
  for (i in seq_len(nrow(x@batches))) cli::cli_text("  {x@batches$batch[i]}  {x@batches$samples[i]}  {format(signif(x@batches$reference_score[i], 4))}  {if (x@batches$selected[i]) 'yes' else 'no'}")
  cli::cli_h2("QL model"); ql <- x@ql[1L, ]
  cli::cli_text("  edgeR pipeline:          v4"); cli::cli_text("  Robust estimation:      yes"); cli::cli_text("  Legacy mode:            no")
  cli::cli_text("  Average QL dispersion:  {format(round(ql$average_ql_dispersion, 3), nsmall = 3)}")
  cli::cli_text("  Median posterior QL:    {format(round(ql$median_posterior_ql_dispersion, 3), nsmall = 3)}")
  cli::cli_h2("Gene outcomes")
  for (i in seq_len(nrow(x@gene_outcomes))) cli::cli_text("  {x@gene_outcomes$outcome[i]}: {x@gene_outcomes$genes[i]} ({format(round(100*x@gene_outcomes$proportion[i], 1), nsmall=1)}%)")
  actions <- x@input_actions[x@input_actions$lossy, , drop = FALSE]
  if (nrow(actions)) { cli::cli_h2("Input actions"); for (message in actions$message) cli::cli_alert_info(message) }
  cli::cli_h2("Timing"); for (i in seq_len(nrow(x@timings))) cli::cli_text("  {x@timings$stage[i]}: {format(round(x@timings$seconds[i], 2), nsmall=2)} s")
  invisible(x)
}
