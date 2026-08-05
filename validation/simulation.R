args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[[1]], c("fast", "full")) else "fast"
devtools::load_all(".", quiet = TRUE)

simulate_counts <- function(seed, mean_pattern, dispersion_ratio, replicates,
                            balance, partial_overlap, weak_reference,
                            genes = 120L) {
  set.seed(seed)
  batch <- factor(rep(c("A", "B"), each = 2 * replicates))
  group <- factor(rep(rep(c("control", "treated"), each = replicates), 2))
  if (balance == "unbalanced") group[batch == "B"][seq_len(replicates %/% 2)] <- "treated"
  if (partial_overlap) group[batch == "B" & group == "control"][1] <- "treated"
  library_factor <- exp(seq(log(0.45), log(2.2), length.out = length(batch)))
  base <- exp(stats::runif(genes, log(2), log(300)))
  biological <- c(rep(log(2), genes %/% 4), rep(0, genes - genes %/% 4))
  delta <- switch(mean_pattern,
    none = rep(0, genes),
    sparse = c(rep(log(1.8), genes %/% 10), rep(0, genes - genes %/% 10)),
    widespread = stats::rnorm(genes, log(1.5), 0.15))
  module <- stats::rnorm(length(batch), 0, 0.15)
  eta <- outer(log(base), rep(1, length(batch))) +
    outer(biological, as.numeric(group == "treated")) +
    outer(delta, as.numeric(batch == "B")) +
    outer(rep(c(-1, 1), length.out = genes), module) +
    matrix(log(library_factor), genes, length(batch), byrow = TRUE)
  phi_a <- 0.08 + 0.35 / sqrt(base + 1)
  if (weak_reference) phi_a <- phi_a * 1.8
  phi <- outer(phi_a, ifelse(batch == "A", 1, dispersion_ratio))
  counts <- matrix(stats::rnbinom(length(eta), mu = exp(eta), size = 1 / phi),
    genes, dimnames = list(paste0("g", seq_len(genes)),
      paste0("s", seq_along(batch))))
  counts[1, 1] <- counts[1, 1] * 20
  list(counts = counts, batch = batch, group = group, delta = delta,
    d = log(dispersion_ratio), biological = biological,
    reference = if (weak_reference) "A" else "A")
}

run_ql <- function(x, mean_shrinkage, dispersion_model) {
  fit <- combat_ref_ql(x$counts, x$batch, x$group, reference = x$reference,
    mean_shrinkage = mean_shrinkage, dispersion_model = dispersion_model,
    verbosity = "quiet")
  contrast <- fit@diagnostics$batch_contrasts
  dispersion <- fit@diagnostics$dispersion
  list(counts = fit@counts, delta = contrast$posterior_delta,
    d = unique(dispersion$d_batch[dispersion$batch == "B"]),
    invariant = identical(fit@counts[, x$batch == x$reference, drop = FALSE],
      x$counts[, x$batch == x$reference, drop = FALSE]))
}

load_external <- function() {
  suppressPackageStartupMessages(library(edgeR))
  directory <- tempfile("combat-ref-")
  dir.create(directory)
  load_pinned <- function(root, main) {
    files <- c("helper_seq.R", main)
    paths <- file.path(directory, paste0(basename(root), "-", files))
    for (i in seq_along(files))
      utils::download.file(paste0(root, files[[i]]), paths[[i]], quiet = TRUE)
    environment <- new.env(parent = globalenv())
    sys.source(paths[[1]], environment)
    sys.source(paths[[2]], environment)
    environment
  }
  combat_ref <- load_pinned(
    "https://raw.githubusercontent.com/xiaoyu12/Combat-ref/f0c6d3f313e9f1f29de2c8c8f12b98113443be55/",
    "Combat_ref.R")
  combat_seq <- load_pinned(
    "https://raw.githubusercontent.com/zhangyuqing/ComBat-seq/2673dd8f23052ea4354ffa5e7600fabcd3f08269/",
    "ComBat_seq.R")
  list(combat_ref = combat_ref$ComBat_ref,
    combat_seq = combat_seq$ComBat_seq,
    combat_ref_commit = "f0c6d3f313e9f1f29de2c8c8f12b98113443be55",
    combat_seq_commit = "2673dd8f23052ea4354ffa5e7600fabcd3f08269")
}

alignment_metrics <- function(adjusted, x) {
  group_batch <- function(b, g) rowMeans(adjusted[, x$batch == b & x$group == g,
    drop = FALSE])
  residual <- rowMeans(cbind(
    log1p(group_batch("B", "control")) - log1p(group_batch("A", "control")),
    log1p(group_batch("B", "treated")) - log1p(group_batch("A", "treated"))))
  lfc <- log1p(rowMeans(adjusted[, x$group == "treated", drop = FALSE])) -
    log1p(rowMeans(adjusted[, x$group == "control", drop = FALSE]))
  variance_ratio <- stats::var(as.numeric(adjusted[, x$batch == "B"])) /
    stats::var(as.numeric(adjusted[, x$batch == "A"]))
  c(mean_alignment_rmse = sqrt(mean(residual^2)),
    variance_alignment = abs(log(variance_ratio)),
    biological_lfc_rmse = sqrt(mean((lfc - x$biological)^2)),
    transformation = mean(abs(adjusted - x$counts)))
}

de_metrics <- function(x) {
  design <- stats::model.matrix(~ x$batch + x$group)
  dge <- edgeR::DGEList(x$counts)
  dge <- edgeR::normLibSizes(dge, method = "TMMwsp")
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE, legacy = FALSE)
  table <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = ncol(design)), n = Inf,
    sort.by = "none")$table
  selected <- table$FDR <= 0.05
  truth <- x$biological != 0
  c(de_fdr = if (any(selected)) mean(!truth[selected]) else 0,
    de_power = mean(selected[truth]))
}

configurations <- data.frame(method = c("shared_none", "shared_ql_eb",
  "batch_scaled_none", "batch_scaled_ql_eb"),
  mean = c("none", "ql_eb", "none", "ql_eb"),
  dispersion = c("shared", "shared", "batch_scaled", "batch_scaled"))
external <- if (mode == "full") load_external() else NULL
grid <- if (mode == "fast") data.frame(
  seed = c(104729L, 130363L), mean_pattern = c("none", "widespread"),
  dispersion_ratio = c(1, 2), replicates = c(3L, 5L),
  balance = c("balanced", "unbalanced"), partial_overlap = c(FALSE, TRUE),
  weak_reference = c(FALSE, TRUE)) else expand.grid(
  seed = c(104729L, 130363L),
  mean_pattern = c("none", "sparse", "widespread"),
  dispersion_ratio = c(0.5, 1, 2), replicates = c(2L, 3L, 5L, 10L),
  balance = c("balanced", "unbalanced"), stringsAsFactors = FALSE)
if (mode == "full") {
  grid$partial_overlap <- rep(c(FALSE, TRUE), length.out = nrow(grid))
  grid$weak_reference <- rep(c(FALSE, FALSE, TRUE), length.out = nrow(grid))
}

rows <- list()
k <- 0L
for (i in seq_len(nrow(grid))) {
  x <- do.call(simulate_counts, as.list(grid[i, ]))
  de <- de_metrics(x)
  methods <- c(configurations$method, "raw_edgeR_design",
    if (mode == "full") c("ComBat-ref", "ComBat-seq"))
  for (method in methods) {
    k <- k + 1L
    started <- proc.time()[["elapsed"]]
    if (method %in% configurations$method) {
      j <- match(method, configurations$method)
      result <- tryCatch(run_ql(x, configurations$mean[j],
        configurations$dispersion[j]), error = identity)
    } else if (method == "raw_edgeR_design") {
      result <- list(counts = x$counts, delta = rep(NA_real_, length(x$delta)),
        d = NA_real_, invariant = TRUE)
    } else {
      output <- tryCatch({
        invisible(utils::capture.output(adjusted <- if (method == "ComBat-ref")
          external$combat_ref(x$counts, x$batch, x$group) else
          external$combat_seq(x$counts, x$batch, x$group)))
        adjusted
      }, error = identity)
      result <- if (inherits(output, "error")) output else list(counts = output,
        delta = rep(NA_real_, length(x$delta)), d = NA_real_,
        invariant = if (method == "ComBat-ref")
          any(vapply(levels(x$batch), function(b) identical(
            output[, x$batch == b, drop = FALSE],
            x$counts[, x$batch == b, drop = FALSE]), logical(1))) else NA)
    }
    if (inherits(result, "error")) {
      rows[[k]] <- cbind(grid[i, ], method = method,
        status = "error", message = conditionMessage(result))
      next
    }
    metrics <- alignment_metrics(result$counts, x)
    rows[[k]] <- data.frame(grid[i, ], method = method, status = "ok",
      delta_rmse = sqrt(mean((result$delta - x$delta)^2, na.rm = TRUE)),
      delta_bias = mean(result$delta - x$delta, na.rm = TRUE),
      d_rmse = abs(result$d - x$d), d_bias = result$d - x$d,
      dispersion_ratio_estimate = exp(result$d), as.list(metrics),
      reference_invariant = result$invariant, de_fdr = de[["de_fdr"]],
      de_power = de[["de_power"]],
      seconds = proc.time()[["elapsed"]] - started)
  }
}
results <- do.call(rbind, rows)
for (column in names(results)[vapply(results, is.numeric, logical(1))])
  results[[column]][is.nan(results[[column]])] <- NA_real_
path <- sprintf("validation/results-%s.tsv", mode)
utils::write.table(results, path, sep = "\t", quote = FALSE, row.names = FALSE,
  na = "NA")
cat(sprintf("%s: %s scenario-method rows; %s errors; reference invariant: %s\n",
  path, nrow(results), sum(results$status != "ok"),
  all(results$reference_invariant[results$status == "ok"])))
