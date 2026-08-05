download_cached <- function(url, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) utils::download.file(url, path, mode = "wb", quiet = TRUE)
  if (!file.exists(path) || file.info(path)$size == 0) stop("Download failed: ", url)
  path
}

load_comparators <- function() {
  suppressPackageStartupMessages(library(edgeR))
  load_method <- function(repository, commit, main) {
    environment <- new.env(parent = globalenv())
    root <- sprintf("https://raw.githubusercontent.com/%s/%s/", repository, commit)
    for (file in c("helper_seq.R", main)) {
      path <- download_cached(paste0(root, file), file.path("validation/.cache/methods",
        paste(repository, commit, file, sep = "-")))
      sys.source(path, environment)
    }
    environment
  }
  ref <- load_method("xiaoyu12/Combat-ref",
    "f0c6d3f313e9f1f29de2c8c8f12b98113443be55", "Combat_ref.R")
  seq <- load_method("zhangyuqing/ComBat-seq",
    "2673dd8f23052ea4354ffa5e7600fabcd3f08269", "ComBat_seq.R")
  list(`ComBat-ref` = ref$ComBat_ref, `ComBat-seq` = seq$ComBat_seq)
}

log_cpm <- function(counts) edgeR::cpm(counts, log = TRUE, prior.count = 1)

pca_scores <- function(counts, batch, group, method, dataset) {
  x <- t(log_cpm(counts))
  keep <- order(apply(x, 2, stats::var), decreasing = TRUE)[seq_len(min(500L, ncol(x)))]
  pca <- stats::prcomp(x[, keep, drop = FALSE], scale. = FALSE)
  data.frame(dataset = dataset, method = method, sample = rownames(x),
    batch = as.character(batch), group = as.character(group),
    PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)
}

conditional_batch_r2 <- function(counts, batch, group) {
  y <- t(log_cpm(counts))
  pcs <- stats::prcomp(y[, order(apply(y, 2, stats::var), decreasing = TRUE)[
    seq_len(min(500L, ncol(y)))], drop = FALSE])
  scores <- pcs$x[, seq_len(min(10L, ncol(pcs$x))), drop = FALSE]
  variance <- pcs$sdev[seq_len(ncol(scores))]^2
  incremental <- vapply(seq_len(ncol(scores)), function(i) {
    reduced <- stats::lm(scores[, i] ~ group)
    full <- stats::lm(scores[, i] ~ group + batch)
    max(0, summary(full)$adj.r.squared - summary(reduced)$adj.r.squared)
  }, numeric(1))
  stats::weighted.mean(incremental, variance)
}

alignment_metrics <- function(counts, batch, group) {
  y <- log_cpm(counts)
  design <- stats::model.matrix(~ group + batch)
  batch_columns <- grep("^batch", colnames(design))
  coefficients <- qr.solve(design, t(y))
  batch_effect <- if (length(batch_columns)) coefficients[batch_columns, , drop = FALSE] else
    matrix(0, 1, nrow(y))
  residuals <- t(y) - design %*% coefficients
  batch_variance <- vapply(levels(batch), function(b)
    mean(apply(residuals[batch == b, , drop = FALSE], 2, stats::var)), numeric(1))
  c(conditional_batch_r2 = conditional_batch_r2(counts, batch, group),
    mean_alignment_rmse = sqrt(mean(batch_effect^2)),
    dispersion_alignment = stats::sd(log(pmax(batch_variance, 1e-12))))
}

save_pca_figure <- function(scores, path, title) {
  plot <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, color = batch, shape = group)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85) +
    ggplot2::facet_wrap(~method, scales = "free") +
    ggplot2::labs(title = title, color = "Batch", shape = "Biology") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom")
  ggplot2::ggsave(path, plot, width = 8, height = 4.8, device = "svg")
}

run_validation <- function(counts, batch, group, dataset, output_dir) {
  batch <- droplevels(factor(batch)); group <- droplevels(factor(group))
  keep <- rowSums(counts) > 0 & vapply(seq_len(nrow(counts)), function(i)
    all(vapply(levels(batch), function(b) any(counts[i, batch == b] > 3), logical(1))),
    logical(1))
  counts <- unname(as.matrix(counts[keep, , drop = FALSE]))
  storage.mode(counts) <- "double"
  dimnames(counts) <- list(paste0("g", seq_len(nrow(counts))), paste0("s", seq_len(ncol(counts))))
  fit <- combat_ref_ql(counts, batch, group, verbosity = "quiet")
  adjusted <- list(Unadjusted = counts, `ComBat-refQL` = fit@counts)
  for (method in names(comparators <- load_comparators())) {
    result <- tryCatch({
      invisible(utils::capture.output(value <- suppressMessages(
        comparators[[method]](counts, batch, group))))
      value
    }, error = identity)
    if (!inherits(result, "error")) adjusted[[method]] <- result
  }
  metrics <- do.call(rbind, lapply(names(adjusted), function(method)
    data.frame(dataset = dataset, method = method,
      as.list(alignment_metrics(adjusted[[method]], batch, group)),
      reference_invariant = if (method == "ComBat-refQL") identical(
        adjusted[[method]][, batch == fit@specification$reference_batch, drop = FALSE],
        counts[, batch == fit@specification$reference_batch, drop = FALSE]) else NA,
      genes = nrow(counts), samples = ncol(counts),
      reference = if (method == "ComBat-refQL") fit@specification$reference_batch else NA,
      stringsAsFactors = FALSE)))
  scores <- do.call(rbind, lapply(names(adjusted), function(method)
    pca_scores(adjusted[[method]], batch, group, method, dataset)))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(metrics, file.path(output_dir, "metrics.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)
  utils::write.table(scores, file.path(output_dir, "pca-scores.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)
  save_pca_figure(scores, file.path(output_dir, "pca.svg"),
    paste(dataset, "reference-batch validation"))
  invisible(metrics)
}
