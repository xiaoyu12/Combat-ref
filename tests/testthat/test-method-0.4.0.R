test_that("normal EB weights respond to QL uncertainty and boundary pooling", {
  precise <- combatrefql:::fit_normal_eb(c(-1, 0, 1), rep(1e-10, 3))
  noisy <- combatrefql:::fit_normal_eb(c(-1, 0, 1), rep(1e6, 3))
  pooled <- combatrefql:::fit_normal_eb(c(2, 2, 2), c(0.1, 0.2, 0.3))

  expect_equal(precise$posterior, c(-1, 0, 1), tolerance = 1e-5)
  expect_equal(noisy$posterior, rep(noisy$mean, 3), tolerance = 1e-5)
  expect_identical(pooled$variance, 0)
  expect_equal(pooled$posterior, rep(pooled$mean, 3))
})

test_that("batch contrasts do not depend on factor level order", {
  z <- combatrefql_test_fixtures()[[2]]
  a <- combat_ref_ql(z$counts, factor(z$batch, c("batch_A", "batch_B")),
    z$group, reference = "batch_A", verbosity = "quiet")
  b <- combat_ref_ql(z$counts, factor(z$batch, c("batch_B", "batch_A")),
    z$group, reference = "batch_A", verbosity = "quiet")
  da <- a@diagnostics$batch_contrasts[order(a@diagnostics$batch_contrasts$gene), ]
  db <- b@diagnostics$batch_contrasts[order(b@diagnostics$batch_contrasts$gene), ]
  expect_equal(da$raw_delta, db$raw_delta, tolerance = 1e-8)
  expect_equal(da$ql_variance, db$ql_variance, tolerance = 1e-8)
})

test_that("QL contrast covariance agrees with the edgeR one-df test", {
  set.seed(42)
  genes <- 200L
  batch <- factor(rep(c("A", "B"), each = 6))
  group <- factor(rep(rep(c("C", "T"), each = 3), 2))
  design <- stats::model.matrix(~ 0 + batch + group)
  means <- outer(stats::runif(genes, 30, 200),
    exp(c(rep(0, 6), rep(0.15, 6))))
  counts <- matrix(stats::rnbinom(genes * 12, mu = means, size = 10), genes,
    dimnames = list(paste0("g", seq_len(genes)), NULL))
  dge <- edgeR::DGEList(counts)
  dge <- edgeR::normLibSizes(dge, method = "TMMwsp")
  dge <- edgeR::estimateDisp(dge, design)
  fit <- combatrefql:::fit_ql_model(dge, design)
  contrast <- combatrefql:::batch_contrast(design, "B", "A")
  covariance <- combatrefql:::contrast_covariance(fit, fit$dispersion, "B", "A")
  qlf <- edgeR::glmQLFTest(fit, contrast = contrast)$table$F
  wald <- covariance$raw_delta^2 / covariance$ql_variance

  expect_gt(stats::cor(qlf, wald), 0.999)
  expect_lt(stats::quantile(abs(qlf - wald) / pmax(qlf, 1e-8), 0.9), 0.03)
})

test_that("batch dispersion is reference-centred and information weighted", {
  z <- combatrefql_test_fixtures()[[1]]
  fit <- combat_ref_ql(z$counts, z$batch, z$group, reference = "1",
    mean_shrinkage = "none", verbosity = "quiet")
  dispersion <- fit@diagnostics$dispersion
  by_batch <- dispersion[!duplicated(dispersion$batch), ]

  expect_identical(by_batch$d_batch[by_batch$batch == "1"], 0)
  expect_equal(dispersion$source_dispersion / dispersion$reference_dispersion,
    exp(dispersion$d_batch), tolerance = 1e-12)
  expect_true(all(dispersion$weighted_residual_information >= 0))
  expect_true(all(dispersion$precision_weight >= 0))
})

test_that("effective information has the required residual total and saturation", {
  batch <- factor(c("A", "A", "B", "B"))
  design <- stats::model.matrix(~ 0 + batch)
  fitted <- matrix(c(10, 20, 30, 40), 1)
  info <- combatrefql:::effective_residual_information(design, batch, fitted,
    matrix(0.1, 1, 4))
  expect_equal(sum(info$design_only$design_residual_information),
    nrow(design) - qr(design)$rank, tolerance = 1e-12)
  expect_equal(sum(info$weighted[1, ]), nrow(design) - qr(design)$rank,
    tolerance = 1e-10)
  expect_gt(info$precision[1, "B"], 0)

  saturated <- combatrefql:::effective_residual_information(diag(4), batch,
    fitted, matrix(0.1, 1, 4))
  expect_equal(saturated$weighted, matrix(0, 1, 2,
    dimnames = list(NULL, c("A", "B"))))
  expect_equal(saturated$precision, matrix(0, 1, 2,
    dimnames = list(NULL, c("A", "B"))))
  expect_lt(1 / trigamma(0.5), 1 / trigamma(2))
})

test_that("zero batch information has explicit source and reference behavior", {
  z <- combatrefql_test_fixtures()[[1]]
  prepared <- combatrefql:::prepare_inputs(z$counts, z$batch, z$group, NULL,
    "1", "round")
  dge <- combatrefql:::normalize_counts(prepared, "TMMwsp")
  keep <- combatrefql:::assess_support(prepared$counts, prepared$batch)
  dge <- dge[keep, , keep.lib.sizes = TRUE]
  design <- combatrefql:::build_design(prepared$batch, prepared$group)$matrix
  preliminary <- combatrefql:::fit_ql_model(dge, design)
  information <- combatrefql:::effective_residual_information(design,
    prepared$batch, preliminary$fitted.values,
    matrix(preliminary$dispersion, nrow(dge), ncol(dge)))

  source_zero <- information
  source_zero$precision[, "2"] <- 0
  estimated <- combatrefql:::estimate_batch_dispersion(dge, design,
    prepared$batch, "1", preliminary, source_zero, "batch_scaled")
  expect_identical(unname(estimated$shifts["2"]), 0)
  expect_true(all(estimated$diagnostics$status[estimated$diagnostics$batch == "2"] ==
    "not_estimable"))

  reference_zero <- information
  reference_zero$precision[, "1"] <- 0
  expect_error(combatrefql:::estimate_batch_dispersion(dge, design,
    prepared$batch, "1", preliminary, reference_zero, "batch_scaled"),
    class = "combatrefql_statistical_error")
})

test_that("mid-P zero policy and distinct dispersions are explicit", {
  y <- matrix(0, 1, 1)
  mu <- matrix(4, 1, 1)
  phi <- matrix(0.2, 1, 1)
  expected <- 0.5 * stats::dnbinom(0, mu = 4, size = 5)
  mapped <- combatrefql:::midp_transport(y, mu, mu, phi, phi, "model")
  expect_equal(mapped$probability, matrix(expected, 1, 1))
  expect_identical(combatrefql:::midp_transport(y, mu, mu, phi, phi,
    "preserve")$counts, matrix(0, 1, 1))

  support <- matrix(0:100, 1)
  equal <- combatrefql:::midp_transport(support, matrix(20, 1, 101),
    matrix(20, 1, 101), matrix(0.2, 1, 101), matrix(0.2, 1, 101))
  changed <- combatrefql:::midp_transport(support, matrix(20, 1, 101),
    matrix(20, 1, 101), matrix(0.5, 1, 101), matrix(0.05, 1, 101))
  expect_identical(equal$counts, support * 1)
  expect_false(identical(changed$counts, equal$counts))
})

test_that("all formal ablations run and label their diagnostics", {
  z <- combatrefql_test_fixtures()[[1]]
  configurations <- expand.grid(mean = c("none", "ql_eb"),
    dispersion = c("shared", "batch_scaled"), stringsAsFactors = FALSE)
  fits <- Map(function(mean, dispersion) combat_ref_ql(z$counts, z$batch,
    z$group, reference = "1", mean_shrinkage = mean,
    dispersion_model = dispersion, verbosity = "quiet"),
    configurations$mean, configurations$dispersion)

  expect_equal(unname(vapply(fits, function(x) x@specification$mean_shrinkage,
    character(1))), configurations$mean)
  expect_equal(unname(vapply(fits, function(x) x@specification$dispersion_model,
    character(1))), configurations$dispersion)
  expect_true(all(vapply(fits, function(x)
    identical(x@counts[, z$batch == 1], z$counts[, z$batch == 1]), logical(1))))
  expect_identical(fits[[4]]@counts,
    combat_ref_ql(z$counts, z$batch, z$group, reference = "1",
      verbosity = "quiet")@counts)
})
