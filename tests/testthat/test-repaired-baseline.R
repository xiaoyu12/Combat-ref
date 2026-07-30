test_that("TMMwsp is default and offsets are fixed before filtering", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  fit <- combat_ref(fixture$counts, fixture$batch, fixture$group, keep_model = TRUE)
  direct <- edgeR::normLibSizes(edgeR::DGEList(fixture$counts), method = "TMMwsp")
  expect_equal(fit$normalization, "TMMwsp")
  expect_equal(unname(fit$diagnostics$normalization_factors), unname(direct$samples$norm.factors))
  expect_equal(fit$diagnostics$effective_library_sizes,
               direct$samples$lib.size * direct$samples$norm.factors)
})

test_that("group NULL retains one baseline coefficient per batch", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  fit <- combat_ref(fixture$counts, fixture$batch, group = NULL, keep_model = TRUE)
  expect_equal(colnames(fit$model$design), paste0("batch", levels(factor(fixture$batch))))
  expect_equal(qr(fit$model$design)$rank, ncol(fit$model$design))
  expect_true(all(is.finite(fit$model$fitted_means)))
  expect_true(all(fit$model$fitted_means > 0))
})

test_that("counterfactual prediction preserves reference samples", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  fit <- combat_ref(fixture$counts, fixture$batch, fixture$group,
                    reference_batch = "1", keep_model = TRUE)
  ref <- factor(fixture$batch) == "1"
  expect_equal(fit$adjusted_counts[, ref], fixture$counts[, ref])
  expect_true(all(fit$model$reference_design[, "batch1"] == 1))
})

test_that("mid-P mapping is finite monotone and identity-safe", {
  y <- matrix(0:100, nrow = 1)
  mu <- matrix(20, 1, length(y))
  mapped <- combatrefql:::match_quantiles(y, mu, 0.2, mu, 0.2, FALSE, "midp")
  expect_true(all(is.finite(mapped)))
  expect_true(all(diff(as.numeric(mapped)) >= 0))
  expect_lte(max(abs(mapped - y)), 1)
})

test_that("QL support keeps NB and QL dispersion separate", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  fit <- combat_ref_ql(fixture$counts, fixture$batch, fixture$group, keep_model = TRUE)
  expect_equal(fit$fit_method, "ql_assisted")
  expect_true(!is.null(fit$model$nb_dispersion))
  expect_true(!is.null(fit$model$ql_dispersion))
  expect_false(identical(fit$model$nb_dispersion, fit$model$ql_dispersion))
})

test_that("unsupported zero-batch genes are visible", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  fit <- combat_ref(fixture$counts, fixture$batch, fixture$group)
  expect_true(all(fit$gene_status$status[fit$gene_status$gene %in% fit$removed_genes] ==
                    "unsupported_zero_batch"))
  expect_equal(dim(fit$adjusted_counts), dim(fixture$counts))
})
