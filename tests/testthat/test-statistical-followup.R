test_that("within-batch designs always include an estimable intercept", {
  group <- factor(c("a", "b", "a", "b", "a", "a", "b", "b"))
  covariates <- data.frame(age = 1:8, site = factor(rep(c("x", "y"), 4)))
  for (idx in list(rep(TRUE, 8), c(TRUE, TRUE, TRUE, TRUE, rep(FALSE, 4)),
                   c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))) {
    x <- combatrefql:::build_within_batch_design(group, covariates, idx)
    expect_equal(unname(x$design[, "(Intercept)"]), rep(1, sum(idx)))
    expect_equal(qr(x$design)$rank, ncol(x$design))
    expect_gt(x$residual_df, 0)
  }
})

test_that("missing local factor levels are dropped and audited", {
  group <- factor(c("a", "a", "a", "b", "b", "b"))
  x <- combatrefql:::build_within_batch_design(group, NULL, seq_along(group) <= 3)
  expect_equal(colnames(x$design), "(Intercept)")
  expect_true(x$fallback)
  expect_match(x$reason, "absent local levels")
})

test_that("automatic reference diagnostics are finite and invariant to level order", {
  f <- combatrefql_test_fixtures()[["numeric_balanced"]]
  first <- combat_ref(f$counts, factor(f$batch, levels = c("1", "2")), f$group)
  second <- combat_ref(f$counts, factor(f$batch, levels = c("2", "1")), f$group)
  expect_true(all(is.finite(first$diagnostics$reference_scores)))
  expect_equal(first$reference_batch, second$reference_batch)
  expect_null(first$diagnostics$reference_selection_frequency)
  expect_true(is.na(first$diagnostics$reference_stable))
  expect_true(all(first$diagnostics$reference_local_designs$columns %in%
                    c("(Intercept),groupb", "(Intercept),grouptrt", "(Intercept),groupctrl")))
})

test_that("QL fit uses one global NB dispersion, never baseline tagwise dispersion", {
  f <- combatrefql_test_fixtures()[["numeric_balanced"]]
  fit <- combat_ref_ql(f$counts, f$batch, f$group, keep_model = TRUE)
  expect_length(fit$model$ql_nb_global_dispersion, 1L)
  expect_true(is.finite(fit$model$ql_nb_global_dispersion))
  expect_false(identical(fit$model$ql_nb_global_dispersion,
                         fit$model$baseline_nb_tagwise_dispersion))
  expect_length(fit$model$ql_dispersion, sum(fit$gene_status$adjusted))
})

test_that("mapping remains numeric above the 32-bit integer range", {
  y <- matrix(c(3e9, 4e9), nrow = 1)
  mu <- matrix(c(3e9, 4e9), nrow = 1)
  mapped <- combatrefql:::match_quantiles(y, mu, 0.01, mu, 0.01,
                                           keep_zero = FALSE, method = "midp")
  expect_type(mapped, "double")
  expect_true(all(is.finite(mapped)))
  expect_true(all(mapped > .Machine$integer.max))
  expect_true(all(mapped == floor(mapped)))
})
