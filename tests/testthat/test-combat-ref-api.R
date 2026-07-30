test_that("combat_ref is deterministic", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  first <- combat_ref(fixture$counts, fixture$batch, fixture$group)
  second <- combat_ref(fixture$counts, fixture$batch, fixture$group)

  expect_equal(second$adjusted_counts, first$adjusted_counts)
  expect_equal(second$reference_batch, first$reference_batch)
})

test_that("combat_ref is silent by default in noninteractive calls and verbose on request", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]

  expect_silent(combat_ref(fixture$counts, fixture$batch, fixture$group))
  expect_silent(combat_ref(fixture$counts, fixture$batch, fixture$group, verbose = FALSE))
  expect_message(
    combat_ref(fixture$counts, fixture$batch, fixture$group, verbose = TRUE),
    "Fitting baseline"
  )
})

test_that("README quick-start example runs", {
  counts <- matrix(
    c(10, 12, 30, 33, 8, 9, 28, 31,
      4, 5, 14, 16, 3, 4, 15, 17,
      50, 55, 120, 130, 48, 52, 125, 135,
      1, 2, 4, 5, 1, 2, 5, 6),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("gene", 1:4), paste0("sample", 1:8))
  )
  metadata <- data.frame(
    batch = rep(c("batch_A", "batch_B"), each = 2, times = 2),
    condition = rep(c("control", "treated"), each = 4),
    row.names = colnames(counts)
  )

  fit <- combat_ref(
    counts = counts,
    batch = metadata$batch,
    group = metadata$condition,
    verbose = FALSE
  )
  adjusted_counts <- fit$adjusted_counts

  expect_s3_class(fit, "combatref_result")
  expect_equal(dim(adjusted_counts), dim(counts))
  expect_equal(dimnames(adjusted_counts), dimnames(counts))
})

test_that("keep_model controls model details", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]

  compact <- combat_ref(fixture$counts, fixture$batch, fixture$group)
  expect_false("model" %in% names(compact))
  expect_false("dispersions" %in% names(compact))
  expect_false("design" %in% names(compact))

  full <- combat_ref(fixture$counts, fixture$batch, fixture$group, keep_model = TRUE)
  expect_true("model" %in% names(full))
  expect_true(all(c("design", "reference_design", "coefficients", "fitted_means",
                    "reference_means", "baseline_nb_tagwise_dispersion") %in% names(full$model)))
})

test_that("explicit reference batch is validated and used", {
  fixture <- combatrefql_test_fixtures()[["character_unbalanced"]]

  fit <- combat_ref(
    counts = fixture$counts,
    batch = fixture$batch,
    group = fixture$group,
    reference_batch = "batch_B"
  )
  expect_equal(fit$reference_batch, "batch_B")
  expect_equal(dimnames(fit$adjusted_counts), dimnames(fixture$counts))

  expect_error(
    combat_ref(fixture$counts, fixture$batch, fixture$group, reference_batch = "missing"),
    "reference_batch"
  )
})

test_that("result print and summary methods are compact", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  fit <- combat_ref(fixture$counts, fixture$batch, fixture$group)

  expect_invisible(print(fit))
  smry <- summary(fit)
  expect_s3_class(smry, "summary.combatref_result")
  expect_invisible(print(smry))
})
