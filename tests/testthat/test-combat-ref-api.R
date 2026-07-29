test_that("combat_ref reproduces legacy fixtures", {
  fixtures <- combatrefql_legacy_fixtures()

  for (nm in names(fixtures)) {
    fixture <- fixtures[[nm]]
    fit <- combat_ref(
      counts = fixture$counts,
      batch = fixture$batch,
      group = fixture$group
    )

    expect_s3_class(fit, "combatref_result")
    expect_equal(fit$adjusted_counts, fixture$expected_adjusted)
    expect_equal(fit$reference_batch, fixture$expected_reference_batch)
    expect_equal(dim(fit$adjusted_counts), dim(fixture$counts))
    expect_equal(dimnames(fit$adjusted_counts), dimnames(fixture$counts))
    expect_true(all(fit$adjusted_counts >= 0))
    expect_true(all(abs(fit$adjusted_counts - round(fit$adjusted_counts)) < 1e-8))
  }
})

test_that("combat_ref is deterministic", {
  fixture <- combatrefql_legacy_fixtures()[["numeric_balanced"]]
  first <- combat_ref(fixture$counts, fixture$batch, fixture$group)
  second <- combat_ref(fixture$counts, fixture$batch, fixture$group)

  expect_equal(second$adjusted_counts, first$adjusted_counts)
  expect_equal(second$reference_batch, first$reference_batch)
})

test_that("combat_ref is silent by default in noninteractive calls and verbose on request", {
  fixture <- combatrefql_legacy_fixtures()[["numeric_balanced"]]

  expect_silent(combat_ref(fixture$counts, fixture$batch, fixture$group, verbose = FALSE))
  expect_message(
    combat_ref(fixture$counts, fixture$batch, fixture$group, verbose = TRUE),
    "Validating input data"
  )
})

test_that("keep_model controls model details", {
  fixture <- combatrefql_legacy_fixtures()[["numeric_balanced"]]

  compact <- combat_ref(fixture$counts, fixture$batch, fixture$group)
  expect_false("model" %in% names(compact))
  expect_false("dispersions" %in% names(compact))
  expect_false("design" %in% names(compact))

  full <- combat_ref(fixture$counts, fixture$batch, fixture$group, keep_model = TRUE)
  expect_true("model" %in% names(full))
  expect_true(all(c("design", "common_dispersions", "coefficients", "fitted_means", "library_sizes") %in% names(full$model)))
})

test_that("explicit reference batch is validated and used", {
  fixture <- combatrefql_legacy_fixtures()[["character_unbalanced"]]

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
  fixture <- combatrefql_legacy_fixtures()[["numeric_balanced"]]
  fit <- combat_ref(fixture$counts, fixture$batch, fixture$group)

  expect_invisible(print(fit))
  smry <- summary(fit)
  expect_s3_class(smry, "summary.combatref_result")
  expect_invisible(print(smry))
})
