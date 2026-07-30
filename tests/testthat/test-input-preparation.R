test_that("matrix input is retained and integer storage is audited", {
  f <- combatrefql_test_fixtures()[["numeric_balanced"]]
  x <- f$counts; storage.mode(x) <- "integer"
  out <- combatrefql:::prepare_counts_input(x, f$batch)
  expect_equal(out$value, f$counts)
  expect_type(out$value, "double")
  expect_true("integer_to_double" %in% out$actions$action)
})

test_that("data frame, gene ID column, and numeric characters are repaired", {
  f <- combatrefql_test_fixtures()[["numeric_balanced"]]
  x <- data.frame(gene_id = rownames(f$counts), f$counts, check.names = FALSE)
  rownames(x) <- NULL
  x[-1] <- lapply(x[-1], as.character)
  out <- combatrefql:::prepare_counts_input(x, f$batch)
  expect_equal(out$value, f$counts)
  expect_true(all(c("first_column_to_gene_ids", "data_frame_to_numeric_matrix") %in% out$actions$action))
})

test_that("one-level groups and method spellings are canonicalized", {
  f <- combatrefql_test_fixtures()[["numeric_balanced"]]
  p <- combatrefql:::prepare_combat_inputs(f$counts, f$batch, rep("one", 8))
  expect_null(p$group)
  expect_true("drop_one_level_group" %in% p$actions$action)
  expect_equal(combatrefql:::canonical_method("TMMWSP", c("TMMwsp", "TMM", "none"), "normalization")$value, "TMMwsp")
  expect_equal(combatrefql:::canonical_method("mid-P", c("midp", "upper"), "quantile_method")$value, "midp")
})

test_that("unambiguous transpose and named sample reorder are audited", {
  f <- combatrefql_test_fixtures()[["numeric_balanced"]]
  x <- t(f$counts)
  batch <- f$batch; names(batch) <- rev(rownames(x))
  out <- combatrefql:::prepare_counts_input(x, batch)
  expect_equal(dim(out$value), dim(f$counts))
  expect_true("transpose" %in% out$actions$action)
})

test_that("lossy count repairs require explicit opt-in and report magnitude", {
  f <- combatrefql_test_fixtures()[["numeric_balanced"]]
  fractional <- f$counts; fractional[1, 1:2] <- c(10.2, 11.6)
  expect_error(combatrefql:::prepare_counts_input(fractional, f$batch), "fractional_counts")
  rounded <- combatrefql:::prepare_counts_input(fractional, f$batch, "round")
  action <- rounded$actions[rounded$actions$action == "round_fractional_counts", ]
  expect_true(action$lossy)
  expect_equal(action$affected, "2")
  expect_match(action$reason, "median.*max")

  duplicated <- rbind(f$counts[1, ], f$counts)
  rownames(duplicated)[1:2] <- "same"
  expect_error(combatrefql:::prepare_counts_input(duplicated, f$batch), "duplicate_genes")
  summed <- combatrefql:::prepare_counts_input(duplicated, f$batch, duplicate_genes = "sum")
  expect_true(any(summed$actions$action == "sum_duplicate_genes" & summed$actions$lossy))
})

test_that("quiet API retains audits while normal API announces conversions", {
  f <- combatrefql_test_fixtures()[["numeric_balanced"]]
  x <- f$counts; storage.mode(x) <- "integer"
  quiet <- expect_silent(combat_ref(x, f$batch, f$group, verbosity = "quiet"))
  expect_true(nrow(quiet$diagnostics$input_actions) > 0)
  expect_message(combat_ref(x, f$batch, f$group, verbosity = "normal"), "edgeR uses numeric doubles")
})

test_that("unsafe inputs remain errors", {
  f <- combatrefql_test_fixtures()[["numeric_balanced"]]
  bad <- f$counts; bad[1] <- Inf
  expect_error(combatrefql:::prepare_counts_input(bad, f$batch), "infinite")
  bad <- f$counts; colnames(bad)[2] <- colnames(bad)[1]
  expect_error(combatrefql:::prepare_counts_input(bad, f$batch), "duplicated sample")
  chars <- matrix("not-a-number", 2, 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
  expect_error(combatrefql:::prepare_counts_input(chars, c("a", "b")), "cannot be parsed")
  expect_error(combatrefql:::prepare_counts_input(f$counts, f$batch, input_scale = "TPM"), "not raw counts")
})
