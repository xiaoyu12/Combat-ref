test_that("counts must be valid", {
  fixture <- combatrefql_legacy_fixtures()[["numeric_balanced"]]
  counts <- fixture$counts

  bad <- counts
  bad[1, 1] <- -1
  expect_error(combat_ref(bad, fixture$batch, fixture$group), "negative")

  bad <- counts
  bad[1, 1] <- NA_real_
  expect_error(combat_ref(bad, fixture$batch, fixture$group), "missing or non-finite")

  bad <- counts
  bad[, 1] <- 0
  expect_error(combat_ref(bad, fixture$batch, fixture$group), "all-zero samples")

  bad <- matrix(as.character(counts), nrow = nrow(counts), dimnames = dimnames(counts))
  expect_error(combat_ref(bad, fixture$batch, fixture$group), "numeric")
})

test_that("metadata must align exactly", {
  fixture <- combatrefql_legacy_fixtures()[["numeric_balanced"]]

  expect_error(combat_ref(fixture$counts, fixture$batch[-1], fixture$group), "length")
  expect_error(combat_ref(fixture$counts, rep("one", ncol(fixture$counts)), fixture$group), "at least two")

  named_batch <- fixture$batch
  names(named_batch) <- rev(colnames(fixture$counts))
  expect_error(combat_ref(fixture$counts, named_batch, fixture$group), "same order")

  named_batch <- fixture$batch
  names(named_batch) <- colnames(fixture$counts)
  expect_no_error(combat_ref(fixture$counts, named_batch, fixture$group))
})

test_that("invalid library sizes fail clearly", {
  fixture <- combatrefql_legacy_fixtures()[["numeric_balanced"]]

  expect_error(combat_ref(fixture$counts, fixture$batch, fixture$group, lib_size = rep(1, 2)), "length")
  expect_error(combat_ref(fixture$counts, fixture$batch, fixture$group, lib_size = rep(0, 8)), "lib_size")
  expect_error(combat_ref(fixture$counts, fixture$batch, fixture$group, lib_size = rep(NA_real_, 8)), "lib_size")
})

test_that("confounded designs fail before fitting", {
  counts <- matrix(
    c(10, 12, 30, 32, 8, 9, 28, 30,
      4, 5, 14, 16, 3, 4, 15, 17,
      50, 55, 120, 130, 48, 52, 125, 135),
    nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("sample", 1:8))
  )
  batch <- rep(c("A", "B"), each = 4)
  group <- batch

  expect_error(combat_ref(counts, batch, group), "rank deficient")
})
