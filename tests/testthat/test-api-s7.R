test_that("sole API returns validated lightweight S7 objects", {
  z <- combatrefql_test_fixtures()[[1]]
  fit <- combat_ref_ql(z$counts, z$batch, z$group, reference = "1", verbosity = "quiet")
  expect_true(S7::S7_inherits(fit, CombatRefQLFit))
  expect_identical(dim(fit@counts), dim(z$counts))
  expect_identical(dimnames(fit@counts), dimnames(z$counts))
  expect_true(is.double(fit@counts))
  expect_true(S7::S7_inherits(summary(fit), CombatRefQLSummary))
  expect_false(any(c("fitted_means", "reference_means", "counts") %in% names(fit@model)))
  expect_error(CombatRefQLFit(counts = matrix(-1, 1, 1), specification = data.frame(reference_batch="x"), samples=data.frame(sample="s",batch="x"), gene_status=data.frame(gene="g"), diagnostics=list(), call=quote(x()), model=list()), "invalid")
})

test_that("old API and plotting are absent", {
  expect_false(exists("combat_ref", asNamespace("combatrefql"), inherits = FALSE))
  expect_false(any(grepl("combatref", methods("plot"))))
  expect_setequal(getNamespaceExports("combatrefql"), c("combat_ref_ql", "CombatRefQLFit", "CombatRefQLSummary"))
})

test_that("fit and summary printing orient without replay", {
  z <- combatrefql_test_fixtures()[[1]]
  fit <- combat_ref_ql(z$counts, z$batch, z$group, reference = "1", verbosity = "quiet")
  out <- testthat::capture_messages(print(fit))
  expect_match(paste(out, collapse="\n"), "(?s)ComBat-refQL fit.*Data.*Model.*Result", perl=TRUE)
  expect_equal(sum(grepl("Reference", out)), 1L)
  sum_out <- testthat::capture_messages(print(summary(fit)))
  expect_match(paste(sum_out, collapse="\n"), "(?s)Design.*Reference.*QL model.*Gene outcomes.*Timing", perl=TRUE)
})
