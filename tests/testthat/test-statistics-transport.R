test_that("v4 QL fit keeps NB and QL dispersion strictly separate", {
  z <- combatrefql_test_fixtures()[[1]]
  fit <- combat_ref_ql(z$counts,z$batch,z$group,reference="1",verbosity="quiet")
  expect_true(isTRUE(fit@diagnostics$ql$robust) && !fit@diagnostics$ql$legacy)
  expect_length(fit@model$nb_dispersion, sum(fit@gene_status$fit_valid))
  expect_length(fit@model$posterior_ql_dispersion, sum(fit@gene_status$fit_valid))
  expect_false(identical(fit@model$nb_dispersion, fit@model$posterior_ql_dispersion))
  code <- paste(deparse(body(combatrefql:::midp_transport)), collapse=" ")
  expect_match(code, "source_dispersion")
  expect_false(grepl("ql_dispersion|s2.post|s2.prior", code))
})

test_that("mid-P transport is deterministic, model consistent, and double safe", {
  y <- matrix(c(0,1,10,100,2^31+10), nrow=1)
  old <- matrix(c(20,20,20,100,2^31+10),nrow=1); new <- old * 1.2
  a <- combatrefql:::midp_transport(y,old,new,0.2,0.2)$counts
  b <- combatrefql:::midp_transport(y,old,new,0.2,0.2)$counts
  expect_identical(a,b); expect_true(is.double(a)); expect_gt(a[5],2^31)
  extremes <- combatrefql:::midp_transport(matrix(c(0, 1, 1e6), 1),
    matrix(c(1e-8, 1e8, 1e6), 1), matrix(c(1e-7, 1e7, 2e6), 1),
    matrix(c(1e-8, 2, 0.01), 1), matrix(c(1e-8, 1, 0.02), 1))$counts
  expect_true(all(is.finite(extremes) & extremes >= 0 & extremes == round(extremes)))
})

test_that("chunking is equivalent and unsupported genes remain unchanged", {
  z <- combatrefql_test_fixtures()[[1]]
  a <- combat_ref_ql(z$counts,z$batch,z$group,reference="1",chunk_size=1L,verbosity="quiet")
  b <- combat_ref_ql(z$counts,z$batch,z$group,reference="1",chunk_size=5000L,verbosity="quiet")
  expect_identical(a@counts,b@counts)
  unsupported <- !a@gene_status$mapping_valid
  expect_equal(a@counts[unsupported,],z$counts[unsupported,])
  expect_identical(a@counts[, z$batch == 1], z$counts[, z$batch == 1])
  expect_true(all(is.finite(a@counts) & a@counts >= 0 & a@counts == round(a@counts)))
})

test_that("genes without support in one batch are identified", {
  batch <- factor(rep(c("A", "B"), each = 3))
  counts <- matrix(1:60, 10, 6)
  counts[1, batch == "B"] <- 0
  keep <- combatrefql:::assess_support(counts, batch)
  expect_false(keep[1])
  expect_true(all(keep[-1]))
})

test_that("mapping failures retain complete original rows and drive outcomes", {
  counts <- matrix(c(10, 12, 20, 22, 5, 6, 7, 8), nrow = 2, byrow = TRUE,
    dimnames = list(c("ok", "failed"), paste0("s", 1:4)))
  batch <- factor(c("A", "A", "B", "B"))
  design <- stats::model.matrix(~ 0 + batch)
  colnames(design) <- c("batchA", "batchB")
  fit <- list(
    coefficients = matrix(c(log(10), log(20), Inf, log(7)), nrow = 2,
      byrow = TRUE),
    offset = matrix(0, 2, 4), dispersion = 0.2,
    fitted.values = matrix(c(10, 12, 20, 22, 5, 6, Inf, 8), nrow = 2,
      byrow = TRUE)
  )
  mapped <- combatrefql:::transport_counts(counts, c(TRUE, TRUE), batch, "A", fit,
    shrinkage = list(posterior = list(B = c(0, 0))),
    dispersion = list(dispersion = matrix(0.2, 2, 4),
      reference = c(0.2, 0.2)), chunk_size = 2L)
  expect_true(mapped$mapping_valid[1])
  expect_false(mapped$mapping_valid[2])
  expect_identical(mapped$counts[2, ], counts[2, ])
  expect_equal(mapped$failures$gene, "failed")

  status <- combatrefql:::make_gene_status(counts, batch, c(TRUE, TRUE),
    mapped$mapping_valid)
  outcomes <- combatrefql:::make_outcomes(status)
  expect_equal(status$status, c("adjusted", "failed"))
  expect_equal(outcomes$genes, c(1L, 0L, 1L))
})

test_that("mapping slices fitted means inside each gene chunk", {
  code <- paste(deparse(body(combatrefql:::transport_counts)), collapse = " ")
  expect_match(code, "fit\\$fitted.values\\[chunk")
  expect_match(code, "dispersion\\$dispersion\\[chunk")
})

test_that("automatic reference and verbosity modes are deterministic", {
  z <- combatrefql_test_fixtures()[[1]]
  quiet <- testthat::capture_messages(a <- combat_ref_ql(z$counts,z$batch,z$group,verbosity="quiet"))
  normal <- testthat::capture_messages(b <- combat_ref_ql(z$counts,z$batch,z$group,verbosity="normal"))
  verbose <- testthat::capture_messages(c <- combat_ref_ql(z$counts,z$batch,z$group,verbosity="verbose"))
  expect_length(quiet,0L); expect_gt(length(normal),0L); expect_gte(length(verbose),length(normal))
  expect_false(any(grepl("reference candidate|gene chunks", normal)))
  expect_true(any(grepl("reference candidate", verbose)) && any(grepl("gene chunks", verbose)))
  expect_identical(a@counts,b@counts); expect_identical(b@counts,c@counts)
  expect_equal(sum(a@diagnostics$batches$selected),1L)
})
