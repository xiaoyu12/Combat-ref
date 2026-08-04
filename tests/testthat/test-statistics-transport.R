test_that("v4 QL fit keeps NB and QL dispersion strictly separate", {
  z <- combatrefql_test_fixtures()[[1]]
  fit <- combat_ref_ql(z$counts,z$batch,z$group,reference="1",verbosity="quiet")
  expect_true(isTRUE(fit@diagnostics$ql$robust) && !fit@diagnostics$ql$legacy)
  expect_length(fit@model$nb_dispersion, 1L)
  expect_length(fit@model$posterior_ql_dispersion, sum(fit@gene_status$fit_valid))
  expect_false(identical(fit@model$nb_dispersion, fit@model$posterior_ql_dispersion))
  code <- paste(deparse(body(combatrefql:::midp_transport)), collapse=" ")
  expect_match(code, "nb_dispersion")
  expect_false(grepl("ql_dispersion|s2.post|s2.prior", code))
})

test_that("mid-P transport is deterministic, zero preserving, and double safe", {
  y <- matrix(c(0,1,10,100,2^31+10), nrow=1)
  old <- matrix(c(20,20,20,100,2^31+10),nrow=1); new <- old * 1.2
  a <- combatrefql:::midp_transport(y,old,new,0.2)$counts
  b <- combatrefql:::midp_transport(y,old,new,0.2)$counts
  expect_identical(a,b); expect_equal(a[1],0); expect_true(is.double(a)); expect_gt(a[5],2^31)
})

test_that("chunking is equivalent and unsupported genes remain unchanged", {
  z <- combatrefql_test_fixtures()[[1]]
  a <- combat_ref_ql(z$counts,z$batch,z$group,reference="1",chunk_size=1L,verbosity="quiet")
  b <- combat_ref_ql(z$counts,z$batch,z$group,reference="1",chunk_size=5000L,verbosity="quiet")
  expect_identical(a@counts,b@counts)
  unsupported <- !a@gene_status$mapping_valid
  expect_equal(a@counts[unsupported,],z$counts[unsupported,])
  expect_equal(a@counts[z$counts==0],rep(0,sum(z$counts==0)))
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
