test_that("fractional counts round by default and are audited", {
  z <- combatrefql_test_fixtures()[[1]]; z$counts[1,1] <- z$counts[1,1] + 0.4
  fit <- combat_ref_ql(z$counts, z$batch, z$group, reference="1", verbosity="quiet")
  action <- fit@diagnostics$input_actions[1,]
  expect_equal(action$affected_values, 1)
  expect_equal(action$affected_genes, 1)
  expect_equal(action$maximum_absolute_change, 0.4)
  expect_true(action$lossy)
  expect_error(combat_ref_ql(z$counts,z$batch,z$group,fractional_counts="error",verbosity="quiet"), class="combatrefql_input_error")
})
test_that("invalid counts and normalized scale have typed errors", {
  z <- combatrefql_test_fixtures()[[1]]
  bad <- z$counts; bad[1,1] <- -1
  expect_error(combat_ref_ql(bad,z$batch,z$group,verbosity="quiet"), class="combatrefql_input_error")
  normalized <- matrix(seq(0.1, 9.6, length.out=length(z$counts)), nrow=nrow(z$counts), dimnames=dimnames(z$counts))
  expect_error(combat_ref_ql(normalized,z$batch,z$group,verbosity="quiet"), class="combatrefql_input_error")
})

test_that("plausible sparse fractional estimated counts remain accepted", {
  z <- combatrefql_test_fixtures()[[1]]
  estimated <- z$counts + ifelse(z$counts > 0, 0.25, 0)
  fit <- combat_ref_ql(estimated, z$batch, z$group, reference = "1",
    verbosity = "quiet")
  expect_s7_class(fit, CombatRefQLFit)
  expect_gt(fit@diagnostics$input_actions$affected_values, 0)
})

test_that("named metadata reorder and bad alignment is typed", {
  z <- combatrefql_test_fixtures()[[1]]
  b <- setNames(z$batch, colnames(z$counts)); b <- b[rev(names(b))]
  fit <- combat_ref_ql(z$counts,b,z$group,reference="1",verbosity="quiet")
  expect_identical(fit@samples$batch, as.character(factor(z$batch)))
  names(b)[1] <- "wrong"
  expect_error(combat_ref_ql(z$counts,b,z$group,verbosity="quiet"), class="combatrefql_alignment_error")
})

test_that("design and reference failures have fields", {
  z <- combatrefql_test_fixtures()[[1]]
  err <- expect_error(combat_ref_ql(z$counts,z$batch,z$batch,verbosity="quiet"), class="combatrefql_design_error")
  expect_true(!is.null(err$rank) && !is.null(err$residual_df))
  err <- expect_error(combat_ref_ql(z$counts,z$batch,z$group,reference="missing",verbosity="quiet"), class="combatrefql_reference_error")
  expect_identical(err$stage, "reference")
})
