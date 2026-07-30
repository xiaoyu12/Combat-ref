test_that("design construction preserves batch and group structure", {
  fixture <- combatrefql_legacy_fixtures()[["character_unbalanced"]]
  design <- combatrefql:::build_combat_design(fixture$batch, fixture$group)

  expect_equal(nrow(design), ncol(fixture$counts))
  expect_equal(qr(design)$rank, ncol(design))
  expect_true(any(grepl("batch", colnames(design))))
  expect_true(any(grepl("group", colnames(design))))
})

test_that("quantile mapping is deterministic and count-like", {
  counts <- matrix(c(0, 1, 2, 100, 1000, 5), nrow = 2)
  old_mu <- matrix(c(1, 2, 3, 100, 1000, 10), nrow = 2)
  new_mu <- old_mu + 1
  mapped <- combatrefql:::match_quantiles(
    counts_sub = counts,
    old_mu = old_mu,
    old_phi = c(0.1, 0.2),
    new_mu = new_mu,
    new_phi = c(0.1, 0.2),
    keep_zero = TRUE, method = "midp"
  )

  expect_equal(mapped, combatrefql:::match_quantiles(counts, old_mu, c(0.1, 0.2), new_mu, c(0.1, 0.2), TRUE, "midp"))
  expect_true(all(mapped >= 0))
  expect_true(all(abs(mapped - round(mapped)) < 1e-8))
  expect_equal(mapped[counts == 0], counts[counts == 0])
})
