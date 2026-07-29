test_that("legacy fixtures remain fixed", {
  fixtures <- combatrefql_legacy_fixtures()

  for (nm in names(fixtures)) {
    fixture <- fixtures[[nm]]
    observed <- legacy_combat_ref(
      counts = fixture$counts,
      batch = fixture$batch,
      group = fixture$group,
      genewise.disp = FALSE
    )

    expect_equal(observed, fixture$expected_adjusted)
    expect_equal(dim(observed), dim(fixture$counts))
    expect_equal(rownames(observed), rownames(fixture$counts))
    expect_equal(colnames(observed), colnames(fixture$counts))
    expect_true(all(observed >= 0))
    expect_true(all(abs(observed - round(observed)) < 1e-8))

    repeated <- legacy_combat_ref(
      counts = fixture$counts,
      batch = fixture$batch,
      group = fixture$group,
      genewise.disp = FALSE
    )
    expect_equal(repeated, observed)
  }
})
