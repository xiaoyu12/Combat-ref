test_that("CLI delegates to package implementation", {
  fixture <- combatrefql_legacy_fixtures()[["numeric_balanced"]]
  tmp <- tempdir()
  count_file <- file.path(tmp, "counts.mat")
  sample_file <- file.path(tmp, "samples.csv")
  output_file <- file.path(tmp, "adjusted.mat")

  utils::write.table(
    fixture$counts,
    file = count_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  samples <- data.frame(
    batch = fixture$batch,
    condition = fixture$group,
    row.names = colnames(fixture$counts)
  )
  utils::write.csv(samples, sample_file, quote = FALSE)

  script <- normalizePath(test_path("../../inst/scripts/combat-ref"))
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(script, count_file, sample_file, output_file),
    stdout = TRUE,
    stderr = TRUE
  )

  expect_null(attr(status, "status"))
  expect_true(file.exists(output_file))

  cli_counts <- as.matrix(utils::read.csv(output_file, sep = "\t", row.names = 1, check.names = FALSE))
  storage.mode(cli_counts) <- "numeric"
  api_counts <- combat_ref(fixture$counts, fixture$batch, fixture$group)$adjusted_counts

  expect_equal(cli_counts, api_counts)
  expect_true(any(grepl("Reference batch", status)))
})
