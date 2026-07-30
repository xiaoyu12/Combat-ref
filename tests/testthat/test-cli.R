write_cli_fixture <- function(tmp, fixture, shuffled = FALSE) {
  count_file <- file.path(tmp, "counts.tsv")
  metadata_file <- file.path(tmp, "metadata.csv")

  utils::write.table(
    fixture$counts,
    file = count_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  metadata <- data.frame(
    batch = fixture$batch,
    condition = fixture$group,
    row.names = colnames(fixture$counts)
  )
  if (isTRUE(shuffled)) {
    metadata <- metadata[rev(rownames(metadata)), , drop = FALSE]
  }
  utils::write.csv(metadata, metadata_file, quote = FALSE)
  list(counts = count_file, metadata = metadata_file)
}

run_cli <- function(args) {
  script <- system.file("scripts", "combat-ref", package = "combatrefql")
  if (!nzchar(script)) {
    script <- normalizePath(test_path("../../inst/scripts/combat-ref"))
  }
  suppressWarnings(
    system2(
      file.path(R.home("bin"), "Rscript"),
      c(script, args),
      stdout = TRUE,
      stderr = TRUE
    )
  )
}

test_that("CLI shows help", {
  out <- run_cli("--help")
  expect_null(attr(out, "status"))
  expect_true(any(grepl("--counts", out, fixed = TRUE)))
  expect_true(any(grepl("--metadata", out, fixed = TRUE)))
})

test_that("CLI shows version", {
  out <- run_cli("--version")
  expect_null(attr(out, "status"))
  expect_true(any(grepl("0.1.0", out, fixed = TRUE)))
})

test_that("CLI validates missing files", {
  out <- run_cli(c(
    "--counts", "missing.tsv",
    "--metadata", "missing.csv",
    "--batch-column", "batch",
    "--output", tempfile()
  ))
  expect_equal(attr(out, "status"), 1L)
  expect_true(any(grepl("Input file not found", out)))
})

test_that("CLI delegates to package implementation", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  tmp <- tempdir()
  files <- write_cli_fixture(tmp, fixture)
  output_file <- file.path(tmp, "adjusted.tsv")

  out <- run_cli(c(
    "--counts", files$counts,
    "--metadata", files$metadata,
    "--batch-column", "batch",
    "--group-column", "condition",
    "--output", output_file,
    "--verbose"
  ))

  expect_null(attr(out, "status"))
  expect_true(file.exists(output_file))
  expect_true(any(grepl("Selected reference batch", out)))

  cli_counts <- as.matrix(utils::read.csv(output_file, sep = "\t", row.names = 1, check.names = FALSE))
  storage.mode(cli_counts) <- "numeric"
  api_counts <- combat_ref(fixture$counts, fixture$batch, fixture$group)$adjusted_counts
  expect_equal(cli_counts, api_counts)
})

test_that("CLI quiet mode suppresses routine output", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  tmp <- tempdir()
  files <- write_cli_fixture(tmp, fixture)
  output_file <- file.path(tmp, "adjusted-quiet.tsv")

  out <- run_cli(c(
    "--counts", files$counts,
    "--metadata", files$metadata,
    "--batch-column", "batch",
    "--group-column", "condition",
    "--output", output_file,
    "--quiet"
  ))

  expect_null(attr(out, "status"))
  expect_length(out, 0)
  expect_true(file.exists(output_file))
})

test_that("CLI reports metadata reordering and rejects silent reordering", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  tmp <- tempdir()
  files <- write_cli_fixture(tmp, fixture, shuffled = TRUE)
  output_file <- file.path(tmp, "adjusted-reordered.tsv")

  out <- run_cli(c(
    "--counts", files$counts,
    "--metadata", files$metadata,
    "--batch-column", "batch",
    "--output", output_file,
    "--verbose"
  ))
  expect_null(attr(out, "status"))
  expect_true(any(grepl("Reordering metadata rows", out)))

  quiet_out <- run_cli(c(
    "--counts", files$counts,
    "--metadata", files$metadata,
    "--batch-column", "batch",
    "--output", tempfile(),
    "--quiet"
  ))
  expect_equal(attr(quiet_out, "status"), 1L)
  expect_true(any(grepl("row order differs", quiet_out)))
})
