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
  expect_true(any(grepl("--validate-only", out, fixed = TRUE)))
  expect_true(any(grepl("Exit codes", out, fixed = TRUE)))
})

test_that("CLI normal output is structured and no-color is plain text", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  files <- write_cli_fixture(tempdir(), fixture)
  out <- run_cli(c("--counts", files$counts, "--metadata", files$metadata,
    "--batch-column", "batch", "--output", tempfile(), "--no-color"))
  expect_null(attr(out, "status"))
  expect_true(all(vapply(c("Input", "Configuration", "Design", "Normalization",
    "Reference selection", "Correction", "Output", "Completed"),
    function(x) any(grepl(x, out, fixed = TRUE)), logical(1))))
  expect_false(any(grepl("\\033", out)))
})

test_that("CLI writes JSON diagnostics, gene status, and log", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  files <- write_cli_fixture(tempdir(), fixture)
  paths <- c(tempfile("out", fileext = ".tsv"), tempfile("diag", fileext = ".json"),
             tempfile("status", fileext = ".tsv"), tempfile("log", fileext = ".log"))
  out <- run_cli(c("--counts", files$counts, "--metadata", files$metadata,
    "--batch-column", "batch", "--output", paths[1],
    "--diagnostics-json", paths[2], "--gene-status", paths[3],
    "--log-file", paths[4], "--no-color"))
  expect_null(attr(out, "status"))
  expect_true(all(file.exists(paths)))
  json <- jsonlite::read_json(paths[2], simplifyVector = TRUE)
  expect_equal(json$reference_selection_frequency, NULL)
  expect_true("input_actions" %in% names(json))
  expect_match(paste(readLines(paths[4]), collapse = "\n"), "Completed")
})

test_that("CLI validate-only performs no correction write", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  files <- write_cli_fixture(tempdir(), fixture)
  output <- tempfile()
  out <- run_cli(c("--counts", files$counts, "--metadata", files$metadata,
    "--batch-column", "batch", "--validate-only", "--no-color"))
  expect_null(attr(out, "status"))
  expect_true(any(grepl("Validation completed", out)))
  expect_false(file.exists(output))
})

test_that("CLI uses stable input and design exit codes", {
  missing <- run_cli(c("--counts", "missing.tsv", "--metadata", "missing.csv",
    "--batch-column", "batch", "--output", tempfile()))
  expect_equal(attr(missing, "status"), 2L)

  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  tmp <- tempdir(); files <- write_cli_fixture(tmp, fixture)
  meta <- utils::read.csv(files$metadata, row.names = 1, check.names = FALSE)
  meta$condition <- meta$batch
  utils::write.csv(meta, files$metadata, quote = FALSE)
  confounded <- run_cli(c("--counts", files$counts, "--metadata", files$metadata,
    "--batch-column", "batch", "--group-column", "condition", "--output", tempfile()))
  expect_equal(attr(confounded, "status"), 3L)
})

test_that("CLI classifies preparation failures as input errors", {
  fixture <- combatrefql_test_fixtures()[["numeric_balanced"]]
  fixture$counts[1, 1] <- fixture$counts[1, 1] + 0.25
  files <- write_cli_fixture(tempdir(), fixture)
  out <- run_cli(c("--counts", files$counts, "--metadata", files$metadata,
    "--batch-column", "batch", "--output", tempfile()))
  expect_equal(attr(out, "status"), 2L)
  expect_true(any(grepl("fractional", out, ignore.case = TRUE)))
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
  expect_equal(attr(out, "status"), 2L)
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
    "--no-color"
  ))

  expect_null(attr(out, "status"))
  expect_true(file.exists(output_file))
  expect_true(any(grepl("Selected reference", out)))

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
  expect_true(any(grepl("different order", out)))

  quiet_out <- run_cli(c(
    "--counts", files$counts,
    "--metadata", files$metadata,
    "--batch-column", "batch",
    "--output", tempfile(),
    "--quiet"
  ))
  expect_null(attr(quiet_out, "status"))
  expect_length(quiet_out, 0)
})
