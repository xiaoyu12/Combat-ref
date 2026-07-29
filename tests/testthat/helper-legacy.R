legacy_env <- new.env(parent = globalenv())
suppressPackageStartupMessages(library(edgeR))

legacy_helper <- test_path("../../helper_seq.R")
legacy_core <- test_path("../../Combat_ref.R")

if (file.exists(legacy_helper) && file.exists(legacy_core)) {
  source(test_path("../../R/legacy-fixtures.R"))
  source(legacy_helper, local = legacy_env)
  source(legacy_core, local = legacy_env)
  legacy_combat_ref <- legacy_env$ComBat_ref
} else {
  legacy_combat_ref <- function(...) {
    out <- combatrefql:::ComBat_ref(...)
    attr(out, "reference_batch") <- NULL
    out
  }
}
