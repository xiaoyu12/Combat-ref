#!/usr/bin/env Rscript

message("`batch_correct.R` is deprecated; use `inst/scripts/combat-ref` or `combatrefql::combat_ref()`.")

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "batch_correct.R"
script <- file.path(dirname(normalizePath(this_file, mustWork = FALSE)), "inst", "scripts", "combat-ref")
if (!file.exists(script)) {
  script <- file.path(getwd(), "inst", "scripts", "combat-ref")
}

source(script, local = new.env(parent = globalenv()))
