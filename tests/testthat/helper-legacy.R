legacy_env <- new.env(parent = globalenv())
suppressPackageStartupMessages(library(edgeR))
source(test_path("../../R/legacy-fixtures.R"))
source(test_path("../../helper_seq.R"), local = legacy_env)
source(test_path("../../Combat_ref.R"), local = legacy_env)

legacy_combat_ref <- legacy_env$ComBat_ref
