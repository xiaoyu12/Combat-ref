devtools::load_all(".", quiet = TRUE)
source("validation/real_data/common.R")

url <- paste0("https://raw.githubusercontent.com/xiaoyu12/Combat-ref/",
  "f0c6d3f313e9f1f29de2c8c8f12b98113443be55/nasa_data/nasa.RData")
path <- download_cached(url, "validation/.cache/nasa/nasa.RData")
environment <- new.env(parent = emptyenv())
load(path, envir = environment)
metadata <- environment$data_list
selected <- metadata$Dataset.Accession..OSD.GLDS. != "GLDS_168"
batch <- interaction(metadata$Dataset.Accession..OSD.GLDS., metadata$Mission,
  drop = TRUE, sep = "_")
counts <- as.matrix(environment$data[, selected, drop = FALSE])
print(run_validation(counts, batch[selected], metadata$Condition[selected], "NASA GeneLab",
  "validation/real_data/nasa"))
