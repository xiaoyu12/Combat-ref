devtools::load_all(".", quiet = TRUE)
source("validation/real_data/common.R")

url <- paste0("https://raw.githubusercontent.com/xiaoyu12/Combat-ref/",
  "f0c6d3f313e9f1f29de2c8c8f12b98113443be55/real_data_application/signature_data.rds")
path <- download_cached(url, "validation/.cache/gfrn/signature_data.rds")
object <- readRDS(path)
metadata <- as.data.frame(SummarizedExperiment::colData(object))
condition <- as.character(metadata$group)
selected <- condition %in% c("gfp_for_egfr", "gfp18", "gfp30", "her2", "egfr", "kraswt")
condition[condition %in% c("gfp_for_egfr", "gfp18", "gfp30")] <- "gfp"
counts <- SummarizedExperiment::assay(object, "counts")[, selected, drop = FALSE]
print(run_validation(counts, metadata$batch[selected], condition[selected], "GFRN",
  "validation/real_data/gfrn"))
