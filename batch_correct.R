# Use Combat-ref to correct batch effects in RNA-seq count data
# This script takes command line arguments for a count matrix and a file describing
# biological conditions and batches of samples
#
# Usage: Rscript batch_correct.R <count_matrix_file> <sample_info_file>

sapply(c("edgeR", "matrixStats"), require, character.only=TRUE)
script_dir <- getwd()
source(file.path(script_dir, "helper_seq.R"))
source(file.path(script_dir, "Combat_ref.R"))

####  Command line arguments
command_args <- commandArgs(trailingOnly=TRUE)
if(length(command_args) < 2) {
  stop("Usage: Rscript batch_correct.R <count_matrix_file> <sample_info_file>")
}
count_file <- command_args[1]   # count matrix file
sample_file <- command_args[2]  # CSV file describing sample conditions and batches

if(!file.exists(count_file)) stop(paste("Count matrix file not found:", count_file))
if(!file.exists(sample_file)) stop(paste("Sample info file not found:", sample_file))

# Read count matrix
cts_mat <- read.csv(count_file, header = TRUE, sep="\t", row.names = 1, as.is=TRUE)
colnames(cts_mat) <- gsub("\\.", "-", colnames(cts_mat))
# Read the sample file
samples <- read.csv(sample_file, row.names=1, header=TRUE)

# Find batch and group columns (exact match first, then partial match)
find_column <- function(df, name) {
  cols <- colnames(df)
  # Exact match
  if(name %in% cols) return(name)
  # Case-insensitive partial match (column name contains the keyword)
  matches <- grep(name, cols, ignore.case=TRUE, value=TRUE)
  if(length(matches) == 1) return(matches)
  if(length(matches) > 1) stop(paste0("Ambiguous: multiple columns match '", name, "': ", paste(matches, collapse=", ")))
  stop(paste0("Sample file must contain a '", name, "' column. Found columns: ", paste(cols, collapse=", ")))
}
batch_col <- find_column(samples, "batch")
group_col <- find_column(samples, "group")
cat("Using batch column:", batch_col, "\n")
cat("Using group column:", group_col, "\n")

## Use the new ComBat-ref to adjust data
start_time <- Sys.time()
new_mat <- ComBat_ref(counts=cts_mat, batch=samples[[batch_col]], group=samples[[group_col]], genewise.disp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)

# Write the corrected matrix to a file with "corrected" added to the name
ext <- tools::file_ext(count_file)
new_filename <- paste0(tools::file_path_sans_ext(count_file), "_corrected.", ext)

write.table(new_mat, file=new_filename, col.names=TRUE, row.names=TRUE, sep="\t")
cat("Corrected data written to:", new_filename, "\n")