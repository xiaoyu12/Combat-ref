# Use Combat-ref to correct batch effects in RNA-seq count data
# This scripts take command line arguments for a count matrix and a file describing
# biological conditions and batches of samples
sapply(c("sva", "dplyr", "DESeq2", "ggplot2", "reshape2", "gridExtra", "scales", 
         "RUVSeq", "ggpubr", "BatchQC"), require, character.only=TRUE)
script_dir <- getwd()
source(file.path(script_dir, "helper_seq.R")) 
source(file.path(script_dir, "Combat_ref.R"))

####  Command line arguments
command_args <- commandArgs(trailingOnly=TRUE)
count_file <- command_args[1]   # count matrix file
sample_file <- command_args[2]   # CSV file describing sample conditions and batches
#count_file <- "gfrn.mat"
#sample_file <- "gfrn_samples.csv"

# Read count matrix
cts_mat <- read.csv(count_file, header = TRUE, sep="\t", row.names = 1, as.is=TRUE)
colnames(cts_mat) <- gsub("\\.", "-", colnames(cts_mat))
# Read the sample file
samples <- read.csv(sample_file, row.names=1, header=TRUE)

## Use the new ComBat-ref to adjust data
start_time <- Sys.time()
new_mat <- ComBat_ref(counts=cts_mat, batch=samples$batch, group=samples$group, genewise.disp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)

# Write the corrected matrix to a file with "corrected" added to the name
# Get the file extension
ext <- tools::file_ext(count_file)
new_filename <- paste0(tools::file_path_sans_ext(count_file), "_corrected.", ext)

write.table(new_mat, file=new_filename, col.names=TRUE, row.names=TRUE, sep="\t")