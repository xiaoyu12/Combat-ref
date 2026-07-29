library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

# Function to read and process the data
process_data <- function(base_dir, file_name) {
  file_path <- paste0(base_dir, "/", file_name)
  data <- read.csv(file_path)
  # Calculate mean FPR and TPR across all repetitions
  means <- colMeans(data[,2:ncol(data)])
  return(means)
}

# Function to plot data for each meanFC and dispFC combination
plot_data <- function(meanFC_levels, dispFC_levels, base_dir_fpr, base_dir_tpr) {
  plot_list <- list()
  for (meanFC in meanFC_levels) {
    for (dispFC in dispFC_levels) {
      fpr_file_name <- sprintf("fpr_sim_N12_meanFC%d_dispFC%d_biofc24_depth5.csv", meanFC, dispFC)
      tpr_file_name <- sprintf("tpr_sim_N12_meanFC%d_dispFC%d_biofc24_depth5.csv", meanFC, dispFC)
      
      fpr_means <- process_data(base_dir_fpr, fpr_file_name)
      tpr_means <- process_data(base_dir_tpr, tpr_file_name)
      
      # Filter to include only '.edgeR' methods
      methods <- names(fpr_means)
      edgeR_methods <- methods[grepl(".edgeR$", methods)| grepl(".lm$", methods)]
      fpr_edgeR <- fpr_means[edgeR_methods]
      tpr_edgeR <- tpr_means[edgeR_methods]
      
      # Create a data frame for plotting
      plot_data <- data.frame(
        Method = edgeR_methods,
        MeanFPR = fpr_edgeR,
        MeanTPR = tpr_edgeR
      )
      methods = c(
        "No batch effect",
        "With batch effect, no adjustment",
        "Treat batch as covariate in DE model",
        "Adjusted by ComBat on logCPM",
        "Adjusted by ComBat-Seq",
        "Adjusted by RUVSeq",
        "Adjusted by SVASeq",
        "Adjusted by new ComBat-ref",
        "Adjusted by NPM"
      )
      plot_data$Method <- factor(methods, levels = methods)
      
      if(meanFC > 10) {
         batch_mean_FC = meanFC/10.0
      } else {
        batch_mean_FC = meanFC/1.0
      }
      if(dispFC > 10) {
        batch_disp_FC = dispFC/10.0
      } else {
        batch_disp_FC = dispFC/1.0
      }
      # Generate plot
      # Generate a color palette with 9 colors
      colors <- brewer.pal(9, "Set1")
      p <- ggplot(plot_data, aes(x = MeanFPR, y = MeanTPR, label = Method)) +
        geom_point(size=3, stroke=1.5, aes(color = Method, shape = Method)) + #geom_text(vjust = 1.5, hjust = 1.5) + scale_color_manual(values = colors) +
        scale_shape_manual(values=0:8) +
        xlim(0, 0.25) +
        ylim(0.4, 1.0) +
        labs(title = paste("mean_FC =", batch_mean_FC, "disp_FC =", batch_disp_FC),
             x = "Mean False Positive Rate (FPR)",
             y = "Mean True Positive Rate (TPR)") +
        theme_bw()  #theme(legend.position = "right")
      
      plot_list[[paste("meanFC", meanFC, "dispFC", dispFC)]] <- p
    }
  }
  return(plot_list)
}

plot_fdr_data <- function(meanFC_levels, dispFC_levels, base_dir_fpr, base_dir_tpr, fdr=0.1) {
  plot_list <- list()
  for (meanFC in meanFC_levels) {
    for (dispFC in dispFC_levels) {
      fpr_file_name <- sprintf("fprADJ_sim_N12_meanFC%d_dispFC%d_biofc24_depth5.csv", meanFC, dispFC)
      tpr_file_name <- sprintf("tprADJ_sim_N12_meanFC%d_dispFC%d_biofc24_depth5.csv", meanFC, dispFC)
      
      fpr_file_path <- paste0(base_dir_fpr, "/", fpr_file_name)
      fpr_data <- read.csv(fpr_file_path)
      idx <- which(fpr_data[, "FDR.cutoff"] == fdr)
      if(length(idx) <= 0) {
        stop(paste("No matching fdr", as.character(fdr)))
      }
      fpr_data <- fpr_data[idx, ]
      tpr_file_path <- paste0(base_dir_tpr, "/", tpr_file_name)
      tpr_data <- read.csv(tpr_file_path)
      tpr_data <- tpr_data[idx, ]
      # Calculate mean FPR and TPR across all repetitions
      fpr_means <- colMeans(fpr_data[, 2:ncol(fpr_data)]) 
      tpr_means <- colMeans(tpr_data[, 2:ncol(tpr_data)])
      
      # Filter to include only '.edgeR' methods
      methods <- names(fpr_means)
      edgeR_methods <- methods[grepl(".edgeR$", methods) | grepl(".lm$", methods)]
      fpr_edgeR <- fpr_means[edgeR_methods]
      tpr_edgeR <- tpr_means[edgeR_methods]
      
      # Create a data frame for plotting
      plot_data <- data.frame(
        Method = edgeR_methods,
        MeanFPR = fpr_edgeR,
        MeanTPR = tpr_edgeR
      )
      methods = c(
        "No batch effect",
        "With batch effect, no adjustment",
        "Treat batch as covariate in DE model",
        "Adjusted by ComBat on logCPM",
        "Adjusted by ComBat-Seq",
        "Adjusted by RUVSeq",
        "Adjusted by SVASeq",
        "Adjusted by new ComBat-ref",
        "Adjusted by NPM"
      )
      plot_data$Method <- factor(methods, levels = methods)
      
      if(meanFC > 10) {
        batch_mean_FC = meanFC/10.0
      } else {
        batch_mean_FC = meanFC/1.0
      }
      if(dispFC > 10) {
        batch_disp_FC = dispFC/10.0
      } else {
        batch_disp_FC = dispFC/1.0
      }
      # Generate plot
      colors <- brewer.pal(9, "Set1")
      p <- ggplot(plot_data, aes(x = MeanFPR, y = MeanTPR, label = Method)) +
        geom_point(size=3, stroke=1.5, aes(color = Method, shape = Method)) + #geom_text(vjust = 1.5, hjust = 1.5) +
        scale_shape_manual(values=0:8) +
        xlim(0, 0.2) +
        ylim(0.0, 1.0) +
        labs(title = paste("mean_FC =", batch_mean_FC, "disp_FC =", batch_disp_FC),
             x = "Mean False Positive Rate (FPR)",
             y = "Mean True Positive Rate (TPR)") +
        theme_bw()  #theme(legend.position = "right")
      
      plot_list[[paste("meanFC", meanFC, "dispFC", dispFC)]] <- p
    }
  }
  return(plot_list)
}

# Define your paths and level settings
base_dir_fpr <- "./simulations"
base_dir_tpr <- "./simulations"
meanFC_levels <- c(1, 15, 2, 24)
dispFC_levels <- c(1, 2, 3, 4)

# Generate plots
plots <- plot_data(meanFC_levels, dispFC_levels, base_dir_fpr, base_dir_tpr)
plots_arr <- ggarrange(plotlist = plots, ncol=4, nrow=4, legend="right", common.legend = TRUE)
ggsave("tpr_fpr_pval_0.05.png", plots_arr, width=18, height = 16, dpi = 1000)

plots_fdr <- plot_fdr_data(meanFC_levels, dispFC_levels, base_dir_fpr, base_dir_tpr, fdr = 0.10)
plots_fdr_arr <- ggarrange(plotlist = plots_fdr, ncol=4, nrow=4, legend="right", common.legend = TRUE)
ggsave("tpr_fpr_fdr_0.1.png", plots_fdr_arr, width=18, height=16, dpi = 1000)

