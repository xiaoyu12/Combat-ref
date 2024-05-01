rm(list=ls())
sapply(c("sva", "dplyr", "edgeR", "DESeq2", "ggplot2", "reshape2", "gridExtra", "scales", 
         "RUVSeq", "ggpubr", "BatchQC"), require, character.only=TRUE)
## Parameters (change paths when necessary)
data_dir <- "~/workspace/Combat-ref/nasa_data"  # path to the signature data (.rds)
source("~/workspace/Combat-ref/real_data_application/gfrn_helpers.R")  # path to gfrn_helpers.R
source("~/workspace/Combat-ref/ComBat_seq.R"); source("~/workspace/Combat-ref/helper_seq.R") 
source("~/workspace/Combat-ref/Combat_ref.R")

load(file.path(data_dir, "nasa.RData"))

group = data_list$Condition
batch = sapply(1:nrow(data_list), function(i) {paste(data_list[i, 2], data_list$Mission[i])})

#idx <- 1:nrow(data_list)
# Remove GLDS_168 RR1_NASA
idx <- which(data_list$Dataset.Accession..OSD.GLDS. != "GLDS_168") # | data_list$Mission != "RR3")
group = as.factor(group[idx])
batch = as.factor(batch[idx])
batch_ids = as.numeric(batch)

cts <- data[, idx]
keep <- apply(cts[, batch_ids==1], 1 ,function(x){!all(x<=3)})
for( i in 2:nlevels(batch)) {
  keep <- keep & apply(cts[, batch_ids==i], 1 ,function(x){!all(x<=3)})
}

cts <- as.matrix(data[keep, idx])

#PCA plot of unadjusted data
col_data <- data.frame(Batch=batch, Group=group) 
rownames(col_data) <- colnames(cts)

plot_Count_PCA <- function(cts, title, ntop=5000, normalization=TRUE, neg_y=1) {
  if (normalization) {
    cts_norm <- apply(cts, 2, function(x){x/sum(x)})
  } else {
    cts_norm = cts
  }
  seobj <- SummarizedExperiment(assays=cts_norm, colData=col_data)
  pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("Batch", "Group"), ntop=ntop)
  plt <- ggplot(pca_obj$data, aes(x=PC1, y=neg_y*PC2, color=Batch, shape=Group)) + 
    geom_point() + 
    labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
         y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
         title=title) 
  return(list("pca_obj"=pca_obj, "plt"=plt))
}
unadjusted <- plot_Count_PCA(cts, title="Unadjusted")
unadjusted$plt

# adjust the data using Combat
cts_combat <- sva::ComBat(cpm(cts, log=TRUE), batch=batch, mod=model.matrix(~group))
combat <- plot_Count_PCA(cts_combat, title="ComBat", normalization=FALSE)
combat$plt

# adjust the data using Combat_seq
cts_combatseq <- ComBat_seq(counts=cts, batch=batch, group=group, shrink=FALSE)
combatseq <- plot_Count_PCA(cts_combatseq, title="ComBat_seq")
combatseq$plt

# adjust the data using combat_ref
cts_combatnew <- ComBat_ref(counts = cts, batch = batch, group = group, genewise.disp = FALSE)
combat_new <- plot_Count_PCA(cts_combatnew, title="New ComBat_ref", neg_y=-1)
combat_new$plt

plt_PCA_full <- ggarrange(unadjusted$plt, combat$plt, combatseq$plt, combat_new$plt, ncol=1, nrow=4, common.legend=TRUE, legend="right")

varexp_full <- list(
  unadjusted=batchqc_explained_variation(cpm(cts, log=TRUE), condition=col_data$Group, batch=col_data$Batch)$explained_variation,
  combatseq=batchqc_explained_variation(cts_combat, condition=col_data$Group, batch=col_data$Batch)$explained_variation,
  combat=batchqc_explained_variation(cpm(cts_combatseq, log=TRUE), condition=col_data$Group, batch=col_data$Batch)$explained_variation,
  combatref=batchqc_explained_variation(cpm(cts_combatnew, log=TRUE), condition=col_data$Group, batch=col_data$Batch)$explained_variation
)
varexp_full_df <- melt(varexp_full)
varexp_full_df$L1 <- factor(varexp_full_df$L1, levels=c("unadjusted", "combat", "combatseq", "combatref"))
varexp_full_df$L1 <- plyr::revalue(varexp_full_df$L1, c("unadjusted"="Unadjusted", "combat"="ComBat",
                                                        "combatseq"="ComBat-Seq", "combatref"="New ComBat-ref"))
varexp_full_df$Var2 <- plyr::revalue(varexp_full_df$Var2, c("Full (Condition+Batch)"="Condition+Batch"))

plt_varexp_full <- ggplot(varexp_full_df, aes(x=Var2, y=value)) +
  geom_boxplot() +
  facet_wrap(~L1, nrow=4, ncol=1) +
  labs(y="Explained variation") +
  theme(axis.title.x = element_blank())

pca_ggar <- ggarrange(plt_PCA_full, plt_varexp_full, ncol=2, widths=c(0.55, 0.45))
ggsave("nasa_pca_variation_aggr.png", plot=pca_ggar, width=8, height=8, dpi=1000)
