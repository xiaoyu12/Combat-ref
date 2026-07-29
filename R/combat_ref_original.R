# Correct batch effects for RNA-seq count data using a negative binomial distribution
# Upstream-derived ComBat-ref statistical implementation.
#
# Adapted from the ComBat-ref research implementation by Xiaoyu Zhang and
# collaborators. Redistribution terms for the upstream-derived implementation
# require manual clarification; see LICENSE_NOTES.md.
#
# Adjust count towards the reference batch with minimum estimated dispersion

#' @param counts Raw count matrix from genomic studies (dimensions gene x sample) 
#' @param batch Batch covariate (only one batch allowed)
#' @param group Vector / factor for condition of interest 
#' @param full_mod Boolean, if TRUE include condition of interest in model
#' @param gene.subset.n Number of genes to use in empirical Bayes estimation, only useful when shrink = TRUE
#' @param genewise.disp Compute dispersion parameter for each gene = FALSE
#' @return data A probe x sample count matrix, adjusted for batch effects.
#' 
#' @examples 
#' #combatseq_new <- ComBat_ref(counts=cts_sub, batch=batch_sub, group=group_sub, genewise.disp=FALSE)
#' #adj_counts <- ComBat_ref(counts=cts, batch=batch, group=group,  genewise.disp=FALSE)
#' 
#' @noRd
#' 

ComBat_ref <- function(counts, batch, group=NULL, full_mod=TRUE,
                       gene.subset.n=NULL, genewise.disp=FALSE,
                       reference_batch=NULL, lib_size=NULL)
{
  ########  Preparation  ########
  counts <- as.matrix(counts)

  ## Does not support 1 sample per batch yet
  batch <- as.factor(batch)
  if(any(table(batch)<=1)){
    stop("ComBat-seq doesn't support 1 sample per batch yet")
  }

  ## Remove genes with only 0 counts in any batch
  keep_lst <- lapply(levels(batch), function(b){
    which(apply(counts[, batch==b], 1, function(x){!all(x==0)}))
  })
  keep <- Reduce(intersect, keep_lst)
  rm <- setdiff(1:nrow(counts), keep)
  countsOri <- counts
  counts <- counts[keep, ]

  # require bioconductor 3.7, edgeR 3.22.1
  dge_obj <- DGEList(counts=counts, lib.size=lib_size)

  ## Prepare characteristics on batches
  n_batch <- nlevels(batch)  # number of batches
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch
  n_batches <- sapply(batches_ind, length)
  n_sample <- sum(n_batches)
  cat("Found",n_batch,'batches\n')

  ## Build covariate model matrix
  group <- as.factor(group)
  if(full_mod & nlevels(group)>1){
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~0+group)
  }else{
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data=as.data.frame(t(counts)))
  }

  ## Helper: build design matrix from current batch factor
  build_design <- function(batch, mod) {
    batchmod <- model.matrix(~batch)
    design <- cbind(batchmod, mod)
    check <- apply(design, 2, function(x) all(x == 1))
    design <- as.matrix(design[,!check])
    list(design=design, batchmod=batchmod)
  }

  dm <- build_design(batch, mod)
  design <- dm$design
  batchmod <- dm$batchmod

  ## Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")}
    if(ncol(design)>(n_batch+1)){
      if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")}}
  }

  cat("Estimating dispersions\n")
  ## Estimate common dispersion within each batch
  disp_common <- sapply(1:n_batch, function(i){
    if((n_batches[i] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){
      # not enough residual degree of freedom
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=NULL, subset=nrow(counts)))
    }else{
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=mod[batches_ind[[i]], ], subset=nrow(counts)))
    }
  })
  for(i in 1:n_batch) {
    cat("Batch", levels(batch)[i], "dispersion =", disp_common[i], "\n")
  }

  # Choose the batch with the smallest dispersion as the reference batch
  ref_batch <- if (is.null(reference_batch)) {
    which.min(disp_common)
  } else {
    match(reference_batch, levels(batch))
  }
  selected_reference_batch <- levels(batch)[ref_batch]
  cat("Reference batch: ", selected_reference_batch, "\n")

  # Reorder factor levels so the reference batch becomes the first level
  if(ref_batch != 1) {
    # Swap dispersion values
    disp_common[c(1, ref_batch)] <- disp_common[c(ref_batch, 1)]
    # Reorder factor levels (not rename) so ref_batch level comes first
    new_order <- seq_along(levels(batch))
    new_order[c(1, ref_batch)] <- new_order[c(ref_batch, 1)]
    batch <- factor(batch, levels = levels(batch)[new_order])
    # Rebuild batch indices and design matrix
    batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])})
    n_batches <- sapply(batches_ind, length)
    dm <- build_design(batch, mod)
    design <- dm$design
    batchmod <- dm$batchmod
  }
  
  if(genewise.disp) {  
    genewise_disp_lst <- lapply(1:n_batch, function(j){
      if((n_batches[j] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
        # not enough residual degrees of freedom - use the common dispersion
        return(rep(disp_common[j], nrow(counts)))
      }else{
        return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], design=mod[batches_ind[[j]], ], 
                                      dispersion=disp_common[j], prior.df=0))
      }
    })
  } else {
    genewise_disp_lst <- lapply(1:n_batch, function(j) {
      return(rep(disp_common[j], nrow(counts)))
    }) 
  }
  
  names(genewise_disp_lst) <- paste0('batch', levels(batch))
  ## construct dispersion matrix
  phi_matrix <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(k in 1:n_batch){
    phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k]) 
  }
  ########  Estimate parameters from NB GLM  ########
  cat("Fitting the GLM model\n")
  glm_f <- glmFit(dge_obj, design=design, dispersion=phi_matrix, prior.count=1e-4)
  
  gamma_hat <- as.matrix(glm_f$coefficients[, 1:(n_batch-1)])
  mu_hat <- glm_f$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)
  
  ########  Adjust the data  ########  
  cat("Adjusting the data\n")
  adjust_counts <- counts
  # adjust batches except for the first one (reference batch)
  for(kk in 2:n_batch) {
    counts_sub <- counts[, batches_ind[[kk]]]
    old_mu <- pmax(mu_hat[, batches_ind[[kk]]], 1e-4)      #numerical stability
    old_phi <- phi_hat[, kk]
    new_mu <- exp(log(old_mu)-vec2mat(gamma_hat[, kk-1], n_batches[kk]))
    # avoid exploding count (mu), new_mu shouldn't increase to more than the ref batch if gamma_hat < 0
    increased_genes <- which(gamma_hat[, kk-1] < -0.2)
    n_cols <- ncol(new_mu)
    ref_max <- vec2mat(matrixStats::rowMaxs(mu_hat[, batches_ind[[1]]]), n_cols)
    new_mu[increased_genes, ] <- pmin(new_mu[increased_genes, ], ref_max[increased_genes, ])
    new_phi <- phi_hat[, 1]
    adjust_counts[, batches_ind[[kk]]] <- match_quantiles(counts_sub=counts_sub,
                                                          old_mu=old_mu, old_phi=old_phi,
                                                          new_mu=new_mu, new_phi=new_phi, keep_zero=FALSE)
  }
  
  ## Add back genes with only 0 counts in any batch (so that dimensions won't change)
  adjust_counts_whole <- matrix(NA, nrow=nrow(countsOri), ncol=ncol(countsOri))
  dimnames(adjust_counts_whole) <- dimnames(countsOri)
  adjust_counts_whole[keep, ] <- adjust_counts
  adjust_counts_whole[rm, ] <- countsOri[rm, ]
  attr(adjust_counts_whole, "reference_batch") <- selected_reference_batch
  attr(adjust_counts_whole, "model") <- list(
    design = design,
    common_dispersions = stats::setNames(disp_common, levels(batch)),
    genewise_dispersions = genewise_disp_lst,
    phi_matrix = phi_matrix,
    coefficients = glm_f$coefficients,
    fitted_means = mu_hat,
    library_sizes = dge_obj$samples$lib.size,
    offset = dge_obj$offset
  )
  return(adjust_counts_whole)
}
