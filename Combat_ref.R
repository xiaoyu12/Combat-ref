# Correct batch effects for RNA-seq count data using a negative binomial distribution
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
#' @export
#' 

ComBat_ref <- function(counts, batch, group=NULL, full_mod=TRUE, 
                       gene.subset.n=NULL, genewise.disp=FALSE) 
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
  dge_obj <- DGEList(counts=counts)
  
  ## Prepare characteristics on batches
  n_batch <- nlevels(batch)  # number of batches
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch  
  n_batches <- sapply(batches_ind, length)
  #if(any(n_batches==1)){mean_only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  n_sample <- sum(n_batches)
  cat("Found",n_batch,'batches\n')
  
  # disp_common = rep(0, n_batch)
  # for (i in 1:n_batch) {
  #   disp_common[i] <- estimate_Batch_CommonDisp(counts[, batches_ind[[i]]], as.character(group[batches_ind[[i]]]))
  #   cat("Batch ",i, " dispersion = ", disp_common[i], "\n")
  # }
  
  #if(any(n_batches==1)){mean_only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  
  ## Make design matrix 
  # batch, use the first batch as the reference
  batchmod <- model.matrix(~batch)  # colnames: levels(batch)
  # covariate
  group <- as.factor(group)
  n_group <- nlevels(group)      # number of groups
  groups_ind <- lapply(1:n_group, function(i) {which(group==levels(group)[i])})  # list of samples in each group
  n_groups <- sapply(groups_ind, length)
  if(full_mod & nlevels(group)>1){
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~0+group)  # model.matrix(~0+group)
  }else{
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data=as.data.frame(t(counts)))
  }
  # combine
  design <- cbind(batchmod, mod)
  ## Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  #if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  ## Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    #if(ncol(design)<=(n_batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
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
    cat("Batch ",i, " dispersion = ", disp_common[i], "\n")
  }
  # Choose the batch with the smallest dispersion as the reference batch
  ref_batch = 1
  for (i in 2:n_batch) {
    if (disp_common[i] < disp_common[ref_batch]) {
      ref_batch = i
    }
  }
  cat("Reference batch: ", ref_batch, "\n")
  # Set reference batch as batch 1
  if(ref_batch != 1) {
    # swap disp_common
    tmp <- disp_common[1]
    disp_common[1] <- disp_common[ref_batch]
    disp_common[ref_batch] <- tmp
    for(i in 1:n_sample) {
      if(batch[i] == ref_batch) {
        batch[i] = 1
      } else if (batch[i] == 1) {
        batch[i] = ref_batch
      }
    }
  }
  # re-compute batches_ind
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch
  n_batches <- sapply(batches_ind, length)
  # Update the design matrix
  batchmod <- model.matrix(~batch)  # colnames: levels(batch)
  # combine
  design <- cbind(batchmod, mod)
  ## Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  #if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  
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
    increased_genes = which(gamma_hat[, kk-1] < -0.2)
    ncol = ncol(new_mu)
    ref_max = vec2mat(rowMaxs(mu_hat[, batches_ind[[1]]]), ncol)
    new_mu[increased_genes, ] = pmin(new_mu[increased_genes, ], ref_max[increased_genes, ])
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
  return(adjust_counts_whole)
}



