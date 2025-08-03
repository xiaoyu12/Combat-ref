####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}


####  Monte Carlo integration functions
monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n){
  weights <- pos_res <- list()
  for(i in 1:nrow(dat)){
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]
    
    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }
    
    #LH <- sapply(1:G_sub, function(j){sum(log2(dnbinom(x, mu=m[j,], size=1/phi_sub[j])+1))})  
    LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/phi_sub[j]))})
    LH[is.nan(LH)]=0; 
    if(sum(LH)==0 | is.na(sum(LH))){
      pos_res[[i]] <- c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i]))
    }else{
      pos_res[[i]] <- c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
    }
    
    weights[[i]] <- as.matrix(LH/sum(LH))
  }
  pos_res <- do.call(rbind, pos_res)
  weights <- do.call(cbind, weights)
  res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"], weights=weights)	
  return(res)
} 


####  Match quantiles
# keep_zero: zero values in the original counts don't change
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi, keep_zero=TRUE){
  new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
  for(a in 1:nrow(counts_sub)){
    for(b in 1:ncol(counts_sub)){
      if(keep_zero) {
        if(counts_sub[a, b] <= 1){
          new_counts_sub[a,b] <- counts_sub[a, b]
        }else{
          tmp_p <- pnbinom(counts_sub[a, b]-1, mu=old_mu[a, b], size=1/old_phi[a])
          if(abs(tmp_p-1)<1e-4){
            new_counts_sub[a,b] <- counts_sub[a, b]  
            # for outlier count, if p==1, will return Inf values -> use original count instead
          }else{
            new_counts_sub[a,b] <- 1+qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
          }
        }
      } else {
        tmp_p = pnbinom(counts_sub[a, b], mu=old_mu[a, b], size=1/old_phi[a])
  
        if(abs(tmp_p-1) < 1e-6) {
          new_counts_sub[a,b] <- counts_sub[a, b]  
          # for outlier count, if p==1, will return Inf values -> use original count instead
        } else {
          new_counts_sub[a,b] <- qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
        }
      }
    }
  }
  return(new_counts_sub)
}



mapDisp <- function(old_mu, new_mu, old_phi, divider){
  new_phi <- matrix(NA, nrow=nrow(old_mu), ncol=ncol(old_mu))
  for(a in 1:nrow(old_mu)){
    for(b in 1:ncol(old_mu)){
      old_var <- old_mu[a, b] + old_mu[a, b]^2 * old_phi[a, b]
      new_var <- old_var / (divider[a, b]^2)
      new_phi[a, b] <- (new_var - new_mu[a, b]) / (new_mu[a, b]^2)
    }
  }
  return(new_phi)
}

#####

# Adjust Batch effect using NPMatch
source(file.path(script_dir, "NPM/R/NPmatch.R"))
source(file.path(script_dir, "NPM/R/normalize.log2CPM.R"))
library(limma)

reverse_log2cpm <- function(log2cpm_matrix, total_counts) {
  # Reverse the log2CPM transformation to get CPM
  cpm_matrix <- 2^(log2cpm_matrix) - 1
  
  total0 <- mean(total_counts)
  total <- ifelse(total0 < 1e6, total0, 1e6)

  # Use per-sample total counts to scale CPM back to raw counts
  count_matrix <- sweep(cpm_matrix, 2, total_counts, FUN = "*") / total
  
  # Convert to integer (rounding the values)
  count_matrix <- round(count_matrix)
  
  # Convert negative values to 0
  count_matrix[count_matrix < 0] <- 0
  
  # Return the reconstructed count matrix
  return(count_matrix)
}

NPM_adjust <- function(counts, batch, group, with.cpm.reverse = FALSE) {
  ## Intra-sample normalization of the raw data.
  ## We use the normalize.log2CPM.R function provided
  nX <- normalize.log2CPM(counts)

    ## Inter-sample normalization by quantile normalization
  nX <- limma::normalizeQuantiles(nX)

  ## Batch correction with NPmatch
  cX <- NPmatch(X=nX, y=group, dist.method="cor", sdtop=5000)

  ## reverse the log2CPM transformation
  if (with.cpm.reverse) {
    cX <- reverse_log2cpm(cX, colSums(counts))
  }

  return(cX)
}

NPM_lm_DEpipe <- function(cts, batch, group, alpha.unadj, alpha.fdr) {
  ## Intra-sample normalization of the raw data.
  ## We use the normalize.log2CPM.R function provided
  nX <- normalize.log2CPM(cts)
  
  ## Inter-sample normalization by quantile normalization
  nX <- limma::normalizeQuantiles(nX)
  
  ## Batch correction with NPmatch
  adj_counts <- NPmatch(X=nX, y=group, dist.method="cor", sdtop=5000)
  
  pval_seq <- apply(adj_counts, 1, function(x, group){
    x_norm <- scale(x, center=TRUE, scale=TRUE)
    fit3 <- lm(x_norm ~ as.factor(group))
    return(summary(fit3)$coefficients[2, 4])
  }, group=group)
  padj_seq <- p.adjust(pval_seq, method="fdr")
  
  de_called <- list(unadj=rownames(cts)[pval_seq < alpha.unadj], fdr=rownames(cts)[padj_seq < alpha.fdr], 
                    de_res=data.frame(PValue=pval_seq, FDR=padj_seq), design=model.matrix(~as.factor(group)))
  return(de_called)
}
