####  Expand a vector into matrix (columns as the original vector)
# Upstream-derived helper functions adapted from the ComBat-ref/ComBat-seq
# research implementation by Xiaoyu Zhang and collaborators. See
# LICENSE_NOTES.md for provisional licensing status.

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
    if(!is.null(gene.subset.n) && is.numeric(gene.subset.n) && length(gene.subset.n)==1){
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


####  Match quantiles (vectorized)
# keep_zero: zero values in the original counts don't change
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi,
                            keep_zero = TRUE, method = c("upper", "midp", "randomized")) {
  method <- match.arg(method)
  nr <- nrow(counts_sub)
  nc <- ncol(counts_sub)
  # Expand per-gene phi vectors into matrices matching counts_sub dimensions
  old_phi_mat <- vec2mat(old_phi, nc)
  new_phi_mat <- vec2mat(new_phi, nc)

  if (method == "midp") {
    lower <- pnbinom(counts_sub - 1, mu = old_mu, size = 1 / old_phi_mat)
    mass <- dnbinom(counts_sub, mu = old_mu, size = 1 / old_phi_mat)
    tmp_p <- lower + 0.5 * mass
    if (keep_zero) tmp_p[counts_sub == 0] <- 0
    new_counts_sub <- qnbinom(pmin(pmax(tmp_p, .Machine$double.xmin), 1 - .Machine$double.eps),
                              mu = new_mu, size = 1 / new_phi_mat)
  } else if (method == "randomized") {
    lower <- pnbinom(counts_sub - 1, mu = old_mu, size = 1 / old_phi_mat)
    tmp_p <- lower + stats::runif(length(counts_sub)) *
      dnbinom(counts_sub, mu = old_mu, size = 1 / old_phi_mat)
    new_counts_sub <- qnbinom(pmin(tmp_p, 1 - .Machine$double.eps), mu = new_mu, size = 1 / new_phi_mat)
  } else if(keep_zero) {
    # Compute CDF for counts-1 (shift for keep_zero mode)
    tmp_p <- pnbinom(counts_sub - 1, mu=old_mu, size=1/old_phi_mat)
    # Map through new distribution, adding 1 back
    new_counts_sub <- 1 + qnbinom(tmp_p, mu=new_mu, size=1/new_phi_mat)
    # Keep original values where count <= 1 or CDF is essentially 1 (outliers)
    keep_mask <- (counts_sub <= 1) | (abs(tmp_p - 1) < 1e-4)
    new_counts_sub[keep_mask] <- counts_sub[keep_mask]
  } else {
    tmp_p <- pnbinom(counts_sub, mu=old_mu, size=1/old_phi_mat)
    new_counts_sub <- qnbinom(tmp_p, mu=new_mu, size=1/new_phi_mat)
    # Keep original for outlier counts where CDF is essentially 1
    outlier_mask <- abs(tmp_p - 1) < 1e-6
    new_counts_sub[outlier_mask] <- counts_sub[outlier_mask]
  }

  storage.mode(new_counts_sub) <- "integer"
  attr(new_counts_sub, "tail_events") <- sum(tmp_p <= .Machine$double.xmin | tmp_p >= 1 - .Machine$double.eps)
  new_counts_sub
}
