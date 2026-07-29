# Compute the GLM model parameters using fully MCMC method
library(rethinking)
# construct the data for the ulam input reshape the data
reshape_data <- function(ctr, nsamples, group, batch) {
  # the first sample
  d_reshaped <- data.frame(ctr[, 1])
  colnames(d_reshaped) <- c("count")
  d_reshaped$sample <- 1
  d_reshaped$group <- group[1]+1
  d_reshaped$batch <- batch[1]
  for (i in 2:nsamples) {
    d_i <- data.frame(ctr[, i])
    colnames(d_i) <- c("count")
    d_i$sample <- i
    d_i$group <- group[i]+1
    d_i$batch <- batch[i]
    d_reshaped <- rbind(d_reshaped, d_i)
  }
  d_reshaped$id = rep(1:nrow(ctr), nsamples)
  return (d_reshaped)
}

d_reshaped <- reshape_data(cts, N_total_sample, group, batch = batch)
num_genes <- max(d_reshaped$id)
dat <- list(
  c = d_reshaped$count,
  s = d_reshaped$sample,
  g = d_reshaped$group,
  b = d_reshaped$batch,
  id = d_reshaped$id,
  ng = max(d_reshaped$id),
  nb = max(d_reshaped$batch),
  nc = max(d_reshaped$group),
  lib_size = log(colSums(cts)/num_genes)
)

m <- ulam(
  alist (
    c ~ dgampois(lambda, phi),
    log(lambda) <- alpha[id] + beta[id, g] + gamma[id, b] + lib_size[s],
    phi <- phi_b[b],
    
    vector[ng]: alpha ~ normal(0, 1),
    matrix[ng, nb]: gamma ~ normal(0, 1),
    matrix[ng, nc]: beta ~ normal(0, 1),
    vector[nb]: phi_b ~ dexp(1)
  ), data = dat, chains = 4, cores = 4 
)

m2 <- ulam(
  alist (
    c ~ dgampois(lambda, phi),
    log(lambda) <- alpha[id] + beta[id, g] + (b-1)*gamma[id] + lib_size[s],
    phi <- phi_b[b],
    
    vector[ng]: alpha ~ normal(0, 1),
    vector[ng]: gamma ~ normal(0, 1),
    matrix[ng, nc]: beta ~ normal(0, 1),
    vector[nb]: phi_b ~ dexp(1)
  ), data = dat, chains = 4, cores = 4 
)

