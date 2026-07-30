match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi,
                            keep_zero = TRUE,
                            method = c("upper", "midp", "randomized")) {
  method <- match.arg(method)
  nc <- ncol(counts_sub)
  old_phi <- matrix(rep(old_phi, nc), ncol = nc)
  new_phi <- matrix(rep(new_phi, nc), ncol = nc)

  lower <- stats::pnbinom(counts_sub - 1, mu = old_mu, size = 1 / old_phi)
  mass <- stats::dnbinom(counts_sub, mu = old_mu, size = 1 / old_phi)
  probability <- switch(method,
    midp = lower + 0.5 * mass,
    randomized = lower + stats::runif(length(counts_sub)) * mass,
    upper = lower + mass
  )
  if (keep_zero) probability[counts_sub == 0] <- 0
  probability <- pmin(pmax(probability, .Machine$double.xmin),
                      1 - .Machine$double.eps)
  mapped <- stats::qnbinom(probability, mu = new_mu, size = 1 / new_phi)
  if (keep_zero) mapped[counts_sub == 0] <- 0
  storage.mode(mapped) <- "double"
  valid <- is.finite(mapped) & mapped >= 0 & mapped == floor(mapped)
  attr(mapped, "tail_events") <- sum(
    probability <= .Machine$double.xmin |
      probability >= 1 - .Machine$double.eps
  )
  attr(mapped, "mapping_valid") <- apply(valid, 1L, all)
  mapped
}
