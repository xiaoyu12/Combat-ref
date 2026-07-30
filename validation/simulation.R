set.seed(104729)
batch <- factor(rep(c("A", "B"), each = 6))
group <- factor(rep(rep(c("control", "treated"), each = 3), 2))
genes <- 200L
base <- rgamma(genes, 2, 0.08)
group_effect <- rep(c(1, 2), each = genes / 2)
batch_effect <- c(A = 1, B = 1.8)
mu <- outer(base, batch_effect[batch]) * outer(group_effect, as.numeric(group == "treated"), function(x, z) x^z)
counts <- matrix(rnbinom(length(mu), mu = mu, size = 10), genes,
                 dimnames = list(paste0("g", seq_len(genes)), paste0("s", seq_along(batch))))
fit <- combatrefql::combat_ref(counts, batch, group, reference_batch = "A", verbose = FALSE)
stopifnot(all(fit$adjusted_counts >= 0), all(fit$adjusted_counts == round(fit$adjusted_counts)))
batch_ratio <- function(x) median((rowMeans(x[, batch == "B"]) + 0.5) /
                                  (rowMeans(x[, batch == "A"]) + 0.5))
result <- data.frame(seed = 104729, genes = genes,
                     raw_batch_ratio = batch_ratio(counts),
                     adjusted_batch_ratio = batch_ratio(fit$adjusted_counts),
                     unsupported_genes = sum(!fit$gene_status$adjusted))
utils::write.table(result, "validation/simulation-summary.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)
