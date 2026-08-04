# Optional reproducible comparison; benchmark-only dependency: sva.
if (!requireNamespace("sva", quietly = TRUE)) stop("Install sva to run this benchmark.")
set.seed(20260804)
genes <- 2000L; samples <- 24L; batches <- 3L
batch <- factor(rep(LETTERS[seq_len(batches)], each = samples / batches))
group <- factor(rep(c("control", "treated"), length.out = samples))
mu <- outer(stats::rgamma(genes, 2, 0.05), rep(c(1, 1.4, 0.8), each = samples / batches))
counts <- matrix(stats::rnbinom(length(mu), mu = mu, size = 8), genes,
  dimnames = list(sprintf("gene_%05d", seq_len(genes)), sprintf("sample_%02d", seq_len(samples))))

ql_time <- system.time(ql <- combatrefql::combat_ref_ql(
  counts, batch, group, reference = "A", verbosity = "quiet"))[["elapsed"]]
sva_time <- system.time(sva_fit <- sva::ComBat_seq(counts, batch, group = group))[["elapsed"]]
result <- data.frame(
  method = c("combatrefql", "sva_ComBat_seq"), elapsed_seconds = c(ql_time, sva_time),
  approximate_result_mb = c(as.numeric(object.size(ql)), as.numeric(object.size(sva_fit))) / 1024^2,
  genes = genes, samples = samples, batches = batches,
  R = as.character(getRversion()), edgeR = as.character(packageVersion("edgeR")),
  sva = as.character(packageVersion("sva")), stringsAsFactors = FALSE)
print(result, row.names = FALSE)

