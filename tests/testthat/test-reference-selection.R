reference_fixture <- function() {
  set.seed(431)
  genes <- 120L
  batch <- factor(rep(c("A", "B"), each = 8L))
  group <- factor(rep(rep(c("control", "treated"), each = 4L), 2L))
  baseline <- stats::rgamma(genes, shape = 3, rate = 0.08)
  group_effect <- rep(c(1, 8), each = genes / 2L)
  mu <- outer(baseline, rep(1, length(batch)))
  mu[, group == "treated"] <- mu[, group == "treated"] * group_effect
  counts <- matrix(stats::rnbinom(length(mu), mu = mu, size = 8), genes,
    dimnames = list(paste0("g", seq_len(genes)), paste0("s", seq_along(batch))))
  list(counts = counts, batch = batch, group = group)
}

test_that("balanced batches use locally estimable biological groups", {
  z <- reference_fixture()
  dge <- edgeR::normLibSizes(edgeR::DGEList(z$counts), method = "TMMwsp")
  selected <- combatrefql:::select_reference(dge, z$batch, z$group)
  expect_true(all(grepl("group", selected$batches$local_formula)))
  expect_true(all(grepl("group", selected$batches$local_columns)))
  expect_false(any(selected$batches$fallback))
  expect_true(all(selected$batches$residual_df > 0))

  idx <- z$batch == "A"
  intercept_fit <- edgeR::glmQLFit(dge[, idx], matrix(1, sum(idx), 1L),
    dispersion = NULL, robust = TRUE, legacy = FALSE)
  expect_lt(selected$batches$reference_score[selected$batches$batch == "A"],
    as.numeric(intercept_fit$dispersion))
})

test_that("unbalanced biology is adjusted where locally estimable", {
  z <- reference_fixture()
  z$group[z$batch == "B"] <- "control"
  z$group <- factor(z$group, levels = c("control", "treated"))
  dge <- edgeR::normLibSizes(edgeR::DGEList(z$counts), method = "TMMwsp")
  selected <- combatrefql:::select_reference(dge, z$batch, z$group)
  batch_a <- selected$batches[selected$batches$batch == "A", ]
  batch_b <- selected$batches[selected$batches$batch == "B", ]
  expect_match(batch_a$local_formula, "group")
  expect_false(batch_a$fallback)
  expect_true(batch_b$fallback)
  expect_match(batch_b$missing_levels, "group:treated")
  expect_match(batch_b$dropped_columns, "group")
})

test_that("missing levels and non-estimable terms are recorded", {
  group <- factor(rep("control", 4), levels = c("control", "treated"))
  covariates <- data.frame(
    varying = 1:4,
    constant = 2,
    category = factor(c("x", "x", "y", "y"), levels = c("x", "y", "z")),
    redundant = 1:4
  )
  design <- combatrefql:::build_local_reference_design(4, group, covariates)
  expect_match(design$missing_levels, "group:treated")
  expect_match(design$missing_levels, "category:z")
  expect_match(design$dropped_columns, "group")
  expect_match(design$dropped_columns, "constant")
  expect_match(design$dropped_columns, "redundant")
  expect_true(design$residual_df > 0)

  fallback <- combatrefql:::build_local_reference_design(3,
    factor(c("a", "b", "c")))
  expect_true(fallback$fallback)
  expect_equal(colnames(fallback$matrix), "(Intercept)")
  expect_match(fallback$fallback_reason, "locally estimable")
})

test_that("explicit and automatic reference selection are deterministic", {
  z <- reference_fixture()
  dge <- edgeR::normLibSizes(edgeR::DGEList(z$counts), method = "TMMwsp")
  first <- combatrefql:::select_reference(dge, z$batch, z$group)
  second <- combatrefql:::select_reference(dge, z$batch, z$group)
  expect_identical(first, second)
  expect_equal(sum(first$batches$second_best), 1L)
  expect_true(is.finite(first$batches$score_margin[1L]))
  expect_gte(first$batches$score_margin[1L], 0)

  other <- setdiff(levels(z$batch), first$reference)
  scored <- 0L
  explicit <- combatrefql:::select_reference(dge, z$batch, z$group,
    explicit = other, progress = function(...) scored <<- scored + 1L)
  expect_identical(explicit$reference, other)
  expect_true(explicit$batches$selected[explicit$batches$batch == other])
  expect_true(all(explicit$batches$selection_method == "explicit"))
  expect_equal(scored, 0L)
  expect_true(all(is.na(explicit$batches$reference_score)))
})
