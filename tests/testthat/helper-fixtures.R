combatrefql_test_fixtures <- function() {
  list(
    numeric_balanced = list(
      counts = structure(
        c(
          10, 0, 1, 6, 3, 60,
          12, 0, 2, 7, 4, 65,
          35, 0, 0, 12, 5, 150,
          38, 0, 0, 14, 6, 160,
          8, 0, 3, 7, 2, 55,
          9, 0, 2, 8, 3, 58,
          31, 0, 0, 15, 7, 140,
          34, 0, 0, 16, 8, 155
        ),
        dim = c(6L, 8L),
        dimnames = list(
          c("gene_stable", "gene_all_zero", "gene_zero_batch2", "gene_low", "gene_small", "gene_large"),
          paste0("sample", 1:8)
        )
      ),
      batch = c(1, 1, 2, 2, 1, 1, 2, 2),
      group = c("ctrl", "ctrl", "ctrl", "ctrl", "trt", "trt", "trt", "trt"),
      expected_reference_batch = "1"
    ),
    character_unbalanced = list(
      counts = structure(
        c(
          5, 1, 0, 20, 2, 90,
          7, 2, 0, 24, 3, 110,
          6, 1, 0, 22, 1, 100,
          20, 0, 4, 50, 8, 190,
          24, 0, 5, 55, 7, 210,
          23, 0, 6, 53, 9, 205,
          25, 0, 4, 58, 8, 220
        ),
        dim = c(6L, 7L),
        dimnames = list(
          c("gene_a", "gene_b", "gene_zero_batch_A", "gene_d", "gene_low", "gene_big"),
          paste0("char_sample", 1:7)
        )
      ),
      batch = c("batch_A", "batch_A", "batch_A", "batch_B", "batch_B", "batch_B", "batch_B"),
      group = c("ctrl", "ctrl", "trt", "ctrl", "ctrl", "trt", "trt"),
      expected_reference_batch = "batch_A"
    )
  )
}
