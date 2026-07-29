build_combat_design <- function(batch, group = NULL) {
  batch <- factor(batch)
  batchmod <- stats::model.matrix(~ -1 + batch)

  if (is.null(group)) {
    mod <- stats::model.matrix(~ 1, data = data.frame(batch = batch))
  } else {
    group <- factor(group)
    mod <- stats::model.matrix(~ group)
  }

  design <- cbind(batchmod, mod[, -1, drop = FALSE])
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  design
}
