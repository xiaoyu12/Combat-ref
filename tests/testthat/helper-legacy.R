legacy_combat_ref <- function(...) {
  out <- combatrefql:::ComBat_ref(...)
  attr(out, "reference_batch") <- NULL
  out
}
