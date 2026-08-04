abort_combatrefql <- function(class, message, ..., stage) {
  dots <- list(...)
  bullet <- names(dots) %in% c("*", "i", "!", "v")
  messages <- c("x" = message, "i" = sprintf("Stage: %s.", stage),
                unlist(dots[bullet], use.names = TRUE))
  args <- c(list(message = messages, class = c(class, "combatrefql_error"),
                 stage = stage, .envir = parent.frame()), dots[!bullet])
  do.call(cli::cli_abort, args)
}

input_error <- function(message, ..., stage = "input") {
  abort_combatrefql("combatrefql_input_error", message, ..., stage = stage)
}
