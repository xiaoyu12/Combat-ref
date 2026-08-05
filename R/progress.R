new_reporter <- function(verbosity) {
  quiet <- verbosity == "quiet"
  completed <- 0L
  total <- 0L

  text <- function(message) if (!quiet) cli::cli_text(message)
  success <- function(message) if (!quiet) cli::cli_alert_success(message)

  list(
    start = function() text(sprintf("ComBat-refQL %s",
      as.character(utils::packageVersion("combatrefql")))),
    success = success,
    info = text,
    warning = text,
    reference = function(candidate) {
      if (verbosity == "verbose" && !interactive())
        text(sprintf("Scored reference candidate: %s", candidate))
    },
    start_mapping = function(steps) {
      total <<- steps
    },
    mapping = function(source, genes) {
      completed <<- completed + 1L
      if (verbosity == "verbose" &&
        (completed == 1L || completed == total ||
          completed %% max(1L, ceiling(total / 5)) == 0L)) {
        text(sprintf("Mapped batch %s: %s/%s gene chunks", source, completed, total))
      }
    },
    finish = function(seconds) text(sprintf("Completed in %.1f seconds", seconds))
  )
}
