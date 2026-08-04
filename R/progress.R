new_reporter <- function(verbosity) {
  quiet <- verbosity == "quiet"
  progress_id <- NULL
  completed <- 0L
  total <- 0L

  text <- function(message) if (!quiet) cli::cli_text(message)
  success <- function(message) if (!quiet) cli::cli_alert_success(message)

  list(
    start = function() text("ComBat-refQL 0.3.0"),
    success = success,
    info = text,
    warning = text,
    reference = function(candidate) {
      if (verbosity == "verbose" && !interactive())
        text(sprintf("Scored reference candidate: %s", candidate))
    },
    start_mapping = function(steps) {
      total <<- steps
      if (!quiet && interactive() && steps > 0L)
        progress_id <<- cli::cli_progress_bar("Mapping batches", total = steps)
    },
    mapping = function(source, genes) {
      completed <<- completed + 1L
      if (!is.null(progress_id)) {
        cli::cli_progress_update(id = progress_id, set = completed)
      } else if (verbosity == "verbose" &&
        (completed == 1L || completed == total ||
          completed %% max(1L, ceiling(total / 5)) == 0L)) {
        text(sprintf("Mapped batch %s: %s/%s gene chunks", source, completed, total))
      }
    },
    finish_mapping = function() {
      if (!is.null(progress_id)) cli::cli_progress_done(progress_id)
    },
    finish = function(seconds) text(sprintf("Completed in %.1f seconds", seconds))
  )
}
