

#' Run simulation and write result to file
#'
#' Wrapper function around some simulation `run` that only run if the results are
#' not already present in the directory defined by `dir_out`
#'
#' @param id (character) id of simulation
#' @param dir_out (character) path to directory where results should be stored
#' @param filename (character) name of path
#' @param run (function) a function that runs a simulation
#' @param ... additional arguments, sent to `run`.
#'
#' @return
#' If `is.null(dir_out):` the output of `run`. Otherwise the output is written to
#' file and the function returns `NULL`.
#'
#' @keywors Internal
sim_and_write_to_file <- function(dir_out, filename, run, ...) {
  if (is.null(dir_out)) return(run(...))

  filepath <-
  if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
  if (file.exists(paste0(dir_out, filename))) return(NULL)

  tic <- Sys.time()
  res <- run(...)

  cat("\nSave results to: ", filepath)
  saveRDS(res, filepath)
  cat("\nRuntime\n")
  print(Sys.time()-tic)

  return(NULL)
}

