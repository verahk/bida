
simulate_and_write_to_file <- function(id, outdir, filename, run, ...) {
  if (is.null(outdir)) return(run(...))

  filepath <- paste0(outdir, filename)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  if (file.exists(filepath)) return(NULL)

  tic <- Sys.time()
  res <- run(...)
  attr(res, "simid") <- id

  cat("\nSave results to: ", filepath)
  saveRDS(res, filepath)
  cat("\nRuntime\n")
  print(Sys.time()-tic)

  return(NULL)
}
