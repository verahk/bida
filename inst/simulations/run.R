

rm(list = ls())

# args ----
args <- as.list(commandArgs(trailingOnly = TRUE))
if (length(args) == 0) {
  args <- list(what = "MCMC",
            test_row = 0,
            nClusters = 4,
            iterStart = 1,
            iterSlutt = 30)
}

stopifnot(length(args) == 5)
names(args) <- c("what", "test_row", "nClusters", "iterStart", "iterSlutt")
args[-1] <- lapply(args[-1], as.numeric)
print(args)
vapply(args, class, character(1))

# load libraries ----
library(doSNOW)

sapply(list.files("./inst/simulations/R", ".R", full.names = T),
       source, echo = T)

# paths ----
branch  <- system("git branch --show-current", intern = TRUE)
subdir  <- switch(args$what, "MCMC" = "MCMCchains/", "results")
outdir <- paste0("./inst/simulations/", branch, "/", subdir)
subdir  <- switch(args$what, "MCMC" = "MCMCchains/")
indir <- paste0("./inst/simulations/", branch, "/", subdir)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# get params ----
pargrid <- sim_load_params("syntethic", args$what, args$iterStart, args$iterSlutt)
sim_run <- switch(args$what == "MCMC", get("sim_run_MCMC"))


# run ----
if (args$test_row > 0) {
  par <- pargrid[args$test_row, ]
  file <- params_to_filename(par)
  path <- paste0(outdir, file)
  sim_run(par, verbose = T)

  file.remove(path)
  sim_and_write_to_file(outdir, file, sim_run, par, verbose = TRUE)
  res <- readRDS(path)
  ls.str(res)
} else if (args$nClusters == 0) {
  for (r in seq_len(nrow(pargrid))) {
    cat(params_to_filename(pargrid[r, ]), "\n")
    sim_and_write_to_file(outdir,
                          params_to_filename(pargrid[r, ]),
                          sim_run,
                          par = pargrid[r, ], verbose = T)
  }
} else {

    # set up cluster
    simId <- format(Sys.time(), "%Y-%m-%d-%H-%m-%S")
    cl <- makeCluster(args$nClusters, type = "SOCK", outfile = paste0(outdir, simId, ".out"))
    export <- ls(pattern = "sim_")
    clusterExport(cl, export)
    registerDoSNOW(cl)


    # run
    foreach (r = seq_len(nrow(pargrid))) %dopar% sim_and_write_to_file(dir_out = outdir,
                                                          params_to_filename(pargrid[r, ]),
                                                          sim_run,
                                                          par = pargrid[r, ],
                                                          verbose = TRUE)
    stopCluster(cl)
}
