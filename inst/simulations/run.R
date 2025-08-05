

rm(list = ls())

# args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args <- c(what = "MCMC",
            doTest = 1,
            nClusters = 4)
}

# load libraries ----
library(doSNOW)

sapply(list.files("./inst/simulations/R", ".R", full.names = T),
       source, echo = T)

# paths ----
branch  <- system("git branch --show-current", intern = TRUE)
subdir  <- switch(args[1], "run_MCMC.R" = "MCMCchains", "results")
outdir <- paste0("./inst/simulations/", branch, "/", subdir)
subdir  <- switch(args[1], "run_MCMC.R" = "MCMCchains")
indir <- paste0("./inst/simulations/", branch, "/", subdir)

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# get params ----
pargrid <- sim_load_params("syntethic", args[1])
sim_run <- switch(args[1] == "MCMC", get("sim_run_MCMC"))

# run ----
if (args[[2]] > 0) {
  par <- pargrid[args[2], ]
  file <- params_to_filename(par)
  path <- paste0(outdir, file)
  sim_run(par, verbose = T)

  file.remove(path)
  sim_and_write_to_file(outdir, file, par, verbose = TRUE)
  res <- readRDS(path)
  ls.str(res)
} else if (nClusters == 0) {
  for (r in seq_len(nrow(pargrid))) {
    sim_and_write_to_file(outdir,
                          params_to_filename(pargrid[r, ]),
                          sim_run,
                          par = pargrid[r, ])
  }
} else {

    # set up cluster
    simId <- format(Sys.time(), "%a %b %d %X %Y")
    cl <- makeCluster(nClusters, type = "SOCK", outfile = paste0(outdir, simId, ".out"))
    export <- ls(pattern = "sim_")
    clusterExport(cl, export)
    registerDoSNOW(cl)
    onExit( stopCluster(cl))

    # run
    foreach (r = seq_len(nrow(pargrid))) %dopar% sim_and_write_to_file(dir_out = outdir,
                                                          params_to_filename(pargrid[r, ]),
                                                          sim_run,
                                                          par = pargrid[r, ])
    stopCluster(cl)
}
