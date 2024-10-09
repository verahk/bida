

rm(list = ls())

# args ----
args <- as.list(commandArgs(trailingOnly = TRUE))
if (length(args) == 0) {
  args <- list("bida", 0, 6)
} else if (length(args) < 3) stop()
names(args) <- c("what", "test_row", "nClusters")
args[-1] <- lapply(args[-1], as.numeric)
print(args)
vapply(args, class, character(1))


# load libraries ----
library(doSNOW)

sapply(list.files("./inst/simulations/R", ".R", full.names = T),
       source, echo = T)

# paths ----
branch  <- system("git branch --show-current", intern = TRUE)
subdir  <- switch(args$what, "MCMC" = "MCMCchains", "bida" = "results")
outdir <- paste0("./inst/simulations/", branch, "/", subdir, "/")
subdir  <- switch(args$what, "bida" = "MCMCchains/")
indir <- paste0("./inst/simulations/", branch, "/", subdir, "/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# source file with specific sim_run fcuntion
filepath <- switch(args$what,
                   "MCMC" =  "./inst/simulations/sim_run_MCMC.R",
                   "bida" = "./inst/simulations/sim_run_bida.R")
source(filepath)

# profile ----
if (FALSE) {
  f <- function(par) {
    dag <- bida:::rand_dag(par$n, 8)

  }

  profvis::profvis(sim_load_bn(pargrid[1, ]))
}

if (args$test_row > 0) {
  # test ----
  par <- pargrid[args$test_row, ]

  file <- params_to_filename(par)
  path <- paste0(outdir, file)
  file.remove(path)

  cat("Run test for", path, "\n")
  sim_and_write_to_file(outdir, file, sim_run, par, verbose = TRUE)
  res <- readRDS(path)
  ls.str(res)
} else if (args$nClusters == 0) {
  # run: for  ----
  for (r in seq_len(nrow(pargrid))) {
    cat(params_to_filename(pargrid[r, ]), "\n")
    sim_and_write_to_file(outdir,
                          params_to_filename(pargrid[r, ]),
                          sim_run,
                          par = pargrid[r, ], verbose = T)
  }
} else {
  # run: foreach  ----
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
