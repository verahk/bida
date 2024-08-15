

#' Simulation: MCMC sampling of DAGs
#'
#' Sample a data set from a bayesian network and draw a sample of DAGs using
#' the MCMC schemes implemented in the `BiDAG`-package.
#'
#'
#' @param bn a bnlearn object
#' @param N sample size
#' @param r iteration (for seed number)
#' @param par
#'  a list that contains the following named elements:
#'  - `init, sample` algorithms for initiating search space and MCMC scheme, see [sample_dags].
#'  - `structure` local structure, see [optimize_partition].
#'  - `ess, edgepf, nlev` bdeu-score parameters, see [define scoreparameter]
#'  - `hardlimit` (integer) hard limit on maximum number of parents, see [sample_dags].
#'  - `regular` (logical) wheter local structures are forced to be regular, see [optimize_partition].
#'  - `N`
#'  - `r`
#' @param outdir (logical constant) if `is.null(outdir)` the results are written to file
#'
#' @return
#' If `is.null(outdir)`, a chain of DAGs (an BiDAG-object which class depend on `sample`).
#' Otherwise, the results is written to file and the functions returns `NULL`.


run_sample_dags <- function(par, verbose = FALSE) {

  bn <- readRDS(paste0("./data/", par$bnname, ".rds"))
  N <- par$N
  r <- par$r

  # draw data
  set.seed(N+r)
  data <- bida:::sample_data_from_bn(bn, N)
  nlev <- sapply(bn, function(x) dim(x$prob)[1])

  # define scorepars
  scoretype <- ifelse(par$struct == "none", "bdecat", par$struct)
  tmp <- c(par, nlev = list(nlev))
  lookup <- rlang::new_environment()
  scorepar  <- bida:::define_scoreparameters(data, scoretype, tmp, lookup = lookup)

  # run MCMC
  smpl <- bida:::sample_dags(scorepar, par$init, par$sample, hardlimit = par$hardlimit, verbose = verbose)

  return(list(smpl = smpl,
              lookup = lookup))
}

params_to_filename <- function(par) {
  tmp <- sprintf("%s_%s_%s_%s_ess%s_epf%s_reg%s_N%s_r%02.0f.rds",
                 par$bnname, par$init, par$struct, par$sample, par$ess, par$edgepf, par$regular*1, par$N, par$r)
  stopifnot(length(tmp) == 1)  # fails if any argument is NULL
  return(tmp)
}


# load libraries and aux functions ----
devtools::install_github("verahk/bida", ref = "dev")
library(doSNOW)
source("./inst/simulations/ldags/simulate_and_write_to_file.R", echo = T)

# specify simulation params ----
par <- list(bnname = "asia",
            init = "pcskel",
            struct = c("none", "ldag", "tree"),
            sample = "order",
            ess = 1,
            edgepf = 2,
            hardlimit = 5,
            regular = TRUE,
            N = 10**c(2:4),
            r = 1:10)
pargrid <- expand.grid(par, stringsAsFactors = FALSE)

outdir <- "./inst/simulations/ldags/MCMCchains/"  # directory for storing res
if (!dir.exists(outdir)) dir.create(outdir)

nClusters <- 6
simId <- format(Sys.time(), "%Y%m%d_%H%M%S")   # name of log file
export  <- c("run_sample_dags")                # objects to export with clusterExport()

# test ----
if (FALSE) {
  # test
  i <- 2
  filename <- "test.rds"
  simulate_and_write_to_file(simId,
                             outdir,
                             filename,
                             run_sample_dags,
                             par = pargrid[i, ],
                             verbose = TRUE)

  res <- readRDS(paste0(outdir, filename))
  ls.str(res$lookup)
  ls.str(res$lookup[[pargrid[i,]$struct]])
  remove.dir(outdir)

}

# run simulation ----
if (nClusters == 1) {
  for (i in seq_len(nrow(pargrid))) simulate_and_write_to_file(simId,
                                                          outdir,
                                                          params_to_filename(pargrid[i, ]),
                                                          run_sample_dags,
                                                          par = pargrid[i, ],
                                                          verbose = TRUE)
} else {
  cl <- makeCluster(nClusters, type = "SOCK", outfile = paste0(outdir, simId, ".out"))
  clusterExport(cl, export)
  registerDoSNOW(cl)
  foreach (i = seq_len(nrow(pargrid))) %dopar% simulate_and_write_to_file(simId,
                                                                    outdir,
                                                                    params_to_filename(pargrid[i, ]),
                                                                    run_sample_dags,
                                                                    par = pargrid[i, ],
                                                                    verbose = TRUE)
  stopCluster(cl)
}


