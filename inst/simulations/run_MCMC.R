

#' Simulation: MCMC sampling of DAGs
#'
#' Sample a data set from a bayesian network and draw a sample of DAGs using
#' the MCMC schemes implemented in the `BiDAG`-package.
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


# load libraries ----
rm(list = ls())
doTest <- FALSE

#devtools::install_github("verahk/bida", ref = "sim_slearn_with_lstruct", upgrade = "never")

library(doSNOW)
sapply(list.files("./inst/simulations/R", ".R", full.names = T),
       source, echo = T)

# paths ----
outdir <- "./inst/simulations/MCMCchains/"  # directory for storing res
if (!dir.exists(outdir)) dir.create(outdir)
simId <- format(Sys.time(), "%Y%m%d_%H%M%S")   # name of log file


# clusters ---
nClusters <- 4

# params ----
par <- list(init = c("pcskel"),
            local_struct = c("ptreereg"),
            sample = "order",
            ess = 1,
            edgepf = c(2),
            hardlimit = 4,
            N = c(300, 1000, 3000),
            n = c(10),
            k = c(2),
            complexity = c(0, .5, 1),
            r = 1:15)

pargrid <- expand.grid(par, stringsAsFactors = FALSE)
indx <- with(pargrid, local_struct == "none" & (edgepf > 2))
indx <- indx | with(pargrid, local_struct == "ptreereg" & (edgepf > 2))
indx <- indx | with(pargrid, local_struct == "pcart" & (k > 2))
#indx <- indx | with(pargrid, r %in% c(2, 3), n == 8, k == 2)
pargrid <- pargrid[!indx, ]


params_to_filename <- function(par) {
  tmp <- sprintf("n%s_k%s_csi%s_%s_%s_%s_ess%s_epf%s_N%s_r%02.0f.rds",
                 par$n, par$k, par$complexity*100, par$init, par$local_struct, par$sample, par$ess, par$edgepf, par$N, par$r)
  stopifnot(length(tmp) == 1)  # fails if any argument is NULL
  return(tmp)
}

sim_run <- function(par, verbose = FALSE) {

  n <- par$n
  k <- par$k
  N <- par$N
  r <- par$r

  nlev <- rep(k, n)

  # draw bn and compute ground truth
  set.seed(r)
  bn <- sim_rand_bn(n, nlev, par$complexity)$bn

  # draw data
  set.seed(N+r)
  data <- bida:::sample_data_from_bn(bn, N)

  # define scorepars
  lookup <- rlang::new_environment()
  if (par$local_struct == "ptreereg") {
    par$local_struct <- "ptree"
    par$regular <- TRUE
  }
  scorepar  <- bida:::define_scoreparameters(data,
                                             scoretype = "bdecat",
                                             par = c(par, nlev = list(nlev)),
                                             lookup = lookup)

  # run MCMC ----
  MCMCchain <- bida:::sample_dags(scorepar, par$init, par$sample, hardlimit = par$hardlimit, verbose = verbose)

  list(par = par,
       lookup = scorepar$lookup,
       MCMCchain = MCMCchain)
}

# test ----
if (doTest) {
  i <- 10
  filename <- params_to_filename(pargrid[i, ])
  file.remove(paste0(outdir, filename))
  sim_and_write_to_file(outdir,
                        params_to_filename(pargrid[i, ]),
                        sim_run,
                        par = pargrid[i, ],
                        verbose = TRUE)

  sim_and_write_to_file(outdir,
                        params_to_filename(pargrid[i, ]),
                        sim_run,
                        par = pargrid[i, ],
                        verbose = TRUE)

  results <- readRDS(paste0(outdir, filename))
  ls.str(results)
  attributes(results)
  stop()
}
# run ----
# set up cluster
cl <- makeCluster(nClusters, type = "SOCK", outfile = paste0(outdir, simId, ".out"))
export <- ls(pattern = "sim_")
clusterExport(cl, export)
registerDoSNOW(cl)

# find first row of pargrid for which no file exists
keepLooking <- TRUE
row <- 0
while (keepLooking) {
  row <- row+1
  keepLooking <- file.exists(paste0(outdir, params_to_filename(pargrid[row, ])))
}

# run
foreach (i = seq(row, nrow(pargrid))) %dopar% sim_and_write_to_file(outdir,
                                                                    params_to_filename(pargrid[i, ]),
                                                                    sim_run,
                                                                    par = pargrid[i, ],
                                                                    verbose = TRUE)
stopCluster(cl)


for (i in seq(row, nrow(pargrid))) sim_and_write_to_file(outdir,
                                                        params_to_filename(pargrid[i, ]),
                                                        sim_run,
                                                        par = pargrid[i, ],
                                                        verbose = TRUE)
