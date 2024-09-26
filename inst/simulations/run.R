

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


# load libraries ----
#devtools::install()

library(doSNOW)

outdir <- "./inst/simulations/MCMCchains/"  # directory for storing res
if (!dir.exists(outdir)) dir.create(outdir)
simId <- format(Sys.time(), "%Y%m%d_%H%M%S")   # name of log file


# specify simulation params ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nClusters <- 4
  bnnames <- c("asia", "sachs", "child")
} else {
  nClusters <- as.numeric(args[1])
  bnnames <- args[-1]
  cat(sprintf("Commandline args: nclusters = %s, bnname(s) = %s\n",
              nClusters, paste(bnnames, collapse = ",")))
}

# additional params
par <- list(init = c("pcskel"),
            local_struct = c("none", "ptree", "pcart"),
            sample = "order",
            ess = 1,
            edgepf = c(2),
            hardlimit = 4,
            N = c(300, 1000, 3000),
            n = c(10),
            k = c(2),
            complexity = c(0, .5, 1),
            r = 1:30)

pargrid <- expand.grid(par, stringsAsFactors = FALSE)
indx <- with(pargrid, local_struct == "none" & (edgepf > 2))
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
  bn <- sim_rand_bn(n, nlev, par$complexity)
  dmat <- bida:::descendants(bn)
  pdo <- bida:::interv_probs_from_bn(bn, "bn")  # ground truth

  # draw data
  set.seed(N+r)
  data <- bida:::sample_data_from_bn(bn, N)
  nlev <- sapply(bn, function(x) dim(x$prob)[1])

  # define scorepars
  lookup <- rlang::new_environment()
  scorepar  <- bida:::define_scoreparameters(data,
                                             scoretype = "bdecat",
                                             par = c(par, nlev = list(nlev)),
                                             lookup = lookup)

  # run MCMC ----
  MCMCchain <- bida:::sample_dags(scorepar, par$init, par$sample, hardlimit = par$hardlimit, verbose = verbose)

  # compute support over unique dags
  dags <- lapply(MCMCchain$traceadd$incidence, as.matrix)
  tmp <- unique(dags)
  support <- bida:::rowsum_fast(rep(1/length(dags), length(dags)), dags, tmp)
  dags <- tmp

  # edge probs  ---
  edgep <- Reduce("+", Map("*", dags, support))
  arp   <- Reduce("+", Map("*", lapply(dags, bida:::descendants), support))

  # compute precision-recall of edges
  compute_avgppv <- function(x, y) {
    indx <- order(x+runif(length(x))/1000, decreasing = TRUE)
    tp <- cumsum(y[indx])
    pp <- seq_along(x)
    mean((tp/pp)[y[indx] == 1])
  }

  dindx <- diag(n) == 1
  tmp <- bida:::rowsum_fast(edgep[!dindx], dag[!dindx], c(0, 1))/tabulate(dag[!dindx]+1, 2)
  avgppv <- c(fpr = tmp[1],
              tpr = tmp[2],
              edgep = compute_avgppv(edgep[!dindx], dag[!dindx]),
              arp   = compute_avgppv(arp[!dindx], dmat[!dindx]))


  # estimate intervention distributions ----
  ## compute support over parent sets
  ps <- bida::parent_support_from_dags(dags)

  ## compute mse of point-estimates (mean) of intervention distribution
  set.seed(r)
  mse_known <- mse_unknown <- matrix(NA, n, n)
  for (x in seq_len(n)) {
    for (y in seq_len(n)[-x]) {
      type <- ifelse(par$local_struct == "none", "cat", par$local_struct)
      pair <- bida::bida_pair(type, data, x, y,
                               sets = ps$sets[[x]],
                               support = ps$support[[x]],
                               hyperpar = c(list(nlev = nlev), par),
                               lookup = NULL)

      mse_unknown[x, y] <- mean( (pdo[[x, y]]-bida::posterior_mean(pair))**2 )

      # estimate intervention prob given true parents
      pa   <- which(dag[x, ] == 1) # true parents
      pair <- bida::bida_pair(type, data, x, y,
                              sets = matrix(pa, nrow = 1),
                              support = 1,
                              hyperpar = c(list(nlev = nlev), par),
                              lookup = NULL)
      mse_known[x, y] <- mean( (pdo[[x, y]]-bida::posterior_mean(pair))**2 )

    }
  }

  mse <- colSums(cbind(unknown = mse_unknown[!dindx], known = mse_known[!dindx]))/(n*(n-1))

  # repeat with known DAG

  list(res = c(avgppv, mse = mse),
       par = par,
       MCMCchain = MCMCchain)
}

# test ----
if (FALSE) {
  # test
  i <- 3
  pargrid[3, ]
  filename <- params_to_filename(pargrid[i, ])
  simulate_and_write_to_file(simId,
                             outdir,
                             filename,
                             sim_run,
                             par = pargrid[i, ],
                             verbose = TRUE)

  res <- readRDS(paste0(outdir, filename))
  str(res, max.level = 2)
  file.remove(paste0(outdir, filename))
}
# run simulation ----
export  <- ls(pattern = "sim_")               # objects to export with clusterExport()

if (nClusters == 1) {
  for (i in seq_len(nrow(pargrid))) simulate_and_write_to_file(simId,
                                                          outdir,
                                                          params_to_filename(pargrid[i, ]),
                                                          sim_run,
                                                          par = pargrid[i, ],
                                                          verbose = TRUE)
} else {

}


