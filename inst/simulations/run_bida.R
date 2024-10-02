

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

on.exit(closeAllConnections())
library(doSNOW)
sapply(list.files("./inst/simulations/R", ".R", full.names = T),
       source, echo = T)

branch <- system("git branch --show-current", intern = TRUE)
indir <- paste0("./inst/simulations/", branch, "/MCMCchains/")
outdir <- paste0("./inst/simulations/" , branch, "/results/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
simId <- format(Sys.time(), "%Y%m%d_%H%M%S")   # name of log file

nClusters <- 4
doTest <- FALSE

sim_run <- function(indir, f, verbose = FALSE) {
  out <- list()

  # read MCMC-results from file
  filepath <- paste0(indir, f)
  res <- readRDS(filepath)
  par <- res$par
  lookup <- res$lookup
  MCMCchain <- res$MCMCchain

  N <- par$N
  r <- par$r

  bn <- readRDS(paste0("inst/data/", par$bnname, ".rds"))
  nlev <- vapply(bn, function(x) dim(x$prob)[1], integer(1))
  n    <- length(bn)
  dag <- bnlearn::amat(bn)
  dmat <- bida:::descendants(dag)
  pdo <- bida:::interv_probs_from_bn(bn, "bn")  # ground truth
  truetau <- matrix(vapply(pdo, bida:::JSD, numeric(1)), n, n)
  dindx <- diag(n) == 1

  # draw data
  set.seed(N+r)
  data <- bida:::sample_data_from_bn(bn, N)
  nlev <- sapply(bn, function(x) dim(x$prob)[1])

  # compute support over unique dags
  dags <- lapply(MCMCchain$traceadd$incidence, as.matrix)
  tmp <- unique(dags)
  support <- bida:::rowsum_fast(rep(1/length(dags), length(dags)), dags, tmp)
  dags <- tmp

  # estimate intervention distributions ----
  ## compute support over parent sets
  ps <- bida::parent_support_from_dags(dags)

  get_size <- function(x) {
    dims <- x$counts$dim
    c(parents = length(dims)-1,
      parts = ifelse(is.null(x$partition), prod(dims[-1]), length(x$partition)))
  }

  ## compute mse of point-estimates (mean) of intervention distribution
  set.seed(r)
  mse <- tau <- parents <- parts <- matrix(list(), n, n)
  cat("Start computing estimates for ", f, "\n")
  for (x in seq_len(n)) {
    cat(" Compute estimates for cause node x", x, "\n")
    for (y in seq_len(n)[-x]) {
      type <- ifelse(par$local_struct == "none", "cat", par$local_struct)
      pa   <- which(dag[, x] == 1) # true parents

      pairs <- list(
        unknown = bida::bida_pair(type, data, x, y,
                                  sets = ps$sets[[x]],
                                  support = ps$support[[x]],
                                  hyperpar = c(list(nlev = nlev), par),
                                  lookup = lookup),
        full = bida::bida_pair("cat", data, x, y,
                               sets = ps$sets[[x]],
                               support = ps$support[[x]],
                               hyperpar = c(list(nlev = nlev), par),
                               lookup = NULL),
        known = bida::bida_pair(type, data, x, y,
                        sets = matrix(pa, nrow = 1),
                        support = 1,
                        hyperpar = c(list(nlev = nlev), par),
                        lookup = lookup)
      )

      pdo_hat     <- lapply(pairs, bida::posterior_mean)
      mse[[x, y]]  <- vapply(pdo_hat, function(p) mean( (p-pdo[[x, y]])**2 ), numeric(1))
      tau[[x, y]]  <- vapply(pdo_hat, bida:::JSD, numeric(1))

      # number of conditioning variables - x + parents
      tmp <- lapply(pairs,
                    function(pair) do.call(rbind, lapply(pair$params, get_size)))
      tmp <- vapply(tmp, colMeans, numeric(2))

      parents[[x, y]] <- tmp[1, ]
      parts[[x, y]]   <- tmp[2, ]
    }
  }

  out$mse_pdo <- colMeans(do.call(rbind, mse[!dindx]))
  out$mse_tau <- colMeans((do.call(rbind, tau[!dindx]) - truetau[!dindx])**2)
  out$parents <- colMeans(do.call(rbind, parents[!dindx]))
  out$parts   <- colMeans(do.call(rbind, parts[!dindx]))

  # edge probs and ranking ---
  edgep <- Reduce("+", Map("*", dags, support))
  arp   <- Reduce("+", Map("*", lapply(dags, bida:::descendants), support))

  # compute average precision-recall
  compute_avgppv <- function(x, y) {
    indx <- order(x+runif(length(x))/1000, decreasing = TRUE)
    tp <- cumsum(y[indx])
    pp <- seq_along(x)
    mean((tp/pp)[y[indx] == 1])
  }
  out$rank <- c(arp = compute_avgppv(arp[!dindx], dmat[!dindx]),
                apply(do.call(rbind, tau[!dindx]), 2, compute_avgppv, y = dmat[!dindx]))

  topmat <- truetau > quantile(truetau[!dindx & truetau > 0], .8)
  out$ranktop <- c(arp = compute_avgppv(arp[!dindx], topmat[!dindx]),
                apply(do.call(rbind, tau[!dindx]), 2, compute_avgppv, y = topmat[!dindx]))
  rates <- rowsum(edgep[!dindx], dag[!dindx])/tabulate(dag[!dindx]+1, 2)
  out$edge <- c(n = sum(edgep[!dindx]),
                fpr = rates[1],
                tpr = rates[2],
                avgppv = compute_avgppv(edgep[!dindx], dag[!dindx]))

  out$par <- par
  out
}

# test ----
if (doTest) {
  # test
  filenames <- list.files(indir, ".rds")
  filename <- filenames[1]
  sim_and_write_to_file(dir_out = outdir,
                        filename = filename,
                        run = sim_run,
                        indir = indir,
                        f = filename)

  res <- readRDS(paste0(outdir, filename))
  str(res, max.level = 2)
  stop()
  file.remove(paste0(outdir, filename))

}

# profile ----
if (FALSE) {
  f <- "sachs_pcskel_ptree_order_ess1_epf2_N300_r16.rds"
  profvis::profvis(sim_run(indir = indir, f))

  # check how many (x, z) sets are in lookup
  ps <- bida:::parent_support_from_dags(dags)

}

# run ----
filenames <- list.files(indir, ".rds")
filenames <- filenames[!grepl("barley", filenames)]
filenames <- filenames[!grepl("alarm", filenames)]
filenames <- filenames[!grepl("insurance", filenames)]
filenames <- sample(filenames)


if (nClusters == 0) {

  for (f in filenames) {
    sim_and_write_to_file(dir_out = outdir, filename = f, run = sim_run, indir = indir, f = f)
  }

  stop()
} else {
  # set up cluster
  cl <- makeCluster(nClusters, type = "SOCK", outfile = paste0(outdir, simId, ".out"))
  export <- ls(pattern = "sim_")
  clusterExport(cl, export)
  registerDoSNOW(cl)

  # run
  foreach (f = filenames) %dopar% sim_and_write_to_file(dir_out = outdir,
                                                        filename = f,
                                                        run = sim_run,
                                                        indir = indir,
                                                        f = f)
  stopCluster(cl)
}



