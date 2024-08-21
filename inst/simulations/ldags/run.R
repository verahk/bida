

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
devtools::install_github("verahk/bida",
                         ref = "sim_backdoor_estimates_with_local_structure")


library(doSNOW)

# specify simulation params ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nClusters <- 6
  bnnames <- c("asia", "child")
} else {
  nClusters <- as.numeric(args[1])
  bnnames <- args[-1]
  cat(sprintf("Commandline args: nclusters = %s, bnname(s) = %s\n",
              nClusters, paste(bnnames, collapse = ",")))
}

# additional params
par <- list(bnname = bnnames,
            init = "pcskel",
            struct = c("tree", "none"),
            sample = "order",
            ess = 1,
            edgepf = 2,
            hardlimit = 5,
            regular = FALSE,
            N = 10**c(2:4),
            r = 1:10)
pargrid <- expand.grid(par, stringsAsFactors = FALSE)


outdir <- "./inst/simulations/ldags/results/"  # directory for storing res
if (!dir.exists(outdir)) dir.create(outdir)
simId <- format(Sys.time(), "%Y%m%d_%H%M%S")   # name of log file
export  <- c("run")                # objects to export with clusterExport()




# define simulation routine -----
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

params_to_filename <- function(par) {
  tmp <- sprintf("%s_%s_%s_%s_ess%s_epf%s_reg%s_N%s_r%02.0f.rds",
                 par$bnname, par$init, par$struct, par$sample, par$ess, par$edgepf, par$regular*1, par$N, par$r)
  stopifnot(length(tmp) == 1)  # fails if any argument is NULL
  return(tmp)
}

run <- function(par, verbose = FALSE) {

  # load bn and compute ground truth ----
  bn <- readRDS(paste0("./data/", par$bnname, ".rds"))
  dag <- bnlearn::amat(bn)
  dmat <- bida:::descendants(dag)
  pdo <- bida:::interv_probs_from_bn(bn, "bn")  # ground truth

  n <- length(bn)
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
  tmp <- bida:::rowsum_fast(edgep[!dindx], dag[!dindx], 0:1)/tabulate(dag[!dindx]+1, 2)
  avgppv <- c(fpr = tmp[1],
              tpr = tmp[2],
              edgep  = compute_avgppv(edgep[!dindx], dag[!dindx]),
              arp    =  compute_avgppv(arp[!dindx], dmat[!dindx]))


  # estimate intervention distributions ----
  ## compute support over parent sets
  ps <- bida::parent_support_from_dags(dags)

  ## compute mse of point-estimates (mean) of intervention distribution
  set.seed(r)
  mse <- matrix(NA, n, n)
  for (x in seq_len(n)) {
    for (y in seq_len(n)[-x]) {
      type <- ifelse(par$struct == "none", "cat", par$struct)
      pair <- bida::bida_pair(type, data, x, y,
                               sets = ps$sets[[x]],
                               support = ps$support[[x]],
                               hyperpar = c(list(nlev = nlev), par),
                               lookup = scorepar$lookup)

      mse[x, y] <- mean( (pdo[[x, y]]-bida::posterior_mean(pair))**2 )
    }
  }

  mse <- sum(mse[!dindx])/(n*(n-1))

  list(res = c(avgppv, mse = mse),
       MCMCchain = MCMCchain)
}

# test ----
if (FALSE) {
  # test
  i <- 1
  filename <- "test.rds"
  simulate_and_write_to_file(simId,
                             outdir,
                             filename,
                             run,
                             par = pargrid[i, ],
                             verbose = TRUE)

  res <- readRDS(paste0(outdir, filename))
  str(res, max.level = 2)
  file.remove(paste0(outdir, filename))
}


# run simulation ----

if (nClusters == 1) {
  for (i in seq_len(nrow(pargrid))) simulate_and_write_to_file(simId,
                                                          outdir,
                                                          params_to_filename(pargrid[i, ]),
                                                          run,
                                                          par = pargrid[i, ],
                                                          verbose = TRUE)
} else {
  cl <- makeCluster(nClusters, type = "SOCK", outfile = paste0(outdir, simId, ".out"))
  clusterExport(cl, export)
  registerDoSNOW(cl)

  keepLooking <- TRUE
  row <- 0
  while (keepLooking) {
    row <- row+1
    keepLooking <- file.exists(paste0(outdir, params_to_filename(pargrid[row, ])))
  }
  foreach (i = seq(row, nrow(pargrid))) %dopar% simulate_and_write_to_file(simId,
                                                                    outdir,
                                                                    params_to_filename(pargrid[i, ]),
                                                                    run,
                                                                    par = pargrid[i, ],
                                                                    verbose = TRUE)
  stopCluster(cl)
}


