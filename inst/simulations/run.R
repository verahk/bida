

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

outdir <- "./inst/simulations/results/"  # directory for storing res
if (!dir.exists(outdir)) dir.create(outdir)
simId <- format(Sys.time(), "%Y%m%d_%H%M%S")   # name of log file


# specify simulation params ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nClusters <- 6
  bnnames <- c("asia", "sachs", "child")
} else {
  nClusters <- as.numeric(args[1])
  bnnames <- args[-1]
  cat(sprintf("Commandline args: nclusters = %s, bnname(s) = %s\n",
              nClusters, paste(bnnames, collapse = ",")))
}

# additional params
par <- list(init = c("pcskel"),
            local_struct = c("none", "ptree", "tree", "pcart"),
            sample = "order",
            ess = 1,
            edgepf = c(2, 4, 16),
            hardlimit = 4,
            N = c(300, 1000), # 3000, 10000),
            n = c(8),
            k = 2,
            complexity = c(.25, .75),
            r = 1:10)

pargrid <- expand.grid(par, stringsAsFactors = FALSE)
indx <- with(pargrid, local_struct == "none" & (edgepf > 2))
pargrid <- pargrid[!indx, ]



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
  tmp <- sprintf("n%s_k%s_csi%s_%s_%s_%s_ess%s_epf%s_N%s_r%02.0f.rds",
                 par$n, par$k, par$complexity*100, par$init, par$local_struct, par$sample, par$ess, par$edgepf, par$N, par$r)
  stopifnot(length(tmp) == 1)  # fails if any argument is NULL
  return(tmp)
}


sim_rand_partitions <- function(dag, nlev, splitprob, depth_ratio) {
  n <- ncol(dag)
  partitions <- vector("list", n)
  for (i in seq_len(n)) {
    pa <- which(dag[, i] == 1)
    if (length(pa) > 1) {
      partitions[[i]] <- sim_rand_partition(nlev[pa], splitprob, maxdepth = ceiling(depth_ratio*length(pa)))
    }
  }
  return(partitions)
}

sim_rand_partition <- function(nlev, splitprob, nextsplitprob = function(x) x, maxdepth = length(nlev)) {
  n <- length(nlev)
  
  
  # define routine for growing a tree
  grow_tree <- function(splitprob, vars, subset) {
    if (length(vars) <= (n-maxdepth) || runif(1) > splitprob) {
      return(list(subset = subset))
    }
    
    # sample a split variable
    pos <- sample.int(length(vars), 1)
    x <- vars[pos]
    
    # split the current subset by values of x
    xval <- (subset%/%stride[x])%%nlev[x]
    new_subsets <- unname(split(subset, xval))
    
    # grow a new tree for each value of the split variable
    list(var = x,
         branches = lapply(new_subsets,
                           function(y) grow_tree(nextsplitprob(splitprob), vars[-pos], y)))
  }
  
  
  # define routine for extracting subsets in each leaf
  unlist_tree <- function(tree) {
    if (is.null(tree$subset)) {
      unlist(lapply(tree$branches, unlist_tree), recursive = FALSE)
    } else {
      unname(tree["subset"])
    }
  }
  
  stride <- c(1, cumprod(nlev[-length(nlev)]))
  vars   <- seq_along(nlev)
  subset <- seq_len(prod(nlev))-1
  tree   <- grow_tree(splitprob, vars, subset)
  unlist_tree(tree)
}

sim_run <- function(par, verbose = FALSE) {

  n <- par$n
  k <- par$k
  N <- par$N
  r <- par$r

  nlev <- rep(k, n)
  
  # draw bn and compute ground truth
  set.seed(r)
  dag <- bida:::randDAG(n, d = 4)
  partitions <- sim_rand_partitions(dag, nlev, 1, par$complexity)
  bn <- bida:::rand_bn(dag, "cat", alpha = 1, nlev = nlev, partitions = partitions)
  dmat <- bida:::descendants(dag)
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
  mse <- matrix(NA, n, n)
  for (x in seq_len(n)) {
    for (y in seq_len(n)[-x]) {
      type <- ifelse(par$local_struct == "none", "cat", local_par$struct)
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
       par = par,
       MCMCchain = MCMCchain)
}

# test ----
if (FALSE) {
  # test
  i <- 2
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
                                                                    sim_run,
                                                                    par = pargrid[i, ],
                                                                    verbose = TRUE)
  stopCluster(cl)
}


