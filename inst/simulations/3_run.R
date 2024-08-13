# prep ----
library(foreach)
library(bida)

files <- list.files("inst/simulations/R", full.names = T)
sapply(files, source)


# param settings
outpath <- "../simres/"
ncores <- 4
par <- expand.grid(N = c(300, 1000, 3000),
                   r = 1:30,
                   method = c("PC05", "MCMCPC", "unadjusted", "KNOWN"),
                   bnname = c("child", "alarm", "win95pts", "sachs", "asia", "water2", "hepar2", "hailfinder", "insurance"))



# check which settings are already runned
filenames <- sprintf("%s_N%0.0f_r%02.0f_%s_sqerr.rds", par$bnname, par$N, par$r, par$method)
indx <- file.exists(paste0(outpath, filenames))
stopifnot(any(!indx))
par <- par[!indx, ]

# function
sim_bida <- function(bnname, N, r, method, outpath = NULL, verbose = T) {

  tic <- Sys.time()
  tag <- sprintf("%s_N%0.0f_r%02.0f", bnname, N, r)
  if (verbose) cat("\n Start simulation", format(Sys.time(), "%a %b %d %H:%M:%S"), "method: ", method, "tag:", tag)


  # import bn and sample data ----
  bn <- readRDS(paste0("../bida/inst/data/", bnname, ".rds"))
  nlev <- vapply(bn, function(node) dim(node$prob)[1], integer(1))
  n <- length(bn)
  dindx <- diag(n) == 1
  set.seed(r+N)
  data <- bida:::sample_data_from_bn(bn, N)

  # compute ground truth ----
  ace_funs <- list(jsd = bida:::avg_jsd_array,
                   abs = bida:::avg_abs_array)
  truth <- compute_ground_truth(bn, ace_funs = ace_funs)

  # compute support over adjustment sets ----
  max_nlev <- max(nlev)
  checksize <- function(x, y, z) length(z) <= 10/log(max_nlev, 2) || sum(log(nlev[z]), 2) < 10
  tic <- c(tic, data = Sys.time())

  if (method == "KNOWN") {
    adjsets <- c("pa", "pa_min", "o", "o_min")
    if (n > 50) adjsets <- adjsets[-3]  # exclude o-set for large networks

    dags <- list(bnlearn::amat(bn))
    ps <- c(bida:::adjsets_support_from_dags(adjsets[-1], dags, checksize = checksize),
            pa = list(bida:::parent_support_from_dags(dags, checksize = checksize)))

  } else if (method == "MCMCPC") {

    filepath <- paste0("../adjsets/", tag, ".rds")
    if (!file.exists(filepath)) return(NULL)
    ps <- readRDS(filepath)

  } else if (method == "PC05") {

    fit <- pcalg::pc(list(dm = data, nlev = truth$nlev, adaptDF = FALSE), u2pd = "retry",
                     indepTest = pcalg::disCItest, alpha = .05, p = ncol(data))
    cpdag <- as(fit@graph, "matrix")
    ps <- bida:::adjsets_support_from_cpdag(c("pa", "o"), cpdag)

  } else if (method == "unadjusted") {

    ps <- list(cond = list(sets =  matrix(list(matrix(NA)), n, n),
                           support = matrix(list(1), n, n)),
               marg = list(sets =  matrix(lapply(seq_len(n), function(i) matrix(i)), n, n, byrow = T),
                           support = matrix(list(1), n, n)))

  }

  tic <- c(tic, slearn = Sys.time())
  tic_prep <- diff(tic)

  # compute posterior intervention distribution and estimate causal params ----
  # prep structures for storing res
  adjsets <- names(ps)

  dims <- list(descendant = c(FALSE, TRUE),
               variable = c("pdo", names(ace_funs)),
               adjset = adjsets)
  sqerr <- array(dim = lengths(dims), dimnames = dims)


  thresholds <- apply(truth$ace, 2, function(x) c(pos = 0, top = quantile(x, .8)))
  n1 <- ceiling(nrow(truth$ace)*c(1, .2))
  npairs <- n*(n-1)
  dims <- list(th = rownames(thresholds),
               variable = names(ace_funs),
               method = c("value", "rank", "arp"),
               adjset = adjsets)

  auc <- array(n1/npairs, dim = lengths(dims), dimnames = dims)

  toc <- setNames(vector("list", length(adjsets)), adjsets)

  # precompute common things
  truth$pdo_unlisted <- data.frame(descendant = rep.int(truth$dmat[!dindx], lengths(truth$pdo[!dindx])),
                                   value = unlist(truth$pdo[!dindx]))


  for (adjset in adjsets) {
    tic <- Sys.time()
    # compute posterior distribution for each cause-effect pair
    bp <- bida:::bida(ps[[adjset]], data, "categorical", list(nlev = nlev, ess = 1))
    tic <- c(tic, "fit" = Sys.time())


    # compute MSE of interventional CPTs
    pdo <- bida::posterior_mean(bp)[!dindx]
    tmp <- unlist(pdo)-truth$pdo_unlisted$value
    sqerr[, "pdo", adjset] <- bida:::rowsum_fast(tmp**2, truth$pdo_unlisted$descendant, c(0, 1)) # return grouped sum also if gr has only one unique value
    tic <- c(tic, "pdo" = Sys.time())

    # posterior mean and rank ----
    if (method %in% c("KNOWN", "MCMCPC")) {
      # enumerate cause-effect pairs with support of non-zero effects
      whichPos <- which(!vapply(bp, function(x) is.null(x) || x$zerosupp == 1, logical(1)))

      if (length(whichPos) == 0) {
        # all predicted effects are zero
        sqerr[1, -1, adjset] <- 0
        sqerr[2, -1, adjset] <- colSums( truth$ace**2 )

      } else {

        # compute mean squared error and auc-PR
        mapTruePos <- match(whichPos, truth$desc, nomatch = 0)

        # sample causal effects
        set.seed(007)
        samplesize <- 10**3
        smpls <- lapply(bp[whichPos], bida:::posterior_sample, n = samplesize, ace_funs = ace_funs)
        tic <- c(tic, "smpl" = Sys.time())

        # compute posterior mean value
        value <- t(vapply(smpls, colSums, numeric(length(ace_funs)))/samplesize)
        sqerr[, names(ace_funs), adjset] <- sqerr_tau_from_pos_pred(value, truth$ace, mapTruePos)

        if (any(mapTruePos > 0)) {
          # if no correct true positives, all positives are given lowest rank and AUCPR equals default n1/n

          # compute posterior mean rank
          #matrix(unlist(list(matrix(1:6, 3, 2), matrix(7:12, 3, 2))), ncol = 2)
          tmp  <- apply(matrix(unlist(smpls), ncol = length(smpls)), 1, rank, ties.method = "random")
          gr    <- rep(names(ace_funs), each = samplesize)
          rank  <- t(rowsum(t(tmp), gr, reorder = F))/samplesize

          # compute ancester relation prob
          arp   <- 1-sapply(bp[whichPos], "[[", "zerosupp")

          true <- matrix(0, nrow = length(mapTruePos), ncol = length(ace_funs))
          true[mapTruePos>0, ] <- truth$ace[mapTruePos]
          colnames(true) <- names(ace_funs)

          for (ace in names(ace_funs)) {
            for (i in seq_len(nrow(thresholds))) {
              y <- true[, ace] > thresholds[i, ace]
              auc[i, ace, "value", adjset] <- aucpr(value[, ace], y, n =  npairs, n1 = n1[i])
              auc[i, ace, "rank", adjset]  <- aucpr(rank[, ace], y, n =  npairs, n1 = n1[i])
              auc[i, ace, "arp", adjset]   <- aucpr(arp, y, n =  npairs, n1 = n1[i])
            }
          }
        }
        tic <- c(tic, "tau" = Sys.time())
      }

    } else if (adjset == "marg") {
      # all predicted effects are zero
      sqerr[1, -1, adjset] <- 0
      sqerr[2, -1, adjset] <- colSums( truth$ace**2 )
      auc[,, c("rank", "arp"), adjset] <- NA
      tic <- c(tic, "tau" = Sys.time())
    } else {


      # compute mean of intervention distribution for each adjustment set
      ida_mean <- function(x, y){
        pdo <- lapply(bp[[x, y]]$params, bida:::posterior_mean.backdoor_params_cat, ess = 1, kx = nlev[x])
        vapply(ace_funs, function(f) mean(vapply(pdo, f, numeric(1))), numeric(1))
      }

      value <- t(mapply(ida_mean,
                      x = .row(c(n, n))[!dindx],
                      y = .col(c(n, n))[!dindx]))

      true <- matrix(0, nrow = n**2, ncol = length(ace_funs))
      true[truth$desc, ] <- truth$ace
      true <- true[!dindx, ]
      colnames(true) <- names(ace_funs)

      sqerr[, -1, adjset] <- rowsum( (value-true)**2, truth$dmat[!dindx])

      for (ace in names(ace_funs)) {
        for (i in seq_len(nrow(thresholds))) {
          y <- true[, ace] > thresholds[i, ace]
          auc[i, ace, "value", adjset] <- aucpr(value[, ace], y, n =  npairs, n1 = n1[i])
          auc[i, ace, "rank", adjset]  <- NA
          auc[i, ace, "arp", adjset]   <- NA
        }
      }
      tic <- c(tic, "tau" = Sys.time())
    }


    toc[[adjset]] <- c(tic_prep, diff(tic))
  }

  if (!is.null(outpath)) {
    saveRDS(sqerr, paste0(outpath, tag, "_", method, "_sqerr.rds"))
    saveRDS(auc, paste0(outpath, tag, "_", method, "_auc.rds"))
    saveRDS(toc, paste0(outpath, tag, "_", method, "_toc.rds"))
  } else {
    return(list(sqerr, auc, toc))
  }
}

# test ----
if (FALSE) {
  bnname <- "hepar2"
  N <- 3000
  r <- 1
  method <- "KNOWN"
  method <- "unadjusted"
  method <- "PC05"

  outpath <- "../simres/"
  files <- sapply(c("sqerr", "auc", "toc"),
                  function(x) sprintf("%s%s_N%0.0f_r%02.0f_%s_%s.rds", outpath, bnname, N, r, method, x))
  file.remove(files)
  sim_bida(bnname, N, r, method, outpath = NULL)

  foreach (f = files) %do% {
    test <- readRDS(f)
    names(test)
    lapply(test, lengths)
    test
  }
}


# run ----
cl <- parallel::makeCluster(ncores, type="SOCK", outfile = "")
doSNOW::registerDoSNOW(cl)

foreach (bnname = par$bnname,
         N = par$N,
         r = par$r,
         method = par$method,
         .packages = c("Matrix", "bida"),
         .verbose = T) %dopar% sim_bida(bnname, N, r, method, outpath = outpath, verbose = T)

on.exit(parallel::stopCluster(cl))
