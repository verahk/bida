

# WHAT: compute support over unique adjustment size
# WHY: intermediate step in bida-simulations, with long runtime.
# HOW:
# For each network, samplesize and iteration do:
# 1. read in chain of sampled DAGs from file,
# 2. define function for checking size of adjset for current bn
# 3. apply bida::parent_support_from_dags for parent sets (identified locally)
# 4. apply bida::adjset_support_from_dags for o-set and minimal sets
# 5. combine output in one list and write to file


library(foreach)
library(bida)

par <- expand.grid(#bnname = c("sachs", "asia", "water2", "hepar2", "hailfinder", "insurance"),
                   bnname =  c("child",  "alarm", "win95pts"),
                   N = c(300, 1000, 3000),
                   r = 1:30)
maxconf <- 2**10
ncores <- 4
outpath <- "../adjsets/"

# check which sim-setting-results that are not already written to file
files <- paste0(outpath, sprintf("%s_N%s_r%02.0f", par$bnname, par$N, par$r), ".rds")
indx  <- file.exists(files)
par <- par[!indx, ]


adjsets_to_file <- function(bnname, N, r, maxconf, outpath) {

  tag <- sprintf("%s_N%0.0f_r%02.0f", bnname, N, r)
  cat("\n", format(Sys.time(), "%a %b %d %H:%M:%S"), "tag:", tag)

  # import MCMC-sample of DAGs
  tmppath <- paste0("../partitionMCMC/", tag, ".rds")
  if (!file.exists(tmppath)) return(NULL)
  partfit <- readRDS(tmppath)

  # collect chain of sampled DAGs
  dags    <- lapply(partfit$traceadd$incidence, as.matrix)
  support <- rep(1, length(dags))/length(dags)

  # max-size function
  bn <- readRDS(paste0("./inst/data/", bnname, ".rds"))
  nlev <- vapply(bn, function(node) dim(node$prob)[1], integer(1))
  checksize <- function(x, y, z) length(z) <= 1 || sum(log(nlev[z], 2)) <= log(maxconf, 2)

  # compute support over unique adjustment size
  ps <- list()

  ## parent sets, locally
  a <- "pa"
  tic <- Sys.time()
  ps[[a]] <- bida:::parent_support_from_dags(dags, support, checksize = checksize)
  attr(ps[[a]], "toc") <- Sys.time()-tic

  ## remaining adjustment sets
  if (bnname %in% c("alarm", "child", "win95pts")) {
    # identify adjset one-by-one, for comparing running times
    for (a in c("o", "o_min", "pa_min")) {
      if (a == "o" && length(nlev) > 50) next
      tic <- Sys.time()
      ps[[a]] <- bida:::adjsets_support_from_dags(a, dags, support, checksize = checksize)[[a]]
      attr(ps[[a]], "toc") <- Sys.time()-tic
      gc()
    }
  } else {
    tic <- Sys.time()
    if (length(nlev) > 50) {
      adjsets <- c("o_min", "pa_min")
    } else {
      adjsets <- c("o", "o_min", "pa_min")
    }
    ps[adjsets] <- bida:::adjsets_support_from_dags(adjsets, dags, support, checksize = checksize)
    attr(ps, "toc") <- Sys.time()-tic
    gc()
  }

  attr(ps, "session_info") <- sessionInfo(package = "bida")
  filepath <- paste0(outpath, tag, ".rds")
  cat("\nwrite to file: ", filepath)
  saveRDS(ps, filepath)
}

# test ----
if (FALSE) {
  bnname <- "child"
  N <- 1000
  r <- 1

  filepath <- paste0(outpath, sprintf("%s_N%0.0f_r%02.0f.rds", bnname, N, r))
  filepath

  file.remove(filepath)
  adjsets_to_file(bnname, N, r, maxconf, outpath)

  test <- readRDS(filepath)
  names(test)
  lapply(test, lengths)
  lapply(test, attr, "toc")
  attr(test, "toc")
}

# run ----
cl <- parallel::makeCluster(ncores, type="SOCK", outfile = "")
doSNOW::registerDoSNOW(cl)
foreach (bnname = par$bnname,
         N = par$N,
         r = par$r,
         .packages = c("Matrix"),
         .verbose = T) %do% adjsets_to_file(bnname, N, r, maxconf, outpat)

on.exit(parallel::stopCluster(cl))
