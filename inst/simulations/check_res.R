
dir_MCMC <- "./inst/simulations/sim_bn2/MCMCchains/"
dir_bida <- "./inst/simulations/sim_bn/results/"
files <- list.files(dir_MCMC, ".rds")

# check if parent sets greater than hardlimit ----
if (FALSE) {
  for (file in files) {
    cat(file, "\n")
    imp <- readRDS(paste0(dir_MCMC, file))
    dags <- lapply(imp$MCMCchain$traceadd$incidence, as.matrix)
    ps1 <- parent_support_from_dags(dags)
    ncols <- vapply(ps1$sets, ncol, integer(1))
    stopifnot(all(ncols <= par$hardlimit+1))
  }
}


# re-compute single sim-setting ----
files <- list.files(dir_MCMC, "asia.*rds")

file <- sample(files, 1)
print(file)
res <- readRDS(paste0(dir_MCMC, file))
par <- res$par

# compute ground truth ----
bn <- readRDS(paste0("inst/data/", par$bnname, ".rds"))
n  <- length(bn)
nlev <- vapply(bn, function(x) dim(x$prob)[1], integer(1))
pdo <- interv_probs_from_bn(bn, "exact")
tau <- matrix(vapply(pdo, bida:::JSD, numeric(1)), n, n)

# simulate data -----
r <- par$r
N <- par$N
set.seed(r+N)
data <- sample_data_from_bn(bn, N)

# check MCMC res ----
scorepar <- define_scoreparameters(data, "bdecat", list(par, nlev = nlev))
tic <- Sys.time()
cat("run MCMC:\n")
smpl <- sample_dags(scorepar, par$init, par$sample, hardlimit = par$hardlimit, verbose = F)

Sys.time()-tic
all.equal(smpl, res$MCMCchain)

# compute parent support ----
ps <- parent_support_from_dags(lapply(smpl$traceadd$incidence, as.matrix))
stopifnot(all(vapply(ps$sets, ncol, integer(1)) < par$hardlimit +1))
# check MSE
type <- ifelse(par$local_struct == "none", "cat", par$local_struct)
hyperpar <- list(ess = par$ess, nlev = nlev, local_struct = par$local_struct)
mse <- array(NA, c(n, n, 2))
tauhat <- matrix(NA, n, n)
for (x in seq_len(n)) {
  cat("estimate interv probs:", x, "\n")
  for (y in seq_len(n)[-x]) {
    tic <- Sys.time()
    nlevx <- dim(pdo[[x, y]])[2]
    pair <- bida_pair(type, data, x, y, ps$sets[[x]], ps$support[[x]], hyperpar)
    est  <- posterior_mean(pair)
    mse[x, y, 1] <- mean( (est-pdo[[x, y]])**2 )

    est  <- tauhat[x, y] <- avg_jsd_array(t(est))
    mse[x, y, 2] <- (est-tau[x, y])**2
  }
  print(Sys.time()-tic)
}

colSums(mse, na.rm = T, dims = 2)/(n*(n-1))
res <- readRDS(paste0(dir_bida, file))
res$mse_pdo
res$mse_tau
res$rank

dindx <- diag(n) == 1
compute_avgppv(tauhat[!dindx], tau[!dindx] > 0)
compute_avgppv(tauhat[!dindx], tau[!dindx] > quantile(tau[!dindx][tau[!dindx] > 0], .8))
