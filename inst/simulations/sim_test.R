

sapply(list.files("./inst/simulations/R", ".R", full.names = T),
       source, echo = T)

indir_MCMC <- "./inst/simulations/sim_slearn_with_lstruct5/MCMCchains/"
indir_bida <- "./inst/simulations/sim_slearn_with_lstruct5/results/"

keep_sampling <- T
while (keep_sampling){
  files <- list.files(indir_MCMC, ".rds")
  filename <- sample(files, 1)
  keep_sampling <- !file.exists(paste0(indir_bida, filename))
}

#filename <- "n20_k2_depth100_pcskel_none_order_ess1_epf2_N300_r20.rds"
print(filename)

#filename <- "MCMC_n20_k4_depth0_pcskel_ptree_order_ess1_epfN_N3000_r01.rds"
#indir <- "./inst/simulations/sim_slearn_with_lstruct5/"
res <- readRDS(paste0(indir_MCMC, filename))
par <- res$par

cat("loading bn..\n")
N <- par$N
r <- par$r
set.seed(r)
bn <- sim_load_bn(par)
nlev <-  vapply(bn, function(x) dim(x$prob)[1], integer(1))
n <- length(bn)

set.seed(N+r)
data <- bida:::sample_data_from_bn(bn, N)

# MCMC ----
cat("BiDAG is running\n")
lookup <- rlang::new_environment()
scorepar <- bida:::define_scoreparameters(data, "bdecat", c(par, list(nlev = nlev)), lookup)
smpl <- bida:::sample_dags(scorepar, par$init, par$sample, par$hardlimit, F)
test <- all.equal(res$MCMCchain$traceadd, smpl$traceadd)
if (!is.logical(test)) {
  print(test)
  print(filename)
  par
  stop()
}


cat("MCMCchain was replicated, compute ground truth intervention probs\n")
set.seed(r)
dag <- bnlearn::amat(bn)
pdo <- bida:::interv_probs_from_bn(bn, "bn")
tau <- matrix(vapply(pdo, bida:::JSD, numeric(1)), n, n)
dindx <- diag(n) == 1


# BIDA ----
cat("Estimate intervention probs\n")
burnin <- 1:200
dags <- lapply(res$MCMCchain$traceadd$incidence[-burnin], as.matrix)

ps <- bida:::parent_support_from_dags(dags)
mse <- tau_hat <- parents <- npar <- list()
compute_mse <- function(x, y) mean( (x-y)**2 )
for (name in c("unknown", "full", "known")) {
  if (name == "known") {
    ps <- bida:::parent_support_from_dags(list(dag), 1)
  }
  params <- list(nlev = nlev,
                 ess = par$ess)
  type <- ifelse(par$local_struct == "none" || name == "full", "cat", par$local_struct)
  bida    <- bida:::bida(ps, data, type, list(nlev = nlev, ess = par$ess), lookup = NULL)
  pdo_hat <- bida:::posterior_mean(bida)
  tau_hat[[name]] <- vapply(pdo_hat[!dindx], bida:::JSD, numeric(1))
  mse[[name]]     <- mapply(compute_mse, pdo_hat[!dindx], pdo[!dindx])

  # compute number of params
  compute_avg_number_of_params <- function(pair){
    tmp <- vapply(pair$params, function(counts) prod(dim(counts)), double(1))
    sum(tmp*pair$support)
  }
  tmp <- vapply(bida[!dindx], compute_avg_number_of_params, numeric(1))
  params[[name]] <- sum(tmp)
}

out <- list()
out$mse_pdo <- colMeans(do.call(cbind, mse))
out$mse_tau <- colMeans((do.call(cbind, tau_hat) - tau[!dindx])**2)
cat("Re-computed MSEs:\n")
print(out$mse_pdo)
print(out$mse_tau)

cat("MSEs loaded from file:\n")
#filename <- "bida_n20_k4_depth0_pcskel_ptree_order_ess1_epfN_N3000_r01.rds"
indir <- "./inst/simulations/sim_slearn_with_lstruct5/results/"
res <- readRDS(paste0(indir_bida, filename))
print(res$mse_pdo)
print(res$mse_tau)

