


filename <- "n20_k4_depth0_pcskel_none_order_ess1_epfN_N300_r01.rds"
filename <- "n20_k4_depth0_pcskel_ptree_order_ess1_epfN_N300_r01.rds"
indir <- "./inst/simulations/sim_slearn_with_lstruct5/MCMCchains/"

#filename <- "MCMC_n20_k4_depth0_pcskel_ptree_order_ess1_epfN_N3000_r01.rds"
#indir <- "./inst/simulations/sim_slearn_with_lstruct5/"
res <- readRDS(paste0(indir, filename))
par <- res$par

N <- par$N
r <- par$r
bn <- sim_load_bn(par)
nlev <-  vapply(bn, function(x) dim(x$prob)[1], integer(1))
n <- length(bn)

set.seed(r)
dag <- bnlearn::amat(bn)
pdo <- interv_probs_from_bn(bn, "bn")
tau <- matrix(vapply(pdo, bida:::JSD, numeric(1)), n, n)
dindx <- diag(n) == 1

set.seed(N+r)
data <- sample_data_from_bn(bn, N)


# MCMC ----
lookup <- rlang::new_environment()
scorepar <- define_scoreparameters(data, "bdecat", c(par, list(nlev = nlev)), lookup)
smpl <- sample_dags(scorepar, par$init, par$sample, par$hardlimit, F)
cat("Compare MCMC chains:")
all.equal(res$MCMCchain$traceadd, smpl$traceadd)
all.equal(res$lookup, scorepar$lookup)

# BIDA ----
mse <- list()
compute_mse <- function(x, y) mean( (x-y)**2 )
burnin <- 1:200
dags <- lapply(res$MCMCchain$traceadd$incidence[-burnin], as.matrix)

# known parents, local structure
ps <- parent_support_from_dags(dags)
params <- list(nlev = nlev, ess = par$ess, local_struct = par$local_struct)
bida    <- bida(ps, data, "cat", params, lookup = FALSE)
pdo_hat <- posterior_mean(bida)
mse$known <- mapply(compute_mse, pdo_hat[!dindx], pdo[!dindx])

ps <- parent_support_from_dags(dags)
mse <- list()
tau_hat <- list()
for (name in c("unknown", "full", "known")) {
  if (name == "known") {
    ps <- parent_support_from_dags(list(dag), 1)
  }

  params <- list(nlev = nlev,
                 ess = par$ess,
                 local_struct = ifelse(name == "full", "none", par$local_struct))
  bida    <- bida(ps, data, "cat", params, lookup = FALSE)
  pdo_hat <- posterior_mean(bida)
  tau_hat[[name]] <- vapply(pdo_hat[!dindx], bida:::JSD, numeric(1))
  mse[[name]] <- mapply(compute_mse, pdo_hat[!dindx], pdo[!dindx])
}

out <- list()
out$mse_pdo <- colMeans(do.call(cbind, mse))
out$mse_tau <- colMeans((do.call(cbind, tau_hat) - tau[!dindx])**2)

cat("Compare MSE:")
#filename <- "bida_n20_k4_depth0_pcskel_ptree_order_ess1_epfN_N3000_r01.rds"
#indir <- "./inst/simulations/sim_slearn_with_lstruct5/"
res <- readRDS(paste0(indir, filename))

