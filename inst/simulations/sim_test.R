


filename <- "n20_k4_depth50_pcskel_none_order_ess1_epflogN_N300_r23.rds"

indir <- "./inst/simulations/sim_slearn_with_lstruct5/MCMCchains/"
res <- readRDS(paste0(indir, filename))
par <- res$par

N <- res$par$N
r <- res$par$r
bn <- sim_load_bn(par)

set.seed(N+r)
data <- sample_data_from_bn(bn, N)


# MCMC ----
smpl <- sample_dags()
