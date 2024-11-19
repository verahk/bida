

n <- 3
nlev <- rep(3, n)
dag <- matrix(0, n, n)
dag[upper.tri(dag)] <- 1
cpdag <- dag + t(dag)

cpts <- rand_dist(dag, type = "cat", nlev = nlev)
pdos <- interv_probs(cpts, "exact")

tmp  <- pcalg::pdag2allDags(as(cpdag, "matrix"))
dags <- apply(tmp$dags, 1, function(x) matrix(x, n, n), simplify = FALSE)

# bdeu mean
sim_run <- function(N) {
  data <- sample_data_from_cpts(N, cpts)
  for (i in 1:n) {
    ps <- parent_support_from_dags(dags, x = i)[[1]]
    bida_posteriors <- bida_posterior_bdeu_local(ps, data, i, 1:n, 1, nlev)
    pdo_hat <- lapply(bida_posteriors, mean.bida_posterior_bdeu)
    for (j in seq_along(bida_posteriors)) mean.bida_posterior_bdeu(bida_posteriors[[j]])
  }
  ayxz <- bida_posterior_bdeu()
}

parent_support <- parent_support_from_dags(dags)
ps <- matrix(list(), n, n)
for (i in seq_len(n)) {
  tmp <- list(ps$sets[[i]], ps$support[[i]])
  ps <- bida_posterior_bdeu_local(tmp, data, i, 1:n, 1, nlev)
}
