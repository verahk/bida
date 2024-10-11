

branch  <- system("git branch --show-current", intern = TRUE)
subdir  <- switch(args$what, "MCMC" = "MCMCchains", "bida" = "results")
indir <- paste0("./inst/simulations/", branch, "/MCMCchains/")

res <- readRDS(list.files(indir, ".rds", full.names = T))
MCMCchain <- res$MCMCchain
par <- res$par

N <- par$N
r <- par$r
n <- par$n
nlev <- rep(par$k, par$n)

set.seed(par$r)
system.time(bn <- sim_load_bn(par))

# compute support over unique dags
dag_posterior <- function(MCMCchain) {
  burnin <- 1:200
  dags <- lapply(MCMCchain$traceadd$incidence[-burnin], as.matrix)
  tmp <- unique(dags)
  support <- bida:::rowsum_fast(rep(1/length(dags), length(dags)), dags, tmp)
  list(dags = dags,
       p = support)
}

f <- function(MCMCchain) {
  tmp <- dag_posterior(MCMCchain)
  ps <- bida::parent_support_from_dags(dags, support)

  type <- ifelse(par$local_struct == "none", "cat", par$local_struct)

}


f <- function(y) {
  type <- ifelse(par$local_struct == "none", "cat", par$local_struct)
  pa   <- which(dag[, x] == 1) # true parents

  pairs <- list(
    unknown = bida::bida_pair(type, data, x, y,
                              sets = ps$sets[[x]],
                              support = ps$support[[x]],
                              hyperpar = c(list(nlev = nlev), par)),
    full = bida::bida_pair("cat", data, x, y,
                           sets = ps$sets[[x]],
                           support = ps$support[[x]],
                           hyperpar = c(list(nlev = nlev), par)),
    known = bida::bida_pair(type, data, x, y,
                            sets = matrix(pa, nrow = 1),
                            support = 1,
                            hyperpar = c(list(nlev = nlev), par))
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
profvis::profvis({
  f(y)
})


type = "ptree"
pairs <- means <- list()

y <- 2
f <- function(y) {
  pair <- bida::bida_pair(type, data, x, y,
                          sets = ps$sets[[x]],
                          support = ps$support[[x]],
                          hyperpar = c(list(nlev = nlev), par))
  mean <- bida:::posterior_mean(pair)
  list(pair, mean)
}
prof <- list()
for (y in 2:20) {
  Rprof()
  f(y)
  Rprof(NULL)
  prof[[y]] <- summaryRprof()
}


profvis::profvis(test <- lapply(seq_len(n)[-x], f))

fprofvis::profvis({
  for ( y in seq_len(n)[-x]) {
    pa   <- which(dag[, x] == 1) # true parents

    pairs <- list(
      unknown = bida::bida_pair(type, data, x, y,
                                sets = ps$sets[[x]],
                                support = ps$support[[x]],
                                hyperpar = c(list(nlev = nlev), par)),
      full = bida::bida_pair("cat", data, x, y,
                             sets = ps$sets[[x]],
                             support = ps$support[[x]],
                             hyperpar = c(list(nlev = nlev), par)),
      known = bida::bida_pair(type, data, x, y,
                              sets = matrix(pa, nrow = 1),
                              support = 1,
                              hyperpar = c(list(nlev = nlev), par))
  )

  pdo_hat     <- lapply(pairs, bida::posterior_mean)
  }
})

# estimate intervention distributions ----
## compute support over parent sets
ps <- bida::parent_support_from_dags(dags, support)
tic <- list("compute parent support" = Sys.time())


