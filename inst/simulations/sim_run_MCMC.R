
#' Simulation: MCMC sampling of DAGs
#'
#' Sample a data set from a bayesian network and draw a sample of DAGs using
#' the MCMC schemes implemented in the `BiDAG`-package.
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


par <- list(local_struct = c("ptree", "none"),
            init = c("pcskel"),
            sample = "order",
            ess = 1,
            edgepf = c("2", "logN"),
            hardlimit = 4,
            N = c(300, 1000, 3000),
            r = seq.int(iterStart, iterStop))

# add params controlling random-cpt generation
par$n <- c(10, 20)
par$k <- c(2, 4)
par$maxdepth <- c(0, .5, 1)

params_to_filename <<- function(par) {
  tmp <- sprintf("n%s_k%s_depth%s_%s_%s_%s_ess%s_epf%s_N%s_r%02.0f.rds",
                 par$n, par$k, par$maxdepth*100, par$init, par$local_struct, par$sample, par$ess, par$edgepf, par$N, par$r)
  stopifnot(length(tmp) == 1)  # fails if any argument is NULL
  tmp
}
expand.grid(par, stringsAsFactors = FALSE)

sim_run <- function(par, verbose = TRUE) {

  N <- par$N
  r <- par$r
  par$edgepf <- ifelse(par$edgepf == "logN", log(N), as.numeric(par$edgepf))

  # import bn
  set.seed(r)
  bn <- sim_load_bn(par)
  nlev <- vapply(bn, function(x) dim(x$prob)[1], integer(1))

  # draw data
  set.seed(N+r)
  data <- bida:::sample_data_from_bn(bn, N)

  # define scorepars
  lookup <- rlang::new_environment()
  scorepar  <- bida:::define_scoreparameters(data,
                                             scoretype = "bdecat",
                                             par = c(par, nlev = list(nlev)),
                                             lookup = lookup)

  # run MCMC ----
  MCMCchain <- bida:::sample_dags(scorepar,
                                  par$init,
                                  par$sample,
                                  hardlimit = par$hardlimit,
                                  verbose = verbose)

  list(par = par,
       lookup = scorepar$lookup,
       MCMCchain = MCMCchain)
}


