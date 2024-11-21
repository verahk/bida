

# WHAT: Small simulation experiment to see that implemented BIDA method
# WHY: See that implemented BIDA-functions works as expected - that the posterior
#      means of the causal parameters converges towards zero given the true DAG.
# HOW: Sample data sets from random networks for small-medium-large sample sizes
#      Apply bida() given the true DAG and compute point-estimates.
#      Plot RMSE-values to check if the error approaches zero.

if (FALSE) {
  n  <- 10
  nlev <- sample.int(3, n, TRUE)+1
  params <- list(nlev = nlev, ess = 1)

  sim_run <- function(r) {
    set.seed(r)
    n    <- 10
    bn   <- rand_bn(n, d = 4, "cat", nlev = nlev)
    dag  <- bnlearn::amat(bn)
    dag0 <- dag*0
    dindx <- diag(n) == 1
    pdos <- interv_probs(bn, "exact")[!dindx]
    tau  <- vapply(pdos, jsd, numeric(1))

    dnames <- list(N = c(100, 300, 1000, 3000, 10**4),
                   var = c("pdo", "tau"))
    res <- array(list(), lengths(dnames), dnames)
    for (N in dnames$N) {
      cat("iter =", r, "N = ", N, "\n")
      set.seed(r+N)
      data <- sample_data_from_bn(bn, N)
      fits <- list(pa = bida(list(dag), data, adjset = "pa", params = params),
                   o  = bida(list(dag), data, adjset = "o", params = params),
                   marg = bida(list(dag0), data, adjset = "o", params = params),
                   cond = bida(list(dag0), data, adjset = "pa", params = params))

      means <- lapply(fits, posterior_mean)
      diff  <- matrix(unlist(means), ncol = length(fits))-unlist(pdos)
      colnames(diff) <- names(fits)
      res[[paste0(N), "pdo"]] <- colMeans(diff**2)

      means <- lapply(fits, posterior_mean, contrasts = list(jsd = jsd))
      diff  <- matrix(unlist(means), ncol = length(fits))-tau
      colnames(diff) <- names(fits)
      res[[paste0(N), "tau"]] <- colMeans(diff**2)
    }
    data.frame(expand.grid(dnames), do.call(rbind, res))
  }


  dfs <- lapply(1:30, function(r) sim_run(r))
  df  <- dplyr::bind_rows(dfs, .id = "r")
  df$N <- factor(df$N)
  library(dplyr)
  library(ggplot2)

  df %>%
    tidyr::pivot_longer(-any_of(c("r", "N", "var"))) %>%
    ggplot(aes(N, value, fill = name)) +
    facet_wrap(var~., scales = "free") +
    geom_boxplot()
}
