

files <- list.files("inst/data/", "alarm*|asia*|child*|sachs*", full.names = TRUE)


## INTERVENTIONAL DISTRIBUTION ----
## compare marginal distaributions of sampled data
for (f in files){
  bn <- readRDS(f)

  pdo_bn     <- interv_probs_from_bn(bn, "bn")
  pdo_exact  <- interv_probs_from_bn(bn, "exact")

  plot(unlist(pdo_exact), unlist(pdo_bn), main = f)
  abline(a = 0, b = 1)
}

## OBSERVATIONAL DISTRIBUTION ----
## compare marginal distribution of sampled data
for (f in files){
  bn <- readRDS(f)

  # compute marginals - exact
  probs <- lapply(bn, "[[", "prob")
  indx <- vapply(probs, function(x) is.null(names(dimnames(x))), logical(1))
  for (i in which(indx)){
    names(dimnames(probs[[i]])) <- names(probs)[i]
  }

  # compute expected values: Monte Carlo estimates of marginal probs
  samplesize <- 10**6
  df <- bnlearn::rbn(bn, samplesize)
  mc <- lapply(df, function(x) tabulate(x, nlevels(x))/samplesize)

  # marginals probs using
  exact <- marginal_probs_exact(probs, descendants(bn))
  testthat::expect_equal(lengths(mc), lengths(exact), ignore_attr = T)

  bnn <- marginal_probs_bn(bn)
  testthat::expect_equal(lengths(mc), lengths(bnn), ignore_attr = T)

  cat("\n", f, `\n`)
  print(range(unlist(mc)-unlist(exact)))
  print(range(unlist(mc)-unlist(bnn)))
}




### benchmarking
bn <- readRDS("inst/data/alarm.rds")

microbenchmark::microbenchmark(interv_probs_from_bn(bn, "exact"),
                               interv_probs_from_bn(bn, "bn"),
                               #interv_prob_from_cpts_factors(factors, dag = bnlearn::amat(bn)),
                               times = 10)


profvis::profvis({
  exact <- interv_probs_from_bn(bn, "exact")
  mc <- interv_probs_from_bn(bn, "bn")
})
