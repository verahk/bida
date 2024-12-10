
# WHAT: Compare posterior samples over intervention probability tables with analytical moments
# WHY: Validate routines for sampling
# HOW:
# - Simulate random counts over a set of variables Y (cause), X (effect), Z (adjustment set)
# - Compute posterior over the assoicated intervention distribution
# - Sample from the posteiror and compare with analytical moments

if (FALSE) {
  # define routine that draws a sample of IPTs and compare theoretical with empirical moments
  sim <- function(Nyxz, mom, n, ess = 1) {
    dim_yx = dim(mom$mean)
    out <- list()
    tic <- Sys.time()
    smpl <- backdoor_sample.bdeu_posterior(Nyxz, n, ess, dim_yx)
    cat("dim =", dim(Nyxz), "n =", n, "\n")
    print(Sys.time()-tic)

    est <- apply(smpl, 1:2, mean)
    out$mean <- mean((est-mom$mean)**2)

    out$cov <- array(NA, dim = dim(mom$cov))
    smpl <- smpl-rep(est, n)  # center
    for (x in seq_len(dim_yx[2])) {
      for (xx in seq_len(dim_yx[2])) {
        if (is.null(mom$cov[[x, xx]])) next
        S <- apply(smpl, 3, function(p) tcrossprod(p[, x], p[, xx]))
        est <- rowSums(S)/(n-1)
        out$cov[x, xx] <- mean((est-mom$cov[[x, xx]])**2)
      }
    }
    return(out)
  }


  samplesizes <- round(10**seq.int(2, 4, length.out = 21))
  nlev <- 2:4
  dim_yx <- nlev[1:2]

  # One adjustment variable -----
  Nyxz <- bida_sparse_array(sample.int(100, size = prod(nlev)), 1:prod(nlev)-1, nlev)
  mom  <- backdoor_moments.bdeu_posterior(Nyxz, 1, dim_yx)
  res <- lapply(samplesizes, function(n) sim(Nyxz, mom, n))

  mse <- sapply(res, "[[", "mean")
  plot(samplesizes, mse, type = "l",  main = "Mean")

  mse <- sapply(res, "[[", "cov")
  matplot(samplesizes, t(mse), type = "l",  main = "Covariance")

  # No adjustment variables -----
  Nyxz <- bida_sparse_array(sample.int(100, size = prod(nlev[-3])), 1:prod(nlev[-3])-1, nlev[-3])
  mom  <- backdoor_moments.bdeu_posterior(Nyxz, 1, dim_yx)
  res <- lapply(samplesizes, function(n) sim(Nyxz, mom, n))

  mse <- sapply(res, "[[", "mean")
  plot(samplesizes, mse, type = "l",  main = "Mean")

  mse <- sapply(res, "[[", "cov")
  matplot(samplesizes, t(mse), type = "l",  main = "Covariance")

  # No effect -----
  Nyxz <- bida_sparse_array(sample.int(100, size = prod(nlev[1])), 1:prod(nlev[1])-1, nlev[1])
  mom  <- backdoor_moments.bdeu_posterior(Nyxz, 1, dim_yx)
  res <- lapply(samplesizes, function(n) sim(Nyxz, mom, n))

  mse <- sapply(res, "[[", "mean")
  plot(samplesizes, mse, type = "l",  main = "Mean")

  mse <- sapply(res, "[[", "cov")
  matplot(samplesizes, t(mse), type = "l",  main = "Covariance")


  # high-dimensional
  Nyxz <- bida_sparse_array(sample.int(100, size = prod(nlev)), 1:prod(nlev)-1, nlev)
  dim(Nyxz)[3] <- 2**10  # add non-observed levels of adjustment variable
  smpl <- backdoor_sample.bdeu_posterior(Nyxz, 10**3, 1, dim_yx)
  mom  <- backdoor_moments.bdeu_posterior(Nyxz, 1, dim_yx)
  res <- lapply(samplesizes, function(n) sim(Nyxz, mom, n))

  mse <- sapply(res, "[[", "mean")
  plot(samplesizes, mse, type = "l",  main = "Mean")

  mse <- sapply(res, "[[", "cov")
  matplot(samplesizes, t(mse), type = "l",  main = "Covariance")

  # very high-dimensional
  Nyxz <- bida_sparse_array(sample.int(100, size = prod(nlev)), 1:prod(nlev)-1, nlev)
  dim(Nyxz)[3] <- 10**4+1  # add non-observed levels of adjustment variable
  mom  <- backdoor_moments.bdeu_posterior(Nyxz, 1, dim_yx)
  res <- lapply(samplesizes, function(n) sim(Nyxz, mom, n))

  mse <- sapply(res, "[[", "mean")
  plot(samplesizes, mse, type = "l",  main = "Mean")

  mse <- sapply(res, "[[", "cov")
  matplot(samplesizes, t(mse), type = "l",  main = "Covariance")

}
