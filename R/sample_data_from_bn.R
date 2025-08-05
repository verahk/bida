


#' Sample data from bayesian networks
#'
#' Wrapper around function \link{bnlearn::rbn},
#' that resamples data sets until all variables has at least 2 observed levels.
#'
#'
#' @param bn an object of class bn.fit.dnet
#' @param samplesize  sample size
#' @param maxiter (integer) maximum numer of iterations
#'
#' @return a data frame with factor variables
#' @export
sample_data_from_bn <- function(bn, samplesize, maxiter = 10**4) {
  nlev_obs <- rep(0, length(bn))
  iter    <- 0
  while (any(nlev_obs<2) && iter <= maxiter) {
    df <- bnlearn::rbn(bn, samplesize)
    nlev_obs <- vapply(df, function(x) length(unique(x)), integer(1))
    iter <- iter +1
  }
  stopifnot(iter < maxiter)

  # convert to numeric matrix
  mat <- vapply(df, as.double, double(samplesize))-1
  return(mat)
}

