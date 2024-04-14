

#' Compute posterior intervention probability
#'
#' @details
#' - [posterior_sample] return samples from the posterior mixture distribution,
#'   First, adjustment sets from the posterior distribution \eqn{p(G_{ij}|D)} are sampled.
#'   Given the adjustment set of each observation, one can sample observational
#'   parameters from the corresponding (local) Dirichlet distributions.
#'   The interventional distributions are then given by applying the back-door-formula.
#'   If specified, [ace_funs] is applied to the sample of interventional CPTs,
#'   in which case the function returns a list of matrices, one matrix for each
#'   function in [ace_funs], where element (i, j) contains the transformed CPTs.
#'
#'
#' @param ps   nested list with parent_sets and parent score for each variable. See X.
#' @param data data frame with observered variables
#' @param nlevels integer vector. Number of categories of each variable
#' @param ess  imaginary sample size for bdeu-prior
#' @param th   threshold value for discarding parent set. Sorted by the
#' @param dirapprox bolean. If TRUE the posterior distribution is approximated as a Dirichlet.
#'
#' @return a list with parameters of the approximate posterior distribution of the
#'         intervention probabilities.

bida_posterior <- function(ps, data, type, hyperparams){


  n <- ncol(data)
  params  <- matrix(list(), n, n)

  if (is.data.frame(data)) data <- as.matrix(data)

  # check of adjustment sets is specific for pairs (x, y)
  isSpecific <- length(ps$sets) == n**2

  # compute mixture distrib parameters for each x, y
  if (isSpecific){
    for (x in seq_len(n)) {
      for (y in seq_len(n)[-x]){
        if (!is.null(ps$sets[[x, y]])) {
          params[[x, y]] <- bida_posterior_pair(data, type, x, y, ps$sets[[x, y]], ps$support[[x, y]], hyperparams)
        }
      }
    }
  } else {
    for (x in seq_len(n)) {
      for (y in seq_len(n)[-x]){
        if (!is.null(ps$sets[[x]])) {
          params[[x, y]] <- bida_posterior_pair(data, type, x, y, ps$sets[[x, y]], ps$support[[x, y]], hyperparams)
        }
      }
    }
  }


  bp <- structure(params,
                  type = type,
                  class  = c("bida_posterior"))
  return(bp)
}





