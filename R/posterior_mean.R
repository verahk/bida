#' Compute posterior mean
#'
#' @param x object of class bida, bida_params, or backdoor_params
#' @param ... additional arguments
#' @return the posterior mean(s) of the distribution(s) specified by `x`
#' @export
posterior_mean <- function(x, ...){
  UseMethod("posterior_mean", x)
}
