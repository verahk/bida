#' Sample from posterior distribution specified by bida-objects
#'
#' @param x object of class bida, bida_params, or backdoor_params
#' @param n (integer) samplesize
#' @param ... additional arguments
#' @return a posterior sample from the distribution(s) specified by `x`
#' @export
posterior_sample <- function(x, n, ...){
  UseMethod("posterior_sample")
}


