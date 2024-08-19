#' Compute posterior mean
#'
#' @param x object of class bida, bida_params, or backdoor_params
#' @param ... additional arguments
#' @return the posterior mean(s) of the distribution(s) specified by `x`
#' @export
posterior_mean <- function(x, ...){
  UseMethod("posterior_mean", x)
}

#' @export
posterior_mean.bida_pair <- function(x, contrasts = NULL) {
  if (is.null(contrast)) {
    Reduce("+", Map("*", lapply(x$params, backdoor_mean), x$support))
  } else {
    NextMethod()
  }
}

#' @export
posterior_mean.bida_pair_bdeu <- function(x, contrasts) {
  smpl <- posterior_sample(x, 10**3)
  colMeans(vapply(contrasts, function(f) f(smpl), numeric(10**3)))
}
