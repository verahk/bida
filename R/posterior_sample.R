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


posterior_sample.bida_pair <- function(x, n, contrasts) {
 if (x$zeroeffect == 1 && !is.null(contrasts)) {
    numeric(n)
  } else {
    NextMethod()
  }
}

#' - `posterior_sample`: if `ace_funs == NULL` a three dimensional array with IPTs.
#'    Otherwise a matrix with the
posterior_sample.bida_pair_bdeu <- function(x, n, contrasts = NULL) {

  # sample adjustment sets
  nG <- stats::rmultinom(n=1, size=n, prob=x$support)
  nPos <- n-nG[1]

  # to avoid sampling zero effects when `contrasts` is specified,
  # sample first only intervention distributions for non-zero effects
  indx <- nG[-1] > 0
  tmp  <- mapply(function(par, n) backdoor_sample(par, n),
                 par = x$params[-1][indx],
                 n = nG[-1][indx],
                 SIMPLIFY = T)

  if (is.null(contrasts)) {
    if (nPos < n) {
      # sample also intervention distributions for zero-effects,
      # i.e. realization of the marginal distributions p(y)
      tmp <- c(tmp, backdoor_sample(x$params[1], nG[1]))
    }
    array(unlist(tmp), c(x$dim, n))
  } else {
    out <- matrix(0, n, length(contrasts))
    colnames(out) <- names(contrasts)
    pdo <- array(unlist(tmp), c(x$dim, nPos))
    out[seq_len(nPos), ] <- vapply(contrasts, function(f) f(pdo), numeric(nPos))

    return(out)
  }
}
