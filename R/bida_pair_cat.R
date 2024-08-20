




#' Class: `bida_pair`
#'
#' Compute parameters defining the posterior mixture distribution over the
#' parameters defining interventional distribution p(y|do(x)).
#'
#' @name bida_pair_cat
#' @param data a N-by-n data matrix
#' @param x (integer)
#' column position of cause variable
#' @param y (integer)
#' column position of effect variable
#' @param sets (integer matrix)
#' each row contains column positions of a valid adjustment set.
#' NA entries in each row is ignored.
#' @param support (numeric vector)
#' a vector with support for each adjustment set. Should be of length `nrow(sets)`
#' and the entries should sum to 1 (not checked).
#' @param nlev (integer vector)
#' cardinality of each variable.
#' @param ess (numeric)
#' imaginary sample size. Defaults to 1.
#' @return
#' An object of class [bida_pair_cat] is a list that contains:
#' - `params`:
#' - `support`:
#' - `nlev`:
#' - `ess`:
#' @examples
#' nlev <- rep(2, 3)
#' lev  <- lapply(nlev-1, seq.int, from = 0)
#' data <- as.matrix(expand.grid(lev))
#' sets <- matrix(c(NA, 2), nrow = 3)
#' fit <- bida_pair_cat(data, 1, 2, sets, rep(1/3, 3), nlev, ess = 0)
#' fit
#' # compute posterior mean
#' posterior_mean(fit)
#' # sample from postrior
#' posterior_sample(fit)
#'
#'
bida_pair_cat <- function(data, x, y, sets, support, nlev, ess = 1) {

  # indicator for zero-effects
  indx     <- rowSums(sets == y, na.rm = T) > 0

  # compute backdoor params for non-zero effects
  params <- vector("list", nrow(sets))
  for (r in seq_along(params)[!indx]) {
    z <- sets[r, ]
    params[[r]] <- backdoor_params_cat(data, x, y, z[!is.na(z)], nlev)
  }

  # compute backdoor params and support for zero-effects
  if (any(indx)) {
    zerosupp <- sum(support[indx])
    tmp <- backdoor_params_cat(data, x, y, y, nlev)
    params <- c(list(tmp), params[!indx])
    support <- c(zerosupp, support[!indx])
  } else {
    zerosupp <- 0
  }

  structure(list(x = x,
                 y = y,
                 params = params,
                 support = support,
                 zerosupp = zerosupp,
                 nlev = nlev[c(x, y)],
                 ess = ess),
            class = c("bida_pair", "bida_pair_cat"))
}







#' @rdname bida_pair
#' @return
#' - `posterior_sample`: A kx-by-ky matrix with posterior mean values.
#' @export
posterior_mean.bida_pair_cat <- function(x, n = 10**3, ace_funs = NULL) {

    if (is.null(ace_funs)) {
      means <- lapply(x$params,
                      function(v) posterior_mean.backdoor_params_cat(v, x$ess, x$nlev[1]))
      Reduce("+", Map("*", means, x$support))
    } else {
      smpl <- posterior_sample(x, n, ace_funs)
      colSums(smpl)/n
    }
}

#' @rdname bida_pair
#' @return
#' - `posterior_sample`: if `ace_funs == NULL` a three dimensional array with IPTs.
#'    Otherwise a matrix with the
#' @export
posterior_sample.bida_pair_cat <- function(x, n = 10**3, ace_funs = NULL) {

  # sample adjustment sets nG
  if (length(n) == 1 && length(x$support) > 1) {
    nG <- stats::rmultinom(n=1,
                           size=n,
                           prob=x$support)
  } else {
    nG <- n
    n <- sum(n)
  }

  if (is.null(ace_funs)) {
    # sample intervention probability tables
    indx <- nG > 0
    tmp  <- mapply(function(par, n) posterior_sample.backdoor_params_cat(par, n, x$ess, x$nlev[1]),
                   par = x$params[indx],
                   n = nG[indx],
                   SIMPLIFY = T)

    array(unlist(tmp), c(x$nlev, n))

  } else {
    # sample causal effect by transforming intervention probability tables
    if (x$zerosupp == 0) {
      tmp <- posterior_sample(x, nG, NULL)
      vapply(ace_funs, function(f) f(tmp), numeric(n))
    } else if (x$zerosupp == 1 || nG[1] == n) {
      tmp <- matrix(0, nrow = n, ncol = length(ace_funs))
      colnames(tmp) <- names(ace_funs)
      tmp
    } else {
      zeros <- rep(0, nG[1])
      tmp <- posterior_sample(x, replace(nG, 1, 0), NULL)
      vapply(ace_funs, function(f) c(zeros, f(tmp)), numeric(n))
    }
  }
}
