

#' Class: `bida_pair`
#'
#' Defines a posterior mixture distribution over the
#' parameters defining intervention distribution p(y|do(x)).
#'
#' @name bida_pair
#' @param type
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
#' @param hyperpar list of parameters, depending on the prior distribution. See details.
#' @param lookup
#' @param nlev
#' @param
#'
#' @details
#' The `hyperpar` argument is a list that contains parameters of the distribution
#' or other variables required to update the parameters.
#'
#' * Categorical data:
#'   - `nlev` (integer vector) cardinality of each variable.
#'   - `ess` (numeric) imaginary sample size.
#'   Also, if type `ldag` or `tree`, parameters for optimizing a partition of the CPT.
#'   See [optimize_partition].
#'
#' @return
#' An object of class [bida_pair] is a list that contains:
#' - `params`:
#' - `support`:
#' - `zerosupp`:
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
bida_pair <- function(type, data, x, y, sets, support, par, lookup = NULL) {

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
    tmp <- backdoor_params(data, x, y, y, nlev)
    params <- c(list(tmp), params[!indx])
    support <- c(zerosupp, support[!indx])
  } else {
    zerosupp <- 0
  }

  new_bida_pair(type, x, y, params, support, zerosupp, list(dim = nlev[c(y, x)]))
}

new_bida_pair <- function(type, x, y, params, support, zerosupp, par) {
  tmp <- list(x = x,
              y = y,
              params = params,
              support = support,
              zerosupp = zerosupp)

  if (type %in% c("cat", "ldag", "tree")) {
    tmp$dim  <- par$dim
    subclass <- "bida_pair_bdeu"
  }

  structure(tmp, class("bida_pair", "bida_pair_bdeu"))
}

print.bida_pair <- function(x) {
  cat(sprintf("\nAn object of class %s:", class(x)[2]),
      sprintf("\nx = %s, y = %s", x$x, x$y),
      sprintf("\nUnique adjustment sets: %s", length(x$params)),
      sprintf("\nSupport for zero-effect: %s", x$support[1]))
}

