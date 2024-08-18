

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
#'
#' x <- 1 # cause node
#' y <- 2 # effect node
#'
#' # construct parent sets
#' sets <- matrix(NA, 4, 2)
#' sets[1, ]          # no parents
#' sets[2, 1] <- 3    # single parent
#' sets[3, ]  <- 3:4  # multiple parents
#' sets[4, 1] <- y    # parent set includes node `y`
#' sets
#'
#' support <- rep(1/nrow(sets), nrow(sets)) # uniform support
#'
#' # categorical data ----
#' nlev <- 1:4 +1
#' lev  <- lapply(nlev-1, seq.int, from = 0)
#' data <- as.matrix(expand.grid(lev))
#' hyperpar <- list(nlev = nlev,  levels = lev, ess = 1)
#'
#' pair <- bida_pair("cat", data, x, y, sets, support, hyperpar)
#' str(pair, max.level = 1)
#'
#'
#' # compute posterior mean
#' posterior_mean(pair)
#'
#' # sample from postrior
#' posterior_sample(pair)
#'
bida_pair <- function(type, data, x, y, sets, support, hyperpar, lookup = NULL) {

  # indicator for zero-effects
  indx     <- rowSums(sets == y, na.rm = T) > 0

  # compute backdoor params for non-zero effects
  params <- vector("list", nrow(sets))
  for (r in seq_along(params)[!indx]) {
    z <- sets[r, ]
    params[[r]] <- backdoor_params(type, data, x, y, z[!is.na(z)], hyperpar, lookup)
  }

  # compute backdoor params and support for zero-effects
  zeropar <- backdoor_params(type, data, x, y, y, hyperpar, lookup)
  zerosupp <- sum(support[indx])

  # store parameters in a list, where the zero-effects are placed first
  params <- c(list(zeropar), params[!indx])
  support <- c(zerosupp, support[!indx])

  if (match(type, c("cat", "ldag", "tree"), 0L) > 0) {
    new_bida_pair_bdeu(x, y, params, support, zerosupp, dim = hyperpar$nlev[c(y, x)])
  }
}

new_bida_pair_bdeu <- function(x, y, params, support, zerosupp, dim){
  structure(list(x = x,
                 y = y,
                 params = params,
                 support = support,
                 zerosupp = zerosupp,
                 dim = dim),
            class = c("bida_pair", "bida_pair_bdeu"))
}


