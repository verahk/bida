

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
#' @details
#' The `hyperpar` argument is a list that contains parameters of the distribution
#' or other variables required to update the parameters.
#'
#' * Categorical data:
#'   - `nlev` (integer vector) cardinality of each variable.
#'   - `ess` (numeric) imaginary sample size.
#'   Also, if type `ldag` or `tree`, parameters for optimizing a partition of the CPT.
#'   See [optimize_partition].
#' @export
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
#' # uniform support over all sets
#' support <- rep(1/nrow(sets), nrow(sets))
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
#' posterior_mean(pair, contrasts = list(jsd = JSD))
#' stopifnot(all(dim(posterior_mean(pair)) == nlev[c(y, x)]))
#'
#' # sample from postrior
#' posterior_sample(pair, n = 10)
#'
bida_pair <- function(type, data, x, y, sets, support, hyperpar, lookup = NULL) {

  # indicator for zero-effects
  indx     <- rowSums(sets == y, na.rm = T) > 0
  zerosupp <- sum(support[indx])

  # compute backdoor params for non-zero effects
  params <- vector("list", nrow(sets))
  for (r in seq_along(params)[!indx]) {
    z <- sets[r, ]
    params[[r]] <- backdoor_params(type, data, x, y, z[!is.na(z)], hyperpar, lookup)
  }

  # add params for zero effect
  if (zerosupp > 0) {
    params <- c(params[!indx], list(backdoor_params(type, data, x, y, y, hyperpar, lookup)))
    support <- c(support[!indx], zerosupp)
  }


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

#' @rdname posterior_sample
#' @export
posterior_sample.bida_pair <- function(x, n, contrasts = NULL) {
  if (x$support[1] == 1 && !is.null(contrasts)) {
    matrix(0, nrow = n, ncol = length(contrasts), dimnames = list(NULL, names(contrasts)))
  } else {
    NextMethod()
  }
}


# Methods for bida_pair_bdeu ----

#' @rdname bida_pair_bdeu
#' @export
posterior_sample.bida_pair_bdeu <- function(x, n, contrasts = NULL) {

  nlevx <- x$dim[2]

  # sample adjustment sets
  nZ <- stats::rmultinom(n=1, size=n, prob=x$support)

  # indicator for sampled parameter sets, excluding that for a zero-effect (last parameter set)
  indx <- replace(c(nZ > 0), length(nZ), x$zerosupp == 0)

  # sample first only intervention distributions for non-zero effects
  tmp  <- mapply(function(bdeu, size) backdoor_sample(bdeu, size, nlevx),
                 bdeu = x$params[indx],
                 size = nZ[indx],
                 USE.NAMES = FALSE,
                 SIMPLIFY = FALSE)

  if (is.null(contrasts)) {
    if (x$zerosupp > 0 && nZ[length(nZ)] > 0) {
      # sample also intervention distributions for zero-effects,
      # i.e. realization of the marginal distributions p(y)
      tmp <- c(tmp, backdoor_sample(x$params[[1]], nZ[1], nlevx))
    }
    array(unlist(tmp), c(x$dim, n))
  } else {

    # store pos-effect samples in array, for computing contrast
    nPos <- n-(x$zerosupp>0)*nZ[length(nZ)]
    pdo <- array(unlist(tmp), c(x$dim, nPos))

    out <- matrix(0, n, length(contrasts))
    colnames(out) <- names(contrasts)
    out[seq_len(nPos), ] <- vapply(contrasts, function(f) f(pdo), numeric(nPos))

    return(out)
  }
}



#' @export
posterior_mean.bida_pair_bdeu <- function(x, contrasts = NULL) {
  if (is.null(contrasts)) {
    Reduce("+", Map("*", lapply(x$params, backdoor_mean, nlevx = x$dim[2]), x$support))
  } else {
    smpl <- posterior_sample(x, 10**3)
    colMeans(vapply(contrasts, function(f) f(smpl), numeric(10**3)))
  }
}
