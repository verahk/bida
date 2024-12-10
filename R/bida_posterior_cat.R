#' Compute the BIDA posterior over categorical intervention distributions
#'
#' Compute the approximate posterior distribution over the interventional probability
#' table (IPT) for a cause-effect pair `(x, y)`, given the posterior over a set
#' of adjustment sets.
#' The local version [bida_posterior_cat_local()] can be used the parent sets of `x`
#' is used for adjustment. Allows for some speed gain by re-using the calculations
#' involving `x` and each adjustment set for each effect variable `y` in `ys`.
#'
#' @param ps (list) posterior support over adjustment sets.
#' @param data (integer matrix) matrix with categorical variables
#' @param x (integer) column position of cause variable
#' @param y (integer) column position of effect variable.
#' @param ess (numeric) imaginary sample size
#' @param nlev (integer vector) cardinality of the variables in `data`.
#' @param Nys (list) list with parameters for the posterior dirichlet
#'  distribution over the marginal distribution of variable.
#'  If `NULL` these parameters are computed as needed (default).
#'  When the function is called for each cause-effect pair,
#'  these can be pre-computed for speed and memory.
#'
#' @return an object with class `posterior_bida` and
#'  subclass `bida_posterior_cat`. The local version returns a
#'  list of such objects.
#' @export
#'
#' @examples
#'
#' data(bida_example_cat)
#'
#' # compute dag posterior from sample of DAGs
#' dags <- bida_example$partitionMCMC$traceadd$incidence[-c(1:200)]
#' dag_supp <- dag_support(dags)
#'
#' # compute adjset posterior
#' x <- match("X", colnames(dag))
#' y <- match("Y", colnames(dag))
#' ps <- adjset_support(dag_supp, x, y, "pa")
#'
#' # compute bida posterior
#' pb <- bida_posterior_cat(ps, data, x, y, ess = 1, nlev = nlev)
#' pb
#'
#' posterior_mean(pb)
#' posterior_sample(pb, 10)
#'
#'
#' @rdname bida_posterior_cat
#' @param ys (integer vector) column position(s) of effect variable(s).
#' @export
bida_posterior_cat <- function(ps, data, x, ys, ess, nlev, Nys = NULL) {

  stopifnot(length(x) == 1)
  if (length(ys) == 1 && ys == x) return(NULL)

  sets <- ps[[1]]
  p    <- ps[[2]]
  if (is.matrix(sets)) sets <- apply(sets, 1, function(z) z[!is.na(z)], simplify = FALSE)

  mode(data) <- mode(nlev) <- "integer"
  kx    <- nlev[[x]]        # cardinality of x
  n     <- ncol(data)
  seqy <- seq_along(ys)

  # for each parent set z,
  # compute parameters for each variable in y that is not in c(x, z)
  # - x has no causal effect on (x, z)
  params <- matrix(list(), length(sets), n)  # init
  for (i in seq_along(sets)) {
    z    <- sets[[i]]
    nz <- length(z)

    if (nz > 0) {
      dim_xz <- c(kx, nlev[z])
      kz   <- prod(dim_xz[-1])  # cardinality of adjustment set
      axz0 <- ess/(kx*kz)       # concentration param of dirichlet prior over P(y|x,z)
      dxz  <- data[, c(x, z)]%*%c(1, cumprod(c(kx, nlev[z[-nz]])))
    } else {
      dim_xz <- kx
      axz0 <- ess/kx
      dxz  <- data[, x]
    }

    poss_descendants <- match(ys, c(x, z), 0L) == 0
    for (y in ys[poss_descendants]) {
      ky <- nlev[[y]]
      dyxz <- data[, y] + ky*dxz  # joint outcome (y, x, z)
      uyxz <- c(unique(dyxz))     # unique outcomes
      Nyxz  <- new_bida_sparse_array(value = tabulate(match(dyxz, uyxz)),
                                     index = uyxz,
                                     dim = c(ky, dim_xz),
                                     dimnames = NULL,
                                     default = 0)
      params[[i, y]] <- Nyxz
    }
  }


  # update params for sets that imply a zero-effect and construct posterior_bida-objects
  out <- vector("list", n)
  for (y in ys[match(ys, x, 0L) == 0]) {
    ky <- nlev[[y]]
    isZeroEffect <- vapply(params[, y], is.null, logical(1))
    if (!any(isZeroEffect)) {
      out[[y]] <- new_bida_posterior_cat(params[, y], p, 0, ess, c(ky, kx))
    } else {
      if (is.null(Nys[[y]])) {
        Ny <- array(tabulate(1+data[, y], ky), ky)
      } else {
        Ny <- Nys[[y]]
      }
      zeroprob     <- sum(p[isZeroEffect])
      out[[y]] <- new_bida_posterior_cat(c(list(Ny), params[!isZeroEffect, y]),
                                          c(zeroprob, p[!isZeroEffect]),
                                          zeroprob,
                                          ess,
                                          c(ky, kx))
    }
  }
  out
}


new_bida_posterior_cat <- function(params, p, zeroprob, ess, dim) {
  structure(list(params = params,
                 p = p,
                 zeroprob = zeroprob,
                 ess = ess,
                 dim = dim),
            class = c("posterior_bida", "bida_posterior_cat"))
}


#' @rdname bida_posterior_cat
#' @section Methods:
#' - `posterior_mean()` returns the posterior mean IPT or the posterior mean
#'   causal effect if `contrasts` is specified. Note that the contrasting functions
#'   are applied to the mean IPTs (before computing the weighted sum), which is
#'   only gives the exact posteiror mean for linear transformations.
#'   For non-linear transformations, a Monte-Carlo esimate can be computed using
#'   posterior_sample().
#' @param n (integer) sample size
#' @param contrasts (list) a list of functions for contrasting IPTs and summarize
#'  the effect of interventions into a single value, such as [jsd()].
#'  The functions are assumed to return 0 when the interventional distributions
#'  are identical.
#' @export
posterior_mean.bida_posterior_cat <- function(obj, contrasts = list()) {
  means <- lapply(obj$params,
                  backdoor_mean.bdeu_posterior, ess = obj$ess, dim_yx = obj$dim)
  if (length(contrasts) == 0) {
    Reduce("+", Map("*", means, obj$p))
  } else {
    if (length(means) == 1) {
      vapply(contrasts, function(f) f(means[[1]]), numeric(1))
    } else {
      # apply causal contrast to each ipt
      taus <- vapply(contrasts,
                     function(f) vapply(means, f, numeric(1)),
                     numeric(length(means)))
      # compute weighted average
      colSums(taus*obj$p)
    }
  }
}

#' @rdname bida_posterior_cat
#' @param n sample size
#' @section Methods:
#' - `posterior_sample()` returns a posterior sample of
#' intervention probability tables (IPTs) or causal effects if `contrasts`
#' is specified.
#' @export
posterior_sample.bida_posterior_cat <- function(obj, n, contrasts = list()) {
  # sample adjustment sets
  nG <- stats::rmultinom(n=1, size=n, prob=obj$p)

  if (length(contrasts) == 0) {
    # draw a sample of IPTs according to the sampled adjustment sets
    indx <- nG > 0
    tmp  <- mapply(backdoor_sample.bdeu_posterior,
                   Nyxz = obj$params[indx],
                   n = nG[indx],
                   MoreArgs = list(ess = obj$ess,
                                   dim_yx = obj$dim),
                   USE.NAMES = FALSE,
                   SIMPLIFY = FALSE)
    array(unlist(tmp), c(obj$dim, n))
  } else {
    # draw a set of IPTs and apply the contrasting functions in `contrasts`

    if (obj$zeroprob > 0) {
      # as the contrasting functions are assumed to evaluate to 0 when applied
      # to IPTs representing a zero-effect, do not sample IPTs with adjustment
      # sets that implies a zero-effect
      # - if obj$zeroprob > 0 the first set of parameters are associated with a zero-effect
      out <- matrix(0, n, length(contrasts))
      colnames(out) <- names(contrasts)
      if (nG[1] < n) {
        obj$p[1] <- 0
        nPos <- n-nG[1]
        pdo  <- posterior_sample.bida_posterior_cat(obj, nPos, contrasts = NULL)
        out[seq_len(nPos), ] <- vapply(contrasts, function(f) f(pdo), numeric(nPos))
      }
      return(out)
    } else {
      pdo  <- posterior_sample.bida_posterior_cat(obj, n, contrasts = NULL)
      vapply(contrasts, function(f) f(pdo), numeric(n))
    }
  }
}
