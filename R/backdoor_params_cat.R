
#' Categorical intervention distributions
#'
#' Compute sample counts defining discrete intervention distributions
#'
#' @param data a N-by-n matrix with integers representing a set of categorical variables.
#'  Each variable is assumed to be coded as 0, 1, ..., nlev-1.
#' @param x (integer)
#'  column position of cause variable
#' @param y (integer)
#'  column position of effect variable
#' @param z (integer vector)
#'  column position(s) of variable(s) in the backdoor set. A numeric/integer vector of
#'  length 0 indicate an empty adjustment set. If the adjustment set includes `y`,
#'  `x` is assumed to have no effect on `y` and the marginal counts of `y` is returned.
#' @param nlev (integer vector)
#'  vector of length n with the cardinality of each variable.
#' @param sparse (integer)
#'  if the number of possible combinations exceeds this value, `prod(dims[c(y, x, z)]) > sparse`,
#'  the sample counts are stored in a sparse structure.
#'  Otherwise, the counts are stored in an standard array.
#' @return
#'  - `backdoor_params_cat`: an array or a [bida::bida_sparse_array] with sample counts.
#'  - `posterior_mean.backdoor_params_cat`: a kx-by-ky matrix with posterior means
#'  - `posterior_sample.backdoor_params_cat`: a kx-by-ky-by-`samplesize` array with interventional CPTs
#' @details
#'  `backdoor_params_cat` computes sample counts for each configuration `(y, x, z)`.
#'  If `z` contains `y` the object contains the marginal counts of `y` (zero effect).
#'  The dimension of the array is then of length 1.
#'  Otherwise, over the triplet `(y, x, z)`. Sufficient stats for the posterior
#'  distribution over categorical intervention distributions.
#'
#'  `posterior_mean.backdoor_params_cat` computes (analytically) the posterior
#'  mean values of the interventional CPT P(y|do(x))
#'
#'  `posterior_sample.backdoor_params_cat` samples a set of interventional CPTs
#'  from the exact posterior, by sampling the relevant observational distributions
#'  from their respective dirichlet posteriors.
#'
#'
#' @export
backdoor_params_cat <- function(data, x, y, z, nlev, min_sparse = 2**10) {

  if (any(z == y)) {

    # compute marginal counts
    nyxz <- array(tabulate(data[, y]+1, nlev[y]))

  } else {

    subset <- c(y, x, z)
    dims   <- nlev[subset]
    n      <- length(subset)

    # enumerate observed joint outcomes
    cump   <- cumprod(dims)
    yxz    <- c(data[, subset]%*%c(1, cump[-n]))

    if (n <= 3) {
      nyxz   <- array(tabulate(yxz+1, cump[n]), dims)
    } else {

      # adjust dimension, for storing as 3-dim array
      dims <- c(dims[1:2], cump[n]/cump[2])

      if (cump[n] <= min_sparse) {
        nyxz <- array(tabulate(yxz+1, cump[n]), dims)
      } else {
        # count observed outcomes, store as sparse array
        uyxz <- unique(yxz)
        tmp  <- tabulate(match(yxz, uyxz), length(uyxz))
        nyxz <- new_bida_sparse_array(tmp,
                                      index = uyxz,
                                      dim  = dims)
      }
    }
  }
  return(nyxz)
}



#' @rdname backdoor_params_cat
#' @param nyxz counts
#' @param ess (integer)
#'        equivalent sample size
#' @param kx (integer)
#'        number of possible interventions on the cause variable,
#'        Used to replicate the marginal distribution when the causal effect is zero
#'        (when \code{isZeroeffect(x) == TRUE}).
#'        Otherwise the dimensions are given by [dim(x)] and this parameter is ignored.
#' @return
#' -  `posterior_mean` a matrix with the posterior means of each intervention distribution,
#'    i.e. a conditional probability table where element (x, y) equals is $P(Y = y| do(X = x))$
#' @export
posterior_mean.backdoor_params_cat <- function(nyxz, ess, kx) {
  if (is.array(nyxz) || length(nyxz$dim) < 3) {
    posterior_mean.backdoor_params_cat_array(nyxz, ess, kx)
  } else {
    posterior_mean.backdoor_params_cat_sparse(nyxz, ess, kx)
  }
}


posterior_mean.backdoor_params_cat_array <- function(nyxz, ess, kx = 1) {

  dims <- dim(nyxz)
  n <- length(dims)
  ayxz <- nyxz + ess/length(nyxz)

  if (n == 1) {
    p <- array(ayxz/sum(ayxz), c(dims, kx))
  } else if (n  == 2) {
    p <- sweep(ayxz, 2, colSums(ayxz), "/")
  } else {
    axz <- colSums(ayxz)
    az  <- colSums(axz)
    a <- sum(az)
    p <- rowSums(sweep(ayxz, 2:3, sweep(axz, 2, az/a, "/"), "/"), dims = 2)
  }

  return(t(p))
}



posterior_mean.backdoor_params_cat_sparse <- function(obj, ess, kx = 1) {

  dims <- obj$dim

  if (length(dims) < 3) {
    # if there is no adjustment variables, convert to array
    posterior_mean.backdoor_params_cat_array(as.array(obj), ess, kx)
  } else {

    # cardinality of (sets of) variables
    kyxz <- prod(dims)
    kyx  <- kyxz%/%dims[3]

    # find observed levels of adjustment set z
    yxz <- nyxz$index
    yx  <- yxz%%kyx     # observed levels of (y, x)
    z   <- yxz%/%kyx    # observed levels of z
    uz  <- unique(z)    # unique levels
    nuz <- length(uz)   # number of unique levels

    ## compute params of posterior for the observed levels of z
    ayxz <- array(ess/kyxz, c(dims[1:2], nuz))  # prior

    # update
    if (nuz < dims[3]) {
      z <- match(z, uz)-1      # relabel z, such that the observed levels are 0, 1, .. nuz.
      upd  <- yx + kyx*z       # create index with relabeled z
    } else {
      upd  <- yxz
    }
    ayxz[upd+1] <- ess/kyxz + obj$counts
    axz <- colSums(ayxz)
    az  <- colSums(axz)

    # average p(y, x, z)/p(x|z) over observed levels of z..
    p <- rowSums(sweep(ayxz, 2:3, sweep(axz, 2, az, "/"), "/"), dims = 2)

    # .. and add sum over unobserved levels (constant for each x, z)
    p <- p + (1-nuz/dims[3])*1/dims[1]

    # normalize
    p <- p/(ess+sum(obj$counts))
    return(t(p))
  }
}


#' @rdname backdoor_params_cat
#' @inheritParams posterior_mean.backdoor_params_cat
#' @param samplesize (integer)
#' sample size
#' @return
#' - `posterior_sample` a 3-dimensional array with samples of the intervention
#'    distribution parameters drawn from the posterior
#' @export
posterior_sample.backdoor_params_cat <- function(nyxz, samplesize, ess, kx) {
  #NextMethod()
  if (is.array(nyxz)) {
    posterior_sample.backdoor_params_cat_array(nyxz, samplesize, ess, kx)
  } else {
    posterior_sample.backdoor_params_cat_sparse(nyxz, samplesize, ess, kx)
  }
}


posterior_sample.backdoor_params_cat_array <- function(nyxz, samplesize, ess, kx) {

  dims <- dim(nyxz)
  n <- length(dims)

  # update params
  ayxz <- nyxz + ess/length(nyxz)


  if (n == 1){
    # replicate samples for each level of x (kx)
    p <- array(rDirichlet(samplesize, ayxz, length(ayxz)),
               dim = c(samplesize, length(ayxz), kx))

  } else if (n == 2) {
    # init matrix
    p  <- array(0, dim = c(samplesize, dims[-3]))

    # draw posterior sample for each level of x
    for (x in seq_len(dims[2])) {
      p[,, x] <- rDirichlet(samplesize, ayxz[, x], dims[1])
    }

  } else {

    # init matrix
    p  <- array(0, dim = c(samplesize, dims[-3]))

    # draw pz
    az <- colSums(ayxz, dims = 2)
    pz <- round(rDirichlet(samplesize, az), 15)
    anyPos <- colSums(pz > 0) > 0

    # apply backdoor formula
    seqx <- seq_len(dims[2])
    py.xz  <- p

    for (z in seq_along(az)[anyPos]) {

      # sample conditional distrib p(y|x, z)
      for (x in seqx) {
        py.xz[,, x] <- rDirichlet(samplesize, ayxz[, x, z], dims[1])
      }

      # update backdoor sum
      p <- p + pz[, z]*py.xz
    }
  }

  aperm(p, c(3, 2, 1))
}




posterior_sample.backdoor_params_cat_sparse <- function(obj, samplesize, ess, kx){

  stopifnot(samplesize > 0)
  dims <- obj$dim

  if (length(dims) < 3) {

    # if there is no adjustment variables, convert to array
    posterior_sample.backdoor_params_cat_array(as.array(obj), samplesize, ess, kx)

  } else {

    # init n-ky-kx dimensional vector for storing do-probs
    p <- array(0, dim = c(samplesize, dims[-3]))

    # cardinality of (sets of) variables
    tmp <- cumprod(dims)
    kyxz <- tmp[3]
    kyx <-  tmp[2]
    kz  <- dims[3]

    # find observed levels of (y, x) and adjustment set z
    yxz <- obj$index
    yx  <- yxz%%kyx +1       # observed outcomes of (y, x) combined
    z   <- yxz%/%kyx +1      # observed outcomes of z

    # compute params of posterior over p(z)
    uz  <- unique(z)                                   # unique observed levels
    az  <- rep(ess/kz, kz)                             # prior hyperparams
    az[uz] <- ess/kz + rowsum_fast(obj$counts, z, uz)  # update by adding counts

    # draw p(z)
    pz <- round(rDirichlet(samplesize, az), 15)
    rm(az)

    isPos <- pz > 0
    nPos  <- colSums(isPos)

    # draw p(y|x,z) and apply backdoor-formula, for each z where p(z) > 0
    # - for reusing arrays, loop over adjsets ordered by number of positive samples

    a0  <- array(ess/kyxz, dim = dims[-3])    # prior hyperparams over p(y|x, z) for each z
    seqx <- seq_len(dims[2])
    nn <- Inf

    for (zz in order(nPos, decreasing = T)) {

      if (nPos[zz] < nn) {
        nn <- nPos[zz]
        if (nn == 0) break
        py.xz <- array(0, c(nn, dims[-3]))
      }

      # update hyperparams
      ayxz <- a0
      tmp.indx <- z == zz
      ayxz[yx[tmp.indx]] <- a0[1] + obj$index[tmp.indx]


      # sample conditional distrib p(y|x, z)
      for (xx in seqx) {
        py.xz[,, xx] <- rDirichlet(nn, ayxz[, xx], dims[1])
      }

      # add to backdoor sum
      indx <- isPos[, zz]
      p[indx,,] <- p[indx,,, drop = FALSE] + pz[indx, zz]*py.xz
    }
  }

  aperm(p, c(3, 2, 1))
}
