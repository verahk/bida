
#' Update BDeu-prior over a conditional probability table (CPT)
#'
#' @param data (integer matrix) data over a set of categorical variables
#' @param y (integer) column position of outcome variable
#' @param x (integer vector) column position(s) of conditioning variable(s)
#' @param ess (numeric) equivalent sample size
#' @param nlev (integer vector) cardinality of the variables in `data`
#' @param dx (integer vector) joint outcome of `x` in data,
#'  enumerated from 0 to `prod(nlev[x])`-1.
#'  Included for settings where these joint outcomes can be reused.
#'  If `NULL` the function computes the joint outcomes as needed (default).
#' @param sparse (logical) If `TRUE`, the counts are stored in a [bida_sparse_array]
#'  object (default). Otherwise, the counts are stored in a standard array.
#'  Otherwise, the function returns an object of class [base::array]
#' @param `ayx` an array of
#'
#'
#' @return
#' * `posterior_bdeu`: an [bida_sparse_array] or [base::array] with dimensions
#'   `nlev[c(y, x)]`, with parameters of the posterior Dirichlet distributions.
#' @export
#'
#' @examples
posterior_bdeu <- function(data, y, x, ess, nlev, dx = NULL, sparse = TRUE) {
  nx <- length(x)
  if (nx == 0) {
    tabulate(data[, y]+1, nlev[y]) + ess/nlev[y]
  } else {
    if (is.null(dx)) {
      if (nx == 1) {
        dx <- data[, x]
      } else {
        stride  <- c(1, cumprod(nlev[x[-nx]]))
        dx <- data[, x]%*%stride
      }
    }

    dims <- nlev[c(y, x)]
    kyx  <- prod(dims)
    dyx  <- data[, y] + nlev[y]*dx

    if (sparse) {
      uyx <- unique(dyx)
      Nyx <- tabulate(match(dyx, uyx))
      bida_sparse_array(Nyx+ess/kyx, uyx, dims, default = ess/kyx)
    } else {
      Nyx <- array(tabulate(dyx+1, kyx), dims)
      Nyx + ess/kyx
    }
  }
}

mean.posterior_bdeu <- function(ayx, dims = dim(ayx)) {
  if (length(dims) == 1) {
    as.array(ayx)/sum(ayx)
  } else {
    p <- as.array(ayxz)
    p/rep(colSums(p), each = dims[2])
  }
}

sample.posterior_bdeu <- function(ayx, n, dims = dim(ayx)) {
  if (length(dims) == 1) {
    t(bida:::rDirichlet(n, as.array(ayx), dims[1]))
  } else {
    ky <- dims[1]
    kx <- prod(dims[-1])
    p <- array(dim = c(n, ky, kx))

    if (inherits(ayx, "array")) {
      dim(ayx) <- c(ky, kx)
      for (xx in seq_len(kx)) {
        p[,, xx] <- bida:::rDirichlet(n, ayx[, xx], ky)
      }
    } else if (inherits(ayx, "bida_sparse_array")) {

      y  <- ayx$indx%%ky +1
      x  <- ayx$indx%/%ky +1
      ux <- unique(x)
      ayx0 <- rep(ayx$default, dims[1])

      # sample from prior distribution for unobserved parents
      p[,,-ux] <- bida:::rDirichlet(n*(kx-length(ux)), ayx0, ky)

      # sample from posterior for observed parents
      for (xx in ux) {
        indx <- x == xx
        tmp  <- replace(ayx0, y[indx], ayx$value[indx])
        p[,, xx] <- bida:::rDirichlet(n, tmp, ky)
      }
    }
    p <- aperm(p, c(2, 3, 1))
    dim(p) <- c(dims, n)
    return(p)
  }
}

backdoor_sample.posterior_bdeu <- function(ayxz, n, dim_yx, th = .999) {
  dims = dim(ayxz)
  if (length(dims) < 3) {
    p <- sample.posterior_bdeu(ayxz, n, dims)
    if (length(dims) == 1) {
      p <- array(rep(p, each = dim_yx[2]), c(dim_yx, n))
    }
    return(p)
  } else {
    dim(axz) <- c(dims[1:2], prod(dims[-c(1, 2)]))
    stopifnot(inherits(ayxz, "bida_sparse_array"))

    kyx <- prod(dims[-3])
    kz  <- dims[3]

    yx  <- ayxz$index%%kyx +1   # observed joint outcomes of (y, x)
    z   <- ayxz$index%/%kyx +1  # observed joint outcomes of (z)
    uz  <- unique(z)

    # update hyperparams over distribution p(z)
    az_obs   <- bida:::rowsum_fast(ayxz$value, z, uz)  # for each observed level
    az_unobs <- ayxz$default*kyx*(kz-length(uz))       # aggregated for ALL unobserved levels

    # sample a distribution over the observed levels of z and all unobserved levels (if any)
    if (az_unobs > 0) {
      pz <- rDirichlet(n, c(az_obs, az_unobs))
    } else {
      pz <- rDirichlet(n, az_obs)
    }

    sample_py.xz <- function(ayxz) {
      for (xx in seqx) {
        py.xz[, , xx] <<- rDirichlet(n, ayxz[, xx], dims[1])
      }
    }
    # for each observed level of z, sample cpts p(y|x,z) and add to backdoor sum
    p <- py.xz <- array(0, c(n, dims[-3]))
    seqx <- seq_len(kx)
    for (zz in uz) {
      indx <- z == zz
      ayx.z  <- matrix(ayxz$default, dims[1], dims[2])
      ayx.z[yx[indx]] <- ayxz$value[indx]
      sample_py.xz(ayx.z)
      p <- p + pz[, zz]*py.xz
    }

    if (az_unobs > 0) {
      # for each unobserved level of z, sample prob pz and cpts p(y|x,z) from prior,
      # until cumulative prob over z reach threshold (for all samples i = 1, ..., n)
      Fz   <- 1-pz[, ncol(pz)]
      az0   <- ayxz$default*kyx
      ayxz0 <- matrix(ayxz$default, dims[1], dims[2])
      iter <- length(uz) + 1
      while (iter < kz && any(Fz < th)) {
        # sample marginal prob pz from conditional beta distrib
        tmp  <- rbeta(n, az, kz_unobs*az)
        pz   <- (1-Fz)*tmp
        Fz   <- Fz + pz
        iter <- iter + 1

        # update backdoor sums
        sample_py.xz(ayxz0)
        p <- p + pz*py.xz
      }
      # ensure that distribution sums to 1
      pz <- 1-Fz
      sample_py.xz(ayxz0)
      p <- p + pz*py.xz
    }

    aperm(p, c(2, 3, 1))
  }
}
backdoor_mean.posterior_bdeu <- function(ayxz, dim_yx, az = NULL) {
  dims <- dim(ayxz)
  if (length(dims) < 3) {
    array(mean.posterior_bdeu(ayxz, dims), dim_yx)
  } else {
    dim(axz) <- c(dims[1:2], prod(dims[-c(1, 2)]))
    axz <- colSums(ayxz)
    az  <- colSums(axz)
    px.z <- axz/rep(az, each = dims[2])
    p <- rowSums(ayxz/rep(px.z, each = dims[1]), dims = 2)
    as.array(p)/sum(az)
  }
}
