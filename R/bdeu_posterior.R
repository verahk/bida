
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
#' * `bdeu_posterior`: an [bida_sparse_array] or [base::array] with dimensions
#'   `nlev[c(y, x)]`, with parameters of the posterior Dirichlet distributions.
#' @export
#'
#' @examples
bdeu_posterior <- function(data, y, x, ess, nlev, dx = NULL, sparse = TRUE) {
  nx <- length(x)
  if (nx == 0) {
    array(tabulate(data[, y]+1, nlev[y]), nlev[y])
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
      uyx <- c(unique(dyx))
      Nyx <- tabulate(match(dyx, uyx))
      new_bida_sparse_array(Nyx, uyx, dims, dimnames = NULL, default = 0)
    } else {
      array(tabulate(dyx+1, kyx), dims)
    }
  }
}

posterior_mean.bdeu_posterior <- function(Nyx, ess) {
  dims <- dim(Nyx)
  if (length(dims) == 1) {
    (as.array(Nyx) + ess/dims[1])/(sum(Nyx) + ess)
  } else {
    ayx0 <- ess/prod(dims)
    p <- as.array(Nyx) + ayx0
    p/rep(colSums(p), each = dims[1])
  }
}

posterior_sample.bdeu_posterior <- function(Nyx, n, ess) {
  dims <- dim(Nyx)
  if (length(dims) == 1) {
    t(bida:::rDirichlet(n, as.array(Nyx)+ess/dims[1], dims[1]))
  } else {
    ky <- dims[1]
    kx <- prod(dims[-1])
    ayx0 <- ess/(ky*kx)

    p <- array(dim = c(n, kx, ky))  # note margin n-x-y
    if (inherits(Nyx, "array")) {
      ayx <- array(Nyx + ayx0, c(ky, kx))
      for (xx in seq_len(kx)) {
        p[, xx, ] <- bida:::rDirichlet(n, ayx[, xx], ky)
      }
    } else if (inherits(Nyx, "bida_sparse_array")) {
      y  <- Nyx$index%%ky +1
      x  <- Nyx$index%/%ky +1
      ux <- unique(x)

      # sample from posterior for observed parents
      ayx0 <- rep(ess/prod(dims), dims[1])
      ayx  <- Nyx$value + ayx0[1]
      for (xx in ux) {
        indx <- x == xx
        p[, xx, ] <- bida:::rDirichlet(n, replace(ayx0, y[indx],  ayx[indx]), ky)
      }
      # sample from prior distribution for unobserved parents
      p[, -ux , ] <- bida:::rDirichlet(n*(kx-length(ux)), ayx0, ky)
    }
    p <- aperm(p, c(3, 2, 1))
    dim(p) <- c(dims, n)
    return(p)
  }
}

backdoor_sample.bdeu_posterior <- function(Nyxz, n, ess, dim_yx, th = .999) {
  dims <- dim(Nyxz)
  if (length(dims) < 3) {
    p <- posterior_sample.bdeu_posterior(Nyxz, n, ess)
    if (length(dims) == 1) {
      p <- array(rep(p, each = dim_yx[2]), c(dim_yx, n))
    }
    return(p)
  } else {
    kyx <- prod(dims[1:2])
    kz  <- prod(dims)/kyx
    p <- py.xz <- array(0, c(n, dims[1:2]))  # init array for storing samples
    seqx <- seq_len(dims[2])
    sample_py.xz <- function(ayxz) {
      # sample CPTs from a two-dimensional array with Dirichlet-counts
      for (xx in seqx) {
        py.xz[,, xx] <<- rDirichlet(n, ayxz[, xx], dims[1])
      }
    }

    if (inherits(Nyxz, "bida_sparse_array")) {
      yx  <- Nyxz$index%%kyx +1   # observed joint outcomes of (y, x)
      z   <- Nyxz$index%/%kyx +1  # observed joint outcomes of (z)
      uz  <- unique(z)

      all_z_observed <- length(uz) == kz

      # sample distributions over the observed levels of z and all unobserved levels (if any)
      az0 <- ess/kz
      az  <- az0 + bida:::rowsum_fast(Nyxz$value, z, uz) # updated counts for observed levels
      if (all_z_observed) {
        pz <- rDirichlet(n, az)
      } else {
        pz <- rDirichlet(n, c(az, az0*(kz-length(uz))))
      }

      # for each observed level of z, sample cpts p(y|x,z) and add to backdoor sum
      ayxz0 <- array(az0/kyx, dims[1:2])  # prior counts for each z
      ayxz  <- Nyxz$value + az0/kyx       # posterior counts for observed outcomes (y, x, z)
      for (zz in seq_along(uz)) {
        indx <- z == uz[zz]
        sample_py.xz(replace(ayxz0, yx[indx], ayxz[indx]))
        p <- p + pz[, zz]*py.xz
      }

      if (!all_z_observed) {
        # for each unobserved level of z, sample prob pz and cpts p(y|x,z) from prior,
        # until cumulative prob over z reach threshold (for all samples i = 1, ..., n)
        Fz   <- 1-pz[, ncol(pz)]
        iter <- kz-length(uz) # number of elements not yet sampled
        while (iter > 1 && any(Fz < th)) {
          # sample marginal prob pz from conditional beta distrib
          iter <- iter-1
          tmp  <- rbeta(n, az0, iter*az0)
          pz   <- (1-Fz)*tmp
          Fz   <- Fz + pz
          # update backdoor sums
          sample_py.xz(ayxz0)
          p <- p + pz*py.xz
        }
        # ensure that distribution sums to 1
        sample_py.xz(ayxz0)
        p <- p + (1-Fz)*py.xz
      }
    } else {
      dim(Nyxz) <- c(dims[1:2], prod(dims[-c(1, 2)]))
      az <- colSums(Nyxz, dims = 2) + ess/kz
      ayxz <- Nyxz + ess/(kz*kyx)

      pz <- round(rDirichlet(n, az), 10)
      for (zz in seq_along(az)) {
        if (!any(pz[, zz] > 0)) next
        sample_py.xz(ayxz[,, zz])
        p <- p + py.xz*pz[, zz]
      }
    }
    aperm(p, c(2, 3, 1))
  }
}

backdoor_mean.bdeu_posterior <- function(Nyxz, ess, dim_yx, az = NULL) {
  dims <- dim(Nyxz)
  if (length(dims) < 3) {
    array(posterior_mean.bdeu_posterior(Nyxz, ess), dim_yx)
  } else {
    kz <- prod(dims[-c(1:2)])
    dim(Nyxz) <- c(dims[1:2], kz)

    Nxz <- colSums(Nyxz)
    Nz  <- colSums(Nxz)

    ayxz <- Nyxz + ess/prod(dims)
    axz  <- Nxz + ess/prod(dims[-1])
    az   <- Nz + ess/kz

    px.z <- (axz)/rep(az, each = dims[2])
    p <- rowSums(ayxz/rep(px.z, each = dims[1]), dims = 2)
    as.array(p)/sum(az)
  }
}
