
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
#'
#' @return
#' * `bdeu_posterior`: an [bida_sparse_array] or [base::array] with dimensions
#'   `nlev[c(y, x)]` with counts for the posterior Dirichlet distributions.
#' @export
#'
#' @examples
#' nlev <- c(3, 3, 3)
#' dim_yx <- nlev[1:2]
#' data <- sapply(nlev, sample.int, size = 10, replace = T)-1
#'
#' Nyxz  <- bdeu_posterior(data, 1, 2:3, 1, nlev)
#' bida:::posterior_mean.bdeu_posterior(Nyxz, ess = 1)
#' mom  <- bida:::posterior_moments.bdeu_posterior(Nyxz, ess = 1)
#'
#' Nyx <- bdeu_posterior(data, 1, 2, 1, nlev)
#' py.dox <- bida:::backdoor_mean.bdeu_posterior(Nyx, ess = 1, dim_yx)
#' all.equal(py.dox, posterior_mean.bdeu_posterior(Nyx, ess = 1))
#'
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

#' @rdname bdeu_posterior
#' @return
#' - `posterior_mean.bdeu_posterior`: posterior means of the CPT, stored in an array with dimensions `dim(Nyx)`
#' @export
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

#' @rdname bdeu_posterior
#' @return
#' - `posterior_mean.bdeu_posterior`: posterior means of the CPT, stored in an array with dimensions `dim(Nyx)`
#' @export
posterior_moments.bdeu_posterior <- function(Nyx, ess) {
  dims <- dim(Nyx)
  my.x <- posterior_mean.bdeu_posterior(Nyx, ess)
  if (length(dims) == 1) {
    list(mean = my.x,
         cov = (diag(my.x) - tcrossprod(my.x))/(sum(Nyx)+ess+1))
  } else {
    ax   <- as.array(colSums(Nyx))+ess/dims[2]
    cov  <- lapply(seq.int(dims[2]),
                   function(x) (diag(my.x[, x])-tcrossprod(my.x[, x]))/(ax[x]+1))
    list(mean = my.x, cov = cov)
  }
}

#' @rdname bdeu_posterior
#' @return
#' - `posterior_sample.bdeu_posterior`: posterior sample of the CPT, stored in an array with dimensions `c(dim(Nyx), n)`
#' @export
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
        p[, xx, ] <- rDirichlet(n, ayx[, xx], ky)
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
#' @rdname bdeu_posterior
#' @return
#' - `backdoor_sample.bdeu_posterior`: posterior sample of the IPT, stored in an array with dimensions `c(nlev[c(y,x)], n)`
#' @export
backdoor_sample.bdeu_posterior <- function(Nyxz, n, ess, dim_yx) {
  dims <- dim(Nyxz)
  if (length(dims) < 3) {
    p <- posterior_sample.bdeu_posterior(Nyxz, n, ess)

    if (length(dims) == 1) {
      indx <- rep(seq.int(dim_yx[1]), n*dim_yx[2]) + dims[1]*rep(seq.int(n)-1, each = dims[1]*dim_yx[2])
      p <- array(p[indx], c(dim_yx, n))
    }
    return(p)
  } else {
    k   <- prod(dims)
    kyx <- prod(dims[1:2])
    kz  <- k/kyx
    p <- py.xz <- array(0, c(n, dims[1:2]))   # init array for storing samples
    seqx <- seq_len(dims[2])
    sample_py.xz <- function(ayxz, py.xz = NULL) {
      # sample CPTs from a two-dimensional array with Dirichlet-counts
      # ayxz: a r x q vector with dirichlet counts for each level of y and x
      # py.xz: a r x q x n vector with samples to be re-used
      for (xx in seqx) {
        py.xz[,, xx] <- rDirichlet(n, ayxz[, xx], dims[1])
      }
      return(py.xz)
    }

    if (inherits(Nyxz, "bida_sparse_array") && log10(kz) > 4) {
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
        py.xz <- sample_py.xz(replace(ayxz0, yx[indx], ayxz[indx]), py.xz)
        p <- p + pz[, zz]*py.xz
      }

      if (!all_z_observed) {
        # approximate sample by sampling from stick-breaking-process
        pz <- round(pz[, ncol(pz)]*rstick(n, min(10**3, kz-length(uz)), ess), 15)
        for (zz in which(colSums(pz)>0)) {
          py.xz <- sample_py.xz(ayxz0, py.xz)
          p <- p + pz[, zz]*py.xz
        }
      }
    } else {

      # collapse dimension of adjustment set, for indexing below
      dim(Nyxz) <- c(dims[1:2], kz)

      # update counts
      ayxz <- as.array(Nyxz) + ess/(kz * kyx)
      az   <- as.array(colSums(Nyxz, dims = 2)) + ess/kz

      # draw distributions over the adjustment set
      pz <- round(rDirichlet(n, az), 15)
      anyPos <- colSums(pz) > 0  # skip evaluation of probabilities equal to zero

      # apply backdoor formula
      for (zz in seq_along(az)[anyPos]) {
        py.xz <- sample_py.xz(ayxz[, , zz], py.xz)
        p <- p + pz[, zz]*py.xz
      }
    }
    aperm(p, c(2, 3, 1))
  }
}

rstick <- function(n, m, alpha0) {
  v <- matrix(c(rbeta(n*(m-1), 1, alpha0), rep(1, n)), n, m)
  cbind(v[, 1],
        exp(log(v[,-1])+t(apply(log(1-v[, -m, drop = FALSE]), 1, cumsum))))
}


#' @rdname bdeu_posterior
#' @return
#' - `backdoor_mean.bdeu_posterior`: posterior mean of the IPT, stored in an array with dimensions `c(nlev[c(y,x)]`
#' @export
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

    px.z <- axz/rep(az, each = dims[2])
    p <- rowSums(ayxz/rep(px.z, each = dims[1]), dims = 2)
    as.array(p)/sum(az)
  }
}


#' @rdname bdeu_posterior
#' @return
#' - `backdoor_moments.bdeu_posterior`: a list with the following elements:
#'    - `mean`: the posterior mean of the IPT
#'    - `cov`: a matrix list where element \eqn{(x, x')} is the covariance between
#'    the associated vector of intervention probabilities, \eqn{\pi_{Y|x}} and \eqn{\pi_{Y|x'}}
#' @export
backdoor_moments.bdeu_posterior <- function(Nyxz, ess, dim_yx) {
  dims <- dim(Nyxz)
  if (length(dims) < 3) {
    tmp <- posterior_moments.bdeu_posterior(Nyxz, ess)
    my.dox <- array(tmp$mean, dim_yx)
    cov <- matrix(list(), dim_yx[2], dim_yx[2])
    if (length(dims) == 2) {
      cov[diag(dim_yx[2]) == 1] <- tmp$cov
    } else {
      cov[diag(dim_yx[2]) == 1] <- list(tmp$cov)
    }
  } else {
    k <- prod(dims)
    ayxz <- array(as.array(Nyxz+ess/k), c(dims[1:2], k/prod(dims[1:2])))

    axz <- colSums(ayxz)
    az  <- colSums(axz)
    a   <- sum(az)

    mz <- az/a
    my.xz <- ayxz/rep(axz, each = dims[1])
    my.dox  <- rowSums(my.xz*rep(mz, each = prod(dims[1:2])), dims = 2)

    cov <- matrix(list(), dims[2], dims[2])
    seqz <- seq_along(mz)
    tmp0 <- tmp <- matrix(0, dims[1], dims[1])

    for (x in seq_len(dims[2])) {
      tmp <- tmp0
      for (z in seqz) {
        tmp <- tmp + mz[z]/(axz[x, z]+1)*((1+az[z])*diag(my.xz[, x, z])+(axz[x, z]-az[z])*tcrossprod(my.xz[ , x, z], my.xz[ , x, z]))
      }
      cov[[x, x]] <-  1/(a+1)*(tmp-tcrossprod(my.dox[, x], my.dox[, x]))
      if (x == 1) next
      for (xx in seq(1, x-1)) {
        tmp <- tmp0
        for (z in seqz) {
          tmp <- tmp + mz[z]*tcrossprod(my.xz[ , x, z], my.xz[ , xx, z])
        }
        cov[[x, xx]] <- 1/(a+1)*(tmp-tcrossprod(my.dox[, x], my.dox[, xx]))
      }
    }
    cov[upper.tri(cov)] <- cov[lower.tri(cov)]
  }
  list(mean = my.dox, cov = cov)
}
