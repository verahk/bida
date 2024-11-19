
bida_local <- function(type, dags, w, data, x = NULL, y = NULL, hyperpar,  checksize = NULL) {
  n   <- ncol(data)
  seqn <- seq_len(n)
  if (is.null(x)) x <- seqn
  if (is.null(y)) y <- seqn
  out <- matrix(list(), n, n)

  for (xx in x) {
    if (method == "exact") {
      stop()
    } else if (method == "exact") {
      ps <- parent_support_from_dags_x(dags, w, x, checksize)
    }
    out[xx, ] <- bida_posterior_local_bdeu(ps, )
  }

}

bida <- function(type = "categorical", dags, data, adjset = "pa", hyperpar, x = 1:ncol(data), y = 1:ncol(data)) {
  n <- ncol(data)
  out <- matrix(list(), n, n)

  toc <- c("adjset" = 0, "params" = 0)
  for (xx in x) {
    if (adjset == "pa") {
      tic <- c(Sys.time(), 0, 0)
      ps <- parent_support_from_dags(dags, xx)
      tic[2] <- Sys.time()

      out[xx, ] <- switch(type,
                         "categorical" = bida_posterior_bdeu(ps, data, xx, y, ess = hyperpar$ess, nlev = hyperpar$nlev))
      tic[3] <- Sys.time()
      toc <- toc + as.numeric(diff(tic), units = "secs")
      next
    }

    for (yy in y) {
      tic <- c(Sys.time(), 0, 0)
      ps <- adjset_support_from_dags(adjset, dags, xx, yy)
      tic[2] <- Sys.time()

      out[[xx, yy]] <- switch(type,
                              "categorical" = bida_posterior_bdeu(ps, data, xx, yy, ess = hyperpar$ess, nlev = hyperpar$nlev))
      tic[3] <- Sys.time()
      toc <- toc + as.numeric(diff(tic), units = "secs")
    }
  }
  structure(out,
            toc   = toc,
            class = c("bida", switch(type,
                                     "categorical" = "bida_cat")))
}

bida_posterior_bdeu <- function(ps, data, x, y, ess, nlev) {
  n <- ncol(data)
  out <- matrix(list(), n, n)
  for (yy in y) {
    ky <- nlev[yy]
    # pre-compute hyperpar for distribution over P(yy)
    ay <- tabulate(data[, yy]+1, ky) + ess/ky
    for (xx in x) {
      if (is.null(ps[[xx, yy]])) next
      sets <- ps[[xx, yy]][[1]]
      w    <- ps[[xx, yy]][[2]]
      if (is.matrix(sets)) sets <- apply(sets, 1, function(z) z[!is.na(z)], simplify = FALSE)

      # identify sets that implies a zero effect
      is_zero_eff <- vapply(sets, function(z) any(z == yy), logical(1))
      if (all(is_zero_eff)) {
        out[[xx, yy]] <- new_bida_posterior_bdeu(ay, 1, 1, nlev[c(yy, xx)])
        next
      }

      # init list for storing counts for each distribution over the CPT P(y|x, z)
      params <- vector("list", length(sets))

      # compute counts for each sets that implies a non-zero effect
      for (l in seq_along(sets)[!is_zero_eff]) {
        z <- sets[[l]]
        nz <- length(z)

        if (nz > 0) {
          cump <- cumprod(nlev[z])
          kz   <- cump[length(z)]

          # compute frequency of observed joint outcomes (y, x, z)
          dyxz <- dyx + data[, z]%*%c(kyx, kyx*cump[-length(z)])
          uyxz <- c(unique(dyxz))             # unique outcomes
          Nyxz <- tabulate(match(dyxz, uyxz)) # (non-zero) counts

          # store updated hyperpar in sparse array
          ayxz0 <- ess/(kyx*kz)
          params[[l]] <- new_bida_sparse_array(value = ayxz0 + Nyxz,
                                               index = uyxz,
                                               dim = c(ky, kx, kz),
                                               default = ayxz0)
        } else {
          Nyx <- tabulate(dyx+1, kyx)
          params[[l]] <-  matrix(Nyx + ess/kyx, nrow = kyx)
        }
      }

      if (any(is_zero_eff)) {
        zeroprob <- sum(w[is_zero_eff])
        w <- c(zeroprob, w[!is_zero_eff])
        params  <- c(list(ay), params[, !is_zero_eff])
      } else {
        zeroprob <- 0
      }
      out[[xx, yy]] <- new_bida_posterior(ayxz, w, zeroprob, dims)
    }
  }
}

new_bida_posterior_bdeu <- function(params, w, zeroprob, dim) {
  structure(list(params = params,
                 w = w,
                 zeroprob = zeroprob,
                 dim = dim),
            class = c("bida_posterior", "bida_posterior_bdeu"))
}

posterior_bdeu <- function(data, y, x, ess, nlev, dx = NULL, sparse = TRUE) {
  if (length(x) == 0) {
    tabulate(data[, y]+1, nlev[y]) + ess/nlev[y]
  } else {
    if (is.null(dx)) {
      if (length(x) == 1) {
        dx <- data[, x]
      } else {
        stride  <- c(1, cumprod(nlev[x[-length(x)]]))
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

bida_posterior_local_bdeu <- function(ps, data, y, ess, nlev, ays = NULL) {

  for (xx in x) {
    sets  <- ps[[1]][[xx]]
    w     <- ps[[2]][[xx]]
    if (is.matrix(sets)) sets <- apply(sets, 1, function(z) z[!is.na(z)], simplify = FALSE)

    mode(data) <- mode(nlev) <- "integer"
    kx    <- nlev[[xx]]        # cardinality of x
    n     <- ncol(data)
    seqy <- seq_along(y)

    # update params for all sets that imply a possible non-zero effect
    params <- matrix(list(), length(sets), length(y))  # init
    for (i in seq_along(sets)) {
      zz    <- sets[[i]]
      nz <- length(zz)

      if (nz > 0) {
        kz   <- prod(nlev[zz])  # cardinality of adjustment set
        axz0 <- ess/(kx*kz)    # concentration param of dirichlet prior over P(y|x,z)
        dxz  <- data[, c(xx, zz)]%*%c(1, cumprod(c(kx, nlev[zz[-nz]])))
      } else {
        kz   <- integer(0)
        axz0 <- ess/kx
        dxz  <- data[, xx]
      }

      poss_descendants <- seqy[match(y, c(xx, zz), 0L) == 0]
      for (j in poss_descendants) {
        yy <- y[j]
        # prior params
        ky <- nlev[[yy]]  # cardinality
        ayxz0 <- axz0/ky

        # update
        dyxz <- data[, yy] + ky*dxz       # joint outcome (y, x, z)
        uyxz <- c(unique(dyxz))             # unique outcomes
        Nyxz <- tabulate(match(dyxz, uyxz)) # (non-zero) counts
        ayxz  <- new_bida_sparse_array(value = ayxz0 + Nyxz,
                                       index = uyxz,
                                       dim = c(ky, kx, kz),
                                       default = ayxz0)
        params[[i, j]] <- ayxz
      }
    }


    # update params for sets that imply a zero-effect and construct bida_posterior-objects
    out <- vector("list", length(y))
    for (j in seqy[match(y, x, 0L)==0]) {
      yy <- y[j]
      ky <- nlev[[yy]]
      isZeroEffect <- vapply(params[, j], is.null, logical(1))
      if (!any(isZeroEffect)) {
        out[[xx, yy]] <- new_bida_posterior_bdeu(params[, j], w, 0, c(ky, kx))
      } else {
        if (is.null(ays[[yy]])) {
          # update hyperparam of prior over P(Y|X, Z)
          ays[[yy]] <- ess/ky + array(tabulate(1+data[, yy], ky), ky)
        }
        zeroprob     <- sum(w[isZeroEffect])
        out[[xx, yy]] <- new_bida_posterior_bdeu(c(ays[yy], params[!isZeroEffect, j]),
                                                 c(zeroprob, w[!isZeroEffect]),
                                                 zeroprob = zeroprob,
                                                 dim = c(ky, kx))
      }
    }
  }
  return(out)
}

mean.bida_posterior_bdeu <- function(obj, contrasts = list()) {
  if (is.null(obj)) return(NULL)
  if (length(contrasts) == 0) {
    means <- lapply(obj$params, backdoor_mean.bdeu_posterior, dim_yx = obj$dim)
    Reduce("+", Map("*", means, obj$w))
  } else {
    smpl <- sample.bida_posterior_bdeu(obj, n = 10**3, contrasts = contrasts)
    colMeans(smpl)
  }
}
backdoor_mean.bdeu_posterior <- function(ayxz, dim_yx, az = NULL) {
  dims <- dim(ayxz)
  if (length(dims) < 3) {
    array(mean.bdeu_posterior(ayxz, dims), dim_yx)
  } else {
    axz <- colSums(ayxz)
    az  <- colSums(axz)
    px.z <- axz/rep(az, each = dims[2])
    p <- rowSums(ayxz/rep(px.z, each = dims[1]), dims = 2)
    as.array(p)/sum(az)
  }
}

mean.bdeu_posterior <- function(ayx, dims = dim(ayx)) {
  if (length(dims) == 1) {
    as.array(ayx)/sum(ayx)
  } else {
    p <- as.array(ayxz)
    p/rep(colSums(p), each = dims[2])
  }
}



sample.bida_posterior_bdeu <- function(obj, n, contrasts = list()) {
  if (is.null(obj)) return(NULL)
  # sample adjustment sets
  nG <- stats::rmultinom(n=1, size=n, prob=obj$w)

  if (length(contrasts) == 0) {
    # draw a sample of IPTs according to the sampled adjustment sets
    indx <- nG > 0
    tmp  <- mapply(backdoor_sample.bdeu_posterior,
                   ayxz = obj$params[indx],
                   n = nG[indx],
                   MoreArgs = list(dim_yx = obj$dim),
                   USE.NAMES = FALSE,
                   SIMPLIFY = FALSE)
    array(unlist(tmp), c(obj$dim, n))
  } else {
    # draw a set of IPTs and apply the contrasting functions in `contrasts`

    if (obj$zeroprob > 0) {
      # avoid sampling IPTs that implies a zero-effect by calling this function
      # with an adjusted sample size and a zero weight on the adjustments assoicated with
      # a zero effect
      # - if obj$zeroprob > 0 the first set of parameters are associated with a zero-effect
      out <- matrix(0, n, length(contrasts))
      colnames(out) <- names(contrasts)
      if (nG == n) {
        return(out)
      } else {
        obj$w[1] <- 0
        nPos <- n-nG[1]
        pdo  <- sample.bida_posterior_bdeu(obj, nPos, contrasts = NULL)
        out[seq_len(nPos), ] <- vapply(contrasts, function(f) f(pdo), numeric(nPos))
      }
    } else {
      pdo  <- sample.bida_posterior_bdeu(obj, n, contrasts = NULL)
      vapply(contrasts, function(f) f(pdo), numeric(n))
    }
  }
}

backdoor_sample.bdeu_posterior <- function(ayxz, n, dim_yx, th = .999) {
  dims = dim(ayxz)
  if (length(dims) < 3) {
    p <- sample.bdeu_posterior(ayxz, n, dims)
    if (length(dims) == 1) {
      p <- array(rep(p, each = dim_yx[2]), c(dim_yx, n))
    }
    return(p)
  } else {

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

sample.bdeu_posterior <- function(ayx, n, dims = dim(ayx)) {
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

# profile ----
if (FALSE) {
  bn <- readRDS("inst/data/alarm.rds")
  data <- sample_data_from_bn(bn, 10**3)
}

