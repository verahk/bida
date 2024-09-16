

backdoor_params_from_counts <- function(counts, partition) {
  dims <- dim(counts)
  kz <- prod(dims[z])
  ky <- prod(dims[y])
  kx <- prod(dims[x])

  perm <- c(z, x, y)
  if (is.unsorted(perm)) counts <- aperm(counts, perm)
  dim(counts) <- c(kz, kx, ky)

}



#' Title
#'
#' @param Nyx
#' @param ess
#' @param partition
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' Nyxz <- bida_sparse_array(1:3, 0:2, c(3, 3, 3))
#' partition <- list(0:2, 3:27)
mean_bdeu <- function(Nyx, ess = 1, partition = NULL) {
  dims <- dim(Nyx)
  r <- dims[1]
  q <- prod(dims[-1])
  if (q == 1) {
    (Nyx + ess/r)/(sum(Nyx)+ess)
  } else if (is.null(partition) || length(partition) == q) {
    a0 <- ess/(r*q)
    (Nyx+a0)/rep(colSums(Nyx) + a0*r, each = r)
  } else {

    if (inherits(Nyx, "bida_sparse_array")) {

      y  <- Nyx$index%%r
      x  <- Nyx$index%/%r
      s  <- get_parts(partition)[x+1]-1
      us <- unique(s)

      group  <- y + r*s                               # enumerate observed outcomes (y, s)
      ugroup <- rep(seq_len(r)-1, length(us)) + r*s   # enumerate possible outcomes (y, s) for observed s
      Nys  <- rowsum_fast(Nyx, group, ugroup)         # zero for levels of group not in ugroup

      ays  <- matrix(Nys, length(us), r, byrow = T) + lengths(partition[us+1])*ess/(r*q)
      tmp  <- ays/rowSums(ays)
      bida_sparse_array(t(tmp), ugroup, c(r, length(partition)), default = 1/r)
    } else if (inherits(Nyx, "array")){

      Nyx <- matrix(Nyx, q, r, byrow = T)
      Nys <- rowsum_fast(Nyx, get_parts(partition), seq_along(partition))
      ays <- Nys + lengths(partition)*ess/(r*q)
      t(ays/rowSums(ays))
    } else stop("Nyx must be an array or an bida_sparse_array.")
  }
}

backdoor_mean_bdeu <- function(Nyxz, ess = 1, partition = NULL, kx = 1) {
  dims <- dim(Nyxz)
  n <- length(dims)
  if (n == 1) {
    p <- rep(mean_bdeu(Nyxz, ess), times = kx)
    dim(p) <- c(dims, kx)
    return(p)
  } else if (n  == 2) {
    mean_bdeu(Nyxz, ess)
  } else if (is.null(partition)) {
    Nxz <- colSums(Nyxz)
    px.z <- mean_bdeu(Nxz+ess/prod(dims[-1]))
    colSums( (Nyxz + ess/prod(dims))/rep(px.z, each = dims[1]) )/sum(Nxz)
  } else {

    Nxz <- colSums(Nyxz)
    Nz  <- colSums(Nxz)
    py.s <- mean_bdeu(Nyxz, ess = ess, partition = partition)

    if (inherits(Nyxz, "bida_sparse_array")) {
      xz  <- Nyxz$index%/%dims[1]
      x   <- xz%%dims[2]
      z   <- xz%/%dims[3]
      s   <- get_parts(partition)[xz+1]

      xs  <- x + dims[2]*(s-1) # enumerate outcomes (x, s)
      uxs <- unique(xs)
      Nz.xs <- rowsum_fast(Nz$value[match(z, Nz$index)], xs, uxs)

      Nys <- bida_sparse_array(rowsum_fast(Nyx, ys, uys), uys, c(r, length(partition)))
      a0ys <- lengths(partition[unique(s)])*ess/(r*q)

    }
    Nz <- colSums(Nyxz, dims = 2)
    py.s <- mean_bdeu(Nyx, ess, partition)
  }
}

rbackdoordist_cat <- function(n, counts) {

  N <- sum(counts$value)
  Nzxy <- counts$value
  zxy  <- counts$index
  dims <- counts$dim

  # cardinality of (sets of) variables
  kz  <- dims[1]
  ky  <- dims[3]
  kx  <- dims[2]
  kxy <- ky*kx

  # prior hyperparams
  a0z <- ess/kz
  a0y.xz <- ess/(kz*kz*ky)

  # map index to observed levels of z, and (x, y) for computation of grouped counts
  xy  <- (zxy%/%kz) +1
  z   <- zxy%%kz +1
  z   <- match(z, unique(z))   # re-index as 1, 2, ..., length(unique(z))
  uz  <- unique(z)


  # sample the backdoor sum over observed levels of z
  # - If z_1, ..., z_j are observed levels, then P(z_1, ..., z_j| sum P(z_i) over i from j+1 to kz) is Dirichlet
  az <- tabulate(z, length(uz)) + a0z                      # posterior hyperparams for observed levels

  # sample the backdoor sum over observed levels of z
  pdo <- array(0, dim = c(n, kx, ky))
  pz <- rDirichlet(n, az, kz)
  seq_kx <- seq_len(kx)
  tmp <- matrix(rep(a0y.xz, kx*ky), nrow = kx, byrow = T)
  for (zz in seq_along(uz)) {
    index <- z == zz
    ay.xz <- replace(tmp, xy[index], Nzxy[index] + a0y.xz[1])
    for (xx in seq_kx) {
      pdo[,xx,] <- pdo[, xx, ] + rDirichlet(n, ay.xz[xx, ], ky)*pz[, zz]
    }
  }

  # sample aggregated prob of z beeing (unobserved)
  pobs <- rbeta(n, N+a0z*length(uz), (kz-length(uz))*a0z)  # prob that z takes one of observed levels
  pdo  <- pdo*1/(1-pobs)

  # add the backdoor sum over unobserved levels of z

  sample_prior <- function(k, cumprob, th = .999) {
    iter <- k
    while (cumprob < th && iter > 1) {
      phi <- rbeta(1, a0z, iter*a0z)
      pz  <- (1-cumprob)*phi

      # add sample to backdoor sum
      py.xz <- rDirichlet(kx, ay.xz)
      pdo <- pdo + py.xz*pz

      cumprob <- cumprob + pz
      iter <- iter-1
    }
    # last iter
    pz  <- (1-cumprob)
    py.xz <- rDirichlet(kx, ay.xz)
    pdo + py.xz*pz

  }

  ay.xz <- rep(a0y.xz, ky)
  tmp <- vapply(pobs,
                function(cumprob) sample_prior(kz-length(uz), cumprob, th = .999),
                matrix(numeric(), kx, ky))

  pdo + tmp

}




# Function to sample from the Dirichlet distribution using conditional Beta distributions
sample_conditional_dirichlet <- function(alpha) {
  n <- length(alpha)
  samples <- numeric(n)
  alpha_sum <- sum(alpha)

  for (i in 1:(n-1)) {
    # Sample from the conditional Beta distribution
    samples[i] <- rbeta(1, alpha[i], alpha_sum - alpha[i])
    # Update alpha sum and remaining alpha
    alpha_sum <- alpha_sum - samples[i]
    alpha <- alpha - samples[i]
  }

  # The last element
  samples[n] <- 1 - sum(samples)

  return(samples)
}

# Example usage
alpha <- c(0.5, 1.5, 2.0)
dirichlet_sample <- sample_conditional_dirichlet(alpha)
print(dirichlet_sample)
