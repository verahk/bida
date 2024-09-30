

#' Class: `bida_bdeu`
#'
#' Represent a distribution over a conditional probability table (CPT)
#'
#' @param data a matrix with categorical data
#' @param j (integer) column position of outcome node
#' @param parentnodes column position(s) of parent node(s)
#' @param ess (integer) imaginary sample size
#' @param nlev (integer vector) cardinality of each variable
#' @param partition a list with partition over the joint outcomes of `parentnodes`.
#' @details
#'
#' The `bida_bdeu` class represent a (posterior) BDeu distribution over the
#' parameters of a conditional probability table (CPT), i.e. a distribution over
#' the distributions `P(X_j|Pa(X_k))`.
#'
#' Methods...
#' - `get_dim`: the dimension of the CPT, i.e. `nlev[c(j, parentnodes)]`.
#' - `posterior_mean`: returns a CPT with the posteriorior mean values.
#'    For `reduced = TRUE` or `is.null(obj$partition)`, the mean values are
#'    stored in an array with `nlev[c(j, parentnodes)]` dimensions. Otherwise,
#'    an `r-by-length(obj$partntion)` matrix.
#' - `posterior_sample` returns an sample of CPTs.
#' - `backdoor_mean`
#' - `backdoor_sample`
#' @return
#' An object of class `bida_bdeu` is a list that contains:
#' - `counts` an array with counts for each joint outcome.
#' - `partition` a list representing a partition over the parent outcomes.
#' - `ess` the imaginary sample size.
#' @export
#' @examples
#'
#' nlev <- 2:4
#' lev  <- lapply(nlev-1, seq.int, from = 0)
#' data <- expand_grid_fast(lev)
#'
#' j <- 1
#' parentnodes <- 2:length(nlev)
#' bdeu <- bida_bdeu(data, j, parentnodes, ess = 1, nlev)
#'
#' # optimize partition of parent space
#' opt <- optimize_bdeu(bdeu, method = "tree", regular = T)
#' opt$partition
#' bdeu_part <- replace(bdeu, "partition", list(opt$partition))
#'
#' # compute score
#' score_bdeu(bdeu)
#' score_bdeu(bdeu_part)
#'
#' # update hyperparams
#' update_bdeu(bdeu)
#' update_bdeu(bdeu_part)
#'
#' # permute parent dimension
#' perm <- c(1, 3, 2)
#' bdeu_perm <- aperm(bdeu_part, perm)
#' stopifnot(all.equal(as.array(bdeu_perm$counts), array(1, nlev[perm])))
#' bdeu_perm
#'
#' # posterior mean
#' posterior_mean(bdeu)
#' posterior_mean(bdeu_part, reduced = T)  # means of reduced CPT
#' posterior_mean(bdeu_part, reduced = F)  # means of full CPT
#'
#' # with local tree-structure defined as string
#' opt <- optimize_partition_from_data(data, j, parentnodes, ess = 1, nlev = nlev, method = "tree")
#' bdeu$partition <- opt
bida_bdeu <- function(data, j, parentnodes, ess, nlev, partition = NULL) {
  if (is.null(colnames(data))) colnames(data) <- paste0("X", seq_len(ncol(data)))
  subset <- c(j, parentnodes)
  counts <- counts_from_data_matrix(data[, subset, drop = FALSE], nlev[subset], sparse = T)
  new_bida_bdeu(colnames(data)[1], colnames(data)[parentnodes], counts, ess, partition)
}

#' @noRd
new_bida_bdeu <- function(node, parents, counts, ess, partition, scope = NULL) {
  structure(list(node = node,
                 parents = parents,
                 counts = counts,
                 partition = partition,
                 ess = ess),
            class = "bida_bdeu")
}

#' S3-methods ----
#' @rdname bida_bdeu
#' @export
aperm.bida_bdeu <- function(obj, perm) {

  if (is.null(obj$partition)) {
    obj$counts <- aperm.bida_sparse_array(obj$counts, perm)
    return(obj)
  }

  stopifnot(perm[1] == 1)

  # create matrix with parent outcomes
  dims_pa <- obj$counts$dim[-1]
  perm_pa <- perm[-1]-1
  lev_pa  <- lapply(dims_pa-1, seq.int, from = 0)

  # order current parent config by permuted
  conf_pa <- expand_grid_fast(lev_pa[perm_pa])
  stride_pa <- c(1, cumprod(dims_pa[-length(dims_pa)]))[perm_pa]
  new_indx_pa <- conf_pa%*%stride_pa

  # adjust counts and partition
  dims <- obj$counts$dim
  obj$counts$dim <- dims[perm]
  indx_y  <- obj$counts$index%%dims[1]
  indx_pa <- obj$counts$index%/%dims[1]
  obj$counts$index <- indx_y + dims[1]*(match(indx_pa, new_indx_pa)-1)
  obj$partition <- relist(match(unlist(obj$partition), new_indx_pa)-1, obj$partition)

  return(obj)
}


#' @rdname bida_bdeu
#' @param reduced (bolean) if TRUE (default), the reduced CPT is returned.
#' @export
posterior_mean.bida_bdeu <- function(obj, reduced = TRUE) {
  dims <- dim(obj$counts)
  if (length(dims) == 1) {
    p <- (obj$counts + obj$ess/dims)/sum(obj$counts + ess)
  } else {
    r <- dims[1]
    q <- prod(dims[-1])
    if (is.null(obj$partition)) {
      p <- (obj$counts+obj$ess/(r*q))/rep(colSums(obj$counts)+obj$ess/q, each = r)
    } else {
      alpha <- update_bdeu(obj)

      p <- alpha/rep(colSums(alpha), each = dims[1])

      if (reduced) {
        return(p)
      } else {
        parts <- attr(alpha, "parts")
        return(array(p[ , parts], dim(obj$counts)))
      }
    }
  }
  return(as.array(p))
}

#' @rdname bida_bdeu
#' @returns
#' - `posterior_sample.bida_bdeu` returns an array of dimensions `c(dim(x), n)`,
#'   a `n` sized sample of CPTs drawn from the posterior distribution.
#' @export
posterior_sample.bida_bdeu <- function(obj, n, reduced = TRUE) {

  alpha <- update_bdeu(obj)

  if (length(dim(alpha)) == 1) {
    t(rDirichlet(n, alpha))
  } else {

    r <- dim(alpha)[1]
    tmp <- seq_along(dim(alpha))
    p <- array(apply(alpha, tmp[-1], rDirichlet, n = n, k = r), c(n, dim(alpha)))
    p <- aperm(p, c(tmp+1, 1))

    if (reduced || is.null(obj$partition)) {
      return(p)
    } else {
      # replicate and return as array
      parts <- attr(alpha, "parts")
      array(p[ , parts, ], c(dim(obj), n))
    }
  }
}

#' @rdname bida_bdeu
#' @param nlevx (integer) cardinality of cause variable, for replicating vectors
#'  with marginal probabilities for each intervention level.
#' @export
backdoor_mean.bida_bdeu <- function(obj, nlevx) {
  dims <- dim(obj$counts)
  ndims <- length(dims)
  if (ndims == 1) {
    # no adjustment
    array(posterior_mean(obj), c(dims, nlevx))
  } else if (ndims == 2) {
    posterior_mean(obj)
  } else {

    # apply backdoor-formula
    ess <- obj$ess

    # compute conditional means
    py.xz <- posterior_mean(obj, reduced = F)

    # compute counts over adjustment sets
    kyx <- prod(dims[1:2])
    kz  <- prod(dims[-c(1:2)])
    z  <- obj$counts$index%/%kyx
    Nz <- rowsum_fast(obj$counts$value, z, seq_len(kz)-1)

    # sum out z
    N <- sum(Nz)
    rowSums(py.xz*rep(Nz + ess/length(Nz), each = kyx), dims = 2)/(N+ess)
  }
}

#' @rdname bida_bdeu
#' @export
backdoor_sample.bida_bdeu <- function(obj, n, nlevx, digits = 16) {
  dims = dim(obj$counts)
  ndims = length(dims)
  if (ndims == 1) {
    tmp <- round(posterior_sample(obj, n), digits)
    array(tmp[rep(seq_len(dims), 2), ], c(dims, nlevx, n))
  } else if (ndims == 2) {
    round(posterior_sample(obj, n), digits)
  } else {

    ess <- obj$ess
    k <- prod(dims)
    kyx <- prod(dims[1:2])
    kz  <- k/kyx

    # compute counts over adjustment sets
    # - rowsum_fast(v, group, ugroup) returns zero for elements of ugroup not found in group
    z  <- obj$counts$index%/%kyx
    az <- rowsum_fast(obj$counts$value, z, seq_len(kz)-1) + ess/kz

    # sample distributions over adjustment set
    pz <- t(round(rDirichlet(n, az, length(az)), digits))

    # collapse dimensions of adjustment set, for aperm below
    obj$counts$dim <- c(dims[1:2], kz)

    # sample conditional means
    py.xz <- round(posterior_sample(obj, n, reduced = F), digits)

    colSums(aperm(py.xz*rep(pz, each = kyx), c(3, 1, 2, 4)))

  }
}


#' @rdname bida_bdeu
#' @param method name of optimization procedure, see `optimize_partition`.
#' @export
optimize_bdeu <- function(obj, method, levels = NULL, ...) {

  dims <- get_dim(obj$counts)
  if (length(dims) < 2) return(NULL)
  counts <- as.array(obj$counts)
  ess <- obj$ess

  # represent counts in q-by-r matrix
  r <- dims[1]
  q <- prod(dims[-1])
  counts <- matrix(counts, q, r, byrow = T)

  # optimize partitition
  if (is.null(levels)) levels <- lapply(dims[-1]-1, seq.int, from = 0)
  optimize_partition(counts, levels, ess, method, ...)
}


#' @rdname bida_bdeu
#' @details
#' - `update_bdeu`: updates BDeu-hyper parameters.
#'    If `is.null(obj$partition)`, returns a full (non-sparse) array with the same
#'    dimensions as `obj$counts`. Otherwise a `r-by-nparts` matrix, where `r` is
#'    the cardinality of the outcome and `nparts` the size of the partition of the
#'    parent space.
update_bdeu <- function(obj, parent_config = NULL) {
  if (is.null(parent_config)) {
    counts <- as.array(obj$counts)
    ess <- obj$ess
    partition <- obj$partition

    if (is.null(partition)) {
      ess/length(counts) + counts
    } else {

      dims <- dim(obj$counts)
      r <- dims[1]         # cardinality of outcome variable
      q <- prod(dims[-1])  # cardinality of parent variables

      # represent counts in q-by-r matrix for rowsum
      counts <- matrix(as.array(obj$counts), q, r, byrow = T)

      # compute posterior hyperparameters
      if (inherits(obj$partition, "partition")) {
        partsize <- obj$partition$size
        tmp <- expand_grid_fast(lapply(dims[-1]-1, seq.int, from = 0))
        colnames(tmp) <- obj$parents
        parts <- predict(obj$partition, tmp)
      } else {
        partsize <- lengths(partition)
        parts <- get_parts(partition)
      }
      alpha <- ess/(r*q)*partsize + rowsum_fast(counts, parts, seq_along(partsize))
      attr(alpha, "parts") <- parts
      t(alpha)
    }
  } else {
    stop()
  }
}

#' @rdname bida_bdeu
#' @details
#' - `score_bdeu` computes the Bayesian score of the CPT
#' @export
score_bdeu <- function(obj) {

  dims <- dim(obj$counts)
  ess <- obj$ess
  r <- dims[1]
  q <- prod(dims[-1])

  # represent counts in q-by-r matrix for scoring functions
  counts <- as.array(obj$counts)
  counts <- matrix(counts, q, r, byrow = T)

  if (is.null(obj$partition)) {
    famscore_bdeu(counts, ess, r, q)
  } else {

    # aggregate counts
    parts <- get_parts(obj$partition)
    agg <- rowsum_fast(counts, parts, seq_along(obj$partition))

    # compute family score
    sum(famscore_bdeu_byrow(agg, ess, r, q, lengths(obj$partition)))
  }
}
