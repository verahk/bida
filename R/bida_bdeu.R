

#' Class: `bida_bdeu`
#'
#' Parameters of a Dirichlet distribution over a conditional probability table (CPT),
#' assuming a BDeu-prior.
#'
#'
#' @name bida_bdeu
#' @param counts
#' @param ess
#' @param partition
#'
#' @return
#' @export
#'
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
#' mean_bdeu(bdeu)
#' mean_bdeu(bdeu_part, reduced = T)  # means of reduced CPT
#' mean_bdeu(bdeu_part, reduced = F)  # means of full CPT
#'
#' # posterior sample
#' dim(sample_bdeu(1000, bdeu_part))
#' dim(sample_bdeu(1000, bdeu_part, reduced = F))
#'
#'
#'

new_bida_bdeu <- function(counts, ess, partition) {
  structure(list(counts = counts,
                 partition = partition,
                 ess = ess),
            class = "bida_bdeu")
}

bida_bdeu <- function(data, j, parentnodes, ess, nlev, partition = NULL) {
  subset <- c(j, parentnodes)
  counts <- counts_from_data_matrix(data[, subset, drop = FALSE], nlev[subset], sparse = T)
  new_bida_bdeu(counts, ess, partition)
}


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


#' @rdname bdeu_par
#' @details
#' - `update_bdeu`: updates BDeu-hyper parameters.
#'    If `is.null(obj$partition)`, returns a full (non-sparse) array with the same
#'    dimensions as `obj$counts`. Otherwise a `r-by-nparts` matrix, where `r` is
#'    the cardinality of the outcome and `nparts` the size of the partition of the
#'    parent space.
#' @export
update_bdeu <- function(obj, parent_config = NULL) {
  if (is.null(parent_config)) {

    counts <- as.array(obj$counts)
    ess <- obj$ess
    partition <- obj$partition

    if (is.null(partition)) {
      ess/length(counts) + counts
    } else {
      r <- dim(counts)[1]         # cardinality of outcome variable
      q <- prod(dim(counts)[-1])  # cardinality of parent variables

      # represent counts in q-by-r matrix for rowsum
      counts <- matrix(counts, q, r, byrow = T)

      # compute posterior hyperparameters
      parts <- get_parts(partition)
      alpha <- ess/(r*q)*lengths(partition) + rowsum_fast(counts, parts, seq_along(partition))
      attr(alpha, "parts") <- parts
      t(alpha)
    }
  } else {
    stop()
  }
}

#' @rdname bdeu_par
#' @details
#' - `score_bdeu` computes the Bayesian score of the CPT
#' @export
score_bdeu <- function(obj) {

  dims <- get_dim(obj)
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

#' @rdname bdeu_par
#' @details
#' - `mean_bdeu` computes the (posterior) mean of the CPT
#' @export
mean_bdeu <- function(obj, reduced = TRUE) {
  alpha <- p <- update_bdeu(obj)
  if (length(dim(alpha)) < 2) {
    return(alpha/sum(alpha))
  } else {
    dims <- get_dim(obj)
    p <- alpha/rep(colSums(alpha), each = dims[1])
  }

  if (reduced || is.null(obj$partition)) {
    return(p)
  } else {
    parts <- attr(alpha, "parts")
    array(p[ , parts], dims)
  }
}

#' @rdname bdeu_par
#' @details
#' `sample_bdeu` samples CPTs from the posterior
#' @export
sample_bdeu <- function(size, obj, reduced = TRUE) {

  alpha <- update_bdeu(obj)

  if (is.null(dim(alpha))) {
    t(rDirichlet(n, alpha))
  } else {
    dims <- get_dim(obj)

    # sample joint probabilities
    p <- array(t(bida:::rDirichlet(size, alpha)), c(dim(alpha), size))

    # compute conditional probabilities
    p <- p/rep(colSums(p), each = dims[1])

    if (reduced || is.null(obj$partition)) {
      return(p)
    } else {
      # replicate and return as array
      parts <- attr(alpha, "parts")
      array(p[ , parts, ], c(dims, size))
    }
  }
}



