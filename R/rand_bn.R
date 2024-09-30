
#' Draw a random Bayesian network
#'
#' Given a DAG, draw a distribution over the DAG and store the BN as a `bnlearn::bn.fit` object.
#'
#' @keywords internal
#' @param dag (integer matrix) adjacency matrix.
#' @param type (character) type of distribution.
#' @param ...
#'
#' @return a `bnlearn::bn.fit` object
#' @keywords internal
#'
#' @examples
#'
#' n <- 3
#' dag <- matrix(0, n, n)
#' dag[upper.tri(dag)] <- 1
#'
#' # categorical distribution with local structure
#' nlev <- rep(2, ncol(dag))
#' partitions <- list(NULL, NULL, X3 = list(0:2, 1))
#' rand_bn(dag, "cat", alpha = 1, nlev = rep(2, 3), partitions)
#'
#'
#' partitions <- rand_partitions_cat(dag, nlev, 1, maxdepth = 2)
#' rand_bn(dag, "cat", alpha = 1, nlev, partitions = partitions)
rand_bn <- function(dag, type = "cat", ...) {
  n <- ncol(dag)

  # init bn-learn object
  varnames <- colnames(dag)
  if (is.null(varnames)) varnames <- paste0("X", seq_len(n))
  g <- bnlearn::empty.graph(varnames)
  bnlearn::amat(g) <- dag

  # draw random distribution
  dist <- switch(type,
                 "cat" = rand_dist_cat(dag, ...))

  # return bn.fit object
  bnlearn::custom.fit(g, dist)
}

rand_partitions_cat <- function(dag, nlev, splitprob, maxdepth) {
  n <- ncol(dag)
  partitions <- vector("list", n)
  for (i in seq_len(n)) {
    pa <- which(dag[, i] == 1)
    if (length(pa) > 1) {
      partitions[[i]] <- rand_partition(nlev[pa], splitprob, method = "tree", maxdepth = maxdepth)
    }
  }
  return(partitions)
}

rand_dist_cat <- function(dag, alpha, nlev, partitions = NULL) {
  n <- ncol(dag)
  dist <- vector("list", n)
  varnames <- colnames(dag)
  if (is.null(varnames)) varnames <- paste0("X", seq_len(n))

  names(nlev) <- varnames # for dimnames
  for (i in seq_len(n)) {
    pa <- which(dag[, i] == 1)
    r <- nlev[i]
    q <- prod(nlev[pa])
    partition <- partitions[[i]]

    if (is.null(partition)) {
      p <- vapply(integer(q),
                  function(x) rDirichlet(1, rep(alpha, r), r), numeric(r))
    } else {

      stopifnot(length(unlist(partition)) == q)
      parts <- get_parts(partition)
      p <- vapply(seq_along(partition),
                  function(x) rDirichlet(1, rep(alpha, r), r), numeric(r))[, parts]

    }

    dims <- nlev[c(i, pa)]
    dimnames <- lapply(dims-1, seq.int, from = 0)
    dist[[i]] <- array(p, dims, dimnames)
  }

  names(dist) <- varnames
  return(dist)
}
