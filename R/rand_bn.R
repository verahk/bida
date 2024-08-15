
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

rand_dist_cat <- function(dag, alpha, nlev, partitions = NULL) {
  n <- ncol(dag)
  dist <- vector("list", n)
  varnames <- colnames(dag)
  if (is.null(varnames)) varnames <- paste0("X", seq_len(n))

  names(nlev) <- varnames # for dimnames
  for (i in seq_len(n)) {
    pa <- which(dag[, i] == 1)
    dims <- nlev[c(i, pa)]
    dimnames <- lapply(dims-1, seq.int, from = 0)

    bdeu <- new_bida_bdeu(new_bida_sparse_array(0, 0, dims),
                          ess = alpha*prod(dims),
                          partition = NULL)
    if (!(is.null(partitions) && is.null(partitions[[i]]))) {
      bdeu$partition <- partitions[[i]]
    }
    tmp <- sample_bdeu(1, bdeu, reduced = FALSE)
    dist[[i]] <- array(tmp, dims, dimnames)
  }

  names(dist) <- varnames
  return(dist)
}
