
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

rand_dist_cat <- function(dag, alpha, nlev, levels = lapply(nlev-1, seq.int, from = 0), partitions = NULL) {
  n <- ncol(dag)
  dist <- vector("list", n)
  varnames <- colnames(dag)
  if (is.null(varnames)) varnames <- paste0("X", seq_len(n))

  names(nlev) <- varnames # for dimnames
  for (i in seq_len(n)) {
    pa <- which(dag[, i] == 1)

    if (is.null(partitions) || is.null(partitions[[i]])) {
      p <- vapply(1:prod(nlev),
                  function(x) rDirichlet(1, rep(alpha, r), r), numeric(r))
    } else {
      parts <- get_parts(partitions[[i]])
      p <- vapply(seq_along(partitions[[i]]),
                  function(x) rDirichlet(1, rep(alpha, r), r), numeric(r))[, parts]

    }

    dims <- nlev[c(i, pa)]
    dimnames <- lapply(dims-1, seq.int, from = 0)
    dist[[i]] <- array(p, dims, dimnames)
  }

  names(dist) <- varnames
  return(dist)
}
