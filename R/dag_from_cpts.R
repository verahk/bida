dag_from_cpts <- function(cpts) {

  # list scope of each CPT
  scopes <- lapply(cpts, function(x) names(dimnames(x)))
  lens   <- lengths(scopes)
  stopifnot(all(lens>0))  # all CPTs do not have named dimnames attribute

  # init adjacency matrix
  n   <- length(cpts)
  dag <- matrix(0, n, n)
  colnames(dag) <- rownames(dag) <- sapply(scopes, "[[", 1)

  # set edges
  for (i in seq_len(n)[lens > 1]) {
    pa <- scopes[[i]][-1]
    dag[pa, i] <- 1
  }
  return(dag)
}
