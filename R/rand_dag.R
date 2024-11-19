#' Draw a random adjacency matrix
#'
#' Wrapper around [pcalg::randDAG] returning a adjacency matrix rather than a graphNEL object
#' @param n (integer) number of nodes
#' @param d (integer) expected number of neighbours
#' @param weighted (boolean) wheter each edge should be assigned weights. Default is FALSE.
#'
#' @return a n-by-n adjacency matrix
rand_dag <- function(n, d, names = paste0("X", seq_len(n)), weighted = FALSE) {
  dag <- as(pcalg::randDAG(n, d, weighted = weighted), "matrix")
  colnames(dag) <- rownames(dag) <- names
  dag
}

