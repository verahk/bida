#' Draw a random adjacency matrix
#'
#' Wrapper around [pcalg::randDAG] returning a adjacency matrix rather than a graphNEL object
#' @param n (integer) number of nodes
#' @param d (integer) expected number of neighbours
#' @param weighted (boolean) wheter each edge should be assigned weights. Default is FALSE.
#'
#' @return a n-by-n adjacency matrix
#' @keywords internal
randDAG <- function(n, d, weighted = FALSE) as(pcalg::randDAG(n, d, weighted = weighted), "matrix")

