#' Draw a random distribution over a DAG
#'
#' Draw a random distribution over a DAG.
#' The distribution is represented in a format matching [bnlearn::bn.fit] class,
#' and is used to construct such objects.
#'
#' @param dag (integer matrix) an adjacency matrix
#' @param type (character) type of distribution
#' @param ... additional arguments for the function producing each local distribution.
#'  See [rand_cpt()] for categorical distributions.
#'
#' @return a list with distributions
#' @export
#'
#' @seealso [rand_cpt()]
#' @examples
#'
#' dag  <- matrix(0, 3, 3, dimnames = list(c("Z", "X", "Y"), c("Z", "X", "Y")))
#' dag[upper.tri(dag)] <- 1
#' rand_dist(dag, "cat", nlev = rep(3, ncol(dag)))
#'
rand_dist <- function(dag, type = "cat", ...) {
  type <- match.arg(type, c("categorical"))

  n <- ncol(dag)
  seqn <- seq_len(n)
  varnames <- colnames(dag)
  if (is.null(varnames)) varnames <- paste0("X", seqn)

  args <- list(...)
  dist <- setNames(vector("list", n), varnames)
  for (j in seqn) {
    parentnodes <- seqn[dag[, j] == 1]
    dist[[j]] <- switch(type,
                        "categorical" = rand_cpt(dim = args$nlev[c(j, parentnodes)],
                                                 scope = varnames[c(j, parentnodes)],
                                                 ...))
  }
  dist
}
