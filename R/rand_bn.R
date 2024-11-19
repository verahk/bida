
#' Draw a random Bayesian network
#'
#' Draw a random DAG and a distribution over the DAG, and store the BN as a `bnlearn::bn.fit` object.
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
#'
#' # categorical ---
#'
#' # draw random CPT with local structure
#' nlev <- c(2, 2, 2)
#' cpt <- rand_cpt(nlev, c("y", "x", "z"), local_structure = "tree", prob = 1, maxdepth = .5)
#'
#' # draw random DAG
#' set.seed(007)
#' dag <- rand_dag(3, 2)
#' dag
#'
#' # draw random distrib over DAG
#' cpts <- rand_dist(
#'   dag,
#'   "cat",
#'   nlev = nlev,
#'   local_structure = "tree",
#'   prob = 1,
#'   maxdepth = .5
#' )
#'
#' # create bn object
#' bn <- custom_bn(dag, cpts)
#'
#' #
#' set.seed(007)
#' bn2 <- rand_bn(3, 2, "cat", nlev = nlev, local_structure = "tree", prob = 1, maxdepth = .5)
#' stopifnot(all.equal(bn, bn2))
#'

rand_bn <- function(n, d, type = "cat", attr_to_keep = character(0), ...) {
  dag <- rand_dag(n, d)
  dist <- rand_dist(dag, type, ...)
  bn <- custom_bn(dag, dist)

  for (a in attr_to_keep) {
    attr(bn, a) <- lapply(dist, attr, which = a)
  }

  return(bn)
}

custom_bn <- function(dag, dist) {

  # init bn-learn object
  varnames <- colnames(dag)
  if (is.null(varnames)) varnames <- paste0("X", seq_len(ncol(dag)))
  g <- bnlearn::empty.graph(varnames)
  bnlearn::amat(g) <- dag

  # create bn.fit object
  bnlearn::custom.fit(g, dist)

}




