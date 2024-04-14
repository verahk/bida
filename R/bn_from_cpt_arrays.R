

#' Construct a bnlearn-network from a list of CPTs
#'
#' @param cpts a list of arrays representing the CPT of each node.
#' @return An object of class `bn.fit.dnet`
#' @keywords internal
bn_from_cpt_arrays <- function(cpts) {
  varnames <- names(cpts)
  bn <- bnlearn::empty.graph(varnames)
  bnlearn::amat(bn) <- dag_from_cpt_arrays(cpts)
  bnlearn::custom.fit(bn, cpts)
}
