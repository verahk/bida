

#' List CPTs of bn network
#'
#' The CPTs associated with each `node` in a `bn.fit.dnet` object is available
#' `node$prob`.
#' To compute exact queries with [cpquery_from_cpt_arrays],
#' the scope of the CPTs must be given by the `dimnames` attribute.
#'
#' @param bn a `bn.fit.dnet` object (see [bnlearn::bn.fit])
#' @details
#'  - `cpt_arrays_from_bn`
#' @return a list with cpts stored as arrays
#' @keywords internal
cpt_arrays_from_bn <- function(bn) {
  cpts <- vector("list", length(bn))
  names(cpts) <- names(bn)
  for (i in seq_along(bn)) {
    node  <- bn[[i]]
    nlev  <- dim(node$prob)
    scope <- c(node$node, node$parents)
    cpts[[i]] <- array(node$prob,
                       nlev,
                       setNames(dimnames(node$prob), scope))
  }
  return(cpts)
}
