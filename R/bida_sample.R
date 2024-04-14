
#' Main function for computing the BIDA posterior using a sample of DAGs.
#'
#' @inheritParams bida_exact
#' @param dags list of DAGs (adjacency matrices, obtained e.g. with BiDAG::partitionMCMC)
#' @param support numeric vector of same length as [dags]. Support of each DAG. Default to [1/length(dags)].
#' @param adjset type of adjustment set to compute in each DAG.
#'
#' @return an object of class [bida]
#' @seealso [bida]
#' @export

bida_sample <- function(data, dags, support = rep(1/length(dags), length(dags)), type = "categorical", adjset = "pa", params = NULL){


  if (tolower(type) == "categorical"){

    if (is.null(params)) {
      params <- list(nlev = apply(data, 2, max)+1,
                     ess  = 1,
                     maxconf = 2**16)
    }

    if (params$maxconf < Inf) {
      checksize <- function(x, y, z) length(z) <= 1 || prod(params$nlev[z]) < params$maxconf
    } else {
      checksize <- NULL
    }
  }

  # compute support adjustment sets
  if (adjset == "local") {
    ps <- parent_support_from_dags(dags, support, checksize = checksize)
  } else {
    ps <- adjsets_support_from_dags(adjset, dags, support, checksize = checksize)[[adjset]]
  }


  bida(ps, data, type, params)
}








