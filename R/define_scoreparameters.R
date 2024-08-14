#' User defined scoreparameters for `BiDAG`
#'
#' Define a [BiDAG::scoreparameters] object with user defined score.
#'
#' @param data (matrix)
#'  each column is assumed to be outcomes of a categorical variable, with
#'  possible outcomes `0, ..., K-1.`
#' @param scoretype
#'  parameters controlling the local-structure-optimization routine. See [optimize_partition()].
#' @param lookup (environment)
#'  an environment in which the optimized structures are saved, initiated by rlang:::new_environment().
#'  If `NULL` the structures are not  saved.
#' @details
#'  For categorical data, the list of parameters `par` have to include:
#'  - `ess`: equivalent sample size
#'  - `edgepf`: a factor `-log(edgepf)*length(parentnodes)` is added to the marginal likelihood score.
#'  - `nlev`: a vector with the cardinalities of all variables
#' @return an object of class [BiDAG::scoreparameters]
#' @export
#'
define_scoreparameters <- function(data, scoretype, par = NULL, lookup = NULL) {

  if (scoretype == "bdecat") {
    # generate data.frame, removing unobserved levels in data, as required by BiDAG
    df <- data.frame(apply(data, 2, factor, exclude = NULL, simplify = FALSE))
    scorepar <- BiDAG::scoreparameters("bdecat",
                                       data = df,
                                       bdecatpar = list(chi = par$ess, edgepf = par$edgepf))

  } else if (FALSE || scoretype %in% c("tree", "ldag", "part")) {

    scorepar <- BiDAG::scoreparameters("usr",
                                       data = data,
                                       usrpar = list(pctesttype = "bdecat",
                                                     chi = par$ess,
                                                     edgepf = par$edgepf))

    # add cardinality of variables
    scorepar$Cvec <- par$nlev  # CVec is used in PC-routine of BiDAG-functions (?)

    # additional parameters for user-specified score params
    scorepar$levels <- lapply(par$nlev-1, seq.int, from = 0)
    scorepar$local_struct <- scoretype
    scorepar$regular <- regular

    # define function for computing scores
    if (is.null(lookup)) {
      usrDAGcorescore <- function(j, parentnodes, n, scorepar) {
        # define function
      }
    } else {
      # define environment in which to store optimized local structures
      tmp <- list()
      assign(local_struct, tmp, lookup)
      scorepar$lookup <- lookup

      usrDAGcorescore <- function(j, parentnodes, n, scorepar) {
       # define function
      }

    }



    # assign function to name-space of BiDAG package
    assignInNamespace("usrDAGcorescore", usrDAGcorescore, ns = "BiDAG")
  } else {
    stop(sprintf("Invalid argument for `scoretype`: %s", scoretype))
  }
  return(scorepar)
}
