#' User defined scoreparameters for `BiDAG`
#'
#' Define a [BiDAG::scoreparameters] object with user defined score.
#'
#' @param data (matrix)
#'  each column is assumed to be outcomes of a categorical variable, with
#'  possible outcomes `0, ..., K-1.`
#' @param scoretype
#'  parameters controlling the local-structure-optimization routine. See [optimize_partition()].
#' @param `par`
#'  For categorical data, the list of parameters `par` have to include:
#'  - `ess`: equivalent sample size
#'  - `edgepf`: a factor `-log(edgepf)*length(parentnodes)` is added to the marginal likelihood score.
#'  - `nlev`: a vector with the cardinalities of all variables.
#' @param lookup (environment)
#'  an environment in which the optimized structures are saved, initiated by rlang:::new_environment().
#'  If `NULL` the structures are not  saved.
#' @return an object of class [BiDAG::scoreparameters]
#' @export
#'
#' @examples
#'
#' lookup <- rlang:::new_environment() # environment for storing scores and CPTs
#' bn  <- readRDS("./data/asia.rds")
#' data <- sample_data_from_bn(bn, 10)
#'
#' par <- list(nlev = apply(data, 2, max)+1, ess = 1, edgepf = 2, regular = T)
#'
#'
#' j <- 1
#' parentnodes <- 2
#'
#' # no partitioning - use BiDAG-score
#' scorepar_bdecat <- define_scoreparameters(data, "bdecat", par, lookup = lookup)
#' score <- BiDAG:::DAGcorescore(j, parentnodes, ncol(data), scorepar_bdecat)
#' score
#'
#' # user-specified function is assigned to BiDAG-namespace
#' scorepar_tree <- define_scoreparameters(data, "tree", par, lookup = lookup)
#' BiDAG:::usrDAGcorescore

#' # with less than 2 parents, the partition space is not partitioned
#' stopifnot(score == BiDAG:::DAGcorescore(j, parentnodes, ncol(data), scorepar_tree))
#'
#' # user-specified functions writes to lookup table
#' ls.str(lookup)    # empty
#' BiDAG:::usrDAGcorescore(1, 2:3, ncol(data), scorepar_tree)
#' ls.str(lookup)    # updated
#' lookup$tree[["1.2.3"]]
define_scoreparameters <- function(data, scoretype, par = NULL, lookup = NULL) {

  if (scoretype == "bdecat") {
    # generate data.frame, removing unobserved levels in data, as required by BiDAG
    df <- data.frame(apply(data, 2, factor, exclude = NULL, simplify = FALSE))
    scorepar <- BiDAG::scoreparameters("bdecat",
                                       data = df,
                                       bdecatpar = list(chi = par$ess, edgepf = par$edgepf))

  } else if (scoretype %in% c("tree", "ldag", "part")) {

    scorepar <- BiDAG::scoreparameters("usr",
                                       data = data,
                                       usrpar = list(pctesttype = "bdecat"))

    # add cardinality of variables
    scorepar$nlev <- par$nlev
    scorepar$Cvec <- par$nlev # CVec is used in PC-routine of BiDAG-functions (?)

    # additional parameters for user-specified score params
    scorepar$ess <- par$ess
    scorepar$edgepf <- par$edgepf
    scorepar$levels <- lapply(scorepar$nlev-1, seq.int, from = 0)
    scorepar$local_struct <- scoretype
    scorepar$regular <- par$regular


    if (!is.null(lookup)) {
      # define environment in which to s  tore optimized local structures
      assign(scoretype, list(), lookup)  # add list with name `scoretype`
      scorepar$lookup <- lookup
    }

    # define function for computing scores
    usrDAGcorescore <- function(j, parentnodes, n, scorepar) {

      npar <- length(parentnodes)
      pG   <- -npar*log(scorepar$edgepf) # penality term for edges
      bdeu <- bida_bdeu(scorepar$data, j, parentnodes, scorepar$ess, scorepar$nlev)

      if (npar < 2) {
        # no partitioning possible, return score
        return(score_bdeu(bdeu)+pG)
      } else if (is.null(scorepar$lookup)) {
        # optimize partitioning of parent space
        opt <- optimize_bdeu(bdeu, method = scorepar$local_struct, levels = scorepar$levels[parentnodes], regular = scorepar$regular)
        return(sum(opt$scores)+pG)
      } else {
        parId <- paste(c(j, sort(parentnodes)), collapse = ".")
        scores <- scorepar$lookup[[scorepar$local_struct]]

        # if score is already computed, return score
        if (length(scores) > 0 && parId %in% names(scores)) return(scores[[parId]]$score)

        # optimize partitioning of parent space
        opt <- optimize_bdeu(bdeu, method = scorepar$local_struct, levels = scorepar$levels[parentnodes], regular = scorepar$regular)
        score <- sum(opt$scores) +pG

        # store score and bdeu-params in lookup
        bdeu$partition  <- opt$partition
        scorepar$lookup[[scorepar$local_struct]][[parId]] <- list(score = score, bdeu = bdeu)

        return(score)
      }
    } else {
      stop(paste0("Invalid scoretype:", scoretype))
    }


    # assign function to name-space of BiDAG package
    assignInNamespace("usrDAGcorescore", usrDAGcorescore, ns = "BiDAG")
  } else {
    stop(sprintf("Invalid argument for `scoretype`: %s", scoretype))
  }
  return(scorepar)
}
