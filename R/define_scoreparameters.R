#' User defined scoreparameters for `BiDAG`
#'
#' Define a [BiDAG::scoreparameters] object with a user defined score for learning
#' DAGs with local structure.
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
    if (is.null(par$local_struct) || par$local_struct == "none") {
      # use BiDAG-package score function

      # generate data.frame, removing unobserved levels in data, as required by BiDAG
      data <- data.frame(apply(data, 2, factor, exclude = NULL, simplify = FALSE))
      scorepar <- BiDAG::scoreparameters("bdecat",
                                         data = data,
                                         bdecatpar = list(chi = par$ess, edgepf = par$edgepf))

    } else if (par$local_struct %in% c("tree", "ldag", "ptree", "pcart")) {

      # init BiDAG::scoreparameters-object
      scorepar <- BiDAG::scoreparameters("usr", data = data, usrpar = list(pctesttype = scoretype))

      # add parameters for user-defined score and local-structure optimiziation
      scorepar <- c(scorepar, par)

      # add cardinality of variables
      scorepar$levels <- lapply(par$nlev-1, seq.int, from = 0)
      scorepar$Cvec   <- par$nlev # CVec is used in PC-routine of BiDAG-functions (?)

      if (par$local_struct == "pcart") {
        #  add data.frame version of data with factor variables, INCLUDING unobserved levels
        #  - because rpcart:::opt.pcart.cat.bdeu takes a data.frame argument

        varnames <- colnames(data)
        if (is.null(varnames)) varnames <- paste0("X", seq_along(scorepar$levels))
        # generate data.frame
        tmp <- lapply(seq_along(scorepar$levels),
                      function(x) factor(data[, x], scorepar$levels[[x]]))
        scorepar$data.frame <- setNames(data.frame(tmp), varnames)
      }

      if (is.environment(lookup)) {
        # add lookup-environment in which to store optimized local structures
        assign(scorepar$local_struct, list(), lookup)  # add list with name `scoretype`
        scorepar$lookup <- lookup
      }

      # define function for computing scores
      usrDAGcorescore <- function(j, parentnodes, n, scorepar) {
        parentnodes <- sort(parentnodes)   # for parID
        npar <- length(parentnodes)
        pG   <- -npar*log(scorepar$edgepf) # penality term for edges

        if (is.environment(scorepar$lookup)) {
          # look for score in lookup-table
          parID <- paste(c(j, parentnodes), collapse = ".")
          scores <- scorepar$lookup[[scorepar$local_struct]]
          if (length(scores) > 0 && parID %in% names(scores)) {
            return(scores[[parID]]$score+pG)
          }
        }

        # construct bida_bdeu-object - for computing score and optimize partition
        bdeu <- bida_bdeu(scorepar$data, j, parentnodes, scorepar$ess, scorepar$nlev)

        # optimize partition of parent space and return score
        if (npar < 2) {
          # no partitioning possible, return score
          # construct bdeu-object with frequency counts
          return(score_bdeu(bdeu)+pG)
        } else if (scorepar$local_struct == "pcart") {
          opt <- optimize_partition_from_df_pcart(scorepar$data.frame, j, parentnodes, scorepar$ess, scorepar$nlev)
        } else {
          # optimize partitioning of parent space
          opt <- optimize_bdeu(bdeu,
                               method = scorepar$local_struct,
                               levels = scorepar$levels[parentnodes],
                               regular = scorepar$regular)
        }

        stopifnot("scores" %in% names(opt))
        score <- sum(opt$scores)

        if (!is.null(scorepar$lookup)) {
          # store score and bdeu-params in lookup
          bdeu$partition  <- opt$partition
          bdeu$score <- score
          scorepar$lookup[[scorepar$local_struct]][[parID]] <- bdeu
        }

        return(score+pG)
      }

      # assign function to name-space of BiDAG package
      assignInNamespace("usrDAGcorescore", usrDAGcorescore, ns = "BiDAG")

    } else {
      stop(cat("Invalid argument for `par$local_struct`: ", par$local_struct))
    }
  } else {
    stop(cat("Invalid argument for `scoretype`:", scoretype))
  }
  return(scorepar)
}
