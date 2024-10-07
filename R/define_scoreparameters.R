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
#' nlev <- sapply(bn, function(x) dim(x$prob)[1])
#' data <- sample_data_from_bn(bn, 10)
#'
#' par <- list(nlev = nlev, ess = 1, edgepf = 2, regular = T)
#'
#' j <- 8
#' parentnodes <- 5:6
#'
#' # no partitioning - use BiDAG-score
#' scorepar <- define_scoreparameters(data, "bdecat", par, lookup = lookup)
#' score <- BiDAG:::DAGcorescore(j, parentnodes, ncol(data), scorepar)
#' score
#'
#' # user-specified function is assigned to BiDAG-namespace
#' par$local_struct = "tree"
#' scorepar <- define_scoreparameters(data, "bdecat", par, lookup = lookup)
#' BiDAG:::usrDAGcorescore
#'
#' # computes score with a call to BiDAG:::usrDAGcorescore
#' score <- BiDAG:::DAGcorescore(j, parentnodes, ncol(data), scorepar)
#' score
#'
#' # user-specified functions writes to lookup table
#' ls.str(lookup)    # updated
#' lookup$tree[["8.5.6"]]
#'
define_scoreparameters <- function(data, scoretype, par = NULL, lookup = NULL) {

  if (scoretype == "bdecat") {
    if (is.null(par$local_struct) || par$local_struct == "none") {
      # use BiDAG-package score function

      # generate data.frame, removing unobserved levels in data, as required by BiDAG
      data <- data.frame(apply(data, 2, factor, exclude = NULL, simplify = FALSE))
      scorepar <- BiDAG::scoreparameters("bdecat",
                                         data = data,
                                         bdecatpar = list(chi = par$ess, edgepf = par$edgepf))

    } else {

      # init BiDAG::scoreparameters-object
      scorepar <- BiDAG::scoreparameters("usr",
                                         data = data,
                                         usrpar = list(chi = par$ess,
                                                       edgepf = par$edgepf,
                                                       pctesttype = scoretype))

      # additional params
      if (is.null(par$nlev)) {
        scorepar$Cvec <- apply(data, 2, max) +1
      } else {
        scorepar$Cvec   <- par$nlev # CVec is used in PC-routine of BiDAG-functions (?)
      }
      scorepar$levels <- lapply(scorepar$Cvec-1, seq.int, from = 0)
      scorepar$local_struct <- par$local_struct

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

        ess <- scorepar$chi
        nlev <- scorepar$Cvec
        local_struct <- scorepar$local_struct

        npar <- length(parentnodes)
        pG   <- -npar*log(scorepar$pf) # penality term for edges

        # optimize partition of parent space and return score
        if (npar == 0) {
          r <- nlev[j]
          q <- 1
          freqtab <- matrix(tabulate(data[, j]+1, r), q, r)
          famscore <- famscore_bdeu(freqtab, ess, r, q)
        } else if (npar == 1) {
          r <- nlev[j]
          q <- nlev[parentnodes]
          joint <- scorepar$data[, c(parentnodes, j), drop = FALSE]%*%c(1, q)
          freqtab <- matrix(tabulate(joint+1, r*q), q, r)
          famscore <- famscore_bdeu(freqtab, ess, r, q)
        } else {
          if (local_struct == "pcart") {
            opt <- optimize_partition_from_df_pcart(scorepar$data.frame[, c(j, parentnodes)], ess)
          } else {
            opt <- optimize_partition_from_data(scorepar$data, j, parentnodes, ess, nlev, local_struct)
          }

          famscore  <- attr(opt, "score")
          # store partition in lookup-table
          if (!is.null(scorepar$lookup)) {
            parentnodes <- sort(parentnodes)
            parID <- sprintf("%s.%s", j, paste0(parentnodes, collapse = "."))
            scorepar$lookup[[scorepar$local_struct]][[parID]] <- opt
          }
        }
        return(famscore+pG)
      }

      # assign function to name-space of BiDAG package
      assignInNamespace("usrDAGcorescore", usrDAGcorescore, ns = "BiDAG")

    }
  } else {
    stop(cat("Invalid argument for `scoretype`:", scoretype))
  }
  return(scorepar)
}
