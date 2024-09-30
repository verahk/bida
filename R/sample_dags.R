
#' Sample DAGs with local structure
#'
#' This function is a wrapper around [BiDAG::sampleBN], and can be used to sample
#' DAGs with local structure using the different MCMC schemes implemented in [BiDAG::sampleBN].
#' The wrapper function allows for different structure learning algorithms to initiate
#' a search space for the MCMC chain, and forces the maximal parent set size
#' (`hardlimit`) to be  respected by re-running the chosen structure learning procedure
#' with stricter add-edge-policies until the criteria is met. If the criteria is
#' not satisfied either under the strictest add-edge-policy, then edges are removed
#' at random until the criteria is met.
#'
#' @param scorepar (object of class [BiDAG::scoreparameters])
#'  See \link{define_scorepar}.
#' @param algo_init (character)
#'  Name of optimization routine for locastructure.
#' @param algo_sample (character)
#'  Name of MCMC implementation. Either `partition` or `order`.
#'  See [BiDAG::sampleBN].
#' @param hardlimit (integer)
#'  Maximanumber of parents allowed. See [BiDAG::learnBN].
#' @param verbose (logical)
#' @return
#' @export
#' @details
#' The `init_search_space` function is a wrapper around structural learning routines
#' ([bnlearn::hc] and [pcalg::pc]), intended to learn a start-space for the
#' [Bidag::iterativeMCMC] procedure.
#' It forces the maximal number of parents (`hardlimit`) to be respected,
#' by re-running the structure learning procedure with stricter add-edge-policies
#' (`maxp` and `alpha`, respectively) until a structure that satifies the `hardlimit`
#' criteria is inferred.
#' Simply running `Bidag::iterativeMCMC` will run an error message if no skeleton
#' that satisfies the `hardlimit` constraint can be found.
#'
#' @examples
#'
#' nlev <- 2:4
#' data <- sapply(nlev, sample, size = 100, replace = T) -1
#' colnames(data) <- paste0("X", seq_along(nlev))
#'
#' scorepar <- define_scoreparameters(data, "bdecat", par = list(chi  = 1, edgepf = 2))
#' smpl <- sample_dags(scorepar, "pcskel", "order", verbose = T)
#' attributes(smpl)
#' smpl
sample_dags <- function(scorepar, algo_init = "pcskel", algo_sample = "order", hardlimit = 5, verbose = F) {

  tic <- Sys.time()

  # init search space
  if (verbose) cat("\nInit search space:")
  startspace <- init_search_space(scorepar, algo_init, hardlimit, verbose = verbose)
  tic <- c(tic, init = Sys.time())

  # expand search space
  if (verbose) cat("\nExpand search space using BiDAG::iterativeMCMC:")
  #iterfit <- BiDAG::iterativeMCMC(scorepar, hardlimit = hardlimit, startspace = startspace, scoreout = T, verbose = verbose)
  iterfit <- BiDAG::learnBN(scorepar, "orderIter", hardlimit = hardlimit, startspace = startspace, scoreout = T, verbose = verbose)
  tic <- c(tic, expand = Sys.time())

  # run MCMC
  if (verbose) cat("\nSample DAGs using BiDAG::sampleBN:\n")
  smpl <- BiDAG::sampleBN(scorepar, algo = algo_sample, hardlimit = hardlimit, scoretable = BiDAG::getSpace(iterfit), verbose = verbose)
  tic  <- c(tic, sample = Sys.time())

  # add time-tracking as an attribute
  attr(smpl, "toc") <- diff(tic)

  return(smpl)

}

#' @rdname sample_dags
#' @param maxp (integer)
#'  Maximanumber of parents in CPDAG learned by [bnlearn::hc].
#'  Ignored if not `algo=="hc"` or `algo == "hcskel"`.
#' @param alpha (numeric)
#'  Significance levefor the conditionaindependencies test in [pcalg::pc].
#'  Ignored if not `algo=="pc"` or `algo == "pcskel"`.
init_search_space <- function(scorepar, algo, hardlimit, maxp = hardlimit, alpha = .05, verbose) {

  if (algo == "hc") {
    df <- data.frame(apply(scorepar$data, 2, factor, exclude = NULL, simplify = FALSE))
    fit  <- bnlearn:::hc(df, maxp = maxp)
    startspace <- bnlearn::amat(fit)
  } else if (algo == "hcskel") {
    df <- data.frame(apply(scorepar$data, 2, factor, exclude = NULL, simplify = FALSE))
    fit  <- bnlearn:::hc(df, maxp = maxp)
    startspace <- bnlearn::amat(fit)
    startspace <- startspace + t(startspace)
  } else if (algo == "pc") {
    suffStat <- list(dm = scorepar$data, nlev = scorepar$Cvec, adaptDF = FALSE)
    fit <- pcalg::pc(suffStat, pcalg::disCItest, alpha = alpha,
                     labels = colnames(scorepar$data), verbose = verbose)
    startspace <- as(fit@graph, "matrix")
  } else if (algo == "pcskel") {
    suffStat <- list(dm = scorepar$data, nlev = scorepar$Cvec, adaptDF = FALSE)
    fit <- pcalg::skeleton(suffStat, pcalg::disCItest, alpha = alpha,
                           labels = colnames(scorepar$data), verbose = verbose)
    startspace <- as(fit@graph, "matrix")
  }

  if (! any(colSums(startspace) > hardlimit)) {
    return(startspace)
  } else if (switch(substr(algo, 1, 2), "hc" = maxp > 1, "pc" = alpha > 10**-10)) {
    init_search_space(scorepar, algo, hardlimit, maxp = maxp-1, alpha = alpha/2, verbose)
  } else {
    cat("\nCould not find startspace satisfying `hardlimit`. Edges are randomly removed from the startspace. \n")
    cat(sprintf("algo: %s. hardlimit: %s. maxp: %s. alpha: %s. \n", algo, hardlimit, maxp, alpha))

    npar <- colSums(startspace)
    nremove <- npar-hardlimit
    for (i in which(nremove > 0)) {
      remove <- sample(which(startspace[, i] == 1), nremove[i])
      startspace[remove, i] <- 0
    }

    return(startspace)
  }
}
