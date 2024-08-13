
#' Optimize local bdeu score
#'
#' @param data
#' @param levels
#' @param nlev
#' @param j (integer)
#'  column index of node
#' @param parentnodes (integer vector)
#'  column indicies of parent nodes. A vector of length 0 (e.g. `integer(0)`)
#'  translates to no parents.
#' @param ess (numeric constant)
#'  equivalent sample size.
#' @param method (character)
#'  algorithm to learn local structure. See \link{optimize_partition}.
#'  If `NULL` (default) no local structure is learned. Ignored if `length(parentnodes) < 2`.
#' @param ... additional parameters sent to [optimize_partition()].
#'
#' @return a numeric constant, the Bdeu-score. If a local structure is inferred,
#'  this is returned as an named attribute (`structure`).
#' @export
#'
#' @examples

optimize_famscore_bdeu <- function(data, levels, nlev = lengths(levels), j, parentnodes, ess = 1, method = NULL, ...) {
  npar <- length(parentnodes)
  r <- nlev[j]
  if (npar == 0) {
    famscore_bdeu_1row(tabulate(data[, j]+1, r), ess, r)
  } else if (npar == 1) {
    q <- nlev[parentnodes]
    joint <- data[, parentnodes] + q*data[, j]
    tab <- matrix(tabulate(joint+1, q*r), q, r)
    famscore_bdeu(tab, ess, r, q)
  } else {

    # enumerate joint outcomes
    stride <- c(1, cumprod(nlev[parentnodes]))
    joint  <- data[, c(parentnodes, j)]%*%stride

    # compute freq table
    r <- nlev[j]
    q <- stride[npar+1]
    tab <- matrix(tabulate(joint+1, q*r), q, r)

    # optimize partition over parent outcomes
    struct <- optimize_partition(tab, levels[parentnodes], ess, method = method, ...)

    # return family score
    score  <- sum(struct$scores)
    attr(score, "structure") <- struct
    return(score)
  }
}
