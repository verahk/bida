#' Compute the JSD of a conditional probability table
#'
#' Compute the Jenson-Shannon-Divergence (JSD) of the categorical distributions
#' described by a conditional probability table (CPT).
#'
#' @param p (numeric array) a CPT stored in a `r-by-q` matrix or a collection of
#' `n` such CPTs stored in a `r-by-q-by-n` array, where `r` is the cardinality of
#'  the outcome and `q` the cardinality of the parent variable(s).
#' @param digits (numeric) precision of output.
#'
#' @return a numeric array with the JSD of each CPT in `p`
#' @export
#'
#' @examples
#'
#' a <- .25
#' p1 <- matrix(c(a, 1-a, 1-a, a), nrow = 2, ncol = 2)
#' all(colSums(p1) == 1)  # sum over outcomes equals 1
#' jsd(p1)
#'
#' p2 <- matrix(c(a, 1-a, a, 1-a), nrow = 2, ncol = 2)
#' jsd(p2) # the JSD beteen identical distributions is 0
#'
#' p3 <- diag(2)
#' jsd(p3) == log(2) # log(q) is the upper bound
#' jsd(p3)
#'
#' p3 <- array(c(p1, p2, p3), c(2, 2, 3))
#' jsd(p3) # compute the JSD of a collection of CPTs stored in a array
#'
jsd <- function(p, digits = 10) {
  dims <- dim(p)
  q <- dims[2]
  if (length(dims) == 2) {
    M  <- colSums(t(p))/q                 # average distribution
    HM <- -sum(M*log(M), na.rm = TRUE)    # entropy of average distribution
    Hx <- -sum(p*log(p), na.rm = TRUE)/q  # average entropy of each distribution
  } else {
    M  <- colSums(aperm(p, c(2, 1, 3)))/q
    HM <- -colSums(M*log(M), na.rm = TRUE)
    Hx <- -colSums(p*log(p), na.rm = TRUE, dims = 2)/q
  }
  round(HM-Hx, digits)
}
