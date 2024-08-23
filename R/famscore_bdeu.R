
#' Compute family score for categorical data (BDeu-score)
#' @param m (matrix) a `q-by-r` matrix with counts (a frequency table)
#' @param x (array) an `1-by-r` matrix or vector of length `r` with counts.
#' @param ess equivalent sample size parameter.
#' @param r cardinality of outcome variable
#' @param q cardinality of parent set
#' @param s number of parent outcomes in each row of `x`. Default to 1.
#' @details
#' - `famscore_bdeu` assumes `x` is a matrix (due to `rowSums`, the functions fails if `dim(x) = NULL`)
#' - `famscore_bdeu_1row` assumes `x` is an 1-by-`r`or array. Faster than `famscore_bdeu` in this case. Used in [optimize_partition].
#' - `famscore_bdeu_byrow` assumes `x` is a matrix and returns a vector with the scores of each row of `x`. Used in [optimize_partition].
#'
#' @returns the bdeu score
#' @keywords internal
#' @example
#'
#' m <- matrix(1:4, 2, 2)
#' ess <- 1
#' famscore_bdeu(m, ess)
#' famscore_bdeu(m, ess) == sum(famscore_bdeu_byrow(m, ess))
#' famscore_bdeu_byrow(m, ess) == apply(m, 1, function(x) famscore_bdeu_1row(x, ess, q = nrow(m)))
#'
#' m <- matrix(0, 2, 2)
#' famscore_bdeu_byrow(m, 1)  # not exact zero
famscore_bdeu <- function(m, ess = 1, r = ncol(m), q = nrow(m)){
  a <- ess/(r*q)
  score <- q*lgamma(r*a) - r*q*lgamma(a) + sum(lgamma(a+m)) - sum(lgamma(r*a+rowSums(m)))
  round(score, 15)
}

#' @rdname famscore_bdeu
famscore_bdeu_byrow <- function(m, ess, r = ncol(m), q = nrow(m), s = 1) {
  ralpha <- ess*s/q
  alpha <- ralpha/r
  scores <- lgamma(ralpha) - r*lgamma(alpha) + rowSums(lgamma(alpha + m)) - lgamma(ralpha + rowSums(m))
  round(scores, 15)
}

#' @rdname famscore_bdeu
famscore_bdeu_1row <- function(x, ess, r = length(x), q = 1, s = 1) {
  ralpha <- ess*s/q
  alpha <- ralpha/r
  score <- lgamma(ralpha) - r*lgamma(alpha) + sum(lgamma(alpha + x)) - lgamma(ralpha + sum(x))
  round(score, 15)
}


