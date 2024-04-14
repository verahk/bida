#' Compute BDeu family score for categorical data
#' @param m a \code{kx}-by-\code{kpa} matrix with counts (a frequency table).
#' @param ess equivalent sample size parameter.
#' @returns family score.
#'
#'
famscore_bdeu <- function(m, ess = 1){

  kx  <- ncol(m)     # cardinality of x
  kpa <- nrow(m)     # cardinality of parent set
  a  <- ess/(kx*kpa) # bdeu prior for each p(x|pa = u)

  kpa*lgamma(kx*a) - sum(lgamma(kx*a+rowSums(m))) + sum(lgamma(a+m)) - kx*kpa*lgamma(a)
}

famscore_bdeu_part <- function(m, s, q = sum(s), ess = 1) {
  r    <- ncol(m)    # cardinality of x
  aji  <- ess/q*s    # bdeu prior for each context p(pa = u)
  ajil <- aji/r      # bdeu prior for each p(x|pa = u)


  sum(lgamma(aji) - lgamma(rowSums(m)+aji) - r*lgamma(ajil)) + sum(lgamma(ajil+m))
}
