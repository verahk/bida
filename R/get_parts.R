

#' Unlist a partition
#'
#' @param P a list that defines a partition of an vector on the form `0:(q-1)`,
#'  e.q. the vector that enumerates the `q` joint outcomes of
#' @param lens (integer vector) number of outcomes in each part of `P`
#' @return an vector of length `unlist(P)` that assigns each outcome two a part in the partition.
#'
#' @examples
#'
#' # list all configurations of two categorical variables
#' levels <- list(0:1, 0:2)
#' conf <- expand.grid(levels)
#'
#' # enumerate joint outcomes
#' nlev <- lengths(levels)
#' stride <- c(1, cumprod(nlev[-length(nlev)]))
#' joint <- as.matrix(conf)%*%stride
#'
#' # define partition of the outcome space
#' P <- list(c(0, 1, 3, 5), 2, 4)
#'
#' cbind(conf, joint, get_parts(P))
get_parts <- function(P, lens = lengths(P)) {
  rep.int(seq_along(P), lens)[order(unlist(P))]
}
