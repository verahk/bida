

#' Check if a partition is CSI consistent
#'
#' @param P (list of integer vectors) 
#'  a partition of the outcome space defined by `levels`.
#'  enumerating the joint outcomes from `0` to `prod(lengths(nlev))-1`. 
#' @param levels (list of integer vectors)
#'  a list with outcomes of each variable.
#' @return a logical constant 
#' @export
#'
#' @examples
#' levels <- list(0:1, 0:1)
#' 
#' # CSI-consistent
#' P <- list(c(0, 1, 3), 2)
#' is_CSI_consistent(P, levels)
#' 
#' # not CSI-consistent partition
#' P <- list(c(0, 3), 2, 1)
#' is_CSI_consistent(P, levels)
#' 
#' P <- list(0, 1:2, 3)
#' is_CSI_consistent(P, levels)
#' 
#' 
#' # CSI-consistent, but not regular
#' P <- list(c(0, 1), c(2, 3))
#' is_CSI_consistent(P, levels)
#' is_regular(P, lengths(levels))
is_CSI_consistent <- function(P, levels, 
                              nlev = lengths(levels), stride = c(1, cumprod(nlev[-length(nlev)]))) {
  for (j in seq_along(P)[lengths(P) > 1]) {
    if (! is_CSI_consistent_part(P[[j]], levels, nlev, stride)) {
      return(FALSE)
    }
  }
  return(TRUE)
}

#' @rdname is_CSI_consistent
#' @param p (integer vector) 
#'  a part of `P`, i.e. a subset of joint outcomes.
#' @param stride (integer vector)
#'  the stride of each variable, corresponding to levels.
is_CSI_consistent_part <- function(p, levels, nlev, stride) {
  seqn <- seq_along(levels)
  K <- p[1]
  j <- 1
  while (j <= length(K) && length(K) < length(p)) {
    pprime  <- K[j]
    pprime0 <- pprime - stride*(pprime%/%stride)%%nlev
    for (i in seqn) {
      tmp <- pprime0[i] + stride[i]*levels[[i]] #stride[i]*(levels[[i]]-pprime%/%stride[i]%%nlev[i])
      if (all(match(tmp, p, 0L) > 0)) {
        K <- c(K, tmp[match(tmp, K, 0L) == 0])
      }
    }
    j <- j+1
  }
  all(match(p, K, 0L) > 0)
}
