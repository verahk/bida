

#' Compute the mean of an intervention distribution
#'
#' Compute the mean of an intervention distribution by applying the backdoor formula.
#'
#' @param x object of class [`bida_bdeu`]
#' @param ... additional arguments
#' @return the posterior mean(s) of the distribution(s) specified by `x`
#' @export
#'

backdoor_sample <- function(x, ...){
  UseMethod("backdoor_sample", x)
}



