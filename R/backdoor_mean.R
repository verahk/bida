

#' Compute the mean of an intervention distribution
#'
#' Compute the mean of an intervention distribution by applying the backdoor formula.
#'
#' @param x object of class [`bida_bdeu`]
#' @return the mean intervention distribution defined by `x`.
#'  * If `class(x) = "bida_bdeu"` this is a conditional probability table.
#' @export
backdoor_mean <- function(x, nlevx){
  UseMethod("backdoor_mean", x)
}


