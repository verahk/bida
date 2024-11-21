

#' Methods for objects of class `posterior_bida`
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
#'
print.posterior_bida <- function(obj) {
  cat("Object of class", class(obj)[[2]], "\n")
  cat("Number of conditional distributions", length(obj$params), "\n")
  cat("Probability mass at zero effect", obj$zeroprob)
}
