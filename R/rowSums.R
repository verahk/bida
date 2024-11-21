#' @export
rowSums <- function(x, na.rm = FALSE, dims = 1, ...) {
  UseMethod("rowSums")
}
#' @export
rowSums.default <- function(x, na.rm = FALSE, dims = 1, ...){
  base:::rowSums(x, na.rm = na.rm, dims = dims, ...)
}
