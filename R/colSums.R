
#' @export
colSums <- function(x, na.rm = FALSE, dims = 1, ...) {
  UseMethod("colSums")
}
#' @export
colSums.default <- function(x, na.rm = FALSE, dims = 1){
  base:::colSums(x, na.rm = na.rm, dims = dims)
}

