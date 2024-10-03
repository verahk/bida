



#' Generic function definitions
#'
#' @export
rowSums <- function(x, na.rm = FALSE, dims = 1, ...) {
  UseMethod("rowSums")
}
#' @export
rowSums.default <- function(x, na.rm = FALSE, dims = 1, ...){
  base:::rowSums(x, na.rm = na.rm, dims = dims, ...)
}
#' @export
colSums <- function(x, na.rm = FALSE, dims = 1, ...) {
  UseMethod("colSums")
}
#' @export
colSums.default <- function(x, na.rm = FALSE, dims = 1){
  base:::colSums(x, na.rm = na.rm, dims = dims)
}

#' @rdname bida_sparse_array
#' @export
as.bida_sparse_array <- function(x, default = 0,...){
  UseMethod("as.bida_sparse_array")
}


#' @export
score <- function(x){
  UseMethod("score")
}

#' @export
asplit <- function(x, MARGIN, ...) {
  UseMethod("asplit")
}
#' @export
asplit.default <- function(x, MARGIN, ...) base:::asplit(x, MARGIN)
