

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_dim <- function(x) {
  UseMethod("get_dim", x)
}
#' @rdname
#' @export
get_dim.default <- function(x) dim(x)
#' @rdname
#' @export
get_dim.bida_sparse_array <- function(x) x$dim

#' @rdname
#' @export
get_dim.bida_bdeu <- function(x) get_dim(x$counts)
