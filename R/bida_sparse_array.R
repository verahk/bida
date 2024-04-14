


#' class: cida sparse array
#'
#' Sparse representation of multi-dimensional arrays.
#'
#' @name bida_sparse_array
#' @param x (numeric vector) values of non-zero elements
#' @param index (integer vector) position of non-zero elements
#' @param dims (integer vector) dimension of array
#'
#' @return an object of class `bida_sparse_array`
#' @keywords internal
#'
#' @examples
#' x <- array(c(0, 0, 1, 0, 0, 1, 0, 0, 1), dim = c(3, 3))
#' y <- new_bida_sparse_array(c(1, 1, 1), c(3, 6, 9), c(3, 3))
#' y
#' all.equal(x, as.array(x))
new_bida_sparse_array <- function(x, index, dims, names = NULL) {
  structure(x,
            index = index,
            dims = dims,
            names = names,
            class = "bida_sparse_array")
}


#' @export
as.array.bida_sparse_array <- function(x) {
  y <- array(0, attr(x, "dims"), attr(x, "names"))
  y[index+1] <- x
  return(y)
}

#' @export
dim.bida_sparse_array <- function(x) attr(x, "dims")
