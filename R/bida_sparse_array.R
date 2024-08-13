


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
#'
#'
#' # compare size of sparse and non-sparse arrays
#' obj_size <- matrix(NA, 10, 2)
#' sparsity <- .25 # share of non-zero elements
#' for (i in seq_len(nrow(obj_size))) {
#'   dims <- rep(2, i)
#'   index  <- sample(prod(dims), round(prod(dims)*sparsity))
#'   counts <- sample(10, size = length(index), replace = T)
#'
#'   sparse <- new_bida_sparse_array(counts, index, dims)
#'   arr    <- as.array(sparse)
#'   obj_size[i, 1] <- pryr::object_size(sparse)
#'   obj_size[i, 2] <- pryr::object_size(arr)
#' }
#'
#' matplot(obj_size, type = "l", log = "y")

new_bida_sparse_array <- function(value, index, dim, dimnames = NULL) {
  structure(list(value = value,
                 index = index,
                 dim = dim,
                 dimnames = dimnames),
            class = "bida_sparse_array")
}

#' @rdname bida_sparse_array
#' @export
as.array.bida_sparse_array <- function(x) {
  y <- array(0, x$dim, x$dimnames)
  y[x$index+1] <- x$value
  return(y)
}

as.numeric.bida_sparse_array <- function(x) {
  y <- numeric(prod(x$dim))
  y[x$index+1] <- x$value
}

is.bida_sparse_array <- function(x) class(x) == "bida_sparse_array"

colSums_bida_sparse_array <- function(x, dims = 1) {
  tmp <- seq_len(dims)
  colIndex  <- x$index%/%prod(x$dim[tmp])
  newIndex <- unique(group)
  new_bida_sparse_array(rowsum_fast(x$value, colIndex, newIndex),
                        newIndex,
                        x$dim[-tmp])
}

#' @export
levels.bida_sparse_array <- function(x) {
  if (is.null(x$dimnames)) {
    setNames(lapply(x$dim-1, seq.int, from = 0), paste0("X", seq_along(x$dim)))
  } else {
    x$dimnames
  }
}



