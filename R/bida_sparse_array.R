


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
#'
#' vals <- seq(3, 27, by = 3)
#' x <- new_bida_sparse_array(vals, vals-1, c(3, 3, 3))
#' y <- as.array(x)
#' y
#'
#' dims <- 2:4
#' x <- new_bida_sparse_array(rep(1, prod(dims)), seq_len(prod(dims))-1, dims)
#' y <- aperm(x, c(3:1))
#' stopifnot(all.equal(as.array(y), aperm(as.array(x), c(3:1))))
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

#' @rdname bida_sparse_array
#' @export
aperm.bida_sparse_array <- function(x, perm) {
  indx <- x$index
  dims <- x$dim
  n <- length(dims)

  stride <- c(1, cumprod(dims[-n]))
  conf <- mapply(function(s, k) (indx%/%s)%%k, stride, dims)

  new_dims <- dims[perm]
  new_indx <- c(conf[, perm, drop = F]%*%c(1, cumprod(new_dims[-n])))

  x$dim    <- new_dims
  x$index  <- new_indx

  return(x)
}

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



