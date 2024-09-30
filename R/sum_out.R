

#' Compute marginal sums
#'
#' @param x an object of class `base::array` or `bida_sparse_array`
#' @param MARGIN (integer vector)
#'
#' @return an array with marginal sums
#' @examples
#'
#' arr <- array(1:27, c(3, 3, 3), list(x = 0:2, y = 0:2, z = 0:2))
#' sum_out(arr, 2)
#' sum_out(arr, 2:3) == colSums(arr, dims = 1)
#' sum_out(arr, 1:2) == rowSums(arr, dims = 2)
#' sum_out(arr, c(1, 3)) == apply(arr, c(1, 3), sum)
#'
#' microbenchmark::microbenchmark(sum_out(arr, c(1, 3)), apply(arr, c(1, 3), sum))
#'
#' # sparse array
#' sarr <- as.bida_sparse_array(arr)

sum_out <- function(x, MARGIN) {
  dims <- dim(x)

  if (is.array(x)) {
    perm <- c(seq_along(dims)[-MARGIN], MARGIN)
    colSums(aperm(x, perm), dims = length(dims[-MARGIN]))
  } else if (inherits(x, "bida_sparse_array")) {
    index <- get_index(x, MARGIN)
    new_index <- c(unique(index))
    new_dims  <- dims[MARGIN]
    new_dimnames <- dimnames(x)[MARGIN]
    # aggregate values for each level of MARGIN
    if (x$default == 0) {
      new_value <- rowsum_fast(x$value, index, new_index)
      new_default <- 0
    } else {
      k <- prod(dims[-MARGIN])    #
      nmiss  <- k-tabulate(match(index, new_index))
      new_value <- rowsum_fast(x$value, index, new_index) + nmiss*x$default
      new_default <- k*x$default
    }
    new_bida_sparse_array(new_value, new_index, new_dims, new_dimnames, new_default)
  } else {
    stop("x is neither an array or an bida_sparse_array")
  }
 }
