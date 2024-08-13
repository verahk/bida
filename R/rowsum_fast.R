#' Compute sums over grouping variable
#'
#' Speedy version of [base::rowsum] that assumes takes the unique elements of
#' grouping variable as an argument, and calls the internal method without any
#' error check.
#'
#' @param x (numeric) vector or matrix
#' @param group grouping variable
#' @param ugroup unique elements of grouping variable. Also defines the order of
#' the group sums returned by the function.
#' @param na.rm See [base::rowsum].
#'
#' @return An array or a matrix, depending on the class of x.
#' @export
#'
#' @examples
#'
#' # grouped sum of vector
#' x <- rep(1:3, 3)
#' ugroup <- letters[3:1]
#' group <- rep(ugroup, 3)
#' rowsum(x, group, reorder = F)
#' rowsum_fast(x, group, ugroup)
#'
#' microbenchmark::microbenchmark(rowsum(x, group, reorder = F), rowsum_fast(x, group, ugroup))
#'
#' # grouped sum of matrix
#' x <- matrix(rep(1:9, 3), nrow = 9, ncol = 3)
#' rowsum(x, group, reorder = F)
#' rowsum_fast(x, group, ugroup)
#'
#' microbenchmark::microbenchmark(rowsum(x, group, reorder = F), rowsum_fast(x, group, ugroup))
#'
#' #' # no error on misspecification
#' rowsum_fast(x, group[-1], ugroup, na.rm = F)
#'
#'
#' # with NAs
#' x <- c(rep(1:3, 3), NA)
#' ugroup <- letters[3:1]
#' group <- c(rep(ugroup, 3), letters[3])
#' rowsum_fast(x, group, ugroup, na.rm = F)
#' rowsum_fast(x, group, ugroup, na.rm = T)
#'
rowsum_fast <- function(x, group, ugroup, na.rm = FALSE) {
  UseMethod("rowsum_fast", x)
}

#' @rdname rowsum_fast
#' @export
rowsum_fast.default <- function(x, group, ugroup, na.rm = FALSE) {
  c(.Internal(rowsum_matrix(x, group, ugroup, na.rm, vector("character", length(ugroup)))))
}
#' @rdname rowsum_fast
#' @export
rowsum_fast.matrix <- function(x, group, ugroup, na.rm = FALSE) {
  out <- .Internal(rowsum_matrix(x, group, ugroup, na.rm, vector("character", length(ugroup))))
  dimnames(out) <- NULL
  return(out)
}

