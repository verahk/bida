
#' Combine vectors of various length
#'
#' This function combines a list of vectors into a matrix (row-wise),
#' padding vectors of less-than-maximal length with NA values.
#' Used to store adjustment sets in a matrix in [bida::adjsets_support_from_dags].
#'
#' @param x (list)
#' a list with numeric vectors
#' @return a `n-by-m` with `n = min(1, length(x))` and `m = min(1, max(lengths(x)))`.
#' @keywords internal
rbind_fill <- function(x){
  if (is.null(x)) {
    return(NULL)
  } else if (is.list(x)) {
    n    <- length(x)
    lens <- vapply(x, length, integer(1))
    if (n == 1){
      if (lens == 0){
        # return a 1 x 1 matrix with NA if input is a single empty vector
        return(matrix(NA, 1, 1))
      } else {
        return(matrix(x[[1]], 1, lens))
      }
    } else {
      max_len  <- max(lens, 1)
      mat <- matrix(NA, ncol = max_len, nrow = length(lens))

      # enumerate pos in matrix with non-NAs: (i, j) = i + ncol*(j-1)
      indx <- rep.int(seq_along(x), lens) + n*(unlist(lapply(lens, seq_len))-1)
      mat[indx] <- unlist(x)
      return(mat)
    }
  } else {
    if (length(x) == 0){
      return(rbind(NA))
    } else{
      return(rbind(x))
    }
  }
}
