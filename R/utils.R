




match_vec <- function(vec, lookup, lens = lengths(lookup), nomatch = NA) {
  for (i in seq_along(lookup)[lens == length(vec)]) {
    if (identical(vec, lookup[[i]])) {
      return(i)
    }
  }
  return(nomatch)
}



#' Faster version of expand.grid for integer values.
#'
#' @param x a list with integer vectors.
#' @param k the length of each vecotr in `x`.
#' @return a `k-by-k` matrix with all combinations of the vectors in [x].
#'   In contrast \link{base}[expand.grid] returns a data.frame.
#'
#' @examples
#' x = list(0:1, 0:1, 0:1)
#' all(expand.grid(x) == bida:::expand_grid_fast(x))
#' microbenchmark::microbenchmark(as.matrix(expand.grid(x)), bida:::expand_grid_fast(x))
#'
expand_grid_fast <- function(x = NULL, k = NULL){

  if (is.null(x)) {
    if (is.null(k)) stop()
    x <- lapply(k-1, seq.int, from = 0)
  } else if (is.null(k)) {
    k <- vapply(x, length, integer(1))
  }

  n <- length(x)

  if (n == 1){
    return(matrix(x[[1]], ncol = 1))
  } else {
    cump  <- cumprod(k)
    each  <- c(1, cump[-n])
    indx  <- seq.int(n, 1)
    times <- c(cumprod(k[indx[-n]])[indx[-1]], 1)

    vapply(seq_along(x),
           function(i) rep(x[[i]], each = each[i], times = times[i]),
           vector(class(x[[1]]), cump[n]))
  }
}













