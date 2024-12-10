#' Find ancestors relations corresponding to a network
#'
#' @param x a n-by-n adjacency matrix or a bn.fit object
#' @return a n-by-n adjacency matrix for the acnestor relations
#'
descendants <- function(x) {
  UseMethod("descendants")
}

#' @rdname descendants
#' @export
descendants.matrix <- function(x) {
  x <- as.matrix(x)
  n <- ncol(x)
  I <- diag(n)
  dmat <- try(sign(round(solve(I-x), 1)), silent = TRUE)

  if (inherits(dmat, "try-error")) {
    old <- array(0, dim(x))
    new <- x
    while (!all(old == new)) {
      old <- new
      new <- sign(old%*%(I+x))
    }
    return(sign(I + new))
  } else {
    return(dmat)
  }
}


#' @rdname descendants
#' @export
descendants.default <- function(x) {
  descendants.matrix(as.matrix(x))
}

#' @rdname descendants
#' @export
descendants.bn.fit <- function(x) {
  x <- bnlearn::amat(x)
  descendants.matrix(x)
}

# #profiling:
# set.seed(007)
# res <- list()
# size <- seq(10, 100, by = 10)
# for (n in size) {
#   x <- bida:::randDAG(n, 3)
#   tmp <- microbenchmark::microbenchmark(descendants(x),
#                                         descendants_while(x),
#                                         check = "equal")
#   tmp <- summary(tmp)
#   res[[paste0(n)]] <- setNames(tmp$mean, tmp$expr)
# }
#
#
# matplot(size, do.call(rbind, res), df, type = "l", log = "y")


