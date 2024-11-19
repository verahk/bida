


#' @rdname adjsets_support_from_dags
#' @inheritParams adjsets_support_from_dags
#' @param maxsize (integer) maximum size of parent sets.
#' The support of parent sets `z` such that `length(z) > maxsize` is sat to zero.
#' If no parent set is smaller than `maxsize`, the empty set is returned as the unique parent set..
#' @return
#' - `parent_support_from_dags`returns a `n` length list with elements:
#'    - `sets`: a matrix with the unique parents sets of each variable.
#'    - `support`: the support over parent sets.
#' @export
parent_support_from_dags <- function(dags, w = NULL, x = NULL, checksize = NULL) {

  n <- ncol(dags[[1]])

  if (!is.matrix(dags[[1]])) {
    dags <- lapply(dags, as.matrix)
  }
  if (is.null(w)) {
    w <- rep(1/length(dags), length(dags))
  }
  if (is.null(x)) {
    x <- seqn
  }

  # aggregate weights over unique dags
  u    <- unique(dags)
  w    <- rowsum_fast(w, dags, u)
  dags <- u

  lapply(x, parent_support_from_dags_x, dags = dags, w = w, checksize = checksize)
}

parent_support_from_dags_x <- function(dags, w, x, checksize) {
  parents <- lapply(dags, function(dag) which(dag[, x] == 1))
  u <- unique(parents)
  w <- c(rowsum_fast(w, parents, u))

  if (!is.null(checksize)) {
    tooLarge <- vapply(u, function(z) !checksize(x, 0L, z), logical(1))
    if (all(tooLarge)) {
      # let empty set be only adjset
      u <- list(integer(0))
      w <- 1
    } else if (any(tooLarge)) {
      # remove too large sets
      u[tooLarge] <- NULL
      w <- w[!tooLarge]/(1-sum(w[tooLarge]))
    }
  }
  list(sets = u, w = w)
}
