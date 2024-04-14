


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
parent_support_from_dags <- function(dags, support = rep(1/length(dags), length(dags)), checksize = NULL) {

  if (!is.matrix(dags[[1]])) {
    dags <- lapply(dags, as.matrix)
  }

  # list unique dags
  unique_dags <- unique(dags)
  support <- rowsum_fast(support, dags, unique_dags)
  dags <- unique_dags

  n <- ncol(dags[[1]])
  seqn <- seq(n)

  # init
  out <- list(sets = vector("list", n),
              support = vector("list", n))

  # find unique sets
  for (x in seqn) {
    parents <- lapply(dags, function(dag) seqn[dag[, x] == 1])
    u <- unique(parents)
    w <- c(rowsum_fast(support, parents, u))

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

    out$sets[[x]]    <- rbind_fill(u)
    out$support[[x]] <- w
  }

  return(out)
}
