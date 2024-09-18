

#' Identify adjustment sets in a dag
#'
#' @param adjsets (character vector)
#' which adjset to identify. A subset of `c("anc", "pa", "pa_min", "o", "o_min")`.
#' @param dag (integer matrix)
#' `n-by-n` adjacency matrix, where `dag[i, j] = 1` indicates an edge from node `i` to node  `j`.
#' @param dmat (integer matrix, optional)
#' adjacenecy matrix indicating ancestor relationships corresponding to `dag`.
#' @param xvars,yvars (integer vector, optional)
#' column position(s) of causes and effects for which adjustment sets should be identified.
#' If `NULL` the full sequence `1:n` is used.
#' @param checksize (function, optional)
#' to evaluate if an adjustment set it to large. Defaults to `NULL`. See details.
#' @return
#' - `adjsets_from_dag`: a 3-dimensional list where element (`x`, `y`, `adjset`) is a
#' is the adjustment set of class `adjset` w.r.t. `x`, `y` and `dag`.
#' For every `x` that is not a cause of `y`, as indicated by `dmat[x, y] == 0`,
#' the element (`x`, `y`, `adjset`) is `y`.
#'
#' @details
#' `adjsets_from_dag` identifies all adjustment sets in `adjsets` for every
#'  cause-effect pair given by `xvars` and `yvars`.
#'
#'
#' The adjustment sets are (w.r.t a cause-effect pair `(x, y)`):
#' - `anc`: Joint ancestors of `x` and `y`, excluding descendants of `x`.
#' - `pa`: Parents of `x`.
#' - `pa_min`: Minimal parent set of `x` w.r.t `y`.
#' - `o`: o-set of `x` w.r.t `y`.
#' - `o_min`: Minimal o-set of `x` w.r.t `y`.
#'
#' Treatment of non-descendants:
#' - if `dmat[x, y] = 0`, `y` is returned for the sets `pa_min, o, o_min` (indicating a zero-effect).
#'
#' Large adjustment sets:
#' If `checksize = NULL` no test for size of the set is applied (default).
#' Otherwise, if `checksize(x, y, z) = FALSE` the adjustment set is replaced with smaller,
#' valid set in this order:
#' `o` > `o_min` > `pa_min` > `NA`
#' If the minimal parent set (usually the smallest set) is too large, the
#' function return `NA`.
#' @seealso [find_nearest_adjset]
#' @example man/examples/adjsets_from_dag.R
#' @export
adjsets_from_dag <- function(adjsets, dag, dmat = NULL, xvars = NULL, yvars = NULL, checksize = NULL, simplify = TRUE){

  n <- ncol(dag)
  seqn <- seq_len(n)
  dag  <- unname(as.matrix(dag, list(NULL, NULL)))

  if (is.null(xvars)) xvars <- seqn
  if (is.null(yvars)) yvars <- seqn
  if (is.null(dmat))  dmat  <- descendants(dag)


  # init a 3-dim array for storing adjsets
  sets <- array(list(), c(n, n, length(adjsets)))
  dimnames(sets) <- list(NULL, NULL, adjsets)

  # adjustment sets for non-zero effects
  for (x in xvars) {

    # backdoor graph
    bdag <- dag
    bdag[x, ] <- 0

    # ancestors of x
    ancx <- dmat[, x] == 1


    for (y in yvars[!yvars == x]) {

      # ancestors of x and/or y
      anc   <- ancx | dmat[, y] == 1

      # find adjsets
      tmp <- sets[x, y, ]
      for (adjset in adjsets) {
        tmp[[adjset]] <- adjset_from_bdag(adjset, bdag, dmat, x, y, anc, tmp)
      }

      # replace large sets
      if (!is.null(checksize)) {
        tmp <- replace_large_adjsets(tmp, checksize, bdag, dmat, x, y, anc)
      }

      sets[x, y, ] <- tmp
    }
  }

  if (simplify) {
    return(sets[xvars, yvars, adjsets])
  } else {
    return(sets)
  }
}



#' @rdname adjsets_from_dag
#' @param bdag (integer matrix)
#' adjacency matrix of backdoor graph w.r.t. a single cause `x`
#' @param dmat (integer matrix)
#' adjacenecy matrix indicating ancestor relationships corresponding to `dag`.
#' Used to identify ancestors of the conditioning set in `find_nearest_adjset`.
#' @param x,y (integer)
#' column position of a cause variable `x` and effect variable `y`
#' @param anc (bolean vector)
#' indicator for the column position(s) of ancestors of `x` and `y`. rowSums(dmat[, c(x, y)]) > 0
#' @param sets (list, optional)
#' named list of valid adjustment sets. For avoiding repeated identification of
#' larger sets.
#' @return (integer vector)
#' - `adjset_from_bdag`: column positions of the variables that consititutes the
#' adjustment set of class `adjset` w.r.t `x`, `y` in `bdag`.
#' @details
#' Identify a single adjustment set given the backdoor graph and joint ancestors.
#' Called from the wrapper [bida::adjsets_from_dag].
adjset_from_bdag <- function(adjset, bdag, dmat, x, y, anc, sets = NULL) {
  if (!is.null(sets[[adjset]])) {
    return(sets[[adjset]])
  } else if (dmat[x, y] == 0 && adjset %in% c("o", "o_min", "pa_min")) {
    return(y)
  } else {

    # find adjustment set
    z <- switch(adjset,
              "o" = {
              z0 <- seq_along(anc)[anc][dmat[x, anc] == 0]
              find_nearest_adjset(bdag, dmat, y-1, anc, z0-1) +1
            }, "o_min" = {
              z0 <- adjset_from_bdag("o",  bdag, dmat, x, y, anc, sets)
              find_nearest_adjset(bdag, dmat, x-1, anc, z0-1) +1
            }, "pa_min" = {
              z0 <- seq_along(anc)[bdag[, x] == 1]
              find_nearest_adjset(bdag, dmat, y-1, anc, z0-1) +1
            }, "anc" = seq_along(anc)[anc][dmat[x, anc] == 0],
              seq_along(anc)[bdag[, x] == 1]) # nomatch: return parent set

    return(as.integer(z))
  }
}

#' @rdname adjsets_from_dag
replace_large_adjsets <- function(sets, checksize, bdag, dmat, x, y, anc) {

  matches <- match(c("pa_min", "pa", "o_min", "o"), names(sets), nomatch = 0)

  # pa_min:
  a <- matches[1]
  if (a > 0) {
    z <- sets[[a]]
    if (!(checksize(x, y, z)) ) {
      sets[[a]] <- NA
    }
  }

  # pa:
  a <- matches[2]
  if (a > 0) {
    z <- sets[[a]]
    if (!(checksize(x, y, z)) ) {
      z <- adjset_from_bdag("pa_min", bdag, dmat, x, y, anc, sets)
      if (matches[1] == 0 && !( checksize(x, y, z)) ) { # check smaller set if not already done
        z <- NA
      }
    }
  }

  # o_min:
  a <- matches[3]
  if (a > 0){
    z <- sets[[a]]
    if (!(checksize(x, y, z)) ) {
      z <- adjset_from_bdag("pa_min", bdag, dmat, x, y, anc, sets)
      if (matches[1] == 0 &&!(checksize(x, y, z)) ) { # check smaller set if not already done
        z <- NA
      }

      sets[matches[3:4]] <- list(z) # replace both o and o_min
      return(sets)                  # no need to check o again
    }
  }

  # o:
  a <- matches[4]
  if (a > 0){
    z <- sets[[a]]
    if (!(checksize(x, y, z)) ) {
      z <- adjset_from_bdag("o_min", bdag, dmat, x, y, anc, sets)
      if (matches[3] == 0 &&!(checksize(x, y, z)) ) { # check smaller set if not already done
        z <- adjset_from_bdag("pa_min", bdag, dmat, x, y, anc, sets)
        if (matches[1] == 0 &&!(checksize(x, y, z)) ) { # check smaller set if not already done
          z <- NA
        }
      }
      sets[[a]] <- z
    }
  }
  return(sets)
}
