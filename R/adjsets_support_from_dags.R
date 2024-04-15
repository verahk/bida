#' Compute support over adjustment sets
#'
#' Identify and compute support over unique adjustment sets in a list of DAGs.
#'
#' @inheritParams adjsets_from_dag
#' @param dags (list)
#' a list of dags, e.g. a sample obtained with [BiDAG::partitionMCMC]
#' @param support (numeric vector)
#' support for each dag in `dags`.
#' @return
#' - `adjsets_support_from_dags` returns a list that includes for each adjset in `adjsets`:
#'    - `sets`: a `n-by-n` list-matrix, where `sets[[x, y]]` is a matrix
#'    listing the unique adjsets w.r.t a cause `x` and effect `y` rowwise
#'    - `support`: a `n-by-n` list-matrix where `support[[x, y]]` is support over
#'    the corresponding sets
#'
#'@details
#' - `adjsets_support_from_dags`:
#'    Computes support over unique adjustment set for each cause-effect pair (`x`, `y`),
#'    where the adjustment eqauls `y` in all DAGs where `y` is not a descendant of `x`.
#'    See also [adjsets_from_dag].
#' - `parent_support_from_dags`:
#'    Computes support over unique parent sets of each variable, without encoding
#'    zero-effects for non-descendants of the variables.
#'@seealso [adjsets_from_dag]
#' @export
adjsets_support_from_dags <- function(adjsets, dags, support = rep(1/length(dags), length(dags)), checksize = NULL){

  # sort adjset for more eff computation / reusing large sets
  indx <- match(c("pa", "pa_min", "o", "o_min"), adjsets, nomatch = 0)
  adjsets_ordered <- adjsets[indx]

  n <- dim(dags[[1]])[1]
  seqn <- seq_len(n)
  seqa <- seq_along(adjsets)

  # ensure that list of dags are matrix (not sparseMatrix)
  if (!is.matrix(dags[[1]])) {
    dags <- lapply(dags, as.matrix)
  }

  # list unique dags
  unique_dags <- unique(dags)
  support <- rowsum_fast(support, dags, unique_dags)
  dags <- unique_dags
  rm(unique_dags)


  # init sets as is full support for zero-effects
  tmp <- list(sets = matrix(lapply(seqn, matrix), n, n, byrow = T),
              support = matrix(list(1), n, n))
  diag(tmp$sets) <- diag(tmp$support) <- list(NULL)
  out <- rep(list(tmp), length(adjsets))
  names(out) <- adjsets_ordered

  ## find unique adjusment sets
  ## - for every dag, x, y:
  ##   if adjustment do not match previous sets, add to list of unique sets
  ##   otherwise update support for that set
  ## - for efficiency: compute support for non-Descendants separately

  # list for storing sets and support
  sets <- supp <- lens <- array(list(), c(n, n, length(adjsets)))

  # matrix for ancestor relations
  arp <- matrix(0, n, n)

  # find unique sets by looping through every g, x, y
  for (g in seq_along(dags)){
    w <- support[g]
    dag  <- as.matrix(dags[[g]])
    dmat <- descendants(dag)

    # update support of anc relations
    arp <- arp + w*dmat

    for (x in seqn[rowSums(dmat) > 1]) {      # for every cause

      # identify descendants
      De <- seqn[-x][dmat[x, -x] == 1]

      # find adjsets
      tmp <- adjsets_from_dag(adjsets_ordered, dag, dmat, xvars = x, yvars = De, checksize = checksize, simplify = FALSE)

      for (y in De) { # for every descendant of x
        for (a in seqa){
          z <- tmp[[x, y, a]]
          if (! (length(z) == 1 && is.na(z)) ) {
            # find position of z in current list of unique sets
            pos  <- match_vec(z, sets[[x, y, a]], lens[[x, y, a]], nomatch = 0)
            if (pos > 0) {
              supp[[x, y, a]][[pos]] <- supp[[x, y, a]][[pos]] + w
            } else {
              pos <- length(sets[[x, y, a]]) + 1
              sets[[x, y, a]][[pos]] <- z
              supp[[x, y, a]][pos] <- w
              lens[[x, y, a]][pos] <- length(z)
            }
          }
        }
      }
    }
  }
  rm(dags, support)
  arp <- round(arp, 10)

  ## store sets in matrix
  possEffects <- seqn[colSums(arp) > 0]
  for (y in possEffects) {

    possCauses <- seqn[-y][arp[-y, y] > 0]

    for (x in possCauses) {

      # add support for zero-effect
      w <- 1-arp[x, y]
      if (w > 0) {
        for (a in seqa) {
          pos <- length(sets[[x, y, a]]) + 1
          sets[[x, y, a]][[pos]] <- y
          supp[[x, y, a]][[pos]] <- w
        }
      }


      for (a in seqa) {

        # if all sets are NA (too large) - return empty set with full support
        if (length(sets[[x, y, a]]) == 0) {
          out[[a]]$sets[[x, y]] <- matrix(NA, 1, 1)
          out[[a]]$support[[x, y]] <- 1
          next
        }

        # store as matrix
        w <- supp[[x, y, a]]
        out[[a]]$support[[x, y]] <- w/sum(w)
        out[[a]]$sets[[x, y]] <- rbind_fill(sets[[x, y, a]])

        sets[x, y, a] <- supp[x, y, a] <- list(NULL)
      }
    }
  }

  rm(sets)
  rm(supp)
  out <- out[adjsets] # re-order
  return(out)
}




adjsets_support_from_dags_2 <- function(adjsets, dags, support = rep(1/length(dags), length(dags)), checksize = NULL){

  # sort adjset for more eff computation / reusing large sets
  indx <- match(c("pa", "pa_min", "o", "o_min"), adjsets, nomatch = 0)
  stopifnot(any(indx>0))
  adjsets_ordered <- adjsets[indx]

  n <- dim(dags[[1]])[1]
  seqn <- seq_len(n)

  # ensure that list of dags are matrix (not sparseMatrix)
  if (!is.matrix(dags[[1]])) {
    dags <- lapply(dags, as.matrix)
  }

  # list unique dags
  unique_dags <- unique(dags)
  support <- rowsum_fast(support, dags, unique_dags)
  dags <- unique_dags
  rm(unique_dags)

  dmats <- lapply(dags, descendants)
  arp   <- round(Reduce("+", mapply("*", dmats, support, SIMPLIFY = FALSE)), 10)

  possCauses  <- seqn[rowSums(arp) > 0]

  # init sets as is full support for zero-effects
  tmp <- list(sets = matrix(lapply(seqn, matrix), n, n, byrow = T),
              support = matrix(list(1), n, n))
  diag(tmp$sets) <- diag(tmp$support) <- list(NULL)

  out <- rep(list(tmp), length(adjsets))
  names(out) <- adjsets_ordered

  seqn_dags <- seq_along(dags)

  for (x in possCauses) {
    possEffects <- seqn[-x][arp[x, -x] > 0]
    sets <- lapply(seqn_dags,
                   function(g) adjsets_from_dag(adjsets_ordered, dags[[g]], dmats[[g]],
                                                xvars = x,
                                                yvars = possEffects,
                                                checksize = checksize,
                                                simplify = T))
    sets <- array(unlist(sets, recursive = F),
                  dim = c(length(possEffects), length(adjsets), length(dags)))

    for (yy in seq_along(possEffects)) {
      y <- possEffects[yy]
      for (a in seq_along(adjsets)) {

        u <- unique(sets[yy, a, ])
        w <- rowsum_fast(support, sets[yy, a, ], u)

        # check for NAs
        indx <- vapply(u, function(z) length(z) == 1 && is.na(z), logical(1))
        if (any(indx)) {
          if (length(u) == 1) {
            out[[a]]$sets[[x, y]]    <- matrix(NA, 1, 1)  # encodes an empty adjustment set
            out[[a]]$support[[x, y]] <- 1
          } else {
            u[indx] <- NULL
            w <- w[indx]/(1-sum(w[indx]))
          }
        }

        out[[a]]$sets[[x, y]] <- rbind_fill(u)
        out[[a]]$support[[x, y]] <- w
      }
    }
  }

  return(out)
}

