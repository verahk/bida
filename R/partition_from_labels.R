#' Convert label-set to partition
#'
#' Convert the set of labels associated with a node $X_j$ in a labeled DAG to a
#' partition of the outcome space of $Pa(X_j)$.
#'
#' @param labels a list with labels on each edge from each parent $X_i$ into $X_j$.
#' Empty labels have to be included as `NULL` elements in the list, as the variables are referred to by
#' their position not names.
#' @param nlev (integer vector) the cardinality of the parent variables.
#' @param levels (optional) a list with integer vectors, zero-indexed levels of each variable.
#' Included as argument for saving computation time if pre-computed. Otherwise, it can be derived from `nlev`.
#' @return a integer vector of length `prod(nlev)` assigning the parent outcomes to a partition.
#'
#' @examples
#' nlev   <- rep(2, 3)
#' levels <- lapply(nlev-1, seq.int, from = 0)
#' labels <- vector("list", 2)
#' labels[[1]] <- rbind(c(0, 0), c(0, 1))
#' cbind(expand.grid(levels), labels_to_parts(labels, nlev, levels))
#'
#' labels[[2]] <- rbind(c(0, 0))  # add label to edge from parent 2
#' cbind(expand.grid(levels), labels_to_parts(labels, nlev, levels))

partition_from_labels <- function(labels, nlev, levels = lapply(nlev-1, seq.int, from = 0)) {

  n <- length(nlev)
  stride <- c(1, cumprod(nlev[-length(levels)]))
  parts <- seq_len(prod(nlev))

  # check that all labels are correctly formatted
  indx <- !vapply(labels, is.null, logical(1))
  stopifnot(all(vapply(labels[indx], is.matrix, logical(1))))
  stopifnot(all(vapply(labels[indx], ncol, integer(1)) == n-1))

  # update partition with labels on edge from i
  for (i in which(indx)) {
    stopifnot(ncol(labels[[i]]) == n-1)

    # map matrix with labels on edge (i, j) to rows in CPT
    rows  <- lapply(labels[[i]]%*%stride[-i], "+", levels[[i]]*stride[i]+1)

    # for each label, collapse the associated parts
    for (x in rows) {
      collapse <- parts[x]                       # parts to be collapsed by label
      isMember <- match(parts, collapse, 0L) > 0 # members of parts to be collapsed
      parts[isMember] <- collapse[1]             # collapse parts
    }
  }
  parts
}


