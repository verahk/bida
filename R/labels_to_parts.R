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
#' levels <- rep(list(0:2), 2)
#' nlev   <- lengths(levels)
#' labels <- vector("list", 2)
#' labels[[1]] <- rbind(0, 1)
#' cbind(expand.grid(levels), labels_to_parts(labels, levels))
#'
#' labels[[2]] <- rbind(0)  # add label to edge from parent 2
#' cbind(expand.grid(levels), labels_to_parts(labels, levels))

labels_to_parts <- function(labels, nlev, levels = lapply(nlev-1, seq.int, from = 0)) {
  stride <- c(1, cumprod(nlev[-length(levels)]))
  parts <- seq_len(prod(nlev))
  for (i in seq_along(labels)) {
    parts <- add_label_to_parts(parts, labels[[i]], i, levels, stride)
  }
  parts
}

add_label_to_parts <- function(parts, label, i, levels, stride) {
  if (is.null(label)) return(parts)
  stopifnot(ncol(label) == length(levels)-1)
  # map labels on edge from node i to rows in CPT
  rows <- lapply(label%*%stride[-i], "+", levels[[i]]*stride[i]+1)

  # for each label, collapse the associated parts
  for (x in rows) {
    collapse <- parts[x]                       # parts collapsed by r'th label
    isMember <- match(parts, collapse, 0L) > 0 # members of collapsed part
    parts[isMember] <- collapse[1]             # update partition
  }
  parts
}



