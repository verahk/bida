
#' Optimize a local tree structure for a frequency table
#'
#' Optimize a regression-tree partitioning of the conditional outcome space.
#'
#' @rdname optimize_partition
optimize_partition_tree <- function(counts, levels, ess = 1, min_score_improv = 0, verbose = verbose) {


  nlev <- lengths(levels)
  stride <- c(1, cumprod(nlev[-length(nlev)]))
  conf  <- as.matrix(expand.grid(levels))
  r <- ncol(counts)
  q <- nrow(counts)

  score <- famscore_bdeu_1row(colSums(counts), ess = ess, r = r, q = q, s = q)
  tree  <- grow_tree(counts, conf, score, stride, ess, r, q, min_score_improv, verbose = verbose)

  # return partition and scores of each part / leaf
  out <- list(partition = list_leaves(tree, "outcomes"),
              scores = unlist(list_leaves(tree, "scores")))


}


list_leaves <- function(tree, name) {
  if (is.null(tree$branches)) {
    unname(tree[name])
  } else {
    unlist(lapply(tree$branches, list_leaves, name = name), recursive = F)
  }
}


grow_tree <- function(counts, conf, score, stride, ess, r, q, min_score_improv, verbose){


  keep_splitting <- FALSE
  if (nrow(counts) > 1) {
    # find best split
    best_split <- find_best_split(counts, conf, ess, r, q, best_score = score+min_score_improv)
    keep_splitting <- length(best_split) > 0  # if FALSE, no split returns score greater than best_score
  }
  if (keep_splitting) {

    if (verbose) cat("\ndiff:", sum(best_split$scores)-score,
                      "Splitvariable:", best_split$var,
                      "Splitvalue:", best_split$vals,
                      "Scores:", best_split$scores)

    # grow a new tree for each value of split variable
    best_split$branches <- vector("list", length(best_split$vals))
    nvals <- length(best_split$vals)
    for (k in seq_len(nvals)) {
      indx <- conf[, best_split$var] == k-1
      best_split$branches[[k]] <- grow_tree(counts = counts[indx, , drop = F],
                                            conf = conf[indx, , drop = F],
                                            score = best_split$scores[k],
                                            stride, ess, r, q,
                                            min_score_improv,
                                            verbose)
    }
  } else {

    # return leaf node
    best_split <- list(outcomes = conf%*%stride,   # enumerate outcomes in part
                       scores = score)             # score contribution of part
  }
  return(best_split)
}

find_best_split <- function(counts, conf, ess, r, q, best_score = -Inf) {
  best_split <- list()
  for (i in 1:ncol(conf)) {
    vals <- unique(conf[, i])
    if (length(vals) > 1) {
      tmp <- unname(rowsum(counts, conf[, i]))
      scores <- famscore_bdeu_byrow(tmp, ess, r, q, s = tabulate(conf[, i]+1))
      score  <- sum(scores)
      if (score-best_score > 10e-16) {
        best_score <- score
        best_split$var <- i
        best_split$vals <- vals
        best_split$scores <- scores
      }
    }
  }
  return(best_split)
}
