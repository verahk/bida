
#' Optimize a local tree structure for a frequency table
#'
#' Optimize a regression-tree partitioning of the conditional outcome space.
#'
#' @rdname optimize_partition
#' @examples
#'
#' levels <- list(0:1, 0:1)
#' counts <- cbind(10, c(10, 100, 1000, 1000))
#' ess <- 1
#' optimize_partition_tree(counts, levels, ess, min_score_improv = 0, verbose = TRUE)
#'
#' # illustrate grow-full-and-then-prune procedure
#' levels <- list(0:1, 0:1, 0:1)
#' counts <- cbind(c(1, 1, 1, 1, 10, 20, 20, 2), c(1, 1, 1, 1, 10, 2, 2, 10))
#'
#' ## stop splitting when additional split do not improve score (default)
#' fit <- optimize_partition_tree(counts, levels, ess, min_score_improv = 0, verbose = TRUE)
#' cbind(counts, get_parts(fit$partition))
#' sum(fit$scores)
#'
#' ##  grow full tree
#' fit <- optimize_partition_tree(counts, levels, ess, min_score_improv = -Inf, verbose = TRUE)
#' cbind(counts, get_parts(fit$partition))
#' sum(fit$scores)
#'
#' #
#' fit <- optimize_partition_tree(counts, levels, ess, min_score_improv = -Inf, TRUE, verbose = TRUE)
#' cbind(counts, get_parts(fit$partition))
#' sum(fit$scores)
#'
optimize_partition_tree <- function(counts, levels, ess = 1, min_score_improv = 0, prune = FALSE, verbose = verbose) {

  nlev <- lengths(levels)
  stride <- c(1, cumprod(nlev[-length(nlev)]))
  conf  <- as.matrix(expand.grid(levels))
  r <- ncol(counts)
  q <- nrow(counts)

  score <- famscore_bdeu_1row(colSums(counts), ess = ess, r = r, q = q, s = q)
  tree  <- grow_tree(counts, conf, score, ess, r, q, min_score_improv, stride, verbose = verbose)

  if (prune) {
    tree <- prune_tree(tree)
  }
  # return partition and scores of each part / leaf
  out <- list(tree = tree,
              partition = list_leaves(tree, "outcomes"),
              scores = unlist(list_leaves(tree, "score")))

  return(out)
}


prune_tree <- function(tree) {
  if (is.null(tree$branches)) {
    return(tree)
  }

  for (b in seq_along(tree$branches)) {
    tree$branches[b] <- list(prune_tree(tree$branches[[b]]))
  }

  indx <- vapply(tree$branches, function(x) is.null(x$branches), logical(1))
  if (all(indx) && (tree$score > sum(unlist(list_leaves(tree, "score"))))) {
      list(score = tree$score,
           outcomes = unlist(list_leaves(tree, "outcomes")))
  } else {
    tree
  }
}

list_leaves <- function(tree, name) {
  if (is.null(tree$branches)) {
    unname(tree[name])
  } else {
    unlist(lapply(tree$branches, list_leaves, name = name), recursive = F)
  }
}


grow_tree <- function(counts, conf, score, ess, r, q, min_score_improv, stride, verbose){

  keep_splitting <- FALSE
  if (nrow(counts) > 1) {
    # find best split
    best_split <- find_best_split(counts, conf, score, ess, r, q, min_score_improv)
    keep_splitting <- length(best_split) > 1  # if FALSE, no split returns score greater than score
  }
  if (keep_splitting) {

    if (verbose) cat("\ndiff:", best_split$score-score,
                      "Splitvariable:", best_split$var,
                      "Splitvalue:", best_split$vals,
                      "Scores:", best_split$leaf_scores)

    # grow a new tree for each value of split variable
    best_split$branches <- setNames(vector("list", length(best_split$vals)), best_split$vals)
    nvals <- length(best_split$vals)
    for (k in seq_len(nvals)) {
      indx <- conf[, best_split$var] == k-1
      best_split$branches[[k]] <- grow_tree(counts = counts[indx, , drop = F],
                                            conf = conf[indx, , drop = F],
                                            score = best_split$leaf_scores[k],
                                            ess, r, q,
                                            min_score_improv,
                                            stride,
                                            verbose)
    }
  } else {

    # return leaf node
    best_split <- list(outcomes = c(conf%*%stride), # enumerate outcomes in part
                       score = score)
  }
  return(best_split)
}

find_best_split <- function(counts, conf, score, ess, r, q, min_score_improv) {
  best_split <- list(score = score)
  best_score <- score+min_score_improv
  for (i in 1:ncol(conf)) {
    vals <- unique(conf[, i])
    if (length(vals) > 1) {
      tmp <- unname(rowsum(counts, conf[, i]))
      leaf_scores <- famscore_bdeu_byrow(tmp, ess, r, q, s = tabulate(conf[, i]+1))
      score  <- sum(leaf_scores)
      if (score-best_score > 10e-16) {
        best_score <- score
        best_split$var <- i
        best_split$vals <- vals
        best_split$leaf_scores <- leaf_scores
      }
    }
  }
  return(best_split)
}
