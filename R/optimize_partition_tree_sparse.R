#' Title
#'
#' @param counts
#' @param score
#' @param coord
#'
#' @return
#' @export
#'
#' @examples
#'
#' # binary split data with missing values of outcome y
#' data <- cbind(z = rep(0:2, 9),
#'               x = rep(0:2, each = 9),
#'               y = c(rep(0, 9), rep(1, 2*9)))
#' nlev <- rep(3, 3)
#'
#' optimize_partition_from_data_tree_sparse(
#'   data,
#'   j = 3,
#'   parentnodes = 1:2,
#'   ess = 1,
#'   nlev = rep(3, 3),
#'   min_score_improv = -Inf,
#'   prune = T,
#'   verbose = TRUE
#' )
#'
optimize_partition_tree_sparse <- function(counts, ess, min_score_improv = 0, prune = F, verbose = FALSE) {

  dims <- dim(counts)
  r <- dims[1]
  q <- prod(dims)/r

  grow_tree <- function(counts, score, coord) {
    # function for growing a tree recursively
    best_split <- find_best_split(counts, score, coord)
    if (length(best_split) > 2) {
      # if best_split is a split, and not a leaf
      if (verbose) {
        cat("Score-diff:", sum(best_split$leaf_scores)-score,
            "Splitvariable:", best_split$var,
            "Splitvalue:", best_split$vals,
            "Leaf-scores:", best_split$leaf_scores,
            "\n")
      }
      # grow new trees for each split value
      add_coord <- function(val, name) c(coord, setNames(val, name))
      best_split$branches <- mapply(grow_tree,
                                    counts = asplit.bida_sparse_array(counts, best_split$var),
                                    score = best_split$leaf_scores,
                                    coord  = lapply(best_split$values, add_coord, name = best_split$var),
                                    SIMPLIFY = FALSE)

    }
    return(best_split)
  }
  find_best_split <- function(counts, score, coord) {
    best_split <- list(score = score, coord = coord)
    dims <- dim(counts)
    if (length(dims) > 1 && length(counts$value) > 0) {
      best_score <- score + min_score_improv      # minimum value for score of new split
      size   <- prod(dims[-1])                    # size of current part
      stride <- get_stride(counts)
      y <- get_coordinates(counts, 1, stride)
      for (i in seq_along(dims)[-1]) {

        # sum observations for each level of y and variable i
        group  <- get_coordinates(counts, i, stride) + dims[i]*y
        ugroup <- seq_len(r*dims[i])-1
        tmp <- matrix(rowsum_fast(counts$value, group, ugroup), ncol = r)

        # compute score of each leaf and compare against current best
        leaf_scores <- famscore_bdeu_byrow(tmp, ess, r, q, s = size/dims[i])
        score <- sum(leaf_scores)
        if (score > best_score) {
          best_score <- score
          best_split$var    <- names(counts$dimnames[i])
          best_split$values <- counts$dimnames[[i]]
          best_split$leaf_scores <- leaf_scores
        }
      }
    }
    return(best_split)
  }
  get_from_leaves <- function(tree, name) {
    # collect object `name` from each leaf in `tree`
    if (is.null(tree$branches)) unname(tree[name])
    else unlist(lapply(tree$branches, get_from_leaves, name = name), recursive = F)
  }

  prune_tree <- function(tree) {
    # recursive function for pruning trees
    if (is.null(tree$branches)) return(tree)
    for (b in seq_along(tree$branches)) {
      tree$branches[b] <- list(prune_tree(tree$branches[[b]]))
    }

    # compare score of pruned tree to root
    score <- sum(unlist(get_from_leaves(tree, "score")))
    if (score > tree$score) {
      tree
    } else {
      if (verbose) cat(sprintf("Collapse split on variable %s. Score-diff: %s \n",
                       tree$var, tree$score-score))
      list(score = tree$score, coord = tree$coord)  # return root as a leaf
    }
  }

  # fit ----
  tree <- grow_tree(counts, score, coord = integer(0))

  if (prune) {
    # prune tree
    tree <- prune_tree(tree)
  }

  # enumerate nodes
  enumerate_leaves <- function(tree) {
    if (is.null(tree$branches)) {
      part <<- part + 1
      tree$part <- part
    } else {
      tree$branches <- lapply(tree$branches, enumerate_leaves)
    }
    return(tree)
  }
  part <- 0
  tree <- enumerate_leaves(tree)
  return(list(coord = get_from_leaves(tree, "coord"),
              scores = get_from_leaves(tree, "score"),
              tree = enumerate_leaves(tree)))
}

optimize_partition_from_data_tree_sparse <- function(data, j, parentnodes, ess, nlev, min_score_improv = 0, prune = F, verbose = FALSE){
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("X", seq_len(ncol(data)))
  }

  r <- nlev[j]
  q <- prod(nlev[parentnodes])

  vars   <- colnames(data)[c(j, parentnodes)]
  counts <- counts_from_data_matrix(data[, vars], nlev, TRUE)
  counts$dimnames <- setNames(lapply(counts$dim-1, seq.int, from = 0), vars)

  tmp   <- tabulate(data[, j]+1, r)
  score <- famscore_bdeu_1row(matrix(tmp, 1, r), ess, r, q)

  optimize_partition_tree_sparse(counts, ess, min_score_improv, prune, verbose)
}

get_part <- function(tree, coordinates) {
  if (is.null(tree$branches)) {
    tree$part
  } else {
    val <- coordinates[tree$var]
    get_part(tree$branches[[val+1]], coordinates)
  }
}

