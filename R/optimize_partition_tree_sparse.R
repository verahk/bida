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
#' dnames <- list(y = 0:1, x = 0:1, z = 0:1)
#' counts <- bida_sparse_array(c(10, 100, 10, 100, 10, 1000, 10, 1000),
#'                             index = 0:7,
#'                             dim = lengths(dnames),
#'                             dimnames = dnames)
#' ess <- 1
#' tree <- optimize_partition_tree_sparse(counts, ess, min_score_improv = 0, prune = F, verbose = TRUE)
#' print(tree)
#' tree <- optimize_partition_tree_sparse(counts, ess, min_score_improv = -Inf, prune = F, verbose = TRUE)
#' print(tree)
#' tree <- optimize_partition_tree_sparse(counts, ess, min_score_improv = -Inf, prune = F, verbose = TRUE)
#' print(tree)
#' summary.tree(tree)
#' predict.tree(tree, expand.grid(dnames))
optimize_partition_tree_sparse <- function(counts, ess, min_score_improv = 0, prune = F, simplify = TRUE, verbose = FALSE) {

  vars <- names(dimnames(counts))
  if (is.null(vars)) stop("Counts must have named dimnames")
  dims <- dim(counts)
  r <- dims[1]
  q <- prod(dims)/r

  grow_tree <- function(counts, score) {
    # function for growing a tree recursively
    best_split <- find_best_split(counts, score)
    if (!is.null(best_split$var)) {
      # if best_split is a split, and not a leaf
      if (verbose) {
        cat("Score-diff:", sum(best_split$leaf_scores)-score,
            "Splitvariable:", best_split$var,
            "Splitvalue:", best_split$vals,
            "Leaf-scores:", best_split$leaf_scores,
            "\n")
      }
      # grow new trees for each split value
      best_split$branches <- vector("list", length(best_split$values))
      names(best_split$branches) <- best_split$values
      tmp <- list(counts = asplit(counts, best_split$var),
                  score = best_split$leaf_scores)
      for (b in seq_along(best_split$branches)) {
        best_split$branches[[b]] <- grow_tree(counts = tmp$counts[[b]], score  = tmp$score[[b]])
      }
    }
    return(best_split)
  }
  find_best_split <- function(counts, score) {
    dims <- dim(counts)
    best_split <- list(score = score, size = prod(dims[-1]))
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


  # fit ----
  score <- famscore_bdeu_1row(as.array(rowSums(counts)), ess)
  tree  <- grow_tree(counts, score)

  get_from_leaves <- function(tree, name) {
    # collect object `name` from each leaf in `tree`
    if (is.null(tree$branches)) unname(tree[name])
    else unlist(lapply(tree$branches, get_from_leaves, name = name), recursive = F)
  }

  # prune ----
  if (prune) {
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
        list(score = tree$score,
             size  = tree$size)  # return root as a leaf
      }
    }

    # prune tree
    tree <- prune_tree(tree)
  }

  # keep only information about split variable and value in non-leaf nodes
  # - kept for making split nodes into leaf during pruning..
  if (simplify) {
    simplify_tree <- function(tree) {
      if (is.null(tree$branches)) return(tree)
      else {
        tree$branches <- lapply(tree$branches, simplify_tree)
        tree[c("var", "values", "branches")]
      }
    }
    tree <- simplify_tree(tree)
  }

  # enumerate nodes ----
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

  structure(tree,
            class = c("partition", "tree"))

}



#' @export
print.tree <- function(tree) {

  tmp <- function(tree, prefix = "") {
    if (is.null(tree$branches)) {
      cat(prefix, "part:", tree$part, "score:", tree$score, "size:", tree$size, "\n")
    } else {
      cat(prefix, "-+", tree$var, ":", tree$values, "\n")
      for (branch in tree$branches) tmp(branch, prefix = paste0(prefix, "  |"))
    }
  }

  tmp(tree)
}


#' @rdname optimize_partition_tree
#' @export
summary.tree <- function(tree, names = c("part", "score", "size"), prettify = TRUE) {
  get_leaf <- function(tree, name) {
    if (is.null(tree$branches)) return(tree[names])
    unname(unlist(lapply(tree$branches, get_leaf, name = name), recursive = FALSE))
  }

  tmp <- get_leaf(tree, name)
  if (prettify) {
    out <- matrix(unlist(tmp, recursive = FALSE), ncol = length(names), byrow = TRUE)
    colnames(out) <- names
    out
  } else {
    tmp
  }
}

#' @rdname optimize_partition_tree
#' @export
score.tree   <- function(tree) {
  sum(unlist(summary.tree(tree, name = "score", prettify = FALSE)))
}


#' @rdname optimize_partition_tree
#' @export
predict.tree <- function(fit, newdata, name = "part") {

  get_leaf <- function(tree, coordinates, name) {
    # collect object `name` from each leaf in `tree`
    if (is.null(tree$branches)) {
      if (length(name) == 0 || nchar(name) == 0) return(tree)
      else return(unname(tree[[name]]))
    }
    val <- as.character(coordinates[tree$var])
    if (is.na(val)) stop("could not find value of split variable ", tree$var)
    get_leaf(tree$branches[[val]], coordinates, name)
  }


  if (is.null(dim(newdata))) {
    get_leaf(tree, newdata, name)
  } else {
    apply(newdata, 1, get_leaf, tree = tree, name = name)
  }
}




