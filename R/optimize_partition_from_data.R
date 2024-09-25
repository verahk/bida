

#' Title
#'
#' @param data
#' @param j
#' @param parentnodes
#' @param nlev
#' @param ess
#'
#' @return
#' @export
#'
#' @examples
#'
#' # binary split data
#' data <- cbind(z = rep(0:2, 9),
#'               x = rep(0:2, each = 9),
#'               y = c(rep(0, 9), rep(1, 2*9)))
#' nlev <- rep(3, 3)
#'
#' # optimize local structure
#' res <- optimize_partition_from_data_pcart(data, 3, 1:2, 1, nlev, verbose = TRUE)
#' parts <- get_parts(res$partition)
#' stopifnot(all(match(parts, unique(parts)) == c(rep(1, 3), rep(2, 6))))
#'
#' # score do not match when there are unobserved levels of outcome variable
#' res[-1]
#'
#' res <- optimize_partition_from_data_pcart(data, 3, 1:2, 1, nlev = NULL, recompute_score = TRUE, verbose = TRUE)
#' res[-1]
#' stopifnot(abs(res$score-res$pcart$score) < 10**-10)
optimize_partition_from_data <- function(data, j, parentnodes, ess, nlev, method, verbose = FALSE, ...) {
  if (method == "pcart") {
    optimize_partition_from_data_pcart(data, j, parentnodes, ess, nlev, verbose = verbose, ...)
  } else {
    stride <- c(1, cumprod(nlev[parentnodes]))
    r <- nlev[j]
    q <- stride[length(parentnodes)+1]
    counts <- matrix(tabulate(data[, c(parentnodes, j)]%*%stride+1, q*r), q, r)
    optimize_partition(counts,
                       levels = lapply(nlev[parentnodes]-1, seq.int, from = 0),
                       ess,
                       method,
                       regular = (method == "ldag"),
                       verbose = verbose,
                       ...)
  }
}
#' Title
#'
#' @param data
#' @param j
#' @param parentnodes
#' @param nlev
#' @param ess
#'
#' @return
#' @export
#'
#' @examples
#'
#' # binary split data
#' data <- cbind(z = rep(0:2, 9),
#'               x = rep(0:2, each = 9),
#'               y = c(rep(0, 9), rep(1, 2*9)))
#' nlev <- rep(3, 3)
#'
#' # optimize local structure
#' res <- optimize_partition_from_data_pcart(data, 3, 1:2, 1, nlev, verbose = TRUE)
#' parts <- get_parts(res$partition)
#' stopifnot(all(match(parts, unique(parts)) == c(rep(1, 3), rep(2, 6))))
#'
#' # score do not match when there are unobserved levels of outcome variable
#' res[-1]
#'
#' res <- optimize_partition_from_data_pcart(data, 3, 1:2, 1, nlev = NULL, recompute_score = TRUE, verbose = TRUE)
#' res[-1]
#' stopifnot(abs(res$score-res$pcart$score) < 10**-10)
optimize_partition_from_data_pcart <- function(data, j, parentnodes, ess, nlev = NULL, recompute_score = !is.null(nlev), verbose = FALSE) {

  if (is.null(nlev)) nlev <- apply(data, 2, function(x) length(unique(x)))

  r <- nlev[j]
  q <- prod(nlev[parentnodes])

  # optimize structure
  df <- as.data.frame(data)
  predictors <- names(df)[parentnodes]
  response   <- names(df)[j]
  if (verbose) cat("\nOptimize local structure using rpcart::opt.pcart.cat.bdeu\n")
  result <- rpcart:::opt.pcart.cat.bdeu(df, predictors, response, ess, use_structure_score = FALSE)
  if (verbose) cat(result$tree)

  # extract partitioning
  tree <- strsplit(result$tree, "\n")[[1]]
  nlev <- nlev[parentnodes]
  stride <- c(1, cumprod(nlev[-length(nlev)]))
  names(nlev) <- names(stride) <- predictors
  unlist_tree <- function(tree, nlev, stride, parts) {
    # recursive function for extracting a partition from a rpcart-tree
    if (startsWith(tree[1], "-+")) {
      root <- tree[1]
      split <- strsplit(substring(root, 4), ": | \\| ")[[1]]
      splitvar <- split[1]
      splitval <- as.numeric(strsplit(split[2], " +")[[1]])

      vals <- (parts%/%stride[splitvar])%%nlev[splitvar]
      subparts <- split(parts, vals %in% splitval)
      subtrees <- split(substring(tree[-1], 3), grepl("^ \\|", tree[-1]))
      if (length(subparts) != length(subtrees)) {
        stop()
      }
      c(unlist_tree(subtrees[[1]], nlev, stride, subparts[[1]]),
        unlist_tree(subtrees[[2]], nlev, stride, subparts[[2]]))
    } else {
      return(list(parts))
    }
  }
  partition <- unlist_tree(tree, nlev, stride, seq_len(prod(nlev))-1)

  # compute score
  if (recompute_score) {
    n_parts <- length(partition)
    data_parts <- get_parts(partition)[data[, parentnodes]%*%stride+1]
    counts <- matrix(tabulate(data_parts + n_parts*data[, j], n_parts*r), n_parts, r)
    score  <- sum(famscore_bdeu_byrow(counts, ess, r, q, lengths(partition)))
  } else {
    score <- result$score
  }

  list(partition = partition,
       score = score,
       pcart = result)
}


if (FALSE) {
  rules_from_tree <- function(tree, rule = "") {
    if (startsWith(tree[1], "-+")) {
      root <- tree[1]
      split <- strsplit(substring(root, 4), ": | \\| ")[[1]]

      splitvals <- lapply(split[-1],
                          function(x) strsplit(x, " +")[[1]])
      add_rules <- lapply(splitvals,
                          function(x) paste0(paste0(split[1], "==", x), collapse = "|"))
      new_rules <- lapply(add_rules,
                          function(x) paste0(rule, x, collapse = "&"))

      subtrees <- split(substring(tree[-1], 3), grepl("^ \\|", tree[-1]))

      c(rules_from_tree(subtrees[[1]], new_rules[[1]]),
        rules_from_tree(subtrees[[2]], new_rules[[2]]))
    } else {
      list(rule)
    }
  }

  rules_chr  <- rules_from_tree(tree, "")
  rules_expr <- lapply(rules_chr, function(x) parse(text = x))
  with(tail(df), eval(rules_expr[[1]]))
  lapply(rules_expr, function(xx) which(with(df, eval(xx))))

}


optimize_partition_from_data_tree <- function(data, j, parentnodes, ess, nlev, min_score_improv = 0){
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


}


optimize_partition_tree_sparse <- function(counts, ess, r, q, min_score_improv) {

  grow_tree_sparse <- function(counts, score, coord) {
    # function for growing a tree recursively
    best_split <- find_best_split(counts, score, coord)
    if (length(best_split) > 2) {
      # if best_split is a split, and not a leaf
      if (verbose) {
        cat("\nScore-diff:", sum(best_split$leaf_scores)-score,
            "Splitvariable:", best_split$var,
            "Splitvalue:", best_split$vals,
            "Leaf-scores:", best_split$leaf_scores)
      }
      # grow new trees for each split value
      add_coord <- function(val, name) c(coord, setNames(val, name))
      best_split$branches <- mapply(grow_tree_sparse,
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
    if (dims > 1 && length(counts$values) > 0) {
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

  prune <- function(tree, verbose = FALSE) {
    # recursive function for pruning trees
    if (is.null(tree$branches)) return(tree)
    for (b in seq_along(tree$branches)) {
      tree$branches[b] <- list(prune_tree(tree$branches[[b]], verbose))
    }

    # compare score of pruned tree to root
    score <- sum(unlist(get_from_leaves(tree, "score")))
    if (score > tree$score) {
      tree
    } else {
      if (verbose) cat(sprintf("Collapse split on variable %s. Score-diff: %s"),
                       tree$var, tree$score-score)
      list(score = tree$score, space = tree$space)  # return root as a leaf
    }
  }

  # fit ----
  tree <- grow_tree(counts, score, coord = integer(0), min_score_improv)

  if (prune) {
    # prune tree
    tree <- prune(tree)
  }

  return(list(space = get_from_leaves(tree, "space"),
              scores = get_from_leaves(tree, "score"),
              tree = tree))
}


