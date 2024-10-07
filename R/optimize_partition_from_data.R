

#' Optimize the local structure of a conditional probability table (CPT)
#'
#' Given data on a response variable and its parents, all assumed categorical,
#' learn a partitioning over the space of the parent variables (the rows of the CPT).
#'
#' @param data (numeric matrix) matrix with integers representing zero-based categorical variables.
#' @param j (integer) column position of response variable.
#' @param parentnodes (integer vector) column positions of parent nodes.
#' @param ess (numeric) equivalent sample size for the BDeu-score.
#' @param nlev (integer vector) cardinality of every variable in `data`
#' @param method (character) name of the algorithm used to learn a partition
#' @param verbose (bolean)
#' @param obj an object of class `partition`
#' @param newdata a matrix or data.frame with values of the `parentnodes`
#' @details
#' The `method` argument must match one of the following algorithms:
#' - `tree`: fit a decision tree using a greedy algorithm, where each
#'   split adds a new branch for each level of the split variable and splits are
#'   added until no new splits improves the score.
#' - `ptree`: fit a decision tree as above, but grow a tree of full depth
#'   (by searcing for the split that gives the least worsening of the score) and
#'   then prune the tree by removing all subtrees that do not improve the score.
#' - `treereg` or `ptreereg`: as the above, yet the `reg` suffix will force all
#'   parents to be present in the tree and the partition to be regular. If there
#'   are parents that are not split variables in the tree, splits on each of
#'   those parents are added to *every* leaf node.
#' - `pcart`: fit the optimial binary decision tree in terms of the BDeu-score,
#'   using the algorithm implemented in `rpcart::opt.pcart.cat.bdeu`.
#'
#' @section Methods for the class `partition`:
#' - [get_rules()]: returns a list with [base::expression()] that defines each
#'   part of the partition in terms of logical rules, e.g. `x == 1 & y == 2`.
#' - [get_predictors()]: returns the subset of `parentnodes` that is relevant
#'   for the partition.
#' - [predict.partition()]: assigns each row in `newdata` to a part in the partition.
#'
#' @return an object of class `partition` with subclass `tree` or `pcart`.
#' @export
#'
#' @examples
#'
#' # Example 1: binary split data with missing levels
#' data <- cbind(z = rep(0:2, 9),
#'               x = rep(0:2, each = 9),
#'               y = c(rep(0, 9), rep(1, 2*9)))
#' nlev <- rep(3, 3)
#' names(nlev) <- colnames(data)
#'
#' j <- 3
#' parentnodes <- 1:2
#' ess <- 1
#' newdata <- expand.grid(lapply(nlev[parentnodes]-1, seq.int, from = 0))
#'
#' # tree: greedy decision tree
#' fit <- optimize_partition_from_data(data, j, parentnodes, ess, nlev, "tree", verbose = TRUE)
#' print(fit)
#' predict(fit, newdata)
#'
#' # treereg: regular greedy decision tree
#' fit <- optimize_partition_from_data(data, j, parentnodes, ess, nlev, "treereg")
#' print(fit)
#' predict(fit, newdata)
#'
#' # pcart: score-maximizing binary decision tree
#' fit <- optimize_partition_from_data(data, j, parentnodes, ess, nlev, "pcart")
#' print(fit)
#' predict(fit, newdata)
#'
#' # Example 2: no best split
#' p <- c(.9, .1, .1, .9)
#' z <- runif(100) > .5
#' x <- runif(100) > .5
#' y <- runif(100) > p[z + 2*x+1]
#' table(z+2*x, y)  # frequency table
#'
#' data <- cbind(z = z, x = x, y = y)*1
#' nlev <- c(2, 2, 2)
#'
#' # tree: greedy decision tree
#' fit <- optimize_partition_from_data(data, j, parentnodes, ess, nlev, "tree", verbose = TRUE)
#' get_rules(fit)
#'
#' # ptree: pruned greedy decision tree
#' fit <- optimize_partition_from_data(data, j, parentnodes, ess, nlev, "ptree", verbose = TRUE)
#' print(fit)
#'
#' # list rules
#' get_rules(fit)
#'
#' #' # Example 3:
#' p <- c(.9, .9, .9, .1, .1, .7)
#' z <- sample(0:2, 1000, replace = T)
#' x <- sample(0:1, 1000, replace = T)
#' joint <- z + 2*x
#' y <- runif(length(z)) > p[joint+1]
#' table(joint, y)  # frequency table
#'
#' data <- cbind(z = z, x = x, y = y)*1
#' nlev <- c(3, 2, 2)
#' j <- 3
#' parentnodes <- 1:2
#'
#' newdata <- expand.grid(list(z = 0:2, x = 0:1))
#' fit <- optimize_partition_from_data(data, j, parentnodes, ess, nlev, "tree", verbose = TRUE)
#' fit
#' cbind(newdata, predict(fit, newdata))
#' fit <- optimize_partition_from_data(data, j, parentnodes, ess, nlev, "ptree", verbose = TRUE)
#' fit
#' cbind(newdata, predict(fit, newdata))
#'
#' fit <- optimize_partition_from_data(data, j, parentnodes, ess, nlev, "pcart", verbose = TRUE)
#' print(fit)
#' cbind(newdata, predict(fit, newdata))
#'
#' # methods
#' methods(class = "partition")
#' methods(class = "tree")
#' methods(class = "pcart")
optimize_partition_from_data <- function(data, j, parentnodes, ess, nlev, method, verbose = FALSE) {


  if (is.null(colnames(data))) colnames(data) <- paste0("X", seq_len(ncol(data)))

  if (method == "pcart") {
    # convert data to data frame with factors with level corresponding to nlev
    tmp <- lapply(c(j, parentnodes), function(i) factor(data[, i], seq.int(0, nlev[i]-1)))
    names(tmp) <- colnames(data)[c(j, parentnodes)]
    df <- data.frame(tmp)
    return(optimize_partition_from_df_pcart(df, ess, use_structure_score = FALSE))
  }

  vars <- c(j, parentnodes)
  data <- data[, vars]
  nlev <- nlev[vars]
  names(nlev) <- colnames(data)

  if (method == "tree" || method == "treereg") {
    optimize_partition_from_data_tree(data, ess, nlev, min_improv = 0, prune = FALSE, regular = endsWith(method, "reg"), verbose = verbose)
  } else if (method == "ptree" || method == "ptreereg") {
    optimize_partition_from_data_tree(data, ess, nlev, min_improv = -Inf, prune = TRUE, regular = endsWith(method, "reg"), verbose = verbose)
  } else {
    stop("invalid `method` argument")
  }
}

#' @rdname optimize_partition_from_data
#' @export
get_predictors <- function(obj, ...) {
  UseMethod("get_predictors")
}
#' @noRd
#' @export
get_predictors.partition <- function(obj) {
  if (inherits(obj, "tree") || inherits(obj, "pcart")) {
    unique(stringr::str_extract_all(fit, "(?<=-\\+ )[^:]*")[[1]])
  } else {
    NextMethod()
  }
}

#' @rdname optimize_partition_from_data
#' @export
get_rules <- function(obj, ...) {
  UseMethod("get_rules")
}


#' @rdname optimize_partition_from_data
#' @export
predict.partition <- function(obj, newdata) {
  if (inherits(obj, "tree") || inherits(obj, "partition")) {
    if (is.null(dim(newdata))) newdata <- rbind(newdata)
    rules <- get_rules(obj)
    if (length(rules) == 1) {
      parts <- rep(1L, nrow(newdata))
    } else {
      newdata <- as.data.frame(newdata)
      parts   <- vector("integer", nrow(newdata))
      for (part in seq_along(rules)) {
        indx <- with(newdata, eval(rules[[part]]))
        parts[indx] <- part
      }
    }
  } else stop("only applicable to objects of class `tree` or `pcart`")

  return(parts)
}

#' @rdname optimize_partition_from_data
#' @export
print.partition <- function(obj) {
  cat("Partition of subclass", class(obj)[-1], "\n")
  cat(obj, "\n")
  cat("score: ", attr(obj, "score"), "\n")
  cat("size of parts: ", attr(obj, "sizes"), "\n")
}


# decision trees ----

#' Greedy optimization of local structure: decision trees
#'
#' This function learn a partition of the parent space, using a greedy recursive
#' decision tree algorithm. At every non-leaf node, one branch is added to the tree
#' for every level of the split variable. The score of the tree is the sum of its
#' leaf-scores. The next split is found by iterating through every variable
#' that is yet not a split variable in the current branch and compute the score
#' of the new leaves.
#'
#' @keywords internal
#' @param data (numeric matrix) a matrix with integers representing zero-based
#'  categorical cariables. The response variable is assumed placed in the first
#'  column, the predictors in the remaining.
#' @param ess ()
#' @param min_improv (numeric) minimum score improvement for the next split.
#'  If 0, splits are added greedily until no split improves the score.
#'  If `-Inf`, full trees are grown, choosing at each iteration
#'  the split that gives the least worsening of the score.
#' @param prune (bolean) if `TRUE` the tree is pruned, that is all subtrees that
#'  do not improve the score is removed. Defaults to `FALSE`.
#' @param regular (bolean) if `TRUE`, the partition is forced regular by adding
#'  to each leaf a split for each variable that is not present in the fitted tree.
#'  Defaults to `FALSE`.
#' @param verbose (bolean) if `TRUE`, messages is printed during the fitting process.
#' @seealso [optimize_partition_from_data()]
#' @return an object of class `partition` with sub-class `tree`
optimize_partition_from_data_tree <- function(data, ess, nlev, min_improv, prune = FALSE, regular = FALSE, verbose = FALSE) {
  r <- nlev[[1]]
  q <- prod(nlev[-1])

  grow_tree <- function(leaf, bins, force_split = FALSE) {
    # function for growing a tree recursively

    dims <- dim(bins)
    find_split <- dims[2] > 0 && (dims[1] > 1 | force_split)
    if (find_split) {
      # find best split, if possible
      best_split <- leaf
      best_score <- leaf$score + min_improv
      for (v in colnames(bins)) {
        # sum observations for each level of y and variable i
        tab <- matrix(tabulate(bins[, v], nbins[v]), ncol = r, byrow = TRUE)

        # compute score of each leaf and compare against current best
        leaf_scores <- famscore_bdeu_byrow(tab, ess, r, q, s = leaf$size/nlev[[v]])
        score <- sum(leaf_scores)
        if (score > best_score) {
          best_score <- score
          best_split$var  <- v
          best_split$leaf_scores <- leaf_scores
          best_split$tab <- tab
        }
      }
    }

    if (!find_split || is_leaf(best_split)) {
      return(leaf)
    } else {
      # if a split was found, keep growing the tree

      if (verbose) {
        sprintf("Score-diff: %1.2f Splitvariable: %s Leaf-scores: %s\n",
                sum(best_split$leaf_scores)-score,
                best_split$var,
                paste(round(best_split$leaf_scores, 2), collapse = " "))
      }

      v <- best_split$var
      col <- match(v, colnames(bins)) # position of col - drop in next iter
      dv  <- (bins[, col]-1)%/%r      # value of v, to split data
      k <- nlev[[v]]                  # cardinality of k - add a branch for each level
      best_split$branches <- vector("list", k)
      for (b in seq_len(k)) {
        new_leaf <- list(score = best_split$leaf_scores[[b]],
                         size  = best_split$size/k,
                         counts = best_split$tab[b, ])
        new_bins <- bins[dv == b-1, -col, drop = FALSE]
        best_split$branches[[b]] <- grow_tree(new_leaf, new_bins, force_split)
      }
      return(best_split[keepinsplit])
    }
  }

  is_leaf <- function(tree) length(tree) == 3
  get_leaf_scores <- function(tree) {
    if (is_leaf(tree)) tree$score
    else lapply(tree$branches, get_leaf_scores)
  }
  get_leaf_sizes <- function(tree) {
    if (is_leaf(tree)) tree$size
    else lapply(tree$branches, get_leaf_sizes)
  }
  # get_leaf_counts <- function(tree) {
  #   if (is_leaf(tree)) tree$counts
  #   else unlist(lapply(tree$branches, get_leaf_counts), recursive = FALSE)
  # }


  # fit ----
  # init root node as leaf
  counts <- tabulate(data[, 1]+1, r)
  leaf <- list(counts = counts, score = famscore_bdeu_1row(counts, ess, r), size = q)

  # enumerate joint outcome (y, x) for each x for faster tabulation
  bins <- data[, 1] + r*data[, -1, drop = FALSE] +1
  nbins <- r*nlev[-1]

  # grow tree
  keepinsplit  <- c("var", "score", "size", "counts", "branches")
  tree  <- grow_tree(leaf, bins)

  # prune ----
  if (prune) {
    prune_tree <- function(tree) {
      # recursive function for pruning trees
      if (is_leaf(tree)) return(tree)
      tree$branches <- lapply(tree$branches, prune_tree)

      # compare score of pruned tree to root
      score <- sum(unlist(get_leaf_scores(tree)))
      if (score > tree$score) {
        tree
      } else {
        if (verbose) cat("Collapse split on variable:", tree$var,
                          "Score-diff:", round(tree$score-score, 2),
                         "\n")
        tree[keepinleaf] # return root as a leaf
      }
    }

    # prune tree
    keepinleaf <- c("score", "size", "counts")
    tree <- prune_tree(tree)
  }


  # make regular ----
  if (regular) {
    get_splitvars <- function(tree) {
      if (!is_leaf(tree)) c(tree$var, unlist(lapply(tree$branches, get_splitvars)))
      else (character(0))
    }

    add_split <- function(tree, bins) {
      if (is.null(tree$var)) {
        pos  <- match(new_splitvars, colnames(bins))
        tree <- grow_tree(bins[, pos, drop = FALSE],
                          tree$score,
                          tree$size,
                          tree$counts,
                          force_split = TRUE)
      } else {
        v <- tree$var
        col <- match(v, colnames(bins))
        dv <- bins[, col]%/%nlev[[1]]
        for (b in seq_along(tree$branches)) {
          tree$branches[[b]] <- add_split(tree = tree$branches[[b]],
                                          bins = bins[dv == b-1, -col, drop = FALSE])
        }
      }
      return(tree)
    }

    min_improv <- -Inf
    new_splitvars <- setdiff(names(nlev)[-1], unique(get_splitvars(tree)))
    if (length(new_splitvars) > 0) {
      if (verbose) cat("Add splits for every non-split variable to all leaves\n")
      tree <- add_split(tree, bins)
    }
  }

  # output ----
  to_string <- function(tree, prefix = "") {
    if (is.null(tree$branches)) {
      sprintf("%s-- score: %1.2f size: %s\n", # counts: %s\n",
              prefix, tree$score, tree$size) # paste(tree$counts, collapse = " "))
    } else {
      split <- sprintf("%s-+ %s:\n", prefix, tree$var)
      c(split, unlist(lapply(tree$branches, to_string, prefix = paste0(prefix, " | "))))
    }
  }

  obj <- paste0(to_string(tree), collapse = "")
  structure(obj,
            score = sum(unlist(get_leaf_scores(tree))),
            sizes = unname(unlist(get_leaf_sizes(tree))),
            class = c("partition", "tree"))

}


#' @rdname optimize_partition_from_data_tree
#' @export
get_rules.tree <- function(obj) {
  rules_from_tree <- function(tree, rule = "") {
    # recurxive function that returns the rules defining each leaf of the tree
    if (grepl("^\\-\\+", tree[1])) {
      root <- tree[1]
      var <- substring(root, 4, nchar(root)-1)

      # collect subtrees
      tmp <- gsub("^ \\| ", "", tree[-1])         # remove pre-fix
      pos <- grep("^-+|^--", tmp)                 # positions of next nodes or leaves
      subtrees <- split(tmp, rep.int(seq_along(pos), diff(c(pos, length(tmp)+1))))

      # update rules
      new_rules <- lapply(seq_along(pos)-1,
                          function(x) paste0(rule, " & ", var, "==", x))

      # add rules of each subtree
      unlist(lapply(seq_along(subtrees),
                    function(b) rules_from_tree(subtrees[[b]], new_rules[[b]])),
             recursive = FALSE)

    } else {
      list(parse(text = gsub("^ & ", "", rule)))
    }
  }
  rules_from_tree(strsplit(obj, "\n")[[1]])
}

#' @rdname optimize_partition_tree
#' @export
summary.tree <- function(obj) {
  tree <- strsplit(obj, "\n")[[1]]
  leaves <- tree[grepl("--", tree)]
  tmp <- list(part = seq_along(leaves),
              score = gsub(".*score: (-?[0-9]+\\.[0-9]+).*", "\\1", leaves),
              size  = gsub(".*size: ([0-9]+).*", "\\1", leaves))
  tmp <- lapply(tmp, as.numeric)
  do.call("cbind", tmp)
}



# pcart ----

#' Greedy optimization of local structure using `rpcart`
#'
#' This function learn a partition of the parent space of an variable,
#' defined by the binary decision tree that maximizes the score.
#' It is a wrapper around `rpcart:::opt.pcart.cat.bdeu` that fits the optimal binary tree.
#' If the variables in `df` are not factor variables, they will be converted to factors
#' with levels corresponding to those observed. Note that the score will depend
#' on wheter or not missing levels are included or not.
#'
#' @keywords internal
#' @param df (data frame) an data frame with factors.
#'  The first column includes the repsonse variable, the remaining predictor variables.
#' @seealso [optimize_partition_from_data()]
#' @return an object of class `partition` with sub-class `pcart`
optimize_partition_from_df_pcart <- function(df, ess, verbose = FALSE, use_structure_score = FALSE) {

  predictors <- names(df)[-1]
  response   <- names(df)[1]
  if (verbose) cat("\nOptimize local structure using rpcart::opt.pcart.cat.bdeu\n")
  fit <- rpcart:::opt.pcart.cat.bdeu(df, predictors, response, ess, use_structure_score = use_structure_score)

  # compute size of each part
  nlev <- vapply(df[-1], nlevels, integer(1))
  tree <- strsplit(fit$tree, "\n")[[1]]
  compute_part_sizes <- function(tree, size = prod(nlev)) {
    if (startsWith(tree[1], "-+")) {
      root <- tree[1]
      split <- strsplit(substring(root, 4), ": | \\| ")[[1]]
      var <- split[1]
      kTRUE  <- length(strsplit(split[2], " ")[[1]])
      kFALSE <- length(strsplit(split[3], " ")[[1]])
      k <- kTRUE + kFALSE
      subtrees <- split(substring(tree[-1], 3), grepl("^ \\|", tree[-1]))
      c(compute_part_sizes(subtrees[[2]], size/k*kTRUE),
        compute_part_sizes(subtrees[[1]], size/k*kFALSE))
    } else return(size)
  }
  compute_part_sizes(tree, prod(nlev))

  structure(fit$tree,
            score = fit$score,
            sizes = compute_part_sizes(tree),
            class = c("partition", "pcart"))
}

#' @rdname optimize_partition_from_df_pcart
#' @export
get_rules.pcart <- function(obj) {
  rules_from_tree <- function(tree, rule = "") {
    if (grepl("^\\-\\+", tree[1])) {
      root <- tree[1]
      split <- strsplit(substring(root, 4), ": | \\| ")[[1]]

      splitvals <- lapply(split[-1],
                          function(x) strsplit(x, " +")[[1]])
      new_rules <- vapply(splitvals,
                          function(x) paste0(paste0(split[1], "==", x), collapse = "|"),
                          character(1))
      upd_rules <- paste(rule, new_rules, sep = " & ")

      subtrees <- split(substring(tree[-1], 3), grepl("^ \\|", tree[-1]))
      c(rules_from_tree(subtrees[[2]], upd_rules[[1]]),
        rules_from_tree(subtrees[[1]], upd_rules[[2]]))
    } else {
      list(parse(text = gsub("^ & ", "", rule)))
    }
  }
  tree <- strsplit(obj, "\n")[[1]]
  rules_from_tree(tree)
}


if (FALSE) {
  # profile ----
  ## compare tree-optimizer from data and from counts

  data <- structure(c(0, 2, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0,
                      1, 2, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0,
                      1, 2, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0,
                      0, 0, 1, 1, 0, 2, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1,
                      1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 2,
                      1, 2, 1, 1, 1, 0, 0, 1, 2, 1, 2, 2, 1, 1, 1, 2, 1, 2, 2, 2, 2,
                      2, 2, 0, 1, 2, 1, 2, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 2, 2, 1,
                      1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 0, 0, 2, 0, 1, 2, 1, 1, 1, 0,
                      1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 1, 1, 1, 1, 2, 2, 1, 1, 2, 0, 1,
                      1, 1, 1, 1, 0, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 2, 0, 1, 1,
                      0, 1, 1, 1, 0, 0, 2, 2, 1, 1, 0, 1, 1, 2, 1, 0, 0, 0, 0, 1, 0,
                      1, 0, 0, 2, 1, 1, 1, 1, 1, 0, 1, 0, 2, 2, 0, 0, 1, 1, 1, 2, 1,
                      0, 1, 1, 1, 2, 1, 1, 1, 2, 0, 1, 0, 0, 0, 2, 0, 0, 2, 0, 0, 2,
                      1, 0, 2, 2, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 2, 1, 1, 1, 2, 1, 0,
                      0, 1, 2, 1, 1, 1, 2, 1, 2, 1, 1, 0, 1, 0, 1, 2, 1, 0, 1, 2, 1,
                      2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      2, 2, 0, 2, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1,
                      1, 2, 2, 1, 2, 1, 1, 1, 2, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
                      1, 1, 0, 0, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 1, 2, 1, 1,
                      1, 1, 1, 2, 1), dim = c(100L, 4L), dimnames = list(NULL, c("Z1",
                                                                                 "Z2", "X", "Y")))

  f <- function(data, j, parentnodes, ess, nlev, method, verbose) {
    counts <- counts_from_data_matrix(data[, c(j, parentnodes)], nlev[c(j, parentnodes)])
    tab <- matrix(counts, ncol = nlev[j], byrow =T)
    levels <- lapply(nlev[parentnodes]-1, seq.int, from = 0)
    names(levels) <- colnames(data)[-j]
    optimize_partition_tree(tab, levels, ess, -Inf, TRUE, verbose)
  }

  j <- 4
  parentnodes <- 1:3
  nlev <- rep(3, 4)
  setNames(nlev) <- colnames(data)

  org <- data
  data <- rbind(data, data, data)
  microbenchmark::microbenchmark(
    optimize_partition_from_data(data, j, parentnodes, ess, nlev, "ptree", FALSE),
    f(data, j, parentnodes, ess, nlev, "ptree", FALSE)
  )

  profvis::profvis(
    replicate(30, optimize_partition_from_data(data, j, parentnodes, ess, nlev, "ptree", FALSE)))
}
