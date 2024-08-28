

#' Draw a random partition of an outcome space
#'
#' Draw a random partition of the outcome space of a set of categorical variables
#'
#' @param nlev cardinality of the variables.
#' @param prob parameter controlling the size of the partition.
#' @param method name of algorithm.
#' @param regular If `TRUE`, the partition is forced to be regular. See [bida::make_regular].
#' @return an vector of length `prod(nlev)` assigning each joint outcome to a
#'  subset of the partition.
#' @export
#'
#' @examples
#'
#' set.seed(007)
#' nlev <- c(2, 2, 2)
#' rand_partition(nlev, prob = .5)
#'
#' set.seed(007)
#' nlev <- c(2, 2, 2)
#' rand_partition(nlev, .5, method = "tree")
#'
#' # force at least one split by manipulating `nextsplitprob` rule
#' set.seed(007)
#' nextsplitprob <- function(x) x/2
#' rand_partition(nlev, 1, "tree", nextsplitprob)
#'
#' stopifnot(length(unique(rand_partition(nlev, 1, "tree", function(x) 0))) == 2)
#'
#'
#' # method = "`dgraph" to merge leaves in decision tree
#' set.seed(007)
#' rand_partition(nlev, 1, "dgraph", nextsplitprob)
#'
#' # high cardinality
#' nlev <- rep(8, 4)
#' parts <- list(labels = rand_partition(nlev, .5, "labels"),
#'               tree = rand_partition(nlev, .5, "tree"),
#'               dgraph = rand_partition(nlev, .5, "dgraph"))
#' lapply(parts, function(x) length(unique(x))/length(x))
#'
#' parts <- rand_partition(nlev, .5, "dgraph")
#' P <- split(seq_along(parts), parts)
#' is_regular(P, nlev)
rand_partition <- function(nlev,
                           prob,
                           method = "tree",
                           regular = FALSE,
                           ...) {
  P <- switch(method,
             "labels" = partition_from_labels(rand_labels(nlev, prob, ...), nlev),
             "tree" = rand_partition_tree(nlev, prob, doMerge = F, ...),
             "dgraph" = rand_partition_tree(nlev, prob, doMerge = T, ...))

  if (regular) {
    P <- split(seq_along(P), P)
    return(get_parts(make_regular(P, nlev)))
  } else {
    # ensure parts are enumerated from 1,..., length(unique(P))
    match(P, unique(P))
  }


}


rand_labels <- function(nlev, lprob) {

  n  <- length(nlev)      # number of nodes
  q  <- prod(nlev)        # number of joint outcomes
  joint  <- seq_len(q)    # joint outcomes

  stride <- c(1, cumprod(nlev[-n]))
  levels <- lapply(nlev-1, seq.int, from = 0)
  conf   <- expand_grid_fast(levels)
  labels <- vector("list", length(nlev))

  for (i in seq_len(n)) {
    if (runif(1) < lprob) {
      nlabels <- rbinom(1, q/nlev[i]-1, lprob)
      contexts <- sample(joint[conf[, i] == 0], nlabels, FALSE)
      labels[[i]] <- conf[contexts, -i, drop = FALSE]
    }
  }
  labels
}



#' @param splitprob (numeric constant) probability of splitting a node
#' @param doMerge (logical constant) randomly merge leaves in tree
#' @param nextsplitprob (function) manipulate `splitprob` for each new split
rand_partition_tree <- function(nlev, splitprob, doMerge = TRUE, nextsplitprob = function(x) x) {

  # define routine for growing a tree
  grow_tree <- function(splitprob, vars, subset) {
    if (length(vars) == 0 || runif(1) > splitprob) {
      return(list(subset = subset))
    }

    # sample a split variable
    pos <- sample.int(length(vars), 1)
    x <- vars[pos]

    # split the current subset by values of x
    xval <- (subset%/%stride[x])%%nlev[x]
    new_subsets <- unname(split(subset, xval))

    # grow a new tree for each value of the split variable
    list(var = x,
         branches = lapply(new_subsets,
                           function(y) grow_tree(nextsplitprob(splitprob), vars[-pos], y)))
  }


  # define routine for extracting subsets in each leaf
  unlist_tree <- function(tree) {
    if (is.null(tree$subset)) {
      unlist(lapply(tree$branches, unlist_tree), recursive = FALSE)
    } else {
      unname(tree["subset"])
    }
  }

  stride <- c(1, cumprod(nlev[-length(nlev)]))
  vars   <- seq_along(nlev)
  subset <- seq_len(prod(nlev))-1
  tree   <- grow_tree(splitprob, vars, subset)
  partition <- unlist_tree(tree)
  nparts <- length(partition)

  if (doMerge && nparts > 1) {
    nparts <- length(partition)
    tmp <- sample.int(nparts, nparts, TRUE)  # assign parts to new parts
    partition <- lapply(split(partition, tmp), do.call, what = "c")
  }

  get_parts(partition)
}




