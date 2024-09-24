

#' Title
#'
#' @param dag
#' @param nlev
#' @param splitprob
#' @param depth_ratio
#'
#' @return
#' @export
#'
#' @examples
#'
#' n <- 3
#' dag <- matrix(0, n, n)
#' dag[upper.tri(dag)] <- 1
#' nlev <- rep(2, n)
#'
#' splitprob = 1
#' nextsplitprob = function(x) x
#'
#' sim_rand_partitions(dag, nlev, depth_ratio = 0)
#' sim_rand_partitions(dag, nlev, depth_ratio = 1)
#' sim_rand_partitions(dag, nlev, depth_ratio = .5)
#'
#' sim_rand_bn(n, nlev, depth_ratio = 0)
NULL

sim_rand_bn <- function(n, nlev, depth_ratio, d = min(n-1, 4)) {
  dag <- bida:::randDAG(n, d)

  partitions <- NULL
  if (depth_ratio < 1) {
    # draw a set of partitions for each dag given its parents
    partitions <- sim_rand_partitions(dag, nlev, depth_ratio)
  }

  # draw a bn given dag and partitions / a set of CPTs consistent with dag and partitions
  bida:::rand_bn(dag, "cat", alpha = 1, nlev = nlev, partitions = partitions)
}

sim_rand_partitions <- function(dag, nlev, depth_ratio) {
  n <- ncol(dag)
  partitions <- vector("list", n)
  for (i in seq_len(n)) {
    pa <- which(dag[, i] == 1)
    if (length(pa) > 1) {
      # draw a tree of given depth - use make_regular argument to add splits until all vars are present
      depth <- max(1, ceiling(depth_ratio*length(pa)))
      partitions[[i]] <- sim_rand_partition(nlev[pa], splitprob = 1, mindepth = depth, maxdepth = 1, make_regular = T)
    }
  }
  return(partitions)
}

sim_rand_partition <- function(nlev, splitprob, nextsplitprob = function(x) x, mindepth = 1, maxdepth = length(nlev), make_regular = T) {
  n <- length(nlev)
  stride <- c(1, cumprod(nlev[-length(nlev)]))

  # define routine for growing a tree
  grow_tree <- function(splitprob, vars, subset) {
    if ( (n-length(vars) >= mindepth && n-length(vars) >= maxdepth) || runif(1) > splitprob) {
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
  unlist_tree <- function(tree, name) {
    if (is.null(tree$subset)) {
      unlist(lapply(tree$branches, unlist_tree), recursive = FALSE)
    } else {
      unname(tree[name])
    }
  }

  get_splitvars <- function(tree, name = "var") {
    if (is.null(tree$subset)) {
      c(tree$var, unlist(lapply(tree$branches, get_splitvars), recursive = FALSE))
    }
  }


  vars   <- seq_along(nlev)
  subset <- seq_len(prod(nlev))-1
  tree   <- grow_tree(splitprob, vars, subset)
  P <- unlist_tree(tree)

  if (make_regular) {
   vars_in <- unique(get_splitvars(tree))
   for (x in setdiff(vars, vars_in)) {

     # sample leaf / region to split
     l <- sample(seq_along(P), 1)

     # make a binary split of split l
     subset <- P[[l]]
     xval <- (subset%/%stride[x])%%nlev[x]
     new_subsets <- unname(split(subset, xval))
     P <- c(P[-l], new_subsets)
   }
  }
  return(P)
}
