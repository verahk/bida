

#' Draw a conditional probability table (CPT)
#'
#' Construct a random conditional probability table (CPT) by sampling independently
#' a vector of probabilities from a Dirichlet distribution for each parent configuration.
#' The function returnes an array with named dimension, in the same format
#' as [bnlearn::bn.fit.dnode] store CPTs.
#'
#' @param dim (integer vector) of cardinality of the outcome variable and the
#'  parent variables.
#' @param scope (character vector) names of the variables.
#' @param alpha (integer)
#' @param local_structure (character) name of algorithm to form a partitioning /
#'  produce parameter restrictions.
#'  Defaults to `"none"`, which returns a CPT with no parameter restrictions.
#' @param ... additional arguments sendt to [rand_partition].
#'
#' @return an array with dimension `nlev` and dimension named by `scope`.
#' @keywords internal
#'
#' @examples
#'
#' # no parents
#' dim <- c(2)
#' scope <- c("Y")
#' rand_cpt(dim, scope)
#'
#' # two parents
#' dim <- c(2, 3, 4)
#' scope <- c("Y", "X", "Z")
#' rand_cpt(dim, scope)
#'
rand_cpt <- function(dim, scope = NULL, alpha = 1, local_structure = "none", ...) {
  r <- dim[1]
  if (length(alpha) == 1) alpha <- rep(alpha, r)
  if (length(dim) < 3 || local_structure == "none") {
    p <- t(rDirichlet(prod(dim[-1]), alpha, r))
  } else {

    args <- list(...)
    args$nlev <- dim[-1]
    args$method <- local_structure
    tmp <- do.call(rand_partition, args)
    stopifnot(length(tmp$parts) ==  q)

    p <- t(rDirichlet(length(tmp$partition), alpha, r))[, tmp$parts]
    attr(p, "parts") <- tmp$parts
  }

  dim(p) <- dim
  dimnames <- lapply(dim-1, seq.int, from = 0)
  names(dimnames) <- scope
  dimnames(p) <- dimnames

  return(p)
}

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
#' set.seed(1)
#' nlev <- c(2, 2, 2)
#' rand_partition(nlev, .75, method = "tree")
#'
#' # tree of given depth
#' rand_partition(nlev, 1, "tree", maxdepth = 1)
#' rand_partition(nlev, 1, "tree", maxdepth = .5)
#'
#' # tree of given depth, forced to be regular
#' rand_partition(nlev, 1, "tree", maxdepth = 1, regular = TRUE)
rand_partition <- function(nlev,
                           prob,
                           method = "tree",
                           regular = FALSE,
                           ...) {
  method <- match.arg(method, c("tree"))
  if (length(nlev) < 2) {
    stop("Need levels of at least two variable to produce a partition.")
  } else {

    P <- switch(method,
                "labels" = partition_from_labels(rand_labels(nlev, prob, ...), nlev),
                "tree" = rand_partition_tree(nlev, splitprob = prob, regular, ...))
  }
  indx <- try(order(unlist(P)))
  if (inherits(indx, "try-error")) {
    stop()
  }
  list(partition = P,
       parts = rep.int(seq_along(P), lengths(P))[order(unlist(P))])
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


#' grow a tree that partition the space of `length(nlev)` variables
#' @param splitprob (numeric constant) probability of splitting a node
#' @param doMerge (logical constant) randomly merge leaves in tree
#' @param nextsplitprob (function) manipulate `splitprob` for each new split
rand_partition_tree <- function(nlev,
                                splitprob,
                                regular,
                                nextsplitprob = function(x) x,
                                mindepth = 1,
                                maxdepth = length(nlev)) {

  n <- length(nlev)
  if (mindepth%/%1 == 0) mindepth <- max(1, round(n*mindepth))
  if (maxdepth%/%1 == 0 || maxdepth == 1) maxdepth <- min(round(n*maxdepth), n)

  # define routine for growing a tree
  grow_tree <- function(splitprob, vars, subset) {

    if ( (n-length(vars) >= mindepth && n-length(vars) >= maxdepth) || runif(1) > splitprob) {
      return(list(subset = subset))
    }

    # sample a split variable
    x <- sample(vars, 1)

    # split the current subset by values of x
    xval <- (subset%/%stride[x])%%nlev[x]
    new_subsets <- unname(split(subset, xval))

    # grow a new tree for each value of the split variable
    list(var = x,
         branches = lapply(new_subsets,
                           function(y) grow_tree(nextsplitprob(splitprob),
                                                 vars[!vars == x],
                                                 y)))
  }



  # define routine for extracting subsets in each leaf
  unlist_tree <- function(tree, name) {
    if (is.null(tree$branches)) unname(tree[name])
    else unlist(lapply(tree$branches, unlist_tree), recursive = FALSE)
  }

  n <- length(nlev)
  stride <- c(1, cumprod(nlev[-length(nlev)]))
  subset <- seq_len(prod(nlev))-1
  vars <- seq_along(nlev)
  tree   <- grow_tree(splitprob, vars, subset)
  P <- unlist_tree(tree)

  if (regular) {
    get_splitvars <- function(tree) {
      if (!is.null(tree$branch)) {
        c(tree$var, unlist(lapply(tree$branches, get_splitvars), recursive = FALSE))
      }
    }
    vars_in <- unique(get_splitvars(tree))
    for (x in setdiff(vars, vars_in)) {
      subset <- P[[1]]
      xval <- (subset%/%stride[[x]])%%nlev[[x]]
      P <- c(P[-1], unname(split(subset, xval)))
    }
  }

  stopifnot(!is.null(P))
  stopifnot(length(unlist(P)) == prod(nlev))
  return(P)
}

