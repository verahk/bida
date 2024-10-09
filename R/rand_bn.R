
#' Draw a random Bayesian network
#'
#' Draw a random DAG and a distribution over the DAG, and store the BN as a `bnlearn::bn.fit` object.
#'
#' @keywords internal
#' @param dag (integer matrix) adjacency matrix.
#' @param type (character) type of distribution.
#' @param ...
#'
#' @return a `bnlearn::bn.fit` object
#' @keywords internal
#'
#' @examples
#'
#' n <- 3
#' dag <- matrix(0, n, n)
#' dag[upper.tri(dag)] <- 1
#'
#'
#' # categorical ---
#'
#' # draw random CPT with local structure
#' nlev <- c(2, 2, 2)
#' cpt <- rand_cpt(nlev, c("y", "x", "z"), local_structure = "tree", prob = 1, maxdepth = .5)
#'
#' # draw random DAG
#' set.seed(007)
#' dag <- rand_dag(3, 2)
#' dag
#'
#' # draw random distrib over DAG
#' cpts <- rand_dist(
#'   dag,
#'   "cat",
#'   nlev = nlev,
#'   local_structure = "tree",
#'   prob = 1,
#'   maxdepth = .5
#' )
#'
#' # create bn object
#' bn <- custom_bn(dag, cpts)
#'
#' #
#' set.seed(007)
#' bn2 <- rand_bn(3, 2, "cat", nlev = nlev, local_structure = "tree", prob = 1, maxdepth = .5)
#' stopifnot(all.equal(bn, bn2))
#'

rand_bn <- function(n, d, type = "cat", attr_to_keep = character(0), ...) {
  dag <- rand_dag(n, d)
  dist <- rand_dist(dag, type, ...)

  # store as bn-objects. Skip tests, assumed correct from rand_dag() and rand_dist()
  bn <- custom_bn(dag, dist, run_checks = FALSE)

  for (a in attr_to_keep) {
    attr(bn, a) <- lapply(dist, attr, which = a)
  }

  return(bn)
}

custom_bn <- function(dag, dist, run_checks = TRUE) {

  # init bn-learn object
  varnames <- colnames(dag)
  if (is.null(varnames)) varnames <- paste0("X", seq_len(ncol(dag)))
  g <- bnlearn::empty.graph(varnames)
  bnlearn::amat(g) <- dag

  # create bn.fit object
  if (run_checks) {
    bnlearn::custom.fit(g, dist)
  } else {
    # skip time-consuming checks of whether or not the dag is valid and so on
    bnlearn:::custom.fit.backend(x = x, dist = dist, ordinal = character(0), debug = FALSE)
  }

}


rand_dist <- function(dag, type = "cat", ...) {
  type <- match.arg(type, c("categorical"))

  varnames <- colnames(dag)
  if (is.null(varnames)) varnames <- paste0("X", seq_len(ncol(dag)))


  draw_dist <- function(node, args) {
    nodes <- c(node,  which(dag[, node] == 1))
    if (type == "categorical") {
      args$nlev <- args$nlev[nodes]
      args$scope <- varnames[nodes]
      do.call(rand_cpt, args)
    }
  }
  args <- list(...)
  setNames(lapply(seq_along(varnames), draw_dist, args = args), varnames)
}


rand_cpt <- function(nlev, scope = NULL, alpha = 1, local_structure = "none", ...) {
  stopifnot(length(nlev) > 0 && class(nlev) == "numeric")
  r <- nlev[1]
  q <- prod(nlev[-1])


  dimnames <- lapply(nlev-1, seq.int, from = 0)
  names(dimnames) <- scope
  p <- array(NA, nlev, dimnames)
  if (length(nlev) < 3 || local_structure == "none") {
    p[] <- vapply(integer(q),
                  function(x) rDirichlet(1, rep(alpha, r), r), numeric(r))
  } else {

    args <- list(...)
    args$nlev <- nlev[-1]
    args$method <- local_structure
    tmp <- do.call(rand_partition, args)
    stopifnot(length(tmp$parts) ==  q)

    p[] <- vapply(seq_along(tmp$partition),
                  function(x) rDirichlet(1, rep(alpha, r), r), numeric(r))[, tmp$parts]

    attr(p, "parts") <- tmp$parts
  }

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
