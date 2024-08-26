

#' Title
#'
#' @param splitprob
#' @param nlev
#' @return
#' @export
#'
#' @examples
#'
#' set.seed(007)
#' nlev <- c(2, 2, 2)
#' rand_partition(nlev, prob = .5)
#'
rand_partition <- function(nlev, prob, method = "labels") {
  P <- switch(method,
             "labels" = partition_from_labels(rand_labels(nlev, prob), nlev),
             "tree" = rand_partition_tree(nlev, prob, doMerge = F, ...),
             "dgraph" = rand_partition_tree(nlev, prob, doMerge = T, ...))

  # ensure parts are enumerated from 1,..., length(unique(P))
  match(P, unique(P))
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
      # draw number of labels
      #p <- lprob/2**seq.int(q/nlev[i])
      #p <- p/sum(p)
      #nlabels <- sample.int(length(p), 1, prob = p)
      nlabels <- sample.int(length(p), 1)

      # draw labels
      contexts <- sample(joint[conf[, i] == 0], nlabels, FALSE)
      labels[[i]] <- conf[contexts, -i, drop = FALSE]
    }
  }
  labels
}



#' Title
#'
#' @param nlev
#' @param splitprob
#' @param doMerge
#' @param nextsplitprob
#' @return
#' @keywords internal
#' @examples
#'
#' set.seed(007)
#' nlev <- c(2, 2, 2)
#' rand_partition_tree(nlev, splitprob = .5, doMerge = F)
#'
#' # force at least one split by manipulating `nextsplitprob` rule
#' set.seed(007)
#' nextsplitprob <- function(x) x/2
#' part_tree <- rand_partition_tree(nlev, splitprob = 1, doMerge = F, nextsplitprob)
#' part_tree
#'
#' # set doMerge = T to merge leaves in decision tree
#' set.seed(007)
#' part_ldag <- rand_partition_tree(nlev, splitprob = 1, doMerge = T, nextsplitprob)
#'
#' # compare
#' levels <- lapply(nlev-1, seq.int, from = 0)
#' cbind(expand.grid(levels), part_tree, part_ldag)
#'
rand_partition_tree <- function(nlev, splitprob, doMerge = T, nextsplitprob = function(x) x) {

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

  if (doMerge && nparts > 1) {
    nparts <- length(partition)
    tmp <- sample.int(nparts, nparts, TRUE)  # assign parts to new parts
    partition <- lapply(split(partition, tmp), do.call, what = "c")
  }

  get_parts(partition)
}




