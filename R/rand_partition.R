

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
#' rand_partition(nlev, splitprob = .5, doMerge = F)
#'
#' # force at least one split by manipulating `nextsplitprob` rule
#' set.seed(007)
#' nextsplitprob <- function(x) x/2
#' rand_partition(nlev, splitprob = 1, doMerge = F, nextsplitprob)
#'
#' set.seed(007)
#' tree <- rand_partition(nlev, splitprob = 1, doMerge = F, nextsplitprob)
#'
#' set.seed(007)
#' ldag <- rand_partition(nlev, splitprob = 1, doMerge = T, nextsplitprob)
#'
#' # compare
#' levels <- lapply(nlev-1, seq.int, from = 0)
#' cbind(expand.grid(levels), tree, ldag)
#'
rand_partition <- function(nlev, splitprob, doMerge = T, nextsplitprob = function(x) x) {

  stride   <- c(1, cumprod(nlev[-length(nlev)]))
  outcomes <- parts <- seq_len(prod(nlev))-1

  grow_tree <- function(splitprob, vars, outcomes) {
    if (length(vars) == 0 || runif(1) > splitprob) {
      return(list(outcomes = outcomes))
    }

    # sample a split variable
    pos <- sample.int(length(vars), 1)
    x <- vars[pos]
    new_vars <- vars[-pos]

    # split the members in the current part by values of x
    new_partition <- unname(split(outcomes, (outcomes%/%stride[x])%%nlev[x]))

    # grow a new tree for each value of the split variable
    list(var = x,
         branches = lapply(new_partition,
                           function(y) grow_tree(nextsplitprob(splitprob), new_vars, y)))
  }

  unlist_tree <- function(tree) {
    if (is.null(tree$outcomes)) {
      unlist(lapply(tree$branches, unlist_tree), recursive = FALSE)
    } else {
      unname(tree["outcomes"])
    }
  }

  tree <- grow_tree(splitprob, vars, outcomes)
  partition <- unlist_tree(tree)
  nparts <- length(partition)
  if (doMerge && nparts > 1) {
    tmp <- sample.int(nparts, nparts, TRUE)  # assign parts to new parts
    partition <- lapply(split(partition, tmp), do.call, what = "c")
  }

  get_parts(c(partition, as.list(outcomes[-(unlist(partition)+1)])))
}



