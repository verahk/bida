#' Optimize partition over parent outcome space
#'
#' @param counts (integer matrix)
#'  a frequency table.
#' @param levels (list of integer vectors)
#'  levels of each conditioning variable in `counts`,
#'  such that `expand.grid(levels)` gives the joint configurations corresponding
#'  to each row in `counts`.
#' @param ess (numeric)
#'  imaginary sample size for the bdeu-prior.
#' @param method (character)
#' @param regular (logical)
#'  if `FALSE` (default) the optimized partition is returned, also if it implies
#'  conditional independencies (is not "regular"). If `TRUE` each part of the
#'  optimized partition is divided by the outcomes of the relevant variables.
#'  See [make_regular()].
#' @param verbose (logical)
#' @return a list with named elements:
#' - `partition`: the partition implied by the tree
#' - `scores`: each part's contribution to the local score
#' -  additional output from the different optimization procedures.n
#' @export
#' @examples
#'
#'
#' # non-binary-split-CSI
#' levels <- list(0:1, 0:1)
#' counts <- cbind(c(100, 10, 10, 10), c(10, 100, 100, 100))
#'
#' tree = optimize_partition(counts, levels, 1, "tree")
#' ldag = optimize_partition(counts, levels, 1, "ldag")
#' part = optimize_partition(counts, levels, 1, "part")
#'
#' cbind(expand.grid(levels),
#'       counts,
#'       get_parts(tree$partition),
#'       get_parts(ldag$partition),
#'       get_parts(part$partition))
#'
#' # mixed cardinality
#' levels <- list(0:1, 0:2, 0:3)
#' r <- 3
#' q <- prod(lengths(levels))
#' counts <- cbind(1, 10, 10*rep.int(0:5, q/5))
#'
#' tree = optimize_partition(counts, levels, 1, "tree")
#' ldag = optimize_partition(counts, levels, 1, "ldag")
#'
#' cbind(expand.grid(levels),
#'       counts,
#'       get_parts(tree$partition),
#'       get_parts(ldag$partition))
#'
#' # lower stopping criteria
#' tree = optimize_partition(counts, levels, 1, "tree", min_score_improv = -10)
#' ldag = optimize_partition(counts, levels, 1, "ldag", min_score_improv = -10)
#'
#' cbind(expand.grid(levels),
#'       counts,
#'       get_parts(tree$partition),
#'       get_parts(ldag$partition))
#'
#' # non-regular partitions
#' levels <- list(0:1, 0:1)
#' counts <- cbind(c(10, 10, 10, 10), c(100, 100, 100, 100))
#'
#' tree = optimize_partition(counts, levels, 1, "tree", regular = F)
#' ldag = optimize_partition(counts, levels, 1, "ldag", regular = F)
#' part = optimize_partition(counts, levels, 1, "part", regular = F)
#'
#' cbind(expand.grid(levels),
#'       counts,
#'       get_parts(tree$partition),
#'       get_parts(ldag$partition),
#'       get_parts(part$partition))
#'
#' # force regular structure
#' tree = optimize_partition(counts, levels, 1, "tree", regular = T)
#' ldag = optimize_partition(counts, levels, 1, "ldag", regular = T)
#' part = optimize_partition(counts, levels, 1, "part", regular = T)
#'
#' cbind(expand.grid(levels),
#'       counts,
#'       get_parts(tree$partition),
#'       get_parts(ldag$partition),
#'       get_parts(part$partition))
#'
#'
optimize_partition <- function(counts, levels, ess, method, regular = F, min_score_improv = 0, verbose = FALSE){
  method <- match.arg(method, c("tree", "ldag", "part"))
  res <- switch(method,
               "tree" = optimize_partition_tree(counts, levels, ess, min_score_improv, verbose = verbose),
               "ldag" = optimize_partition_ldag(counts, levels, ess, regular, min_score_improv, verbose = verbose),
               "part" = optimize_partition_part(counts, levels, ess, regular, min_score_improv, verbose = verbose))

  if (regular && method == "tree") {
    # ensure that partition is regular
    new_P <- make_regular(res$partition, lengths(levels))
    if (!length(new_P) == length(res$partition)) {
      new_counts <- rowsum(counts, get_parts(new_P))
      new_scores <- famscore_bdeu_byrow(new_counts, ess,
                                        r = ncol(counts), q = nrow(counts), s = lengths(new_P))
      res$partition <- new_P
      res$counts <- new_counts
      res$scores <- new_scores
    }
  }
  return(res)
}

#' @rdname optimize_partition
optimize_partition_from_cpt <- function(cpt, ess, method, regular = T, return_cpt = TRUE) {
  n <- length(cpt$levels)
  if (n > 2) {
    res <- optimize_partition(cpt$counts, cpt$levels[-n], ess, method, regular, verbose = FALSE)
    if (return_cpt) {
      cpt$partition <- res$partition
      cpt$counts <- res$counts
      cpt$score  <- sum(res$scores)
      return(cpt)
    } else {
      res$partition
    }
  } else if (return_cpt) {
    return(cpt)
  }
}
