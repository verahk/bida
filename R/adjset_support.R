
#' Compute support of adjustment set
#'
#' Identify the set of unique adjustment sets in a sample of DAGs.
#'
#' @param dag_supp (list) posterior support over DAGs.
#' @param x,y (integer vector) column position(s) of cause and effect variable(s), respectively
#' @param adjset (character) name of adjustment set. See [adjset::find_adjset()].
#' @param dmats (optional) list of matrices representing the ancestral relationships
#'  in each dag in `pdags[[1]]`. If not specified, [descendants()] is called for each
#'  dag.
#'
#' @examples
#'
#' data(bida_example_cat)
#' dags <- partitionMCMC$traceadd$incidence[-c(1:200)]
#'
#' x <- match("X", colnames(dag))
#' y <- match("Y", colnames(dag))
#' adjset_support(dag_support(dags), x, y)
#'
#' \dontrun{
#' library(adjset)
#' adjset_support(dag_support(dags), x, y, "o")
#' }
adjset_support <- function(dag_supp, x, y = 0L, adjset = "pa", dmats = NULL, replace_large_adjset = NULL) {

  dags <- dag_supp[[1]]
  p    <- dag_supp[[2]]
  n    <- ncol(dags[[1]])

  if (adjset == "pa") {
    sets <- lapply(dags, function(dag) seq_len(n)[dag[, x] == 1])
    zeroprob <- NULL
  } else {
    stopifnot(length(y) == 1 && y > 0)
    if (is.null(dmats)) dmats <- lapply(dags, descendants)
    yIsDescendant <- vapply(dmats, function(m) m[[x+n*(y-1)]], numeric(1)) == 1
    if (any(yIsDescendant)) {
      # list adjustment sets in DAGs phere y is descendant of x
      tmp  <- mapply(function(dag, dmat) adjset::find_adjset(dag, x-1, y-1, adjset, dmat)+1,
                     dag = dags[yIsDescendant],
                     dmat = dmats[yIsDescendant],
                     SIMPLIFY = FALSE)
      if (all(yIsDescendant)) {
        sets <- tmp
        zeroprob <- 0
      } else {
        sets <- c(list(y), tmp)
        p    <- c(sum(p[!yIsDescendant]), p[yIsDescendant])
        zeroprob <- p[1]
      }
    } else {
      # set (unique) adjustment set equal to y, indicating no effect from x to y
      sets <- list(y)
      p    <- 1
      zeroprob <- p[1]
    }
  }

  # aggregate probability over unique sets
  ps  <- frequency_table(sets, p)

  # replace large sets
  if (!is.null(replace_large_adjset)) {
    sets  <- lapply(ps[[1]], function(z) replace_large_adjset(x, y, z, adjset))
    ps <- frequency_table(sets, ps[[2]])
  }

  out <- setNames(ps, c("sets", "p"))
  out$zeroprob <- zeroprob
  attr(out, "adjset") <- adjset
  return(out)
}

