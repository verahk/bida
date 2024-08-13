


#' Make a CPT partition regular
#' 
#' If the partition implies any conditional independencies, return a finer partitioning
#' that does not, by splitting all parts by the outcomes of the relevant variable(s).
#'
#' @param P (list) a partition of the parent space
#' @param nlev (integer vector) the cardinality of the variables
#'
#' @return a regular partition 
#' @export
#' @seealso [is_regular()]
#' 
#' @examples
#' nlev <- c(2, 2)
#' 
#' P <- list(0:3)
#' make_regular(P, nlev)
#' 
#' P <- list(c(0, 2), c(1, 3))
#' make_regular(P, nlev)
#' 
#' P <- list(0, 1:3)
#' make_regular(P, nlev)
make_regular <- function(P, nlev) {
  if ((min(lengths(P)) >= min(nlev))) {
    stride <- c(1, cumprod(nlev[-length(nlev)]))
    new_P <- P 
    for (i in seq_along(nlev)) {
      if (!is_regular_varwise(P, nlev[i], stride[i])) {
        tmp   <- lapply(new_P, function(p) split(p, (p%/%stride[i])%%nlev[i]))
        new_P <- unlist(tmp, recursive = F, use.names = FALSE)
      }
    }
    return(new_P)
  } else {
    return(P)
  }
}
