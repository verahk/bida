#' Check if partition is regular.
#'
#' A partition is called regular if it does not encode any conditional independencies.
#' The tests if there is any variable for which all parts in the partition is
#' independent of that variable's value.
#' More precisely, it iterates through the variables and parts and checks if all
#' values of the variable is repeated for every configuration of the remainder
#' variables in that part.
#'
#' @param P (list of integer vectors) a partitioning of a zero-indexed outcome space
#'  `0, 1, ..., prod(nlev)-1, see examples.
#' @param nlev (integer vector) cardinality of the variables
#' @param stride (integer vector) stride corresponding to `nlev`. Optional.
#' @seealso [make_regular()]
#' @return a logical constant; `TRUE` if partition is regular, `FALSE` otherwise.
#' @examples
#'
#' # specify cardinality of variables
#' nlev   <- c(2, 2)
#'
#' ## regular
#' P <- list(0, 1:3)
#' stopifnot(is_regular(P, nlev) == TRUE)
#'
#' P <- list(c(0, 3), c(1, 2))
#' stopifnot(is_regular(P, nlev) == TRUE)
#'
#' ## not regular
#' P <- list(0:3)
#' stopifnot(is_regular(P, nlev) == FALSE)
#'
#' P <- list(0:1, 2:3)
#' stopifnot(is_regular(P, nlev) == FALSE)

is_regular <- function(P, nlev, stride = NULL) {
  if (any(lengths(P) < min(nlev))) return(TRUE)
  if (is.null(stride)) stride <- c(1, cumprod(nlev[-length(nlev)]))

  is_regular <- T
  i <- 0
  while (is_regular && i < length(nlev)) {
    i <- i+1
    is_regular <- is_regular_varwise(P, nlev, stride, i = i)
  }
  return(is_regular)
}

#'@rdname is_regular
#'@param i position of node
#'@details
#' - `is_regular_varwise`: returns `TRUE` if `P` is regular with respect to variable `i`
is_regular_varwise <- function(P, nlev, stride, i = 1) {
  # check if all values of i is repeated for every configuration of the other variables
  if (any(lengths(P)%%nlev[i] > 0)) return(TRUE)
  for (p in P) {
    x <- (p%/%stride[i])%%nlev[i]    # values of variable i
    if (all(x[1] == x) || any(table(p-stride[i]*x, x) == 0)) return(TRUE)
  }
  return(FALSE)
}
