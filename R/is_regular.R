#' Check if partition is regular.
#'
#' A partition is called regular if it do not encode any conditional independencies.
#' That is, the partition is independent of one or more conditioning variable(s).
#' The function iterates through each conditioning variable until a variable and 
#' return `FALSE` such a variable is met.
#' 
#' @param P a partitioning of a outcome space. 
#'  See class documentation [partition()].
#' @param nlev 
#' @param stride 
#' @seealso [make_regular()]
#' @return a logical constant; `TRUE` if partition is regular, `FALSE` otherwise.
#' @export
#' @examples
#' 
#' # check if partition is regular 
#' nlev   <- c(2, 2)
#' 
#' ## regular 
#' P <- list(0, 1:3)
#' is_regular(P, nlev)
#' 
#' P <- list(c(0, 3), c(1, 2)) 
#' is_regular(P, nlev)
#' 
#' ## not regular
#' P <- list(0:3)
#' is_regular(P, nlev)
#' 
#' P <- list(0:1, 2:3)
#' is_regular(P, nlev)
is_regular <- function(P, nlev, stride = c(1, cumprod(nlev[-length(nlev)]))) {
  is_regular <- TRUE
  if (all(lengths(P) >= min(nlev))) {
    n <- length(nlev)
    i <- 0 
    while (is_regular && i < n) {
      i <- i+1
      is_regular <- is_regular_varwise(P, nlev, stride, i)
    }
  } 
  return(is_regular)
}

is_regular <- function(P, nlev, stride = c(1, cumprod(nlev[-length(nlev)]))) {
  is_regular <- TRUE
  if (all(lengths(P) >= min(nlev))) {
    i <- 0
    while (is_regular && i < length(nlev)) {
      i <- i+1
      is_regular <- is_regular_varwise(P, nlev[i], stride[i])
    }
  } 
  return(is_regular)
}

is_regular_varwise <- function(P, k, s) {
  is_regular <- TRUE
  j <- 0
  while (is_regular && j < length(P)) {
    j <- j+1
    is_regular <- part_depends(P[[j]], k, s)
  }
  return(is_regular)
}

part_depends <- function(p, k, s) {
  if (length(p)%%k > 0) {
    return(TRUE)
  } else {
    x <- (p%/%s)%%k  
    all(x[1] == x) || any(tabulate(x+1, k) != length(x)%/%k) || any(table(x, p-s*x) == 0)
  }
}
