

#' Compute the mean of an intervention distribution
#'
#' Compute the mean of an intervention distribution by applying the backdoor formula.
#'
#' @param x object of class [`bida_bdeu`]
#' @param ... additional arguments
#' @return the posterior mean(s) of the distribution(s) specified by `x`
#' @export
backdoor_mean <- function(x, ...){
  UseMethod("backdoor_mean", x)
}

backdoor_mean.bida_bdeu <- function(x) {
  dims <- get_dim(x)
  if (length(dims) < 2) {
    mean_bdeu(x)
  } else {

    ess <- x$ess

    # compute conditional means
    py.xz <- mean_bdeu(x, reduced = F)

    # compute counts over adjustment sets
    Nyxz <- as.array(x$counts)
    az <- colSums(Nyxz, dims = 2) + ess/prod(dims[-c(1:2)])

    # sum out z
    rowSums(py.xz*rep(az, each = prod(dims[1:2])), dims = 2)/sum(az)
  }
}
