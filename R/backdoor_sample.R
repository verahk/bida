

#' Compute the mean of an intervention distribution
#'
#' Compute the mean of an intervention distribution by applying the backdoor formula.
#'
#' @param x object of class [`bida_bdeu`]
#' @param ... additional arguments
#' @return the posterior mean(s) of the distribution(s) specified by `x`
#' @export
#'

backdoor_sample <- function(x, ...){
  UseMethod("backdoor_mean", x)
}


#' @rdname bdeu_par
#' @export
backdoor_sample.bida_bdeu <- function(size, obj, dims_ipt = get_dim(obj), digits = 10e-16) {

  dims <- get_dim(obj)
  ndims <- length(dims)
  if (ndims < 2) {
    tmp <- round(sample_bdeu(size, obj), digits)
    if (ndims == 2 || length(dims_ipt) == 1)  return(tmp)
    # produce a sample of CPTs by replicating the marginal distribution
    array(rep(tmp, each = dims_ipt[2]), c(dims_ipt, size))
  } else if (dims == 2) {
    round(sample_bdeu(size, obj), digits)
  } else {

    ess <- obj$ess
    k <- prod(dims)
    kyx <- prod(dims[1:2])
    kz  <- k/kyx

    # force 3-dim array
    obj$counts$dim <- c(dims[1:2], kz)

    # compute counts over adjustment sets
    # - rowsum_fast(v, group, ugroup) returns zero for elements of ugroup not found in group
    z  <- obj$counts$index%/%kyx
    az <- rowsum_fast(obj$counts$value, z, seq_len(kz)-1) + ess/kz

    # sample distributions over adjustment set
    pz <- t(round(rDirichlet(size, az, length(az)), digits))

    # sample conditional means
    py.xz <- round(sample_bdeu(size, obj, reduced = F), digits)

    colSums(rep(pz, each = kyx)*aperm(py.xz, c(3, 1, 2, 4)))
  }
}
