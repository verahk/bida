

#' Class: bida
#'
#' An object of class bida is a matrix list with objects of class [bida::bida_pair],
#' that describes the posterior mixture distribution over the causal effect between every variable pair.
#' The `posterior_mean` and `posterior_sample` methods calls the respective methods
#' for each [bida::bida_pair] and collect the output in a matrix list.
#'
#' @param ps (list)
#' Support over adjustment sets, stored in a list with two named elements `sets` and `support`.
#' See [adjsets_support_from_dags] [parent_support_from_dags] [parent_support_exact]
#' @param data (matrix)
#' @param params additional params
#' @param pairs (matrix) a n-by-n matrix indicating which pairs should be considered.
#'  If NULL (default), all pairs are considered and `pairs[i, j] = 1` for all but `i == j`.
#' @param x an object of class [bida]
#' @param samplesize (integer) samplesize
#' @seealso [bida_pair_cat]
#' @return
#' `bida`: a n-by-n list where each element is an object of class `bida_pair`
#' @examples
#' dag  <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
#' colnames(dag) <- rownames(dag) <- letters[1:3]
#'
#' # Categorical data
#' nlev <- 2:4
#' bn <- rand_bn(dag, "cat", alpha = 1, nlev = nlev)
#' data <- sample_data_from_bn(bn, 10**4)
#'
#' # list parent sets in true DAG
#' ps <- parent_support_from_dags(list(dag))
#'
#' params <- list(nlev = nlev, ess = 1)
#' fit <- bida(ps, data, "cat", params)
#' means <- posterior_mean(fit)
#' means
#'
#' # compare with true distributions
#' pdo <- interv_probs_from_bn(bn, "exact")
#' dindx <- diag(ncol(dag)) == 1
#' plot(unlist(pdo[!dindx]), unlist(means[!dindx]))
#'
bida <- function(ps, data, type, params, lookup = NULL, pairs = NULL) {

  data <- as.matrix(data)
  n <- ncol(data)
  seqn <- seq_len(n)
  if (is.null(pairs)) pairs <- matrix(1, n, n)-diag(n)

  # store sets and support in matrix for for-loop below
  sets <- ps$sets
  supp <- ps$support
  if (length(ps$sets) == n) {
    sets <- matrix(sets, n, n)
    supp <- matrix(supp, n, n)
  }

  # construct a bida-pair object for each cause-effect pair
  out <- matrix(list(), n, n)
  for (x in seqn) {
    for (y in seqn) {
      if (pairs[x, y] == 1) {
        out[[x, y]] <- bida_pair(type, data, x, y, sets[[x, y]], supp[[x, y]], params, lookup)
      }
    }
  }

  colnames(out) <- rownames(out) <- colnames(data)

  bp <- structure(out,
                  params = params,
                  type   = type,
                  class  = c("bida"))
  return(bp)
}


#' @export
print.bida <- function(x) {
  cat("\nbida-object of type", attr(x, "type"))
  cat("\nNumber of variables: ", dim(x)[1])
}


#' @rdname bida
#' @export
posterior_mean.bida <- function(x, ...) {
  out  <- array(list(), dim(x), dimnames(x))
  indx <- sapply(x, is.null)
  out[!indx] <- lapply(x[!indx], function(v) posterior_mean(v, ...))
  return(out)
}

#' @rdname bida
#' @export
posterior_sample.bida <- function(x, n, ...) {
  out  <- array(list(), dim(x), dimnames(x))
  indx <- sapply(x, is.null)
  out[!indx] <- lapply(x[!indx], function(v) posterior_sample(v, n, ...))
  return(out)
}
