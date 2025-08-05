#' Draw a conditional probability table (CPT)
#'
#' Construct a random conditional probability table (CPT) by sampling independently
#' a vector of probabilities from a Dirichlet distribution for each parent configuration.
#' If the scope of the CPT is specified, the function returns an array with
#' named dimensions, in the same format as [bnlearn::bn.fit.dnode] store CPTs.
#'
#' @param nlev (integer vector) of cardinality of the outcome variable and the
#'  parent variables (the dimensions of the CPT).
#' @param scope (character vector) names of the variables. Defaults to `names(nlev)`.
#' @param alpha (integer)
#' @param local_struct (character) name of algorithm to form a partitioning /
#'  produce parameter restrictions.
#'  Defaults to `"none"`, which returns a CPT with no parameter restrictions.
#' @param ... additional arguments sendt to [rand_partition].
#'
#' @return an array with dimension `nlev` and dimension named by `scope`.
#' @keywords internal
#'
#' @examples
#'
#' # no parents
#' cpt <- rand_cpt(nlev = 2, scope = "Y")
#' cpt
#' dimnames(cpt)
#'
#' # two parents
#' nlev <- c(Y = 2, X = 3, Z = 4)
#' cpt <- rand_cpt(nlev)
#' cpt
#'
#' # with local structure
#' \dontrun{
#' cpt <- rand_cpt(nlev, local_struct = "tree", prob = 1, maxdepth = .5)
#' }
rand_cpt <- function(nlev,
                     scope = names(nlev),
                     alpha = 1,
                     j = 1,
                     parentnodes = seq_along(nlev)[-1],
                     local_struct = "none",
                     ...) {

  nlev <- nlev[c(j, parentnodes)]
  r <- nlev[1]
  q <- prod(nlev[-1])
  if (length(alpha) == 1) alpha <- rep(alpha, r) else stopifnot(length(alpha) == r)
  if (length(nlev) < 3 || local_struct == "none") {
    p <- t(rDirichlet(q, alpha, r))
  } else {

    args <- list(...)
    args$nlev <- nlev[-1]
    args$method <- local_struct
    print(args)
    tmp <- do.call(rand_partition, args)
    print(tmp)
    stopifnot(length(tmp$parts) ==  q)

    p <- t(rDirichlet(length(tmp$partition), alpha, r))[, tmp$parts]
    attr(p, "parts") <- tmp$parts
  }

  dim(p) <- nlev
  dimnames <- lapply(nlev-1, seq.int, from = 0)
  names(dimnames) <- scope
  dimnames(p) <- dimnames

  return(p)
}
