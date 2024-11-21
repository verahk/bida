

#' Sample categorical data
#'
#' @param n (integer) sample size
#' @param cpts (list of arrays) a list of array-CPTs, stored in the same format as
#'  in [bnlearn::bn.fit.dnode], i.e. with each column corresponding to a parent config
#'  and with named dimensions
#' @param top_order
#'
#' @return
#' @export
#' @details
#' Note that the index of the categorical variables are 0-based in [sample_data_from_cpts()]
#' and 1-based in [sample_data_from_cpts()].
#'
#' @examples
#'
#' dag  <- matrix(0, 3, 3, dimnames = list(c("Z", "X", "Y"), c("Z", "X", "Y")))
#' dag[upper.tri(dag)] <- 1
#' cpts <- rand_dist(dag, "cat", nlev = rep(3, ncol(dag)))
#'
#' data <- sample_data_from_cpts(10**4, cpts)
#'
#'
#' pz <- cpts[["Z"]]
#' pz
#' tabulate(data[, 1] +1, nlev[1])/nrow(data)
#'
#' pxz <- sweep(cpts[["X"]], "Z", pz, "*")
#' rowSums(pxz)
#' tabulate(data[, 2] +1, nlev[2])/nrow(data)
#'
#' pyxz <- sweep(cpts[["Y"]], c("X", "Z"), pxz, "*")
#' rowSums(pyxz)
#' tabulate(data[, 3] +1, nlev[3])/nrow(data)
#'
sample_data_from_cpts <- function(n, cpts, top_order = NULL) {
  if (is.null(top_order)) top_order <- top_order(dag_from_cpts(cpts))

  # collect cardinality of each variable
  nlev <- vapply(cpts, function(x) dim(x)[1], integer(1))

  data <- matrix(integer(), n, length(cpts)) # init
  for (i in top_order) {
    pa <- which(dag[, i] == 1)
    npa <- length(pa)
    if (npa == 0) {
      data[, i] <- sample.int(nlev[i], n, replace = T, prob = cpts[[i]])-1
    } else {
      if (npa == 1) {
        pa_obs <- data[, pa]
      } else {
        pa_obs <- data[, pa]%*%c(1, cumprod(nlev[pa[-npa]]))
      }
      data[, i] <- sample_data_from_cpt(cpts[[i]], pa_obs + 1)-1
    }
  }
  return(data)
}

top_order <- function(dag){
  n <- ncol(dag)
  marked <- vector("logical",n)
  order <- vector("numeric",n)
  for (i in 1:n){
    for (j in 1:n){
      found <- FALSE
      if (marked[j] == FALSE && all(marked[dag[,j] == 1]) == TRUE){
        marked[j] <- TRUE
        order[i] <- j
        found <- TRUE
        break
      }
    }
    if (found == FALSE){
      stop("No topological order exists!")
    }
  }
  return(order)
}

#' Draw a categorical variable given its parents
#'
#' @param cpt (numeric array) a CPT
#' @param parents (integer vector) parent configuration of each sample.
#' @return an integer vector of length `NROW(parents)` sampled from the CPT
#'
#' @keywords internal
#' @examples
#'
#' nlev <- c(3, 3, 3)
#' r <- nlev[1]        # cardinality of outcome variable
#' q <- prod(nlev[-1]) # cardinality of parent variables
#' cpt <- t(rDirichlet(q, rep(1, r)))
#' dim(cpt) <- nlev
#'
#' # sample parents uniformily
#' parents <- sample.int(q, 10**4, TRUE)
#'
#' # sample from cpt
#' y <- sample_from_cpt(cpt, parents)
#'
#' # compare means
#' tabulate(y + 1, r)/length(y)
#' rowMeans(cpt)
sample_data_from_cpt <- function(cpt, parents) {
  dims <- dim(cpt)
  r <- dims[1]                    # cardinality of outcome
  y <- integer(length(parents))   # init data vector
  if (length(dims) == 1) {
    y <- sample.int(r, NROW(parents), replace = T, prob = cpt)
  } else {
    q <- prod(dims[-1])
    dim(cpt) <- c(r, q)
    for (k in unique(parents)) {
      indx <- parents == k
      y[indx] <- sample.int(r, sum(indx), replace = T, prob = cpt[, k])
    }
  }
  return(y)
}

