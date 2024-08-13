


#' Compute Bayesian score of DAG
#'
#' @param data (matrix) data matrix
#' @param j (integer) column position of node
#' @param parentnodes (integer vector) column position of parents of node `j`
#' @param type (character) name of score
#' @param params (list) list with additional parameters
#' @return (numeric constant) the score of node `j` given `parentnodes`
#' @param ess (integer) equivalent sample size
#' @param nlev (integer vector) cardinality of each variable
#' @param partitions (list) list with partitioning of the parent space of node `j`
#' @export
#'
#' @examples
#'
#' n   <- 3
#' dag <- matrix(0, n, n)
#' dag[upper.tri(dag)] <- 1
#'
#' # categorical
#' N <- 100
#' params <- list(nlev = rep(3, n), ess = 1)
#' data <- sapply(params$nlev, sample, size = N, replace = T) -1
#' score <- score_dag(dag, data, type = "cat", params)
#' score == score_dag(t(dag), data, type = "cat", params) # score equiv
#'
#' # labeled DAG - add list of partitions of each CPT
#' params$partition <- list(NULL, NULL, list(c(0, 1), 2, 3))
#' score_dag(dag, data, type = "cat", params)
score_dag <- function(dag, data, type = "cat", params) {
  n <- ncol(dag)
  scores <- numeric(n)

  for (j in seq_len(n)) {
    pa <- which(dag[, j] == 1)
    scores[j] <- score_fam(data, j, pa, type, params)
  }
  sum(scores)
}

score_fam <- function(data, j, parentnodes, type, params) {
  switch(type,
         "cat" = score_fam_cat(data, j, parentnodes, params$ess, params$nlev, params$partitions))
}

score_fam_cat <- function(data, j, parentnodes, ess, nlev, partitions = NULL) {
  r <- nlev[j]
  npar <- length(parentnodes)
  if (npar == 0) {
    tab <- tabulate(data[, j]+1, r)
    famscore_bdeu_1row(tab, ess)
  } else if (is.null(partitions) || is.null(partitions[[j]])) {

    # enumerate joint outcomes
    stride <- c(1, cumprod(nlev[parentnodes]))
    joint  <- data[, c(parentnodes, j), drop = FALSE]%*%stride

    # compute frequency table
    q   <- stride[npar+1]
    tab <- matrix(tabulate(joint + 1, q*r), q, r, byrow = F)

    # compute family score
    famscore_bdeu(tab, ess, r, q)

  } else {

    q <- prod(nlev[parentnodes])
    stopifnot(length(unlist(partitions[[j]])) == q)

    # enumerate joint parent outcomes
    stride <- c(1, cumprod(nlev[parentnodes[-npar]]))
    pa     <- data[, parentnodes, drop = FALSE]%*%stride

    # map to parts
    parts <- get_parts(partitions[[j]])[pa+1]

    # enumerate joint outcomes and compute freq table
    nparts <- length(partitions[[j]])
    joint  <- parts + nparts*data[, j]
    tab <- matrix(tabulate(joint, nparts*r), nparts, r)

    # compute family score
    sum(famscore_bdeu_byrow(tab, ess, r, q, lengths(partitions[[j]])))
  }
}

