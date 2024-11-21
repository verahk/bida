#' Compute Bayesian score of a DAG
#'
#' @param data (matrix) data matrix
#' @param j (integer) column position of node
#' @param parentnodes (integer vector) column position of parents of node `j`
#' @param type (character) name of score
#' @param params (list) list with additional parameters
#' @return (numeric constant) the score of node `j` given `parentnodes`
#' @param ess (integer) equivalent sample size
#' @param nlev (integer vector) cardinality of each variable
#' @export
#'
#' @examples
#'
#' n   <- 3
#' dag <- matrix(0, n, n)
#' dag[upper.tri(dag)] <- 1
#'
#' # Categorical data ---
#' params <- list(nlev = 2:4, ess = 1)
#'
#' ## draw random data set
#' N <- 100
#' data <- sapply(params$nlev, sample, size = N, replace = T) -1
#'
#' ## test score equivalence
#' score1 <- score_dag(data, dag, type = "cat", params)
#' score2 <- score_dag(data, t(dag), type = "cat", params)
#' stopifnot(score1 == score2)
#'
score_dag <- function(data, dag, type = "cat", params) {
  n <- ncol(dag)
  scores <- numeric(n)

  for (j in seq_len(n)) {
    parentnodes <- which(dag[, j] == 1)
    scores[j] <- score_fam(data, j, parentnodes, type, params)
  }
  sum(scores)
}

#' @rdname score_dag
#' @export
score_fam <- function(data, j, parentnodes, type = c("cat"), params) {
  type = match.arg(type, c("categorical"))
  switch(type,
         "categorical" = score_fam_cat(data, j, parentnodes, params$ess, params$nlev))
}

#' @rdname score_dag
#' @export
score_fam_cat <- function(data, j, parentnodes, ess, nlev) {
  r <- nlev[j]
  q <- prod(nlev[parentnodes])
  subset <- c(parentnodes, j)
  tab <- counts_from_data_matrix(data[, subset, drop = FALSE], nlev[subset])
  dim(tab) <- c(q, r)
  return(famscore_bdeu(tab, ess, r, q))
}
