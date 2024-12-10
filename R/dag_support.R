
#' Compute support of DAG posterior
#'
#' @param dags list of adjacencey matrices
#' @param p vector of weights. Uniform by default.
#'
#' @return a list representing with the unique DAGs in `dags` and their relative
#'  freuency
#' @export
dag_support <- function(dags, p = rep(1/length(dags), length(dags))) {
  setNames(frequency_table(lapply(dags, as.matrix), p), c("dags", "p"))
}

