
dag_support <- function(dags, p = rep(1/length(dags), length(dags))) {
  setNames(frequency_table(lapply(dags, as.matrix), p), c("dags", "p"))
}

