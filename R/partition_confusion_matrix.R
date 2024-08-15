
#  Compare partitions
#  Compute paired-confusion-matrix https://stats.stackexchange.com/questions/548778/how-to-compute-a-pair-confusion-matrix.
compute_confusion_matrix_for_partitions <- function(x, y) {
  stopifnot(length(x) == length(y))
  res <- matrix(character(), length(x), length(y))
  # loop over all pairs
  for (i in seq_along(x)[-1]) {
    for (j in seq.int(1, i-1)) {
      if (x[i] == x[j]) {
        if (y[i] == y[j]) {
          res[i, j] <- "a"
        } else {
          res[i, j] <- "b"
        }
      } else {
        if (y[i] == y[j]) {
          res[i, j] <- "c"
        } else {
          res[i, j] <- "d"
        }
      }
    }
  }
  table(factor(res, c("a", "b", "c", "d")))
}
