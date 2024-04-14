
#' Compute BDeu family score for categorical data.
#'
#' Computes BDeu family score for all possible parent-child combinations. The scores
#' are written to file in GOBNILP format.
#'
#' @param data an n-by-d data frame or data matrix.
#' @param max_parent_size limit on the maximum parent set size.
#' @param file_out name of the file to which the scores are written (".score"
#           is added to the file name).
#' @param ess equivalent sample size parameter.
#'
#' @return a file with the parent scores in GOBNILP format.
#'

famscore_to_file_bdeu <- function(data, max_parent_size=ncol(data)-1, file_out, ess = 1, nlev){

  if (is.data.frame(data)) data <- as.matrix(data)

  d <- ncol(data)

  if (max_parent_size == Inf | max_parent_size >= d){
    max_parent_size <- d-1
  }
  nps <- sum(choose(d-1, seq(0, max_parent_size)))
  if (nps > 1e6){
    stop("The number of parent sets is > 10^6. Comment this out if you still want to proceed.")
  }


  # create .score file
  fid <- file(file_out,"wt")
  on.exit(close(fid))
  writeLines(toString(d), con = fid, sep = "\n")

  # for each node, compute score of all possible parent_sets and write to file
  for (x in 1:d) {
    writeLines(paste(x, nps), con = fid, sep = "\n")
    for (size in 0:max_parent_size) {
      if (size == 0){
        parent_sets <- list(NULL)
      } else {
        parent_sets <- utils::combn(seq_len(d)[-x], size, simplify = FALSE)
      }
      for (pa in parent_sets) {
        counts <- counts_from_data_matrix(data[, c(x, pa), drop = FALSE], nlev[c(x, pa)])
        scr <- famscore_bdeu(matrix(counts, ncol = nlev[x], byrow = TRUE), ess)

        writeLines(paste(trimws(format(round(scr, 6), nsmall=6)),
                         paste(c(size, pa), collapse = " "),
                         sep = " "),
                   con = fid, sep = "\n")
      }
    }
  }
}







