

#' Fast cross tabulation of categorical data
#'
#' Faster version of `table(as.data.frame(data))` for numeric matrices with
#' integer values.
#'
#' @param data (integer matrix)
#'  each column represent observed values of a categorical variable.
#' @param nlev (integer vector)
#'  number of categories of each variable.
#' @param sparse (boolean)
#'  if TRUE, a [bida_sparse_array] with counts is returned.
#'  If false, the counts are returned in a [base::array] with dimensions given by `nlev`.
#'  Default to false.
#' @details
#'  The ith variable (column i of `data`) is assumed to take values from `0` to `nlev[i]-1`.
#'  This is not checked.
#' @return An array with counts over each joint combination in `data`
#' @keywords internal
counts_from_data_matrix <- function (data, nlev, sparse = F) {

  n <- ncol(data)
  cump <- cumprod(nlev)
  group  <- c(data%*%c(1, cump[-n]))

  if (sparse) {

    # compute counts over observed groups
    ugroup <- unique(group)
    tmp <- tabulate(match(group, ugroup), length(ugroup))

    # return sparse array with counts
    new_bida_sparse_array(tmp, index = ugroup, dims = nlev)
  } else {
    tmp <- tabulate(group + 1, cump[n])
    array(tmp, nlev)                    # return array with counts
  }
}










