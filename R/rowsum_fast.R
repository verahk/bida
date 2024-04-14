


#' Compute sum conditional on a grouping variable
#'
#' Faster version of [base::rowsum], to be used when the unique levels of the grouping variable is known.
#' Note that the sanity-checks of rowsum is excluded.
#'
#' @inheritParams base::rowsum
#' @param ugroup unique groups
#' @keywords internal
rowsum_fast <- function(x, group, ugroup, na.rm = FALSE) {
  if (is.null(dim(x))) {
    c(.Internal(rowsum_matrix(x, group, ugroup, na.rm, character(length(ugroup)))))
  } else {
    .Internal(rowsum_matrix(x, group, ugroup, na.rm, character(length(ugroup))))
  }
}
