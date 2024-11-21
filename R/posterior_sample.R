
#' @noRd
#' @export
posterior_sample <- function(x, n, ...) {
  UseMethod("posterior_sample")
}

#' @noRd
#' @export
posterior_sample.NULL <- function(x, n, ...) {
  NULL
}
