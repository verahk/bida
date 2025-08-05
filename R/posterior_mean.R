
#' @noRd
#' @export
posterior_mean <- function(x, ...) {
  UseMethod("posterior_mean")
}

#' @noRd
#' @export
posterior_mean.NULL <- function(x, ...) {
  NULL
}

