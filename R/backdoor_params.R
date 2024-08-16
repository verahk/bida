

#' Compute parameters for direct backdoor estimation of an intervention distribution

#' @inheritParams bida_pair
#' @param z (integer vector) column position of the adjustment variables.
#'   If `any(y == z)`, the marginal distribution of `y` is returned.
#'   If `length(z) = 0`, the conditional distribution of `y` given `x`.
#' @return
#' - If `type` is in `c("cat", "ldag", "tree")`, an object of class [bida_bdeu].
#' @export
#'
backdoor_params <- function(type, data, x, y, z, hyperpar, lookup) {
  if (type %in% c("cat", "ldag", "tree")) {
    if (y == z) {
      bida_bdeu(data, y, integer(0), hyperpar$ess, hyperpar$nlev)
    } else if (type == "cat" || length(z) == 0) {
      bida_bdeu(data, y, c(x, z), hyperpar$ess, hyperpar$nlev)
    } else {
      if (is.null(lookup)) {
        parID <- paste(c(y, sort(c(x, z))), collapse = ".")
        if (parID %in% names_lookup) {
          return(lookup[[parId]]$bdeu)
        }
      }

      bdeu <- bida_bdeu(data, x, c(x, z), hyperpar$ess, hyperpar$nlev, NULL)
      opt  <- optimize_bdeu(bdeu,
                            method = type,
                            ess    = hyperpar$ess,
                            levels = hyperpar$levels,
                            regular = hyperpar$regular)
      bdeu$partition <- opt$partition
      return(bdeu)
    }
  } else {
    stop(paste0("type, ", type, " is not supported."))
  }
}
