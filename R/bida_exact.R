
#' Main function for computing the BIDA posterior using dynamic programming.
#'
#' @param data an n-by-d data matrix
#' @param type data type ("Categorical" or "Gaussian")
#' @param max_parent_size max size of parent sets (default: 4)
#' @param temp_id name for file in which scores are temporarily saved
#' @param filepath_apssolver filepath directing to the apssolver
#' @param ess equivalent sample for BDeu prior (default: 1). Ignored if [type != "Categorical"]
#' @param dirapprox indicator for whether parameters of the approximate dirichlet posterior of the intervention probabilities
#'  should be computed. If FALSE, the posterior parameters of the dirichlet distribution of the joint probability p(x, y, z) is returned.
#'
#' @return bida object containing the parameters that specify the
#                      BIDA posterior for all cause-effect pairs.
#'
#' @export
#'


bida_exact <- function(data, type="", max_parent_size=min(5, ncol(data)-1), temp_id = "", filepath_apssolver = "./aps-0.9.1/aps", ess = 1, nlev, dirapprox = FALSE){

  if (tolower(type) == "categorical") {
    ps <- parent_support_exact(data, type, max_parent_size, temp_id, filepath_apssolver, th = .999, nlev = nlev, ess = ess)
    bida_post <- bida_posterior_cat(ps, data, nlev, ess = ess, dirapprox = dirapprox)
  } else if (tolower(type) == "gaussian") {
    data <- apply(data, 2, function(x) x-mean(x))
    ps <- parent_support_exact(data, type, max_parent_size, temp_id, filepath_apssolver, ess)
    bida_post <- bida_posterior_gauss(ps, data, th=0.999)
  } else {
    stop("Data type must be either Cateogorical or Gaussian!")
  }
  return(bida_post)
}


