

#' Compute backdoor parameters
#'
#' Compute parameters for direct backdoor estimation of an intervention distribution.
#'
#' @inheritParams bida_pair
#' @param z (integer vector) column position of the adjustment variables.
#'   If `any(y == z)`, the marginal distribution of `y` is returned.
#'   If `length(z) = 0`, the conditional distribution of `y` given `x`.
#' @return
#' - If `type` is in `c("cat", "ldag", "tree")`, an object of class [bida_bdeu].
#' @export
#' @examples
#'
#' # Categorical data ----
#' nlev <- 2:4
#' lev  <- lapply(nlev-1, seq.int, from = 0)
#' data <- as.matrix(expand.grid(lev))
#' y <- 3
#' x <- 2
#' z <- 1
#'
#' type <- "tree"
#' hyperpar <- list(nlev = nlev, ess = 1)
#'
#' # empty adjustment set
#' backdoor_params(type, data, x, y, integer(0), hyperpar, NULL)
#'
#' # no effect, y is in adjustment set
#' backdoor_params(type, data, x, y, c(y, z), hyperpar, NULL)
#'
#' # return params for reduced CPT
#' bdeu_tree <- backdoor_params(type, data,  x, y, z, hyperpar, NULL)
#'
#' # from lookup table
#' ## compute score with user-defined score-function in BiDAG::scoreparameter
#' lookup   <- rlang:::new_environment()
#' scorepar <- define_scoreparameters(data, "tree", c(hyperpar, list(edgepf = 1)), lookup = lookup)
#' BiDAG:::usrDAGcorescore(y, c(x, z), n = 3, scorepar)
#' ls.str(lookup)
#' bdeu_lookup <- backdoor_params(type, data,  x, y, z, hyperpar, lookup)
#' stopifnot(all.equal(bdeu_tree$counts, bdeu_lookup$counts))
#' all.equal(bdeu_tree$partition, bdeu_lookup$partition)      # ordering differ
#' table(get_parts(bdeu_tree$partition), get_parts(bdeu_lookup$partition))
#'
backdoor_params <- function(type, data, x, y, z, hyperpar, lookup) {
  if (type %in% c("cat", "ldag", "tree")) {
    if (length(z) > 0 && any(y == z)) {
      bida_bdeu(data, y, integer(0), hyperpar$ess, hyperpar$nlev)
    } else if (length(z) == 0 || type == "cat") {
      bida_bdeu(data, y, c(x, z), hyperpar$ess, hyperpar$nlev)
    } else {

      if (!is.null(lookup[[type]])) {
        parentnodes <- sort.int(c(x, z))
        parID <- paste(c(y, parentnodes$x), collapse = ".")
        if (parID %in% names(lookup[[type]])) {
          bdeu <- lookup[[type]][[parID]]
          if (parentnodes[1] == x) {
            return(bdeu)
          } else {
            # permute dimension of bdeu object
            # the intervention variable x is assumed to be second dim
            pos_x <- match(x, parentnodes, 0L)
            perm <- c(0, pos_x, seq_along(parentnodes)[-pos_x]) +1
            return(aperm.bida_bdeu(bdeu, perm))
          }
        }
      }

      bdeu <- bida_bdeu(data, y, c(x, z), hyperpar$ess, hyperpar$nlev, NULL)
      opt  <- optimize_bdeu(bdeu,
                            method = type,
                            levels = hyperpar$levels[c(x, z)],
                            regular = hyperpar$regular)
      bdeu$partition <- opt$partition
      return(bdeu)
    }
  } else {
    stop(paste0("type, ", type, " is not supported."))
  }
}

