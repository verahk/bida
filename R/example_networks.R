


#' Title
#'
#' @param name (character) the name of an example network
#'
#' @return a bn.fit object
#' @export
#'
#' @examples
#' bn <- example_networks("illustrate_posteiror")
#' g <- as
example_networks <- function(name) {
  if (name == "illustrate_posterior") {

    # Network designed for illustrating the multi-modality of the posterior
    # distributions over causal effects

    dag <- cbind(x1 = c(0, 0, 0, 0),
                 x2 = c(1, 0, 0, 0),
                 x3 = c(1, 1, 0, 1),
                 x4 = c(0, 0, 0, 0))
    rownames(dag) <- colnames(dag)

    n <- ncol(dag)
    nlev <- setNames(rep(2, n), colnames(dag))
    lev  <- lapply(nlev-1, seq.int, from = 0)

    cpts <- list()
    cpts[[1]] <- array(c(.5, .5), nlev[1], lev[1]) # adjust distribution of root node

    pa <- 1
    cpts[[2]] <- array(c(.1, .9, .8, .2), nlev[c(2, pa)], lev[c(2, pa)])


    pa <- which(dag[, 3] == 1)
    tmp <- c(c(.1, .9, .7, .3),
             c(.05, .95, .95, .05),
             c(.15, .85,  .85, .15),
             c(.3, .7, .9, .1))
    cpts[[3]] <- array(tmp, nlev[c(3, pa)], lev[c(3, pa)])

    cpts[[4]] <- array(c(.5, .5), nlev[4], lev[4])
    names(cpts) <- colnames(dag)
    return(bn_from_cpt_arrays(cpts))
  }
}
