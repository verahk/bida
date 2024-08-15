
# WHAT: compute causal parameters of a discrete network
# WHY:  ground truth in BIDA-simulations
# HOW:


#' @param bn bn.fit.dnet object
#' @param ace_funs list of functions that evaluate the causal effect given a rx-by-ry interventional CPT
#' @return a list with causal parameter of the network `bn`
compute_ground_truth <- function(bn,
                                 ace_funs = list(jsd = bida:::avg_jsd_array,
                                                 abs = bida:::avg_abs_array)) {
  n    <- length(bn)
  nlev <- vapply(bn, function(x) dim(x$prob)[1], integer(1))
  dag  <- bnlearn::amat(bn)
  cpdag <- bnlearn::amat(bnlearn::cpdag(bn))

  dmat <- bida:::descendants(dag)
  desc <- which(dmat == 1 & diag(n) == 0)

  pdmat <- bida:::descendants(cpdag)
  pdesc <- which(pdmat == 1 & diag(n) == 0)

  pdo <- bida:::interv_probs_from_bn(bn, method = "bn")
  ace <- vapply(ace_funs, function(f) vapply(pdo[desc], f, numeric(1)), numeric(length(desc)))

  list(n = n,
       nlev = nlev,
       dag = dag,
       cpdag = bnlearn::amat(bnlearn::cpdag(bn)),
       dmat = dmat,
       desc = desc,
       pdmat = pdmat,
       pdesc = pdesc,
       pdo = pdo,
       ace = ace,
       ace_funs = ace_funs)
}
