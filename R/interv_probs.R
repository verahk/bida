

#' Compute intervention probabilities
#'
#' @name interv_probs
#' @param cpts a list of cpt_arrays of length `n`, where `names(dimnames(cpts[[i]]))` gives
#'  scope of the variable `i`.
#' @param bn a object of class `bnlearn::bn.fit.dnet` of length `n`
#' @param method either `"exact"` or `"bn"`
#' @return a `n`-by-`n` matrix with intervention probabilities
#' @details
#' These functions computes interventional CPTs for every cause-effect pair.
#'
#' - If `method == "exact"` the exact intervention probabilities is computed.
#'   The `interv_probs_from_bn` function then first extract the cpt_arrays from
#'   the bn-network.
#' - If `method == "bn"`, the parameters are estimated using [bnlearn::cpdist].
#'   The `interv_probs_from_cpt_arrays` then first construct a `bn.fit.dnet` object.
#' `interv_probs` loops over every node, identifies the true descendants of
#'  the node and compute the intervention distributions. For non-descendants,
#'  the interventional CPTs are filled with the marginal distribution of the
#'  pre-intervention distribution.
NULL

#' @rdname interv_probs
#' @param obj a object defining the distribution, see `cpts` or `bn`
#' @param dag adjacency matrix of `obj`
#' @export
interv_probs <- function(obj, dmat, method = c("exact", "bn")) {

  # precompute marginal probs
  margs <- switch(method,
                  "exact" = marginal_probs_exact(obj, dmat),
                  "bn" = marginal_probs_bn(obj))

  # prep
  nlev  <- vapply(margs, length, integer(1))
  n <- length(nlev)
  seqn <- seq_len(n)


  # compute intervention probs for all nodes
  pdo  <- matrix(list(), n, n)
  for (i in seqn) {

    k <- nlev[i]
    pdo[[i, i]] <- diag(k)

    # for non-descendants
    nonDe <- which(dmat[i, ] == 0)
    if (length(nonDe) > 0) {
      pdo[i, nonDe] <- lapply(margs[nonDe],
                              function(p) matrix(rep(p, each = k), nrow = k))
    }

    # for descendants
    if (length(nonDe) < n-1) {
      De <- seqn[-c(i, nonDe)]
      pdo[i, De] <- switch(method,
                           "exact" = interv_probs_x_exact(obj, i, De, dmat),
                           "bn" = interv_probs_x_bn(obj, i, De, dmat))
    }
  }
  return(pdo)
}

#' @rdname interv_probs
#' @export
interv_probs_from_cpt_arrays <- function(cpts, method =  c("exact", "bn")) {
  dag  <- dag_from_cpt_arrays(cpts)
  dmat <- bida::descendants(dag)
  if (method == "bn") {
    bn <- bn_from_cpt_arrays(cpts, dag)
    interv_probs(bn, dmat, method)
  } else {
    interv_probs(cpts, dmat, method)
  }
}


#' @rdname interv_probs
#' @export
interv_probs_from_bn <- function(bn, method = c("exact", "bn")) {
  if (method == "exact") {
    cpts <- cpt_arrays_from_bn(bn)
    interv_probs(cpts, descendants(bn), method)
  } else {
    interv_probs(bn, descendants(bn), method)
  }
}


#' Compute exact intervention probabilities under interventions on `x`
#' @rdname interv_probs
#' @inheritParams interv_probs
#' @return a list of `length(y)` with interventional cpts
interv_probs_x_exact <- function(cpts, x, y, dmat) {

  varnames <- names(cpts)
  k <- dim(cpts[[x]])[1]
  seqk <- seq_len(k)

  out <- vector("list", length(y))
  # compute marginals of intervention distirbution
  for (j in seq_along(y)) {

    yy   <- y[j]
    anc  <- dmat[, yy] == 1
    elim <- replace(anc, c(x, yy), FALSE)

    if (!any(elim)) {

      # x is only ancestor of yy
      out[[j]] <- t(unname(cpts[[yy]]))

    } else {

      #  for each intervention level
      evidence <- array(0, k, setNames(list(NULL), varnames[x]))
      probs <- vector("list", k)
      for (kk in seqk) {

        # replace cpt of intervention variable
        do_cpts <- cpts
        do_cpts[[x]] <- replace(evidence, kk, 1)

        # compute marginal probs
        tmp <- sum_product_ve(do_cpts[anc], varnames[elim])
        tmp <- aperm(tmp, varnames[c(x, yy)])
        probs[[kk]] <- tmp[kk, ]
      }

      out[[j]] <- do.call("rbind", probs)
    }
  }
  return(out)
}

#' Compute approximate intervention probabilities under interventions on `x`,
#' using [bnlearn::mutilated] and [bnlearn::cpdist].
#' @inheritParams interv_probs
#' @return a list of `length(y)` with interventional cpts
interv_probs_x_bn <- function(bn, x, y, dmat) {

  varnames <- names(bn)
  evidence <- as.list(dimnames(bn[[x]]$prob)[[1]])
  k   <- length(evidence)
  names(evidence) <- rep(varnames[x], k)

  # for each intervention level,
  # sample from the multilated network and compute marginal probs
  probs <- array(list(), c(k, length(y)))
  for (kk in seq_len(k)){
    bn_do  <- bnlearn::mutilated(bn, evidence[kk])
    data   <- bnlearn::cpdist(bn_do, varnames[y], evidence = TRUE)
    probs[kk, ]  <- lapply(data, function(x) tabulate(x, nlevels(x))/nrow(data))
  }

  apply(probs, 2, function(x) do.call("rbind", x), simplify = FALSE)
}


interv_prob_from_cpts_mc <- function(cpts, x, y, oc, top_ordering = NULL, samplesize = 10**4){

  if (is.null(top_ordering)) top_ordering <- top_order(dag_from_cpt_arryays(cpts))


  # enumerate nodes preceding all nodes in y, i.e. nodes that need to be sampled (don't sample descendants of y not in y)
  indx <- seq_len(max(match(y, top_ordering)))

  # init
  evidence <- list(pos = x, values = NULL)
  out  <- lapply(y, function(j) array(dim = noc[c(x, j)], dimnames = oc[c(x, j)]))
  for (k in seq_len(noc[x])){
    # sample data from cpts, given values of x, and compute freq for each variable in y
    evidence$values <- rep.int(k-1, samplesize)
    data <- sample_data_from_cpts(cpts, samplesize, top_ordering[indx], list(evidence)) +1
    for (j in seq_along(y)) {
      out[[j]][k, ] <- tabulate(data[, y[j]], noc[y[j]])/samplesize
    }
  }
  return(out)
}

# compute marginal probabilities
marginal_probs_exact <- function(cpts, dmat) {
  lapply(seq_along(cpts),
         function(j) cpquery_from_cpt_arrays(cpts, j, NULL, dmat[, j] == 1))
}
marginal_probs_bn  <- function(bn) {
  data  <- bnlearn::cpdist(bn, names(bn), evidence = TRUE)
  lapply(data, function(x) tabulate(x, nlevels(x))/nrow(data))
}



## OLD

interv_probs_from_bn_old <- function(bn, samplesize = 10**4){

  n <- length(bn)
  oc  <- lapply(bn, function(x) dimnames(x$prob)[[1]])
  noc <- lengths(oc)

  dag <- bnlearn::amat(bn)
  dmat <- bida:::descendants(dag)

  # pre-compute marginal probabilities
  data  <- bnlearn::rbn(bn, samplesize)
  margs <- lapply(data, function(x)  tabulate(x, nlevels(x))/samplesize)


  pdo <- matrix(list(), n, n)
  for (i in seq_len(n)){
    nonDesc <- which(dmat[i, ] == 0)
    Desc <- seq_len(n)[-c(i, nonDesc)]

    # interv probs of node i
    pdo[[i, i]] <- array(diag(noc[i]), noc[c(i, i)], oc[c(i, i)])

    # intervention probs of non-descendants are the marginal distrib
    if (length(nonDesc) > 0){
      pdo[i, nonDesc] <- lapply(nonDesc, function(j) array(rep(margs[[j]], each = noc[i]),
                                                           noc[c(i, j)],
                                                           oc[c(i, j)]))
    }

    # compute MC-estimate of intervention probs for descendants
    if (length(Desc) > 0){
      pdo[i, Desc] <- lapply(Desc, function(j) array(NA, noc[c(i, j)], oc[c(i, j)]))
      for (k in seq_len(noc[i])){
        bn_do  <- bnlearn::mutilated(bn, setNames(list(oc[[i]][k]), names(oc[i])))
        data   <- bnlearn::rbn(bn_do, samplesize)
        for (j in Desc){
          pdo[[i, j]][k, ] <- tabulate(data[, j], noc[j])/samplesize
        }
      }
    }
  }
  return(pdo)
}


