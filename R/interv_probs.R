

#' Compute discrete intervention distributions
#'
#' For each variable in a discrete Bayesian network and every possible intervention
#' on that variable, compute the marginal intervention distribution of all other
#' variables in the network.
#' The default version is based on [bnlearn::mutilated] and [bnlearn::cpdist],
#' and estimates the interventional distributions by sampling from the network
#' under all possible interventions.
#' An exact version is via [SumProductVE::cpquery] is applicable to smaller networks.
#'
#' @param obj a object of class `bnlearn::bn.fit.dnet` of length `n` or
#'  a list of CPTs stored as arrays.
#' @param method either `"bnlearn"` (default) or `"exact"`.
#' @return a `n`-by-`n` matrix with intervention probabilities
#' @examples
#' bn <- readRDS("./data/alarm.rds")
#'
#' pdo <- interv_probs(bn)
#' pdo
#'
#' system.time(pdo2  <- interv_probs(bn))
#' system.time(exact <- interv_probs(bn, method = "exact"))
#'
#' plot(unlist(pdo), unlist(exact))
#' abline(a = 0, b = 1, col = "red")
#'
#' plot(unlist(pdo), unlist(pdo2))
#' abline(a = 0, b = 1, col = "red")
#'
#' dag <- matrix(0, 3, 3)
#' dag[upper.tri(dag)] <- 1
#' cpts <- rand_dist(dag, "cat", nlev = rep(3, 3))
#' pdo_bn <- interv_probs(cpts)
#' pdo_exact <- interv_probs(cpts, method = "exact")
#' plot(unlist(pdo_exact), unlist(pdo_bn))
#' abline(a = 0, b = 1, col = "red")
interv_probs <- function(obj, method = "bnlearn") {
  stopifnot(method == "bnlearn" || nchar(system.file(package = "SumProductVE")) > 0)
  if (inherits(obj, "bn.fit.dnet")) {
    dmat <- descendants.bn.fit(obj)
    if (method == "exact") {
      # extract CPTs from bn
      obj <- cpts_from_bn(obj)
    }
  } else if (all(vapply(obj, is.array, logical(1)))) {
    if (method == "bnlearn") {
      obj  <- custom_bn(dag_from_cpts(cpts), cpts)
      dmat <- descendants.bn.fit(obj)
    } else {
      if (!all(lengths(lapply(obj, dimnames)) > 0)) {
        stop("all CPTs must have named dimnames")
      }
      dmat <- descendants.matrix(dag_from_cpts(obj))
    }
  } else {
    stop("obj must be either an object of class bn.fit.dnet or a list of arrays")
  }

  n <- length(obj)
  seqn <- seq_len(n)

  # pre-compute marginal probs for zero-effets
  margs <- switch(method,
                  "exact" = marginal_probs_from_cpts(obj, seqn, dmat),
                  "bnlearn" = marginal_probs_from_bn(obj, seqn))
  nlev  <- vapply(margs, length, integer(1))


  # compute intervention probs for all nodes
  pdo  <- matrix(list(), n, n)
  for (i in seqn) {

    k <- nlev[i]

    # for intervention node, IPT equal diagonal matrix
    pdo[[i, i]] <- diag(k)

    # for non-descendants, IPTs equal marginal probs
    nonDe <- seqn[dmat[i, ] == 0]
    if (length(nonDe) > 0) {
      pdo[i, nonDe] <- lapply(margs[nonDe],
                              function(p) matrix(rep(p, k), ncol = k))
    }

    # for descendants, compute IPTs in mutilated network
    if (length(nonDe) < n-1) {
      De <- seqn[-c(i, nonDe)]
      pdo[i, De] <- switch(method,
                           "exact" = interv_probs_from_cpts(obj, i, De, dmat),
                           "bnlearn" = interv_probs_from_bn(obj, i, De))
    }
  }
  return(pdo)
}





#' Compute exact intervention probabilities
#'
#' Compute exact intervention probabilities under interventions on single cause variable `x`
#'
#' @param cpts a list of of length `n` with CPTs stored in arrays, where `names(dimnames(cpts[[i]]))` gives
#'  the scope of the `i`th CPT. Note that
#' @return a list of length `length(y)` with the interventional probability tables
#' @keywords internal
interv_probs_from_cpts <- function(cpts, x, y, dmat) {

  varnames <- names(cpts)
  k <- dim(cpts[[x]])[1]
  seqk <- seq_len(k)
  px.do   <- array(0, k, setNames(list(NULL), varnames[x]))

  out <- vector("list", length(y))
  for (j in seq_along(y)) {
    yy   <- y[j]
    anc  <- dmat[, yy] == 1
    if (anc[x]) {
      # indicator for variables to eliminate / marginalize out
      elim <- replace(anc, c(x, yy), FALSE)
      if (any(elim)) {
        p <- array(NA, dim = c(dim(cpts[[yy]])[1], k))
        for (kk in seqk) {
          # replace cpt of intervention variable
          cpts[[x]] <- replace(px.do, kk, 1)
          p[, kk] <- SumProductVE::cpquery(cpts, yy, anc = anc)
        }
        out[[j]] <- p
      } else {
        out[[j]] <- cpts[[yy]]
      }
    } else {
      p <- SumProductVE::cpquery(cpts, yy, anc = anc)
      out[[j]] <- array(rep(p, k), c(length(p), k))
    }
  }
  return(out)
}

#' Approximate intervention probabilities
#'
#' Approximate intervention probabilties under interventions on single`x`,
#' using [bnlearn::mutilated] and [bnlearn::cpdist].
#' @param bn an object of class bn.fit.dnet
#' @return a list of length `length(y)` with the interventional probability tables
#' @keywords internal
interv_probs_from_bn <- function(bn, x, y) {

  varnames <- names(bn)
  evidence <- as.list(dimnames(bn[[x]]$prob)[[1]])
  k   <- length(evidence)
  names(evidence) <- rep(varnames[x], k)

  # for each intervention level,
  # sample from the multilated network and compute marginal probs
  probs <- array(list(), c(k, length(y)))
  for (kk in seq_len(k)){
    bn_do  <- bnlearn::mutilated(bn, evidence[kk])
    df  <- bnlearn::cpdist(bn_do, varnames[y], evidence = TRUE)
    probs[kk, ]  <- lapply(df, function(x) tabulate(x, nlevels(x))/nrow(df))
  }

  apply(probs, 2, function(x) do.call("cbind", x), simplify = FALSE)
}

# compute marginal probabilities
marginal_probs_from_cpts <- function(cpts, y = seq_along(cpts), dmat) {
  lapply(y, function(j) SumProductVE::cpquery(cpts, j, NULL, dmat[, j] == 1))
}
marginal_probs_from_bn  <- function(bn, y = seq_along(bn)) {
  data  <- bnlearn::cpdist(bn, names(bn)[y], evidence = TRUE)
  lapply(data, function(v) tabulate(v, nlevels(v))/nrow(data))
}



cpts_from_bn <- function(bn) {
  lapply(bn, function(x) {
    if (length(x$parents) == 0) {
      array(x$prob, length(x$prob), setNames(vector("list", 1), x$node))
    } else {
      x$prob
    }
  })
}
dag_from_cpts <- function(cpts) {

  # list scope of each CPT
  scopes <- lapply(cpts, function(x) names(dimnames(x)))
  lens   <- lengths(scopes)
  stopifnot(all(lens>0))  # all CPTs do not have named dimnames attribute

  # init adjacency matrix
  n   <- length(cpts)
  dag <- matrix(0, n, n)
  colnames(dag) <- rownames(dag) <- sapply(scopes, "[[", 1)

  # set edges
  for (i in seq_len(n)[lens > 1]) {
    pa <- scopes[[i]][-1]
    dag[pa, i] <- 1
  }
  return(dag)
}

