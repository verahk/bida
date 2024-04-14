

#' Exact computation of conditional probability distributions
#'
#' Compute conditional probability distributions from a set of CPTs associated
#' associated with a Bayesian network, using the variable elimination algorithm.
#'
#' @param cpts (list)
#'  list of [base::array]s that specifies a conditional probability table $P(X_i| Pa(X_i))$
#'  For each array the `dimnames` attribute must be a named list, indicating the
#'  scope of the vector - i.e. the names of X_i and Pa(X_i).
#' @param y (integer vector)
#'  position of the (set of) variable(s) in `cpts`
#' @param x (integer vector)
#'  position of variables that constitutes the conditioning set.
#'  Defaults to the `integer(0)`
#' @param anc (logical vector)
#'  indicator of ancestors of `y` and `x`.
#'  If `NULL`, the ancestors is computed via the DAG implied by `cpts`.
#'
#' @return an array with the conditional distribution of `y` given `x`.
#'
#' @keywords internal
cpquery_from_cpt_arrays <- function(cpts, y, x = integer(0), anc = NULL) {

  if (is.null(anc)) {
    dag  <- dag_from_cpt_arrays(cpts)
    dmat <- descendants(dag)
    anc  <- rowSums(dmat[, c(y, x), drop = FALSE]) > 0
  }

  # compute marg joint prob
  varnames <- names(cpts)
  keep <- c(y, x)
  stopifnot(all(anc[keep]))
  pyx  <- sum_product_ve(cpts[anc], varnames[replace(anc, keep, FALSE)])

  if (length(keep) > 1) {

    # permute array
    pyx  <- aperm(pyx, varnames[keep])

    if (!is.null(x)) {
      # compute cond prob p(y|x)
      ny <- length(y)
      ky <- prod(dim(pyx)[seq_len(ny)])
      px <- colSums(pyx, ny)
      pyx <- pyx/rep(px, each = ky)
    }
  }
  return(pyx)
}


dag_from_cpt_arrays <- function(cpts) {
  n   <- length(cpts)
  scopes <- lapply(cpts, scope)
  varnames <- names(cpts)

  dag <- matrix(0, n, n)
  colnames(dag) <- rownames(dag) <- varnames

  for (i in seq_len(n)[lengths(scopes) > 1]) {
    pa <- scopes[[i]][-1]
    dag[pa, i] <- 1
  }
  return(dag)
}

stride <- function(x) {
  nlev <- dim(x)
  c(1, cumprod(nlev[-length(nlev)]))
}
stride.bn.fit <- function(x) {
  stride.default(x$prob)
}

scope <- function(x) names(dimnames(x))
scope.bn.fit.dnode <- function(x) c(x$node, x$parents)

sum_out <- function(x, margin) {

  nlev <- dim(x)
  perm <- c(margin, seq_along(nlev)[-margin])
  tmp <- colSums(aperm(x, perm), dims = length(margin))

  # return array with dim and dimnames
  array(tmp, nlev[-margin], dimnames(x)[-margin])
}

sum_product_ve <- function(factors, vars_to_elim){
  if (length(vars_to_elim) == 0) {
    if (length(factors) == 1) {
      factors[[1]]
    } else {
      factors_product_c(factors)
    }
  } else {
    for (z in vars_to_elim){

      factors <- sum_product_elim_var(factors, z)

      # factors_ <- sum_product_elim_var(factors, z)
      # for (f in factors_) stopifnot(!is.null(dim(f)))
      # factors <- factors_
    }
    factors_product_c(factors)
  }
}



sum_product_elim_var <- function(factors, z){
  stopifnot(is.list(factors))
  z_in_scope <- vapply(factors, function(f) z %in% scope(f), logical(1))

  if (!any(z_in_scope)) {
    return(factors)
  } else {

    # compute factor product of factors with z in scope
    # f <- Reduce(factor_product_c, factors[z_in_scope])
    f <- factors_product_c(factors[z_in_scope])
    scope_f <- scope(f)
    nlev_f  <- dim(f)

    if (length(nlev_f) == 1) {
      stopifnot(scope_f == z)
      sum_f <- sum(f)

      if (all(z_in_scope)) {
        return(sum(f))
      } else if (sum_f == 1) {
        return(factors[!z_in_scope])
      } else {
        return(lapply(factors[!z_in_scope], function(x) x*sum_f))
      }

    } else {

      # sum out z
      new_factor <- sum_out(f, match(z, scope_f))

      # return list with factors
      return(c(factors[!z_in_scope], list(new_factor)))
    }
  }
}


factor_product <- function(x, y) {

  if (length(x) == 1) {
    if (length(y) == 1) {
        return(x[1]*y[1])
    } else {
      return(y*x[1])
    }
  } else if (length(y) == 1){
    return(x*y[1])
  }

  scope_x <- scope(x)
  scope_y <- scope(y)

  # find which variables in y that are not in x
  yinx    <- match(scope_y, scope_x)
  nomatch <- is.na(yinx)

  if (all(nomatch)) {
    return(outer(x, y))
  } else {

    scope <- c(scope_x, scope_y[nomatch])
    nlev  <- c(dim(x), dim(y)[nomatch])
    value <- numeric(prod(nlev))

    stride_x <- c(stride(x), rep(0, sum(nomatch)))
    stride_y <- replace(nlev*0, match(scope_y, scope), stride(y))

    seqn <- seq_along(scope)
    ass  <- seqn*0
    j <- 1
    k <- 1

    for (i in seq_along(value)){

      value[i] = x[j]*y[k];

      for (v in seqn){
        if (ass[v] == nlev[v]-1){
          ass[v] = 0;
          j = j - (nlev[v]-1)*stride_x[v];
          k = k - (nlev[v]-1)*stride_y[v];
        } else {
          ass[v] = ass[v]+1;
          j = j + stride_x[v];
          k = k + stride_y[v];
          break;
        }
      }
    }
  }

  dim(value) <- nlev
  dimnames(value) <- setNames(vector("list", length(nlev)), scope)

  return(value)
}


