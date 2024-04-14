




rand_cpt_arrays <- function(dag, nlev = 2, alpha = 1) {
  d <- ncol(dag)
  if (length(nlev)==1){
    nlev <- rep(nlev, d)
  }
  oc   <- lapply(nlev-1, seq, from = 0)
  cpts <- vector("list", d)
  names(cpts) <- names(oc) <- names(nlev) <- colnames(dag)
  for (i in seq_len(d)){
    pa <- which(dag[,i]==1)

    # normalize
    if (length(pa) == 0){
      u <- array(rgamma(prod(nlev[i]), alpha), nlev[i], oc[i])
      cpts[[i]] <- u/sum(u)
    } else {
      npaconf <- prod(nlev[pa])
      u <- array(rgamma(nlev[i]*npaconf, alpha), nlev[c(i, pa)], oc[c(i, pa)])
      cpts[[i]] <- u/rep.int(colSums(u), rep.int(nlev[i], npaconf))
    }
  }
  return(cpts)
}


cpt_arrays_from_mat <- function(cpts){

  stopifnot(!is.null(names(cpts)))
  oc  <- lapply(cpts, function(x) colnames(x$cpd))
  noc <- vapply(oc, length, integer(1))
  out <- setNames(vector("list", length(cpts)), names(cpts))

  for (x in names(cpts)){
    pa <- colnames(cpts[[x]]$par_conf)
    if (is.null(pa)) {
      out[[x]] <- array(cpts[[x]]$cpd, noc[x], oc[x])
    } else {
      out[[x]] <- array(t(cpts[[x]]$cpd), noc[c(x, pa)], oc[c(x, pa)])
    }
  }
  return(out)
}




