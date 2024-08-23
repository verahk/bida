

sample_data_from_cpts <- function(size, cpts, top_order) {
  n <- length(cpts)
  seq_n <- seq_len(n)
  nlev <- vapply(cpts, function(x) dim(x)[1], integer(1))

  data <- matrix(integer(), size, n)
  for (i in top_order) {
    pa <- seqn[dag[, i] == 1]
    npa <- length(pa)
    if (npa == 0) {
      data[, i] <- sample.int(nlev[i], size, replace = T, prob = cpts[[i]])
    } else if (npa == 1) {
      data[, i] <- sample_from_cpt(cpts[[i]], data[, pa], nlev[i])
    } else {
      stride <- c(1, cumprod(nlev[pa[-npa]]))
      pa_obs <- data[, pa]%*%stride
      data[, i] <- sample_from_cpt(cpts[[i]], data[, pa]%*%stride, nlev[i])
    }
  }
  return(data)
}

# r <- 3
# q <- 3
# cpt <- sapply(1:q, function(x) rDirichlet(1, rep(1, r), r))
#
# # sample parents uniformily
# parents <- sample.int(q, 10**4, TRUE)
#
# # sample from cpt
# x <- sample_from_cpt(cpt, parents-1, r)


sample_from_cpt <- function(cpt, parents, r = dim(cpt)[1]) {
  y <- integer(length(parents))
  for (k in unique(parents)) {
    indx <- parents == k
    y[indx] <- sample.int(r, sum(indx), replace = T, prob = cpt[, k+1])
  }
  return(y-1)
}
