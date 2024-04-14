

avg_abs_array <- function(p, k = dim(p), digits = 15){

  if (length(k) == 2) dim(p) <- c(dim(p), 1)
  indx <- combn(1:k[1], 2)
  diff <- abs(p[indx[1, ], , , drop = FALSE]-p[indx[2, ], , ,drop = FALSE])
  round(1/2*1/ncol(indx)*colSums(diff,  dims = 2), digits)
}


avg_akl_array <- function(p, k = dim(p), digits = 15){

  if (length(k) == 2) dim(p) <- c(dim(p), 1)
  logp <- log(p)
  x  <- apply(p, 3, function(pp) tcrossprod(pp, log(pp)))
  indx <- which(diag(k[1])==1)
  round(1/k[1]*colSums(x[indx,])-1/k[1]**2*colSums(x), digits)
}


avg_jsd_array  <- function(p, k = dim(p), digits = 15){

  M <- 1/k[1]*colSums(p)
  logM <- log(M)
  logM[M == 0] <- 0

  logp <- log(p)
  logp[p == 0] <- 0

  if (length(k) > 2) {
    round(1/k[1]*colSums(p*logp, dims = 2)-colSums(M*logM), digits)
  } else {
    round(1/k[1]*sum(p*logp)-sum(M*logM), digits)
  }
}



abs_dev <- function(p, k = dim(p), M) {
  if (length(k) == 2) {
    mean(abs(p-rep(M, each = k[1])))
  } else {
    colSums(abs(p-rep(M, each = k[1])), dims = 2)/prod(k[1:2])
  }
}
