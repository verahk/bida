  # Computes Gaussian BIDA posterior for some given parent support.
# ------------------------------------------------------------------------------
# INPUT:  - data: an n-by-d data matrix
#         - ps: list containing parent supports for each node.
#         - th: threshold parameter; approximate the parent posterior by the
#               largest scores that together execeed the threshold.
# OUTPUT: - bida_post: bida object containing the parameters that specify the
#                      BIDA posterior for all cause-effect pairs.
# ------------------------------------------------------------------------------
bida_posterior_gauss <- function(ps, data, th=0.999){
  n <- nrow(data)
  d <- ncol(data)
  S <- t(data)%*%data
  a0 <- 1
  b0 <- 1
  c_mu0 <- 0
  diag_V0 <- 1
  for (node in 1:d){
    tmp <- sort(x=ps$support[[node]],
                index.return=TRUE,
                decreasing=TRUE)
    ind <- 1:which(cumsum(tmp$x)>=th)[1]
    ps$support[[node]] <- tmp$x[ind]
    ps$support[[node]] <- ps$support[[node]]/sum(ps$support[[node]])
    ps$sets[[node]] <- ps$sets[[node]][tmp$ix[ind],,drop=FALSE]
  }
  parameters <- matrix(list(),d,d)
  parents <- vector("list",d)

  for (x in 1:d){
    par_mat <- ps$sets[[x]]
    temp <- lapply(1:d, function(j) if(j != x){cbind(ps$support[[x]],matrix(0,nrow(par_mat),3))})
    for (i in 1:nrow(par_mat)){
      z <- par_mat[i,!is.na(par_mat[i,])]
      k <- length(z)
      mu0 <- matrix(c_mu0,k+1,1)
      A0 <- diag(k+1)/diag_V0
      xzTxz <- S[c(x,z),c(x,z)]
      A1 <- A0+xzTxz
      V1 <- solve(A1)
      a1 <- a0+n/2
      for (y in setdiff(1:d,x)){
        if (y %in% z) {
          temp[[y]][i,2:4] <- rep(NA,3)
        } else {
          xzTy <- S[c(x,z),y]
          mu1 <- V1%*%(A0%*%mu0+xzTy)
          b1 <- as.vector(b0+(t(mu0)%*%A0%*%mu0+S[y,y]-t(mu1)%*%A1%*%mu1)/2)
          v <- 2*a1
          V2 <- (b1/a1)*V1
          temp[[y]][i,2:4] <- c(mu1[1],V2[1,1],v)
        }
      }
    }
    for (y in setdiff(1:d,x)){
      para_mat <- temp[[y]]
      colnames(para_mat) <- c("w","mu","sigma","ny")
      parameters[[x,y]] <- para_mat
    }
    parents[[x]] <- par_mat
  }
  bp <- list("parameters"=parameters, "sets"=parents, "type"="Gaussian")
  class(bp) <- "bida"
  return(bp)
}
