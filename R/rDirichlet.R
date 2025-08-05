#' Wrapper around rBeta2009::rdirichlet that always return an n-by-k matrix
rDirichlet <- function(n, alpha, k = length(alpha)){
  if (k == 2){
    p <- matrix(numeric(), n, 2)
    p[, 1] <- rBeta2009::rbeta(n, alpha[1], alpha[2])
    p[, 2] <- 1-p[, 1]
  } else {
    p <- rBeta2009::rdirichlet(n, alpha)
  }
  return(p)
}
