
compute_avgppv <- function(x, y, method = "noise") {
  switch(method,
         "noise" = compute_avgppv_noise(x, y),
         "dups"  = compute_avgppv_dups(x, y))
}
compute_avgppv_dups <- function(x, y) {

  indx <- order(x, decreasing = TRUE)
  tp   <- cumsum(y[indx])
  pp   <- seq_along(x)

  dups <- duplicated(x[indx], fromLast = TRUE)
  #cbind(dups, x[indx], c(y[indx]), tp, pp)

  w <- diff(c(0, tp[!dups]))
  sum(w*tp[!dups]/pp[!dups])/tp[length(x)]
}

compute_avgppv_noise <- function(x, y) {
  indx <- order(x+runif(length(x))/1000, decreasing = TRUE)
  tp <- cumsum(y[indx])
  pp <- seq_along(x)
  mean((tp/pp)[y[indx] == 1])
}
