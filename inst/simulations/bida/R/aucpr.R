


auc_trapzoid <- function(x, y){
  dx <- base::diff(x)
  dy <- base::diff(y)

  sum(dx*y[-1]-dx*dy/2)
}

agg_pred_over_ties <- function(x, y) {
  if (anyDuplicated(x)) {
    # sum number of predictions at each unique value of x
    # - return a matrix with rows sorted from lowest to highest value of x
    tmp  <- unname(rowsum(cbind(1, y), x, reorder = T))
    indx <- seq.int(nrow(tmp), 1)
    n   <- tmp[indx, 1]
    y   <- tmp[indx, 2]

    # compute true and false positive rates
    pp  <- cumsum(n)
    tp  <- cumsum(y)
  } else {
    n   <- rep(1, length(x))
    y   <- y[order(x, decreasing = T)]
    pp  <- seq_along(x)
    tp  <- cumsum(y)
  }
  return(data.frame(n = n, y = y, pp = pp, tp = tp))
}

auc_roc <- function(x, y, method = c("rank", "trapzoid")){
  n = length(x)
  n1 = sum(y)
  n0 <- n-n1
  if (method == "rank") {
    # compute AUC as prob ranking a random non-event higher than event
    r <- rank(x)
    return((sum(r[y == 1])-n1*(n1+1)/2)/(n0*n1))
  } else if (method == "trapzoid") {
    tmp <- agg_pred_over_ties(x, y)
    tpr <- tmp$tp/n1
    fpr <- (tmp$pp-tmp$tp)/n0
    return(auc_trapzoid(fpr, tpr))
  }
}


aucpr <- function(x, y, n = length(x), n1 = sum(y), method = "avgppv") {

  if (length(x) == 0 || all(y == 0)) {
    switch(method,
           "avgppv" = n1/n,
           "trapezoidal" = 1/2*(1 + n1/n))
  } else {

    # aggregate predictions
    tmp <- agg_pred_over_ties(x, y)

    if (method == "avgppv") {

      # maximal number of positive prediction given x
      s <- tmp$tp[length(tmp$tp)]

      # return average precision
      # - adjustment term for zero-effects:
      #   (n1-s) = number of true positives assigned zero score
      #   (n1-s)/n1 * n1/n = weight * precision
      wppv <- with(tmp[tmp$y > 0, ], tp/pp*y)
      sum(wppv)/n1 + (n1-s)/n

    } else if (method == "avgppvr") {
      rank <- rank(x, ties.method = "average")
      ord  <- order(rank, decreasing = T)

      ppv  <- seq_len(sum(y))/rank[ord[y>0]]
    } else if (method == "interpolate") {
      stop()
      n1 <- sum(y)
      n0 <- length(y)-n1
      if (any(tmp$y) > 2) {
        tmp <- tmp[tmp$y>0,]
        fp <- (tmp$pp-tmp$tp)
        fpr <- fp/n0
        tpr <- tmp$tp/n1
        slope <- fpr/tpr

        # interpolate false positives
        fphat <- numeric(n1+1)
        k  <- 2
        for (i in seq_along(tmp$y)) {
          for (j in seq_along(tmp$y)) {
            fphat[k] <- fphat[k-1] + slope[i]*j
            k <- k +1
          }
        }
      }
    } else if (method == "trapezoidal") {
      # aggregate true positives
      # - aggregate
      if (length(x) > length(unique(x))) {
        tmp <- unname(rowsum(cbind(1, y), x, reorder = T))
        yy  <- rev(tmp[, 2])   # number of true postives, ordered by x
        np  <- rev(tmp[, 1])   # number of predictions, ordered by x
      } else {
        yy  <- y[order(x, decreasing = T)]
        np  <- rep(1, length(x))
      }

      tp <- cumsum(yy)         # cumulative true positives
      pp <- cumsum(np)         # cumulative positive predictions
      ppv  <- tp/pp
      s    <- tp[length(tp)]
      1/(2*n1)*(sum((yy>0)*(c(1, ppv[-length(ppv)])+ppv)) + (n1-s)*(s/length(x)+n1/n))
    }
  }
}

if (FALSE) {
  # perfectly correct
  x <- 1:3
  auc_pr(x, y = c(0, 0, 1), method = "avgppv") == 1
  auc_pr(x, y = c(0, 0, 1), method = "trapezoidal") == 1


  # the "ties at zero" option do not affect auc if there are no true positives among scores
  auc_pr(x, y = c(0, 0, 1), n = 10, n1 = 1, method = "avgppv") == 1
  auc_pr(x, y = c(0, 0, 1), n = 10, n1 = 1, method = "trapezoidal") == 1

  auc_pr(x, y = c(0, 0, 1), n = 10, n1 = 2, method = "avgppv") == 1/2 + (2-1)/10
  auc_pr(x, y = c(0, 0, 1), n = 10, n1 = 2, method = "trapezoidal")

  # all ties
  x <- rep(0, 3)
  auc_pr(x, y = c(0, 0, 1), method = "avgppv") == 1/3
  auc_pr(x, y = c(0, 0, 1), method = "trapezoidal") == 1

  # perfectly incorrect
  x <- 3:1
  auc_pr(x, y = c(0, 0, 1), method = "avgppv") == 1/3
  auc_pr(x, y = c(0, 0, 1), method = "trapezoidal") == 1/2*(1+1/3)

  # all false positive estimates
  x <- 1:3
  auc_pr(x, y = c(0, 0, 0), n = 10, n1 = 2, method = "avgppv") ==  (2-1)/10
  auc_pr(x, y = c(0, 0, 0), n = 10, n1 = 2, method = "trapezoidal")


  aucpr(x, c(0, 0, 1, 1)) == 1
  aucpr(x, c(0, 0, 1, 1), n = 10, n1 = 2) == 1
  aucpr(x, c(0, 0, 1, 0)) == .75
  aucpr(1:4, c(0, 0, 1, 0)) == (1/3+2/3*1/2) + .5*(1/4+(1/3-1/4)/2)

  avgppv(x, c(0, 0, 1, 1)) == 1
  avgppv(x, c(0, 0, 1, 1), n = 10, n1 = 2) == 1
  avgppv(x, c(0, 0, 1, 0)) == .5
  avgppv(1:4, c(0, 0, 1, 0)) == 1*1/4

  x <- c(1, 2, 3, 3) + runif(4)/100
  aucpr(x, c(0, 0, 1, 0))
  avgppv(x, c(0, 0, 1, 0))

  n <- 10
  n1 <- 2
  aucpr(x, c(0, 0, 1, 0), n, n1) == .5*(1+.5)/2+.5*(1/4+n1/n)/2

  x <- c(1, 0, 0, 0)
  y <- c(0, 1, 0, 0)
  aucpr(x, y)
}


