

#' Compute precision recall
#'
#' @param x score
#' @param y outcome
#' @param method which method to use
#'
#' @return
#' @keywords internal
#' @examples
#'
#' # no ties
#' compare_prec_recall(1:0, 1:0, niter = 10)
#' compare_prec_recall(0:1, 1:0, niter = 10)
#'
#' # perfectly right
#' x <- c(2, 2, 1, 0)
#' y <- c(1, 1, 1, 0)
#' compare_prec_recall(x, y, niter = 20)
#'
#' # perfectly wrong
#' x <- c(0, 0, 1, 1)
#' y <- c(1, 1, 0, 0)
#' compare_prec_recall(x, y, niter = 20)
#'
#'
#' #' # get first, miss last
#' x <- c(1, 1, 0)
#' y <- c(1, 0, 1)
#' compare_prec_recall(x, y, niter = 20)
#'
#' # get first, miss last, no ties
#' x <- c(2, 1, 0)
#' y <- c(1, 0, 1)
#' compare_prec_recall(x, y, niter = 20)
#'
#' x <- c(8:1, rep(0, 2))
#' y <- c(x[x > 0]/10 > runif(8), .2 > runif(2))
#' compare_prec_recall(x, y, niter = 20)
#'
#' x <- c(80:1, rep(0, 20))
#' y <- c(x[x > 0]/100 > runif(80), .2 > runif(20))
#' compare_prec_recall(x, y, niter = 20)
#'
#'
#'


compute_prec_recall <- function(x, y, method) {
  switch(method,
         "noise" = compute_prec_recall_noise(x, y),
         "dups"  = compute_prec_recall_dups(x, y),
         "interp" = compute_prec_recall_interp(x, y))
}

compute_prec_recall_noise <- function(x, y) {
  indx <- order(x+runif(x)/10**5, decreasing = TRUE)
  tp   <- cumsum(y[indx])
  list(df = data.frame(x = x[indx], TPR = tp/sum(y), PPV = tp/seq_along(x)),
       avgppv = mean((tp/seq_along(x))[y[indx]>0]))
}

compute_prec_recall_dups <- function(x, y) {

  indx <- order(x, decreasing = TRUE)
  tp   <- cumsum(y[indx])
  pp   <- seq_along(x)
  n1   <- tp[length(y)]

  dups <- duplicated(x[indx], fromLast = TRUE)
  #cbind(x = x[indx], dups, y[indx], tp, pp, ppv = tp/pp)

  PPV  <- tp[!dups]/pp[!dups]
  w    <- diff(c(0, tp[!dups]))
  list(df = data.frame(x = x[indx][!dups],
                       TPR = tp[!dups]/n1,
                       PPV = PPV),
       avgppv = 1/n1*sum(w*PPV))
}

compute_prec_recall_interp <- function(x, y) {

  if (all(x[y == 1]  == max(x))) {
    # -approx- fails if there are not more than two unique values
    # hardcode the special case where the highest score is given to all positives
    # e.g. x = c(1, 1, 0) = y

    n1  <- sum(y)
    TPR <- PPV <- rep(1, n1)
  } else {

    indx <- order(x, decreasing = TRUE)

    # compute cumulative predictions, aggregated over duplicated scores
    dups <- duplicated(x[indx], fromLast = TRUE)
    pp   <- seq_along(x)[!dups]     # number of positive predictions
    tp   <- cumsum(y[indx])[!dups]  # true positives
    fp   <- pp-tp
    n1   <- sum(y)
    # cbind(x = x[indx], y = y[indx], dups, tp = cumsum(y[indx]), fp = cumsum(!y[indx]), pp = seq_along(x))

    # create vectors for storing the cumulative TP and FP at each true positive
    # - for computing average PPV, only the positives contributes.
    new_tp <- new_fp <- seq_len(n1)

    # fill in false positives at each
    indx   <- match(new_tp, tp, 0L) # identify which counts are in tp
    new_fp[tp[indx]] <- fp[indx]        # include false positives at each observed TP

    if (any(indx > 0)) {
      # interpolate false positives at "missing" true positives
      # - that is where the cumulative count of true positives
      #   increases by more than 1 due to duplicated scores
      # - pad with zeros for interpolation at tp < min(tp)
      tmp <- approx(c(0, tp), c(0, fp),
                    xout = new_tp[indx == 0],
                    ties = "max")
      new_fp[indx == 0] <- tmp$y
    }


    #plot(tp, fp, xlim = c(0, n1))
    #lines(new_tp, new_fp)

    TPR <- new_tp/n1
    PPV <- new_tp/(new_tp+new_fp)
  }

  list(df = data.frame(TPR = TPR, PPV = PPV),
       avgppv = 1/n1*sum(PPV))
}

compare_prec_recall <- function(x, y, niter = 10) {
  avgppvs <- numeric(niter)

  res <- compute_prec_recall(x, y, method = "dups")
  plot(res$df$TPR, res$df$PPV,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "recall", ylab = "precision")
  avgppvs[1] <- res$avgppv

  res <- compute_prec_recall(x, y, method = "interp")
  points(res$df$TPR, res$df$PPV, col = "blue", pch = 4)
  lines(res$df$TPR, res$df$PPV, col = "blue", lty = "dotted")
  avgppvs[2] <- res$avgppv


  for (i in seq_len(niter)) {
    res <- compute_prec_recall(x, y, method = "noise")
    lines(res$df$TPR, res$df$PPV, col = "red", type = "s")
    avgppvs[i+2] <- res$avgppv
    cat("\navgppv, noise, iteration:", i, ":", res$avgppv)

  }

  cat("\navgppv, dups:", avgppvs[1])
  cat("\navgppv, interp:", avgppvs[2])
  cat("\navgppv, noise, mean:",   mean(avgppvs[-c(1:2)]))
}

