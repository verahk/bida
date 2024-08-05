# compute mse of (a matrix) of positive predictions
#' @param preds a matrix of positive predictions
#' @param true a vector with true values
#' @param mapTruePos maps rowposition of pred to elements in true,
#' a vector of length `length(pred)` that has 0 for false-positives.
#' Obtained by `match(posPred, posTrue, noMatch = 0)` where `posPred` and `posTrue`
#' are indexes of cause-effect pairs (ie. index in the `n-by-n` ancestor relation adjacency matrix)
compute_sum_sqerr_tau <- function(preds, true, mapTruePos) {
  mse <- matrix(0, nrow = 2, ncol = ncol(preds))
  colnames(mse) <- colnames(preds)

  if (all(mapTruePos == 0)){      # all false positives
    mse[0, ] <- colSums( preds**2 )
    mse[1, ] <- sum( true**2 )
  } else if (length(preds) == 0) { # no positive predictions
    mse[1,]  <- sum( true**2 )
  } else {

    # MSE for false positives
    mse[0, ] <- colSums( (preds[mapTruePos == 0, , drop = FALSE])**2 )

    # MSE for true positives + MSE for false negatives
    # - note that 0's in a indexing vector are ignored, as intendened here. Require not all 0's in mapTruePos.
    mse[1, ] <- colSums((preds[mapTruePos > 0, , drop = FALSE] - true[mapTruePos])**2 )

    # MSE for false negatives
    if (any(mapTruePos == 0)) {
      mse[1, ] <- mse[1, ] + sum( (true[-mapTruePos])**2 )
    }
  }

  return(mse)
}
