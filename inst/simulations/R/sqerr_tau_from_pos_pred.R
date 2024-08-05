
sqerr_tau_from_pos_pred <- function(x, true, mapTruePos) {
  # x a matrix with positive predictions
  # true a matrix with true vals
  # mapTruePos a integer vector mapping rows in x to rows in true
  sqerr <- matrix(0, nrow = 2, ncol = ncol(true))

  if (length(mapTruePos) == 0) { # no positive predictions
    # error of false negatives
    sqerr[2, ] <-  colSums( true**2 )
  } else if (all(mapTruePos == 0)) { # no correct positive predictions

    # error of false positives
    sqerr[1, ] <- colSums( x*2 )

    # error of false negatives
    sqerr[2, ] <-  colSums( true**2 )

  } else {

    # error of true positives
    # - note that 0's in a indexing vector are ignored, as intendened here. Require not all 0's in mapTruePos.
    sqerr[2, ] <- colSums((x[mapTruePos > 0, , drop = FALSE] - true[mapTruePos,,drop = FALSE])**2 )

    # error of false positives
    sqerr[1, ] <- colSums( (x[mapTruePos == 0, , drop = FALSE])**2 )

    # error of false negatives
    if (any(mapTruePos == 0)) {
      sqerr[2, ] <-  sqerr[2, ] + colSums( true[-mapTruePos, , drop = FALSE]**2 )
    }
  }
  return(sqerr)
}
