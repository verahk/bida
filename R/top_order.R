top_order <- function(dag){
  n <- ncol(dag)
  marked <- vector("logical",n)
  order <- vector("numeric",n)
  for (i in 1:n){
    for (j in 1:n){
      found <- FALSE
      if (marked[j] == FALSE && all(marked[dag[,j] == 1]) == TRUE){
        marked[j] <- TRUE
        order[i] <- j
        found <- TRUE
        break
      }
    }
    if (found == FALSE){
      stop("No topological order exists!")
    }
  }
  return(order)
}
