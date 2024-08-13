
#' Title
#'
#' @param counts 
#' @param levels 
#' @param nlev 
#' @param P 
#' @param ess 
#' @param lkappa 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
optimize_partition_part <- function(counts, levels,  ess, 
                                    regular, min_score_improv = 0,
                                    P = as.list(1:nrow(counts)-1), 
                                    ldag_consistent = T, 
                                    verbose = FALSE){
  
  q   <- nrow(counts)
  r   <- ncol(counts)

  # compute seach part's contribution to the local score
  part_size   <- lengths(P)
  part_counts <- rowsum(counts, get_parts(P, part_size), reorder = T)
  part_scores <- famscore_bdeu_byrow(part_counts, ess, r, q, part_size)
  
  # find the two current parts that increase the score most if collapsed
  keep_climb <- FALSE
  if (length(P) > 1) {
    best_diff <- min_score_improv
    nlev <- lengths(levels)
    stopifnot(all(nlev == 2))
    stride <- c(1, cumprod(nlev[-length(nlev)]))
 
    
    for (i in seq.int(1, length(P)-1)) {
      for (j in seq.int(i+1, length(P))) {
        
        # check if collapsed part can be represented by labels
        new_part <- c(P[[i]], P[[j]])
   
        # compute score of collapsed parts
        tmp_count <- part_counts[i,, drop = FALSE] + part_counts[j, , drop = FALSE]
        tmp_score <- famscore_bdeu_1row(tmp_count, ess, r, q, length(new_part))
        
        # compare score with score of non-collapsed regions
        diff <- tmp_score - (part_scores[i] + part_scores[j])
        
        # compare score with current best partition
        if (diff > best_diff) {
          
          if (ldag_consistent && !is_CSI_consistent_part(new_part, levels, nlev, stride)) next 
          
          # candidate partition
          PP  <- c(P[-c(i, j)], list(new_part))
          
          # check if collapsing part i and j result in regular partition
          if (regular && !is_regular(PP, nlev, stride))  next
       
          best_diff <- diff
          best_partition <- PP
          best_parts <- c(i, j)
          keep_climb <- TRUE
        }
      }
    }
  }
  
  if (keep_climb) {
    if (verbose) cat(sprintf("\ndiff = %.5f, parts = %s, new partition: %s",
                             best_diff,
                             paste(best_parts, collapse = ","),
                             paste(get_parts(best_partition), collapse = " ")))
    
    optimize_partition_part(counts, levels,  ess, regular, min_score_improv = 0,
                            best_partition, ldag_consistent, verbose)
  } else {
    return(list(partition = P, counts = part_counts, scores = part_scores))
  }
}
