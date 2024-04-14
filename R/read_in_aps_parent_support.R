# Computes binary BIDA posterior for some given parent support.
# ------------------------------------------------------------------------------
# INPUT:  - file_path: path to a parent support file in APS format.
#         - max_parent_size: limit on the maximum parent size used when
#                            computing the family scores
# OUTPUT: ps: list containing parent supports for each node.
# ------------------------------------------------------------------------------
read_in_aps_parent_support <- function(file_path, max_parent_size){
  temp <- utils::read.table(file_path,
                            header = FALSE,
                            sep = " ",
                            col.names = paste0("V",seq_len(max_parent_size+2)),
                            fill = TRUE)
  d <- temp[1,1]
  pos <- 2
  parent <- vector("list", length = d)
  parent_support <- vector("list", length = d)
  for (i in 1:d){
    node <- temp[pos,1]
    npar <- temp[pos,2]
    pos <- pos+1
    parent[[node]] <- as.matrix(temp[pos:(pos+npar-1),3:(max_parent_size+2)])
    s <- temp[pos:(pos+npar-1),1]
    s <- exp(s-max(s))
    parent_support[[node]] <- s/sum(s)
    pos <- pos+npar
  }
  ps <- list()
  ps$sets <- parent
  ps$support <- parent_support
  return(ps)
}
