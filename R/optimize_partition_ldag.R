

#' Optimize a local LDAG structure for a frequency table
#'
#' @rdname optimize_partition
#' @examples
optimize_partition_ldag <- function(counts, levels, ess, regular,
                                    min_score_improv = 0,
                                    P = as.list(1:nrow(counts)-1),
                                    labels = rep(list(integer()), length(levels)),
                                    conf = as.matrix(expand.grid(levels)),
                                    verbose = FALSE){

  r <- ncol(counts)
  q <- nrow(counts)

  nlev <- lengths(levels)
  stride <- c(1, cumprod(nlev[-length(nlev)]))
  joint  <- seq_len(nrow(conf))-1 # enumerates joint outcome
  #conf_enum_par <- conf_enum - sweep(conf, 2, stride, "*")

  # compute seach part's contribution to the local score
  partition   <- get_parts(P)
  part_size   <- lengths(P)
  part_counts <- rowsum(counts, partition, reorder = T)
  part_scores <- famscore_bdeu_byrow(part_counts, ess, r, q, part_size)

  # search for label that improves score the most
  keep_climb <- FALSE
  best_diff <- min_score_improv
  best_lab <- list()

  # check if all co-parent configs is already added to label
  for (i in seq_along(nlev)) {

    # find candidate  labels /  co-parent configs to add to label
    # - enumerate these contexts by the rows in the CPT where node i is zero
    contexts <- joint[conf[, i] == 0]

    # remove contexts already in labels on edge from i
    contexts <- contexts[match(contexts, labels[[i]], nomatch = 0) == 0]

    # rows in the CPT satisfied by each candidate label
    rows  <- outer(contexts, levels[[i]]*stride[i]+1, "+")

    # correspononding parts in the current partition
    parts <- array(partition[rows], dim(rows))
    # check
    # apply(rows, 1, function(x) cbind(conf[x, ], p = partition[x]), simplify = F)

    # find labels that implies no new independencies
    if (nlev[i] == 2) {
      redundant_labels <- parts[, 1] == parts[, 2]
    } else {
      redundant_labels <- rowSums(parts[, -1, drop = F] == parts[, 1]) == nlev[i]-1
    }


    # add redundant labels to current set of labels
    if (any(redundant_labels)) {
      if (verbose) cat("\nAdd implicit labels: node", i,
                       "contexts = ", paste(contexts[redundant_labels], collapse = ","))
      labels[[i]] <- append(labels[[i]], contexts[redundant_labels])
      if (all(redundant_labels)) next
      parts <- parts[!redundant_labels, , drop = F]
      contexts <- contexts[!redundant_labels]
    }

    # if no additional label can be added to edge from i, go to next node
    if (regular && !(length(labels[[i]]) < q/nlev[i]-1) ) next


    # find labels that implies the same independencies
    # - evaluate only the first as candidate context
    # - add redudant labels in next call to the function
    # dups  <- duplicated(lapply(split(parts, .row(dim(parts))), tabulate))

    # do not evaluate collapse parts with zero-counts
    # - collapsing rows with zero counts do not change the score
    part_zero_counts <- rowSums(part_counts) == 0
    indx <- rowSums(array(!part_zero_counts[parts], dim(parts))) > 0

    for (j in seq_along(contexts)[indx]) {

      collapse <- unique(parts[j, ])

      # compare score of collapsed vs separated parts
      tmp_counts <- colSums(part_counts[collapse,])
      tmp_size   <- sum(part_size[collapse])
      tmp_score  <- famscore_bdeu_1row(tmp_counts, ess, r, q, tmp_size)
      diff  <-  tmp_score - sum(part_scores[collapse])

      if (diff - best_diff > 10e-16) {

        # check if regular
        new_part <- list(unlist(P[collapse]))
        PP <- c(P[-collapse], new_part)

        if (regular && !is_regular(PP,  nlev, stride))  next

        best_partition <- PP
        best_diff <- diff
        best_lab  <- list(node = i, context = contexts[j], parts = collapse)
        keep_climb <- TRUE
      }
    }
  }



  if (keep_climb) {

    # update labels and partition
    labels[[best_lab$node]] <- c(labels[[best_lab$node]], best_lab$context)
    P <- best_partition

    if (verbose) cat(sprintf("\ndiff = %.5f, node %s, context = %s, parts = %s, new partition: %s",
                             best_diff,
                             best_lab$node,
                             paste(best_lab$context, collapse = ","),
                             paste(best_lab$parts, collapse = ","),
                             paste(get_parts(P), collapse = " ")))


    return(optimize_partition_ldag(counts, levels, ess,
                                   regular, min_score_improv,
                                   P, labels, conf, verbose))

  } else {

    return(list(partition = P,
                labels = labels,
                counts = part_counts,
                scores = part_scores))
  }
}
