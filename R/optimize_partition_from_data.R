

#' Title
#'
#' @param data
#' @param j
#' @param parentnodes
#' @param nlev
#' @param ess
#'
#' @return
#' @export
#'
#' @examples
#'
#' # binary split data with missing levels
#' data <- cbind(z = rep(0:2, 2),
#'               x = rep(c(0, 2), each = 3),
#'               y = c(rep(0, 3), rep(1, 3)))
#' nlev <- rep(3, 3)
#' pcart <- optimize_partition_from_data(data, 3, 1:2, 1, nlev, "pcart", verbose = TRUE)
#' tree  <- optimize_partition_from_data(data, 3, 1:2, 1, nlev, "tree", regular = FALSE, verbose = TRUE)
#' ptree <- optimize_partition_from_data(data, 3, 1:2, 1, nlev, "ptree", verbose = TRUE)
#'
#' counts <- matrix(table(data.frame(data)), 6, 2)
optimize_partition_from_data <- function(data, j, parentnodes, ess, nlev, method, verbose = FALSE, ...) {
  if (method == "pcart") {
    # construct data frame with factor variables
    vars <- c(j, parentnodes)
    varnames <- colnames(data)[vars]
    if (is.null(varnames)) varnames <- paste0("X", vars)
    df <- data.frame(lapply(vars,
                            function(x) factor(data[, x], seq.int(0, nlev[x]-1))))
    names(df) <- varnames
    optimize_partition_from_df_pcart(df, 1, seq_along(parentnodes)+1, ess, nlev[vars], verbose = verbose, ...)
  } else {
    stride <- c(1, cumprod(nlev[parentnodes]))
    r <- nlev[j]
    q <- stride[length(parentnodes)+1]
    counts <- matrix(tabulate(data[, c(parentnodes, j)]%*%stride+1, q*r), q, r)
    optimize_partition(counts,
                       levels = lapply(nlev[parentnodes]-1, seq.int, from = 0),
                       ess,
                       method,
                       verbose = verbose,
                       ...)
  }
}
#' @rdname optimize_partition_from_data
#' @param df a data.frame of factor variables, with levels corresponding to nlev.
#' @return
#' @export
#'
#' @examples

optimize_partition_from_df_pcart <- function(df, j, parentnodes, ess, nlev, verbose = FALSE) {

  # optimize structure
  predictors <- names(df)[parentnodes]
  response   <- names(df)[j]
  if (verbose) cat("\nOptimize local structure using rpcart::opt.pcart.cat.bdeu\n")
  result <- rpcart:::opt.pcart.cat.bdeu(df, predictors, response, ess, use_structure_score = FALSE)
  if (verbose) cat(result$tree)

  # extract partitioning
  tree <- strsplit(result$tree, "\n")[[1]]
  nlev <- nlev[parentnodes]
  stride <- c(1, cumprod(nlev[-length(nlev)]))
  names(nlev) <- names(stride) <- predictors
  unlist_tree <- function(tree, nlev, stride, parts) {
    # recursive function for extracting a partition from a rpcart-tree
    if (startsWith(tree[1], "-+")) {
      root <- tree[1]
      split <- strsplit(substring(root, 4), ": | \\| ")[[1]]
      splitvar <- split[1]
      splitval <- as.numeric(strsplit(split[2], " +")[[1]])

      vals <- (parts%/%stride[splitvar])%%nlev[splitvar]
      subparts <- split(parts, vals %in% splitval)
      subtrees <- split(substring(tree[-1], 3), grepl("^ \\|", tree[-1]))
      if (length(subparts) != length(subtrees)) {
        stop()
      }
      c(unlist_tree(subtrees[[1]], nlev, stride, subparts[[1]]),
        unlist_tree(subtrees[[2]], nlev, stride, subparts[[2]]))
    } else {
      return(list(parts))
    }
  }
  partition <- unlist_tree(tree, nlev, stride, seq_len(prod(nlev))-1)
  c(partition = list(partition),
    scores = result$score,
    tree = result$tree)
}


if (FALSE) {
  rules_from_tree <- function(tree, rule = "") {
    if (startsWith(tree[1], "-+")) {
      root <- tree[1]
      split <- strsplit(substring(root, 4), ": | \\| ")[[1]]

      splitvals <- lapply(split[-1],
                          function(x) strsplit(x, " +")[[1]])
      add_rules <- lapply(splitvals,
                          function(x) paste0(paste0(split[1], "==", x), collapse = "|"))
      new_rules <- lapply(add_rules,
                          function(x) paste0(rule, x, collapse = "&"))

      subtrees <- split(substring(tree[-1], 3), grepl("^ \\|", tree[-1]))

      c(rules_from_tree(subtrees[[1]], new_rules[[1]]),
        rules_from_tree(subtrees[[2]], new_rules[[2]]))
    } else {
      list(rule)
    }
  }

  rules_chr  <- rules_from_tree(tree, "")
  rules_expr <- lapply(rules_chr, function(x) parse(text = x))
  with(tail(df), eval(rules_expr[[1]]))
  lapply(rules_expr, function(xx) which(with(df, eval(xx))))

}
