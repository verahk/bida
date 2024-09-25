

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
#' # binary split data
#' data <- cbind(z = rep(0:2, 9),
#'               x = rep(0:2, each = 9),
#'               y = c(rep(0, 9), rep(1, 2*9)))
#' nlev <- rep(3, 3)
#'
#' # optimize local structure
#' res <- optimize_partition_from_data_pcart(data, 3, 1:2, 1, nlev, verbose = TRUE)
#' parts <- get_parts(res$partition)
#' stopifnot(all(match(parts, unique(parts)) == c(rep(1, 3), rep(2, 6))))
#'
#' # score do not match when there are unobserved levels of outcome variable
#' res[-1]
#'
#' res <- optimize_partition_from_data_pcart(data, 3, 1:2, 1, nlev = NULL, recompute_score = TRUE, verbose = TRUE)
#' res[-1]
#' stopifnot(abs(res$score-res$pcart$score) < 10**-10)
optimize_partition_from_data <- function(data, j, parentnodes, ess, nlev, method, verbose = FALSE, ...) {
  if (method == "pcart") {
    optimize_partition_from_data_pcart(data, j, parentnodes, ess, nlev, verbose = verbose, ...)
  } else {
    stride <- c(1, cumprod(nlev[parentnodes]))
    r <- nlev[j]
    q <- stride[length(parentnodes)+1]
    counts <- matrix(tabulate(data[, c(parentnodes, j)]%*%stride+1, q*r), q, r)
    optimize_partition(counts,
                       levels = lapply(nlev[parentnodes]-1, seq.int, from = 0),
                       ess,
                       method,
                       regular = (method == "ldag"),
                       verbose = verbose,
                       ...)
  }
}
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
#' # binary split data
#' data <- cbind(z = rep(0:2, 9),
#'               x = rep(0:2, each = 9),
#'               y = c(rep(0, 9), rep(1, 2*9)))
#' nlev <- rep(3, 3)
#'
#' # optimize local structure
#' res <- optimize_partition_from_data_pcart(data, 3, 1:2, 1, nlev, verbose = TRUE)
#' parts <- get_parts(res$partition)
#' stopifnot(all(match(parts, unique(parts)) == c(rep(1, 3), rep(2, 6))))
#'
#' # score do not match when there are unobserved levels of outcome variable
#' res[-1]
#'
#' res <- optimize_partition_from_data_pcart(data, 3, 1:2, 1, nlev = NULL, recompute_score = TRUE, verbose = TRUE)
#' res[-1]
#' stopifnot(abs(res$score-res$pcart$score) < 10**-10)
optimize_partition_from_data_pcart <- function(data, j, parentnodes, ess, nlev = NULL, recompute_score = !is.null(nlev), verbose = FALSE) {

  if (is.null(nlev)) nlev <- apply(data, 2, function(x) length(unique(x)))

  r <- nlev[j]
  q <- prod(nlev[parentnodes])

  # optimize structure
  df <- as.data.frame(data)
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

  # compute score
  if (recompute_score) {
    n_parts <- length(partition)
    data_parts <- get_parts(partition)[data[, parentnodes]%*%stride+1]
    counts <- matrix(tabulate(data_parts + n_parts*data[, j], n_parts*r), n_parts, r)
    score  <- sum(famscore_bdeu_byrow(counts, ess, r, q, lengths(partition)))
  } else {
    score <- result$score
  }

  list(partition = partition,
       score = score,
       pcart = result)
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

