

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
#'
#'
#' #' # binary split data with missing values of outcome y
#' data <- cbind(z = rep(0:2, 9),
#'               x = rep(0:2, each = 9),
#'               y = c(rep(0, 9), rep(1, 2*9)))
#' nlev <- rep(3, 3)
#'
#' optimize_partition_from_data_tree_sparse(
#'   data,
#'   j = 3,
#'   parentnodes = 1:2,
#'   ess = 1,
#'   nlev = rep(3, 3),
#'   min_score_improv = -Inf,
#'   prune = T,
#'   verbose = TRUE
#' )
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
#'
#' levels <- list(x = 0:2, y = 0:2, z = 0:2)
#' data <- sapply(levels, sample, size = 10, replace = T)
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
  structure(result, class = "pcart")

}

extract_partitioning <- function(fit, nlev, parentnodes) {
  # extract partitioning

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
print.pcart(fit) {
  cat(fit$tree, "\n")
  cat("score: ", fit$score, "\n")
}

predict.pcart <- function(fit, newdata) {

  tree <- strsplit(fit$tree, "\n")[[1]]
  index <- grep("--", tree)
  for (part in seq_along(index)) tree[index[part]] <- gsub("--", paste0("-- ", part), tree[index[part]])


  tmp <- function(tree, newdata) {
    while (startsWith(tree[1], "-+")) {
      root <- tree[1]
      split <- strsplit(substring(root, 4), ": | \\| ")[[1]]
      var <- split[1]
      val <- split[2]
      index <- grepl("^ `", tree[-1]) == !(newdata[var] == val)
      tree <- substring(tree[-1], 3)[index]
    }
    as.numeric(substr(tree, 4, 4))
  }

  if (dim(newdata) == 0) [
    tmp(tree, newdata)
  ] else {
    apply(newdata, 1, function(x) tmp(tree, newdata = x))
  }
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

      if (nchar(rule) == 0) {
        new_rules <- add_rules
      } else {
        new_rules <- lapply(add_rules, function(x) sprintf("%s & %s", rule, x))
      }

      subtrees <- split(substring(tree[-1], 3), grepl("^ `", tree[-1]))
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


optimize_partition_from_data_tree_sparse <- function(data, j, parentnodes, ess, nlev, min_score_improv = 0, prune = F, verbose = FALSE){
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("X", seq_len(ncol(data)))
  }

  r <- nlev[j]
  q <- prod(nlev[parentnodes])

  vars   <- colnames(data)[c(j, parentnodes)]
  counts <- counts_from_data_matrix(data[, vars], nlev, TRUE)
  counts$dimnames <- setNames(lapply(counts$dim-1, seq.int, from = 0), vars)

  tmp   <- tabulate(data[, j]+1, r)
  score <- famscore_bdeu_1row(matrix(tmp, 1, r), ess, r, q)

  optimize_partition_tree_sparse(counts, ess, min_score_improv, prune, verbose)
}
