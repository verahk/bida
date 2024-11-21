

#' Class: bida
#'
#' An object of class bida is a matrix list with objects of class [bida::bida_posterior],
#' that describes the posterior mixture distribution over the causal effect between every variable pair.
#' The `posterior_mean` and `posterior_sample` methods calls the respective methods
#' for each [bida::bida_posterior] and collect the output in a matrix list.
#'
#' @param ps (list)
#' Support over adjustment sets, stored in a list with two named elements `sets` and `support`.
#' See [adjsets_posterior]
#' @param data (matrix)
#' @param params additional params
#' @param pairs (matrix) a n-by-n matrix indicating which pairs should be considered.
#'  If NULL (default), all pairs are considered and `pairs[i, j] = 1` for all but `i == j`.
#' @param x an object of class [bida]
#' @param samplesize (integer) samplesize
#' @seealso [bida_pair_cat]
#' @return
#' `bida`: a n-by-n list where each element is an object of class `bida_pair`
#' @examples
#'
#' dag  <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
#' colnames(dag) <- rownames(dag) <- letters[1:3]
#' dags <- list(dag, t(dag))
#'
#' # Categorical data
#' nlev <- 2:4
#' cpts <- rand_dist(dag, "cat", nlev = nlev)
#' bn   <- custom_bn(dag, cpts)
#' data <- sample_data_from_bn(bn, 10**4)
#' head(data)
#'
#'
#' params <- list(nlev = nlev, ess = 1)
#' bida <- bida(list(dag, t(dag)), data, "cat", "pa", params)
#' bida
#'
#' # posterior mean
#' means <- posterior_mean(bida)
#' means
#'
#' # posterior mean of intervention effect (jsd)
#' means <- posterior_mean(bida, contrasts = list(jsd = jsd))
#' means
bida <- function(dags,
                 data,
                 type = "cat",
                 adjset = "pa",
                 params,
                 x = 1:ncol(data),
                 y = 1:ncol(data),
                 replace_large_adjset = NULL,
                 verbose = FALSE) {

  type <- match.arg(type, c("categorical"))
  n <- ncol(data)

  out <- matrix(list(), n, n)
  colnames(out) <- rownames(out) <- colnames(data)

  # prep
  if (type == "categorical") {
    if (is.null(params$ess)) ess <- 1 else ess <- params$ess
    if (is.null(params$nlev)) nlev <- apply(data, 2, max+1) else nlev <- nlev

    tmp <- lapply(y, function(yy) bdeu_posterior(data, yy, integer(0), ess, nlev))
    ays <- replace(list(), y, tmp)
  }

  if (adjset == "pa") {
    tic <- c()
    tic[1] <- Sys.time()
    if (missing(dags)) {
      stopifnot(adjset == "pa")
      stop("exact method for computing parent posterior is not available.")
    } else {
      dag_supp <- dag_support(dags)
      ps <- lapply(x,
                   function(xx) adjset_support(dag_supp, xx, 0L, "pa",
                                               replace_large_adjset = replace_large_adjset))
    }
    tic["adjset"] <- Sys.time()
    for (xx in x) {
      out[xx, ] <- bida_posterior_cat_local(ps[[xx]], data, xx, y, ess, nlev, ays)
    }
    tic["params"] <- Sys.time()
    toc <- diff(tic)
  } else {
    dag_supp <- dag_support(dags)
    dmats <- lapply(dag_supp[[1]], descendants)
    tic <- array(NA, c(n, n, 3),
                 dimnames = list(NULL, NULL, c("start", "adjset", "params")))
    for (xx in x) {
      for (yy in y[match(y, xx, 0L) == 0]) {
        tic[[xx, yy, 1]] <- Sys.time()

        ps <- adjset_support(dag_supp, xx, yy, adjset, dmats, replace_large_adjset)
        tic[[xx, yy, 2]] <- Sys.time()

        out[[xx, yy]] <- switch(type,
                                "categorical" = bida_posterior_cat(ps, data, xx, yy, ess, nlev, ays))
        tic[[xx, yy, 3]] <- Sys.time()
      }
    }
    toc <- colSums(tic[,,-1]-tic[,,-3], dims = 2, na.rm = TRUE)

  }
  structure(out,
            toc   = toc,
            type = type,
            class = c("bida"))
}

#' @export
print.bida <- function(obj) {
  cat("\nbida-object of type", attr(obj, "type"))
  cat("\nVariables: ", colnames(obj), "\n")
  cat("Runtimes in seconds:\n")
  print(attr(obj, "toc"))
}


#' @rdname bida
#' @export
posterior_mean.bida <- function(obj, ...) {
  out   <- array(list(), dim(obj), dimnames(obj))
  out[] <- lapply(obj, function(v) posterior_mean(v, ...))
  return(out)
}

#' @rdname bida
#' @export
posterior_sample.bida <- function(obj, n, ...) {
  out  <- array(list(), dim(obj), dimnames(obj))
  out[] <- lapply(obj, function(v) posterior_sample(v, n, ...))
  return(out)
}
