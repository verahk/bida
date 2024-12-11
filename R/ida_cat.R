
#' Title
#'
#' @param cpdag (integer matrix) an adjacency matrix representing a CPDAG, where
#'
#' - `cpdag[i, j] = 1` indicates an edge from node `i` to node `j`,
#' - `cpdag[i, j] = 1` indicates absence of an edge from node `i` to node `j`,
#' - `cpdag[i, j] = cpdag[j, i] = 1` indicates an undirected edge from node `i` to node `j`.
#' @param data (integer matrix) data matrix representing categorical variables. Assumed to be zero-based.
#' @param adjset (character) which adjustment set to use in backdoor formula.
#' @param params (list) list with hyperparameters for the local bdeu-priors.
#' @param x,y (integer vector) position of cause and effect variables, respectively.
#' @param verbose (logical) print output to follow progression.
#'
#' @return a bida-object.
#' @export
#'
#' @examples
#'
#' data(bida_example_cat)
#' attach(bida_example_cat)
#'
#' # infer a CPDAG with pcalg::pc
#' suffStat <- list(dm = data,
#'                  nlev = apply(data, 2, max)+1,
#'                  adaptDF = FALSE)
#' fit <- pcalg::pc(suffStat,
#'                  pcalg::disCItest,
#'                  alpha = .05,
#'                  labels = colnames(data),
#'                  u2pd = c("retry"),
#'                  verbose = FALSE)
#' cpdag <- as(fit@graph, "matrix")
#' # Rgraphviz::plot(fit@graph)
#' cpdag
#'
#' # estimate intervention probabilities
#' params <- list(ess = 1, nlev = apply(data, 2, max)+1)
#' fit <- ida_cat(cpdag, data, "pa", params, contrasts = list(jsd = jsd))
#' posterior_mean()
#'
#' # optimal-IDA, set effects to zero exactly
#' res <- ida_cat(cpdag, data, "o", params, contrasts = list(jsd = jsd))
#' res
ida_cat <- function(cpdag,
                data,
                adjset = c("pa", "o"),
                params,
                x = 1:ncol(data),
                y = 1:ncol(data),
                contrasts = NULL,
                verbose = FALSE) {

  adjset <- match.arg(adjset, c("pa", "o"))
  mode(data) <- mode(params$nlev) <- "integer"

  n <- ncol(data)
  out <- matrix(list(), n, n)
  for (xx in x) {
    if (verbose) cat("Estimating unique intervention distributions for cause", xx, "\n")
    tmp <- orient_siblings_in_cpdag(cpdag, xx)
    if (adjset == "pa") {
      ps <- list(tmp$parents, tmp$p)
      out[xx, ] <- bida_posterior_cat(ps, data, xx, y, ess = params$ess, nlev = params$nlev)
    } else {
      dag_supp <- dag_support(tmp$pdags)
      posteriors <- vector("list", ncol(data))
      for (yy in y[-match(xx, y)]) {
        ps <- adjset_support(dag_supp, xx, yy, adjset)
        out[[xx, yy]] <- bida_posterior_cat(ps, data, xx, yy, ess = params$ess, nlev = params$nlev)[[yy]]
      }
    }
  }
  out
}

orient_siblings_in_cpdag <- function(cpdag, x){

  tcpdag <- t(cpdag)
  stopifnot(pcalg::isValidGraph(tcpdag, type = "pdag"))

  # find siblings of x
  posspa  <- which(cpdag[, x] == 1)
  indx <- cpdag[x, posspa] == 1
  sib  <- posspa[indx]              # subset with undirected edges into x
  pa   <- posspa[!indx]             # subset with directed edges into x
  nsib <- length(sib)

  if (nsib == 0){
    parents  <- list(pa)
    pdags <- list(cpdag)
  } else {
    ## list all subsets of siblings = possible additional parents
    if (length(sib) == 1) {
      subsib <- list(vector("integer", 0), sib)
      parents <- pdags <- vector("list", 2)
    } else {
      subsib  <- unlist(lapply(c(0, seq_along(sib)), combn, x = sib, simplify = FALSE), recursive = FALSE)
      parents <- pdags <- vector("list", length(subsib))
    }

    # make graphNEL-object for pcalg::addBgKnowledge
    colnames(cpdag) <- rownames(cpdag) <- seq_len(ncol(cpdag))

    isInvalid <- vector("logical", length = length(subsib))
    for (g in seq_along(subsib)){
      # direct edges of subset s of siblings
      s    <- subsib[[g]]
      ns   <- length(s)
      from <- c(s, rep(x, nsib-ns))
      to   <- c(rep(x, ns), setdiff(sib, s))
      gprime <- pcalg::addBgKnowledge(tcpdag, from, to)
      if (is.null(gprime)) {
        isInvalid[g] <- TRUE
      } else {
        parents[[g]]  <- c(s, pa)
        pdags[[g]] <- t(gprime)
      }
    }
    if (all(isInvalid)){
      stop("Edges can not be directed to a valid PDAG")
    }

    # remove list elements for invalid parent sets
    parents[isInvalid] <- NULL
    pdags[isInvalid] <- NULL
  }
  return(list(pdags = pdags, parents  = parents, p = rep(1/length(parents), length(parents))))
}

