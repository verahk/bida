#' List unique adjustment sets in a CPDAG
#'
#' @rdname adjsets_support_from_dags
#' @param cpdag (matrix)
#' @details
#' - `adjsets_support_from_cpdags`:
#'    Identifies unique adjustment sets (parent or o-set only) locally in a CPDAG,
#'    by directing siblings of every node. The same set of variables can consititute
#'    the o-set in multiple of the resulting (P)DAGs, reflected in the returned "support".
#' @export
adjsets_support_from_cpdag <- function(adjsets, cpdag) {
  n <- ncol(cpdag)

  if (all(cpdag == 0)) {
    ps <- list()
    if ("pa" %in% adjsets){
      ps$pa <- list(sets = rep(list(matrix(NA)), n),
                    support = rep(list(1), n))
    }
    if ("o" %in% adjsets){
      ps$o <-  list(sets = rep(lapply(seq_len(n), function(y) matrix(y)), each = n),
                    support = matrix(list(1), n, n))
    }
  } else {

    seqn <- seq_len(n)

    # init list for storing sets and support
    tmp <- matrix(list(), n, n)   # adjustment set dep on y
    ps <- list(pa = list(sets = tmp[, 1], support = tmp[, 1]),
               o  = list(sets = tmp, support = tmp))
    ps <- ps[adjsets]

    # find which variables that are possible causes and effects
    posspa  <- rowSums(cpdag)>0  # might be a cause
    possch  <- colSums(cpdag)>0  # might be an effect

    # o-sets for non-effects and non-causes
    if ("o" %in% adjsets) {

      tmp <- lapply(seqn, matrix)
      for (x in seqn[!posspa]) { # non-causes
        ps$o$sets[x, ] <- tmp
      }
      for (y in seqn[!possch]) {
        ps$o$sets[, y] <- tmp[y] # non-effects
      }
      ps$o$support[!posspa, ] <- ps$o$support[, !possch] <- list(1)
    }


    # adjsets for possible causes
    for (x in seqn[posspa]) {
      tmp <-  orient_siblings_in_cpdag(cpdag, x)
      pdags <- tmp$pdags
      supp  <- rep(1/length(pdags), length(pdags))

      if ("pa" %in% adjsets) {
        ps$pa$sets[x] <- list(rbind_fill(tmp$parents))
        ps$pa$support[x] <- list(supp)
      }
      if ("o" %in% adjsets) {

        sets <- array(list(), c(n, length(pdags)))
        for (g in seq_along(pdags)) {
          pdag <- pdags[[g]]
          dmat <- descendants(pdag)
          sets[, g] <- adjsets_from_dag("o", pdag, dmat = dmat, xvars = x, yvars = seqn[possch])
        }

        for (y in seqn[possch]) {
          if (y == x) next
          u <- unique(sets[y, ])
          ps$o$sets[[x, y]] <- rbind_fill(u)
          ps$o$support[[x, y]] <- c(rowsum_fast(supp, sets[y, ], u))
        }
      }
    }

    # parent sets for non-causes
    if ("pa" %in% adjsets) {
      for (x in seqn[!posspa]) {
        tmp <-  orient_siblings_in_cpdag(cpdag, x)
        pdags <- tmp$pdags
        supp  <- rep(1/length(pdags), length(pdags))

        ps$pa$sets[x, ] <- list(rbind_fill(tmp$parents))
        ps$pa$support[x, ] <- list(supp)
      }
    }
  }

  return(ps)
}



orient_siblings_in_cpdag <- function(cpdag, x){

  # replace names with colnumbers for using position rather than names with pcalg::addBgKnowledge
  colnames(cpdag) <- rownames(cpdag) <- seq(ncol(cpdag))
  gcpdag <- as(cpdag, "graphNEL")

  # apply Meek's rules and check if input is a valid CPDAG
  gcpdag <- pcalg::addBgKnowledge(gcpdag, verbose = T, checkInput = T)
  stopifnot(!is.null(gcpdag))

  # find siblings of x
  posspa  <- which(cpdag[, x] == 1)
  indx <- cpdag[x, posspa] == 1
  pa   <- posspa[!indx]             # subset with directed edges into x
  sib  <- posspa[indx]              # subset with undirected edges into x
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

    isInvalid <- vector("logical", length = length(subsib))
    for (g in seq_along(subsib)){
      # direct edges of subset s of siblings
      s    <- subsib[[g]]
      ns   <- length(s)
      from <- c(s, rep(x, nsib-ns))
      to   <- c(rep(x, ns), setdiff(sib, s))
      gprime <- pcalg::addBgKnowledge(gcpdag, from, to)
      if (is.null(gprime)) {
        isInvalid[g] <- TRUE
      } else {
        parents[[g]]  <- c(s, pa)
        pdags[[g]] <- as(gprime, "matrix")
      }
    }
    if (all(isInvalid)){
      stop("Edges can not be directed to a valid PDAG")
    }

    # remove list elements for invalid parent sets
    parents[isInvalid] <- NULL
    pdags[isInvalid] <- NULL
  }
  return(list(pdags = pdags, parents  = parents))
}

