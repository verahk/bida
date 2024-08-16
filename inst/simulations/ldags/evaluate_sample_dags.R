

dir_in <- "./inst/simulations/ldags/MCMCchains/"

# import files
files <- list.files(dir_in, "*.rds", full.names = T)[1:100]
imp <- lapply(files, readRDS)

compute_prec_recall <- function(x, y) {
  indx <- order(x+runif(x)/10**5, decreasing = TRUE)
  tp   <- cumsum(y[indx])
  list(df = data.frame(x = x[indx], TPR = tp/sum(y), PPV = tp/seq_along(x)),
       avgppv = mean((tp/seq_along(x))[y[indx]>0]))
}

bn <- readRDS("./data/asia.rds")
ground_truth <- list(dag = bnlearn::amat(bn),
                     dmat = descendants(bn))

eval <- function(MCMCchain, lookup, burnin = seq_len(200)) {

  dag <- ground_truth$dag
  dmat <- ground_truth$dmat
  n <- ncol(dag)
  seqn <- seq_len(n)             # for iterating over every node
  dindx <- diag(n) == 1          # for excluding diagonal elements in adj matrices

  dags <- lapply(MCMCchain$traceadd$incidence[-burnin], as.matrix)
  support <- rep(1/length(dags), length(dags))

  # list unique dags, to not duplicate operations
  tmp <- unique(dags)
  support <- rowsum_fast(support, dags, tmp)
  dags <- unique(dags)

  # POSTERIOR EDGE PROBABILITIES ----
  dmats <- lapply(dags, descendants)

  # compute posterior edge prob for DAG and ancestor relations
  compute_edgeps <- function(amat, support) Reduce("+", Map("*", amat, support))
  edgeps <- list(dag = compute_edgeps(dags, support)[!dindx],
                 dmat = compute_edgeps(dmats, support)[!dindx])

  # compute prec_recall
  pr <- Map(compute_prec_recall, edgeps, list(dag[!dindx], dmat[!dindx]))

  # INTERVENTION DISTRIBUTION ----

  #

  # list all unique parent set of every variable
  adjsets <- parent_support_from_dags(dags, support, NULL)
  names_lookup <- names(lookup)
  for (x in seqn) {
    means

    }
  }
}


bida_pairs_global <- function(data, x, y, sets, support, par) {

  out <- matrix(list(), nrow = length(support), length(y))
  for (yy in y) {
  # indicator for zero effects
  isZeroEffect <- lapply(y, function(yy) rowSums(yy == sets, na.rm = TRUE))

  for (r in seq_along(support)) {
    z <- sets[r, ]
    z <- z[is.na(z)]

      out[[r, y]] <- backdoor_params(data, x, y, z, par)
    }
  }
  return
}
bida_bdeu_lookup <- function(data, x, y, sets, support, nlev, ess = 1, lookup, method, ...) {

    out <- matrix(list(), nrow = length(support), length(y))
    for (r in seq_along(support)) {
      z <- sets[r, ]
      z <- z[is.na(z)]
      for (yy in y) {
        if (any(yy == z)) {
          out[[r]] <- bida_bdeu(data, j, integer(0L), ess, nlev)
        } else {
          parID <- paste(c(y, sort(c(x, z))), collapse = ".")
          if (parID %in% names_lookup) {
            bdeu <- lookup[[parID]]$bdeu
          } else {
            # !! req data and params as arguments
            bdeu <- bida_bdeu(data, j, parentnodes, ess, nlev, NULL)
            opt  <- optimize_bdeu(bdeu, method, levels)
            bdeu$partition <- opt$partition
          }

          pdo_mean <- backdoor_mean(bdeu,)
        }
      }
    }

