


# test code to see what takes so long ...
indir <- "./inst/simulations/sim_bn/MCMCchains/"
f <- "asia_pcskel_ptree_order_ess1_epf2_N300_r16.rds"

# import DAGs
filepath <- paste0(indir, f)
res <- readRDS(filepath)
par <- res$par
lookup <- res$lookup
MCMCchain <- res$MCMCchain
dags <- MCMCchain$traceadd$incidence

# check how many (x, z) sets are in lookup
ps <- bida:::parent_support_from_dags(dags)

n <- length(ps$sets)
res <- matrix(NA, n, 3, dimnames = list(node = 1:n, metric = c("sets", "size", "inLookup")))
for (x in seq_len(n)) {

  eval <- function(z) {
    z <- z[!is.na(z)]
    parIDs <- paste(seq_len(n)[-x], paste0(sort(c(x, z)), collapse = "."), sep = ".")
    inLookup <- !is.na(match(parIDs, names(lookup$ptree)))
    c(length(z), sum(inLookup))
  }

  nsets <- nrow(ps$sets[[x]])
  res[x, ] <- c(nsets, rowSums(apply(ps$sets[[x]], 1, eval))/nsets)
}
