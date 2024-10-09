sim_load_bn <- function(par) {
  if (is.null(par$bnname)) {
    n <- par$n
    k <- par$k
    nlev <- rep(k, n)
    bida:::rand_bn(n = par$n,
            d = 8,
            type = "cat",
            nlev = nlev,
            local_structure = "tree",
            prob = 1,
            regular = TRUE,
            maxdepth = par$maxdepth)
  } else {
    readRDS(paste0("inst/data/", par$bnname, ".rds"))
  }
}


