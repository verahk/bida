sim_load_bn <- function(par) {
  if (is.null(par$bnname)) {
    n <- par$n
    k <- par$k
    nlev <- rep(k, n)
    if (par$maxdepth < 1) {
      bida:::rand_bn(n = par$n,
                     d = 6,
                     type = "cat",
                     nlev = nlev,
                     local_structure = "tree",
                     prob = 1,
                     regular = TRUE,
                     maxdepth = par$maxdepth)
    } else {
      bida:::rand_bn(n = par$n,
                     d = 6,
                     type = "cat",
                     nlev = nlev)
    }


    # profvis::profvis({
    #   set.seed(par$r)
    #   dag <- bida:::rand_dag(n, d = 6)
    #   dist <- bida:::rand_dist(dag, "cat", nlev = nlev, local_structure = "tree", prob = 1, regular = TRUE, maxdepth = 0)
    # })

    # profvis::profvis({
    #   bida:::rand_bn(n = par$n,
    #                  d = 6,
    #                  type = "cat",
    #                  nlev = nlev)
    # })
    #bida:::custom_bn(dag, dist, FALSE)
  } else {
    readRDS(paste0("inst/data/", par$bnname, ".rds"))
  }
}


